import igraph
import os
import pandas as pd
import matplotlib.pyplot as plt
import itertools
from igraph import plot
from icecream import ic
import numpy as np
import scipy.stats as stats

def read_num_results() -> pd.DataFrame:
  #read the numerical results from the results folder:
  directory = "../results/"
  ic("reading files from " + directory)
  data = pd.DataFrame()
  for filename in os.listdir(directory):
    ic("reading file " + filename)
    if filename.endswith(".num"):
      with open(os.path.join(directory, filename), "r") as file:
        # The name of the file is of the form name_XX.num. Extract XX, which is an interger:
        name = filename.split(".")[0]
        name = name.split("_")[1]
        X = int(name)
        # Read the file into a temporary dataframe:
        temp = pd.read_csv(file, header=None)
        # X is the index of the dataframe:
        temp["X"] = X
        # Add the temporary dataframe to the data:
        data = data.append(temp)
  # Rename the columns:
  data.columns = ["num_detached", "epoi_cut", "epoi_cut_base", "epoi_greedy", "epoi_greedy_base", "X"]
  # Set the index to X:
  data.set_index("X", inplace=True)
  # Sort the data by X:
  data.sort_index(inplace=True)
  # Calculate the absolute value of the difference between epoi_cut_base and epoi_greedy_base, and name that into base_epoi_dev:
  data["base_epoi_dev"] = abs(data["epoi_cut_base"] - data["epoi_greedy_base"])
  # Combine epoi_cut_base and epoi_greedy_base into one column, as the averge of the two:
  data["epoi_base"] = (data["epoi_cut_base"] + data["epoi_greedy_base"])/2
  # Remove the columns epoi_cut_base and epoi_greedy_base:
  data.drop(columns=["epoi_cut_base", "epoi_greedy_base"], inplace=True)
  # return the dataframe:
  return data
        
## Create a comparitive plot using the dataframe: 
def comp_plot(data:pd.DataFrame) -> None:
  ## The dataframe has X as the index and the columns are num_detached, epoi_cut, epoi_greedy, base_epoi_dev, epoi_base
  ## We want to plot the following:
  # The epoi_base should be plotted as a function of X, this should be a continuous line.
  # The line plotting epoi_base should be surrounded by a shaded area, which extends in both directions.
  # The epoi_cut and epoi_greedy should be plotted as a function of X, as continuous lines with different colors. 
  # The num_detached should be plotted as a function of X, as a bar plot behind all the other plots.

  # all in the same plot.
  # The plot should have a legend, and the axes should be labeled.
  # The plot should be saved in the results folder.
  plt.figure(figsize=(10,10))
  # Plot the epoi_base as a continuous line:
  plt.plot(data.index, data["epoi_base"], color="black", label="epoi_base")
  # Plot the epoi_cut and epoi_greedy as continuous lines:
  plt.plot(data.index, data["epoi_cut"], color="red", label="epoi_cut")
  plt.plot(data.index, data["epoi_greedy"], color="blue", label="epoi_greedy")
  # Plot a deviation plot as folows: plot a line above and below the epoi_base line, with a line width of 0.1, and a color of
  # a shade of gray that is ligher

  # Set the labels:
  plt.xlabel("Number of nodes")
  plt.ylabel("EPOI")
  # Set the legend:
  plt.legend()
  # Save the plot:
  plt.savefig("../results/comp_plot.png")
  plt.show()


## Take the dataframe and create a latex table from it. 
def create_latex_table(data:pd.DataFrame) -> None:
  ## Create a latex table with the following columns:
  # X, EPOI, EPOI after cut, EPOI after greedy. 
  # four rows: averages when X < 500, X between 500 and 100, X between 1000 and 1500, X > 1500
  # The averages should be rounded to two decimals.
  # The table should be saved in the results folder.
  # The table should be saved as a .tex file.
  # The table should have a caption and a label.
  # The table should be centered.
  table_data = []
  # Create the table data:
  # Create the header:
  row = []
  # Add the X:
  row.append("X")
  # Add the EPOI:
  row.append("EPOI")
  # Add the EPOI after cut:
  row.append("EPOI after cut")
  # Add the EPOI after greedy:
  row.append("EPOI after greedy")
  # Add the row to the table data:
  table_data.append(row)
  X_old = 0
  for X in [500, 1000, 1500, 2000]:
    # Create the row:
    row = []
    # Add the X:
    row.append(str(X_old) + "< X \leq " + str(X))
    # Add the EPOI:
    row.append(round(data.loc[X_old:X, "epoi_base"].mean(), 2))
    # Add the EPOI after cut:
    row.append(round(data.loc[X_old:X, "epoi_cut"].mean(), 2))
    # Add the EPOI after greedy:
    row.append(round(data.loc[X_old:X, "epoi_greedy"].mean(), 2))
    # Add the row to the table data:
    table_data.append(row)
    # Set the X_old to X:
    X_old = X
  # Create the latex table:
  latex_table = pd.DataFrame(table_data)
  # Save the latex table:
  latex_table.to_latex("../results/latex_table.tex", index=False, caption="Comparison of EPOI values for different algorithms", label="tab:epoi")

# Read all the graphs that are given in csv-format in the results folder
# and create a list of igraph objects

def read_graphs():
  graphs = {}
  #iterate over all .csv files in the ../results and the ../data folder
  directory = "../results/"
  for filename in os.listdir(directory):
    if filename.endswith(".csv"):
      with open(os.path.join(directory, filename), "r") as file:
        #read the file and create a graph object
        data = pd.read_csv(file)
        #create a graph object from the data
        graph = igraph.Graph.TupleList(data.itertuples(index=False), directed=False, weights=True)
        #add the graph to the dictionary, but use the name that has .csv removed, also some
        # have an ending ..csv, so remove that as well, but the name is of the form N_C_p, where
        # p is a probability and it has a dot in it, so split at the dot and take the first two
        # parts
        name =  filename.split(".")[0] + "." + filename.split(".")[1]
        graphs[name] = graph

  return graphs

# Create triplets of graphs, where each triplet contains the base graph, the cut graph and the greedy graph:
def triplets(graphs):
  #create a list of all the names of the graphs
  names = list(graphs.keys())
  #find the base names: these are the ones that have the form N_C_p, and the greedy and cut graphs
  # have names greedy_N_C_p and cut_N_C_p. Group these into triplets.
  base_names = set([])
  for name in names:
    if name.startswith("greedy"):
      base_names.add(name[7:])
    if name.startswith("cut"):
      base_names.add(name[4:])
    if name.startswith("original"):
      base_names.add(name[9:])
  triplet_dict = {}
  for name in base_names:
    try:
      triplet_dict[name] = (graphs["original_"+name], graphs["cut_" + name], graphs["greedy_" + name])
    except:
      print("Warning: The graph " + name + " does not have a cut or greedy graph, skipping")
      continue
  return triplet_dict

# create a layout for visualization of a given graph:

def create_layout(graph):
  #create a layout for the graph
  layout = graph.layout(layout = "fruchterman_reingold")
  #set the size of the layout
  layout.size = (800,800)
  #set the vertex size
  visual_style = {}
  visual_style["vertex_size"] = 20
  visual_style["layout"] = layout
  return visual_style



# create a visualization of graphs

def visualize():
  graphs = read_graphs()
  print("Graphs read")
  name_prefix = ["original", "cut", "greedy"]
  triplet_s = triplets(graphs)
  # create a gnuplot for each triplet
  for name in triplet_s:
    fig, ax = plt.subplots(1,3)
    g = triplet_s[name][1]
    visual_style = create_layout(g)
    Name_coord = {}
    for vertex_index,coord in enumerate(visual_style["layout"].coords):
      Name_coord[g.vs[vertex_index]["name"]] = coord
    print(Name_coord)
    for i,h in enumerate(triplet_s[name]):
      # adjust the visual style in such a way that the layout is the same for all three graphs, but
      # vertex names are laid out in the same places as in the original graph
      new_coords = []
      for vertex in h.vs:
        new_coords.append(Name_coord[vertex["name"]])
      modified_layout = igraph.Layout(new_coords)
      modified_visual_style = visual_style.copy()
      modified_visual_style["layout"] = modified_layout
      plot(h, target = ax[i], **modified_visual_style)
      ax[i].set_title(name + f"-{name_prefix[i]}")
      #print(g)
    plt.show()

if __name__ == "__main__":
  ic("starting visualization")
  data = read_num_results()
  ic(data)
  create_latex_table(data)
  comp_plot(data)

