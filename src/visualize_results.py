import igraph
import os
import pandas as pd
import matplotlib.pyplot as plt
import itertools
from igraph import plot

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

if __name__ == '__main__':
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

    