# Input-output class for the bridge-block network
import igraph as ig
import matplotlib.pyplot as plt


class IO:
  def __init__(self, filenames):
    if type(filenames) == str:
      filenames = [filenames]
    self.graphs = []
    self.nblocks = []
    self.blockdegrees = []
    self.nnodes = []
    for filename in filenames: 
      # Read the text file and create the briges and blocks
      # a single line of the file describes a block: b;[a1, ... , an]
      # where b is the block number and a1, ... , an are the vertices in the block
      inputfile = open(filename, 'r')
      nodes = set([])
      blocks = {}
      for line in inputfile.readlines():
        line = line.strip()
        if line == "":
          continue
        if ';' not in line:
          line = line.split('[')
          line[1] = "[" + line[1]
        else:
          line = line.split(';')
        block = int(line[0])
        vertices = line[1].strip('[]').split(',')
        vertices = [int(v) for v in vertices]
        blocks[block] = vertices
        nodes.update(vertices)
      inputfile.close()
      # Create the graph. Each bridge is a vertex that belongs to more than one block
      # The graph is bipartite: bridges are connected to blocks, blocks are connected to bridges
      # First create a dictionary that lists which blocks each vertex belongs to:
      vertex_blocks = {}
      for block in blocks:
        for vertex in blocks[block]:
          if vertex not in vertex_blocks:
            vertex_blocks[vertex] = []
          vertex_blocks[vertex].append(block)
      # now recognize the bridges:
      bridges = []
      for vertex in vertex_blocks:
        if len(vertex_blocks[vertex]) > 1:
          bridges.append(vertex)
      # Then create the edges:
      edges = []
      active_blocks = set([])
      for vertex in bridges:
        for block in vertex_blocks[vertex]:
          edges.append((vertex, block))
          active_blocks.add(block)
      
      # Now the problem is that block and bridge numbers overlap
      # so we need to renumber the bridges and blocks. Blocks will be numbered
      # from 0 to len(self.blocks)-1 and bridges from len(self.blocks) to len(self.blocks)+len(self.bridges)-1
      bridge_renumbering = {}
      block_renumbering = {}
      i = 0
      for b in blocks:
        if b not in active_blocks:
          continue
        block_renumbering[b] = i
        i+=1
      for b in bridges:
        bridge_renumbering[b] = i
        i+=1
      # Now we can renumber the edges:
      try:
        edges = [(bridge_renumbering[edge[0]], block_renumbering[edge[1]]) for edge in edges]
      except Exception as e:
        print("Error in renumbering the edges")
        print([edge for edge in edges if edge[0] not in bridge_renumbering or edge[1] not in block_renumbering])
        raise e

      # Create the graph:
      graph = ig.Graph()
      graph.add_vertices(len(active_blocks) + len(bridges))
      graph.add_edges(edges)
      # Set the names of the vertices:
      graph.vs['name'] = ["block_" + str(i) for i in range(len(active_blocks))] + ["bridge_" + str(i) for i in range(len(bridges))]
      # Set the names of the edges:
      graph.es['name'] = ["bridge_" + str(edge[0]) + "_block_" + str(edge[1]) for edge in edges]
      self.graphs.append(graph)
      self.nblocks.append({block_renumbering[block]:len(blocks[block]) for block in active_blocks})
      self.blockdegrees.append({block_renumbering[block]:0 for block in active_blocks})
      for edge in edges:
        self.blockdegrees[-1][edge[1]] += 1
      self.nnodes.append(len(nodes))
      statstring = "Graph " + filename + " has " + str(len(active_blocks)) + " blocks," + str(len(bridges)) + " bridges, and " + str(len(nodes)) + " nodes" 
      print (statstring)  

  ## Calculate statistics of the graph:
  def plot_statistics(self):
    nplot = 140
    # create a subplot for each graph:
    for i,graph in enumerate(self.graphs):
      #create a histogram of block degrees and bridge degrees and plot them side by side:
      vertex_degrees = graph.degree()
      plt.subplot(nplot+1)
      nblocks = len(self.nblocks[i])
      plt.hist(vertex_degrees[:nblocks], alpha = 0.5, bins=range(max(vertex_degrees[:nblocks])+2))
      plt.xlabel('Block degree')
      plt.ylabel('Frequency')
      plt.subplot(nplot+2)
      plt.hist(vertex_degrees[nblocks:], alpha = 0.5, bins=range(max(vertex_degrees[nblocks:])+2))
      plt.xlabel('Bridge degree')
      plt.ylabel('Frequency')
      # create a histogram of block sizes:
      plt.subplot(nplot+3)
      plt.hist(self.nblocks[i].values(), alpha = 0.5, bins=range(max(self.nblocks[i].values())+2))
      plt.xlabel('Block size')
      plt.ylabel('Frequency')
      print("Graph " + str(i) + " has " + str(len(graph.es)) + " edges and " + str(len(graph.vs)) + " vertices")
      # Plot the block sizes, sorted by size:
      plt.subplot(nplot+4)
      plt.plot(sorted(self.nblocks[i].values(),  reverse=True), alpha = 0.5)
      plt.xlabel('Block index')
      plt.ylabel('Block size')

    plt.show()  

if __name__ == "__main__":
  # Test the class:
  io = IO(["../data/empirical.txt","../data/randomgraph_1000.txt"])
  io.plot_statistics()

    

