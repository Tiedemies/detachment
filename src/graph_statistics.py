# Input-output class for the bridge-block network
import igraph as ig
import matplotlib.pyplot as plt


class IO:
  def __init__(self, filename):
    self.filename = filename
    # Read the text file and create the briges and blocks
    # a single line of the file describes a block: b;[a1, ... , an]
    # where b is the block number and a1, ... , an are the vertices in the block
    inputfile = open(filename, 'r')
    self.nodes = set([])
    self.blocks = {}
    for line in inputfile.readlines():
      line = line.strip()
      line = line.split(';')
      block = int(line[0])
      vertices = line[1].strip('[]').split(',')
      vertices = [int(v) for v in vertices]
      self.blocks[block] = vertices
      self.nodes.update(vertices)
    inputfile.close()
    # Create the graph. Each bridge is a vertex that belongs to more than one block
    # The graph is bipartite: bridges are connected to blocks, blocks are connected to bridges
    # First create a dictionary that lists which blocks each vertex belongs to:
    self.vertex_blocks = {}
    for block in self.blocks:
      for vertex in self.blocks[block]:
        if vertex not in self.vertex_blocks:
          self.vertex_blocks[vertex] = []
        self.vertex_blocks[vertex].append(block)
    # now recognize the bridges:
    self.bridges = []
    for vertex in self.vertex_blocks:
      if len(self.vertex_blocks[vertex]) > 1:
        self.bridges.append(vertex)
    # Then create the edges:
    self.edges = []
    for vertex in self.bridges:
      for block in self.vertex_blocks[vertex]: 
        self.edges.append((vertex, block))
    # Now the problem is that block and bridge numbers overlap
    # so we need to renumber the bridges and blocks. Blocks will be numbered
    # from 0 to len(self.blocks)-1 and bridges from len(self.blocks) to len(self.blocks)+len(self.bridges)-1
    self.bridge_renumbering = {}
    self.block_renumbering = {}
    i = 0
    for b in self.blocks:
      self.block_renumbering[b] = i
      i+=1
    for b in self.bridges:
      self.bridge_renumbering[b] = i
      i+=1
    # Now we can renumber the edges:
    try:
      self.edges = [(self.bridge_renumbering[edge[0]], self.block_renumbering[edge[1]]) for edge in self.edges]
    except Exception as e:
      print("Error in renumbering the edges")
      print([edge for edge in self.edges if edge[0] not in self.bridge_renumbering or edge[1] not in self.block_renumbering])
      raise e

    # Create the graph:
    self.graph = ig.Graph()
    self.graph.add_vertices(len(self.blocks) + len(self.bridges))
    self.graph.add_edges(self.edges)
    # Set the names of the vertices:
    self.graph.vs['name'] = ["block_" + str(i) for i in range(len(self.blocks))] + ["bridge_" + str(i) for i in range(len(self.bridges))]
    # Set the names of the edges:
    self.graph.es['name'] = ["bridge_" + str(edge[0]) + "_block_" + str(edge[1]) for edge in self.edges]
    ## Calculate statistics of the graph:
    # Degree of each vertex:
    self.vertex_degrees = self.graph.degree()
    print("Vertex degrees: ", self.vertex_degrees)
    # Number of blocks:
    self.nblocks = len(self.blocks)
    self.nbridges = len(self.bridges)
    print("Number of nodes: ", len(self.nodes))
    print("Number of blocks: ", self.nblocks)
    print("Number of bridges: ", self.nbridges)
    # Number of bridges:
  def plot_statistics(self):
    # create a histogram of block sizes and bridge degrees and plot them side by side:
    plt.subplot(121)
    plt.hist(self.vertex_degrees[:self.nblocks], bins=range(max(self.vertex_degrees[:self.nblocks])+2))
    plt.xlabel('Block size')
    plt.ylabel('Frequency')
    plt.subplot(122)
    plt.hist(self.vertex_degrees[self.nblocks:], bins=range(max(self.vertex_degrees[self.nblocks:])+2))
    plt.xlabel('Bridge degree')
    plt.ylabel('Frequency')
    plt.show()  
    # print some statistics of the graph to the screen:
    print("Number of blocks: ", self.nblocks)
    print("Number of bridges: ", len(self.bridges))
    # Calculate the tree width of the block-bridge graph:
    print("Tree width: ", self.graph.tree_width())

if __name__ == "__main__":
  # Test the class:
  io = IO("../data/company_dict_insiders.txt")
  io.plot_statistics()

    

