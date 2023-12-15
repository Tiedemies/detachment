# Random generator for graphs that exhibit the bridge-block structure

import random
import sys
from icecream import ic

class RandomGraphGenerator:
  def __init__(self, n: int, m: int, k: int, e: int, p: float, q: float) -> None:
    # n: number of nodes
    self.n = n
    self.nodes = 0
    # m: number of blocks
    self.m = m
    self.blocknumber = 0
    # k: number of bridges
    self.k = k
    self.brigenumber = 0
    # e: number of edges
    self.e = e
    self.edges = []

    # p: prerential attachment parameter for the bridges
    self.p = p
    # q: the preferential attachment parameter for the blocks
    self.q = q
    self.weights = []
    self.bridge_weights = []
    self.blocks =  []
    self.bridges = []
    self.block_of_node = {}
  """ Update the weights of the blocks """
  def update_weights(self, q: float = None):
    if q is None:
      q = self.q
    if not self.weights:
      for i in range(self.m):
        self.weights.append(1/self.m)
      return
    for i in range(self.m):
      self.weights[i] = (1-q)*(1/self.m) + q*len(self.blocks[i])/self.nodes

  """ Update the weights of the bridges """
  def update_bridge_weights(self, q: float = None):
    if q is None:
      q = self.p
    if not self.bridge_weights:
      for i in range(self.k):
        self.bridge_weights.append(1/self.k)
      return
    for i in range(self.k):
      self.bridge_weights[i] = (1-q)*(1/self.k) + q*(len(self.block_of_node[self.bridges[i]])-1)/len(self.edges)
  
  """ Select a block with probability proportional to its weight """
  ## Exclude any block in the list given
  def select_block(self, exclude: list = []) -> int:
    r = random.random()
    blocks_to_select = list(set(range(self.m)) - set(exclude))
    blocks_to_select.sort()
    weight_sum = sum([self.weights[i] for i in blocks_to_select])
    s = 0
    for b in blocks_to_select:
      s += self.weights[b]/weight_sum
      if r < s:
        return b
    return blocks_to_select[-1]
    
  
  def select_bridge(self):
    if len(self.bridge_weights) < self.k:
      self.bridge_weights = [1/self.k for i in range(self.k)]
    r = random.random()
    s = 0
    for i in range(self.k):
      s += self.bridge_weights[i]
      if r < s:
        return self.bridges[i]
    return self.bridges[-1]
  
  """ Generate the graph """
  def generate(self)-> None:
    # The blocks are represented as a list of lists. The only guarantee is that each block has at least one node.
    for i in range(self.m):
      self.blocks.append([i])
      self.block_of_node[i] = [i]
      self.nodes += 1
    self.update_weights()
    for j in range(self.m, self.n):
      # select a block
      block = self.select_block()
      self.blocks[block].append(j)
      self.block_of_node[j] = [block]
      self.nodes += 1
      self.update_weights(q = 0.02)
    # First we make sure the graph is connected, by generating bridges between the blocks.   
    # To do this we create a "tree" of blocks. Select a random node:
    node = random.randint(0, self.n-1)
    block = self.block_of_node[node][0]
    connected_nodes = set(self.blocks[block])
    unconnected_nodes = set(range(self.n)) - connected_nodes
    while unconnected_nodes:
      # select a random node from the unconnected nodes
      node =  random.choice(list(unconnected_nodes))
      block = self.block_of_node[node][0]
      # select a random connected node:
      connected_node = random.choice(list(connected_nodes))
      # pick the block of the connected node
      connected_block = self.block_of_node[connected_node][0]
      # add node to the connected block
      self.blocks[connected_block].append(node)
      # update the block of the node
      self.block_of_node[node].append(connected_block)
      # update the sets
      connected_nodes.update(self.blocks[block])
      unconnected_nodes = unconnected_nodes - connected_nodes
      self.edges.append((node, connected_block))
    # Detecte bridges:
    self.bridges = [node for node in range(self.n) if len(self.block_of_node[node]) > 1]
    # Now we add the remaing edges, by selecting random nodes and assigning them to random extra blocks.
    # Random nodes are chosen with a mix of uniform and preferential attachment, with parameter p. 
    for i in range(self.k-len(self.bridges)):
      # select a brigde:
      bridges_to_select = list(set(range(self.k)) - set(self.bridges))
      bridge = random.choice(bridges_to_select)
      # select a block in which the bridge is not present:
      block = self.select_block(exclude=self.block_of_node[bridge])
      # add the bridge to the block:
      # add the block to the bridge:
      self.blocks[block].append(bridge)
      self.block_of_node[bridge].append(block)
      self.bridges.append(bridge)
      self.edges.append((bridge, block))
    assert (len(self.bridges) == self.k)
    # Add the remaining edges by randomly picking bridges:
    ic(len(self.edges))
    ic(self.e)
    ic(len(self.blocks))
    ic(len(self.bridges))
    n = self.e - 2*len(self.edges)
    for i in range(n):
      bridge = self.select_bridge()
      block = self.select_block(exclude=self.block_of_node[bridge])
      self.blocks[block].append(bridge)
      self.block_of_node[bridge].append(block)
      self.edges.append((bridge, block))
      self.update_bridge_weights()
    
  def print(self, filename:str) -> None:
    # print all blocks into the file
    with open(filename, "w") as f:
      for i in range(self.m):
        f.write(str(i) + ";")
        f.write(str(self.blocks[i]) + "\n")
    

if __name__ == "__main__":
  base_n = 967
  base_m = 106
  base_k = 140
  base_e = 330
  base_p = 0.3
  base_q = 0.3
  base_mn_ratio = base_m/base_n
  base_kn_ratio = base_k/base_n
  base_em_ratio = base_e/base_m
  for n in range (100,2000,20):
    m = int(base_mn_ratio*n)
    k = int(base_kn_ratio*n)
    e = int(base_em_ratio*m)
    p = base_p
    q = base_q
    g = RandomGraphGenerator(n=n, m=m, k=k, e=e, p=p, q=q)
    g.generate()
    g.print("../data/randomgraph_" + str(n) + ".txt") 

