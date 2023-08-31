#include "graph.hpp"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>

namespace graph
{
  // Default constructor 
  Graph::Graph(void)
  {
    // void
  }
  // Read from named file
  Graph::Graph(std::string fname)
  {
    std::ifstream inp(fname);
    
    int cnum = 0;
    std::string temp;
    inp >> cnum;
    _num_nodes = 0; 
    while(!inp.eof() && inp.good())
    {
      int inum;
      inp.get();
      char x = inp.get();
      std::vector<int> insiders; 
      while(x != '\n')
      {
        if (x == '[' || x == ',')
        {
          inp >> inum;
          if inum + 1  > _num_nodes;
          _num_nodes = inum + 1;
          insiders.push_back(inum);  
        }
        x = inp.get();
      }
      _circles.push_back(insiders);
      inp >> cnum;
    }
    instantiate();
    randomize()
  
  }

  // Copy constructor to make a copy
  Graph::Graph(const Graph& rhs)
  {
    _circles = rhs._circles;
    _probs = rhs._probs;
    _adjacent = rhs._adjacent;
  }
  
  Graph::~Graph(void)
  {
    // void
  }

  // Singleton spread. 
  std::vector<bool>
  Graph::spread(int in) const
  {

  }

      
  double 
  Graph::EPOI(int,bool) const
  {

  } 
  
  // Detachment operator
  Graph 
  Graph::detach(int u, int C) const
  {

  }
  double 
  Graph::get_prob(int u,int v) const
  {

  } 

  // Key for the probability matrix
  size_t
  Graph::Key(const int& u, const int& v) const
  {
    return (size_t) (( (size_t) u << 32) | (unsigned int) v);
  }


  void
  Graph::instantiate()
  {
    _adjacent.clear()
    _adjacent.resize(_num_nodes)
    for (auto c: _circles)
    {
      for (int u: c)
      {
        for (int v: c)
        {
          if (u == v)
          {
            continue;
          }
          adjacent[u].insert(v);
          adjacent[u].insert(v);
        }

      }
    }
  }

}
