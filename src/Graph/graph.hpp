#ifndef GRAPH_HPP
#define GRAPH_HPP

#include<vector>
#include<iostream>
#include<string>
#include<unordered_map>
#include<unordered_set>

namespace graph
{
  // Graph class supports at most 32 bits
  class Graph
  {
    public:
      // Constructors
      Graph(void);
      Graph(std::string);
       Graph(const Graph&);
      ~Graph(void);

      // Calculate a single spread pattern
      std::vector<bool> spread(int) const;
      // Calculate EPOI
      double EPOI(int,bool) const; 
      Graph detach(int, int) const;
      double get_prob(int,int) const; 

    private:
      int _num_nodes;
      size_t Key(const int&, const int&) const; 
      std::vector<std::vector<int>> _circles;
      std::vector<std::unordered_set<int>> _adjacent; 
      std::vector<double> _initial_probabilities; 
      std::unordered_map<size_t, double> _probs; 
      // make circles into connections parameter is max number
      void instantiate();
      // Randomize the connections
      void randomize(double);
  };
}

#endif