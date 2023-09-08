#ifndef GRAPH_HPP
#define GRAPH_HPP

#include<vector>
#include<iostream>
#include<string>
#include<unordered_map>
#include<unordered_set>

// Forward declaration
namespace algo
{
  class GreedyOptimizer;
}


namespace graph
{
  // Graph class supports at most 32 bits
  class Graph
  {
    public:
      // Constructors: Default (empty graph)
      Graph(void);
      // Filename (read)
      Graph(std::string);
      // Random (modified stochastic block model) n vertex, r circles, p redundancy parameter
      Graph(int,int,double);

      // Copy constructor
      Graph(const Graph&);
      ~Graph(void);

      // Calculate a single spread pattern
      std::vector<bool> spread(int, std::unordered_map<size_t,double>&) const;
      //  Calculate activation probabilities for single circle
      std::vector<double> activation_probs(int,int=1000) const;

      //  Calculate activation heuristic probabilities for single circle
      std::vector<double> numeric_spread(int,std::unordered_map<size_t,double>&) const;
      // Calculate EPOI, (number of simulations, use initial_probabilities distribtion)
      double EPOI(int=1000,bool=false) const; 
      Graph detach(int, int) const;
      double get_prob(int,int) const;
      // DMP Epoi
      double DMP_EPOI(bool=false) const;

      // Calculate EPOI, (number of simulations, use initial_probabilities distribtion)
      double numeric_EPOI(int=1000,bool=false) const; 

      // Print to file name 
      void print(std::string);
      void printcsv(std::string,bool=false);

      // Friends and family
      friend class algo::GreedyOptimizer;

    private:
      int _num_nodes;
      size_t Key(const int&, const int&) const; 
      std::vector<std::vector<int>> _circles;
      std::vector<std::vector<int>> _circle_of_node;
      std::vector<std::unordered_set<int>> _adjacent; 
      std::vector<double> _initial_probabilities; 
      std::unordered_map<size_t, double> _probs; 
      // make circles into connections parameter is max number
      void instantiate();
      // Randomize the connections, i.e, give each existing connection a coefficient from beta distribution.
      // Default values are mu=0.2, sample size=100 
      void randomize(double=20,double=80);
      std::vector<double> _activations; 
      void clean_up();
  };
}

#endif