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
  static std::vector<double> DEFAULT_VECTOR;
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
      // std::vector<bool> spread(int, std::unordered_map<size_t,double>&) const;

      // Calculate a single spread pattern
      int spread_num(int, std::vector<double>&) const;
      //  Calculate activation probabilities for single circle
      // std::vector<double> activation_probs(int,int=1000) const;

      //  Calculate activation heuristic probabilities for single circle
      // std::vector<double> numeric_spread(int,std::unordered_map<size_t,double>&) const;
      // Calculate EPOI, (number of simulations, use initial_probabilities distribtion)
      // double EPOI(int=1000,bool=false,std::vector<double>& =DEFAULT_VECTOR) const; 
      
      // Detach and undo. 
      void detach(int, int);
      void undo(); 

      // Resolve the detach stack and create a copy that is clean and delete. 
      Graph clean_copy() const; 
      
      double get_prob(int,int) const;
      // DMP Epoi
      double DMP_EPOI(bool=false) const;

      // Calculate EPOI, (number of simulations, use initial_probabilities distribtion, verbose [i.e., count deviation])
      double EPOI_num(int=1000,bool=false, std::vector<double>& =DEFAULT_VECTOR) const; 

      // Print to file name 
      void print(std::string);
      void printcsv(std::string,bool=false);

      // Friends and family
      friend class algo::GreedyOptimizer;

    private:
      
      size_t Key(const int&, const int&) const; 
      
      // Copyable attributes:
      int _num_nodes;
      int _num_edges;
      std::vector<std::vector<int>> _circles;
      std::vector<std::vector<int>> _circle_of_node;
      std::vector<double> _initial_probabilities; 
      std::vector<double> _prob_vector; 
      std::vector<int> _adjacency_vector;
      std::vector<int> _node_places;
      std::vector<std::pair<int,int>> _detachment_stack;
      // Total 6.


      // make circles into connections parameter is max number
      void instantiate();
      // Randomize the connections, i.e, give each existing connection a coefficient from beta distribution.
      // Default values are mu=0.2, sample size=100 
      void randomize(double=20,double=80);
      std::vector<double> _activations; 
      void clean_up();
      void validate();
  };
}

#endif