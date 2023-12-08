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
  class Optimizer;
  class GreedyOptimizer;
  class Cutter; 
}


namespace graph
{
  static std::vector<double> DEFAULT_VECTOR;

  // Graph class for detacment problem on the stochastic block graphs. 
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

      /// @brief a single Monte carlo spread to find a single cascade number of activatiaons
      /// @param  int: index of the source circle 
      /// @param  vector<double>& :  A stored vector reference to implement antithetic distribution
      /// @return int: number of activations
      int spread_num(int, std::vector<double>&) const;

      /// @brief Monte Carlo approximation of a single circle epoi upon infection
      /// @param  int: The circle index
      /// @param  int: The number of simulations
      /// @return double: estimated EPOI. 
      double circle_epoi(int, int = 10000, std::vector<double>& = DEFAULT_VECTOR) const; 

      
      /// @brief Detach a vertex - circle pair
      /// @param  int: vertex
      /// @param  int: circle
      void detach(int, int);

      /// @brief Undo the last detachment
      void undo(); 

      /// @brief Resolve the detachment stack and create a copy of the implemented detachments
      /// @return A fresh graph with the detachments implemented and a clear cache.
      Graph clean_copy() const; 
      
      /// @brief Get the probability of transmission u -> v
      /// @param  int u: vertex index
      /// @param  int v: vertex index
      /// @return double p: probability of transmission
      double get_prob(int,int) const;
      
      /// @brief DMP EPOI calculation
      /// @param  bool: whether to use a predetermined distribution (false: uniform)
      /// @param  vector: A vector to store the result vector of the simulations, if none is given, it is omitted
      /// @return double: the EPOI estimate
      double DMP_EPOI(bool=false) const;
      std::vector<double> RunDMP(const std::vector<int>&) const;

      /// @brief Numeric (Monte Carlo) calculation of EPOI
      /// @param  int: Number of simulations per circle
      /// @param  bool: Whether to use predetermined distribution (false = uniform)
      /// @param  vector: A vector to store the result vector of the simulations, if none is given, it is omitted
      /// @return double: The EPOI estimate 
      double EPOI_num(int=1000,bool=false, std::vector<double>& =DEFAULT_VECTOR) const; 

      // Print to file name
      void print(std::string);
      void printcsv(std::string,bool=false) const;

      // Friends and family
      friend class algo::Optimizer;
      friend class algo::GreedyOptimizer;
      friend class algo::Cutter;

      /// @brief Check if two vertices are adjacent
      /// @param  int u: integer vertex number
      /// @param  int v: integer vertex number
      /// @return boolean that tells you if (u,v) is an edge 
      bool is_adjacent(int,int) const;

      /// @brief Get an iterable of vertices adjacent to u
      /// @param  int u: integer vertex
      /// @return const std::vector<int> reference to adjacent vertices
      std::vector<int> adjacent_to(int) const; 
      int Get_num_nodes() const;

    private:
      
      size_t Key(const int&, const int&) const; 
      // inverse Key
      std::pair<int,int> InvKey(const size_t&) const;

      int EdgeIndex(const int&, const int&) const;
      
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