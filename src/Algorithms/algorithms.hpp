#ifndef ALGO_HH
#define ALGO_HH

#include<vector>
#include<numeric>
#include<set>
#include<unordered_map>
#include<unordered_set>
#include "../Graph/graph.hpp"


namespace algo
{
  
  class Optimizer
  {
    public:
      Optimizer() = delete;
      Optimizer(const graph::Graph&);
      void optimize(int=1) = delete;
      const graph::Graph& get_result();
      double _result_epoi; 
      double _base_epoi;
      std::set<std::pair<int,int>> _detached; 
      void set_iterations();
    protected:
      void pre_process();
      const graph::Graph& _parent; 
      graph::Graph _result; 
       
  };
  
  class GreedyOptimizer : public Optimizer
  {
    public:
      GreedyOptimizer() = delete;
      GreedyOptimizer(const graph::Graph&);
      void optimize(int=1);
  };
  
  class Cutter : public Optimizer
  {
    public:
      Cutter() = delete;
      Cutter(const graph::Graph&);
      void optimize(int=1);
    protected:
      /// @brief Create an alternative internal representation using block structure
      void pre_process();
    
      std::unordered_map<int,std::vector<int>> _bridge_block_edges;
      std::unordered_map<int,std::vector<int>> _block_bridge_edges; 

      /// @brief Not using atm 
      // std::unordered_map<std::pair<int,int>, double> _bridge_block_probs;
      // std::unordered_map<std::pair<int,int>, double> _block_bridge_probs;
  };


  /// Estimate mean
  inline double mean(const std::vector<double>& input);
  // Estimate variance
  inline double variance(const std::vector<double>& input);
  inline double variance(const std::vector<double>&, const double&);
  // 
  std::pair<double,double> confidence_intervals(const std::vector<double>&, const double&);
  std::pair<double,double> confidence_intervals(const std::vector<double>&, const double&, const double&);

}

#endif