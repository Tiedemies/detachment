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
      const graph::Graph& get_result() const;
      double _result_epoi; 
      double _base_epoi;
      std::set<std::pair<int,int>> _detached; 
      void set_iterations(int);
      int _sims; 
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
      void initialize_circle_epois();
      void initialize_division(int k = 1);


      double _max_weight; 
    
      std::unordered_map<int,std::vector<int>> _bridge_block_edges;
      std::unordered_map<int,std::vector<int>> _block_bridge_edges; 
      std::unordered_map<size_t,double> _bridge_block_probs;
      std::unordered_map<size_t,double> _block_bridge_probs;

      std::vector<double> _circle_epois;

      // SOurces and sinks
      std::vector<int> _S_blocks;
      std::vector<int> _T_blocks;

      // @brief Not using atm 
      // std::unordered_map<std::pair<int,int>, double> _bridge_block_probs;
      // std::unordered_map<std::pair<int,int>, double> _block_bridge_probs;

      /// @brief Minimum cut Edmonds Karp
      /// @return detachment vector
      std::vector<std::pair<int,int>> min_cut() const;
      std::vector<std::pair<int,int>> extract_cut(std::unordered_map<size_t, int>&, std::unordered_map<size_t,int>&) const;
      void initialize_residual(std::unordered_map<size_t, int>&, std::unordered_map<size_t,int>&) const;
      void update_residual(std::unordered_map<size_t, int>&, std::unordered_map<size_t,int>&,const std::vector<int>&) const;

      std::vector<int> find_aug_path(std::unordered_map<size_t, int>&, std::unordered_map<size_t,int>&) const;
      size_t Key(int,int) const;


      // Debug functions
      bool verify_connectedness() const;
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