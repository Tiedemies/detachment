#ifndef ALGO_HH
#define ALGO_HH

#include<vector>
#include<numeric>
#include<unordered_map>
#include "../Graph/graph.hpp"


namespace algo
{
  class GreedyOptimizer
  {
    public:
      GreedyOptimizer() = delete;
      GreedyOptimizer(graph::Graph);
      void optimize(int=1);
      const graph::Graph& get_result();

      void set_iterations();
    private:
      const graph::Graph& _parent; 
      graph::Graph _result;
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