#ifndef ALGO_HH
#define ALGO_HH

#include<vector>
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
}

#endif