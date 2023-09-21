#include "algorithms.hpp"
#include<numeric>
#include<math.h>
#include<set>
#include<iostream>
#include "../Graph/graph.hpp"
#include "../Utils/utils.hpp"
#include<boost/math/distributions/normal.hpp>

namespace algo
{

  static double previous_mu = 0.0;
  static double previous_sigma = 0.0;

  Optimizer::Optimizer(const graph::Graph& g): _parent(g)
  {
    // void
  }

  GreedyOptimizer::GreedyOptimizer(const graph::Graph& g) : Optimizer(g)
  {
    // void
  }

  Cutter::Cutter(const graph::Graph& g) : Optimizer(g)
  {
    // void
  }


  inline double mean(const std::vector<double>& input)
  {
    if (input.size() == 0)
    {
      return 0.0;
    }
    return std::accumulate(input.begin(),input.end(),0.0)/static_cast<double>(input.size()); 
  }


  inline double variance(const std::vector<double>& input)
  {
    double mu = mean(input);
    return variance(input,mu);
  }

  inline double variance(const std::vector<double>& input, const double& mu)
  {
    size_t n = input.size();
    if (n <= 1)
    {
      return 0;
    }
    auto var_func = [&mu, &n](double accumulator, const double& val)
    { 
      return accumulator + (val - mu)*(val - mu)/(n-1); 
    };
    return std::accumulate(input.begin(),input.end(),0.0,var_func); 
  }

  std::pair<double,double> confidence_intervals(const std::vector<double>& input, const double& p)
  {
    const double mu = mean(input);
    DEBUG("average: " << mu);
    return confidence_intervals(input,p,mu);
  }

  std::pair<double,double> confidence_intervals(const std::vector<double>& input, const double& p, const double& mu)
  {
    size_t n = input.size();
    double s2 = variance(input,mu);
    double sigma = std::sqrt(s2);
    DEBUG("mu: " << mu << " sigma: " << sigma);
    DEBUG("mu-difference: " << mu - previous_mu << ", sigma-difference: " << sigma - previous_sigma);
    previous_mu = mu; 
    previous_sigma =  sigma;
    boost::math::normal_distribution nd;
    double z = boost::math::quantile(nd, 1.0 - ( (1-p)/2.0)); 
    double delta = std::sqrt(s2 / static_cast<double>(n))  * z;
    DEBUG("delta: " << delta << "\n");
    return std::make_pair(mu-delta,mu+delta);
  }

  /// @brief Optimize by finding the best k-detachment with a greedy algorithm
  /// @param k 
  void
  GreedyOptimizer::optimize(int k)
  {
    const int base_sim = 10000;
    // Initialize a copy.
    const int n = _parent._num_nodes;
    DEBUG("working with " << n << " nodes");
    std::vector<double> result(2*base_sim,0.0);
    DEBUG("establishing baseline");
    _base_epoi = _parent.EPOI_num(2*base_sim, false, result);
    DEBUG("calculating intervals");
    const auto ref_interval = confidence_intervals(result,0.99);
    double cur_min = _base_epoi; 
    DEBUG("base interval established " << ref_interval.first << " -- " << _base_epoi << " -- " << ref_interval.second);
    _result = _parent;
    graph::Graph temporary_min = _parent;
    _detached.clear(); 
    for (int i = 0; i < k; ++i)
    {
      std::pair<int,int> rpair; 
      for (int j = 0; j < n; ++ j)
      {
        // Never try to remove one that is in only one circle.
        if (_result._circle_of_node[j].size() <= 1)
        {
          continue;
        }
        for (int c: _result._circle_of_node[j])
        {
          _result.detach(j,c);
          double comp_epoi = _result.EPOI_num(base_sim,false);
          if (comp_epoi < cur_min)
          {
            cur_min = comp_epoi;
            rpair.first = j; rpair.second = c; 
          }
          _result.undo(); 
        }
      }
      _result.detach(rpair.first,rpair.second);
      _result_epoi = cur_min; 
      _detached.insert(rpair); 
      _result = _result.clean_copy(); 
    }
  }

  /// @brief Optimize to make an n-cut
  /// @param n 
  void 
  Cutter::optimize(int n)
  {
    pre_process();
    
  }

  void 
  Cutter::pre_process()
  {
    int bridge_to_block = 0;
    int block_to_bridge = 0;
    int num_bridges = 0;  
    for (int c = 0; c < _parent._circles.size(); ++ c)
    {
      for (int u: _parent._circles.at(c))
      {
        // If the node is in more than one circle, then it is a bridge
        if (_parent._circle_of_node.at(u).size() > 1)
        {
          auto br_bl_n =_bridge_block_edges[u];
          br_bl_n.push_back(c);
          auto bl_br_n = _block_bridge_edges[c];
          bl_br_n.push_back(u);
        }
      }
    }
  }

}