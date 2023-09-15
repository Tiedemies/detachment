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

  GreedyOptimizer::GreedyOptimizer(graph::Graph& g): _parent(g)
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
          graph::Graph temp = _result.detach(j,c);
          double comp_epoi = temp.EPOI_num(base_sim,false);
          if (comp_epoi < cur_min)
          {
            temporary_min = temp;
            cur_min = comp_epoi;
            rpair.first = j; rpair.second = c; 
          }
        }
      }
      _result = temporary_min;
      _result_epoi = cur_min; 
      _detached.insert(rpair); 
    }
  }

}