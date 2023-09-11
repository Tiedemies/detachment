#include "algorithms.hpp"
#include<numeric>
#include<math.h>
#include<iostream>
#include "../Graph/graph.hpp"
#include<boost/math/distributions/normal.hpp>

namespace algo
{

  GreedyOptimizer::GreedyOptimizer(graph::Graph g): _parent(g)
  {
    //void
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
    std::cerr << "average: " << mu << "\n";
    return confidence_intervals(input,p,mu);
  }

  std::pair<double,double> confidence_intervals(const std::vector<double>& input, const double& p, const double& mu)
  {
    size_t n = input.size();
    double s2 = variance(input,mu);
    std::cerr << "mu: " << mu << "\n";
    boost::math::normal_distribution nd;
    double z = boost::math::quantile(nd, 1.0 - ( (1-p)/2.0)); 
    double delta = std::sqrt(s2 / static_cast<double>(n))  * z;
    std::cerr << "delta: " << delta << "\n";
    return std::make_pair(mu-delta,mu+delta);


  }

}