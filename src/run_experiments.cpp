#include "Graph/graph.hpp"
#include<iostream>
#include <chrono>
#include "Algorithms/algorithms.hpp"
#include "Utils/utils.hpp"


void run_test()
{
  graph::Graph test1("../data/company_dict_insiders.txt");
  //graph::Graph test2(105, 10, 0.5);
  // test1.printcsv("big.csv");
  algo::GreedyOptimizer optim(test1);
  int tests = 5;
  for (int k = 1; k < 8; ++k)
  {
    std::cout << k << ": \n";
    for (int i = 0; i < tests; ++i)
    {
      DEBUG("Start optimize for " << k);
      optim.optimize(k);
      DEBUG("Optimization complete \n");
      std::cout << optim._base_epoi << ", " << optim._result_epoi << "\n";
      for (auto cpair: optim._detached)
      {
        std::cout << "(" << cpair.first << "," << cpair.second << ")";
      }
      std::cout << "\n **** \n";
    }
    std::cout << "************** \n";
  }
}

int main()
{
  run_test(); 
  return 0;
}