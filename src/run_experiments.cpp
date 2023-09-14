#include "Graph/graph.hpp"
#include<iostream>
#include <chrono>
#include "Algorithms/algorithms.hpp"


int main()
{
  std::cerr << "foo";
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
      optim.optimize(k);
      std::cout << optim._base_epoi << ", " << optim._result_epoi << "\n";
      for (auto cpair: optim._detached)
      {
        std::cout << "(" << cpair.first << "," << cpair.second << ")";
      }
      std::cout << "\n **** \n";
    }
    std::cout << "************** \n";
  }
  return 0;
}