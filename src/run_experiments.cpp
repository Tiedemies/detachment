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

void run_cut_test()
{
  graph::Graph test1("../data/company_dict_insiders.txt");
  //graph::Graph test2(105, 10, 0.5);
  // test1.printcsv("big.csv");
  algo::Cutter optim(test1);
  DEBUG("Start max cut optimizer ");
  optim.optimize();
  DEBUG("Optimization complete \n");
  std::cout << "Base: " << optim._base_epoi << ", Optimized: " << optim._result_epoi << "\n";
  for (auto cpair: optim._detached)
  {
      std::cout << "(" << cpair.first << "," << cpair.second << ")";
  }
  optim.get_result().printcsv("cut_result.csv");
}

void speed_test()
{
  for (int n =10; n < 10001;n=n*10)
  {
    std::vector<double> result(n,0.0);
    graph::Graph test1("../data/company_dict_insiders.txt");

    auto start_time = std::chrono::system_clock::now();
    double ep = test1.EPOI_num(n,false,result);
    auto end_time = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_time = end_time - start_time;
    std::cout << "epoi " << n << " Elapsed time: " << elapsed_time.count() << "s" << std::endl;
    auto interval = algo::confidence_intervals(result,0.99);
    std::cout << interval.first << " -- " << ep << " -- " << interval.second << std::endl;
  }
}



int main()
{
  run_cut_test(); 
  return 0;
}