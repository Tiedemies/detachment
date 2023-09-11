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
  std::cerr << "bar";
  double avg = 0.0;
  double div = 0;
  for (int n = 100; n < 20000; n = n + 100)
  {
    auto clock_start = std::chrono::system_clock::now();
    std::vector<double> result(n,0.0);
    double ep = test1.EPOI(n,false,result); 
    auto clock_now = std::chrono::system_clock::now();
    double ct = double(std::chrono::duration_cast <std::chrono::microseconds>(clock_now - clock_start).count());
    //std::cerr << n << ": " << ct << "\n";
  
    auto interval = algo::confidence_intervals(result,0.95);
    std::cout << n << ", " << interval.first << "," << ep << "," << interval.second << "\n";
    std::cout.flush();
  } 
  return 0;
}