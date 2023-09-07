#include "Graph/graph.hpp"
#include<iostream>
#include <chrono>
#include "Algorithms/algorithms.hpp"


int main()
{
  graph::Graph test1("../data/company_dict_insiders.txt");
  //graph::Graph test2(105, 10, 0.5);
  // test1.printcsv("big.csv");

  graph::Graph bar = algo::greedy(5, test1);

  for (int n = 1000; n < 200000; n = n + 200)
  {
    auto clock_start = std::chrono::system_clock::now();
    double ep = test1.EPOI(n); 
    auto clock_now = std::chrono::system_clock::now();
    double ct = double(std::chrono::duration_cast <std::chrono::microseconds>(clock_now - clock_start).count());
    std::cerr << n << ": " << ct << "\n";
    std::cout << n << ", " << ep << "\n";
    std::cout.flush();
  } 
  return 0;
}