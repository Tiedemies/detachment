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
  for (int n = 1000; n < 200000; n = n + 200)
  {
    auto clock_start = std::chrono::system_clock::now();
    double ep = test1.numeric_EPOI(n); 
    auto clock_now = std::chrono::system_clock::now();
    double ct = double(std::chrono::duration_cast <std::chrono::microseconds>(clock_now - clock_start).count());
    std::cerr << n << ": " << ct << "\n";

    avg = (div*avg + ep) / (div + 1);
    ++div; 
    std::cout << n << ", " << ep << "," << avg << "\n";
    std::cout.flush();
  } 
  return 0;
}