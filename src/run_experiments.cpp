#include "Graph/graph.hpp"
#include<iostream>

int main()
{
  graph::Graph test1("../data/company_dict_insiders.txt");
  //graph::Graph test2(105, 10, 0.5);
  // test1.printcsv("big.csv");
  for (int n = 1000; n < 200000; n = n + 200)
  {
    std::cout << n << ", " << test1.EPOI(n) << "\n";
    std::cout.flush();
  } 
  return 0;
}