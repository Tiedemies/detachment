#include "Graph/graph.hpp"
#include<iostream>
#include <chrono>
#include "Algorithms/algorithms.hpp"
#include "Utils/utils.hpp"
#include "Tests/tests.hpp"
#include <sstream>
#include <iomanip>
#include <string>
#include <fstream>





int main()
{
  //Run the comparison test
  // test::Run_comparative_tests(300);
  test::Run_test("../data/empirical.txt","../data/empirical", 10000);
  return 0;
}