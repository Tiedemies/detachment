#ifndef TESTS_HPP
#define TESTS_HPP
#include<string>

// The header file for the tests, for each function in tests.cpp there is a declaration here
// The functions are defined in tests.cpp
namespace test
{
  // Generate the graphs for the experiments
  void generate_graphs();

  /// Run experiments for different algorithms
  void Run_comparative_tests(int max_nodes = -1);

  /// Run a single test for the company data
  void Run_company_test();

  /// Test the cut algorithm
  void run_cut_test();

  /// Test the greedy algorithm
  void run_test();

  /// Test speed
  void speed_test();

  // A debugging test
  void run_test2();

  void Run_test(std::string filename, std::string result_handle, int iterations = 100);
}

#endif // TESTS_HPP