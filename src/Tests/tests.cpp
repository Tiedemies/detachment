#include"tests.hpp"
#include"../Graph/graph.hpp"
#include<iostream>
#include <chrono>
#include"../Algorithms/algorithms.hpp"
#include"../Utils/utils.hpp"
#include <sstream>
#include <iomanip>
#include <string>
#include <fstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <random>
#include <map>
#include <set>

namespace test
{
  // Run single test for the company data
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
    graph::Graph test1("../data/100_20_0.30.txt");
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
      graph::Graph test1("..");

      auto start_time = std::chrono::system_clock::now();
      double ep = test1.EPOI_num(n,false,result);
      auto end_time = std::chrono::system_clock::now();
      std::chrono::duration<double> elapsed_time = end_time - start_time;
      std::cout << "epoi " << n << " Elapsed time: " << elapsed_time.count() << "s" << std::endl;
      auto interval = algo::confidence_intervals(result,0.99);
      std::cout << interval.first << " -- " << ep << " -- " << interval.second << std::endl;
    }
  }

  // Generate the graphs for the experiments
  void generate_graphs()
  {
    std::vector<int> node_n = {50, 100, 300, 1000};
    std::vector<int> circle_n = {50,10,5};
    std::vector<double> p_s = {0.0, 0.3, 0.5, 0.8};
    DEBUG("Init...");
    for (int nodes: node_n)
    {
      for (int c_coeff: circle_n)
      {
        int circles = static_cast<int>( ((double) nodes) / ((double) c_coeff));
        for (double p: p_s)
        {
          std::stringstream stream;
          stream << std::fixed << std::setprecision(2) << p;
          std::string s = stream.str();
          std::string filename = "../data/" + std::to_string(nodes) + "_" + std::to_string(circles) + "_" + s + ".csv";
          std::string filename2 = "../data/" + std::to_string(nodes) + "_" + std::to_string(circles) + "_" + s + ".txt";
          
          // file_exists is a flag that is true if the file already exists
          bool file_exists = util::file_exists(filename);
          if (file_exists)
          {
            DEBUG("File " << filename << " already exists, skipping...");
            continue;
          }

          DEBUG("Generate " << nodes << "," << circles << "," << p);
          graph::Graph nn(nodes,circles,p);
          DEBUG("Generated, now saving...");
          nn.printcsv(filename);
          nn.print(filename2);
        }
      }
    }
  }

  /// Run experiments for different algorithms
  void Run_comparative_tests(int max_nodes)
  {
    // Iterate over all files with .txt extension in data folder
    std::vector<std::string> files = util::get_files("../data/", ".txt");
    for (auto file: files)
    {
      // Load graph
      graph::Graph g("../data/" + file);
      if (max_nodes > 0 && g.Get_num_nodes() > max_nodes)
      {
        continue;
      }
      DEBUG("Running test for " << file);
      // Run cut optimizer
      algo::Cutter cut(g);
      cut.set_iterations(10000);
      cut.optimize();
      int n = cut._detached.size();     

      // Run greedy optimizer to cut n connections  
      algo::GreedyOptimizer greedy(g);
      greedy.optimize(n);


      // Save results, but preface the filename with the algorithm name and replace the extension with .csv
      std::string cutfilename = "../results/cut_" + file.substr(8, file.size() - 11) + "csv";
      std::string greedyfilename = "../results/greedy_" + file.substr(8, file.size() - 11) + "csv";
      std::string originalfilename = "../results/original_" + file.substr(8, file.size() - 11) + "csv";
      // Print the results to the files
      g.printcsv(originalfilename);
      cut.get_result().printcsv(cutfilename);
      greedy.get_result().printcsv(greedyfilename);
      // Save the result epoi and the base epoi to a file with the same name but with .num extension
      std::string numfilename = "../results/" + file.substr(8, file.size() - 11) + "num";
      std::ofstream numfile(numfilename);
      // Prevace the result with n, the number of connections cut
      numfile << n << "," << cut._result_epoi << "," << cut._base_epoi << "," << greedy._result_epoi << "," << greedy._base_epoi << "\n";
    
    }
  }

  // Run a debug test
  void run_test2()
  {
    DEBUG("Run the 100 vertex test");
    // Iterate over all the graph files in the data folder that have 100 vertices:
    std::vector<std::string> files = util::get_files("../data/", ".txt");
    for (auto file: files)
    {
      // Skip the file unless it starts with "100_" 
      if (file.substr(8,7) != "100_2_0")
      {
        DEBUG("Skipping " << file);
        continue;
      }
      DEBUG("***** Running test for " << file);
      // Load the graph
      graph::Graph g("../data/" + file);
      // Run the cut optimizer  
      algo::Cutter cut(g);
      cut.optimize();
      DEBUG("Cut optimizer complete, result epoi: " << cut._result_epoi);
      DEBUG("The cut is: ");
      #ifndef NDEBUG
      for (auto cpair: cut._detached)
      {
        std::cerr << "(" << cpair.first << "," << cpair.second << ")";
      }
      std::cerr << std::endl;
      #endif

      int n = cut._detached.size();
      // Run the greedy optimizer
      algo::GreedyOptimizer greedy(g);
      greedy.optimize(n);

      // output the resuts:

    }
  }
}