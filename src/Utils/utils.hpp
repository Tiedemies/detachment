#ifndef UTILS_HPP
#define UTILS_HPP 


#include <filesystem>
#include <iostream>

namespace util
{
  // A routine that iterates over all the files in a directory with a given extension and returns a vector of filenames
  std::vector<std::string> get_files(std::string path, std::string extension);

  // A routine that checks if a file exists
  bool file_exists(std::string filename);

}

#endif

#ifdef NDEBUG
#define DEBUG(x)
#else
#define DEBUG(x) do {std::cerr << x << std::endl << std::flush; } while(0)
#endif

