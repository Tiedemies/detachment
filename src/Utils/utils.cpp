#include "utils.hpp"
#include <filesystem>
#include <iostream>
#include <fstream>

namespace util
{
  // A routine that iterates over all the files in a directory with a given extension and returns a vector of filenames
  std::vector<std::string> get_files(std::string path, std::string extension)
  {
    std::vector<std::string> files;
    for (const auto & entry : std::filesystem::directory_iterator(path))
    {
      std::string file = entry.path();
      if (file.find(extension) != std::string::npos)
      {
        files.push_back(file);
      }
    }
    return files;
  }

  // A routine that checks if a file exists
  bool file_exists(std::string filename)
  {
    std::ifstream infile(filename);
    return infile.good();
  }


}


