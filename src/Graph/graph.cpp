#include "graph.hpp"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>
#include<set>
#include<algorithm>
#include "../Utils/grandom.h"

namespace graph
{
  // Default constructor 
  Graph::Graph(void)
  {
    // void
  }
  // Read from named file
  Graph::Graph(std::string fname)
  {
    std::ifstream inp(fname);
    
    int cnum = 0;
    std::string temp;
    inp >> cnum;
    _num_nodes = 0; 
    while(!inp.eof() && inp.good())
    {
      int inum;
      inp.get();
      char x = inp.get();
      std::vector<int> insiders; 
      while(x != '\n')
      {
        if (x == '[' || x == ',')
        {
          inp >> inum;
          if(inum + 1  > _num_nodes)
          {
            _num_nodes = inum + 1;
          }
          insiders.push_back(inum);  
        }
        x = inp.get();
      }
      _circles.push_back(insiders);
      inp >> cnum;
    }
    instantiate();
    clean_up();
    instantiate(); 
    randomize();
  }

  // Random graph generation. 
  // Model: n vertices, r circles, uniform random (one vertex per circle guaranteed). 
  // 0 <= p < 1 is a connection coefficient. p == 0 results in a block graph, theoretically p=1 is fully connected. 
  Graph::Graph(int n, int r, double p)
  {
    BOOST_ASSERT(n > r);
    BOOST_ASSERT(0 <= p && p < 0.95);
     _circles.resize(r);
    [[maybe_unused]] int left = n;
    util::Random rnd_circle(r);
    util::Random rnd_seed;  
    // Every circle has a vertex
    for (int i = 0; i < r; ++i)
    {
      _circles[i].push_back(i);
    }
    // Rest are evenly distributed
    for(int i = r; i < n; ++i)
    {
      int j = (int) rnd_circle.get();
      _circles[j].push_back(i);
    }
    // Then we generate a tree:
    for (int i = 1; i < r; ++i)
    {
      int j = ((int) (i*rnd_seed.get()));
      int index = ((int) (_circles[i].size()*rnd_seed.get()));
      int k = _circles[i][index];
      _circles[j].push_back(k);
    }
    for (int i = 0; i < r; ++i)
    {
      while (rnd_seed.get() < p)
      {
        int j = (int) rnd_circle.get();
        if (i == j)
        {
          continue;
        }
        int index = ((int) (_circles[i].size()*rnd_seed.get()));
        int k = _circles[i][index];
        if (find(_circles[j].begin(), _circles[j].end(),k) != _circles[j].end())
        {
          continue;
        }
        _circles[j].push_back(k);
        std::cerr << "created 1\n";
      }
    }
    instantiate();
    randomize();
  }

  // Copy constructor to make a copy
  Graph::Graph(const Graph& rhs)
  {
    _circles = rhs._circles;
    _probs = rhs._probs;
    _adjacent = rhs._adjacent;
  }
  
  Graph::~Graph(void)
  {
    // void
  }

  // Singleton spread. 
  std::vector<bool>
  Graph::spread(int in) const
  {
    // Precondition
    BOOST_ASSERT(in >= 0);
    BOOST_ASSERT((size_t) in < _circles.size());
  
    std::vector<bool> output(_num_nodes,false);
    if (_circles[in].empty())
    {
      return output; 
    }

    std::vector<int> infected;
    infected.reserve(_num_nodes); 
    std::copy(_circles.at(in).begin(), _circles.at(in).end(), std::back_inserter(infected));
    for (int i: infected)
    {
      output[i] = true;
    }
    util::Random rnd;
    while(infected.size() > 0)
    {
      int u = infected.back();
      infected.pop_back();
      for (int v: _adjacent.at(u))
      {
        if (output[v]) continue;
        if (rnd.get() < _probs.at(Key(u,v)))
        {
          output[v] = true;
          infected.push_back(v);
        }
      }
    }
    return output; 
  }

  // Calculate activation probabilities for circle i
  std::vector<double> 
  Graph::activation_probs(int i, int n) const
  {
    BOOST_ASSERT(0 <= i);
    BOOST_ASSERT((size_t) i < _circles.size());
    std::vector<double> activations(_num_nodes,0.0);
    #pragma omp parallel for shared(activations)
    for(int j = 0; j < n; ++j)
    {
      auto res = spread(i);
      for (int k = 0; k < _num_nodes; ++k)
      {
        if (res[k])
        {
          #pragma omp critical
          {
            activations.at(k) += 1.0/n;   
          }
        }
      }
    }
    return activations; 
  }
      
  // Calculate epoi 
  double 
  Graph::EPOI(int n,bool probs) const
  {
    int total_sims = 0;
    double epoi = 0.0;
    for(int i =0; i < (int) _circles.size(); ++i)
    {
      int sims = (int) (probs?(_initial_probabilities.at(i) * n): (n/_circles.size()));
      total_sims += sims; 
      double poi = 0.0;
      #pragma omp parallel for shared(epoi)
      for (int j = 0; j < sims; ++j)
      {
        auto res = spread(i);
        int nn = (std::count(res.begin(), res.end(), true)) - _circles.at(i).size();
        #pragma omp critical
        {
          poi += ((double) nn)/((double) (_num_nodes - _circles.at(i).size())); 
        }
      }
      epoi = epoi+poi;
    }
    epoi = epoi/total_sims;
    return epoi;  
  } 
  
  // Detachment operator
  Graph 
  Graph::detach(int u, int C) const
  {
    // void 
  }

  // Probability (defunct now)
  double 
  Graph::get_prob(int u,int v) const
  {
    double foo = 0;
    return foo;
  } 

  // Key for the probability matrix
  size_t
  Graph::Key(const int& u, const int& v) const
  {
    return (size_t) (( (size_t) u << 32) | (unsigned int) v);
  }


  void
  Graph::instantiate()
  {
    _adjacent.clear();
    _adjacent.resize(_num_nodes);
    for (auto c: _circles)
    {
      for (int u: c)
      {
        for (int v: c)
        {
          if (u == v)
          {
            continue;
          }
          _adjacent[u].insert(v);
          _adjacent[u].insert(v);
        }

      }
    }
  }

  void
  Graph::randomize(double alpha, double beta)
  {
    util::BetaRandom rnd(alpha,beta);
    for (int i = 0; i < _num_nodes; ++i)
    {
      for (int j: _adjacent[i])
      {
        size_t key = Key(i,j);
        _probs[key] = rnd.get();
      }
    }
  }

  void 
  Graph::print(std::string fname)
  {
    std::ofstream ofs(fname);
    for (unsigned int i = 0; i < _circles.size(); ++ i)
    {
      ofs << i << "[";
      auto citer = _circles[i].begin();
      for (;citer != _circles[i].end()-1; ++citer)
      {
        ofs << *citer << ",";
      }
      ofs << *citer << "]\n";
    }
  }
  void 
  Graph::printcsv(std::string fname, bool directed)
  {
    std::ofstream ofs(fname);
    ofs << "source, target, weight \n";
    for (int u = 0; u < _num_nodes; ++u)
    {
      for(auto v: _adjacent[u])
      {
        if (!directed && v > u)
        {
          continue;
        }
        double p = 0.1;
        size_t key = Key(u,v);
        auto c_it = _probs.find(key);
        if (c_it != _probs.end())
        {
          p = c_it->second; 
        }
        ofs << u << ", " << v << ", " << p << "\n";
      }
    }
  }

  // Make it all a one big component.
  void 
  Graph::clean_up()
  {
    std::set<int> found;
    std::vector<int> main_component;
    std::cerr << "Number of nodes: " << _num_nodes << "\n";
    for (int i= 0; i < _num_nodes; ++i)
    {
      if (found.find(i) != found.end()) continue;
      std::vector<int> comp = {i};
      std::vector<int> Stack = {i};
      while (!Stack.empty())
      {
        int s = Stack.back();
        Stack.pop_back();
        for (int t: _adjacent[s])
        {
          auto t_it = found.find(t);
          if (t_it != found.end())
          {
            continue;
          }
          Stack.push_back(t);
          comp.push_back(t);
          found.insert(t);
        }
      }
      if (comp.size() > main_component.size())
      {
        std::cerr << "found component of size " << comp.size() << "\n"; 
        main_component = std::move(comp);
      }
    }
    // main_component now contains all. 
    std::unordered_map<int,int> num_map;
    std::set<int> main_cset;
    for (int i = 0; i <= main_component.size(); ++i)
    {
      num_map[main_component[i]] = i;
    }
    std::vector<std::vector<int>> new_circles; 
    for (auto c: _circles)
    {
      if (c.empty() ) continue;
      int s = c.front();
      auto s_it = num_map.find(s);
      if (s_it == num_map.end()) continue;
      std::vector<int> circle;
      circle.reserve(c.size());
      for (int u: c)
      {
        int s = num_map[u];
        circle.push_back(s); 
      }
      new_circles.push_back(circle);
    }
    _num_nodes = main_component.size(); 
    _circles = std::move(new_circles);
  }
}
