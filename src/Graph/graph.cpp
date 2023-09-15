#include "graph.hpp"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>
#include<set>
#include<algorithm>
#include "../Utils/grandom.h"
#include "../Utils/utils.hpp"

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
    DEBUG("instantiate");
    instantiate();
    DEBUG("clean up");
    clean_up();
    DEBUG("instantiate again");
    instantiate(); 
    DEBUG("randomize");
    randomize();
    #ifndef NDEBUG
    validate();
    #endif
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
        DEBUG("created 1\n");
      }
    }
    instantiate();
    randomize();
  }

  // Copy constructor to make a copy
  Graph::Graph(const Graph& rhs): _num_nodes(rhs._num_nodes), _circles(rhs._circles), _adjacent(rhs._adjacent), _circle_of_node(rhs._circle_of_node),_probs(rhs._probs),_adjacency_vector(rhs._adjacency_vector), _node_places(rhs._node_places),_prob_vector(rhs._prob_vector)
  {
    // Void
  }
  
  Graph::~Graph(void)
  {
    // void
  }

  // Singleton spread. 
  std::vector<bool>
  Graph::spread(int in, std::unordered_map<size_t,double>& previous_results) const
  {
    bool clear = previous_results.empty();
    // Precondition
    BOOST_ASSERT(in >= 0);
    BOOST_ASSERT((size_t) in < _circles.size());
  
    std::vector<bool> output(_num_nodes,false);
    if (_circles.at(in).empty())
    {
      return output; 
    }

    std::vector<int> infected;
    infected.reserve(_num_nodes); 
    std::copy(_circles.at(in).begin(), _circles.at(in).end(), std::back_inserter(infected));
    for (int i: infected)
    {
      try
      {
        output.at(i) = true;
      }
      catch(const std::exception& e)
      {
        DEBUG(e.what() << " at infection loop");
        DEBUG("in: " << in);
        DEBUG(_circles.at(in).size());
        DEBUG(infected.size());
        throw e;
      }
      
    }
    util::Random rnd;
    while(infected.size() > 0)
    {
      int u = infected.back();
      infected.pop_back();
      for (int v: _adjacent.at(u))
      {
        try
        {
          if (output.at(v)) continue;
        }
        catch(const std::exception& e)
        {
          DEBUG(e.what() << " at adjacency loop");
          DEBUG( "u:" << u);
          DEBUG( "v:" << v);
          DEBUG( "output size " << output.size());
;          
          throw e; 
        }
        
        size_t k = Key(u,v); 
        double rndn =0.0;
        if (!clear)
        {
          auto pr_it = previous_results.find(k);
          if (pr_it != previous_results.end())
          {
            rndn = pr_it->second;
          }
          else
          {
            rndn = rnd.get();
          }  
        }
        else 
        {
          rndn = rnd.get();
          previous_results[k] = 1-rndn; 
        }
        if (rndn < _probs.at(Key(u,v)))
        {
          output.at(v) = true;
          infected.push_back(v);
        }
      }
    }
    if (!clear)
    {
      previous_results.clear();
    }
    return output; 
  }

  // Singleton spread. 
  int
  Graph::spread_num(int in, std::vector<double>& previous) const
  {
    bool clear = previous.empty();
    // Precondition
    BOOST_ASSERT(in >= 0);
    BOOST_ASSERT((size_t) in < _circles.size());
    int output = 0;
    if (_circles.at(in).empty())
    {
      return output; 
    }
    std::vector<bool> found(_num_nodes,false);
    std::vector<int> infected;
    infected.reserve(_num_nodes); 
    std::copy(_circles.at(in).begin(), _circles.at(in).end(), std::back_inserter(infected));
    output = (int) infected.size();
    for (int i: infected)
    {
      found[i] = true;
    }
    util::Random rnd;
    if (clear)
    {
      previous.resize(_num_edges,std::nan(""));
    }
    while(infected.size() > 0)
    {
      int u = infected.back();
      infected.pop_back();
       for (int j = _node_places[u];j <_node_places[u+1];++j)
      {
        int v = _adjacency_vector[j];
        if (found[v]) continue;
        double die = (clear || std::isnan(previous[j]))?rnd.get():previous[j];
        if (die < _prob_vector[j])
        {    
          ++output;
          found[v] = true;
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
      std::unordered_map<size_t,double> dies; 
      auto res = spread(i,dies);
      auto res2 = spread(i,dies);
      for (int k = 0; k < _num_nodes; ++k)
      {
        if (res[k])
        {
          #pragma omp critical
          {
            activations.at(k) += 0.5/n;   
          }
        }
        if (res2[k])
        {
           #pragma omp critical
          {
            activations.at(k) += 0.5/n;   
          }

        }
      }
    }
    return activations; 
  }
      
  // Calculate epoi 
  double 
  Graph::EPOI(int n,bool probs, std::vector<double>& num_activated) const
  {
    int total_sims = 0;
    double epoi = 0.0;
    
    // If verbose, we use numbers activated 
    bool verbose = (num_activated.size() == (size_t) n); 
    
    DEBUG("start base simulation, verbose: " << verbose);
    for(int i =0; i < (int) _circles.size(); ++i)
    {
      int sims = n;
      total_sims += sims;
      double poi = 0.0;
      #pragma omp parallel for shared(epoi,num_activated)
      for (int j = 0; j < sims; ++j)
      {
        std::unordered_map<size_t,double> results;
        auto res = spread(i,results);
        double nn1 = (double) (std::count(res.begin(), res.end(), true)) - _circles.at(i).size();
        res = spread(i,results);
        double nn2 = (double) (std::count(res.begin(), res.end(), true)) - _circles.at(i).size();
        double nn = nn1+nn2;
        double marginal = (nn/2)/((double) (_num_nodes - _circles.at(i).size())); 
        #pragma omp critical
        {
          poi +=  marginal;
          //DEBUG(nn);
        }
        if (verbose)
        {
          #pragma omp critical
          {
            num_activated.at(j) += marginal/_circles.size();
          }
        }
      }
      epoi += poi;
    }
    epoi = epoi/total_sims;
    return epoi;  
  }

 // Calculate epoi, numbers only 
  double 
  Graph::EPOI_num(int n,bool probs, std::vector<double>& num_activated) const
  {
    int total_sims = 0;
    double epoi = 0.0;
    
    // If verbose, we use numbers activated 
    bool verbose = (num_activated.size() == (size_t) n); 
    
    DEBUG("start base simulation, verbose: " << verbose);
    for(int i =0; i < (int) _circles.size(); ++i)
    {
      int sims = n;
      total_sims += sims;
      double poi = 0.0;
      #pragma omp parallel for shared(epoi,num_activated)
      for (int j = 0; j < sims; ++j)
      {
        std::vector<double> results;
        double res = spread_num(i,results);
        double nn1 = res - _circles.at(i).size();
        res = spread_num(i,results);
        double nn2 = res - _circles.at(i).size();
        double nn = nn1+nn2;
        double marginal = (nn/2)/((double) (_num_nodes - _circles.at(i).size())); 
        #pragma omp critical
        {
          poi +=  marginal;
          // DEBUG(nn);
        }
        if (verbose)
        {
          #pragma omp critical
          {
            num_activated.at(j) += marginal/_circles.size();
          }
        }
      }
      epoi += poi;
    }
    epoi = epoi/total_sims;
    return epoi;  
  }

    
  // Detachment operator
  Graph 
  Graph::detach(int u, int C) const
  {
    Graph other = *this;
    auto conit = std::find(other._circle_of_node[u].begin(), other._circle_of_node[u].end(), C);
    if (conit != other._circle_of_node[u].end()) //the element was not found
    {
      other._circle_of_node[u].erase(conit);
    }
    auto cinit = std::find(other._circles[C].begin(), other._circles[C].end(), u);
    if (cinit != other._circles[C].end()) //the element was not found
    {
      other._circles[C].erase(cinit);
    }
    for (int v: other._circles[C])
    {
      auto nodit = std::find(other._adjacent[v].begin(), other._adjacent[v].end(),u);
      if (nodit != other._adjacent[v].end())
      {
        other._adjacent[v].erase(nodit);
      }
      auto adit = std::find(other._adjacent[u].begin(), other._adjacent[u].end(), v);
      if (adit != other._adjacent[u].end())
      {
        other._adjacent[u].erase(adit);
      }
    }
    return other;
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
    _circle_of_node.clear();
    _circle_of_node.resize(_num_nodes);
    for (int i = 0; i < _circles.size(); ++i)
    {
      auto c = _circles[i];
      for (int u: c)
      {
        // std::cerr << "u;" << u << ", i;" << i << "\n";
        _circle_of_node[u].push_back(i);
        // std::cerr << "done\n";
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
    _adjacency_vector.clear();
    _node_places.clear();
    _num_edges = 0;
    for (int u = 0; u < _num_nodes; ++u)
    {
      _node_places.push_back(_num_edges);
      for (int v: _adjacent[u])
      {
        _adjacency_vector.push_back(v);
        ++_num_edges;
      }
    }
    _node_places.push_back(_num_edges);
  }

  void
  Graph::randomize(double alpha, double beta)
  {
    util::BetaRandom rnd(alpha,beta);
    for (int i = 0; i < _num_nodes; ++i)
    {
      for (int j: _adjacent[i])
      {
        double p = rnd.get();
        size_t key = Key(i,j);
        _probs[key] = p;
      }
    }
    _prob_vector.resize(_num_edges,std::nan(""));
    for (int u = 0; u < _num_nodes;++u)
    {
      size_t end = (u < _num_nodes - 1)?_node_places[u+1]:_num_edges;
      for(int j = _node_places[u]; j < end; ++j)
      {
        int v = _adjacency_vector[j];
        size_t key = Key(u,v);
        _prob_vector[j] = _probs[key];
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
    for (int i = 0; i < main_component.size(); ++i)
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

  // Singleton spread. 
  std::vector<double>
  Graph::numeric_spread(int in, std::unordered_map<size_t,double>& previous_results) const
  {
    bool clear = previous_results.empty();
    // Precondition
    BOOST_ASSERT(in >= 0);
    BOOST_ASSERT((size_t) in < _circles.size());
  
    std::vector<double> output(_num_nodes,0.0);
    std::vector<bool> found(_num_nodes,false); 
    if (_circles[in].empty())
    {
      return output; 
    }

    std::vector<int> infected;
    infected.reserve(_num_nodes); 
    std::copy(_circles.at(in).begin(), _circles.at(in).end(), std::back_inserter(infected));
    for (int i: infected)
    {
      output[i] = 1.0;
      found[i] = 0.0;
    }
    util::Random rnd;
    while(infected.size() > 0)
    {
      int u = infected.back();
      infected.pop_back();
      for (int v: _adjacent.at(u))
      {
        size_t k = Key(u,v); 
        double rndn;
        if (!clear)
        {
          auto pr_it = previous_results.find(k);
          if (pr_it != previous_results.end())
          {
            rndn = pr_it->second;
          }
          else
          {
            rndn = rnd.get();
          }  
        }
        else 
        {
          rndn = rnd.get();
          previous_results[k] = 1-rndn; 
        }
        output[v] = std::min(1.0, output[v] + (1-output[v])*(_probs.at(k)*output[u]));
        if (rndn < _probs.at(k))
        {
          if (!found[v])
          {
            infected.push_back(v);
            found[v] = true; 
          }
        }
      }
    }
    if (clear)
    {
      previous_results.clear();
    }
    return output; 
  }

  void
  Graph::validate()
  { 
    for (int i = 0; i < _num_nodes;++i)
    {
      int u = i;
      std::set<int> neighbours_a;
      for (int v: _adjacent[u])
      {
        neighbours_a.insert(v);
      }
      std::set<int> neighbours_b;
      for (int j = _node_places[i]; j < _node_places[i+1]; ++j)
      {
        int v = _adjacency_vector[j];
        double prob1 = _probs[Key(u,v)];
        double prob2 = _prob_vector[j];
        BOOST_ASSERT(std::abs(prob1-prob2) < 0.000001);
        neighbours_b.insert(v);
      }
      BOOST_ASSERT(neighbours_a == neighbours_b);
    }
    DEBUG("Validation passed");

  }
}
