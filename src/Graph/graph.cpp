#include "graph.hpp"
#include<vector>
#include<iostream>
#include<fstream>
#include<string>
#include<unordered_map>
#include<set>
#include<algorithm>
#include<queue>
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
  Graph::Graph(int n, int r, double p): _num_nodes(n)
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
    DEBUG("Vertices initialized");
    // Rest are evenly distributed
    for(int i = r; i < n; ++i)
    {
      int j = static_cast<int>(rnd_circle.get());
      _circles[j].push_back(i);
    }
    DEBUG("Vertices distributed");
    // Then we generate a tree:
    for (int i = 1; i < r; ++i)
    {
      int j = ((int) (i*rnd_seed.get()));
      int index = ((int) (_circles[i].size()*rnd_seed.get()));
      int k = _circles[i][index];
      _circles[j].push_back(k);
    }
    DEBUG("Tree Generated");
    if (p > 0.001)
    {
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
          //DEBUG("created a connection \n");
        }
      }
    }
    instantiate();
    randomize();
  }

  // Copy constructor to make a copy
  Graph::Graph(const Graph& rhs): _num_nodes(rhs._num_nodes), _circles(rhs._circles), _circle_of_node(rhs._circle_of_node), _adjacency_vector(rhs._adjacency_vector), _node_places(rhs._node_places), _prob_vector(rhs._prob_vector)
  {
    // Void
  }
  
  Graph::~Graph(void)
  {
    // void
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
      // Iterate over the adjacent vertices
      for (int j = _node_places[u];j <_node_places[u+1];++j)
      {
        int v = _adjacency_vector[j];
        if (v < 0 || found[v]) continue;
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
      
  double 
  Graph::circle_epoi(int i, int sims, std::vector<double>& num_activated) const
  {
    bool verbose = (num_activated.size() == (size_t) sims); 
    BOOST_ASSERT(0 <= i);
    BOOST_ASSERT(i < _circles.size());
    double poi = 0.0;
    #pragma omp parallel for shared(poi,num_activated)
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
    return poi;
  } 

 // Calculate epoi, numbers only 
  double 
  Graph::EPOI_num(int n,bool probs, std::vector<double>& num_activated) const
  {
    int total_sims = 0;
    double epoi = 0.0;
    
    // If verbose, we use numbers activated 
    [[maybe_unused]] bool verbose = (num_activated.size() == (size_t) n); 
    
    DEBUG("start base simulation, verbose: " << verbose);
    for(int i =0; i < (int) _circles.size(); ++i)
    {
      /*
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
      */
      double poi = circle_epoi(i,n,num_activated);
      total_sims+=n;
      epoi += poi;
    }
    epoi = epoi/total_sims;
    return epoi;  
  }



std::vector<double> 
Graph::RunDMP(const std::vector<int>& inside) const
{
  const double epsilon = 0.0000001;
  const int n = _num_nodes;
  const double err_thres = n*epsilon;

  // Pi_j is as follows: Pi_j[j] contains (i, p') when p(i->j) = p' I.e. Pj_i[j] contains *incoming* 
  // messages from i to j 
  std::vector<std::vector<int>> Pi_j(n);
  std::vector<std::vector<int>> Pi_j_prev(n);
  std::unordered_map<size_t, double> messages;
  std::unordered_map<size_t, double> messages_prev;
  // Marginal probabilities
  std::vector<double> P0(n,0.0);
  std::vector<double> Marginals(n,0.0);
  std::vector<double> Marginals_prev(n,0.0);
  std::queue<int> Waiting;
  std::vector<bool> Wait(n,false);
  for (int u = 0; u < n; ++u)
  {
    // Instantiate inside
    if(std::find(inside.begin(), inside.end(), u) != inside.end())
    {
      P0[u] = 1.0; 
      Marginals_prev[u] = 1.0;
      // Iterate over adjacent vertices
      for (int j = _node_places[u];j <_node_places[u+1];++j)
      {
        int v = _adjacency_vector[j];
        if (v < 0) continue;
        messages_prev[Key(u,v)] = 1.0;
      } 
    }
    // Instantiate other
    else
    {
      // Iterate over adjacent vertices:
      for (int j = _node_places[u];j <_node_places[u+1];++j)
      {
        int v = _adjacency_vector[j];
        if (v < 0) continue;     
        messages_prev[Key(u,v)] = 0.0;
      }
    }
  }
  double err = (double) n;
  int iter = 0;
  while (err > err_thres && iter < n)
  {
    err = 0.0;
    ++iter; 
    for (int u = 0; u < n; ++u)
    {
      double p_co = 1-P0[u];
      // Iterate over adjacent vertices:
      for (int j = _node_places[u];j <_node_places[u+1];++j)
      {
        int v = _adjacency_vector[j];
        auto msp = messages_prev.find(Key(v,u));
        if (msp == messages_prev.end())
        {
          continue;
        }
        p_co*=(1.0 - (_prob_vector.at(j)*msp->second));
      }
      Marginals[u] = std::max(Marginals_prev[u], 1-p_co);
      err += std::abs(Marginals[u] - Marginals_prev[u]);
      Marginals_prev[u] = Marginals[u];
    }
    for (auto k_pair: messages_prev)
    {
      auto nodes = InvKey(k_pair.first);
      int u = nodes.first;
      int v = nodes.second;
      double m_p = 0.0;
      double temp_p = 0.0;
      int j = EdgeIndex(u,v); 
      BOOST_ASSERT(j >= 0);
      try
      {
        temp_p = 1.0 - _prob_vector[j]*messages_prev[Key(v,u)];
      }
      catch(const std::exception& e)
      {
        std::cerr << "Pmap or message fail" << '\n';
        std::cerr << "Kpair; " << k_pair.first << "\n";
        std::cerr << "v: " << v << ", u: " << u << "\n";
        throw e;
      }
      if (temp_p < epsilon)
      {
        m_p = 1.0 - P0[u];  
        for (int l = _node_places[v]; l < _node_places[v+1];++l)
        {
          int w = _adjacency_vector[l];
          // Skip backward edges
          if (w == u)
          {
            continue;
          }
          // If there is no connection or, if there is no message
          auto l_it = messages_prev.find(Key(w,v));
          int m = EdgeIndex(w,v);
          BOOST_ASSERT(m >= 0);
          if (_prob_vector[m] < epsilon ||  l_it == messages_prev.end())
          {
            continue;
          }
          m_p *= (1 - (_prob_vector[m]*l_it->second)); 
        }
        m_p = 1.0 - m_p; 
      }
      else
      {
        m_p = 1.0 - (1.0 - Marginals[u])/temp_p; 
      }
      messages[k_pair.first] = std::max(m_p, k_pair.second);
      err += std::abs(messages[k_pair.first] - k_pair.second);
    }
    messages_prev = messages;
  } 
  return Marginals;
}

  double
  Graph::DMP_EPOI(bool probs) const
  {
    // Use the DMP algorithm
    // Iterate over all circles:
    double epoi = 0.0;
    for (int i = 0; i < (int) _circles.size(); ++i)
    {
      // Create a vector of inside vertices:
      std::vector<int> inside = _circles.at(i);
      const auto marginals = RunDMP(inside);
      // Add together all marginals:
      double current_epoi = std::accumulate(marginals.begin(), marginals.end(), 0.0);
      // Subtract the marginals of the inside set:
      for (int j: inside)
      {
        current_epoi -= marginals[j];
      }
      // Divide by the number of outside vertices:
      current_epoi = current_epoi/((double) (_num_nodes - inside.size()));
      // Add to epoi:
      epoi += current_epoi;
    }
    return epoi/((double) _circles.size()); 
  }

  int 
  Graph::EdgeIndex(const int& u, const int& v) const
  {
    int index = -1;
    for (int i = _node_places[u]; i < _node_places[u+1]; ++i)
    {
      if (_adjacency_vector[i] == v)
      {
        index = i;
        break;
      }
    }
    return index;
  }
    
  // Detachment operator
  void 
  Graph::detach(int u, int C)
  {
    DEBUG("Detaching " << u << " from " << C);
    DEBUG("Detacment stack size: " << _detachment_stack.size());
    BOOST_ASSERT(C >= 0);
    BOOST_ASSERT(u >= 0);
    BOOST_ASSERT(C < _circles.size());
    BOOST_ASSERT(u < _num_nodes);
    /*
    auto c_iter = std::find(_circles.at(C).begin(),_circles.at(C).end(),u);
    BOOST_ASSERT(c_iter != _circles.at(C).end());
    auto cn_iter = std::find(_circle_of_node.aEPOI_numt(u).begin(),_circle_of_node.at(u).end(),C);
    BOOST_ASSERT(cn_iter != _circle_of_node.at(u).end());
    _circles.at(C).erase(c_iter);
    _circle_of_node.at(u).erase(cn_iter);
    */
    _detachment_stack.push_back(std::make_pair(u,C));
    for (int v: _circles[C])
    {
      for (int i = _node_places[v]; i < _node_places[v+1];++i)
      {
        if (_adjacency_vector[i] == u)
        {
          _adjacency_vector[i] = -u;
          break;
        }
      }
      for (int i = _node_places[u]; i < _node_places[u+1]; ++i)
      {
        if (_adjacency_vector[i] == v)
        {
          _adjacency_vector[i] = -v;
        }
      
      }
    }
  }

  void 
  Graph::undo()
  {
    if (_detachment_stack.empty())
    {
      return;
    }
    auto fpair = _detachment_stack.back();
    _detachment_stack.pop_back();
    int u = fpair.first;
    int C = fpair.second;
    for (int v: _circles[C])
    {
      if (v == u) continue;
      // If v is a detached vertex, then don't back down!
      if (std::find(_detachment_stack.begin(), _detachment_stack.end(), std::make_pair(v,C)) != _detachment_stack.end())
      {
        continue;
      }
      [[maybe_unused]] bool found = false;
      for (int i = _node_places[v]; i < _node_places[v+1];++i)
      {
        if (_adjacency_vector[i] == -u)
        {
          found = true;
          _adjacency_vector[i] = u;
          break;
        }
      }
      found = false;
      for (int i = _node_places[u]; i < _node_places[u+1]; ++i)
      {
        if (_adjacency_vector[i] == -v)
        {
          found = true;
          _adjacency_vector[i] = v;
          break;
        }
      }
    }
    /*
    _circles[C].push_back(u);
    _circle_of_node[u].push_back(C);
    */
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

  // Inverse key for the probability matrix
  std::pair<int,int>
  Graph::InvKey(const size_t& k) const
  {
    int u = (int) (k >> 32);
    int v = (int) (k & 0xFFFFFFFF);
    return std::make_pair(u,v);
  }

  void
  Graph::instantiate()
  {
    DEBUG("Instantiationg adjacency lists for " << _num_nodes << " nodes");
    _circle_of_node.clear();
    _circle_of_node.resize(_num_nodes);
    std::vector<std::vector<int>> adjacent(_num_nodes); 
    for (size_t i = 0; i < _circles.size(); ++i)
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
          adjacent[u].push_back(v);
        }
      }
    }
    _adjacency_vector.clear();
    _node_places.clear();
    _num_edges = 0;
    for (int u = 0; u < _num_nodes; ++u)
    {
      _node_places.push_back(_num_edges);
      for (int v: adjacent[u])
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
    _prob_vector.resize(_num_edges,std::nan(""));
    for (int u = 0; u < _num_nodes;++u)
    {
      size_t end = (u < _num_nodes - 1)?_node_places[u+1]:_num_edges;
      for(int j = _node_places[u]; j < (int) end; ++j)
      {
        // int v = _adjacency_vector[j];
        _prob_vector[j] = rnd.get();
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
  Graph::printcsv(std::string fname, bool directed) const
  {
    std::ofstream ofs(fname);
    ofs << "source, target, weight \n";
    for (int u = 0; u < _num_nodes; ++u)
    {
      for(int i = _node_places[u];i < _node_places[u+1]; ++i)
      {        
        int v = _adjacency_vector[i];
        if (!directed && v > u)
        {
          continue;
        }
        double p = _prob_vector[i];
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
        for (int i = _node_places[s]; i < _node_places[s+1]; ++i)
        {
          int t = _adjacency_vector[i];
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
    for (size_t i = 0; i < main_component.size(); ++i)
    {
      num_map[main_component[i]] = (int) i;
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
    DEBUG("Circles after clean up:" <<_circles.size());
  }

  Graph
  Graph::clean_copy() const
  {
    Graph newg;
    newg._num_nodes = _num_nodes;
    newg._num_edges = 0;
    newg._adjacency_vector.clear();
    newg._prob_vector.clear();
    newg._node_places.clear();
    newg._circles = _circles;
    newg._circle_of_node = _circle_of_node;
    for (auto d_pair: _detachment_stack)
    {
      int u = d_pair.first;
      int C = d_pair.second;
      auto c_iter = std::find(newg._circles.at(C).begin(),newg._circles.at(C).end(),u);
      if (c_iter == newg._circles.at(C).end())
      {
        std::cerr << "Circle: " << C << ", removing " << u << "\n";
        for (auto v: newg._circles.at(C))
        {
            std::cerr << v << " ";
        }
        std::cerr << "\n";
      }
      BOOST_ASSERT(c_iter != newg._circles.at(C).end());
      auto cn_iter = std::find(newg._circle_of_node.at(u).begin(),newg._circle_of_node.at(u).end(),C);
      BOOST_ASSERT(cn_iter != newg._circle_of_node.at(u).end());
      newg._circles.at(C).erase(c_iter);
      newg._circle_of_node.at(u).erase(cn_iter);
    }

    for (int u = 0; u < _num_nodes; ++u)
    {
      newg._node_places.push_back(newg._num_edges);
      for (int i = _node_places[u]; i < _node_places[u+1]; ++ i)
      {
        int v = _adjacency_vector[i];
        if (v < 0)
        {
          continue;
        }
        newg._adjacency_vector.push_back(v);
        newg._prob_vector.push_back(_prob_vector[i]);
        ++newg._num_edges;
      }
    }
    newg._node_places.push_back(newg._num_edges);
    return newg; 
  }

  void
  Graph::validate()
  { 
    DEBUG("Validation passed");
  }

  bool Graph::is_adjacent(int u, int v) const
  {
    for (int b = _node_places[u]; b < _node_places[v]; ++b)
    {
      if (_adjacency_vector[b] == v)
      {
        return true;
      }
    }
    return false;
  }

  int Graph::Get_num_nodes() const
  {
    return _num_nodes;
  }

}
