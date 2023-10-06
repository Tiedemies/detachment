#include "algorithms.hpp"
#include<numeric>
#include<math.h>
#include<queue>
#include<set>
#include<iostream>
#include "../Graph/graph.hpp"
#include "../Utils/utils.hpp"
#include<boost/math/distributions/normal.hpp>

namespace algo
{

  static double previous_mu = 0.0;
  static double previous_sigma = 0.0;

  const graph::Graph& 
  Optimizer::get_result() const
  {
    return _result;
  }

  Optimizer::Optimizer(const graph::Graph& g): _parent(g)
  {
    // void
  }

  void
  Optimizer::set_iterations(int sims)
  {
    _sims = sims;
  }

  GreedyOptimizer::GreedyOptimizer(const graph::Graph& g) : Optimizer(g)
  {
    // void
  }

  Cutter::Cutter(const graph::Graph& g) : Optimizer(g)
  {
    // void
  }


  inline double mean(const std::vector<double>& input)
  {
    if (input.size() == 0)
    {
      return 0.0;
    }
    return std::accumulate(input.begin(),input.end(),0.0)/static_cast<double>(input.size()); 
  }


  inline double variance(const std::vector<double>& input)
  {
    double mu = mean(input);
    return variance(input,mu);
  }

  inline double variance(const std::vector<double>& input, const double& mu)
  {
    size_t n = input.size();
    if (n <= 1)
    {
      return 0;
    }
    auto var_func = [&mu, &n](double accumulator, const double& val)
    { 
      return accumulator + (val - mu)*(val - mu)/(n-1); 
    };
    return std::accumulate(input.begin(),input.end(),0.0,var_func); 
  }

  std::pair<double,double> confidence_intervals(const std::vector<double>& input, const double& p)
  {
    const double mu = mean(input);
    DEBUG("average: " << mu);
    return confidence_intervals(input,p,mu);
  }

  std::pair<double,double> confidence_intervals(const std::vector<double>& input, const double& p, const double& mu)
  {
    size_t n = input.size();
    double s2 = variance(input,mu);
    double sigma = std::sqrt(s2);
    DEBUG("mu: " << mu << " sigma: " << sigma);
    DEBUG("mu-difference: " << mu - previous_mu << ", sigma-difference: " << sigma - previous_sigma);
    previous_mu = mu; 
    previous_sigma =  sigma;
    boost::math::normal_distribution nd;
    double z = boost::math::quantile(nd, 1.0 - ( (1-p)/2.0)); 
    double delta = std::sqrt(s2 / static_cast<double>(n))  * z;
    DEBUG("delta: " << delta << "\n");
    return std::make_pair(mu-delta,mu+delta);
  }

  /// @brief Optimize by finding the best k-detachment with a greedy algorithm
  /// @param k 
  void
  GreedyOptimizer::optimize(int k)
  {
    const int base_sim = 10000;
    // Initialize a copy.
    const int n = _parent._num_nodes;
    DEBUG("working with " << n << " nodes");
    std::vector<double> result(2*base_sim,0.0);
    DEBUG("establishing baseline");
    _base_epoi = _parent.EPOI_num(2*base_sim, false, result);
    DEBUG("calculating intervals");
    const auto ref_interval = confidence_intervals(result,0.99);
    double cur_min = _base_epoi; 
    DEBUG("base interval established " << ref_interval.first << " -- " << _base_epoi << " -- " << ref_interval.second);
    _result = _parent;
    graph::Graph temporary_min = _parent;
    _detached.clear(); 
    for (int i = 0; i < k; ++i)
    {
      std::pair<int,int> rpair; 
      for (int j = 0; j < n; ++ j)
      {
        // Never try to remove one that is in only one circle.
        if (_result._circle_of_node[j].size() <= 1)
        {
          continue;
        }
        for (int c: _result._circle_of_node[j])
        {
          _result.detach(j,c);
          double comp_epoi = _result.EPOI_num(base_sim,false);
          if (comp_epoi < cur_min)
          {
            cur_min = comp_epoi;
            rpair.first = j; rpair.second = c; 
          }
          _result.undo(); 
        }
      }
      _result.detach(rpair.first,rpair.second);
      _result_epoi = cur_min; 
      _detached.insert(rpair); 
      _result = _result.clean_copy(); 
    }
  }

  /// @brief Optimize to make an n-cut
  /// @param n 
  void 
  Cutter::optimize(int n)
  {
    DEBUG("Preprocess");
    pre_process();
    DEBUG("initialize epois");
    initialize_circle_epois(); 
    initialize_division();
    DEBUG("Base epoi calculation");
    _base_epoi = _parent.EPOI_num(_sims);
    DEBUG(_base_epoi);
    DEBUG("Min Cut optimization");
    std::vector<std::pair<int,int>> detacment = min_cut();
    DEBUG("Cutting " << detacment.size() << " detachment");
    _result = _parent;
    _detached.clear();
    for (auto pp: detacment)
    {
      _result.detach(pp.first,pp.second);
      _detached.insert(pp);
    }
    _result = _result.clean_copy();
    DEBUG("Post cut epoi calculation");
    _result_epoi = _result.EPOI_num(_sims);
    DEBUG(_result_epoi);
  }


  std::vector<std::pair<int,int>> 
  Cutter::min_cut() const
  {
    std::unordered_map<size_t, int> LR_res;
    std::unordered_map<size_t, int> RL_res;
    bool ff = verify_connectedness();
    BOOST_ASSERT(ff);
    auto path = find_aug_path(LR_res, RL_res);
    int count = 0;
    DEBUG("first path found:" << path.size());
    while (!path.empty())
    {
      ++count;
      update_residual(LR_res,RL_res,path);
      try
      {
        /* code */
        path = find_aug_path(LR_res,RL_res);
      }
      catch(const std::exception& e)
      {
        DEBUG("Problem finding path");
        std::cerr << e.what() << '\n';
        throw e;
      }
      DEBUG("Path " << count << ", size " << path.size());
    }
    return extract_cut(LR_res,RL_res);
  }

  std::vector<int> 
  Cutter::find_aug_path(std::unordered_map<size_t,int>& LR_res, std::unordered_map<size_t,int>& RL_res) const
  {
    if (LR_res.empty() && RL_res.empty())
    {
      initialize_residual(LR_res,RL_res);
    }
    DEBUG("Res sizes " << LR_res.size() << "," << RL_res.size());
    // BOOST_ASSERT(LR_res.size() == RL_res.size());
    const int s = _S_blocks[0];
    const int t = _T_blocks[0];
    DEBUG("s:" << s << ", t:" << t);
    // Queue for blocks.
    std::queue<int> L_queue;
    L_queue.push(s);
    // Queue for bridges
    std::queue<int> R_queue ;
    std::unordered_map<int,int> RL_pred;
    std::unordered_map<int,int> LR_pred;
    LR_pred[s] = -1;
    bool L_turn = true;
    bool found = false;
    while (!found && !(L_queue.empty() && L_turn) && !(R_queue.empty() && !L_turn))
    {
      auto& Q = L_turn?L_queue:R_queue;
      auto& Q2 = L_turn?R_queue:L_queue; 
      int u = Q.front();
      Q.pop();
      // DEBUG( (L_turn?"Circle:":"Bridge:") << u);
      const auto& neighbour = L_turn?_block_bridge_edges:_bridge_block_edges;
      auto& pred = L_turn?RL_pred:LR_pred;
      const auto res = L_turn?LR_res:RL_res; 
      
      for (auto v: neighbour.at(u))
      {
        //DEBUG("Node v:" << v);
        size_t k = Key(u,v); 
        // Skip edges with zero residual capacity
        BOOST_ASSERT(res.find(k) != res.end());
        if (res.at(k) == 0)
        {
          // DEBUG("Skipped: " << u << "->" << v);
          continue;
        }
        if (pred.find(v) != pred.end())
        {
          continue;
        }
        pred[v] = u;
        // Did we take t out of the L-queue?
        if (L_turn && v == t)
        {
          found = true;
          break;
        } 
        Q2.push(v);
      }
      // if the current queue is empty, flip. 
      if (Q.empty())
      {
        L_turn = !L_turn;
      }
    }
    DEBUG("visited: " << LR_pred.size() + RL_pred.size() << " nodes ");

    /// Extract the result from the predecessor
    std::vector<int> result;
    result.clear();
    if (LR_pred.find(t) == LR_pred.end())
    {
      return result;
    }
    int u = t;
    L_turn = true;
    while (u != -1)
    {
      result.push_back(u);
      auto& pred = L_turn?LR_pred:RL_pred;
      auto pred_it = pred.find(u);
      BOOST_ASSERT(pred_it != pred.end());
      u = pred_it->second; 
      L_turn = !L_turn;
    }
    std::reverse(result.begin(), result.end());
    // Postcondition
    BOOST_ASSERT(result.front() == s);
    BOOST_ASSERT(result.back() == t);
    return result; 
  }

  void 
  Cutter::initialize_residual(std::unordered_map<size_t,int>& LR_res, std::unordered_map<size_t,int>& RL_res) const
  {
    DEBUG("Residual initialization");
    for (auto cu: _block_bridge_edges)
    {
      int c = cu.first;
      for (int u: cu.second)
      {
        LR_res[Key(c,u)] =  1;//static_cast<int>(std::round(std::max(_block_bridge_probs.at(Key(c,u)), 1.0) + 0.1));
        RL_res[Key(u,c)] =  1;//static_cast<int>(std::round(std::max(_bridge_block_probs.at(Key(u,c)), 1.0) + 0.1));
      }
    }
    DEBUG("Residual sizes " << LR_res.size() << "," << RL_res.size());
  }
    
  void
  Cutter::initialize_circle_epois()
  {
    _circle_epois.clear();
    _circle_epois.resize(_parent._circles.size(),0.0);
    for (int c = 0; c < _parent._circles.size();++c)
    {
      _circle_epois[c] = _parent.circle_epoi(c, _sims);
    }
  }

  void 
  Cutter::initialize_division()
  {
    double max_epoi = 0.0;
    int max_index = -1;
    int second_index = -1;
    double second_epoi = 0.0;
 
    for (int c = 0; c < _parent._circles.size(); ++c)
    {
      if (_circle_epois[c] >= max_epoi)
      {
        second_index = max_index;
        max_index = c;
        max_epoi = _circle_epois[c];
      }
      else if (_circle_epois[c] > second_epoi)
      {
        second_index = c;
        second_epoi = _circle_epois[c];
      }
    }
    BOOST_ASSERT(second_index >= 0);
    // Currently just s = max, t = second max
    _S_blocks.clear();
    _T_blocks.clear();
    _S_blocks.push_back(max_index);
    _T_blocks.push_back(second_index);

  }

  void 
  Cutter::update_residual(std::unordered_map<size_t, int>& LR_res, std::unordered_map<size_t,int>& RL_res,const std::vector<int>& path) const
  {
    DEBUG("Updating");
    bool L_turn = true;
    int min_cap = 2; 
    for (int i = 0; i < path.size()-1; ++i)
    {
      const auto& res = L_turn?LR_res:RL_res;
      int u = path[i];
      int v = path[i+1];
      size_t k = Key(u,v);
      BOOST_ASSERT(res.find(k) != res.end());
      double cap = res.at(k);
      if (cap < min_cap)
      {
        min_cap = cap;
      }
      L_turn = !L_turn;
    }
    DEBUG("min cap" << min_cap);
    BOOST_ASSERT(min_cap > 0);
    for (int i = 0; i < path.size()-1; ++i)
    {
      auto& res = L_turn?LR_res:RL_res;
      auto& other = L_turn?RL_res:LR_res;
      int u = path[i];
      int v = path[i+1];
      res[Key(u,v)]-= min_cap;
      res[Key(v,u)]+= min_cap;
      L_turn = !L_turn;
    }
    DEBUG("Update done");
  }

  size_t
  Cutter::Key(int u, int v) const
  {
    return _parent.Key(u,v);
  }

  void 
  Cutter::pre_process()
  {
    _max_weight = 0.0;
    int num_bridges = 0;
    int max = 1;
    for(auto c:_parent._circles)
    {
      max = std::max(max, (int)c.size());
    }
    ++max; 

    for (int c = 0; c < _parent._circles.size(); ++c)
    {
      _block_bridge_edges[c].clear();
      for (int u: _parent._circles.at(c))
      {
        BOOST_ASSERT(0 <= u);
        BOOST_ASSERT(u <= _parent._num_nodes);
        // If the node is in more than one circle, then it is a bridge
        if (_parent._circle_of_node.at(u).size() > 1)
        {
          auto& br_bl_n =_bridge_block_edges[u];
          if (br_bl_n.empty())
          {
            ++num_bridges;
          }
          br_bl_n.push_back(c);
          auto& bl_br_n = _block_bridge_edges[c];
          bl_br_n.push_back(u);
          double lwcu = max - _parent._circles.at(c).size();
          //DEBUG(lwcu);
          _max_weight = std::max(_max_weight,lwcu);
          _bridge_block_probs[Key(u,c)] = lwcu;
          _block_bridge_probs[Key(c,u)] = lwcu;         
        }
      }
    }
    DEBUG("Bridges: " << num_bridges << ", Blocks: " << _parent._circles.size());
    DEBUG("Initialized probs " << _bridge_block_probs.size());
  }

  std::vector<std::pair<int,int>> 
  Cutter::extract_cut(std::unordered_map<size_t, int>& LR_res, std::unordered_map<size_t,int>& RL_res) const
  {
    const int s = _S_blocks[0];
    const int t = _T_blocks[0];
    std::queue<int> L_queue;
    L_queue.push(s);
    std::queue<int> R_queue ;
    std::unordered_map<int,int> RL_pred;
    std::unordered_map<int,int> LR_pred;
    LR_pred[s] = -1;
    bool L_turn = true;

    while (!(L_queue.empty() && L_turn) && !(R_queue.empty() && !L_turn))
    {
      auto& Q = L_turn?L_queue:R_queue;
      auto& Q2 = L_turn?R_queue:L_queue; 
      int u = Q.front();
      Q.pop();
      const auto& neighbour = L_turn?_block_bridge_edges:_bridge_block_edges;
      auto& pred = L_turn?RL_pred:LR_pred;
      const auto res = L_turn?LR_res:RL_res; 
      BOOST_ASSERT(neighbour.find(u) != neighbour.end());
      for (auto v: neighbour.at(u))
      {
        // Skip edges with zero residual capacity
        if (res.at(Key(u,v)) == 0)
        {
          continue;
        }
        if (pred.find(v) != pred.end())
        {
          continue;
        }
        pred[v] = u;
        Q2.push(v);
      }
      // if the current queue is empty, flip. 
      if (Q.empty())
      {
        L_turn = !L_turn;
      }
    }
    // We assume t was not found.
    BOOST_ASSERT(LR_pred.find(t) == LR_pred.end());
    std::vector<std::pair<int,int>> result;
    for (auto nbr: _bridge_block_edges)
    {
      int c = nbr.first;
      bool c_found = (RL_pred.find(c) != RL_pred.end()); 
      for (int u: nbr.second)
      {
        bool u_found = (LR_pred.find(u) != LR_pred.end());
        // Is the c--> u edge part of the cut or
        // Is the u--> c edge part of the cut?
        if ( c_found  != u_found )
        {
          result.push_back(std::make_pair(c,u));
        }
      }
    }
    return result; 
  }
  bool
  Cutter::verify_connectedness() const
  {
    int s = 0;
    int num_found = 1;
    std::vector<bool> found_C(_parent._circles.size(), false);
    std::vector<bool> found_R(_parent._num_nodes, false);
    std::vector<int> Stack_C;
    std::vector<int> Stack_R;
    bool is_C = true;
    found_C[s] = true;
    Stack_C.push_back(s);
    while ((is_C && !Stack_C.empty()) || (!is_C && !Stack_R.empty()))
    {
      auto & Stack = is_C?Stack_C:Stack_R;
      auto & Stack2 = is_C?Stack_R:Stack_C;
      auto & found = is_C?found_R:found_C;
      int u = Stack.back();
      Stack.pop_back();
      const auto& neighbours = is_C?_block_bridge_edges:_bridge_block_edges;
      for (int v: neighbours.at(u))
      {
        if (found[v]) continue;
        Stack2.push_back(v);
        found[v] = true;
      }
      if (Stack.empty())
      {
        is_C = !is_C;
      }
    }
    for (int i = 0; i < _parent._circles.size();++i)
    {
      if(!found_C[i])
      {
        DEBUG(" Circle: " << i << " was not found");
        return false;
      }
    }
    DEBUG("All circles found");
    return true;
  }
}