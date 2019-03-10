/* (c) 2019 Jithin K. Sreedharan
The following program implements all the main functions in the TIMES project
*/
#include "times.hpp"

void write_dot(PNEANet G, const char *FName_t, const char *Desc) {
  // Output node IDs as numbers
  TIntStrH NIdLabelH;
  // Generate labels for random graph
  for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
    NIdLabelH.AddDat(NI.GetId(), TStr::Fmt("Node%d", NI.GetId()));
  }
  TSnap::SaveGViz(G, FName_t, Desc, NIdLabelH);
}

void print_TIntV(TIntV v, char *S) {
  printf("%s: ", S);
  for (int i = 0; i < v.Len(); i++) {
    printf("%d ", v[i].Val);
  }
  printf("\n");
}

void print_TUInt64V(TUInt64V v, char *S) {
  printf("%s: ", S);
  for (int i = 0; i < v.Len(); i++) {
    printf("%" PRIu64 " ", v[i].Val);
  }
  printf("\n");
}

int rangeRandom(int min, int max) {
  int n = max - min + 1;
  int remainder = RAND_MAX % n;
  int x;
  do {
    x = rand();
  } while (x >= RAND_MAX - remainder);
  return min + x % n;
}

PNEANet GenPrefAttachGeneral_undirected(
    const int &time_n, const float &pr_alpha, const int &vec_p_1,
    const int &vec_p_2, const float &pr_beta, const int &vec_q_1,
    const int &vec_q_2, const float &pr_delta, const float &pr_gamma,
    const bool &self_loops_allowed) {
  PNEANet G = PNEANet::New();
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_int_distribution<> rv_m_new(vec_p_1,
                                           vec_p_2);  // unif distbn is a class
  std::uniform_int_distribution<> rv_m_old(vec_q_1, vec_q_2);
  std::bernoulli_distribution rv_alpha(pr_alpha);
  std::bernoulli_distribution rv_beta(pr_beta);
  std::bernoulli_distribution rv_delta(pr_delta);
  std::bernoulli_distribution rv_gamma(pr_gamma);

  std::vector<int> target_list;
  std::vector<int> target_list_temp;
  int target;
  uint16_t m;

  G->AddNode(0);
  m = rv_m_new(gen);
  if (self_loops_allowed) {
    for (int i = 0; i < m; i++) {
      G->AddEdge(0, 0);
      target_list.push_back(0);
      target_list.push_back(0);
    }
  }
  G->AddNode(1);
  m = rv_m_new(gen);
  for (int i = 0; i < m; i++) {
    G->AddEdge(1, 0);
    G->AddEdge(0, 1);
    target_list.push_back(0);
    target_list.push_back(1);
  }
  int source = 2;
  int source_exist;

  TNEANet::TNodeI NI_temp;

  for (int time_i = 2; time_i < time_n; time_i++) {
    if (rv_alpha(gen)) {
      target_list_temp.clear();
      G->AddNode(source);
      m = rv_m_new(gen);
      for (int i = 0; i < m; i++) {
        if (rv_beta(gen)) {
          target = *select_randomly(target_list.begin(), target_list.end());
          G->AddEdge(source, target);
          G->AddEdge(target, source);
          target_list_temp.push_back(source);
          target_list_temp.push_back(target);
        } else {
          std::uniform_int_distribution<int> random_node(
              0, (source - 1));  // guaranteed unbiased
          target = random_node(gen);
          G->AddEdge(source, target);
          G->AddEdge(target, source);
          target_list_temp.push_back(source);
          target_list_temp.push_back(target);
        }
      }
      source++;
      target_list.insert(target_list.end(), target_list_temp.begin(),
                         target_list_temp.end());
    } else {
      target_list_temp.clear();
      m = rv_m_old(gen);
      for (int i = 0; i < m; i++) {
        if (rv_delta(gen)) {
          source_exist =
              *select_randomly(target_list.begin(), target_list.end());
          //* is added because it returns a pointer Iter
        } else {
          std::uniform_int_distribution<int> random_node(
              0, (source - 1));  // guaranteed unbiased
          source_exist = random_node(gen);
        }
        if (rv_gamma(gen)) {
          target = *select_randomly(target_list.begin(), target_list.end());
        } else {
          std::uniform_int_distribution<int> random_node(
              0, (source - 1));  // guaranteed unbiased
          target = random_node(gen);
        }
        G->AddEdge(source_exist, target);
        G->AddEdge(target, source_exist);
        target_list_temp.push_back(source_exist);
        target_list_temp.push_back(target);
      }
      target_list.insert(target_list.end(), target_list_temp.begin(),
                         target_list_temp.end());
    }
  }
  return G;
}

std::tuple<double, double, double, double> count_perfect_correct_pairs_DAG(
    PNEANet G_DAG, TUInt64UInt64VH rank_new, TUInt64UInt64H rank_orig,
    uint64_t G_no_nodes) {
  std::pair<float, float> eta;

  uint64_t rank_new_len = rank_new.Len();
  uint64_t len_rank_bin_i;
  uint64_t len_rank_bin_j;

  double beta_tilde = 0;
  double alpha_tilde = 0;
  double beta_ij = 0;
  double alpha_ij = 0;
  uint64_t node_v, node_u, node_u_rank;
  TUInt64V bin_i, bin_j;
  double tt = ((double)(G_no_nodes - 1) * G_no_nodes) / 2.0;
  TIntV NIdV, NIdV_temp;
  TIntIntVH NId_reachb;

  PNGraph G_tree;

  double count_perfect_pairs_nlzd = 0;
  for (PNEANet::TObj::TNodeI NI = G_DAG->BegNI(); NI < G_DAG->EndNI(); NI++) {
    G_tree = TSnap::GetBfsTree(G_DAG, NI.GetId(), true, false);
    G_tree->GetNIdV(NIdV);
    NId_reachb.AddDat(NI.GetId(), NIdV);
    count_perfect_pairs_nlzd += (NIdV.Len() - 1) / tt;
    NIdV.Clr();
  }
  // printf("No. of true perfect pairs: %f\n", count_perfect_pairs_nlzd);
  for (uint64_t i = 0; i < (rank_new_len - 1); i++) {
    bin_i = TUInt64V(rank_new.GetDat(i));
    len_rank_bin_i = bin_i.Len();
    for (uint64_t j = (i + 1); j < rank_new_len; j++) {
      bin_j = TUInt64V(rank_new.GetDat(j));
      len_rank_bin_j = bin_j.Len();
      beta_ij = 0;
      alpha_ij = 0;
      for (uint64_t u_i = 0; u_i < len_rank_bin_i; u_i++) {
        node_u = bin_i[u_i].Val;
        node_u_rank = rank_orig.GetDat(node_u);
        for (uint64_t v_i = 0; v_i < len_rank_bin_j; v_i++) {
          node_v = bin_j[v_i].Val;
          if (node_u_rank < rank_orig.GetDat(node_v)) {
            beta_ij += 1.0;
            // TIntV ttt = TIntV(NId_reachb.GetDat(node_v));
            // bool checkk = ttt.IsIn(node_u);
            if ((NId_reachb.GetDat(node_v)).IsIn(node_u)) {
              alpha_ij += 1.0;
            }
          }
        }
      }
      beta_tilde += beta_ij / tt;
      alpha_tilde += alpha_ij / tt;
    }
  }
  // alpha_tilde is proportion of no. of perfect pairs that are correct
  // beta_tilde is proportion of no. of correct pairs
  double temp_theta_rperf = alpha_tilde / count_perfect_pairs_nlzd;
  uint64_t no_pairs = no_pairs_DAG(rank_new);
  double temp_theta_rpeel = beta_tilde * tt / no_pairs;
  return std::make_tuple(alpha_tilde, temp_theta_rperf, beta_tilde,
                         temp_theta_rpeel);
}

/* Given a DAG, find the average and maximum size of the tree starting from a
 * node in the youngest bin (top bin). The size of the tree is defined as the
 * number of nodes in the tree. */
std::pair<double, double> count_tree_size_DAG(PNEANet G_DAG,
                                              TUInt64UInt64VH rank_new,
                                              uint64_t G_no_nodes) {
  double count_tree_len = 0;
  double count_tree_len_tt;
  uint64_t rank_new_len = rank_new.Len();
  uint64_t node_u;
  TIntV NIdV;
  PNGraph G_tree;

  TUInt64V bin_l = TUInt64V(rank_new.GetDat(rank_new_len - 1));
  uint64_t len_rank_bin_l = bin_l.Len();
  double max_s = 1;
  for (uint64_t u_i = 0; u_i < len_rank_bin_l; u_i++) {
    node_u = bin_l[u_i].Val;
    G_tree = TSnap::GetBfsTree(G_DAG, node_u, true, false);
    G_tree->GetNIdV(NIdV);
    count_tree_len_tt = (NIdV.Len() - 1);
    count_tree_len += count_tree_len_tt;
    if (count_tree_len_tt > max_s) {
      max_s = count_tree_len_tt;
    }
    NIdV.Clr();
  }
  count_tree_len /= len_rank_bin_l;
  return std::make_pair(count_tree_len, max_s);
}

/* Given a DAG and a node, find the length of the longest path to any other
 * node. */
uint64_t longest_distance(const PNEANet &G_DAG, const uint64_t &G_no_nodes,
                          uint64_t node_u) {
  uint64_t node_v = 0;
  uint64_t len_long_path = 0, len_long_path_t = 0;
  uint64_t max;
  TNEANet::TNodeI NI_temp;
  if (node_u == 0) {
    len_long_path = 1;
    return len_long_path;
  } else {
    NI_temp = G_DAG->GetNI(node_u);
    max = 0;
    for (int i = 0; i < NI_temp.GetOutDeg(); i++) {
      node_v = (uint64_t)NI_temp.GetOutNId(i);
      len_long_path_t = longest_distance(G_DAG, G_no_nodes, node_v) + 1;
      if (max < len_long_path_t) max = len_long_path_t;
    }
    len_long_path += max;
    return len_long_path;
  }
}

/* Given a DAG, find average and maximum height of a tree starting from any node
 * in the youngest bin (top bin). The height is defined as the length of the
 * longest path between a node in the youngest bin to any other reachable node
 */
std::pair<double, double> count_tree_height_DAG(PNEANet G_DAG,
                                                TUInt64UInt64VH rank_new,
                                                uint64_t G_no_nodes) {
  uint64_t rank_new_len = rank_new.Len();
  uint64_t node_u;
  TIntV NIdV;
  PNGraph G_tree;
  TUInt64V bin_l = TUInt64V(rank_new.GetDat(rank_new_len - 1));
  uint64_t len_rank_bin_l = bin_l.Len();
  double max_l = 0, avg_l = 0;
  int TreeHeight_t;

  for (uint64_t u_i = 0; u_i < (len_rank_bin_l); u_i++) {
    node_u = bin_l[u_i].Val;
    TreeHeight_t = longest_distance(G_DAG, G_no_nodes, node_u) - 1;
    if (TreeHeight_t > max_l) {
      max_l = (double)TreeHeight_t;
    }
    avg_l += (double)TreeHeight_t;
  }
  avg_l /= len_rank_bin_l;
  return std::make_pair(avg_l, max_l);
}

uint64_t no_pairs_DAG(TUInt64UInt64VH rank) {
  uint64_t rank_len = rank.Len();
  std::vector<uint64_t> len_bin(rank_len);
  for (uint64_t i = 0; i < (rank_len); i++) {
    len_bin[i] = TUInt64V(rank.GetDat(i)).Len();
  }
  uint64_t no_pairs = 0;
  uint64_t len_bin_t = 0;
  for (unsigned i = len_bin.size(); i-- > 1;) {
    len_bin_t = len_bin[i];
    for (unsigned j = i; j-- > 0;) {
      no_pairs += (len_bin[j] * len_bin_t);
    }
  }
  return no_pairs;
}

float FindThetaVH(TUInt64UInt64VH rank_new, TUInt64UInt64H rank_orig,
                  uint64_t G_no_nodes) {
  std::pair<float, float> eta;

  uint64_t rank_new_len = rank_new.Len();
  uint64_t len_rank_bin_i;
  uint64_t len_rank_bin_j;

  double beta_tilde = 0;
  double beta_ij = 0;
  uint64_t node_v, node_u, node_u_rank;
  TUInt64V bin_i, bin_j;
  double tt = ((float)(rank_new_len - 1) * rank_new_len) / 2.0;
  // double tt = ((float)(G_no_nodes-1)*G_no_nodes)/2.0;
  for (uint64_t i = 0; i < (rank_new_len - 1); i++) {
    bin_i = TUInt64V(rank_new.GetDat(i));
    len_rank_bin_i = bin_i.Len();
    for (uint64_t j = (i + 1); j < rank_new_len; j++) {
      bin_j = TUInt64V(rank_new.GetDat(j));
      len_rank_bin_j = bin_j.Len();
      beta_ij = 0;
      for (uint64_t u_i = 0; u_i < len_rank_bin_i; u_i++) {
        node_u = bin_i[u_i].Val;
        node_u_rank = rank_orig.GetDat(node_u);
        for (uint64_t v_i = 0; v_i < len_rank_bin_j; v_i++) {
          node_v = bin_j[v_i].Val;
          if (node_u_rank < rank_orig.GetDat(node_v)) {
            beta_ij += 1.0 / ((float)len_rank_bin_i * len_rank_bin_j);
            // beta_ij += 1.0;
          }
        }
      }
      beta_tilde += beta_ij / tt;
    }
  }
  eta.second = beta_tilde;
  return eta.second;
}

float FindThetaH_kendaul_tau(TUInt64UInt64H rank_new, TUInt64UInt64H rank_orig,
                             uint64_t G_no_nodes, const float p) {
  double beta_tilde = 0;
  double eta = 0;

  uint64_t rank_new_len = rank_new.Len();
  uint64_t node_u, node_v, node_u_rank, node_v_rank;
  float tt = (G_no_nodes * (G_no_nodes - 1) / (float)2);
  for (uint64_t i = 0; i < (rank_new_len - 1); i++) {
    node_u = rank_new.GetDat(i);
    node_u_rank = rank_orig.GetDat(node_u);
    beta_tilde = 0;
    for (uint64_t j = (i + 1); j < rank_new_len; j++) {
      node_v = rank_new.GetDat(j);
      node_v_rank = rank_orig.GetDat(node_v);
      if (node_u_rank < node_v_rank) {
        beta_tilde += 1;
        // printf("beta_tilde: %g \n", beta_tilde);
      } else if (node_u_rank == node_v_rank) {
        beta_tilde += p;
      }
    }
    eta += beta_tilde / tt;
    // printf("eta: %g \n", eta);
  }
  // eta = beta_tilde /(G_no_nodes* (G_no_nodes-1)/(float)2);
  return eta;
}

float FindThetaVH_kendaul_tau(TUInt64UInt64VH rank_new,
                              TUInt64UInt64H rank_orig, uint64_t G_no_nodes,
                              const float p) {
  double eta;
  uint64_t rank_new_len = rank_new.Len();
  uint64_t len_rank_bin_i;
  uint64_t len_rank_bin_j;
  double beta_tilde_1 = 0, beta_tilde_2 = 0, beta_tilde_3 = 0;
  double beta_ij = 0;
  uint64_t node_v, node_u, node_u_rank, node_v_rank, node_uu, node_uu_rank;
  TUInt64V bin_i, bin_j;
  float tt = (G_no_nodes * (G_no_nodes - 1) / (float)2);
  for (uint64_t i = 0; i < (rank_new_len - 1); i++) {
    bin_i = TUInt64V(rank_new.GetDat(i));
    len_rank_bin_i = bin_i.Len();
    for (uint64_t u_i = 0; u_i < (len_rank_bin_i - 1); u_i++) {
      node_u = bin_i[u_i].Val;
      node_u_rank = rank_orig.GetDat(node_u);
      beta_ij = 0;
      for (uint64_t u_ii = u_i + 1; u_ii < len_rank_bin_i; u_ii++) {
        node_uu = bin_i[u_ii].Val;
        node_uu_rank = rank_orig.GetDat(node_uu);
        if (node_u_rank == node_uu_rank) {
          beta_ij += 1;
        }
      }
      beta_tilde_1 += beta_ij / tt;
    }
    for (uint64_t j = (i + 1); j < rank_new_len; j++) {
      bin_j = TUInt64V(rank_new.GetDat(j));
      len_rank_bin_j = bin_j.Len();
      beta_ij = 0;
      for (uint64_t u_i = 0; u_i < len_rank_bin_i; u_i++) {
        node_u = bin_i[u_i].Val;
        node_u_rank = rank_orig.GetDat(node_u);
        for (uint64_t v_i = 0; v_i < len_rank_bin_j; v_i++) {
          node_v = bin_j[v_i].Val;
          node_v_rank = rank_orig.GetDat(node_v);
          if (node_u_rank < node_v_rank) {
            beta_ij++;
          } else if (node_u_rank == node_v_rank) {
            beta_ij += p;
          }
        }
      }
      beta_tilde_2 += beta_ij / tt;
    }
  }

  bin_i = TUInt64V(rank_new.GetDat((rank_new_len - 1)));
  len_rank_bin_i = bin_i.Len();
  for (uint64_t u_i = 0; u_i < (len_rank_bin_i - 1); u_i++) {
    node_u = bin_i[u_i].Val;
    node_u_rank = rank_orig.GetDat(node_u);
    beta_ij = 0;
    for (uint64_t u_ii = u_i + 1; u_ii < len_rank_bin_i; u_ii++) {
      node_uu = bin_i[u_ii].Val;
      node_uu_rank = rank_orig.GetDat(node_uu);
      if (node_u_rank == node_uu_rank) {
        beta_ij += 1;
      }
    }
    beta_tilde_3 += beta_ij / tt;
  }
  eta = beta_tilde_1 + beta_tilde_2 + beta_tilde_3;
  return eta;
}

// This calculates the rho and theta of real networks when there are very few
// bins in the original ranking  Meaning the original ranking is not at all a
// total order.  Uses Peeling algorithm
// BECAREFUL WHEN DEALING WITH LARGE GRAPHS; OVERFLOW MAY OCCUR
std::pair<uint64_t, uint64_t> find_prec_recall_real_networks(
    TUInt64UInt64VH rank_parallel, TUInt64UInt64H rank_orig,
    uint64_t G_no_nodes) {
  uint64_t rho = 0;
  uint64_t rank_new_len = rank_parallel.Len();
  uint64_t len_rank_bin_i;
  uint64_t len_rank_bin_j;
  uint64_t rho_temp = 0;
  uint64_t node_v, node_u, node_u_rank, node_v_rank;
  TUInt64V bin_i, bin_j;
  uint64_t count_K = 0;

  for (int64_t i = (rank_new_len - 1); i >= 1; --i) {
    bin_i = TUInt64V(rank_parallel.GetDat(i));
    len_rank_bin_i = bin_i.Len();
    for (int64_t j = (i - 1); j >= 0; --j) {
      bin_j = TUInt64V(rank_parallel.GetDat(j));
      len_rank_bin_j = bin_j.Len();
      rho_temp = 0;
      for (uint64_t u_i = 0; u_i < len_rank_bin_i; ++u_i) {
        node_u = bin_i[u_i].Val;
        node_u_rank = rank_orig.GetDat(node_u);
        for (uint64_t v_i = 0; v_i < len_rank_bin_j; ++v_i) {
          node_v = bin_j[v_i].Val;
          node_v_rank = rank_orig.GetDat(node_v);
          if (node_u_rank != node_v_rank) {
            count_K += 1;
            if (node_u_rank > node_v_rank) rho_temp += 1.0;
          }
        }
      }
      rho += rho_temp;
    }
  }
  return std::make_pair(rho, count_K);
}

std::pair<uint64_t, uint64_t> find_prec_recall_real_networks_STL(
    Tuint64uint64VM rank_parallel, Tuint64uint64M rank_orig,
    const uint64_t &G_no_nodes, const int &no_threads) {
  uint64_t rho = 0;
  uint64_t rank_new_len = rank_parallel.size();
  uint64_t count_K = 0;

  omp_set_num_threads(no_threads);
#pragma omp parallel
  {
    std::vector<uint64_t> bin_i, bin_j;
    uint64_t node_v, node_u, node_u_rank, node_v_rank;
    uint64_t len_rank_bin_i;
    uint64_t len_rank_bin_j;
    uint64_t i, j, u_i, v_i;
    for (i = (rank_new_len - 1); i >= 1; --i) {
      bin_i = rank_parallel.at(i);
      len_rank_bin_i = bin_i.size();
      for (j = i; j-- > 0;) {  // Loop from i-1 to 0
        bin_j = rank_parallel.at(j);
        len_rank_bin_j = bin_j.size();
#pragma omp for reduction(+ : count_K, rho)
        for (u_i = 0; u_i < len_rank_bin_i; u_i++) {
          node_u = bin_i[u_i];
          node_u_rank = rank_orig.at(node_u);
          for (v_i = 0; v_i < len_rank_bin_j; v_i++) {
            node_v = bin_j[v_i];
            node_v_rank = rank_orig.at(node_v);
            // MAKES SURE THAT WE COUNT THE PAIRS IN TO RHO ONLY WHEN THEIR
            // ORIGINAL RANKING DIFFER
            if (node_u_rank != node_v_rank) {
              count_K += 1;
              if (node_u_rank > node_v_rank) rho += 1;
            }
          }
        }
      }
    }
  }
  return std::make_pair(rho, count_K);
}

uint64_t find_max_neighb_level(PNEANet G, const uint64_t &node,
                               TUInt64UInt64H rank_node) {
  const PNEANet::TObj::TNodeI NodeI = G->GetNI(node);
  uint64_t max_level = 0;
  uint64_t temp;
  uint64_t DstNId;
  for (int v = 0; v < NodeI.GetDeg(); v++) {
    // DstNId = NodeI.GetOutNId(v);
    DstNId = NodeI.GetNbrNId(v);
    temp = rank_node[DstNId].Val;
    // printf("temp:%d\n", temp);
    if (max_level < temp) max_level = temp;
  }
  return max_level;
}

double find_avg_neighb_level(PNEANet G, const uint64_t &node,
                             TUInt64UInt64H &rank_node) {
  const PNEANet::TObj::TNodeI NodeI = G->GetNI(node);
  double temp = 0;
  uint64_t DstNId;
  // uint64_t deg = NodeI.GetODeg();
  double deg = NodeI.GetOutDeg();
  for (uint v = 0; v < deg; v++) {
    // DstNId = NodeI.GetNbrNId(v);
    DstNId = NodeI.GetOutNId(v);
    // temp += (float)rank_node[DstNId].Val;
    double temp11 = (float)rank_node.GetDat(DstNId);
    temp += temp11;
    // printf("node_u: %d, DstNId:%d, rank_node[DstNId]: %f\n",node,
    // DstNId,temp11);
  }
  return (temp / deg);
}

// Same bin pairs
std::pair<double, uint64_t> Count_extra_pairs(PNEANet G, uint64_t G_no_nodes,
                                              TUInt64UInt64H rank_orig,
                                              TUInt64UInt64VH rank_new,
                                              TUInt64UInt64H rank_node) {
  uint64_t rank_new_len = rank_new.Len();
  uint64_t len_rank_bin_i;
  double eta = 0;
  double beta_ij = 0;
  uint64_t node_u, node_u_rank, node_uu, node_uu_rank;
  TUInt64V bin_i;
  double tt = (G_no_nodes * (G_no_nodes - 1) / (double)2);
  uint64_t count_extra_pairs = 0;
  // uint64_t node_u_level, node_uu_level;
  double node_u_level, node_uu_level;
  // for(uint64_t i = (rank_new_len-1); i >= std::ceil(rank_new_len/2.0); --i)
  // for (uint64_t i = (rank_new_len - 1); i >= (rank_new_len - 1); --i) {
  uint64_t no_orig_pairs_same_bin = 0;
  for (uint64_t i = (rank_new_len - 1); i >= 1; --i) {
    bin_i = TUInt64V(rank_new.GetDat(i));
    len_rank_bin_i = bin_i.Len();
    for (uint64_t u_i = 0; u_i < (len_rank_bin_i - 1); u_i++) {
      node_u = bin_i[u_i].Val;
      node_u_rank = rank_orig.GetDat(node_u);
      beta_ij = 0;
      node_u_level = find_avg_neighb_level(G, node_u, rank_node);
      for (uint64_t u_ii = u_i + 1; u_ii < len_rank_bin_i; u_ii++) {
        node_uu = bin_i[u_ii].Val;
        node_uu_rank = rank_orig.GetDat(node_uu);
        node_uu_level = find_avg_neighb_level(G, node_uu, rank_node);
        if (node_uu_rank != node_u_rank) {
          no_orig_pairs_same_bin += 1;
          if (node_uu_level != node_u_level) {
            //  We say if (node_uu_level > node_u_level) node_uu is younger;
            count_extra_pairs += 1;
            if (((node_uu_level > node_u_level) &&
                 (node_uu_rank > node_u_rank)) ||
                ((node_uu_level < node_u_level) &&
                 (node_uu_rank < node_u_rank))) {
              beta_ij += 1.0;
            }
          }
        }
      }
      eta += beta_ij / tt;
    }
  }
  return std::make_pair(eta, count_extra_pairs);
}

// Same bin pairs
std::tuple<uint64_t, uint64_t, uint64_t> Count_extra_pairs_real_networks(
    PNEANet G, uint64_t G_no_nodes, TUInt64UInt64H rank_orig,
    TUInt64UInt64VH rank_new, TUInt64UInt64H rank_node) {
  uint64_t rank_new_len = rank_new.Len();
  uint64_t len_rank_bin_i;
  uint64_t eta = 0;
  uint64_t beta_ij = 0;
  uint64_t node_u, node_u_rank, node_uu, node_uu_rank;
  TUInt64V bin_i;
  uint64_t count_extra_pairs = 0;
  double node_u_level, node_uu_level;
  uint64_t no_orig_pairs_same_bin = 0;
  for (uint64_t i = (rank_new_len - 1); i >= 1; --i) {
    bin_i = TUInt64V(rank_new.GetDat(i));
    len_rank_bin_i = bin_i.Len();
    for (uint64_t u_i = 0; u_i < (len_rank_bin_i - 1); u_i++) {
      node_u = bin_i[u_i].Val;
      node_u_rank = rank_orig.GetDat(node_u);
      beta_ij = 0;
      node_u_level = find_avg_neighb_level(G, node_u, rank_node);
      for (uint64_t u_ii = u_i + 1; u_ii < len_rank_bin_i; u_ii++) {
        node_uu = bin_i[u_ii].Val;
        node_uu_rank = rank_orig.GetDat(node_uu);
        node_uu_level = find_avg_neighb_level(G, node_uu, rank_node);
        if (node_uu_rank != node_u_rank) {
          no_orig_pairs_same_bin += 1;
          if (node_uu_level != node_u_level) {
            //  We say if (node_uu_level > node_u_level) node_uu is younger;
            count_extra_pairs += 1;
            if (((node_uu_level > node_u_level) &&
                 (node_uu_rank > node_u_rank)) ||
                ((node_uu_level < node_u_level) &&
                 (node_uu_rank < node_u_rank))) {
              beta_ij += 1;
            }
          }
        }
      }
      eta += beta_ij;
    }
  }
  return std::make_tuple(eta, count_extra_pairs, no_orig_pairs_same_bin);
}

double find_avg_neighb_level_STL(const PNEANet::TObj::TNodeI NodeI,
                                 const Tuint64uint64M &rank_node,
                                 const int &no_threads) {
  uint64_t temp = 0;
  uint64_t DstNId;
  uint64_t deg = NodeI.GetOutDeg();
  for (uint v = 0; v < deg; v++) {
    DstNId = NodeI.GetOutNId(v);
    temp += rank_node.at(DstNId);
  }
  return (temp / (double)deg);
}

// Same bin pairs
std::tuple<uint64_t, uint64_t, uint64_t> Count_extra_pairs_real_networks_STL(
    PNEANet G, const uint64_t &G_no_nodes, const Tuint64uint64M &rank_orig,
    const Tuint64uint64VM &rank_new, const Tuint64uint64M &rank_node,
    const int &no_threads) {
  uint64_t rank_new_len = rank_new.size();
  uint64_t eta = 0;
  uint64_t count_extra_pairs = 0;
  uint64_t no_orig_pairs_same_bin = 0;
  omp_set_num_threads(no_threads);
#pragma omp parallel
  {
    std::vector<uint64_t> bin_i;
    uint64_t len_rank_bin_i;
    uint64_t node_u, node_u_rank, node_uu, node_uu_rank;
    uint64_t u_i, u_ii;
    for (uint64_t i = (rank_new_len - 1); i >= 1; --i) {
      bin_i = rank_new.at(i);
      len_rank_bin_i = bin_i.size();
#pragma omp for reduction(+ : eta, no_orig_pairs_same_bin, count_extra_pairs)
      for (u_i = 0; u_i < (len_rank_bin_i - 1); u_i++) {
        node_u = bin_i[u_i];
        node_u_rank = rank_orig.at(node_u);
        const PNEANet::TObj::TNodeI NodeI = G->GetNI(node_u);
        double node_u_level =
            find_avg_neighb_level_STL(NodeI, rank_node, no_threads);
        for (u_ii = u_i + 1; u_ii < len_rank_bin_i; u_ii++) {
          node_uu = bin_i[u_ii];
          node_uu_rank = rank_orig.at(node_uu);
          const PNEANet::TObj::TNodeI NodeII = G->GetNI(node_uu);
          double node_uu_level =
              find_avg_neighb_level_STL(NodeII, rank_node, no_threads);
          if (node_uu_rank != node_u_rank) {
            no_orig_pairs_same_bin += 1;
            if (node_uu_level != node_u_level) {
              count_extra_pairs += 1;
              if (((node_uu_level > node_u_level) &&
                   (node_uu_rank > node_u_rank)) ||
                  ((node_uu_level < node_u_level) &&
                   (node_uu_rank < node_u_rank))) {
                eta += 1;
              }
            }
          }
        }
      }
    }
  }
  return std::make_tuple(eta, count_extra_pairs, no_orig_pairs_same_bin);
}

// Bin_0: Oldest nodes's bin, Bin_k: Newst nodes's bin, k being the depth of DAG
TUInt64UInt64VH find_parallel_rank(PNEANet G) {
  uint64_t no_nodes = G->GetNodes();

  PNEANet G_n = TNEANet::New();
  *G_n = *G;

  TUInt64UInt64VH rank_new;
  TUInt64UInt64VH rank_new_1;

  TUInt64V MnDegV;
  uint64_t MnDeg;
  uint64_t depth_dag;
  uint64_t i = 0;

  while (((uint64_t)G_n->GetNodes()) > 0) {
    MnDeg = 100 * no_nodes;

    MnDegV.Clr();
    for (PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
      if (MnDeg > (uint64_t)NI.GetInDeg()) {
        MnDegV.Clr();
        MnDeg = (uint64_t)NI.GetInDeg();
      }
      if (MnDeg == (uint64_t)NI.GetInDeg()) {
        MnDegV.Add((uint64_t)NI.GetId());
      }
    }
    rank_new.AddDat(i, MnDegV);

    for (uint64_t k = 0; k < (uint64_t)MnDegV.Len(); k++) {
      G_n->DelNode(MnDegV[k].Val);
    }
    i++;
  }
  depth_dag = i;
  for (uint64_t k = 0; k < depth_dag; k++) {
    rank_new_1.AddDat((depth_dag - 1 - k), rank_new.GetDat(k));
  }
  return rank_new_1;
}

// Bin_0: Oldest nodes's bin, Bin_k: Newst nodes's bin, k being the depth of DAG
// Returns node rank too, rank 0 for the oldest node
std::pair<TUInt64UInt64VH, TUInt64UInt64H> find_parallel_rank_returnnoderank(
    PNEANet G) {
  uint64_t no_nodes = G->GetNodes();
  PNEANet G_n = TNEANet::New();
  *G_n = *G;
  TUInt64UInt64VH rank_new;
  TUInt64UInt64VH rank_new_1;
  TUInt64UInt64H rank_node;
  TUInt64V MnDegV;
  TUInt64V temp_vec;
  uint64_t MnDeg;
  uint64_t depth_dag;
  uint64_t i = 0;

  while (((uint64_t)G_n->GetNodes()) > 0) {
    MnDeg = 100 * no_nodes;

    MnDegV.Clr();
    for (PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
      if (MnDeg > (uint64_t)NI.GetInDeg()) {
        MnDegV.Clr();
        MnDeg = (uint64_t)NI.GetInDeg();
      }
      if (MnDeg == (uint64_t)NI.GetInDeg()) {
        MnDegV.Add((uint64_t)NI.GetId());
      }
    }
    rank_new.AddDat(i, MnDegV);
    for (uint64_t k = 0; k < (uint64_t)MnDegV.Len(); k++) {
      G_n->DelNode(MnDegV[k].Val);
    }
    i++;
  }
  depth_dag = i;
  for (uint64_t k = 0; k < depth_dag; k++) {
    temp_vec = TUInt64V(rank_new.GetDat(k));
    // print_TUInt64V(temp_vec, "temp_vec");
    rank_new_1.AddDat((depth_dag - 1 - k), temp_vec);
    for (uint64_t kk = 0, kk_max = temp_vec.Len(); kk < kk_max; kk++) {
      rank_node.AddDat(temp_vec[kk].Val, (depth_dag - 1 - k));
      // int hhh = temp_vec[kk].Val;
      // int kkk = rank_node.GetDat(hhh);
      // printf("node: %d, rank: %d\n",hhh,kkk);
      // cout<<temp_vec[kk].Val<<","<<rank_node.GetDat(hhh)<<endl;
    }
  }
  return std::make_pair(rank_new_1, rank_node);
}

std::pair<Tuint64uint64VM, Tuint64uint64M>
find_parallel_rank_returnnoderank_STL(PNEANet G) {
  uint64_t no_nodes = G->GetNodes();
  PNEANet G_n = TNEANet::New();
  *G_n = *G;
  Tuint64uint64VM rank_new;
  Tuint64uint64VM rank_new_1;
  Tuint64uint64M rank_node;
  std::vector<uint64_t> MnDegV;
  std::vector<uint64_t> temp_vec;
  uint64_t MnDeg;
  uint64_t depth_dag;
  uint64_t i = 0;

  while (((uint64_t)G_n->GetNodes()) > 0) {
    MnDeg = 100 * no_nodes;

    MnDegV.clear();
    for (PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
      if (MnDeg > (uint64_t)NI.GetInDeg()) {
        MnDegV.clear();
        MnDeg = (uint64_t)NI.GetInDeg();
      }
      if (MnDeg == (uint64_t)NI.GetInDeg()) {
        MnDegV.push_back((uint64_t)NI.GetId());
      }
    }
    rank_new.insert({i, MnDegV});
    for (uint64_t k = 0; k < (uint64_t)MnDegV.size(); k++) {
      G_n->DelNode(MnDegV.at(k));
    }
    i++;
  }

  depth_dag = i;
  for (uint64_t k = 0; k < depth_dag; k++) {
    temp_vec = rank_new.at(k);
    rank_new_1.insert({(depth_dag - 1 - k), temp_vec});
    for (uint64_t kk = 0, kk_max = temp_vec.size(); kk < kk_max; kk++) {
      rank_node.insert({temp_vec.at(kk), (depth_dag - 1 - k)});
    }
  }
  return std::make_pair(rank_new_1, rank_node);
}

std::pair<TUInt64UInt64VH, TUInt64UInt64H>
find_parallel_rank_returnnoderank_PUN(PUNGraph G) {
  uint64_t no_nodes = G->GetNodes();
  PUNGraph G_n = TUNGraph::New();
  *G_n = *G;
  TUInt64UInt64VH rank_new;
  TUInt64UInt64VH rank_new_1;
  TUInt64UInt64H rank_node;
  TUInt64V MnDegV;
  TUInt64V temp_vec;
  uint64_t MnDeg;
  uint64_t depth_dag;
  uint64_t i = 0;

  while (((uint64_t)G_n->GetNodes()) > 0) {
    MnDeg = 100 * no_nodes;

    MnDegV.Clr();
    for (PUNGraph::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
      if (MnDeg > (uint64_t)NI.GetInDeg()) {
        MnDegV.Clr();
        MnDeg = (uint64_t)NI.GetInDeg();
      }
      if (MnDeg == (uint64_t)NI.GetInDeg()) {
        MnDegV.Add((uint64_t)NI.GetId());
      }
    }
    rank_new.AddDat(i, MnDegV);
    for (uint64_t k = 0; k < (uint64_t)MnDegV.Len(); k++) {
      G_n->DelNode(MnDegV[k].Val);
    }
    i++;
  }
  depth_dag = i;
  for (uint64_t k = 0; k < depth_dag; k++) {
    temp_vec = TUInt64V(rank_new.GetDat(k));
    rank_new_1.AddDat((depth_dag - 1 - k), temp_vec);
    for (uint64_t kk = 0, kk_max = temp_vec.Len(); kk < kk_max; kk++) {
      rank_node.AddDat(temp_vec[kk].Val, (depth_dag - 1 - k));
    }
  }
  return std::make_pair(rank_new_1, rank_node);
}

// Bin_0: Oldest nodes's bin, Bin_k: Newst nodes's bin, k being the depth of DAG
TUInt64UInt64VH find_parallel_rank_unif_sampling(PNEANet G, const int &vec_p_1,
                                                 const int &vec_p_2) {
  uint64_t no_nodes = G->GetNodes();

  PNEANet G_n = TNEANet::New();
  ;
  *G_n = *G;

  TUInt64UInt64VH rank_new;
  TUInt64UInt64VH rank_new_1;

  TUInt64V MnDegV;
  uint64_t MnDeg;
  uint64_t depth_dag;
  uint64_t i = 0;

  while (((uint64_t)G_n->GetNodes()) > 0) {
    MnDeg = 100 * no_nodes;

    MnDegV.Clr();
    for (PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
      if (MnDeg > (uint64_t)NI.GetInDeg()) {
        MnDegV.Clr();
        MnDeg = (uint64_t)NI.GetInDeg();
      }
      if (MnDeg == (uint64_t)NI.GetInDeg()) {
        MnDegV.Add((uint64_t)NI.GetId());
      }
    }
    rank_new.AddDat(i, MnDegV);

    for (uint64_t k = 0; k < (uint64_t)MnDegV.Len(); k++) {
      G_n->DelNode(MnDegV[k].Val);
    }
    i++;
  }
  depth_dag = i;
  for (uint64_t k = 0; k < depth_dag; k++) {
    rank_new_1.AddDat((depth_dag - 1 - k), rank_new.GetDat(k));
  }
  return rank_new_1;
}

PNEANet create_DAG(PNEANet G, TUInt64UInt64VH &rank_parallel) {
  PNEANet G_DAG = TNEANet::New();
  *G_DAG = *G;
  TUInt64V bin_i, bin_j;
  uint64_t node_v, node_u;
  uint64_t rank_new_len_1 = rank_parallel.Len();
  uint64_t len_rank_bin_i;
  uint64_t len_rank_bin_j;
  for (uint64_t i = 0; i < (rank_new_len_1 - 1); i++) {
    bin_i = TUInt64V(rank_parallel.GetDat(i));
    // bin_0: oldest bin
    len_rank_bin_i = bin_i.Len();
    for (uint64_t j = (i + 1); j < rank_new_len_1; j++) {
      bin_j = TUInt64V(rank_parallel.GetDat(j));
      len_rank_bin_j = bin_j.Len();
      for (uint64_t u_i = 0; u_i < len_rank_bin_i; u_i++) {
        node_u = bin_i[u_i].Val;
        for (uint64_t v_i = 0; v_i < len_rank_bin_j; v_i++) {
          node_v = bin_j[v_i].Val;
          G_DAG->DelEdge(node_u, node_v);
        }
      }
    }
  }
  return G_DAG;
}

std::pair<PNEANet, TUInt64UInt64VH> find_parallel_rank_DAG(PNEANet G) {
  uint64_t no_nodes = G->GetNodes();

  PNEANet G_n = TNEANet::New();
  ;
  *G_n = *G;
  TUInt64UInt64VH rank_new;
  TUInt64UInt64VH rank_new_1;
  TUInt64V MnDegV;
  uint64_t MnDeg;
  uint64_t depth_dag;
  uint64_t i = 0;

  while (((uint64_t)G_n->GetNodes()) > 0) {
    MnDeg = 100 * no_nodes;

    MnDegV.Clr();
    for (PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
      if (MnDeg > (uint64_t)NI.GetInDeg()) {
        MnDegV.Clr();
        MnDeg = (uint64_t)NI.GetInDeg();
      }
      if (MnDeg == (uint64_t)NI.GetInDeg()) {
        MnDegV.Add((uint64_t)NI.GetId());
      }
    }
    rank_new.AddDat(i, MnDegV);

    for (uint64_t k = 0; k < (uint64_t)MnDegV.Len(); k++) {
      G_n->DelNode(MnDegV[k].Val);
    }
    i++;
  }
  depth_dag = i;
  for (uint64_t k = 0; k < depth_dag; k++) {
    rank_new_1.AddDat((depth_dag - 1 - k), rank_new.GetDat(k));
  }

  PNEANet G_DAG = TNEANet::New();
  ;
  *G_DAG = *G;
  TUInt64V bin_i, bin_j;
  uint64_t node_v, node_u;
  uint64_t rank_new_len_1 = rank_new_1.Len();
  uint64_t len_rank_bin_i;
  uint64_t len_rank_bin_j;
  for (uint64_t i = 0; i < (rank_new_len_1 - 1); i++) {
    bin_i = TUInt64V(rank_new_1.GetDat(i)); // bin_0: oldest bin
    len_rank_bin_i = bin_i.Len();
    for (uint64_t j = (i + 1); j < rank_new_len_1; j++) {
      bin_j = TUInt64V(rank_new_1.GetDat(j));
      len_rank_bin_j = bin_j.Len();
      for (uint64_t u_i = 0; u_i < len_rank_bin_i; u_i++) {
        node_u = bin_i[u_i].Val;
        for (uint64_t v_i = 0; v_i < len_rank_bin_j; v_i++) {
          node_v = bin_j[v_i].Val;
          G_DAG->DelEdge(node_u, node_v);
        }
      }
    }
  }
  return std::make_pair(G_DAG, rank_new_1);
}


// MODIFICATION OF find_parallel_rank WHEN MIN AND MAX OF DEGRESS TO BE PEELED
// GIVEN
TIntIntVH find_parallel_rank_minmax_deg(PNEANet G, const int min_indeg,
                                        const int max_indeg) {
  PNEANet G_n = TNEANet::New();
  ;
  *G_n = *G;

  TIntIntVH rank_new;
  TIntIntVH rank_new_1;

  TIntV MnDegV;
  // int MnDeg;
  int depth_dag;
  int i = 0;

  while ((G_n->GetNodes()) > 0) {
    // MnDeg = 100*no_nodes;

    MnDegV.Clr();
    for (PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
      if ((NI.GetInDeg() >= min_indeg) && (NI.GetInDeg() <= max_indeg)) {
        MnDegV.Add(NI.GetId());
        // printf("check\n");
      }
    }
    rank_new.AddDat(i, MnDegV);

    for (int k = 0; k < MnDegV.Len(); k++) {
      G_n->DelNode(MnDegV[k].Val);
    }
    i++;
    // printf("%d\n", i);
  }
  printf("exited\n");
  depth_dag = i;
  for (int k = 0; k < depth_dag; k++) {
    rank_new_1.AddDat((depth_dag - 1 - k), rank_new.GetDat(k));
  }
  return rank_new_1;
}

std::vector<std::vector<double>> find_recall_plot(TUInt64UInt64H rank_node,
                                                  uint64_t width_box,
                                                  uint64_t total_boxes) {
  std::vector<std::vector<double>> recall_matrix(
      total_boxes, std::vector<double>(total_boxes));
  for (uint64_t i = 0; i < total_boxes - 1; ++i) {
    for (uint64_t j = i + 1; j < total_boxes; ++j) {
      std::vector<int> pairs(width_box);
      std::iota(std::begin(pairs), std::end(pairs),
                i * width_box);  // i * width_box is the starting number
      std::vector<int> pairs_temp(width_box);
      std::iota(std::begin(pairs_temp), std::end(pairs_temp), j * width_box);
      pairs.insert(pairs.end(), pairs_temp.begin(), pairs_temp.end());
      int temp = 0;
      for (auto box_i = pairs.begin(); box_i != pairs.end(); ++box_i) {
        int rank_box_i = rank_node.GetDat(*box_i);
        for (auto box_j = box_i; ++box_j != pairs.end(); /**/) {
          int rank_box_j = rank_node.GetDat(*box_j);
          if (((*box_i < *box_j) && (rank_box_i < rank_box_j)) ||
              ((*box_i > *box_j) && (rank_box_i > rank_box_j))) {
            temp += 1;
          }
        }
      }
      recall_matrix[i][j] =
          temp / (double)(2 * nCk(width_box, 2) + std::pow(width_box, 2));
      recall_matrix[j][i] = recall_matrix[i][j];
    }
  }
  return recall_matrix;
}

std::vector<std::vector<double>> find_recall_plot(Tuint64uint64M rank_node,
                                                  uint64_t width_box,
                                                  uint64_t total_boxes) {
  std::vector<std::vector<double>> recall_matrix(
      total_boxes, std::vector<double>(total_boxes));
  for (uint64_t i = 0; i < total_boxes - 1; ++i) {
    for (uint64_t j = i + 1; j < total_boxes; ++j) {
      std::vector<int> pairs(width_box);
      std::iota(std::begin(pairs), std::end(pairs),
                i * width_box);  // i * width_box is the starting number
      std::vector<int> pairs_temp(width_box);
      std::iota(std::begin(pairs_temp), std::end(pairs_temp), j * width_box);
      pairs.insert(pairs.end(), pairs_temp.begin(), pairs_temp.end());
      int temp = 0;
      for (auto box_i = pairs.begin(); box_i != pairs.end(); ++box_i) {
        int rank_box_i = rank_node.at(*box_i);
        for (auto box_j = box_i; ++box_j != pairs.end(); /**/) {
          int rank_box_j = rank_node.at(*box_j);
          if (((*box_i < *box_j) && (rank_box_i < rank_box_j)) ||
              ((*box_i > *box_j) && (rank_box_i > rank_box_j))) {
            temp += 1;
          }
        }
      }
      recall_matrix[i][j] =
          temp / (double)(2 * nCk(width_box, 2) + std::pow(width_box, 2));
      recall_matrix[j][i] = recall_matrix[i][j];
    }
  }
  return recall_matrix;
}

// Implements binomial coefficient
double nCk(int n, int k) {
  if (k > n) std::cerr << "k > n in nCk calculation" << std::endl;
  if (k == 0) return 1;
  if (k > n / 2) return nCk(n, n - k);
  return (n * nCk(n - 1, k - 1) / (double)k);
}
