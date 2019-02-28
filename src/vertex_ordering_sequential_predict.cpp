// This is an improvement over Vertex Ordering with precition of non-perfect
// pairs using Lemma 1 expression with Sequential technique.

#include "times.hpp"

// print vector contents:
void print_vector(std::vector<int> &temp, char *s) {
  printf("%s: ", s);
  for (auto &i : temp) {
    cout << i << " ";
  }
  cout << endl;
}

void print_vector(std::vector<std::pair<int, int>> &temp, char *s) {
  printf("%s: ", s);
  for (auto &i : temp) {
    cout << "(" << i.first << "," << i.second << ")"
         << " ";
  }
  cout << endl;
}

void print_linear_extension_from_rank_pair(
    std::vector<std::pair<int, int>> &temp, char *s) {
  printf("%s:\n", s);
  for (auto it = temp.rbegin(); it != temp.rend(); ++it) {
    if (it != temp.rbegin()) {
      std::cout << ' ';
      std::cout << " " << (*it).second;
    } else {
      std::cout << (*it).first << " " << (*it).second;
    }
  }
  cout << endl;
}
// print set contents:
void print_set(std::set<int> &temp, char *s) {
  printf("%s: ", s);
  for (auto &i : temp) {
    cout << i << " ";
  }
  cout << endl;
}

// print map<int,pair<int,int>> contents:
void print_map(std::map<int, std::pair<int, int>> &temp, char *s) {
  printf("%s:\n", s);
  for (auto &i : temp) {
    cout << i.first << ": "
         << "(" << (i.second).first << "," << (i.second).second << ")" << endl;
  }
}

// Sequential algorithm. Returns rank in pairwaise format
// RANK format rank_id: Node,Node; rank 0: (second oldest,oldest)
std::map<int, std::pair<int, int>> find_seq_rank_unif_mappair(PNEANet G) {
  uint64_t no_nodes = G->GetNodes();
  PNEANet G_n = TNEANet::New();
  *G_n = *G;
  std::map<int, std::pair<int, int>> rank_new;
  uint64_t i = no_nodes - 2;
  TUInt64V MnDegV;
  uint64_t MnDeg;
  uint64_t node_sel;
  uint64_t temp_i;
  bool flag_i = 1;

  std::random_device rd;
  std::mt19937 gen(rd());

  while (i > 0) {
    MnDeg = 100 * no_nodes;
    for (PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
      if (MnDeg > (uint64_t)NI.GetInDeg()) {
        MnDegV.Clr();
        MnDeg = (uint64_t)NI.GetInDeg();
      }
      if (MnDeg == (uint64_t)NI.GetInDeg()) {
        MnDegV.Add((uint64_t)NI.GetId());
      }
    }
    std::uniform_int_distribution<int> min_deg(
        0, (MnDegV.Len() - 1)); // guaranteed unbiased
    node_sel = MnDegV[min_deg(gen)].Val;
    // node_sel = MnDegV[TInt::Rnd.GetUniDevInt(MnDegV.Len())].Val;
    if (flag_i) {
      temp_i = node_sel;
      G_n->DelNode(node_sel);
      flag_i = 0;
      continue;
    }
    rank_new.insert(std::make_pair(i, std::make_pair(temp_i, node_sel)));
    temp_i = node_sel;
    G_n->DelNode(node_sel);
    i--;
  }
  TNEANet::TNodeI NI = G_n->BegNI();
  node_sel = NI.GetId();
  rank_new.insert(std::make_pair(i, std::make_pair(temp_i, node_sel)));
  return rank_new;
}

// Returns perfect pairs, non-perfect pairs (in distinct bins), same bin pairs
std::tuple<double, double, double, double,
           std::vector<std::pair<uint64_t, uint64_t>>,
           std::vector<std::pair<uint64_t, uint64_t>>,
           std::vector<std::pair<uint64_t, uint64_t>>>
FindRhoCountPerfectPeel(PNEANet G_DAG, TUInt64UInt64VH &rank_parallel,
                        TUInt64UInt64H &rank_orig, uint64_t &G_no_nodes) {
  uint64_t rank_parallel_len = rank_parallel.Len();
  uint64_t len_rank_bin_i;
  uint64_t len_rank_bin_j;
  double beta_tilde = 0;
  double alpha_tilde = 0;
  double beta_ij = 0;
  double alpha_ij = 0;
  uint64_t node_v, node_u, node_u_rank, node_v_rank, node_uu;
  TUInt64V bin_i, bin_j;
  double tt = ((double)(G_no_nodes - 1) * G_no_nodes) / 2.0;
  TIntV NIdV, NIdV_temp, NIdV_temp_1;
  TIntIntVH NId_reachb;
  PNGraph G_tree;
  std::vector<std::pair<uint64_t, uint64_t>> nonperfectpairs;
  std::vector<std::pair<uint64_t, uint64_t>> perfectpairs;
  std::vector<std::pair<uint64_t, uint64_t>> samebinpairs;

  double count_perfect_pairs_nlzd = 0;
  for (PNEANet::TObj::TNodeI NI = G_DAG->BegNI(); NI < G_DAG->EndNI(); NI++) {
    G_tree = TSnap::GetBfsTree(G_DAG, NI.GetId(), true, false);
    G_tree->GetNIdV(NIdV);
    NId_reachb.AddDat(NI.GetId(), NIdV);
    count_perfect_pairs_nlzd += (NIdV.Len() - 1) / tt;
    NIdV.Clr();
  }

  for (int64_t i = (rank_parallel_len - 1); i >= 1; --i) {
    bin_i = TUInt64V(rank_parallel.GetDat(i));
    len_rank_bin_i = bin_i.Len();
    for (int64_t j = (i - 1); j >= 0; --j) {
      bin_j = TUInt64V(rank_parallel.GetDat(j));
      len_rank_bin_j = bin_j.Len();
      beta_ij = 0;
      alpha_ij = 0;
      for (uint64_t u_i = 0; u_i < len_rank_bin_i; u_i++) {
        node_u = bin_i[u_i].Val;
        node_u_rank = rank_orig.GetDat(node_u);
        NIdV_temp = TIntV(NId_reachb.GetDat(node_u));
        for (uint64_t v_i = 0; v_i < len_rank_bin_j; v_i++) {
          node_v = bin_j[v_i].Val;
          node_v_rank = rank_orig.GetDat(node_v);
          if (node_u_rank > node_v_rank)
            beta_ij += 1.0;
          if ((NIdV_temp).IsIn(node_v)) {
            // cout<<"Perfect pair detected: "<< node_u <<"->"<<node_v<<endl;
            perfectpairs.push_back(std::make_pair(node_u, node_v));
            if (node_u_rank > node_v_rank)
              alpha_ij += 1.0;
          } else {
            nonperfectpairs.push_back(std::make_pair(node_u, node_v));
          }
        }
      }
      beta_tilde += beta_ij / tt;
      alpha_tilde += alpha_ij / tt;
    }
    for (uint64_t u_i = 0; u_i < (len_rank_bin_i - 1); u_i++) {
      node_u = bin_i[u_i].Val;
      for (uint64_t u_ii = u_i + 1; u_ii < len_rank_bin_i; u_ii++) {
        node_uu = bin_i[u_ii].Val;
        samebinpairs.push_back(std::make_pair(node_u, node_uu));
      }
    }
  }
  // alpha_tilde is proportion of no. of perfect pairs that are correct
  // beta_tilde is proportion of no. of correct pairs
  uint64_t no_pairs = no_pairs_DAG(rank_parallel);
  float temp_theta_rpeel = beta_tilde * tt / no_pairs;
  return std::make_tuple(alpha_tilde, count_perfect_pairs_nlzd, beta_tilde,
                         temp_theta_rpeel, nonperfectpairs, perfectpairs,
                         samebinpairs);
}

std::tuple<double, double, double, double, std::vector<std::pair<int, int>>,
           std::vector<std::pair<int, int>>, std::vector<std::pair<int, int>>>
FindRhoCountPerfectPeel_int(PNEANet G_DAG, TUInt64UInt64VH &rank_parallel,
                            TUInt64UInt64H &rank_orig, uint64_t &G_no_nodes) {
  uint64_t rank_new_len = rank_parallel.Len();
  uint64_t len_rank_bin_i;
  uint64_t len_rank_bin_j;
  double beta_tilde = 0;
  double alpha_tilde = 0;
  double beta_ij = 0;
  double alpha_ij = 0;
  uint64_t node_v, node_u, node_u_rank, node_v_rank, node_uu;
  TUInt64V bin_i, bin_j;
  double tt = ((double)(G_no_nodes - 1) * G_no_nodes) / 2.0;
  TIntV NIdV, NIdV_temp, NIdV_temp_1;
  TIntIntVH NId_reachb;
  PNGraph G_tree;
  std::vector<std::pair<int, int>> nonperfectpairs;
  std::vector<std::pair<int, int>> perfectpairs;
  std::vector<std::pair<int, int>> samebinpairs;

  double count_perfect_pairs_nlzd = 0;
  for (PNEANet::TObj::TNodeI NI = G_DAG->BegNI(); NI < G_DAG->EndNI(); NI++) {
    G_tree = TSnap::GetBfsTree(G_DAG, NI.GetId(), true, false);
    G_tree->GetNIdV(NIdV);
    NId_reachb.AddDat(NI.GetId(), NIdV);
    count_perfect_pairs_nlzd += (NIdV.Len() - 1) / tt;
    NIdV.Clr();
  }

  for (int64_t i = (rank_new_len - 1); i >= 1; --i) {
    bin_i = TUInt64V(rank_parallel.GetDat(i));
    len_rank_bin_i = bin_i.Len();
    for (int64_t j = (i - 1); j >= 0; --j) {
      bin_j = TUInt64V(rank_parallel.GetDat(j));
      len_rank_bin_j = bin_j.Len();
      beta_ij = 0;
      alpha_ij = 0;
      for (uint64_t u_i = 0; u_i < len_rank_bin_i; u_i++) {
        node_u = bin_i[u_i].Val;
        node_u_rank = rank_orig.GetDat(node_u);
        NIdV_temp = TIntV(NId_reachb.GetDat(node_u));
        for (uint64_t v_i = 0; v_i < len_rank_bin_j; v_i++) {
          node_v = bin_j[v_i].Val;
          node_v_rank = rank_orig.GetDat(node_v);
          if (node_u_rank > node_v_rank)
            beta_ij += 1.0;
          if ((NIdV_temp).IsIn(node_v)) {
            // cout<<"Perfect pair detected: "<< node_u <<"->"<<node_v<<endl;
            perfectpairs.push_back(std::make_pair(node_u, node_v));
            if (node_u_rank > node_v_rank)
              alpha_ij += 1.0;
          } else {
            nonperfectpairs.push_back(std::make_pair(node_u, node_v));
          }
        }
      }
      beta_tilde += beta_ij / tt;
      alpha_tilde += alpha_ij / tt;
    }
    for (uint64_t u_i = 0; u_i < (len_rank_bin_i - 1); u_i++) {
      node_u = bin_i[u_i].Val;
      for (uint64_t u_ii = u_i + 1; u_ii < len_rank_bin_i; u_ii++) {
        node_uu = bin_i[u_ii].Val;
        samebinpairs.push_back(std::make_pair(node_u, node_uu));
      }
    }
  }
  // alpha_tilde is proportion of no. of perfect pairs that are correct
  // beta_tilde is proportion of no. of correct pairs
  uint64_t no_pairs = no_pairs_DAG(rank_parallel);
  float temp_theta_rpeel = beta_tilde * tt / no_pairs;
  return std::make_tuple(alpha_tilde, count_perfect_pairs_nlzd, beta_tilde,
                         temp_theta_rpeel, nonperfectpairs, perfectpairs,
                         samebinpairs);
}

// Compute Rho and Theta from multiple runs of Sequential algorithm
void with_Sequential_uniform(const int &TimeN, const float &pr_alpha,
                             const int &vec_p_1, const int &vec_p_2,
                             const float &pr_beta, const int &vec_q_1,
                             const int &vec_q_2, const float &pr_delta,
                             const float &pr_gamma, const uint &no_runs_seq,
                             const uint &no_runs) {
  PNEANet G, G_DAG;
  TUInt64UInt64VH rank_parallel;
  TUInt64UInt64H rank_seq_uni;
  uint64_t G_no_nodes;
  TUInt64V bin_sl;

  // G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
  // pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
  // TSnap::SaveEdgeList(G, FName, Desc);
  // // //----------------------------------
  // // G_no_nodes = (uint64_t)G->GetNodes();
  // // TUInt64UInt64H rank_orig;
  // // for(uint64_t i = 0; i < G_no_nodes; i++) {
  // //  rank_orig.AddDat(i) = i;
  // // }
  // TUInt64UInt64H rank_node;
  // std::tie(rank_parallel,rank_node) = find_parallel_rank_returnnoderank(G);
  // G_DAG = create_DAG(G,rank_parallel);
  // // float
  // rho_perf_tt,theta_perf_tt,rho_guess_tt,theta_guess_tt,rho_peel_tt,theta_peel_tt;
  // //
  // std::tie(rho_perf_tt,theta_perf_tt,rho_guess_tt,theta_guess_tt,rho_peel_tt,theta_peel_tt)
  // = FindRhoThetaPerfectGuess(G_DAG, rank_parallel, rank_orig, rank_node,
  // G_no_nodes);
  // //----------------------------------
  // write_dot(G,"original_G.dot","G");
  // write_dot(G_DAG,"DAG_G.dot","DAG(G)");

  // //----------------------------------
  // Load edge list
  // const char *FName = "PA_graph.dat";
  // const char *FName = "PA_graph_n5_m2.dat";
  // // const char *Desc = "PA Graph generated via Jithin's technique";
  // G = TSnap::LoadEdgeList<sPNEANet>(FName);
  // G_no_nodes = (uint64_t)G->GetNodes();
  // TUInt64UInt64H rank_orig;
  // for(uint64_t i = 0; i < G_no_nodes; i++) {
  //      rank_orig.AddDat(i) = i;
  // }
  // //----------------------------------
  double rho_perf_tt, count_perfect_pairs_nlzd, rho_peel_tt, theta_peel_tt;
  double rho_prob_lemma_seq = 0.0, rho_peel = 0.0, rho_perf = 0.0;
  double theta_prob_lemma_seq = 0.0, theta_peel = 0.0;

  // Temporary variables
  for (uint kk = 0; kk < no_runs; ++kk) {
    const bool self_loops_allowed = 1;
    // G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
    // pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
    const char *FName = "PA_graph_n5_m2.dat";
    G = TSnap::LoadEdgeList<PNEANet>(FName);
    G_no_nodes = (uint64_t)G->GetNodes();
    TUInt64UInt64H rank_orig;
    for (uint64_t i = 0; i < G_no_nodes; i++) {
      rank_orig.AddDat(i) = i;
    }
    // Find rho and theta of perfect pairs and nonperfect pairs across bins
    std::vector<std::pair<uint64_t, uint64_t>> nonperfectpairs;
    std::vector<std::pair<uint64_t, uint64_t>> perfectpairs;
    std::vector<std::pair<uint64_t, uint64_t>> samebinpairs;
    rank_parallel = find_parallel_rank(G);
    G_DAG = create_DAG(G, rank_parallel);
    std::tie(rho_perf_tt, count_perfect_pairs_nlzd, rho_peel_tt, theta_peel_tt,
             nonperfectpairs, perfectpairs, samebinpairs) =
        FindRhoCountPerfectPeel(G_DAG, rank_parallel, rank_orig, G_no_nodes);
    uint64_t nonperfectpairs_len = nonperfectpairs.size();
    uint64_t samebinpairs_len = samebinpairs.size();
    double rho_prob_lemma_seq_tt = 0.0, rho_prob_lemma_seq_temp = 0.0;
    double theta_prob_lemma_seq_tt = 0.0;
    uint64_t count_prob_lemma_seq = 0;
    double nC2 = ((double)(G_no_nodes - 1) * G_no_nodes) / 2.0;

    // Find rho and theta of pairs guessed by sequential technique

    std::vector<double> emp_average(nonperfectpairs_len, 0.0);
    std::vector<double> emp_average_samebin(nonperfectpairs_len, 0.0);
    uint64_t node_u, node_v;
    for (uint ii = 0; ii < no_runs_seq; ii++) {
      std::cout << kk << "th run of Sequential technique: " << ii << std::endl;
      rank_seq_uni = find_seq_rank_unif_Hout_node_rank(G);
      for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
        node_u = nonperfectpairs[jj].first;
        node_v = nonperfectpairs[jj].second;
        if (rank_seq_uni.GetDat(node_u) > rank_seq_uni.GetDat(node_v))
          emp_average[jj] += 1.0 / no_runs_seq;
      }
      for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
        node_u = samebinpairs[jj].first;
        node_v = samebinpairs[jj].second;
        if (rank_seq_uni.GetDat(node_u) > rank_seq_uni.GetDat(node_v))
          emp_average_samebin[jj] += 1.0 / no_runs_seq;
      }
    }

    for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
      cout << "(" << nonperfectpairs[jj].first << ","
           << nonperfectpairs[jj].second << "): " << emp_average[jj] << endl;
      if (emp_average[jj] >= 0.8) {
        count_prob_lemma_seq += 1;
        node_u = nonperfectpairs[jj].first;
        node_v = nonperfectpairs[jj].second;
        if (rank_orig.GetDat(node_u) > rank_orig.GetDat(node_v))
          rho_prob_lemma_seq_temp += 1.0 / nC2;
      }
      // This will change the Peeling bins, and may not give consistent linear
      // order
      // if(emp_average[jj] <= 0.2) {
      //  count_prob_lemma_seq += 1;
      //  node_u = nonperfectpairs[jj].first;
      //  node_v = nonperfectpairs[jj].second;
      //  if(rank_orig.GetDat(node_u) < rank_orig.GetDat(node_v))
      //      rho_prob_lemma_seq_temp +=1.0/nC2;
      // }
    }

    for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
      cout << "(" << samebinpairs[jj].first << "," << samebinpairs[jj].second
           << "): " << emp_average_samebin[jj] << endl;
      if (emp_average_samebin[jj] >= 0.8) {
        count_prob_lemma_seq += 1;
        node_u = samebinpairs[jj].first;
        node_v = samebinpairs[jj].second;
        if (rank_orig.GetDat(node_u) > rank_orig.GetDat(node_v))
          rho_prob_lemma_seq_temp += 1.0 / nC2;
      }
      if (emp_average_samebin[jj] <= 0.2) {
        count_prob_lemma_seq += 1;
        node_u = samebinpairs[jj].first;
        node_v = samebinpairs[jj].second;
        if (rank_orig.GetDat(node_u) < rank_orig.GetDat(node_v))
          rho_prob_lemma_seq_temp += 1.0 / nC2;
      }
    }

    rho_prob_lemma_seq_tt = rho_prob_lemma_seq_temp + rho_perf_tt;
    theta_prob_lemma_seq_tt =
        (rho_prob_lemma_seq_temp + rho_perf_tt) /
        (((double)count_prob_lemma_seq / nC2) + count_perfect_pairs_nlzd);
    rho_peel += rho_peel_tt / no_runs;
    rho_perf += rho_perf_tt / no_runs;
    rho_prob_lemma_seq += rho_prob_lemma_seq_tt / no_runs;
    theta_peel += theta_peel_tt / no_runs;
    theta_prob_lemma_seq += theta_prob_lemma_seq_tt / no_runs;
  }
  TStr FuncStr = TStr::Fmt("%s:graph", __func__);
  PrintGStats(FuncStr.CStr(), G);
  printf("rho(Rperf): %f, theta(Rperf): %f \n", rho_perf, 1.0);
  printf("rho(Rguess): %f, theta(Rguess): %f \n", rho_prob_lemma_seq,
         theta_prob_lemma_seq);
  printf("rho(Rpeel): %f, theta(Rpeel): %f \n", rho_peel, theta_peel);
}

// Ratio form of importance sampling
// a. Write P matrix
// b. Find Rho and Theta by putting threshold for estimated p_{u,v}
void with_importance_sampling(const int &TimeN, const float &pr_alpha,
                              const int &vec_p_1, const int &vec_p_2,
                              const float &pr_beta, const int &vec_q_1,
                              const int &vec_q_2, const float &pr_delta,
                              const float &pr_gamma, const uint &no_runs_seq,
                              const uint &no_runs) {
  PNEANet G, G_DAG;
  TUInt64UInt64VH rank_parallel;
  TUInt64UInt64H rank_seq_uni;
  uint64_t G_no_nodes;
  TUInt64V bin_sl;

  double rho_perf_tt = 0, count_perfect_pairs_nlzd = 0, rho_peel_tt = 0,
         theta_peel_tt = 0;
  double rho_prob_lemma_seq = 0.0, rho_peel = 0.0, rho_perf = 0.0;
  double theta_prob_lemma_seq = 0.0, theta_peel = 0.0;
  double log_prob_seq_rank;
  for (uint kk = 0; kk < no_runs; ++kk) {
    // const bool self_loops_allowed = 1;
    // G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
    // pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
    // // ***SAVE GRAPH
    // const char *FName = "PA_graph.dat";
    // const char *Desc = "PA Graph generated via Jithin's technique";
    // TSnap::SaveEdgeList(G, FName, Desc);
    // *** LOAD GRAPH
    // const char *FName = "PA_graph_n50_m3_good.dat";
    // const char *FName = "PA_graph_n10_m2.dat";
    const char *FName = "PA_graph_n5_m2.dat";
    G = TSnap::LoadEdgeList<PNEANet>(FName);

    G_no_nodes = (uint64_t)G->GetNodes();
    TUInt64UInt64H rank_orig;
    for (uint64_t i = 0; i < G_no_nodes; i++) {
      rank_orig.AddDat(i) = i;
    }

    // Find rho and theta of perfect pairs and nonperfect pairs across bins
    std::vector<std::pair<uint64_t, uint64_t>> nonperfectpairs;
    std::vector<std::pair<uint64_t, uint64_t>> perfectpairs;
    std::vector<std::pair<uint64_t, uint64_t>> samebinpairs;
    rank_parallel = find_parallel_rank(G);
    G_DAG = create_DAG(G, rank_parallel);
    std::tie(rho_perf_tt, count_perfect_pairs_nlzd, rho_peel_tt, theta_peel_tt,
             nonperfectpairs, perfectpairs, samebinpairs) =
        FindRhoCountPerfectPeel(G_DAG, rank_parallel, rank_orig, G_no_nodes);
    uint64_t nonperfectpairs_len = nonperfectpairs.size();
    uint64_t perfectpairs_len = perfectpairs.size();
    uint64_t samebinpairs_len = samebinpairs.size();
    double rho_prob_lemma_seq_tt = 0.0, rho_prob_lemma_seq_temp = 0.0;
    double theta_prob_lemma_seq_tt = 0.0;
    uint64_t count_prob_lemma_seq = 0;
    double nC2 = ((double)(G_no_nodes - 1) * G_no_nodes) / 2.0;

    // Find rho and theta of pairs guessed by sequential technique
    std::vector<double> emp_average(nonperfectpairs_len, 0.0);
    std::vector<double> emp_average_samebin(samebinpairs_len, 0.0);
    double norm_average = 0.0;
    uint64_t node_u, node_v;
    for (uint ii = 0; ii < no_runs_seq; ii++) {
      std::cout << kk << "th run of Sequential technique: " << ii << std::endl;
      std::tie(rank_seq_uni, log_prob_seq_rank) =
          find_seq_rank_unif_n_prob_Hout_node_rank(G);
      cout << "prob_seq_rank: " << exp(log_prob_seq_rank) << endl;
      for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
        node_u = nonperfectpairs[jj].first;
        node_v = nonperfectpairs[jj].second;
        if (rank_seq_uni.GetDat(node_u) > rank_seq_uni.GetDat(node_v)) {
          emp_average[jj] += (1.0 / exp(log_prob_seq_rank));
        }
      }
      norm_average += (1.0 / exp(log_prob_seq_rank));
      for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
        node_u = samebinpairs[jj].first;
        node_v = samebinpairs[jj].second;
        if (rank_seq_uni.GetDat(node_u) > rank_seq_uni.GetDat(node_v))
          emp_average_samebin[jj] += (1.0 / exp(log_prob_seq_rank));
      }
    }
    std::transform(emp_average.begin(), emp_average.end(), emp_average.begin(),
                   std::bind1st(std::multiplies<double>(), 1.0 / norm_average));
    std::transform(emp_average_samebin.begin(), emp_average_samebin.end(),
                   emp_average_samebin.begin(),
                   std::bind1st(std::multiplies<double>(), 1.0 / norm_average));

    // Writing results to a file
    const char *outfile = "P_matrix_importance_sampling_sequential.txt";
    std::ofstream ost{outfile};
    for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
      node_u = nonperfectpairs[jj].first;
      node_v = nonperfectpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t" << emp_average[jj] << "\n";
      ost << node_v << "\t" << node_u << "\t" << 1 - emp_average[jj] << "\n";
    }
    for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
      node_u = samebinpairs[jj].first;
      node_v = samebinpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t" << emp_average_samebin[jj]
          << "\n";
      ost << node_v << "\t" << node_u << "\t" << 1 - emp_average_samebin[jj]
          << "\n";
    }
    for (uint64_t jj = 0; jj < perfectpairs_len; jj++) {
      node_u = perfectpairs[jj].first;
      node_v = perfectpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t"
          << "1"
          << "\n";
      ost << node_v << "\t" << node_u << "\t"
          << "0"
          << "\n";
    }
    ost.close();
    //-------------------------

    for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
      cout << "(" << nonperfectpairs[jj].first << ","
           << nonperfectpairs[jj].second << "): " << emp_average[jj] << endl;
      if (emp_average[jj] >= 0.8) {
        count_prob_lemma_seq += 1;
        node_u = nonperfectpairs[jj].first;
        node_v = nonperfectpairs[jj].second;
        if (rank_orig.GetDat(node_u) > rank_orig.GetDat(node_v))
          rho_prob_lemma_seq_temp += 1.0 / nC2;
      }
      // This will change the Peeling bins, and may not give consistent linear
      // order
      // if(emp_average[jj] <= 0.2) {
      //  count_prob_lemma_seq += 1;
      //  node_u = nonperfectpairs[jj].first;
      //  node_v = nonperfectpairs[jj].second;
      //  if(rank_orig.GetDat(node_u) < rank_orig.GetDat(node_v))
      //      rho_prob_lemma_seq_temp +=1.0/nC2;
      // }
    }

    // SAME BIN PAIRS
    for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
      cout << "(" << samebinpairs[jj].first << "," << samebinpairs[jj].second
           << "): " << emp_average_samebin[jj] << endl;
      if (emp_average_samebin[jj] >= 0.8) {
        count_prob_lemma_seq += 1;
        node_u = samebinpairs[jj].first;
        node_v = samebinpairs[jj].second;
        if (rank_orig.GetDat(node_u) > rank_orig.GetDat(node_v))
          rho_prob_lemma_seq_temp += 1.0 / nC2;
      }
      //  if(emp_average_samebin[jj]<=0.2) {
      //          count_prob_lemma_seq += 1;
      //          node_u = samebinpairs[jj].first;
      //          node_v = samebinpairs[jj].second;
      //          if(rank_orig.GetDat(node_u) < rank_orig.GetDat(node_v))
      //                  rho_prob_lemma_seq_temp +=1.0/nC2;
      //  }
    }

    rho_prob_lemma_seq_tt = rho_prob_lemma_seq_temp + rho_perf_tt;
    theta_prob_lemma_seq_tt =
        (rho_prob_lemma_seq_temp + rho_perf_tt) /
        (((double)count_prob_lemma_seq / nC2) + count_perfect_pairs_nlzd);
    rho_peel += rho_peel_tt / no_runs;
    rho_perf += rho_perf_tt / no_runs;
    rho_prob_lemma_seq += rho_prob_lemma_seq_tt / no_runs;
    theta_peel += theta_peel_tt / no_runs;
    theta_prob_lemma_seq += theta_prob_lemma_seq_tt / no_runs;
  }
  TStr FuncStr = TStr::Fmt("%s:graph", __func__);
  PrintGStats(FuncStr.CStr(), G);
  printf("rho(Rperf): %f, theta(Rperf): %f \n", rho_perf, 1.0);
  printf("rho(Rguess): %f, theta(Rguess): %f \n", rho_prob_lemma_seq,
         theta_prob_lemma_seq);
  printf("rho(Rpeel): %f, theta(Rpeel): %f \n", rho_peel, theta_peel);
  printf("Density(Rpeel): %f \n", rho_peel_tt / theta_peel_tt);
}

// Jithin's algorithm
void estimate_puv_KK_algorithm(const int &TimeN, const float &pr_alpha,
                               const int &vec_p_1, const int &vec_p_2,
                               const float &pr_beta, const int &vec_q_1,
                               const int &vec_q_2, const float &pr_delta,
                               const float &pr_gamma, const uint &no_runs_MC,
                               const uint &no_runs) {
  PNEANet G, G_DAG;
  TUInt64UInt64VH rank_parallel;
  TUInt64UInt64H rank_seq_uni;
  uint64_t G_no_nodes;
  TUInt64V bin_sl;

  double rho_perf_tt = 0, count_perfect_pairs_nlzd = 0, rho_peel_tt = 0,
         theta_peel_tt = 0;
  std::map<int, std::pair<int, int>> rank_pair;
  const char *outfile_theta_delta = "theta_delta.txt";
  std::ofstream ost1{outfile_theta_delta};

  for (uint kk = 0; kk < no_runs; ++kk) {
    // const bool self_loops_allowed = 1;
    // G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
    // pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
    // // ***SAVE GRAPH
    // const char *FName = "PA_graph.dat";
    // const char *Desc = "PA Graph generated via Jithin's technique";
    // TSnap::SaveEdgeList(G, FName, Desc);
    // *** LOAD GRAPH
    // const char *FName = "PA_graph_n50_m3.dat";
    const char *FName = "PA_graph_n10_m2.dat";
    // const char *FName = "PA_graph_n50_m3_good.dat";
    G = TSnap::LoadEdgeList<PNEANet>(FName);

    G_no_nodes = (uint64_t)G->GetNodes();
    TUInt64UInt64H rank_orig;
    for (uint64_t i = 0; i < G_no_nodes; i++) {
      rank_orig.AddDat(i) = i;
    }
    std::vector<std::pair<int, int>> nonperfectpairs;
    std::vector<std::pair<int, int>> perfectpairs;
    std::vector<std::pair<int, int>> samebinpairs;
    rank_parallel = find_parallel_rank(G);
    G_DAG = create_DAG(G, rank_parallel);
    std::tie(rho_perf_tt, count_perfect_pairs_nlzd, rho_peel_tt, theta_peel_tt,
             nonperfectpairs, perfectpairs, samebinpairs) =
        FindRhoCountPerfectPeel_int(G_DAG, rank_parallel, rank_orig,
                                    G_no_nodes);
    rank_pair = find_seq_rank_unif_mappair(G);
    std::set<int> interchangable_pairs;
    uint64_t perfectpairs_len = perfectpairs.size();
    uint64_t nonperfectpairs_len = nonperfectpairs.size();
    uint64_t samebinpairs_len = samebinpairs.size();
    for (auto &x : rank_pair) {
      if (std::find(perfectpairs.begin(), perfectpairs.end(), x.second) ==
          perfectpairs.end()) {
        interchangable_pairs.insert(x.first);
      }
    }
    std::vector<double> emp_average(nonperfectpairs_len, 0.0);
    std::vector<double> emp_average_samebin(samebinpairs_len, 0.0);
    double normalization = 0.0;
    // uint mixing_time = std::pow(G_no_nodes,3)*log(G_no_nodes);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution rv_alpha(0.8);

    // print_vector(interchangable_pairs,"Original interchangable_pairs");
    // print_set(interchangable_pairs,"Original interchangable_pairs");
    // print_map(rank_pair,"Original rank_pair");
    for (uint64_t kkkk = 0; kkkk < no_runs_MC; kkkk++) {
      cout << "run MC:" << kkkk << endl;
      // if (rv_alpha(gen)) {
      uint sel_element = *select_randomly(interchangable_pairs.begin(),
                                          interchangable_pairs.end());
      // uint64_t deg_present_node = interchangable_pairs.size();
      // cout<<"sel_element: "<<sel_element<<endl;
      int temp_1, temp_2, temp_11, temp_22;
      std::pair<int, int> new_next_sel_element, new_prev_sel_element;
      std::tie(temp_1, temp_2) = rank_pair[sel_element];
      std::swap(temp_1, temp_2);
      rank_pair[sel_element] = std::make_pair(temp_1, temp_2);
      if (sel_element > 0) {
        std::tie(temp_11, temp_22) = rank_pair[sel_element - 1];
        new_prev_sel_element = std::make_pair(temp_2, temp_22);
        rank_pair[sel_element - 1] = new_prev_sel_element;
        if (std::find(perfectpairs.begin(), perfectpairs.end(),
                      new_prev_sel_element) != perfectpairs.end()) {
          // interchangable_pairs.erase(std::remove(interchangable_pairs.begin(),
          // interchangable_pairs.end(), (sel_element-1)),
          // interchangable_pairs.end() ); //Assuming (sel_element-1) already
          // present
          interchangable_pairs.erase((sel_element - 1));
        } else {
          interchangable_pairs.insert((sel_element - 1));
        }
      }
      if (sel_element < (G_no_nodes - 2)) {
        std::tie(temp_11, temp_22) = rank_pair[sel_element + 1];
        new_next_sel_element = std::make_pair(temp_11, temp_1);
        rank_pair[sel_element + 1] = new_next_sel_element;
        if (std::find(perfectpairs.begin(), perfectpairs.end(),
                      new_next_sel_element) != perfectpairs.end()) {
          // interchangable_pairs.erase( std::remove(
          // interchangable_pairs.begin(), interchangable_pairs.end(),
          // (sel_element+1) ), interchangable_pairs.end() );
          interchangable_pairs.erase((sel_element + 1));
        } else {
          interchangable_pairs.insert((sel_element + 1));
        }
      }
      // print_map(rank_pair,"rank_pair");
      // print_vector(interchangable_pairs,"interchangable_pairs");
      // print_set(interchangable_pairs,"interchangable_pairs");
      // }
      std::vector<std::pair<int, int>> rank_pair_second;
      for (auto &s : rank_pair)
        rank_pair_second.push_back(s.second);

      // // CHECK: WHETHER ZERP PROB PERF PAIRS NEED TO BE SEARCHED
      // for (uint64_t jj=0; jj< perfectpairs_len; jj++) {
      //      // cout<<"Entered 0 perfect pair check"<<endl;
      //      int aa = perfectpairs[jj].first;
      //      int bb = perfectpairs[jj].second;
      //      std::pair<int,int> pair_temp = std::make_pair(bb,aa);
      //      if (std::find(rank_pair_second.begin(),
      //      rank_pair_second.end(),pair_temp) != rank_pair_second.end()) {
      //              cout<<"Add this constraint!!"<<endl;
      //      }
      // }
      // // *****************************************************

      // //******  Taking into account mixing time *******
      // if (kkkk > mixing_time) {
      //      for (uint64_t jj=0; jj<nonperfectpairs_len; jj++) {
      //              if( std::find(rank_pair_second.begin(),
      //              rank_pair_second.end(), nonperfectpairs[jj]) !=
      //              rank_pair_second.end()) {
      //                      emp_average[jj] += (1.0/(no_runs_MC -
      //                      mixing_time));
      //              }
      //      }
      //      for (uint64_t jj=0; jj<samebinpairs_len; jj++) {
      //              if( std::find(rank_pair_second.begin(),
      //              rank_pair_second.end(), samebinpairs[jj]) !=
      //              rank_pair_second.end()) {
      //                      emp_average_samebin[jj] += (1.0/(no_runs_MC -
      //                      mixing_time));
      //              }
      //      }
      // }
      // // *******
      uint64_t deg_present_node = interchangable_pairs.size();
      for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
        if (std::find(rank_pair_second.begin(), rank_pair_second.end(),
                      nonperfectpairs[jj]) != rank_pair_second.end()) {
          emp_average[jj] += 1.0 / deg_present_node;
        }
      }
      for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
        if (std::find(rank_pair_second.begin(), rank_pair_second.end(),
                      samebinpairs[jj]) != rank_pair_second.end()) {
          emp_average_samebin[jj] += 1.0 / deg_present_node;
        }
      }
      rank_pair_second.clear();
      normalization += 1.0 / deg_present_node;
      // printf("Press Return to continue\n");
      // std::cin.get();
    }
    std::transform(
        emp_average.begin(), emp_average.end(), emp_average.begin(),
        std::bind1st(std::multiplies<double>(), 1.0 / normalization));
    std::transform(
        emp_average_samebin.begin(), emp_average_samebin.end(),
        emp_average_samebin.begin(),
        std::bind1st(std::multiplies<double>(), 1.0 / normalization));

    // Writing results to a file
    const char *outfile = "P_matrix.txt";
    std::ofstream ost{outfile};
    uint node_u, node_v;
    for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
      node_u = nonperfectpairs[jj].first;
      node_v = nonperfectpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t" << emp_average[jj] << "\n";
      ost << node_v << "\t" << node_u << "\t" << 1 - emp_average[jj] << "\n";
    }
    for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
      node_u = samebinpairs[jj].first;
      node_v = samebinpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t" << emp_average_samebin[jj]
          << "\n";
      ost << node_v << "\t" << node_u << "\t" << 1 - emp_average_samebin[jj]
          << "\n";
    }
    for (uint64_t jj = 0; jj < perfectpairs_len; jj++) {
      node_u = perfectpairs[jj].first;
      node_v = perfectpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t"
          << "1"
          << "\n";
      ost << node_v << "\t" << node_u << "\t"
          << "0"
          << "\n";
    }
    ost.close();
    cout << "Written P matrix to file" << endl;
    //-------------------------
    // theta and density into file:
    ost1 << theta_peel_tt << "\t" << rho_peel_tt / theta_peel_tt << "\n";
  }
  ost1.close();
  printf("rho(Rperf): %f, theta(Rperf): %f \n", rho_perf_tt, 1.0);
  printf("rho(Rpeel): %f, theta(Rpeel): %f \n", rho_peel_tt, theta_peel_tt);
  printf("Density(Rpeel): %f \n", rho_peel_tt / theta_peel_tt);
}

// Using my MCMC technique for estimating p_uv.
// It also contains the code to generate the points of following estimators:
// Peeling, Peeling+, Precision-1, degree based binning
void with_puv_jks_algorithm(const int &TimeN, const float &pr_alpha,
                            const int &vec_p_1, const int &vec_p_2,
                            const float &pr_beta, const int &vec_q_1,
                            const int &vec_q_2, const float &pr_delta,
                            const float &pr_gamma, const uint &no_runs_MC,
                            const uint &no_runs, const uint &width_box) {
  PNEANet G, G_DAG;
  TUInt64UInt64VH rank_parallel;
  TUInt64UInt64H rank_seq_uni;
  uint64_t G_no_nodes, G_no_edges;
  TUInt64V bin_sl;

  std::map<int, std::pair<int, int>> rank_pair;
  const char *outfile_theta_delta_peel = "./temp_model/theta_delta_peel.txt";
  const char *outfile_theta_delta_peelp = "./temp_model/theta_delta_peel+.txt";
  const char *outfile_theta_delta_perf = "./temp_model/theta_delta_perf.txt";
  const char *outfile_theta_delta_degree = "./temp_model/theta_delta_deg.txt";
  const char *outfile_theta_delta_ML = "./temp_model/theta_delta_ML.txt";

  std::ofstream ost_peel{outfile_theta_delta_peel};
  std::ofstream ost_peelp{outfile_theta_delta_peelp};
  std::ofstream ost_perf{outfile_theta_delta_perf};
  std::ofstream ost_deg{outfile_theta_delta_degree};
  std::ofstream ost_ML{outfile_theta_delta_ML};

  double rho_perf_tt = 0, count_perfect_pairs_nlzd = 0, rho_peel_tt = 0,
         theta_peel_tt = 0;
  double rho_prob_lemma_seq = 0.0, rho_peel = 0.0, rho_perf = 0.0,
         density_peel = 0, density_guess = 0;
  double theta_prob_lemma_seq = 0.0, theta_peel = 0.0;

  double rho_deg = 0, theta_deg = 0;
  double rho_ML = 0, theta_ML = 0;
  double rho_peelp = 0, theta_peelp = 0, density_peelp = 0;
  // double nC2 = ((double)(G_no_nodes - 1) * G_no_nodes) / 2.0;

  const int no_threads = 1;

  // Forf finding Recall Matrix we define G_no_nodes here. Note that only the
  // classical PA graph model is allowed here.
  // create boxes of equal width_box
  G_no_nodes = TimeN;
  uint64_t total_boxes = std::ceil(G_no_nodes / (double)width_box);
  std::vector<std::vector<double>> recall_matrix(
      total_boxes, std::vector<double>(total_boxes, 0.0));
  std::vector<std::vector<double>> recall_matrix_deg(
      total_boxes, std::vector<double>(total_boxes, 0.0));

  for (uint kk = 0; kk < no_runs; ++kk) {
    double rho_prob_lemma_seq_tt = 0.0, rho_prob_lemma_seq_temp = 0.0;
    double theta_prob_lemma_seq_tt = 0.0;

    double theta_peelp_temp = 0, rho_peelp_temp = 0;
    uint64_t count_prob_lemma_seq = 0;

    const bool self_loops_allowed = 1;
    G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
                                        pr_beta, vec_q_1, vec_q_2, pr_delta,
                                        pr_gamma, self_loops_allowed);
    //  G = GenPrefAttachGeneral(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta,
    //  vec_q_1, vec_q_2, pr_delta, pr_gamma);
    // // ***SAVE GRAPH
    // const char *FName = "PA_graph.dat";
    // const char *Desc = "PA Graph generated via Jithin's technique";
    // TSnap::SaveEdgeList(G, FName, Desc);
    //  write_dot(G,"G.dot","G");
    // *** LOAD GRAPH
    // const char *FName = "PA_graph_n50_m3.dat";
    // const char *FName = "PA_graph_n10_m2.dat";
    // const char *FName = "PA_graph_n5_m2.dat";
    // const char *FName = "PA_graph_n50_m3_good.dat";
    // G = TSnap::LoadEdgeList<PNEANet>(FName);
    // write_dot(G, "PA_graph_n5_m2.dot", "G");

    G_no_nodes = (uint64_t)G->GetNodes();
    G_no_edges = (uint64_t)G->GetEdges();

    double nC2 = ((double)(G_no_nodes - 1) * G_no_nodes) / 2.0;

    TUInt64UInt64H rank_orig;
    for (uint64_t i = 0; i < G_no_nodes; i++) {
      rank_orig.AddDat(i) = i;
    }

    Tuint64uint64M rank_orig_STL;
    for (uint64_t i = 0; i < G_no_nodes; i++) {
      rank_orig_STL.insert({i, i});
    }

    std::vector<std::pair<int, int>> nonperfectpairs;
    std::vector<std::pair<int, int>> perfectpairs;
    std::vector<std::pair<int, int>> samebinpairs;

    // Degree based estimator
    uint64_t rho_unnormalized_deg, nopairsDAG_deg;
    double theta_deg_temp = 0, rho_deg_temp = 0;
    Tuint64uint64VM rank_degree_bins;
    Tuint64uint64M rank_node_deg;
    std::tie(rank_degree_bins, rank_node_deg) = find_degree_rank_STL(G);
    std::tie(rho_unnormalized_deg, nopairsDAG_deg) =
        find_prec_recall_real_networks_STL(rank_degree_bins, rank_orig_STL,
                                           G_no_nodes, no_threads);
    rho_deg_temp =
        2 * rho_unnormalized_deg / ((double)G_no_nodes * (G_no_nodes - 1));
    theta_deg_temp = rho_unnormalized_deg / (double)nopairsDAG_deg;

    // ML estimator
    uint64_t rho_unnormalized_ML, nopairsDAG_ML;
    TUInt64UInt64VH rank_ML = find_seq_rank_unif_VHout(G);
    std::tie(rho_unnormalized_ML, nopairsDAG_ML) =
        find_prec_recall_real_networks(rank_ML, rank_orig, G_no_nodes);
    double rho_ML_temp =
        2 * rho_unnormalized_ML / ((double)G_no_nodes * (G_no_nodes - 1));
    double theta_ML_temp = rho_unnormalized_ML / (double)nopairsDAG_ML;

    // Peel estimator
    rank_parallel = find_parallel_rank(G);
    G_DAG = create_DAG(G, rank_parallel);
    std::tie(rho_perf_tt, count_perfect_pairs_nlzd, rho_peel_tt, theta_peel_tt,
             nonperfectpairs, perfectpairs, samebinpairs) =
        FindRhoCountPerfectPeel_int(G_DAG, rank_parallel, rank_orig,
                                    G_no_nodes);
    rank_pair = find_seq_rank_unif_mappair(G);
    // print_vector(perfectpairs, "perfect pairs");

    // Peel+ estimator
    TUInt64UInt64H rank_node;
    uint64_t rho_extra_unnnormalized, nopairsDAGextra;
    uint64_t no_orig_pairs_same_bin = 0, no_orig_pairs_across_bin = 0,
             no_orig_pairs = 0;
    uint64_t rho_unnormalized, nopairsDAG;
    std::tie(rank_parallel, rank_node) = find_parallel_rank_returnnoderank(G);
    std::tie(rho_unnormalized, nopairsDAG) =
        find_prec_recall_real_networks(rank_parallel, rank_orig, G_no_nodes);
    std::tie(rho_extra_unnnormalized, nopairsDAGextra, no_orig_pairs_same_bin) =
        Count_extra_pairs_real_networks(G, G_no_nodes, rank_orig, rank_parallel,
                                        rank_node);
    no_orig_pairs_across_bin = nopairsDAG;
    no_orig_pairs = no_orig_pairs_same_bin + no_orig_pairs_across_bin;
    rho_peelp_temp =
        (rho_extra_unnnormalized + rho_unnormalized) / (double)no_orig_pairs;
    theta_peelp_temp = (rho_extra_unnnormalized + rho_unnormalized) /
                       (double)(no_orig_pairs_across_bin + nopairsDAGextra);

    // Recall matrix for Peeling
    std::vector<std::vector<double>> recall_matrix_temp =
        find_recall_plot(rank_node, width_box, total_boxes);
    for (uint64_t i = 0; i < total_boxes; i++)
      for (uint64_t j = 0; j < total_boxes; j++) {
        if (kk == 0)
          recall_matrix[i][j] = recall_matrix_temp[i][j];
        else
          recall_matrix[i][j] = ((kk - 1) / (double)kk) * recall_matrix[i][j] +
                                (recall_matrix_temp[i][j] / (double)kk);
      }

    // // Recall matrix for Peeling
    // std::vector<std::vector<double>> recall_matrix_temp =
    //     find_recall_plot(rank_node, width_box, total_boxes);
    // for (uint64_t i = 0; i < total_boxes; i++)
    //   for (uint64_t j = 0; j < total_boxes; j++) {
    //     if (kk == 1)
    //       recall_matrix[i][j] = recall_matrix_temp[i][j];
    //     else
    //       recall_matrix[i][j] = ((kk - 1) / (double)kk) * recall_matrix[i][j]
    //       +
    //                             (recall_matrix_temp[i][j] / (double)kk);
    //   }

    // Recall matrix for Degree
    std::vector<std::vector<double>> recall_matrix_deg_temp =
        find_recall_plot(rank_node_deg, width_box, total_boxes);
    for (uint64_t i = 0; i < total_boxes; i++)
      for (uint64_t j = 0; j < total_boxes; j++) {
        if (kk == 0)
          recall_matrix_deg[i][j] = recall_matrix_deg_temp[i][j];
        else
          recall_matrix_deg[i][j] =
              ((kk - 1) / (double)kk) * recall_matrix_deg[i][j] +
              (recall_matrix_deg_temp[i][j] / (double)kk);
      }

    // MCMC CALCULATION STARTS HERE
    std::set<int> interchangable_pairs;
    uint64_t perfectpairs_len = perfectpairs.size();
    uint64_t nonperfectpairs_len = nonperfectpairs.size();
    uint64_t samebinpairs_len = samebinpairs.size();
    // print_map(rank_pair,"rank_pair");
    // print_vector(perfectpairs,"perfectpairs");
    // print_vector(nonperfectpairs,"nonperfectpairs");
    // print_vector(samebinpairs,"samebinpairs");
    for (auto &x : rank_pair) {
      if (std::find(perfectpairs.begin(), perfectpairs.end(), x.second) ==
          perfectpairs.end()) {
        interchangable_pairs.insert(x.first);
      }
    }
    std::vector<double> emp_average(nonperfectpairs_len, 0.0);
    std::vector<double> emp_average_samebin(samebinpairs_len, 0.0);
    double normalization = 0.0;

    // std::random_device rd;
    // std::mt19937 gen(rd());
    // std::bernoulli_distribution rv_alpha(0.8);

    for (uint64_t kkkk = 0; kkkk < no_runs_MC; kkkk++) {
      cout << "run MC:" << kkkk << endl;
      uint sel_element = *select_randomly(interchangable_pairs.begin(),
                                          interchangable_pairs.end());
      int temp_1, temp_2, temp_11, temp_22;
      std::pair<int, int> new_next_sel_element, new_prev_sel_element,
          new_sel_element;
      std::tie(temp_1, temp_2) = rank_pair[sel_element];
      std::swap(temp_1, temp_2);
      new_sel_element = std::make_pair(temp_1, temp_2);
      rank_pair[sel_element] = new_sel_element;

      if (std::find(perfectpairs.begin(), perfectpairs.end(),
                    new_sel_element) != perfectpairs.end()) {
        interchangable_pairs.erase(sel_element);
      }
      if (sel_element > 0) {
        std::tie(temp_11, temp_22) = rank_pair[sel_element - 1];
        new_prev_sel_element = std::make_pair(temp_2, temp_22);
        rank_pair[sel_element - 1] = new_prev_sel_element;
        if (std::find(perfectpairs.begin(), perfectpairs.end(),
                      new_prev_sel_element) != perfectpairs.end()) {
          interchangable_pairs.erase((sel_element - 1));
        } else {
          interchangable_pairs.insert((sel_element - 1));
        }
      }
      if (sel_element < (G_no_nodes - 2)) {
        std::tie(temp_11, temp_22) = rank_pair[sel_element + 1];
        new_next_sel_element = std::make_pair(temp_11, temp_1);
        rank_pair[sel_element + 1] = new_next_sel_element;
        if (std::find(perfectpairs.begin(), perfectpairs.end(),
                      new_next_sel_element) != perfectpairs.end()) {
          interchangable_pairs.erase((sel_element + 1));
        } else {
          interchangable_pairs.insert((sel_element + 1));
        }
      }
      // std::vector<std::pair<int,int> > rank_pair_second;
      // for (auto &s : rank_pair)
      //      rank_pair_second.push_back(s.second);

      std::vector<int> present_linear_extension;
      for (auto it = rank_pair.begin(); it != rank_pair.end(); ++it) {
        if (it == rank_pair.begin()) {
          present_linear_extension.push_back(((*it).second).second);
          present_linear_extension.push_back(((*it).second).first);
        } else {
          present_linear_extension.push_back(((*it).second).first);
        }
      }

      // cout<<"sel_element: "<<sel_element<<endl;
      // print_map(rank_pair,"rank_pair");
      // print_set(interchangable_pairs,"interchangable_pairs");
      // print_vector(rank_pair_second,"rank_pair_second");

      uint64_t deg_present_node = interchangable_pairs.size();

      for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
        auto pos_newer =
            std::distance(present_linear_extension.begin(),
                          std::find(present_linear_extension.begin(),
                                    present_linear_extension.end(),
                                    nonperfectpairs[jj].first));
        auto pos_older =
            std::distance(present_linear_extension.begin(),
                          std::find(present_linear_extension.begin(),
                                    present_linear_extension.end(),
                                    nonperfectpairs[jj].second));
        // cout<<"pos_newer: "<<pos_newer<<" value:
        // "<<nonperfectpairs[jj].first<<endl; cout<<"pos_older: "<<pos_older<<"
        // value: "<<nonperfectpairs[jj].second<<endl;

        if (pos_newer > pos_older) {
          // cout<<"Entered: nonperfectpairs
          // "<<nonperfectpairs[jj].first<<","<<nonperfectpairs[jj].second<<endl;
          emp_average[jj] += 1.0 / deg_present_node;
        }
      }
      for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
        auto pos_newer = std::distance(
            present_linear_extension.begin(),
            std::find(present_linear_extension.begin(),
                      present_linear_extension.end(), samebinpairs[jj].first));
        auto pos_older = std::distance(
            present_linear_extension.begin(),
            std::find(present_linear_extension.begin(),
                      present_linear_extension.end(), samebinpairs[jj].second));
        // cout<<"pos_newer: "<<pos_newer<<" value:
        // "<<samebinpairs[jj].first<<endl; cout<<"pos_older: "<<pos_older<<"
        // value: "<<samebinpairs[jj].second<<endl;

        if (pos_newer > pos_older) {
          // cout<<"Entered: nonperfectpairs
          // "<<nonperfectpairs[jj].first<<","<<nonperfectpairs[jj].second<<endl;
          emp_average_samebin[jj] += 1.0 / deg_present_node;
        }
      }
      present_linear_extension.clear();
      normalization += 1.0 / deg_present_node;

      // for (uint64_t jj=0; jj<nonperfectpairs_len; jj++) {
      //      if (std::find(rank_pair_second.begin(), rank_pair_second.end(),
      //      nonperfectpairs[jj]) != rank_pair_second.end()) {
      //              // cout<<"Entered: nonperfectpairs
      //              "<<nonperfectpairs[jj].first<<","<<nonperfectpairs[jj].second<<endl;
      //              emp_average[jj] += 1.0/deg_present_node;
      //              // emp_average[jj] += 1.0;
      //      }
      // }
      // for (uint64_t jj=0; jj<samebinpairs_len; jj++) {
      //      if( std::find(rank_pair_second.begin(), rank_pair_second.end(),
      //      samebinpairs[jj]) != rank_pair_second.end()) {
      //              // cout<<"Entered: samebin
      //              "<<samebinpairs[jj].first<<","<<samebinpairs[jj].second<<endl;
      //              emp_average_samebin[jj] += 1.0/deg_present_node;
      //              // emp_average_samebin[jj] += 1.0;
      //      }
      // }
      // rank_pair_second.clear();
      // normalization += 1.0/deg_present_node;
    }
    std::transform(
        emp_average.begin(), emp_average.end(), emp_average.begin(),
        std::bind1st(std::multiplies<double>(), 1.0 / normalization));
    std::transform(
        emp_average_samebin.begin(), emp_average_samebin.end(),
        emp_average_samebin.begin(),
        std::bind1st(std::multiplies<double>(), 1.0 / normalization));
    uint64_t node_u, node_v;

    //// GUESS PAIRS
    //  // This comments are temporary, just for writing the P matrix
    //  for (uint64_t jj=0; jj<nonperfectpairs_len; jj++) {
    //         //
    //         cout<<"("<<nonperfectpairs[jj].first<<","<<nonperfectpairs[jj].second<<"):
    //         "<<emp_average[jj]<<endl; if(emp_average[jj] >= 0.55) {
    //                 count_prob_lemma_seq += 1;
    //                 node_u = nonperfectpairs[jj].first;
    //                 node_v = nonperfectpairs[jj].second;
    //                 if(rank_orig.GetDat(node_u) > rank_orig.GetDat(node_v))
    //                         rho_prob_lemma_seq_temp +=1.0/nC2;
    //         }
    //         // This will change the Peeling bins, and may not give consistent
    //         linear order if(emp_average[jj] <= 0.45) {
    //                 count_prob_lemma_seq += 1;
    //                 node_u = nonperfectpairs[jj].first;
    //                 node_v = nonperfectpairs[jj].second;
    //                 if(rank_orig.GetDat(node_u) < rank_orig.GetDat(node_v))
    //                         rho_prob_lemma_seq_temp +=1.0/nC2;
    //         }
    //  }
    //
    //  // This comments are temporary, just for writing the P matrix
    //  //SAME BIN PAIRS
    //  for (uint64_t jj=0; jj<samebinpairs_len; jj++) {
    //         //
    //         cout<<"("<<samebinpairs[jj].first<<","<<samebinpairs[jj].second<<"):
    //         "<<emp_average_samebin[jj]<<endl;
    //         if(emp_average_samebin[jj]>=0.55) {
    //                 count_prob_lemma_seq += 1;
    //                 node_u = samebinpairs[jj].first;
    //                 node_v = samebinpairs[jj].second;
    //                 if(rank_orig.GetDat(node_u) > rank_orig.GetDat(node_v))
    //                         rho_prob_lemma_seq_temp +=1.0/nC2;
    //         }
    //         if(emp_average_samebin[jj]<=0.45) {
    //                 count_prob_lemma_seq += 1;
    //                 node_u = samebinpairs[jj].first;
    //                 node_v = samebinpairs[jj].second;
    //                 if(rank_orig.GetDat(node_u) < rank_orig.GetDat(node_v))
    //                         rho_prob_lemma_seq_temp +=1.0/nC2;
    //         }
    //  }

    // This comments are temporary, just for writing the P matrix
    cout << "rho_prob_lemma_seq_temp: " << rho_prob_lemma_seq_temp << endl;
    //  rho_prob_lemma_seq_tt = rho_prob_lemma_seq_temp + rho_perf_tt;
    //  theta_prob_lemma_seq_tt = (rho_prob_lemma_seq_temp + rho_perf_tt) /
    //                           (((double)count_prob_lemma_seq/nC2) +
    //                           count_perfect_pairs_nlzd);

    //  rho_prob_lemma_seq += rho_prob_lemma_seq_tt/no_runs;
    //  theta_prob_lemma_seq += theta_prob_lemma_seq_tt/no_runs;
    //  density_guess +=
    //  (rho_prob_lemma_seq_tt/theta_prob_lemma_seq_tt)/no_runs;

    // Writing results to a file
    //  std::string outfile_temp =
    //  "./P_matrices_n50_cooper_frieze_model/P_matrix_RW_ratio_";
    std::string outfile_temp = "./temp_model/P_matrix_RW_ratio_";
    std::string outfile = outfile_temp + std::to_string(kk) + ".txt";
    std::ofstream ost{outfile};
    ost << "#Estimated p_uv according to Lemma 1 (WWW paper)"
        << "\n";
    ost << "#First line is (no_of_nodes,no_of_edges)"
        << "\n";
    ost << G_no_nodes << "\t" << G_no_edges << "\n";
    for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
      node_u = nonperfectpairs[jj].first;
      node_v = nonperfectpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t" << emp_average[jj] << "\n";
      ost << node_v << "\t" << node_u << "\t" << 1 - emp_average[jj] << "\n";
    }
    for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
      node_u = samebinpairs[jj].first;
      node_v = samebinpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t" << emp_average_samebin[jj]
          << "\n";
      ost << node_v << "\t" << node_u << "\t" << 1 - emp_average_samebin[jj]
          << "\n";
    }
    for (uint64_t jj = 0; jj < perfectpairs_len; jj++) {
      node_u = perfectpairs[jj].first;
      node_v = perfectpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t"
          << "1"
          << "\n";
      ost << node_v << "\t" << node_u << "\t"
          << "0"
          << "\n";
    }
    ost.close();
    cout << "Written " << outfile << " to file" << endl;

    // theta and density into file:
    rho_perf += rho_perf_tt / no_runs;

    rho_peel += rho_peel_tt / no_runs;
    theta_peel += theta_peel_tt / no_runs;
    density_peel += (rho_peel_tt / theta_peel_tt) / no_runs;

    rho_peelp += rho_peelp_temp / no_runs;
    theta_peelp += theta_peelp_temp / no_runs;
    density_peelp += (rho_peelp_temp / theta_peelp_temp) / no_runs;

    rho_deg += rho_deg_temp / no_runs;
    theta_deg += theta_deg_temp / no_runs;

    rho_ML += rho_ML_temp / no_runs;
    theta_ML += theta_ML_temp / no_runs;

    ost_peel << theta_peel_tt << "\t" << rho_peel_tt / theta_peel_tt << "\n";
    ost_peelp << theta_peelp_temp << "\t" << rho_peelp_temp / theta_peelp_temp
              << "\n";
    ost_perf << "1.0"
             << "\t" << rho_perf_tt << "\n";
    ost_deg << theta_deg_temp << "\t" << rho_deg_temp / theta_deg_temp << "\n";
    ost_ML << theta_ML_temp << "\t" << rho_ML_temp / theta_ML_temp << "\n";
  }
  ost_peel.close();
  ost_peelp.close();
  ost_perf.close();
  ost_deg.close();
  ost_ML.close();
  printf("rho(Rperf): %f, theta(Rperf): %f, density(Rperf): %f \n", rho_perf,
         1.0, rho_perf);
  printf("rho(Rpeel): %f, theta(Rpeel): %f, density(Rpeel): %f \n", rho_peel,
         theta_peel, density_peel);
  printf("rho(Rpeelp): %f, theta(Rpeelp): %f, density(Rpeelp): %f \n",
         rho_peelp, theta_peelp, density_peelp);
  printf("rho(Rdeg): %f, theta(Rdeg): %f, density(Rdeg): %f \n", rho_deg,
         theta_deg, rho_deg / theta_deg);
  printf("rho(Rseq): %f, theta(Rseq): %f, density(Rseq): %f", rho_ML, theta_ML,
         rho_ML / theta_ML);

  // std::transform(recall_matrix.begin(), recall_matrix.end(),
  //                recall_matrix.begin(),
  //                std::bind1st(std::multiplies<double>(), 1.0 / no_runs));

  std::string outfile_recall = "./temp_model/Recall_matrix.txt";
  std::ofstream ost_recall{outfile_recall};
  for (int i = 0; i < total_boxes; i++)
    for (int j = 0; j < total_boxes; j++) {
      ost_recall << i << "\t" << j << "\t" << recall_matrix[i][j] << "\n";
    }
  ost_recall.close();

  std::string outfile_recall_deg = "./temp_model/Recall_deg_matrix.txt";
  std::ofstream ost_recall_deg{outfile_recall_deg};
  for (int i = 0; i < total_boxes; i++)
    for (int j = 0; j < total_boxes; j++) {
      ost_recall_deg << i << "\t" << j << "\t" << recall_matrix_deg[i][j]
                     << "\n";
    }
  ost_recall_deg.close();
}

// Using K-K technique
void with_puv_KK_algorithm(const int &TimeN, const float &pr_alpha,
                           const int &vec_p_1, const int &vec_p_2,
                           const float &pr_beta, const int &vec_q_1,
                           const int &vec_q_2, const float &pr_delta,
                           const float &pr_gamma, const uint &no_runs_MC,
                           const uint &no_runs) {
  PNEANet G, G_DAG;
  TUInt64UInt64VH rank_parallel;
  TUInt64UInt64H rank_seq_uni;
  uint64_t G_no_nodes;
  TUInt64V bin_sl;

  std::map<int, std::pair<int, int>> rank_pair;
  const char *outfile_theta_delta = "theta_delta.txt";
  std::ofstream ost1{outfile_theta_delta};

  double rho_perf_tt = 0, count_perfect_pairs_nlzd = 0, rho_peel_tt = 0,
         theta_peel_tt = 0;
  double rho_prob_lemma_seq = 0.0, rho_peel = 0.0, rho_perf = 0.0;
  double theta_prob_lemma_seq = 0.0, theta_peel = 0.0;

  std::random_device rd;
  std::mt19937 gen(rd());
  std::bernoulli_distribution rv_select(0.5);

  for (uint kk = 0; kk < no_runs; ++kk) {
    double rho_prob_lemma_seq_tt = 0.0, rho_prob_lemma_seq_temp = 0.0;
    double theta_prob_lemma_seq_tt = 0.0;
    uint64_t count_prob_lemma_seq = 0;
    double nC2 = ((double)(G_no_nodes - 1) * G_no_nodes) / 2.0;

    // const bool self_loops_allowed = 1;
    // G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
    // pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
    // // ***SAVE GRAPH
    // const char *FName = "PA_graph.dat";
    // const char *Desc = "PA Graph generated via Jithin's technique";
    // TSnap::SaveEdgeList(G, FName, Desc);
    // *** LOAD GRAPH
    // const char *FName = "PA_graph_n50_m3.dat";
    // const char *FName = "PA_graph_n10_m2.dat";
    const char *FName = "PA_graph_n5_m2.dat";
    // const char *FName = "PA_graph_n50_m3_good.dat";
    G = TSnap::LoadEdgeList<PNEANet>(FName);

    G_no_nodes = (uint64_t)G->GetNodes();
    TUInt64UInt64H rank_orig;
    for (uint64_t i = 0; i < G_no_nodes; i++) {
      rank_orig.AddDat(i) = i;
    }
    std::vector<std::pair<int, int>> nonperfectpairs;
    std::vector<std::pair<int, int>> perfectpairs;
    std::vector<std::pair<int, int>> samebinpairs;
    rank_parallel = find_parallel_rank(G);
    G_DAG = create_DAG(G, rank_parallel);
    std::tie(rho_perf_tt, count_perfect_pairs_nlzd, rho_peel_tt, theta_peel_tt,
             nonperfectpairs, perfectpairs, samebinpairs) =
        FindRhoCountPerfectPeel_int(G_DAG, rank_parallel, rank_orig,
                                    G_no_nodes);
    rank_pair = find_seq_rank_unif_mappair(G);
    std::set<int> interchangable_pairs;
    uint64_t perfectpairs_len = perfectpairs.size();
    uint64_t nonperfectpairs_len = nonperfectpairs.size();
    uint64_t samebinpairs_len = samebinpairs.size();
    //  print_map(rank_pair,"rank_pair");
    //  print_vector(perfectpairs,"perfectpairs");
    //  print_vector(nonperfectpairs,"nonperfectpairs");
    //  print_vector(samebinpairs,"samebinpairs");
    for (auto &x : rank_pair) {
      if (std::find(perfectpairs.begin(), perfectpairs.end(), x.second) ==
          perfectpairs.end()) {
        interchangable_pairs.insert(x.first);
      }
    }
    std::vector<double> emp_average(nonperfectpairs_len, 0.0);
    std::vector<double> emp_average_samebin(samebinpairs_len, 0.0);

    // USED DYER Method
    for (uint64_t kkkk = 0; kkkk < no_runs_MC; kkkk++) {
      cout << "run MC:" << kkkk << endl;
      std::uniform_int_distribution<int> random_pair(
          0, (G_no_nodes - 2)); // guaranteed unbiased
      uint target_pair = random_pair(gen);
      cout << "target_pair: " << target_pair << endl;
      uint rv_select_instance = rv_select(gen);
      cout << "rv_select_instance: " << rv_select_instance << endl;
      if (rv_select_instance &&
          (std::find(interchangable_pairs.begin(), interchangable_pairs.end(),
                     target_pair) != interchangable_pairs.end())) {
        // uint sel_element = *select_randomly(interchangable_pairs.begin(),
        // interchangable_pairs.end());
        uint sel_element = target_pair;
        int temp_1, temp_2, temp_11, temp_22;
        std::pair<int, int> new_next_sel_element, new_prev_sel_element,
            new_sel_element;
        std::tie(temp_1, temp_2) = rank_pair[sel_element];
        std::swap(temp_1, temp_2);
        new_sel_element = std::make_pair(temp_1, temp_2);
        rank_pair[sel_element] = new_sel_element;

        if (std::find(perfectpairs.begin(), perfectpairs.end(),
                      new_sel_element) != perfectpairs.end()) {
          interchangable_pairs.erase(sel_element);
        }
        if (sel_element > 0) {
          std::tie(temp_11, temp_22) = rank_pair[sel_element - 1];
          new_prev_sel_element = std::make_pair(temp_2, temp_22);
          rank_pair[sel_element - 1] = new_prev_sel_element;
          if (std::find(perfectpairs.begin(), perfectpairs.end(),
                        new_prev_sel_element) != perfectpairs.end()) {
            interchangable_pairs.erase((sel_element - 1));
          } else {
            interchangable_pairs.insert((sel_element - 1));
          }
        }
        if (sel_element < (G_no_nodes - 2)) {
          std::tie(temp_11, temp_22) = rank_pair[sel_element + 1];
          new_next_sel_element = std::make_pair(temp_11, temp_1);
          rank_pair[sel_element + 1] = new_next_sel_element;
          if (std::find(perfectpairs.begin(), perfectpairs.end(),
                        new_next_sel_element) != perfectpairs.end()) {
            interchangable_pairs.erase((sel_element + 1));
          } else {
            interchangable_pairs.insert((sel_element + 1));
          }
        }
        cout << "sel_element: " << sel_element << endl;
      }

      // std::vector<std::pair<int,int> > rank_pair_second;
      // for (auto &s : rank_pair)
      //      rank_pair_second.push_back(s.second);

      std::vector<int> present_linear_extension;
      for (auto it = rank_pair.begin(); it != rank_pair.end(); ++it) {
        if (it == rank_pair.begin()) {
          present_linear_extension.push_back(((*it).second).second);
          present_linear_extension.push_back(((*it).second).first);
        } else {
          present_linear_extension.push_back(((*it).second).first);
        }
      }

      // print_map(rank_pair,"rank_pair");
      // print_set(interchangable_pairs,"interchangable_pairs");
      // print_vector(rank_pair_second,"rank_pair_second");
      // print_vector(present_linear_extension, "linear_extenstion");

      // uint64_t deg_present_node = interchangable_pairs.size();
      uint mixing_time = 5;
      if (kkkk > mixing_time) {
        for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
          auto pos_newer =
              std::distance(present_linear_extension.begin(),
                            std::find(present_linear_extension.begin(),
                                      present_linear_extension.end(),
                                      nonperfectpairs[jj].first));
          auto pos_older =
              std::distance(present_linear_extension.begin(),
                            std::find(present_linear_extension.begin(),
                                      present_linear_extension.end(),
                                      nonperfectpairs[jj].second));
          // cout<<"pos_newer: "<<pos_newer<<" value:
          // "<<nonperfectpairs[jj].first<<endl; cout<<"pos_older:
          // "<<pos_older<<" value: "<<nonperfectpairs[jj].second<<endl;

          if (pos_newer > pos_older) {
            // cout<<"Entered: nonperfectpairs
            // "<<nonperfectpairs[jj].first<<","<<nonperfectpairs[jj].second<<endl;
            // emp_average[jj] += 1.0/deg_present_node;
            emp_average[jj] += 1.0 / (no_runs_MC - mixing_time);
          }
        }
        for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
          auto pos_newer =
              std::distance(present_linear_extension.begin(),
                            std::find(present_linear_extension.begin(),
                                      present_linear_extension.end(),
                                      samebinpairs[jj].first));
          auto pos_older =
              std::distance(present_linear_extension.begin(),
                            std::find(present_linear_extension.begin(),
                                      present_linear_extension.end(),
                                      samebinpairs[jj].second));
          // cout<<"pos_newer: "<<pos_newer<<" value:
          // "<<samebinpairs[jj].first<<endl; cout<<"pos_older: "<<pos_older<<"
          // value: "<<samebinpairs[jj].second<<endl;

          if (pos_newer > pos_older) {
            // cout<<"Entered: nonperfectpairs
            // "<<nonperfectpairs[jj].first<<","<<nonperfectpairs[jj].second<<endl;
            // emp_average[jj] += 1.0/deg_present_node;
            emp_average_samebin[jj] += 1.0 / (no_runs_MC - mixing_time);
          }
        }
        present_linear_extension.clear();
      }
    }
    uint64_t node_u, node_v;
    for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
      cout << "(" << nonperfectpairs[jj].first << ","
           << nonperfectpairs[jj].second << "): " << emp_average[jj] << endl;
      if (emp_average[jj] >= 0.8) {
        count_prob_lemma_seq += 1;
        node_u = nonperfectpairs[jj].first;
        node_v = nonperfectpairs[jj].second;
        if (rank_orig.GetDat(node_u) > rank_orig.GetDat(node_v))
          rho_prob_lemma_seq_temp += 1.0 / nC2;
      }
      // This will change the Peeling bins, and may not give consistent linear
      // order
      if (emp_average[jj] <= 0.2) {
        count_prob_lemma_seq += 1;
        node_u = nonperfectpairs[jj].first;
        node_v = nonperfectpairs[jj].second;
        if (rank_orig.GetDat(node_u) < rank_orig.GetDat(node_v))
          rho_prob_lemma_seq_temp += 1.0 / nC2;
      }
    }

    // SAME BIN PAIRS
    for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
      cout << "(" << samebinpairs[jj].first << "," << samebinpairs[jj].second
           << "): " << emp_average_samebin[jj] << endl;
      if (emp_average_samebin[jj] >= 0.8) {
        count_prob_lemma_seq += 1;
        node_u = samebinpairs[jj].first;
        node_v = samebinpairs[jj].second;
        if (rank_orig.GetDat(node_u) > rank_orig.GetDat(node_v))
          rho_prob_lemma_seq_temp += 1.0 / nC2;
      }
      if (emp_average_samebin[jj] <= 0.2) {
        count_prob_lemma_seq += 1;
        node_u = samebinpairs[jj].first;
        node_v = samebinpairs[jj].second;
        if (rank_orig.GetDat(node_u) < rank_orig.GetDat(node_v))
          rho_prob_lemma_seq_temp += 1.0 / nC2;
      }
    }

    cout << "rho_prob_lemma_seq_temp: " << rho_prob_lemma_seq_temp << endl;
    rho_prob_lemma_seq_tt = rho_prob_lemma_seq_temp + rho_perf_tt;
    theta_prob_lemma_seq_tt =
        (rho_prob_lemma_seq_temp + rho_perf_tt) /
        (((double)count_prob_lemma_seq / nC2) + count_perfect_pairs_nlzd);
    rho_peel += rho_peel_tt / no_runs;
    rho_perf += rho_perf_tt / no_runs;
    rho_prob_lemma_seq += rho_prob_lemma_seq_tt / no_runs;
    theta_peel += theta_peel_tt / no_runs;
    theta_prob_lemma_seq += theta_prob_lemma_seq_tt / no_runs;

    // Writing results to a file
    const char *outfile = "P_matrix_RW_KK.txt";
    std::ofstream ost{outfile};
    for (uint64_t jj = 0; jj < nonperfectpairs_len; jj++) {
      node_u = nonperfectpairs[jj].first;
      node_v = nonperfectpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t" << emp_average[jj] << "\n";
      ost << node_v << "\t" << node_u << "\t" << 1 - emp_average[jj] << "\n";
    }
    for (uint64_t jj = 0; jj < samebinpairs_len; jj++) {
      node_u = samebinpairs[jj].first;
      node_v = samebinpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t" << emp_average_samebin[jj]
          << "\n";
      ost << node_v << "\t" << node_u << "\t" << 1 - emp_average_samebin[jj]
          << "\n";
    }
    for (uint64_t jj = 0; jj < perfectpairs_len; jj++) {
      node_u = perfectpairs[jj].first;
      node_v = perfectpairs[jj].second;
      ost << node_u << "\t" << node_v << "\t"
          << "1"
          << "\n";
      ost << node_v << "\t" << node_u << "\t"
          << "0"
          << "\n";
    }
    ost.close();
    cout << "Written P matrix to file" << endl;
    //-------------------------
    // theta and density into file:
    ost1 << theta_prob_lemma_seq_tt << "\t"
         << rho_prob_lemma_seq_tt / theta_prob_lemma_seq_tt << "\n";
  }
  ost1.close();
  printf("rho(Rperf): %f, theta(Rperf): %f \n", rho_perf_tt, 1.0);
  printf("rho(Rpeel): %f, theta(Rpeel): %f \n", rho_peel_tt, theta_peel_tt);
  printf("Density(Rpeel): %f \n", rho_peel_tt / theta_peel_tt);
  printf("rho(Rguess): %f, theta(Rguess): %f \n", rho_prob_lemma_seq,
         theta_prob_lemma_seq);
  printf("Density(Rguess): %f \n", rho_prob_lemma_seq / theta_prob_lemma_seq);
}

int main(int argc, char *argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("cpm. build: %s, %s. Time: %s", __TIME__, __DATE__,
                         TExeTm::GetCurTm()));
  TExeTm ExeTm;

  const int TimeN = Env.GetIfArgPrefixInt("-timen:", 100, "time_steps");
  const float pr_alpha =
      Env.GetIfArgPrefixFlt("-pralpha:", 1, "add node or not");
  const int vec_p_1 =
      Env.GetIfArgPrefixInt("-vecp1:", 5, "add node: lb for uniform dist of m");
  const int vec_p_2 =
      Env.GetIfArgPrefixInt("-vecp2:", 5, "add node: ub for uniform dist of m");
  const float pr_beta = Env.GetIfArgPrefixFlt(
      "-prbeta:", 1, "add node: preferentially or uniformly");
  const int vec_q_1 = Env.GetIfArgPrefixInt(
      "-vecq1:", 1, "edge b/w existing nodes: lb for uniform dist of m");
  const int vec_q_2 = Env.GetIfArgPrefixInt(
      "-vecq2:", 1, "edge b/w existing nodes: ub for uniform dist of m");
  const float pr_delta = Env.GetIfArgPrefixFlt(
      "-prdelta:", 1,
      "edge b/w existing nodes: source node, preferentially or uniformly");
  const float pr_gamma = Env.GetIfArgPrefixFlt(
      "-prgamma:", 1,
      "edge b/w existing nodes: terminal node, preferentially or uniformly");
  const uint no_runs = Env.GetIfArgPrefixInt("-noruns:", 1, "No. of runs");
  const uint width_box =
      Env.GetIfArgPrefixInt("-widthbox:", 10, "Box width of recall matrix");
  const uint no_runs_seq = Env.GetIfArgPrefixInt(
      "-norunsseq:", 100, "No. of runs for Sequential Importance Sampling");
  const uint no_runs_MC =
      Env.GetIfArgPrefixInt("-norunsMC:", 100, "No. of runs for Sequential");
  const int choice = Env.GetIfArgPrefixInt(
      "-choice:", 0,
      "0:via Sequential technique uniformly sampling m-degree nodes");
  printf("\n \n");
  switch (choice) {
  case 0:
    with_Sequential_uniform(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1,
                            vec_q_2, pr_delta, pr_gamma, no_runs_seq, no_runs);
    break;
  case 1:
    with_importance_sampling(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta,
                             vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs_seq,
                             no_runs);
    break;
  case 2:
    estimate_puv_KK_algorithm(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta,
                              vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs_MC,
                              no_runs);
    break;
  case 3:
    with_puv_jks_algorithm(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1,
                           vec_q_2, pr_delta, pr_gamma, no_runs_MC, no_runs,
                           width_box);
    break;
  case 4:
    with_puv_KK_algorithm(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1,
                          vec_q_2, pr_delta, pr_gamma, no_runs_MC, no_runs);
    break;
  default:
    printf("Can not determine what to execute!\n");
  }
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());
  return 0;
}
