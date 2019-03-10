/* (c) 2019 Jithin K. Sreedharan
Program to estimate the probability, for any given pair (u,v), u is older than
v. This probability is denoted by p_uv.

Available choices:
  0: Using the Markov chain Monte Carlo technique developed in our paper for
counting linear extension
  1: Use the ideas proposed in the following two papers:
    - Bubley, R., & Dyer, M. (1999). Faster random generation of linear
extensions. Discrete mathematics, 201(1-3), 81-88.
    - Karzanov, A., & Khachiyan, L. (1991). On the conductance of order Markov
chains. Order, 8(1), 7-15.

Both the choices require graph generation and below are the parameters of the
generalized preferential attachment graph:
- `timen`: Total time-steps of the procedure
- `pralpha`: With this probability a new node will be added; with probability
`(1-pralpha)`, new edges will be added between existing nodes.
- `vecp1`: Lower end of uniform distribution for `m` (number of edges each new
node brings into the graph) when a new node is added.
- `vecp2`: Upper end of uniform distribution for `m` (number of edges each new
node brings into the graph) when a new node is added.
- `prbeta`: With this probability end points of edges of the new node will be
selected preferentially; with probability `(1-prbeta)` endpoints of edges of new
node will be chosen uniformly at random.
- `vecq1`: Lower end of uniform distribution for `m` (number of edges each new
node brings into the graph) when edges between existing nodes are added;
- `vecq2`: Upper end of uniform distribution for `m` (number of edges each new
node brings into the graph) when edges between existing nodes are added.
- `prdelta`: With this probability,When adding edges between existing nodes,
slource node be selected preferentially.
- `prgamma`: With this probability, when adding edges between existing nodes,
Terminal node will be selected preferentially

Note that even though the programs are for generalized preferential attachment
model, the theory of optimization and estimation is guaranteed only for standard
preferential attachment model.

Sample usage: /times_randomgraph_estimate_puv -timen:50 -pralpha:1 -prbeta:1
-vecp1:3 -vecp2:3 -noruns:100 -norunsMC:400000 -choice:0
 */

#include "times.hpp"

/* Sequential algorithm; removes only one m-degree node at any step, and is
selected uniformly at random.
Output: rank in pairwise format like "rank_id: Node,Node"; e.g.: rank 0: (second
oldest,oldest) */

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
        0, (MnDegV.Len() - 1));  // guaranteed unbiased
    node_sel = MnDegV[min_deg(gen)].Val;
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

/* A function which returns many results together:
  1. Recall of perfect pair estimator
  2. Number of prefect pairs (normalized by nC2)
  3. Recall of Peeling estimator
  4. Precision of Peeling estimator
  5. Vector of non-perfect pairs in the output of Peeling estimator
  6. Vector of perfect pairs
  7. Vector of same-bin pairs in the putput of Peeling estimator */

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
          if (node_u_rank > node_v_rank) beta_ij += 1.0;
          if ((NIdV_temp).IsIn(node_v)) {
            perfectpairs.push_back(std::make_pair(node_u, node_v));
            if (node_u_rank > node_v_rank) alpha_ij += 1.0;
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

/* The following function carry out multiple tasks:
  1. Using MCMC technique proposed in our paper for estimating p_uv and writing
  to a file.
  2. Generate the points of following estimators: Peeling, Peeling+,
  Precision-1. Write them separately to different files
  3. Write recall matrix */

void estimate_puv_jks_algorithm(const int &TimeN, const float &pr_alpha,
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

  std::ofstream ost_peel{outfile_theta_delta_peel};
  std::ofstream ost_peelp{outfile_theta_delta_peelp};
  std::ofstream ost_perf{outfile_theta_delta_perf};

  double rho_perf_tt = 0, count_perfect_pairs_nlzd = 0, rho_peel_tt = 0,
         theta_peel_tt = 0;
  double rho_prob_lemma_seq = 0.0, rho_peel = 0.0, rho_perf = 0.0,
         density_peel = 0, density_guess = 0;
  double theta_prob_lemma_seq = 0.0, theta_peel = 0.0;

  double rho_peelp = 0, theta_peelp = 0, density_peelp = 0;
  const int no_threads = 1;

  /* For finding Recall Matrix we define G_no_nodes here. Note that only the
  classical PA graph model is allowed. Create boxes of equal width_box */
  G_no_nodes = TimeN;
  uint64_t total_boxes = std::ceil(G_no_nodes / (double)width_box);
  std::vector<std::vector<double>> recall_matrix(
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

    // Peel estimator
    rank_parallel = find_parallel_rank(G);
    G_DAG = create_DAG(G, rank_parallel);
    std::tie(rho_perf_tt, count_perfect_pairs_nlzd, rho_peel_tt, theta_peel_tt,
             nonperfectpairs, perfectpairs, samebinpairs) =
        FindRhoCountPerfectPeel_int(G_DAG, rank_parallel, rank_orig,
                                    G_no_nodes);
    rank_pair = find_seq_rank_unif_mappair(G);

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

    // MCMC CALCULATION STARTS HERE
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
      std::vector<int> present_linear_extension;
      for (auto it = rank_pair.begin(); it != rank_pair.end(); ++it) {
        if (it == rank_pair.begin()) {
          present_linear_extension.push_back(((*it).second).second);
          present_linear_extension.push_back(((*it).second).first);
        } else {
          present_linear_extension.push_back(((*it).second).first);
        }
      }
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

        if (pos_newer > pos_older) {
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
        if (pos_newer > pos_older) {
          emp_average_samebin[jj] += 1.0 / deg_present_node;
        }
      }
      present_linear_extension.clear();
      normalization += 1.0 / deg_present_node;
    }
    std::transform(
        emp_average.begin(), emp_average.end(), emp_average.begin(),
        std::bind1st(std::multiplies<double>(), 1.0 / normalization));
    std::transform(
        emp_average_samebin.begin(), emp_average_samebin.end(),
        emp_average_samebin.begin(),
        std::bind1st(std::multiplies<double>(), 1.0 / normalization));
    uint64_t node_u, node_v;

    cout << "rho_prob_lemma_seq_temp: " << rho_prob_lemma_seq_temp << endl;
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

    ost_peel << theta_peel_tt << "\t" << rho_peel_tt / theta_peel_tt << "\n";
    ost_peelp << theta_peelp_temp << "\t" << rho_peelp_temp / theta_peelp_temp
              << "\n";
    ost_perf << "1.0"
             << "\t" << rho_perf_tt << "\n";
  }
  ost_peel.close();
  ost_peelp.close();
  ost_perf.close();
  printf("rho(Rperf): %f, theta(Rperf): %f, density(Rperf): %f \n", rho_perf,
         1.0, rho_perf);
  printf("rho(Rpeel): %f, theta(Rpeel): %f, density(Rpeel): %f \n", rho_peel,
         theta_peel, density_peel);
  printf("rho(Rpeelp): %f, theta(Rpeelp): %f, density(Rpeelp): %f \n",
         rho_peelp, theta_peelp, density_peelp);

  std::string outfile_recall = "./temp_model/Recall_matrix.txt";
  std::ofstream ost_recall{outfile_recall};
  for (int i = 0; i < total_boxes; i++)
    for (int j = 0; j < total_boxes; j++) {
      ost_recall << i << "\t" << j << "\t" << recall_matrix[i][j] << "\n";
    }
  ost_recall.close();
}

// Using K-K technique ro estimate p_uv and write to file
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

    const bool self_loops_allowed = 1;
    G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
                                        pr_beta, vec_q_1, vec_q_2, pr_delta,
                                        pr_gamma, self_loops_allowed);
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

    // USED DYER Method
    for (uint64_t kkkk = 0; kkkk < no_runs_MC; kkkk++) {
      cout << "run MC:" << kkkk << endl;
      std::uniform_int_distribution<int> random_pair(
          0, (G_no_nodes - 2));  // guaranteed unbiased
      uint target_pair = random_pair(gen);
      cout << "target_pair: " << target_pair << endl;
      uint rv_select_instance = rv_select(gen);
      cout << "rv_select_instance: " << rv_select_instance << endl;
      if (rv_select_instance &&
          (std::find(interchangable_pairs.begin(), interchangable_pairs.end(),
                     target_pair) != interchangable_pairs.end())) {
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

      std::vector<int> present_linear_extension;
      for (auto it = rank_pair.begin(); it != rank_pair.end(); ++it) {
        if (it == rank_pair.begin()) {
          present_linear_extension.push_back(((*it).second).second);
          present_linear_extension.push_back(((*it).second).first);
        } else {
          present_linear_extension.push_back(((*it).second).first);
        }
      }

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

          if (pos_newer > pos_older) {
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
          if (pos_newer > pos_older) {
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
  const int choice = Env.GetIfArgPrefixInt(
      "-choice:", 0,
      "0: estimate p_uv - our algorithm; 1: estimate p_uv - K-K algorithm ");
  const uint width_box =
      Env.GetIfArgPrefixInt("-widthbox:", 10, "Box width of recall matrix");
  const uint no_runs_MC =
      Env.GetIfArgPrefixInt("-norunsMC:", 100, "No. of runs for Sequential");

  printf("\n \n");
  switch (choice) {
    case 0:
      /* Estimate probabilities p_uv and write to a file. Use our Markov chain
       * Monte Carlo technique for counting the linear extension. */
      estimate_puv_jks_algorithm(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta,
                                 vec_q_1, vec_q_2, pr_delta, pr_gamma,
                                 no_runs_MC, no_runs, width_box);
      break;
    case 1:
      /* Estimate probabilities p_uv and write to a file. Use the method
       * proposed
       * by Karzanov and Khachiyan, and the modification by Bubley and Dyer. */
      estimate_puv_KK_algorithm(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta,
                                vec_q_1, vec_q_2, pr_delta, pr_gamma,
                                no_runs_MC, no_runs);
      break;
    default:
      printf("Can not determine what to execute!\n");
  }
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());
  return 0;
}
