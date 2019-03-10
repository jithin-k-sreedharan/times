/* (c) 2019 Jithin K. Sreedharan
Available choices:
        0: probability of success of finding the oldest node in bin 0
        1: find the recall and precision of Peeling, Peeling+ and Perfect pairs
algorithm 2: find maximum and average size of the tree starting from a node in
the youngest bin (top bin). 3: find maximum and average height of the tree
starting from a node in the youngest bin (top bin). 4: find average number of
levels See inline comments for more details.

All the choices require graph generation and below are the parameters of the
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

Sample usage: ./times_randomgraph -timen:5000 -pralpha:0.75 -vecp1:5 -vecp2:50
-prbeta:0.5 -vecq1:5 -vecq2:50  -prdelta:0.5 -prgamma:0.5 -noruns:1000 -choice:0
*/
#include "times.hpp"

void calculate_rho_theta_with_extra_pairs(
    const int& TimeN, const float& pr_alpha, const int& vec_p_1,
    const int& vec_p_2, const float& pr_beta, const int& vec_q_1,
    const int& vec_q_2, const float& pr_delta, const float& pr_gamma,
    const int& no_runs, const float& p) {
  PNEANet G;

  TUInt64UInt64VH rank_parallel;
  TUInt64UInt64H rank_node;
  TUInt64UInt64H rank_orig;
  float rho_temp = 0;
  float theta_temp = 0;
  float rho = 0;
  float theta = 0;
  uint64_t G_no_nodes;
  float depth_parallel = 0;
  float normln_temp;
  uint64_t nopairsDAG;

  float rho_extra_temp = 0;
  float theta_extra_temp = 0;
  float rho_extra = 0;
  float theta_extra = 0;
  uint64_t nopairsDAGextra;

  uint64_t rank_parallel_len;
  const bool self_loops_allowed = 1;

  const char* outfile_theta_delta = "theta_delta.txt";
  std::ofstream ost1{outfile_theta_delta};

  // For Perfect pair calculation
  PNEANet G_DAG;
  TUInt64UInt64VH rank_parallel_1;
  double rho_perf_temp = 0, theta_perf_temp = 0, temp_aaa, temp_bbb,
         rho_perf = 0, theta_perf = 0;

  for (int ii = 0; ii < no_runs; ii++) {
    std::cout << " run index: " << ii << std::endl;

    G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
                                        pr_beta, vec_q_1, vec_q_2, pr_delta,
                                        pr_gamma, self_loops_allowed);

    G_no_nodes = (uint64_t)G->GetNodes();
    for (uint64_t i = 0; i < G_no_nodes; i++) {
      rank_orig.AddDat(i) = i;
    }
    // Structure of rank_parallel = {rank:list of nodes},rank_node = {node:rank}
    std::tie(rank_parallel, rank_node) = find_parallel_rank_returnnoderank(G);

    rho_temp = FindThetaVH_kendaul_tau(rank_parallel, rank_orig, G_no_nodes, p);
    normln_temp = (G_no_nodes * (G_no_nodes - 1) / (float)2);
    nopairsDAG = no_pairs_DAG(rank_parallel);
    rho += rho_temp;
    theta_temp = rho_temp * normln_temp / (float)nopairsDAG;
    theta += theta_temp;

    std::tie(rho_extra_temp, nopairsDAGextra) =
        Count_extra_pairs(G, G_no_nodes, rank_orig, rank_parallel, rank_node);
    rho_extra += rho_temp + rho_extra_temp;
    theta_extra_temp = (rho_temp + rho_extra_temp) * normln_temp /
                       (float)(nopairsDAG + nopairsDAGextra);
    theta_extra += theta_extra_temp;

    // // Perfect-Precision algorithms to find perfect pairs
    // cout<<"entering perf-precision algorithm"<<endl;
    // std::tie(G_DAG,rank_parallel_1) = find_parallel_rank_DAG(G);
    // std::tie(rho_perf_temp,theta_perf_temp, temp_aaa, temp_bbb) =
    // count_perfect_correct_pairs_DAG(G_DAG, rank_parallel, rank_orig,
    // G_no_nodes); rho_perf += rho_perf_temp; theta_perf += theta_perf_temp;
    // cout<<"leaving perf-precision algorithm"<<endl;

    rank_parallel_len = rank_parallel.Len();
    depth_parallel = depth_parallel + rank_parallel_len;
    ost1 << theta_perf_temp << "\t" << (rho_perf_temp / theta_perf_temp)
         << "\n";
  }
  ost1.close();
  rho = rho / no_runs;
  theta = theta / no_runs;
  rho_extra = rho_extra / no_runs;
  theta_extra = theta_extra / no_runs;
  depth_parallel = depth_parallel / no_runs;
  rho_perf /= no_runs;
  theta_perf /= no_runs;

  TStr FuncStr = TStr::Fmt("%s:graph", __func__);
  PrintGStats(FuncStr.CStr(), G);
  std::cout << " rho_peel: " << rho << ", theta_peel: " << theta
            << ", depth_bins: " << depth_parallel
            << ", density_peel: " << (rho / theta) << std::endl;
  std::cout << " rho_peel+: " << rho_extra << " theta_peel+: " << theta_extra
            << ", density_peel+: " << (rho_extra / theta_extra) << std::endl;
  std::cout << " --- Pefect pair algorithm disabled by default ---"
            << std::endl;
  std::cout << " rho_perf " << rho_perf << " theta_perf: " << theta_perf
            << ", density_perf: " << (rho_perf / theta_perf) << std::endl;
}

void pr_success_node_0_in_bin_0(const int& TimeN, const float& pr_alpha,
                                const int& vec_p_1, const int& vec_p_2,
                                const float& pr_beta, const int& vec_q_1,
                                const int& vec_q_2, const float& pr_delta,
                                const float& pr_gamma, const int& no_runs) {
  PNEANet G;
  TUInt64UInt64VH rank_parallel;
  TUInt64V bin_0;
  uint64_t theta_par_pref = 0;
  uint64_t bin_0_len = 0;
  const bool self_loops_allowed = 0;
  for (int ii = 0; ii < no_runs; ii++) {
    std::cout << " run index: " << ii << std::endl;

    G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
                                        pr_beta, vec_q_1, vec_q_2, pr_delta,
                                        pr_gamma, self_loops_allowed);

    rank_parallel = find_parallel_rank(G);
    // Oldest node bin
    bin_0 = rank_parallel.GetDat(0);

    if (bin_0.IsIn(0)) {
      theta_par_pref += 1;
      bin_0_len += bin_0.Len();
      // print_TUInt64V(bin_0,"bin_0");
    }
  }
  double prob_success = (double)theta_par_pref / no_runs;
  double avg_bin_0_len = (double)bin_0_len / theta_par_pref;

  printf("Prob of success: %0.3f; Avg size of bin 0: %0.3f\n", prob_success,
         avg_bin_0_len);
}

// Study the size of the directed tree rooted at nodes in the youngest bin
std::pair<double, double> max_avg_tree_size_DAG(
    const int& TimeN, const float& pr_alpha, const int& vec_p_1,
    const int& vec_p_2, const float& pr_beta, const int& vec_q_1,
    const int& vec_q_2, const float& pr_delta, const float& pr_gamma,
    const int& no_runs) {
  PNEANet G;

  TUInt64UInt64H rank_seq_uni;
  TUInt64UInt64VH rank_parallel;
  TUInt64UInt64VH rank_parallel_new;
  TUInt64V rank_parallel_temp;

  TUInt64UInt64H rank_orig;
  uint64_t G_no_nodes;
  TUInt64V bb;
  const bool self_loops_allowed = 1;

  double avg_tr = 0;
  double max_tr = 0;
  double avg_tt, max_tt;
  PNEANet G_DAG;
  for (int ii = 0; ii < no_runs; ii++) {
    std::cout << " run index: " << ii << std::endl;
    G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
                                        pr_beta, vec_q_1, vec_q_2, pr_delta,
                                        pr_gamma, self_loops_allowed);

    G_no_nodes = (uint64_t)G->GetNodes();
    for (uint64_t i = 0; i < G_no_nodes; i++) {
      rank_orig.AddDat(i) = i;
    }
    std::tie(G_DAG, rank_parallel) = find_parallel_rank_DAG(G);
    write_dot(G, "original_G.dot.dat", "G");
    write_dot(G_DAG, "DAG_G.dot.dat", "DAG(G)");

    std::tie(avg_tt, max_tt) =
        count_tree_size_DAG(G_DAG, rank_parallel, G_no_nodes);

    avg_tr = avg_tr + avg_tt;
    max_tr = max_tr + max_tt;
  }
  avg_tr = avg_tr / no_runs;
  max_tr = max_tr / no_runs;

  return std::make_pair(avg_tr, max_tr);
}

// Study the height of the directed tree rooted at nodes in the youngest bin
std::pair<double, double> max_avg_tree_height_DAG(
    const int& TimeN, const float& pr_alpha, const int& vec_p_1,
    const int& vec_p_2, const float& pr_beta, const int& vec_q_1,
    const int& vec_q_2, const float& pr_delta, const float& pr_gamma,
    const int& no_runs) {
  PNEANet G;
  PNEANet G_DAG;
  TUInt64UInt64VH rank_parallel;
  uint64_t G_no_nodes;
  const bool self_loops_allowed = 1;
  double avg_tr = 0;
  double max_tr = 0;
  double avg_tt, max_tt;
  for (int ii = 0; ii < no_runs; ii++) {
    std::cout << " run index: " << ii << std::endl;
    G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
                                        pr_beta, vec_q_1, vec_q_2, pr_delta,
                                        pr_gamma, self_loops_allowed);
    G_no_nodes = (uint64_t)G->GetNodes();
    std::tie(G_DAG, rank_parallel) = find_parallel_rank_DAG(G);
    std::tie(avg_tt, max_tt) =
        count_tree_height_DAG(G_DAG, rank_parallel, G_no_nodes);
    avg_tr = avg_tr + avg_tt;
    max_tr = max_tr + max_tt;
  }
  avg_tr = avg_tr / no_runs;
  max_tr = max_tr / no_runs;
  return std::make_pair(avg_tr, max_tr);
}

// Study number of bins/levels from the Peeling algorithm
double study_levels(const int& TimeN, const float& pr_alpha, const int& vec_p_1,
                    const int& vec_p_2, const float& pr_beta,
                    const int& vec_q_1, const int& vec_q_2,
                    const float& pr_delta, const float& pr_gamma,
                    const int& no_runs) {
  PNEANet G;
  TUInt64UInt64VH rank_parallel;
  double depth_parallel = 0;
  uint64_t rank_parallel_len;
  const bool self_loops_allowed = 1;
  for (int ii = 0; ii < no_runs; ii++) {
    std::cout << " run index: " << ii << std::endl;

    G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2,
                                        pr_beta, vec_q_1, vec_q_2, pr_delta,
                                        pr_gamma, self_loops_allowed);
    rank_parallel = find_parallel_rank(G);
    rank_parallel_len = rank_parallel.Len();
    depth_parallel = depth_parallel + rank_parallel_len;
  }
  depth_parallel = depth_parallel / no_runs;
  return depth_parallel;
}

int main(int argc, char* argv[]) {
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
  const int no_runs = Env.GetIfArgPrefixInt("-noruns:", 5000, "No. of runs");
  const float p = Env.GetIfArgPrefixFlt(
      "-p:", 0, "Parameter for partial Kendaul-Tau distance");
  const int choice = Env.GetIfArgPrefixInt(
      "-choice:", 2,
      "0: Calculate precision and recall for Peeling and Peeling+"
      "1: probability of finding oldest node in bin 0, 2: find maximum and "
      "avergae tree size, "
      "3: find maxim and avergae tree height"
      "4: find average number of bins");
  const bool write2file = Env.GetIfArgPrefixInt(
      "-write2file:", false,
      "false: display results, true :write results for various set "
      "of parameters to a file");
  printf("\n \n");
  double avg_tr, max_tr;

  switch (choice) {
    case 0:
      /* 1. Compute recall and precision using Peeling method.
       * 2. Compute recall and precision using Peeling+ method (differentiate
       * nodes in the same bin too!)
       * 3. Compute recall and precision of pairs that can be orderd with
       * probability 1*/

      calculate_rho_theta_with_extra_pairs(
          TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2,
          pr_delta, pr_gamma, no_runs, p);
      break;

    case 1:
      // Calculate the probability of finding oldest node in the bin 0 of the
      // Peeling method

      pr_success_node_0_in_bin_0(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta,
                                 vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs);
      break;

    case 2:
      /* Find the average and maximum size of the tree starting from a node in
       * the youngest bin (top bin). The size of the tree is defined as the
       * number of nodes in the tree. */

      if (write2file == 0) {
        std::tie(avg_tr, max_tr) = max_avg_tree_size_DAG(
            TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2,
            pr_delta, pr_gamma, no_runs);
        std::cout << " avg_tr: " << avg_tr << ", max_tr: " << max_tr
                  << std::endl;
      } else {
        std::vector<int> TimeN_v{100, 500, 1000, 2000, 3000, 4000, 5000};
        std::vector<int> vec_p_1_v{5, 10, 15, 20, 25};
        const char* outfile = "results_rg_max_avg_tree.csv";
        std::ofstream ost{outfile};
        ost << "n"
            << ","
            << "m"
            << ","
            << "avg_tr"
            << ","
            << "max_tr"
            << "\n";
        for (auto& i : vec_p_1_v) {
          for (auto& j : TimeN_v) {
            std::tie(avg_tr, max_tr) =
                max_avg_tree_size_DAG(j, pr_alpha, i, i, pr_beta, vec_q_1,
                                      vec_q_2, pr_delta, pr_gamma, no_runs);
            ost << j << "," << i << "," << avg_tr << "," << max_tr << "\n";
          }
        }
        ost.close();
      }
      break;

    case 3:
      /* find average and maximum height of a tree starting from any node in the
       * youngest bin (top bin). The height is defined as the length of the
       * longest path between a node in the youngest bin to any other reachable
       * node */

      if (write2file == 0) {
        std::tie(avg_tr, max_tr) = max_avg_tree_height_DAG(
            TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2,
            pr_delta, pr_gamma, no_runs);
        std::cout << " avg_tr: " << avg_tr << ", max_tr: " << max_tr
                  << std::endl;
      }
      break;

    case 4:
      // Study number of bins/levels from the Peeling algorithm

      double depth_parallel;
      if (write2file == 0) {
        depth_parallel =
            study_levels(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1,
                         vec_q_2, pr_delta, pr_gamma, no_runs);
        std::cout << "depth_bins: " << depth_parallel << std::endl;
      } else {
        std::vector<int> TimeN_v{1000,  2000,  5000,   10000,  25000,
                                 50000, 75000, 100000, 500000, 1000000};
        std::vector<int> vec_p_1_v{5, 25, 50, 100};
        const char* outfile = "results_level_size.csv";
        std::ofstream ost{outfile};
        ost << "n"
            << ","
            << "m"
            << ","
            << "depth_parallel"
            << "\n";
        for (auto& i : vec_p_1_v) {
          for (auto& j : TimeN_v) {
            std::cout << "(n,m): " << j << "," << i << std::endl;
            depth_parallel = study_levels(j, pr_alpha, i, i, pr_beta, vec_q_1,
                                          vec_q_2, pr_delta, pr_gamma, no_runs);
            ost << j << "," << i << "," << depth_parallel << "\n";
          }
        }
        ost.close();
      }
      break;

    default:
      printf("Can not determine what to execute!\n");
  }

  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());
  return 0;
}
