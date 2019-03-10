/* (c) 2019 Jithin K. Sreedharan
Program to find the node ranking based on their age.

Available choices:
  0: Find rank using Peeling algorithm, and write to a file
  1: Find rank using Peeling & Peeling+ algorithms, and calculate precision and
recall

Sample usage: ./times_givengraph i:"../data/dynamic-simplewiki.txt"
-ior:"../data/dynamic-simplewiki_data.csv" -choice:1
*/
#include "times.hpp"

void calculate_rho_theta_density_peeling_n_peelingplus_STL(
    const TStr InFNm, const char *InFNm_or, const int no_threads) {
  PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm, 0, 1);

  Tuint64uint64VM rank_parallel;
  Tuint64uint64M rank_node;
  Tuint64uint64M rank_orig;

  uint64_t G_no_nodes;

  uint64_t rho_unnormalized, nopairsDAG;
  uint64_t rho_extra_unnnormalized, nopairsDAGextra;
  uint64_t no_orig_pairs_same_bin = 0, no_orig_pairs_across_bin = 0,
           no_orig_pairs = 0;

  double theta = 0, rho = 0;
  double rho_extra = 0, theta_extra = 0;
  uint64_t depth_parallel = 0;

  uint64_t kk, vv;  // CAREFUL HERE, we use integer, it may be string as well
  std::ifstream infile(InFNm_or);
  while (infile >> kk >> vv) {
    rank_orig.insert({kk, vv});
  }

  G_no_nodes = (uint64_t)G->GetNodes();

  // Structure of rank_parallel = {rank:list of nodes},rank_node = {node:rank}
  printf("Peeling ranking: started \n");
  std::tie(rank_parallel, rank_node) = find_parallel_rank_returnnoderank_STL(G);
  printf("Peeling ranking: finished \n");

  double time_temp = omp_get_wtime();
  printf("Metric calculation: Peeling - started \n");
  std::tie(rho_unnormalized, nopairsDAG) = find_prec_recall_real_networks_STL(
      rank_parallel, rank_orig, G_no_nodes, no_threads);
  no_orig_pairs_across_bin = nopairsDAG;
  printf("Metric calculation: Peeling - ended \n");
  time_temp = omp_get_wtime() - time_temp;
  printf("Time taken for parallel implementation: %.2fs\n", time_temp);
  std::cout << "rho_unnormalized: " << rho_unnormalized
            << " nopairsDAG: " << nopairsDAG << " G_no_nodes: " << G_no_nodes
            << std::endl;

  time_temp = omp_get_wtime();
  printf("Peeling with same bin pairs ranking and metric started \n");
  std::tie(rho_extra_unnnormalized, nopairsDAGextra, no_orig_pairs_same_bin) =
      Count_extra_pairs_real_networks_STL(G, G_no_nodes, rank_orig,
                                          rank_parallel, rank_node, no_threads);
  printf("Peeling with same bin pairs ranking and metric ended \n");
  time_temp = omp_get_wtime() - time_temp;
  printf("Time taken for parallel implementation of Peeling+: %.2fs\n",
         time_temp);

  no_orig_pairs = no_orig_pairs_same_bin + no_orig_pairs_across_bin;
  rho = rho_unnormalized / (double)no_orig_pairs;
  rho_extra =
      (rho_extra_unnnormalized + rho_unnormalized) / (double)no_orig_pairs;

  theta = rho_unnormalized / (double)nopairsDAG;
  theta_extra = (rho_extra_unnnormalized + rho_unnormalized) /
                (double)(nopairsDAG + nopairsDAGextra);

  depth_parallel = rank_parallel.size();

  TStr FuncStr = TStr::Fmt("%s:graph", __func__);
  PrintGStats(FuncStr.CStr(), G);
  std::cout << " rho_peel: " << rho << ", theta_peel: " << theta
            << ", depth_bins: " << depth_parallel
            << ", density_peel: " << rho / theta << std::endl;
  std::cout << " rho_peel+: " << rho_extra << " theta_peel+: " << theta_extra
            << ", density_peel+: " << rho_extra / theta_extra << std::endl;
}

int main(int argc, char *argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("cpm. build: %s, %s. Time: %s", __TIME__, __DATE__,
                         TExeTm::GetCurTm()));
  TExeTm ExeTm;
  const TStr InFNm = Env.GetIfArgPrefixStr(
      "-i:", "./data/cit-HepPh_connected.txt", "Input edgelist file name");
  const TStr InFNm_or_t =
      Env.GetIfArgPrefixStr("-ior:", "./data/cit-HepPh_connected.txt",
                            "Input original rank file name");
  const char *InFNm_or = InFNm_or_t.CStr();
  const int choice = Env.GetIfArgPrefixInt(
      "-choice:", 0,
      "0: Write predicted rank by Peeling algorithm to a file, 1: Find "
      "predicted rank using Peeling & Peeling+ algorithms and utput "
      "performance measures");
  const int no_threads =
      Env.GetIfArgPrefixInt("-nothreads:", 1, "No of threads (default = 1)");
  cout << "\n";

  switch (choice) {
    case 0: {
      // Write ranking to a file
      TUInt64UInt64VH rank_parallel;
      TUInt64UInt64H rank_node;
      PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm, 0, 1);
      TStr FuncStr = TStr::Fmt("%s:graph", __func__);
      PrintGStats(FuncStr.CStr(), G);
      cout << "\n";
      std::tie(rank_parallel, rank_node) = find_parallel_rank_returnnoderank(G);
      const char *outfile = "node_rank.txt";
      std::ofstream ost{outfile};

      for (TUInt64UInt64H::TIter It = rank_node.BegI(); It < rank_node.EndI();
           It++) {
        // get the key
        TUInt64 Key = It.GetKey();
        TUInt64 Value = It.GetDat();
        // write
        ost << Key << "\t" << Value << "\n";
      }
      ost.close();
      cout << "Written results to node_rank.txt";
      break;
    }
    case 1:
      /* Find rank using Peeling, Peeling+ algorithms. Calculate precision and
       * recall. STL and parallel implementation */
      calculate_rho_theta_density_peeling_n_peelingplus_STL(InFNm, InFNm_or,
                                                            no_threads);
      break;
    default:
      printf("Invalid choice\n");
  }
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());
  return 0;
}
