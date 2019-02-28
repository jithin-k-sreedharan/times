
#include "vertex_ordering.hpp"

void certain_correct_pairs_DAG(const TStr InFNm, const char *InFNm_or,
                               const int &no_runs) {
  // PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm,0,1);
  // TUInt64UInt64H rank_seq_uni;
  // TUInt64UInt64VH rank_parallel;
  // TUInt64UInt64VH rank_parallel_new;
  // TUInt64V rank_parallel_temp;
  //
  // TUInt64UInt64H rank_orig;
  // uint64_t G_no_nodes;
  // float depth_parallel = 0;
  // TUInt64V bb;
  // uint64_t rank_parallel_len;
  // const bool self_loops_allowed = 1;
  //
  // //Counting the certain and deductible pairs
  // std::pair<float,float> eta;
  // float eta_cert = 0;
  // float eta_ded = 0;
  // PNEANet G_DAG;
  //
  // uint64_t kk, vv; //CAREFUL HERE, we use integer, it may be string as well
  // // std::ifstream
  // infile("./data/BioGrid_full_human/BioGRID_full_human_net_orig_rank.csv");
  // std::ifstream infile(InFNm_or);
  // while (infile >> kk >> vv) {
  //      rank_orig.AddDat(kk) = vv;
  // }
  //
  // G_no_nodes = (uint64_t)G->GetNodes();
  //
  // printf("Reached 0\n");
  // std::tie(G_DAG,rank_parallel) = find_parallel_rank_DAG(G);
  // printf("Reached 00\n");
  //
  // eta = count_perfect_correct_pairs_DAG(G_DAG, rank_parallel, rank_orig,
  // G_no_nodes); eta_cert = eta_cert + eta.first; eta_ded = eta_ded +
  // eta.second;
  //
  // rank_parallel_len = rank_parallel.Len();
  // depth_parallel = depth_parallel + rank_parallel_len;
  //
  //
  // // eta_cert = eta_cert/no_runs;
  // // eta_ded = eta_ded/no_runs;
  // // depth_parallel = depth_parallel/no_runs;
  //
  // TStr FuncStr = TStr::Fmt("%s:graph", __func__);
  // PrintGStats(FuncStr.CStr(), G);
  // std::cout<<" eta_cert: "<<eta_cert<<", eta_correct: "<<eta_ded<<std::endl;
  // std::cout<<" depth_bins: "<< depth_parallel <<std::endl;
}

void calculate_distance_metric(const TStr InFNm, const char *InFNm_or,
                               const int &no_runs, const float &p,
                               const int &distance_metric) {
  PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm, 0, 1);

  TUInt64UInt64H rank_seq_uni;
  TUInt64UInt64VH rank_parallel;
  TUInt64UInt64VH rank_parallel_new;
  TUInt64V rank_parallel_temp;

  TUInt64UInt64H rank_orig;
  double eta_seq = 0;
  double eta_par = 0;

  uint64_t G_no_nodes = G->GetNodes();
  double eta_t_seq = 0; // redeclared for Kendaul-Tau distance
  double eta_t_par = 0; // redeclared for Kendaul-Tau distance

  uint64_t kk, vv; // CAREFUL HERE, we use integer, it may be string as well
  // std::ifstream
  // infile("./data/BioGrid_full_human/BioGRID_full_human_net_orig_rank.csv");
  std::ifstream infile(InFNm_or);
  while (infile >> kk >> vv) {
    rank_orig.AddDat(kk) = vv;
  }

  double depth_parallel = 0;
  int i_max = 0;

  double normln_temp;
  uint64_t no_pairs;
  double eta_t_par_presn = 0;
  double eta_par_presn = 0;
  for (int ii = 0; ii < no_runs; ii++) {
    std::cout << " run index: " << ii << std::endl;

    // printf("Sequential ranking: started \n");
    // rank_seq_uni = find_seq_rank_unif_Hout(G);
    // printf("Sequential ranking: finished \n");
    printf("Peeling ranking: started \n");
    rank_parallel = find_parallel_rank(G);
    printf("Peeling ranking: finished \n");

    const bool OR_BINS_CONSIDER = 0;

    // PEELING TECHNIQUE: IF ORIGINAL RANK HAS BINS
    //======================================================================
    if (OR_BINS_CONSIDER == 1) {
      // rank_parallel_new = find_parallel_rank_minmax_deg(G,(int)1, (int)133);
      // COMBINE BINS IN RANK_NEW TO MAKE THE NUMBER OF BINS SMALL

      for (int j = 0; j < 16; j++) {
        rank_parallel_temp.Clr();
        if (j == 15) {
          i_max = 247;
        } else {
          i_max = (j * 15) + 15;
        }
        for (int i = (j * 15); i < i_max; i++) {
          rank_parallel_temp.AddV(rank_parallel.GetDat(i));
        }
        rank_parallel_new.AddDat(j, rank_parallel_temp);
      }
      // for (int j = 0; j < 10; j++)
      // {
      //     rank_parallel_temp.Clr();
      //     if (j==9)
      //     {
      //         i_max = 247;
      //     }
      //     else
      //     {
      //         i_max = (j*24)+24;
      //     }
      //     for (int i = (j*24); i < i_max; i++)
      //     {
      //         rank_parallel_temp.AddV(rank_parallel.GetDat(i));
      //     }
      //     rank_parallel_new.AddDat(j,rank_parallel_temp);
      // }
      std::cout << "Length of new rank_parallel: " << rank_parallel_new.Len()
                << std::endl;
    }
    //======================================================================

    if (distance_metric == 1) {
      // printf("Eta: Sequential - started \n");
      // eta_t_seq = FindThetaH_kendaul_tau(rank_seq_uni,
      // rank_orig,G_no_nodes,p); printf("Eta: Sequential - finished \n");
      printf("Eta: Peeling - started \n");
      if (OR_BINS_CONSIDER == 1) {
        eta_t_par = FindThetaVH_kendaul_tau(rank_parallel_new, rank_orig,
                                            G_no_nodes, p);
      } else {
        eta_t_par =
            FindThetaVH_kendaul_tau(rank_parallel, rank_orig, G_no_nodes, p);
      }
      printf("Eta: Peeling - finished \n");
      normln_temp = (G_no_nodes * (G_no_nodes - 1) / (float)2);
      no_pairs = no_pairs_DAG(rank_parallel);
      eta_t_par_presn = eta_t_par * normln_temp / (float)no_pairs;
    } else {
      printf("Eta: Sequential - started \n");
      eta_t_seq = FindThetaH(rank_seq_uni, rank_orig, G_no_nodes);
      printf("Eta: Sequential - finished \n");
      printf("Eta: Peeling - started \n");
      if (OR_BINS_CONSIDER == 1) {
        eta_t_par = FindThetaVH(rank_parallel_new, rank_orig, G_no_nodes);
      } else {
        eta_t_par = FindThetaVH(rank_parallel, rank_orig, G_no_nodes);
      }
      printf("Eta: Peeling - finished \n");
    }
    eta_seq = eta_seq + eta_t_seq;
    eta_par = eta_par + eta_t_par;
    eta_par_presn = eta_par_presn + eta_t_par_presn;
    depth_parallel = depth_parallel + rank_parallel.Len();
  }

  // // get elements by an iterator
  // for (TUInt64UInt64VH::TIter It = rank_parallel.BegI(); It <
  // rank_parallel.EndI(); It++) {
  //      // get the key
  //      int64_t Key = It.GetKey();
  //      // get the value
  //      TUInt64V Value = It.GetDat();
  //
  //      cout<<"Key: "<<Key<<endl;
  //      print_TUInt64V(Value,"value:");
  // }
  // // get elements by an iterator
  // for (TUInt64UInt64H::TIter It = rank_orig.BegI(); It < rank_orig.EndI();
  // It++) {
  //      // get the key
  //      int64_t Key = It.GetKey();
  //      // get the value
  //      int64_t Value = It.GetDat();
  //
  //      cout<<"Key: "<<Key<<", Value: "<<Value<<endl;
  //      // print_TUInt64V(Value,"value:");
  // }

  eta_seq = eta_seq / no_runs;
  eta_par = eta_par / no_runs;
  eta_par_presn = eta_par_presn / no_runs;
  depth_parallel = depth_parallel / no_runs;

  TStr FuncStr = TStr::Fmt("%s:graph", __func__);
  PrintGStats(FuncStr.CStr(), G);
  if (distance_metric == 1) {
    // std::cout<<" rho_seq: "<<eta_seq<<std::endl;
    std::cout << " rho_par: " << eta_par << ", depth_bins: " << depth_parallel
              << std::endl;
    std::cout << " theta_par: " << eta_par_presn << std::endl;
  } else {
    std::cout << " theta_seq: " << eta_seq << std::endl;
    std::cout << " theta_par: " << eta_par << ", depth_bins: " << depth_parallel
              << std::endl;
  }
}

void calculate_rho_theta_density_peeling_n_peelingplus(const TStr InFNm,
                                                       const char *InFNm_or) {
  PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm, 0, 1);

  TUInt64UInt64VH rank_parallel;
  TUInt64UInt64H rank_node;
  TUInt64UInt64H rank_orig;

  uint64_t G_no_nodes;

  uint64_t rho_unnormalized, nopairsDAG;
  uint64_t rho_extra_unnnormalized, nopairsDAGextra;
  uint64_t no_orig_pairs_same_bin = 0, no_orig_pairs_across_bin = 0,
           no_orig_pairs = 0;

  double theta = 0, rho = 0;
  double rho_extra = 0, theta_extra = 0;
  uint64_t depth_parallel = 0;

  uint64_t kk, vv; // CAREFUL HERE, we use integer, it may be string as well
  std::ifstream infile(InFNm_or);
  while (infile >> kk >> vv) {
    rank_orig.AddDat(kk) = vv;
  }

  G_no_nodes = (uint64_t)G->GetNodes();

  // Structure of rank_parallel = {rank:list of nodes},rank_node = {node:rank}
  printf("Peeling ranking: started \n");
  std::tie(rank_parallel, rank_node) = find_parallel_rank_returnnoderank(G);
  printf("Peeling ranking: finished \n");

  printf("Metric calculation: Peeling - started \n");
  // rho = FindThetaVH_kendaul_tau(rank_parallel,rank_orig,G_no_nodes,p);
  // normln = (G_no_nodes* (G_no_nodes-1)/(double)2);
  // nopairsDAG = no_pairs_DAG(rank_parallel);
  // theta = rho * normln/(double)nopairsDAG;

  // uint64_t no_orig_pairs = 85263786;
  // double no_orig_pairs = 0;
  // for (PNEANet::TObj::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
  //      for (PNEANet::TObj::TNodeI NJ = G->BegNI(); ((NJ < G->EndNI())); NJ++)
  //      {
  //              uint64_t node_u = NI.GetId();
  //              uint64_t node_v = NJ.GetId();
  //              if((node_u!=node_v ) && (rank_orig.GetDat(node_u) !=
  //              rank_orig.GetDat(node_v)))
  //                      no_orig_pairs += 0.5;
  //      }
  // }
  // cout<<"no_orig_pairs: "<<no_orig_pairs<<endl;
  double time_temp = omp_get_wtime();
  std::tie(rho_unnormalized, nopairsDAG) =
      find_prec_recall_real_networks(rank_parallel, rank_orig, G_no_nodes);
  no_orig_pairs_across_bin = nopairsDAG;
  printf("Metric calculation: Peeling - ended \n");
  time_temp = omp_get_wtime() - time_temp;
  printf("Time taken for sequential implementation: %.2fs\n", time_temp);

  time_temp = omp_get_wtime();
  printf("Peeling with same bin pairs ranking and metric started \n");
  std::tie(rho_extra_unnnormalized, nopairsDAGextra, no_orig_pairs_same_bin) =
      Count_extra_pairs_real_networks(G, G_no_nodes, rank_orig, rank_parallel,
                                      rank_node);
  printf("Peeling with same bin pairs ranking and metric ended \n");
  time_temp = omp_get_wtime() - time_temp;
  printf("Time taken for sequential implementation of Peeling+: %.2fs\n",
         time_temp);

  no_orig_pairs = no_orig_pairs_same_bin + no_orig_pairs_across_bin;
  rho = rho_unnormalized / (double)no_orig_pairs;
  rho_extra =
      (rho_extra_unnnormalized + rho_unnormalized) / (double)no_orig_pairs;

  theta = rho_unnormalized / (double)nopairsDAG;
  theta_extra = (rho_extra_unnnormalized + rho_unnormalized) /
                (double)(nopairsDAG + nopairsDAGextra);

  depth_parallel = rank_parallel.Len();

  TStr FuncStr = TStr::Fmt("%s:graph", __func__);
  PrintGStats(FuncStr.CStr(), G);
  std::cout << " rho_peel: " << rho << ", theta_peel: " << theta
            << ", depth_bins: " << depth_parallel
            << ", density_peel: " << rho / theta << std::endl;
  std::cout << " rho_peel+: " << rho_extra << " theta_peel+: " << theta_extra
            << ", density_peel+: " << rho_extra / theta_extra << std::endl;
}

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

  uint64_t kk, vv; // CAREFUL HERE, we use integer, it may be string as well
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
  // p = 0 means nodes in the same bin can not be distinguished
  const float p = Env.GetIfArgPrefixFlt(
      "-p:", 0, "Parameter for partial Kendaul-Tau distance");
  const int no_runs = Env.GetIfArgPrefixInt("-noruns:", 1, "No. of runs");
  const int distance_metric = Env.GetIfArgPrefixInt(
      "-dm:", 1, "Distance metric; 1-Generalized KT, 0-Probabilistic");
  const int choice = Env.GetIfArgPrefixInt(
      "-choice:", 1,
      "0:count_certain_correct_pairs_DAG, 1:calculate_distance_metric");
  const int no_threads =
      Env.GetIfArgPrefixInt("-nothreads:", 4, "No of threads (default = 4)");
  cout << "\n";

  switch (choice) {
  case 0:
    certain_correct_pairs_DAG(InFNm, InFNm_or, no_runs);
    break;

  case 1:
    calculate_distance_metric(InFNm, InFNm_or, no_runs, p, distance_metric);
    break;
  case 2: {
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

    // for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
    //      int aa = NI.GetId();
    //      if (aa == 1597787111) {
    //              printf("Found at node level!!\n");
    //      }
    // }

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
  case 3: {
    // Write rankings multiple related graphs to files
    // Data from the paper "Simple models for human brain fMRI"
    TUInt64UInt64VH rank_parallel;
    TUInt64UInt64H rank_node;
    std::string InFm_1_temp = "../data/brain_network/Simple_Network_Models/"
                              "NV_matrices/brain_nw_edgelist_healthy_";
    std::string InFm_1;
    std::string OutFm_1_temp = "./NV_matrices/node_rank_";
    std::string OutFm_1;
    for (int i = 1; i <= 20; i++) {
      InFm_1 = InFm_1_temp + std::to_string(i) + ".txt";
      PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(TStr(InFm_1.c_str()), 0, 1);
      TStr FuncStr = TStr::Fmt("%s:graph", __func__);
      PrintGStats(FuncStr.CStr(), G);
      cout << "\n";
      std::tie(rank_parallel, rank_node) =
          find_parallel_rank_returnnoderank_PUN(G);
      // const char* outfile = "node_rank.txt";
      OutFm_1 = OutFm_1_temp + std::to_string(i) + ".txt";
      std::ofstream ost{OutFm_1};
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
    }

    std::string InFm_2_temp = "../data/brain_network/Simple_Network_Models/"
                              "COS_matrices/brain_nw_edgelist_COS_";
    std::string InFm_2;
    std::string OutFm_2_temp = "./COS_matrices/node_rank_";
    std::string OutFm_2;
    for (int i = 1; i <= 19; i++) {
      InFm_2 = InFm_2_temp + std::to_string(i) + ".txt";
      PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(InFm_2.c_str(), 0, 1);
      TStr FuncStr = TStr::Fmt("%s:graph", __func__);
      PrintGStats(FuncStr.CStr(), G);
      cout << "\n";
      std::tie(rank_parallel, rank_node) =
          find_parallel_rank_returnnoderank_PUN(G);
      // const char* outfile = "node_rank.txt";
      OutFm_2 = OutFm_2_temp + std::to_string(i) + ".txt";
      std::ofstream ost{OutFm_2};
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
    }
    break;
  }

  case 4: {
    // Human connectome project - Cortex, 100 subjects (4 files or scans per
    // subject) Write rankings multiple related graphs to files
    TUInt64UInt64VH rank_parallel;
    TUInt64UInt64H rank_node;
    std::string InFm_1_temp = "../data/brain_network/Human_connectome_project/"
                              "cortex/graphs_spanningtree_180/"
                              "graph_connectome_";
    std::string InFm_1;
    std::string OutFm_1_temp = "./Human_connectome_result/node_rank_";
    std::string OutFm_1;
    for (int i = 0; i <= 399; i++) {
      InFm_1 = InFm_1_temp + std::to_string(i) + ".txt";
      PUNGraph G = TSnap::LoadEdgeList<PUNGraph>(TStr(InFm_1.c_str()), 0, 1);
      TStr FuncStr = TStr::Fmt("%s:graph", __func__);
      PrintGStats(FuncStr.CStr(), G);
      cout << "\n";
      std::tie(rank_parallel, rank_node) =
          find_parallel_rank_returnnoderank_PUN(G);
      // const char* outfile = "node_rank.txt";
      OutFm_1 = OutFm_1_temp + std::to_string(i) + ".txt";
      std::ofstream ost{OutFm_1};
      for (TUInt64UInt64H::TIter It = rank_node.BegI(); It < rank_node.EndI();
           It++) {
        // get the key
        TUInt64 Key = It.GetKey();
        TUInt64 Value = It.GetDat();
        // write
        ost << Key << "\t" << Value << "\n";
      }
      ost.close();
      cout << "Written results to node_rank.txt" << endl;
      ;
    }
    cout << "file used:" << InFm_1_temp << endl;
    break;
  }

  case 5:
    calculate_rho_theta_density_peeling_n_peelingplus(InFNm, InFNm_or);
    break;
  case 6:
    // STL and parallel implementation
    calculate_rho_theta_density_peeling_n_peelingplus_STL(InFNm, InFNm_or,
                                                          no_threads);
    break;
  default:
    printf("Can not determine what to execute!\n");
  }
  printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());
  return 0;
}
