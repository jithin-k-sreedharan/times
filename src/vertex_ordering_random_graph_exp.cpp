// FIDNING EXPECTED BIN SIZES

#include "vertex_ordering.hpp"

void write_dot(PNEANet G,const char *FName_t,const char *Desc)
{
    // Output node IDs as numbers
    TIntStrH NIdLabelH;
    // Generate labels for random graph
    for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
      NIdLabelH.AddDat(NI.GetId(), TStr::Fmt("Node%d", NI.GetId()));
    }
    TSnap::SaveGViz(G, FName_t, Desc,NIdLabelH);
}


int main(int argc, char* argv[])
{
    Env = TEnv(argc, argv, TNotify::StdNotify);
    Env.PrepArgs(TStr::Fmt("cpm. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
    TExeTm ExeTm;

    const int TimeN = Env.GetIfArgPrefixInt("-timen:", 100, "time_steps");
    const float pr_alpha=  Env.GetIfArgPrefixFlt("-pralpha:", 1, "add node or not");
    const int vec_p_1 = Env.GetIfArgPrefixInt("-vecp1:", 5, "add node: lb for uniform dist of m");
    const int vec_p_2 = Env.GetIfArgPrefixInt("-vecp2:", 5, "add node: ub for uniform dist of m");
    const float pr_beta=  Env.GetIfArgPrefixFlt("-prbeta:", 1, "add node: preferentially or uniformly");
    const int vec_q_1 = Env.GetIfArgPrefixInt("-vecq1:", 5, "edge b/w existing nodes: lb for uniform dist of m");
    const int vec_q_2 = Env.GetIfArgPrefixInt("-vecq2:", 5, "edge b/w existing nodes: ub for uniform dist of m");
    const float pr_delta=  Env.GetIfArgPrefixFlt("-prdelta:", 1, "edge b/w existing nodes: source node, preferentially or uniformly");
    const float pr_gamma=  Env.GetIfArgPrefixFlt("-prgamma:", 1, "edge b/w existing nodes: terminal node, preferentially or uniformly");
    const int no_runs = Env.GetIfArgPrefixInt("-noruns:", 1000, "No. of runs");
    const float p =  Env.GetIfArgPrefixFlt("-p:", 0.5, "Parameter for partial Kendaul-Tau distance");
    const int distance_metric = Env.GetIfArgPrefixInt("-dm:", 1, "Distance metric; 1-Generalized KT, 0-Probabilistic");


    PNEANet G;
    PNEANet G_DAG;
    PUNGraph G_1;

    TUInt64UInt64H rank_seq_uni;
    TUInt64UInt64VH rank_parallel;
    TUInt64UInt64VH rank_parallel_new;
    TUInt64V rank_parallel_temp;


    TUInt64UInt64H rank_orig;
    float eta_seq = 0;
    float eta_par = 0;

    uint64_t G_no_nodes;
    float eta_t_seq = 0;//redeclared for Kendaul-Tau distance
    float eta_t_par = 0;//redeclared for Kendaul-Tau distance

    float depth_parallel = 0;

    //---------------------------------------------
    std::vector<float> a;
    // TFltV a;
    uint64_t n;
    double bin_s_t;
    n = TimeN;
    int m = vec_p_1;
    bin_s_t = 0;

    int iii = 0;
    while(n > 0){
        bin_s_t = (n*2.0/(float)(m+2));
        n = n - bin_s_t;
        a.push_back(bin_s_t);
        iii++;
    }
    printf("Total expected bins: %d\n", a.size());
    //---------------------------------------------

    std::vector<float> bb;
    std::vector<int> bb_count;

    uint64_t rank_parallel_len;
    for(int ii = 0; ii<no_runs; ii++)
    {
        std::cout<<" run index: "<<ii<<std::endl;

        // G = GenPrefAttachGeneral(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma);
        G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma);
        std::tie(G_DAG,rank_parallel) = find_parallel_rank_DAG(G);
        // G_1 = TSnap::ConvertGraph<PUNGraph>(G);

        G_no_nodes = (uint64_t)G->GetNodes();
        for(uint64_t i =0; i < G_no_nodes; i++){
            rank_orig.AddDat(i) = i;
        }

        G_no_nodes = (uint64_t)G->GetNodes();
        // rank_parallel = find_parallel_rank(G);

        // rank_parallel_len = rank_parallel.Len();
        // depth_parallel = depth_parallel + (float)rank_parallel_len;
        //
        // for (int i = 0; i < rank_parallel_len; i++)
        // {
        //     TUInt64V bb_temp = rank_parallel.GetDat(i);
        //     if (i >= (bb.size()))
        //     {
        //         bb.push_back(bb_temp.Len());
        //         bb_count.push_back(1);
        //     }
        //     else
        //     {
        //         bb[i] = bb[i] + bb_temp.Len();
        //         bb_count[i] = bb_count[i] + 1;
        //     }
        // }
    }
//-------- Writing Dot file
// TSnap::GetShortPath(G, )

// write_dot(G_1, "undirected_G.dot.dat", "undirected(G)");
write_dot(G,"original_G.dot.dat","G");
write_dot(G_DAG,"DAG_G.dot.dat","DAG(G)");

for(int i = 0; i < rank_parallel.Len(); i++)
{
    TUInt64V bin_ttt = rank_parallel.GetDat(i);
    print_TUInt64V(bin_ttt, "Bin");
}
count_certain_correct_pairs_DAG(G_DAG, rank_parallel, rank_orig, G_no_nodes);

// TUInt64H NIdToDistH;
// TUInt64V KeyV;
// TUInt64UInt64VH Rchble_nodes;
// for(PNEANet::TObj::TNodeI NI = G_DAG->BegNI(); NI < G_DAG->EndNI(); NI++)
// {
//     TSnap::GetShortPath(G_DAG, NI.GetId(), NIdToDistH, IsDir=true);
//     NIdToDistH.GetKeyV(KeyV);
//     Rchble_nodes(NI.GetId(),KeyV);
// }
//
// TUInt64V bin_tt;
// for(int i=0; i <G_no_nodes;i++)
// {
//     printf("%d\n",i);
//     bin_tt = TUInt64V(Rchble_nodes.GetDat(i));
//     print_TUInt64V(Rchble_nodes, "reachable");
// }
//------------------------------


depth_parallel = depth_parallel/no_runs;

uint64_t bb_Len = bb.size();
for (int i = 0; i < bb_Len; i++)
{
    bb[i] = bb[i] / bb_count[i];
}
std::reverse(bb.begin(), bb.end());

FILE *F = fopen("./compare_bins.txt", "w");
fprintf(F, "No. of nodes: %d, m: %d\n", TimeN, vec_p_1);
fprintf(F, "total-empirically-averaged-simulated-bins: %f, total-theoretical-bins: %d\n", depth_parallel, a.size());
// fprintf(F, "simulated-bin-size/n theoretical-bin-size/n\n");
fprintf(F, "simulated-bin-size  theoretical-bin-size\n");
fprintf(F, "========================================\n");

int min_bin_size  = std::min((int)a.size(),(int)bb_Len);
bool flag_a = 0, flag_b = 0;
if(a.size() <= bb_Len){
    flag_a = 1;
}
else if (a.size()>bb_Len){
    flag_b = 1;
}

//BIN SIZE WRITING...
//-------------------
for (int i= 0; i < min_bin_size; i++) {
    // fprintf(F, "%0.4f %0.4f\n", bb[i]/TimeN, a[i]/TimeN);
    fprintf(F, "%0.4f %0.4f\n", bb[i], a[i]);
}
if(flag_a ==1){
    for (int i = (bb_Len - min_bin_size);  i < bb_Len ; i++)
    {
        // fprintf(F, "%0.4f\n", bb[i]/TimeN);
        fprintf(F, "%0.4f\n", bb[i]);
    }
}
if(flag_b ==1){
    for (uint64_t i = (a.size() - min_bin_size);  i < a.size(); i++){
        // fprintf(F, "       %0.4f\n", a[i]/TimeN);
        fprintf(F, "       %0.4f\n", a[i]);
    }
}
fclose(F);

    TStr FuncStr = TStr::Fmt("%s:graph", __func__);
    PrintGStats(FuncStr.CStr(), G);
    std::cout<<" eta_seq: "<<eta_seq<<std::endl;
    std::cout<<" eta_par: "<<eta_par<< ", depth_bins: "<< depth_parallel <<std::endl;


    printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());

    return 0;
}
