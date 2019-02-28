//
//  vertex_ordering.cpp
//  vertex ordering
//
//  Created by Jithin K Sreedharan on 2/20/17.
//  Copyright © 2017 Jithin K Sreedharan. All rights reserved.
//

#include "Snap.h"
#undef max
#undef min

#include "vertex_ordering.hpp"
#include <iostream>
#include <vector>
#include <random>
#include <iterator>
#include <cstdlib>

template<typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator& g) {
    std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
    std::advance(start, dis(g));
    return start;
}

template<typename Iter>
Iter select_randomly(Iter start, Iter end) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    return select_randomly(start, end, gen);
}

template <class PGraph>
void PrintGStats(const char s[], PGraph Graph) {
    printf("graph %s, nodes %d, edges %d, empty %s\n",
           s, Graph->GetNodes(), Graph->GetEdges(),
           Graph->Empty() ? "yes" : "no");
}


int rangeRandom(int min, int max){
    int n = max - min + 1;
    int remainder = RAND_MAX % n;
    int x;
    do{
        x = rand();
    }while (x >= RAND_MAX - remainder);
    return min + x % n;
}


PNEANet GenPrefAttachGeneral(const int& time_n, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma)
{
    PNEANet GraphPt = PNEANet::New();
    TNEANet& G = *GraphPt;
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> rv_m_new(vec_p_1,vec_p_2); //unif distbn is a class
    std::uniform_int_distribution<> rv_m_old(vec_q_1,vec_q_2);
    std::bernoulli_distribution rv_alpha(pr_alpha);
    std::bernoulli_distribution rv_beta(pr_beta);
    std::bernoulli_distribution rv_delta(pr_delta);
    std::bernoulli_distribution rv_gamma(pr_gamma);
    
    
    std::vector<int> target_list;
    std::vector<int> target_list_temp;
    int target;
    
//    TIntV NIdV(time_n, 0);
    // first edge
    G.AddNode(0);
    for(int i = 0; i < rv_m_new(gen); i++){
        G.AddEdge(0,0);
        target_list.push_back(0);
        target_list.push_back(0);
    }
    G.AddNode(1);
    for(int i = 0; i < rv_m_new(gen); i++){
        G.AddEdge(1,0);
        target_list.push_back(0);
        target_list.push_back(1);
    }
    int source = 2;
    int source_exist;
    
    TNEANet::TNodeI NI_temp;

    
    for (int time_i = 2; time_i < time_n; time_i++){
        if(rv_alpha(gen)){
            target_list_temp.clear();
            G.AddNode(source);
            for(int i =0; i<rv_m_new(gen); i++){
                if(rv_beta(gen)){
                    target = *select_randomly(target_list.begin(), target_list.end());
                    G.AddEdge(source,target);
                    target_list_temp.push_back(source);
                    target_list_temp.push_back(target);
                }
                else{
//                    target = random_node(gen);
//                    target = rangeRandom(0,(source-1));
                    std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                    target = random_node(gen);
                    G.AddEdge(source,target);
                    target_list_temp.push_back(source);
                    target_list_temp.push_back(target);
                }
            }
            source ++;
            target_list.insert(target_list.end(), target_list_temp.begin(), target_list_temp.end());
        }
        else{
            target_list_temp.clear();
            for(int i =0 ; i<rv_m_old(gen); i++){
                if (rv_delta(gen)){
                    source_exist = *select_randomly(target_list.begin(), target_list.end());
                }
                else{
//                    source_exist = rangeRandom(0,(source-1));
                    std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                    source_exist = random_node(gen);
                }
                if (rv_gamma(gen)){
                    target = *select_randomly(target_list.begin(), target_list.end());
                }
                else{
//                    target = rangeRandom(0,(source-1));
                    std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                    target = random_node(gen);
                }
                G.AddEdge(source_exist,target);
                target_list_temp.push_back(source_exist);
                target_list_temp.push_back(target);
            }
            target_list.insert(target_list.end(), target_list_temp.begin(), target_list_temp.end());
        }
    }
    return GraphPt;
}
    
std::pair<float,float> FindTheta(TIntIntH rank_new, TIntIntH rank_orig, int G_no_nodes){
    TIntIntVH rank_bin;
    TIntV tempV;
    float eta_t0 = 0;
    float eta_t1;

    for(TIntIntH::TIter It = rank_new.BegI(); It < rank_new.EndI(); It++){
        tempV.Clr();
        if(rank_bin.IsKey(It.GetDat())){
            tempV = rank_bin.GetDat(It.GetDat());
        }
        tempV.Add(It.GetKey());
        rank_bin.AddDat(It.GetDat(),tempV);
//        rank_bin[It.GetDat()].Add(It.GetKey());
    }
 
   
    int rank_bin_len = rank_bin.Len();
    int len_rank_bin_i;
    int len_rank_bin_j;

    TIntV temp_i, temp_j;
    
    std::pair<float,float> eta;
    for(int i=0; i< (rank_bin_len-1); i++){
        temp_i = TIntV(rank_bin.GetDat(i));
        for(TIntV::TIter It = temp_i.BegI(); It < temp_i.EndI(); It++){
            eta_t1 = 0;
            len_rank_bin_i = temp_i.Len();
            for(int j = (i+1); j< rank_bin_len; j++){
                temp_j = TIntV(rank_bin.GetDat(j));
                len_rank_bin_j = temp_j.Len();
                for(int k = 0; k < len_rank_bin_j; k++){
                    if (rank_orig.GetDat(It->Val) < rank_orig.GetDat(temp_j[k])){
                        eta_t1 = eta_t1 + 1/(len_rank_bin_i*len_rank_bin_j);
                    }
                }
            }
            eta_t0 = eta_t0 + eta_t1;
        }
    }
//    printf("bin length: %d \t ",rank_bin_len);
//    printf("G_no_nodes: %d \t ",G_no_nodes);
//    printf("eta_t0: %g \t ",eta_t0);
    eta.second = eta_t0 *2/((rank_bin_len-1)*rank_bin_len);
    eta.first = eta_t0 * 2/((rank_bin_len-1)*G_no_nodes);
    return eta;
}

TIntIntH find_seq_rank_unif(PNEANet G){
    int no_nodes = G->GetNodes();

    PNEANet G_n = PNEANet(G);
//    int min_t;
    TIntIntH rank_new;
    
    TIntV rank_bin_0;
    TNEANet::TNodeI NI_temp;
    int ttt;
    
    TIntV tempV;
    TIntIntVH degree_list_ard;
    
    int i = no_nodes-1;
    
    TIntV MnDegV;
    int MnDeg;
    int node_sel;
    
    while(i >0){
//        degree_list_ard.Clr();
//        for(TNEANet::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++){
//            tempV.Clr();
//            if(degree_list_ard.IsKey(NI.GetDeg())){
//                tempV = degree_list_ard.GetDat(NI.GetDeg());
//            }
//            tempV.Add(NI.GetId());
//            degree_list_ard.AddDat(NI.GetDeg(),tempV);
//        }
//        
//        //Find minimum key in degree_list_ard
//        min_t = 2*no_nodes;
//        for(TIntIntVH::TIter It = degree_list_ard.BegI(); It < degree_list_ard.EndI(); It++){
//            if(It.GetKey() < min_t){
//                min_t = It.GetKey();}
//        }
//        
//
//        rank_bin_0 = TIntV(degree_list_ard.GetDat(min_t));
//        node_sel = *select_randomly(rank_bin_0.BegI(), rank_bin_0.EndI());
////        if(G_n->IsNode(node_sel)){
////            G_n->DelNode(node_sel);
////        }
////        else{
////            printf("success\n");
////            std::cout<<"node_sel: "<<node_sel<<std::endl;
////        }
        
        MnDeg = 3*no_nodes;
        
        for(PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++){
            if (MnDeg > NI.GetDeg()){MnDegV.Clr(); MnDeg = NI.GetDeg();}
            if (MnDeg == NI.GetDeg()){MnDegV.Add(NI.GetId());}
        }
        
        node_sel = MnDegV[TInt::Rnd.GetUniDevInt(MnDegV.Len())];
        rank_new.AddDat(node_sel,i);

        
//        144   TIntV MxDegV;
//        145   int MxDeg=-1;
//        146   for (typename PGraph::TObj::TNodeI NI = Graph->BegNI(); NI < Graph->EndNI(); NI++) {
//            147     if (MxDeg < NI.GetDeg()) { MxDegV.Clr(); MxDeg = NI.GetDeg(); }
//            148     if (MxDeg == NI.GetDeg()) { MxDegV.Add(NI.GetId()); }
//            149   }
//        150   EAssertR(! MxDegV.Empty(), "Input graph is empty!");
//        151   return MxDegV[TInt::Rnd.GetUniDevInt(MxDegV.Len())];
        
        
//        if(G_n->IsNode(node_sel))
//        {G_n->DelNode(node_sel);}
//        else{printf("entha?");}
        G_n->DelNode(node_sel);
        i--;
    }
    for (TNEANet::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
        ttt = NI.GetId();}
    rank_new.AddDat(ttt,i);
    return rank_new;
}



int main(int argc, char* argv[]){
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

//
//    const int TimeN = 100;
//    const float pr_alpha = 1;
//    const int vec_p_1 = 5;
//    const int vec_p_2 = 5;
//    const float pr_beta = 1;
//    const int vec_q_1 = 2;
//    const int vec_q_2 = 5;
//    const float pr_delta = 0.5;
//    const float pr_gamma = 0.5;
//    const int no_runs = 10000;
    
    PNEANet G;
    TIntIntH rank_seq_uni;
    TIntIntH rank_orig;
    float eta = 0;
    float eta_1 = 0;
    
    int G_no_nodes;
    std::pair<float,float> eta_t;
    
    bool flag_no_orig_rank =0;
    
    if((pr_alpha==1) && (pr_beta==1)){
        flag_no_orig_rank = 1;
        G_no_nodes = TimeN;
        for(int i =0; i < G_no_nodes; i++){
            rank_orig.AddDat(i) = i;
        }
    }

    for(int ii = 0; ii<no_runs; ii++){
        std::cout<<" run index: "<<ii<<std::endl;
        // Generate undirected graph
        G = GenPrefAttachGeneral(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma);

        if(~flag_no_orig_rank){
            G_no_nodes = G->GetNodes();
            for(int i =0; i < G_no_nodes; i++){
                rank_orig.AddDat(i) = i;
            }
        }
        
        rank_seq_uni = find_seq_rank_unif(G);
        
//        for(TIntIntH::TIter It = rank_seq_uni.BegI(); It < rank_seq_uni.EndI(); It++){
//            std::cout<<" "<<It.GetKey()<<": "<<It.GetDat()<<")"<<std::endl;
//        }
//
//
//        for(TIntIntH::TIter It = rank_orig.BegI(); It < rank_orig.EndI(); It++){
//            std::cout<<" "<<It.GetKey()<<": "<<It.GetDat()<<")"<<std::endl;
//        }

        eta_t = FindTheta(rank_seq_uni, rank_orig, G_no_nodes);
        eta = eta + eta_t.first;
        eta_1 = eta_1 + eta_t.second;
    }
    
    std::cout<<"G_no_nodes: "<<G_no_nodes<<std::endl;
    eta = eta/no_runs;
    eta_1 = eta_1/no_runs;
    
//    TStr FuncStr = TStr::Fmt("%s:graph", __func__);
//    PrintGStats(FuncStr.CStr(), G);
    std::cout<<" eta: "<<eta<< ",   eta': "<<eta_1<<std::endl;

    printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());

    return 0;
}


// ARCHIVE *************************************
// *********************************************

// // TO PRINT A HASH TABLE

//    for(TIntIntH::TIter It = rank_seq_uni.BegI(); It < rank_seq_uni.EndI(); It++){
//        std::cout<<" "<<It.GetKey()<<": "<<It.GetDat()<<")"<<std::endl;
//    }

// // TRAVERSING THE GRAPH

//    for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
//        printf("node id %d with out-degree %d and in-degree %d\n",
//               NI.GetId(), NI.GetOutDeg(), NI.GetInDeg());
//    }
//

//    // traverse the edges
//    for (TNEANet::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
//        printf("edge (%d, %d)\n", EI.GetSrcNId(), EI.GetDstNId());
//    }


//
//    for (int node = 2; node < Nodes; node++) {
//        NodeSet.Clr(false);
//        while (NodeSet.Len() < NodeOutDeg && NodeSet.Len() < node) {
//            NodeSet.AddKey(NIdV[TInt::Rnd.GetUniDevInt(NIdV.Len())]);
//        }
//        const int N = G.AddNode();
//        for (int i = 0; i < NodeSet.Len(); i++) {
//            Graph.AddEdge(N, NodeSet[i]);
//            NIdV.Add(N);
//            NIdV.Add(NodeSet[i]);
//        }
//    }
//    return GraphPt;
//}


//// Generate PrefAttach graph
//void GeneratePrefAttachGraph() {
//    const int NNodes = 100;
//    const int NodeOutDeg = 15;
//
//    PUNGraph UNGraph;
//
//    // Generate undirected graph
//    UNGraph = TSnap::GenPrefAttach(NNodes, NodeOutDeg);
//
//    TStr FuncStr = TStr::Fmt("%s:UNgraph", __func__);
//    PrintGStats(FuncStr.CStr(), UNGraph);
//}

//PUNGraph GenPrefAttach(const int& Nodes, const int& NodeOutDeg, TRnd& Rnd) {
//    PUNGraph GraphPt = PUNGraph::New();
//    TUNGraph& Graph = *GraphPt;
//    Graph.Reserve(Nodes, NodeOutDeg*Nodes);
//    TIntV NIdV(NodeOutDeg*Nodes, 0);
//    // first edge
//    Graph.AddNode(0);  Graph.AddNode(1);
//    NIdV.Add(0);  NIdV.Add(1);
//    Graph.AddEdge(0, 1);
//    TIntSet NodeSet;
//    for (int node = 2; node < Nodes; node++) {
//        NodeSet.Clr(false);
//        while (NodeSet.Len() < NodeOutDeg && NodeSet.Len() < node) {
//            NodeSet.AddKey(NIdV[TInt::Rnd.GetUniDevInt(NIdV.Len())]);
//        }
//        const int N = Graph.AddNode();
//        for (int i = 0; i < NodeSet.Len(); i++) {
//            Graph.AddEdge(N, NodeSet[i]);
//            NIdV.Add(N);
//            NIdV.Add(NodeSet[i]);
//        }
//    }
//    return GraphPt;
//}

//PNEGraph GenPrefAttachGeneral(const int& time_n, const int& NodeOutDeg, TRnd& Rnd)
