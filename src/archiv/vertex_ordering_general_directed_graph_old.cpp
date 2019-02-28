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
#include <fstream>

typedef THash<TUInt64, TUInt64>  TUInt64UInt64H;

void print_TIntV(TIntV v, char* S){
    printf("%s: ", S);
    for (int i =0; i < v.Len(); i++){
        printf("%d ", v[i].Val);
    }
    printf("\n");
}


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
    PNEANet G = PNEANet::New();
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

    G->AddNode(0);
    for(int i = 0; i < rv_m_new(gen); i++){
        G->AddEdge(0,0);
        target_list.push_back(0);
        target_list.push_back(0);
    }
    G->AddNode(1);
    for(int i = 0; i < rv_m_new(gen); i++){
        G->AddEdge(1,0);
        target_list.push_back(0);
        target_list.push_back(1);
    }
    int source = 2;
    int source_exist;

    TNEANet::TNodeI NI_temp;


    for (int time_i = 2; time_i < time_n; time_i++){
        if(rv_alpha(gen)){
            target_list_temp.clear();
            G->AddNode(source);
            for(int i =0; i<rv_m_new(gen); i++){
                if(rv_beta(gen)){
                    target = *select_randomly(target_list.begin(), target_list.end());
                    G->AddEdge(source,target);
                    target_list_temp.push_back(source);
                    target_list_temp.push_back(target);
                }
                else{
                    //                    target = rangeRandom(0,(source-1));
                    std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                    target = random_node(gen);
                    G->AddEdge(source,target);
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
                    std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                    source_exist = random_node(gen);
                }
                if (rv_gamma(gen)){
                    target = *select_randomly(target_list.begin(), target_list.end());
                }
                else{
                    std::uniform_int_distribution<int> random_node(0,(source-1)); // guaranteed unbiased
                    target = random_node(gen);
                }
                G->AddEdge(source_exist,target);
                target_list_temp.push_back(source_exist);
                target_list_temp.push_back(target);
            }
            target_list.insert(target_list.end(), target_list_temp.begin(), target_list_temp.end());
        }
    }
    return G;
}

// std::pair<float,float> FindThetaVH(TIntIntVH rank_new, TIntIntH rank_orig, uint64_t G_no_nodes){
float FindThetaVH(TIntIntVH rank_new, TIntIntH rank_orig, uint64_t G_no_nodes){

    std::pair<float,float> eta;

    uint64_t rank_new_len = rank_new.Len();

    uint64_t len_rank_bin_i;
    uint64_t len_rank_bin_j;


    float beta_tilde = 0;
    float beta_ij = 0;
    uint64_t node_v, node_u, node_u_rank;
    TIntV bin_i, bin_j;
    for(uint64_t i=0; i< (rank_new_len-1); i++){
        bin_i = TIntV(rank_new.GetDat(i));
        len_rank_bin_i = bin_i.Len();
        for(uint64_t j = (i+1); j< rank_new_len; j++){
            bin_j = TIntV(rank_new.GetDat(j));
            len_rank_bin_j = bin_j.Len();
            beta_ij = 0;
            for(uint64_t u_i = 0; u_i < len_rank_bin_i; u_i++){
                node_u = bin_i[u_i].Val;
				node_u_rank = rank_orig.GetDat(node_u);
                for(uint64_t v_i = 0; v_i < len_rank_bin_j; v_i++){
                    node_v = bin_j[v_i].Val;
                    if (node_u_rank < rank_orig.GetDat(node_v)){
                        beta_ij++;}
                }
            }
            beta_ij = beta_ij/(len_rank_bin_i*len_rank_bin_j);
            beta_tilde = beta_tilde + beta_ij;
        }
    }
    eta.first = beta_tilde * 2/((rank_new_len-1)*G_no_nodes);
    eta.second = beta_tilde *2/((rank_new_len-1)*rank_new_len);
    return eta.second;
}

// std::pair<float,float> FindThetaH(TIntIntH rank_new, TIntIntH rank_orig, uint64_t G_no_nodes){
float FindThetaH(TIntIntH rank_new, TIntIntH rank_orig, uint64_t G_no_nodes){

    float beta_tilde = 0;
    std::pair<float,float> eta;

    uint64_t rank_new_len = rank_new.Len();
    uint64_t node_u,node_v,node_u_rank;
    for(uint64_t i=0; i<(rank_new_len-1); i++){
        node_u = rank_new.GetDat(i);
		node_u_rank = rank_orig.GetDat(node_u);
        for(uint64_t j = (i+1); j< rank_new_len; j++){
            node_v = rank_new.GetDat(j);
            if (node_u_rank < rank_orig.GetDat(node_v)){beta_tilde = beta_tilde + 1;}
        }
    }
	eta.first = beta_tilde * 2/((rank_new_len-1)*G_no_nodes);
    eta.second = beta_tilde *2/((rank_new_len-1)*rank_new_len);
	return eta.second;
}

float FindThetaH_kendaul_tau(TIntIntH rank_new, TIntIntH rank_orig, uint64_t G_no_nodes, const float p){

    float beta_tilde = 0;
    float eta;

    uint64_t rank_new_len = rank_new.Len();
    uint64_t node_u,node_v,node_u_rank, node_v_rank;
    for(uint64_t i=0; i<(rank_new_len-1); i++){
        node_u = rank_new.GetDat(i);
        // printf("%d\n", node_u);
	    node_u_rank = rank_orig.GetDat(node_u);
        // printf("%d\n", node_u);
        for(uint64_t j = (i+1); j< rank_new_len; j++){
            node_v = rank_new.GetDat(j);
            // printf("%d\n", node_u);
            node_v_rank = rank_orig.GetDat(node_v);
            if (node_u_rank > node_v_rank){
                beta_tilde += 1;}
            else if(node_u_rank == node_v_rank){
                beta_tilde += p;
            }
        }
    }
	eta = beta_tilde /(G_no_nodes* (G_no_nodes-1)/(float)2);
	return eta;
}


float FindThetaVH_kendaul_tau(TIntIntVH rank_new, TIntIntH rank_orig, uint64_t G_no_nodes, const float p){

    float eta;
    uint64_t rank_new_len = rank_new.Len();
    uint64_t len_rank_bin_i;
    uint64_t len_rank_bin_j;
    float beta_tilde = 0;
    float beta_ij = 0;
    uint64_t node_v, node_u, node_u_rank, node_v_rank, node_uu, node_uu_rank;
    TIntV bin_i, bin_j;
    for(uint64_t i=0; i< (rank_new_len-1); i++){
        bin_i = TIntV(rank_new.GetDat(i));
        len_rank_bin_i = bin_i.Len();
        for(uint64_t u_i = 0; u_i < (len_rank_bin_i-1); u_i++){
            node_u = bin_i[u_i].Val;
            node_u_rank = rank_orig.GetDat(node_u);
            for(uint64_t u_ii = u_i+1; u_ii < len_rank_bin_i; u_ii++){
                node_uu = bin_i[u_ii].Val;
                node_uu_rank = rank_orig.GetDat(node_uu);
                if (node_u_rank != node_uu_rank){
                    beta_tilde += p;
                }
            }
        }
        for(uint64_t j = (i+1); j< rank_new_len; j++){
            bin_j = TIntV(rank_new.GetDat(j));
            len_rank_bin_j = bin_j.Len();
            beta_ij = 0;
            for(uint64_t u_i = 0; u_i < len_rank_bin_i; u_i++){
                node_u = bin_i[u_i].Val;
				node_u_rank = rank_orig.GetDat(node_u);
                for(uint64_t v_i = 0; v_i < len_rank_bin_j; v_i++){
                    node_v = bin_j[v_i].Val;
                    node_v_rank = rank_orig.GetDat(node_v);
                    if (node_u_rank > node_v_rank){
                        beta_ij++;
                    }
                    else if(node_u_rank == node_v_rank){
                        beta_ij += p;
                    }
                }
            }
            beta_tilde += beta_ij;
        }
    }

    bin_i = TIntV(rank_new.GetDat((rank_new_len-1)));
    len_rank_bin_i = bin_i.Len();
    for(uint64_t u_i = 0; u_i < (len_rank_bin_i-1); u_i++){
        node_u = bin_i[u_i].Val;
        node_u_rank = rank_orig.GetDat(node_u);
        for(uint64_t u_ii = u_i+1; u_ii < len_rank_bin_i; u_ii++){
            node_uu = bin_i[u_ii].Val;
            node_uu_rank = rank_orig.GetDat(node_uu);
            if (node_u_rank != node_uu_rank){
                beta_tilde += p;
            }
        }
    }
    // printf("%d\n",G_no_nodes);
    eta = beta_tilde / (G_no_nodes* (G_no_nodes-1)/(float)2);
    // eta = beta_tilde;
    // float k = (G_no_nodes * (G_no_nodes-1)/(float)2);
    // float k1 = 96047*96046*0.5;
    // printf("%f\n",k );
    // printf("%f\n",k1);
    return eta;
}

TIntIntH find_seq_rank_unif_Hout(PNEANet G){
    int no_nodes = G->GetNodes();

    PNEANet G_n = TNEANet::New();;
    *G_n = *G;

    TIntIntH rank_new;
    int ttt;

    int i = no_nodes-1;

    TIntV MnDegV;
    int MnDeg;
    int node_sel;

    while(i >0){
        MnDeg = 100*no_nodes;

        for(PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++){
            if (MnDeg > NI.GetInDeg()){MnDegV.Clr(); MnDeg = NI.GetInDeg();}
            if (MnDeg == NI.GetInDeg()){MnDegV.Add(NI.GetId());}
        }

        node_sel = MnDegV[TInt::Rnd.GetUniDevInt(MnDegV.Len())].Val;
        rank_new.AddDat(i,node_sel);

        G_n->DelNode(node_sel);
        i--;
    }
    for (TNEANet::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
        ttt = NI.GetId();}
    rank_new.AddDat(i,ttt);
    return rank_new;
}

TIntIntVH find_seq_rank_unif_VHout(PNEANet G){
    //Returns TIntIntVH, bins
    int no_nodes = G->GetNodes();

    PNEANet G_n = TNEANet::New();;
    *G_n = *G;

    TIntIntVH rank_new;
    int i = no_nodes-1;

    TIntV MnDegV;
    int MnDeg;
    int node_sel;
    TIntV node_selV;

    while(i >0){

        MnDeg = 100*no_nodes;

        for(PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++){
            if (MnDeg > NI.GetInDeg()){MnDegV.Clr(); MnDeg = NI.GetInDeg();}
            if (MnDeg == NI.GetInDeg()){MnDegV.Add(NI.GetId());}
        }

        node_sel = MnDegV[TInt::Rnd.GetUniDevInt(MnDegV.Len())].Val;
        node_selV.Clr();
        node_selV.Add(node_sel);
        rank_new.AddDat(i,node_selV);
        G_n->DelNode(node_sel);
        i--;
    }
    for (TNEANet::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++) {
        node_sel = NI.GetId();
        node_selV.Clr();
        node_selV.Add(node_sel);
        }

    rank_new.AddDat(i,node_selV);
    return rank_new;
}


TIntIntVH find_parallel_rank(PNEANet G){
    int no_nodes = G->GetNodes();

    PNEANet G_n = TNEANet::New();;
    *G_n = *G;

    TIntIntVH rank_new;
    TIntIntVH rank_new_1;

    TIntV MnDegV;
    int MnDeg;
    int depth_dag;
    int i = 0;

    while((G_n->GetNodes()) >0){
        MnDeg = 100*no_nodes;

        MnDegV.Clr();
        for(PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++){

            if (MnDeg > NI.GetInDeg()){MnDegV.Clr(); MnDeg = NI.GetInDeg();}
            if (MnDeg == NI.GetInDeg()){
				MnDegV.Add(NI.GetId());
            }
        }
        rank_new.AddDat(i,MnDegV);

        for (int k = 0; k < MnDegV.Len(); k++) {
                G_n->DelNode(MnDegV[k].Val);
	}
        i++;
    }
    depth_dag = i;
    for(int k =0; k<depth_dag; k++){
        rank_new_1.AddDat((depth_dag-1-k),rank_new.GetDat(k));
    }
    return rank_new_1;
}

//MODIFICATION OF find_parallel_rank WHEN MIN AND MAX OF DEGRESS TO BE PEELED GIVEN
TIntIntVH find_parallel_rank_minmax_deg(PNEANet G, const int min_indeg, const int max_indeg){
    int no_nodes = G->GetNodes();

    PNEANet G_n = TNEANet::New();;
    *G_n = *G;

    TIntIntVH rank_new;
    TIntIntVH rank_new_1;

    TIntV MnDegV;
    int MnDeg;
    int depth_dag;
    int i = 0;

    while((G_n->GetNodes()) >0){
        MnDeg = 100*no_nodes;

        MnDegV.Clr();
        for(PNEANet::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++){
            if ((NI.GetInDeg()>= min_indeg) && (NI.GetInDeg()<= max_indeg)){
				MnDegV.Add(NI.GetId());
                // printf("check\n");
            }
        }
        rank_new.AddDat(i,MnDegV);

        for (int k = 0; k < MnDegV.Len(); k++) {
                G_n->DelNode(MnDegV[k].Val);
	}
        i++;
        // printf("%d\n", i);
    }
    printf("exited\n" );
    depth_dag = i;
    for(int k =0; k<depth_dag; k++){
        rank_new_1.AddDat((depth_dag-1-k),rank_new.GetDat(k));
    }
    return rank_new_1;
}


int main(int argc, char* argv[]){
    Env = TEnv(argc, argv, TNotify::StdNotify);
    Env.PrepArgs(TStr::Fmt("cpm. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
    TExeTm ExeTm;

    const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "./data/cit-HepPh_connected.txt", "Input edgelist file name");

    const TStr InFNm_or_t = Env.GetIfArgPrefixStr("-ior:", "./data/cit-HepPh_connected.txt", "Input original rank file name");

    const char* InFNm_or = InFNm_or_t.CStr();
    const float p=  Env.GetIfArgPrefixFlt("-p:", 0.5, "Parameter for partial Kendaul-Tau distance");
    const int no_runs = Env.GetIfArgPrefixInt("-noruns:", 1, "No. of runs");
    const int distance_metric = Env.GetIfArgPrefixInt("-dm:", 1, "Distance metric; 1-Generalized KT, 0-Probabilistic");


    PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm,0,1);
    TIntIntH rank_seq_uni;
    TIntIntVH rank_parallel;
    TIntIntVH rank_parallel_new;
    TIntV rank_parallel_temp;

    TIntIntH rank_orig;
    float eta_seq = 0;
    // float eta_1_seq = 0;
    float eta_par = 0;
    // float eta_1_par = 0;

    uint64_t G_no_nodes = G->GetNodes();
    // std::pair<float,float> eta_t_seq;
    // std::pair<float,float> eta_t_par;
    float eta_t_seq = 0;//redeclared for Kendaul-Tau distance
    float eta_t_par = 0;//redeclared for Kendaul-Tau distance

    uint64_t kk, vv; //CAREFUL HERE, we use integer, it may be string as well
    // std::ifstream infile("./data/BioGrid_full_human/BioGRID_full_human_net_orig_rank.csv");
    std::ifstream infile(InFNm_or);
    while (infile >> kk >> vv){
        rank_orig.AddDat(kk) = vv;
    }

    float depth_parallel = 0;
    uint64_t i_max = 0;
    for(int ii = 0; ii<no_runs; ii++){
        std::cout<<" run index: "<<ii<<std::endl;

        printf("Sequential ranking: started \n");
        rank_seq_uni = find_seq_rank_unif_Hout(G);
        printf("Sequential ranking: finished \n");
        printf("Peeling ranking: started \n");
        rank_parallel = find_parallel_rank(G);
        printf("Peeling ranking: finished \n");

        const bool OR_BINS_CONSIDER = 0;

        // PEELING TECHNIQUE: IF ORIGINAL RANK HAS BINS
        //======================================================================
        if (OR_BINS_CONSIDER == 1){
            // rank_parallel_new = find_parallel_rank_minmax_deg(G,(int)1, (int)133);
            //COMBINE BINS IN RANK_NEW TO MAKE THE NUMBER OF BINS SMALL
            for (int j = 0; j < 16; j++){
                rank_parallel_temp.Clr();
                if (j==15){
                    i_max = 247;}
                else{
                    i_max = (j*15)+15;}
                for (int i = (j*15); i < i_max; i++){
                    rank_parallel_temp.AddV(rank_parallel.GetDat(i));
                }
                rank_parallel_new.AddDat(j,rank_parallel_temp);
            }
            std::cout<<"Length of new rank_parallel: "<<rank_parallel_new.Len()<<std::endl;
        }
        //======================================================================

        if(distance_metric == 1){
            printf("Eta: Sequential - started \n");
            eta_t_seq = FindThetaH_kendaul_tau(rank_seq_uni, rank_orig,G_no_nodes,p);
            printf("Eta: Sequential - finished \n");
            printf("Eta: Peeling - started \n");
            eta_t_par = FindThetaVH_kendaul_tau(rank_parallel, rank_orig,G_no_nodes,p);
            if(OR_BINS_CONSIDER == 1){
                eta_t_par = FindThetaVH_kendaul_tau(rank_parallel_new, rank_orig,G_no_nodes,p);}
            printf("Eta: Peeling - finished \n");
        }
        else{
            eta_t_seq = FindThetaH(rank_seq_uni, rank_orig,G_no_nodes);
            eta_t_par = FindThetaVH(rank_parallel, rank_orig,G_no_nodes);
            if(OR_BINS_CONSIDER == 1){
                eta_t_par = FindThetaVH(rank_parallel_new, rank_orig,G_no_nodes);}
        }
        eta_seq = eta_seq + eta_t_seq;
        eta_par = eta_par + eta_t_par;
        depth_parallel = depth_parallel + rank_parallel.Len();


        // rank_seq_uni = find_seq_rank_unif_VHout(G);
        // eta_t_seq = FindThetaVH(rank_seq_uni, rank_orig,G_no_nodes);
        // eta_seq = eta_seq + eta_t_seq.first;
        // eta_1_seq = eta_1_seq + eta_t_seq.second;

        // eta_t_par = FindThetaVH(rank_parallel, rank_orig,G_no_nodes);
        // eta_par = eta_par + eta_t_par.first;
        // eta_1_par = eta_1_par + eta_t_par.second;
    }

    eta_seq = eta_seq/no_runs;
    eta_par = eta_par/no_runs;
    depth_parallel = depth_parallel/no_runs;
    // eta_1_seq = eta_1_seq/no_runs;
    // eta_1_par = eta_1_par/no_runs;

    TStr FuncStr = TStr::Fmt("%s:graph", __func__);
    PrintGStats(FuncStr.CStr(), G);
    // std::cout<<" eta_seq: "<<eta_seq<< ",   eta'_seq: "<<eta_1_seq<<std::endl;
    // std::cout<<" eta_par: "<<eta_par<< ",   eta'_par: "<<eta_1_par<<", depth_bins: "<< depth_parallel <<std::endl;

    std::cout<<" eta_seq: "<<eta_seq<<std::endl;
    std::cout<<" eta_par: "<<eta_par<< ", depth_bins: "<< depth_parallel <<std::endl;


    printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());

    return 0;
}
