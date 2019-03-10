// (c) 2019 Jithin K. Sreedharan

#ifndef times_hpp
#define times_hpp

#include <inttypes.h>  //To print properly uint64_t
#include <omp.h>
#include <stdio.h>
#include <algorithm>  // for transform
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>  //for bind1st
#include <iomanip>
#include <iostream>
#include <iterator>
#include <map>
#include <numeric>
#include <random>  //C++11 specific
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>
#include "Snap.h"

#undef max
#undef min

typedef THash<TUInt64, TUInt64> TUInt64UInt64H;
typedef THash<TUInt64, TUInt64V> TUInt64UInt64VH;

using std::cout;
using std::endl;
typedef std::unordered_map<uint64_t, uint64_t> Tuint64uint64M;
typedef std::unordered_map<uint64_t, std::vector<uint64_t>> Tuint64uint64VM;

struct HASH {
  size_t operator()(const std::pair<int, int> &x) const {
    return std::hash<long long>()(((long long)x.first) ^
                                  (((long long)x.second) << 32));
  }
};

void print_TIntV(TIntV v, char *S);
void print_TUInt64V(TUInt64V v, char *S);
template <typename Iter, typename RandomGenerator>
Iter select_randomly(Iter start, Iter end, RandomGenerator &g);
template <typename Iter>
Iter select_randomly(Iter start, Iter end);
int rangeRandom(int min, int max);
PNEANet GenPrefAttachGeneral_undirected(
    const int &time_n, const float &pr_alpha, const int &vec_p_1,
    const int &vec_p_2, const float &pr_beta, const int &vec_q_1,
    const int &vec_q_2, const float &pr_delta, const float &pr_gamma,
    const bool &self_loops_allowed);
std::tuple<double, double, double, double> count_perfect_correct_pairs_DAG(
    PNEANet G_DAG, TUInt64UInt64VH rank_new, TUInt64UInt64H rank_orig,
    uint64_t G_no_nodes);
std::pair<double, double> count_tree_size_DAG(PNEANet G_DAG,
                                              TUInt64UInt64VH rank_new,
                                              uint64_t G_no_nodes);
std::pair<double, double> count_tree_height_DAG(PNEANet G_DAG,
                                                TUInt64UInt64VH rank_new,
                                                uint64_t G_no_nodes);
uint64_t no_pairs_DAG(TUInt64UInt64VH rank);
float FindThetaVH(TUInt64UInt64VH rank_new, TUInt64UInt64H rank_orig,
                  uint64_t G_no_nodes);
float FindThetaH_kendaul_tau(TUInt64UInt64H rank_new, TUInt64UInt64H rank_orig,
                             uint64_t G_no_nodes, const float p);
float FindThetaVH_kendaul_tau(TUInt64UInt64VH rank_new,
                              TUInt64UInt64H rank_orig, uint64_t G_no_nodes,
                              const float p);
TUInt64UInt64VH find_parallel_rank(PNEANet G);
std::pair<PNEANet, TUInt64UInt64VH> find_parallel_rank_DAG(PNEANet G);
TIntIntVH find_parallel_rank_minmax_deg(PNEANet G, const int min_indeg,
                                        const int max_indeg);
std::pair<TUInt64UInt64VH, TUInt64UInt64H> find_parallel_rank_returnnoderank(
    PNEANet G);
std::pair<Tuint64uint64VM, Tuint64uint64M>
find_parallel_rank_returnnoderank_STL(PNEANet G);
std::pair<double, uint64_t> Count_extra_pairs(PNEANet G, uint64_t G_no_nodes,
                                              TUInt64UInt64H rank_orig,
                                              TUInt64UInt64VH rank_new,
                                              TUInt64UInt64H rank_node);

PNEANet create_DAG(PNEANet G, TUInt64UInt64VH &rank_parallel);
void write_dot(PNEANet G, const char *FName_t, const char *Desc);

std::pair<uint64_t, uint64_t> find_prec_recall_real_networks(
    TUInt64UInt64VH rank_new, TUInt64UInt64H rank_orig, uint64_t G_no_nodes);

std::tuple<uint64_t, uint64_t, uint64_t> Count_extra_pairs_real_networks(
    PNEANet G, uint64_t G_no_nodes, TUInt64UInt64H rank_orig,
    TUInt64UInt64VH rank_new, TUInt64UInt64H rank_node);

std::pair<TUInt64UInt64VH, TUInt64UInt64H>
find_parallel_rank_returnnoderank_PUN(PUNGraph G);

std::pair<uint64_t, uint64_t> find_prec_recall_real_networks_STL(
    Tuint64uint64VM rank_parallel, Tuint64uint64M rank_orig,
    const uint64_t &G_no_nodes, const int &no_threads);

std::tuple<uint64_t, uint64_t, uint64_t> Count_extra_pairs_real_networks_STL(
    PNEANet G, const uint64_t &G_no_nodes, const Tuint64uint64M &rank_orig,
    const Tuint64uint64VM &rank_new, const Tuint64uint64M &rank_node,
    const int &no_threads);

std::vector<std::vector<double>> find_recall_plot(TUInt64UInt64H rank_node,
                                                  uint64_t width_box,
                                                  uint64_t total_boxes);
std::vector<std::vector<double>> find_recall_plot(Tuint64uint64M rank_node,
                                                  uint64_t width_box,
                                                  uint64_t total_boxes);

double nCk(int n, int k);

// Template definitions
template <class PGraph>
inline void PrintGStats(const char s[], PGraph Graph) {
  printf("graph %s, nodes %d, edges %d, empty %s\n", s, Graph->GetNodes(),
         Graph->GetEdges(), Graph->Empty() ? "yes" : "no");
}

template <typename Iter, typename RandomGenerator>
inline Iter select_randomly(Iter start, Iter end, RandomGenerator &g) {
  std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
  std::advance(start, dis(g));
  return start;
}

template <typename Iter>
inline Iter select_randomly(Iter start, Iter end) {
  static std::random_device rd;
  static std::mt19937 gen(rd());
  return select_randomly(start, end, gen);
}

#endif /* times_hpp */
