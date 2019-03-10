/*(c) 2019 Jithin K. Sreedharan

Useful functions for a given graph in edgelist format with a time stamp
i.e., entries of the form (u,v,t), where t is the unix timestamp of the edge.

The following are the choices available:
  0. The program finds the largest weak connected component of a directed graph
  and saves its edgelist with time stamp
  1. Rank the nodes of the graph based on the time stamp.

Sample usage: ./process_temporal_graph -i:"graph_file.txt" -choice:1 */

#include "times.hpp"
int main(int argc, char *argv[]) {
  Env = TEnv(argc, argv, TNotify::StdNotify);
  Env.PrepArgs(TStr::Fmt("cpm. build: %s, %s. Time: %s", __TIME__, __DATE__,
                         TExeTm::GetCurTm()));
  TExeTm ExeTm;
  // GetIfArgPrefixStr is used to read the commandline arguments
  // GetIfArgPrefixStr format: (prefix to be added, default value, comment)
  const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "./data/temp_file.txt",
                                           "Input edgelist file name");
  const int choice = Env.GetIfArgPrefixInt(
      "-choice:", 0,
      "0:Extract WCC for directed graph of the form (u,v,t), 1:Predict the "
      "rank of nodes of a directed graph of the form (u,v,t)");
  const char *InFNm_cstr = InFNm.CStr();

  switch (choice) {
    case 0: {
      // Extract and write largest weakly connected component of the given graph
      // Considers the graph as DIRECTED.

      std::multimap<std::pair<uint64_t, uint64_t>, uint64_t> Attr;
      std::ifstream edgetime_stream(InFNm_cstr);
      uint64_t node_u, node_v, unixtime;
      while (edgetime_stream >> node_u >> node_v >> unixtime) {
        Attr.insert({std::make_pair(node_u, node_v), unixtime});
      }

      PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm, 0, 1);
      printf("Nodes = %d, Edges = %d\n", G->GetNodes(), G->GetEdges());
      printf("IsWeaklyConnected(G) = %d\n", TSnap::IsWeaklyConn(G));

      // Get graph with largest WCC without time stamps
      PNEANet GWMx = TSnap::GetMxWcc(G);
      printf("GWMx->GetNodes() = %d, GWMx->GetEdges() = %d\n", GWMx->GetNodes(),
             GWMx->GetEdges());

      // Write the largest WCC "with" time stamps
      // The idea is to create another map wih unique edges as the keys.
      // If a key is already present then skip the rest of the steps. If not,
      // create the key, and write all the edges with time stamps of this node
      // pair.

      std::string OutFNm = "WCC_temp.txt";
      std::ofstream ost{OutFNm};
      std::unordered_map<std::pair<uint64_t, uint64_t>, bool, HASH> edge_status;
      for (PNEANet::TObj::TEdgeI EI = GWMx->BegEI(); EI < GWMx->EndEI(); EI++) {
        uint64_t node_v = EI.GetDstNId();
        uint64_t node_u = EI.GetSrcNId();
        if (node_u != node_v) {
          auto search = edge_status.find(std::make_pair(node_u, node_v));
          if (search == edge_status.end()) {
            edge_status[std::make_pair(node_u, node_v)] = 1;
            std::pair<std::multimap<std::pair<uint64_t, uint64_t>,
                                    uint64_t>::iterator,
                      std::multimap<std::pair<uint64_t, uint64_t>,
                                    uint64_t>::iterator>
                ret;
            ret = Attr.equal_range(std::make_pair(
                node_u,
                node_v));  // To find an iterator for all the values of the pair
            for (auto it = ret.first; it != ret.second; ++it) {
              ost << node_u << "\t" << node_v << "\t" << it->second << "\n";
            }
          }
        }
      }
      cout << "\n"
           << "Done!" << endl;
      break;
    }
    case 1: {
      // Find the rank of vertices from the time stamp - if multiple time stamps
      // present for an edge, choose the smallest value. Considers the graph as
      // UNDIRECTED.
      std::ifstream edgetime_stream(InFNm_cstr);
      std::unordered_map<uint64_t, uint64_t> node_time;
      std::set<uint64_t> unixtime_unique;
      uint64_t node_u, node_v, unixtime;
      while (edgetime_stream >> node_u >> node_v >> unixtime) {
        auto search = node_time.find(node_u);
        if (search == node_time.end()) {
          node_time.insert({node_u, unixtime});
        } else {
          if (unixtime < search->second) {
            search->second = unixtime;
          }
        }
        auto search1 = node_time.find(node_v);
        if (search1 == node_time.end()) {
          node_time.insert({node_v, unixtime});
        } else {
          if (unixtime < search1->second) {
            search1->second = unixtime;
          }
        }
        unixtime_unique.insert(unixtime);
      }

      // STL set is automatically sorted
      uint64_t i = 0;
      std::unordered_map<uint64_t, uint64_t> unixtime_rank;
      for (auto const &val : unixtime_unique) {
        unixtime_rank.insert({val, i});
        i++;
      }
      std::string OutFm_1 = "predicted_rank.txt";
      std::ofstream ost{OutFm_1};
      for (auto &val : node_time) {
        uint64_t rank = unixtime_rank.at(val.second);
        ost << val.first << "\t" << rank << "\n";
      }
      cout << "\n"
           << "Done!" << endl;
      break;
    }
    default:
      cout << "Invalid choice" << endl;
  }
  return 0;
}
