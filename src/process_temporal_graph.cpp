// The aim of this program is to process a directed graph
// with entries of the form (u,v,t), where t is the unix timestamp of the edge
// 1. The program finds the largest weak connected component and saves its edgelist
// 2. Saves the rank of each node in the graph of the form 'node rank'
// The code has to be MODIFIED when attributes are not UNSIGNED INTEGER

#include "times.hpp"

int main(int argc, char* argv[])
{
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("cpm. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;

	const TStr InFNm = Env.GetIfArgPrefixStr("-i:", "./data/cit-HepPh_connected.txt", "Input edgelist file name");
	const int choice = Env.GetIfArgPrefixInt("-choice:", 0, "0:Extract WCC for directed graph of the form (u,v,t), 1:Predict the rank of nodes of a directed graph of the form (u,v,t)");
	const TStr InFNm_nl_t = Env.GetIfArgPrefixStr("-inl:", "./data/cit-HepPh_connected.txt", "Input original rank file name");
	const char* InFNm_nl = InFNm_nl_t.CStr();

	switch(choice) {
	case 0: {
		const char* InFNm_cstr = InFNm.CStr();
		std::multimap<std::pair<uint64_t,uint64_t>,uint64_t> Attr;
		std::ifstream edgetime_stream(InFNm_cstr);
		uint64_t node_u, node_v, unixtime;
		while (edgetime_stream >> node_u >> node_v >> unixtime) {
			Attr.insert({std::make_pair(node_u,node_v),unixtime});
		}

		PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm,0,1);
		printf("Nodes = %d, Edges = %d\n", G->GetNodes(), G->GetEdges());
		printf("IsWeaklyConnected(G) = %d\n", TSnap::IsWeaklyConn(G));

		// Get graph with largest WCC
		PNEANet GWMx = TSnap::GetMxWcc(G);
		printf("GWMx->GetNodes() = %d, GWMx->GetEdges() = %d\n", GWMx->GetNodes(), GWMx->GetEdges());

		// PNGraph H = PNGraph::New();
		// H = TSnap::ConvertGraph<PNGraph>(GWMx); //Should be very resource taking
		//** Another idea without converting to H is to create another map wih unique edges as the keys.
		//If a key is already present then skip the rest of the steps. If not, create the key, and write all the edges with this node pair and different attributes.

		std::string OutFNm = "WCC.txt";
		std::ofstream ost {OutFNm};

		// for (PNGraph::TObj::TEdgeI EI = H->BegEI(); EI < H->EndEI(); EI++) {
		//      uint64_t node_v = EI.GetDstNId();
		//      uint64_t node_u = EI.GetSrcNId();
		//      std::pair <std::multimap<std::pair<uint64_t,uint64_t>,float>::iterator, std::multimap<std::pair<uint64_t,uint64_t>,float>::iterator> ret;
		//      ret = Attr.equal_range(std::make_pair(node_u,node_v));
		//      for (auto it=ret.first; it!=ret.second; ++it) {
		//              ost << node_u << "\t" << node_v << "\t" << it->second <<"\n";
		//      }
		// }

		std::unordered_map<std::pair<uint64_t,uint64_t>,bool,HASH> edge_status;
		for (PNEANet::TObj::TEdgeI EI = GWMx->BegEI(); EI < GWMx->EndEI(); EI++) {
			uint64_t node_v = EI.GetDstNId();
			uint64_t node_u = EI.GetSrcNId();
			if (node_u != node_v ) {
				auto search = edge_status.find(std::make_pair(node_u,node_v));
				if (search == edge_status.end()) {
					edge_status[std::make_pair(node_u,node_v)] = 1;
					std::pair <std::multimap<std::pair<uint64_t,uint64_t>,uint64_t>::iterator, std::multimap<std::pair<uint64_t,uint64_t>,uint64_t>::iterator> ret;
					ret = Attr.equal_range(std::make_pair(node_u,node_v));
					for (auto it=ret.first; it!=ret.second; ++it) {
						ost << node_u << "\t" << node_v << "\t" << it->second <<"\n";
					}
				}
			}
		}
		cout<<"\n"<<"Done!"<<endl;
		break;
	}
	case 1: {
		const char* InFNm_cstr = InFNm.CStr();
		std::ifstream edgetime_stream(InFNm_cstr);
		std::unordered_map<uint64_t, uint64_t > node_time;
		std::set<uint64_t> unixtime_unique;
		uint64_t node_u, node_v, unixtime;
		while (edgetime_stream >> node_u >> node_v >> unixtime) {
			// cout<<"node_u: "<<node_u<<", node_v: "<<node_v<<", unixtime: "<<unixtime<<endl;
			auto search = node_time.find(node_u);
			if (search == node_time.end()) {
				node_time.insert({node_u,unixtime});
			}
			else {
				if (unixtime < search->second) {
					search->second = unixtime;
				}
			}
			auto search1 = node_time.find(node_v);
			if (search1 == node_time.end()) {
				node_time.insert({node_v,unixtime});
			}
			else {
				if (unixtime < search1->second) {
					search1->second = unixtime;
				}
			}
			unixtime_unique.insert(unixtime);
		}

		// Automatically sorted
		uint64_t i = 0;
		std::unordered_map<uint64_t,uint64_t> unixtime_rank;
		for(auto const &val: unixtime_unique) {
			unixtime_rank.insert({val,i});
			i++;
		}

		std::string OutFm_1 = "predicted_rank.txt";
		std::ofstream ost {OutFm_1};
		for(auto &val: node_time) {
			uint64_t rank = unixtime_rank.at(val.second);
			ost << val.first << "\t" << rank<< "\n";
		}
		cout<<"\n"<<"Done!"<<endl;
		break;
	}
	case 2: {
		// Extract the induced subgraph with a node list
		PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm,0,1);
		std::ifstream infile(InFNm_nl);
		TIntV induced_subgraph_nodes;
		uint64_t kk, vv;
		while (infile >> kk >> vv) {
			induced_subgraph_nodes.Add(kk);
		}
		PNEANet G_induced_subgraph;
		G_induced_subgraph = TSnap::GetSubGraph(G,induced_subgraph_nodes);
		const char *FName = "edge_time_induced_rank_temp.txt";
		const char *Desc = "edge_time_induced_rank_temp.txt";
		TSnap::SaveEdgeList(G_induced_subgraph, FName, Desc);
		break;
	}
	case 3: {
		// Get the WCC of the directed graph
		PNEANet G = TSnap::LoadEdgeList<PNEANet>(InFNm,0,1);
		PNEANet GWMx = TSnap::GetMxWcc(G);
		printf("G->GetNodes() = %d, G->GetEdges() = %d\n", G->GetNodes(), G->GetEdges());
		printf("GWMx->GetNodes() = %d, GWMx->GetEdges() = %d\n", GWMx->GetNodes(), GWMx->GetEdges());
		const char *FName = "WCC_temp.txt";
		const char *Desc = "WCC_temp.txt";
		TSnap::SaveEdgeList(GWMx, FName, Desc);
		break;
	}
	}
	return 0;
}

// // Save the graph of largest WCC into a file
// const char *FName = "largest_WCC.txt";
// const char *Desc = "largest_WCC graph";
//
// // TSnap::SaveEdgeList(GWMx, FName, Desc);
