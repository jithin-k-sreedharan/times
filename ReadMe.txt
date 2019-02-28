./vertex_ordering_random_graph -timen:5000 -pralpha:1 -vecp1:50 -vecp2:50 -prbeta:1 -vecq1:10 -vecq2:10  -prdelta:1 -prgamma:1 -noruns:1000 -p:0.5 -dm:1 -choice:2

./vertex_ordering_random_graph -timen:5000 -pralpha:0.75 -vecp1:5 -vecp2:50 -prbeta:0.5 -vecq1:5 -vecq2:50  -prdelta:0.5 -prgamma:0.5 -noruns:1000

./vertex_ordering_random_graph -timen:100 -pralpha:1 -vecp1:5 -vecp2:5 -prbeta:1 -vecq1:5 -vecq2:5  -prdelta:0.5 -prgamma:0.5 -noruns:1000

./vertex_ordering_random_graph -timen:5000 -pralpha:1 -vecp1:5 -vecp2:5 -prbeta:1 -vecq1:5 -vecq2:5  -prdelta:1 -prgamma:1 -noruns:1000 -p:0 -dm:1 -choice:2

./vertex_ordering_general_directed_graph  -i:"../data/ArXiv-HEP-TH_from_KDD03/hep-th-citations_weakly_connected.txt" -ior:"../data/ArXiv-HEP-TH_from_KDD03/hep-th-citations_orig_rank.csv" -p:0 -dm:1

./vertex_ordering_general_directed_graph -i:"../data/dynamic-simplewiki/dynamic-simplewiki_cleaned.txt" -ior:"../data/dynamic-simplewiki/dynamic-simplewiki_data_cleaned.csv" -p:0.5 -dm:1

./vertex_ordering_general_directed_graph -i:"../data/dnc_emailleak/dnctemporalGraph_cleaned.txt" -ior:"../data/dnc_emailleak/dnctemporalGraph_data_cleaned.csv" -p:0 -dm:0

./vertex_ordering_general_directed_graph -i:"../data/digg_friends/digg_friends_cleaned.txt" -ior:"../data/digg_friends/digg_friends_data_cleaned.csv" -p:0 -dm:0

./vertex_ordering_general_directed_graph -i:"../data/BioGrid_full_human/BioGRID_full_human_net_cleaned.txt" -ior:"../data/BioGrid_full_human/BioGRID_full_human_net_cleaned_rank.csv" -p:0 -dm:1 -choice:1

./vertex_ordering_random_graph -timen:100 -pralpha:1 -prbeta:1 -vecp1:5 -vecp2:5  -noruns:1000 -p:0 -choice:2

./vertex_ordering_random_graph -timen:10 -pralpha:1 -prbeta:1 -vecp1:2 -vecp2:2  -noruns:1 -choice:4

./vertex_ordering_general_directed_graph -i:"../data/internet_topology/topology_cleaned.txt" -ior:"../data/internet_topology/topology_data_cleaned.csv" -choice:1 -p:0 -dm:1

./vertex_ordering_random_graph -timen:100 -pralpha:1 -prbeta:1 -vecp1:5 -vecp2:5  -noruns:1 -p:0 -choice:1


For WWW 2018:
./vertex_ordering_sequential_predict -timen:50 -pralpha:1 -prbeta:1 -vecp1:3 -vecp2:3 -noruns:1 -norunsMC:400000 -choice:3
./vertex_ordering_sequential_predict -timen:200 -pralpha:1 -prbeta:1 -vecp1:1 -vecp2:20 -noruns:100 -norunsMC:1 -choice:3

 ./vertex_ordering_sequential_predict -timen:50 -pralpha:1 -prbeta:1 -vecp1:3 -vecp2:3 -noruns:100 -norunsMC:1 -choice:3

./vertex_ordering_sequential_predict -noruns:1 -norunsMC:400000 -choice:3
./vertex_ordering_sequential_predict -timen:50 -pralpha:0.5 -vecp1:1 -vecp2:8 -prbeta:0.5 -vecq1:1 -vecq2:8  -prdelta:0.5 -prgamma:0.5 -noruns:100 -norunsMC:100 -choice:3


./vertex_ordering_general_directed_graph -i:"../data/BioGrid_full_human/BioGRID_full_human_net_cleaned.txt" -ior:"../data/BioGrid_full_human/BioGRID_full_human_net_cleaned_rank.csv" -p:0 -choice:3

./vertex_ordering_general_directed_graph -i:"../data/digg_friends/digg_friends_cleaned.txt" -ior:"../data/digg_friends/digg_friends_data_cleaned.csv" -p:0 -choice:3

./vertex_ordering_general_directed_graph -i:"../data/dynamic-simplewiki/dynamic-simplewiki_cleaned.txt" -ior:"../data/dynamic-simplewiki/dynamic-simplewiki_data_cleaned.csv" -p:0 -choice:3

./vertex_ordering_general_directed_graph -i:"../data/dnc_emailleak/dnctemporalGraph_cleaned.txt" -ior:"../data/dnc_emailleak/dnctemporalGraph_data_cleaned.csv" -p:0 -choice:3

./vertex_ordering_general_directed_graph  -i:"../data/ArXiv-HEP-TH_from_KDD03/hep-th-citations_weakly_connected.txt" -ior:"../data/ArXiv-HEP-TH_from_KDD03/hep-th-citations_orig_rank.csv" -p:0 -choice:3

./vertex_ordering_random_graph -timen:5000 -pralpha:1 -prbeta:1 -vecp1:5 -vecp2:5  -noruns:1 -p:0 -choice:6

./vertex_ordering_general_directed_graph -i:"../data/internet_topology/topology_cleaned.txt" -ior:"../data/internet_topology/topology_data_cleaned.csv" -p:0 -choice:3

./vertex_ordering_random_graph -timen:5000 -pralpha:0.75 -vecp1:5 -vecp2:50 -prbeta:0.5 -vecq1:5 -vecq2:50  -prdelta:0.5 -prgamma:0.5 -noruns:1000 -choice:6
./vertex_ordering_random_graph -timen:5000 -pralpha:1 -vecp1:5 -vecp2:50 -prbeta:1 -noruns:1000 -choice:6
./vertex_ordering_random_graph -timen:5000 -pralpha:1 -vecp1:25 -vecp2:25 -prbeta:0 -noruns:1000 -choice:6
----------

./vertex_ordering_general_directed_graph -i:"../data/brain_network/brain_nw_edgelist_temp.txt" -p:0 -choice:3

Bitcoin graph:
./vertex_ordering_general_directed_graph -i:"../data/Bitcoin/graphs_njp/lt_graph_digraph.txt" -ior:"../data/Bitcoin/graphs_njp/lt_graph_node_rank.txt" -p:0 -choice:4
./vertex_ordering_general_directed_graph -i:"../data/Bitcoin/graphs_njp/lt_graph_multidigraph.txt" -ior:"../data/Bitcoin/graphs_njp/lt_graph_node_rank.txt" -p:0 -choice:4
./vertex_ordering_general_directed_graph -i:"../data/Bitcoin/graphs_njp/au_graph_digraph.txt" -ior:"../data/Bitcoin/graphs_njp/au_graph_node_rank.txt" -p:0 -choice:4
./vertex_ordering_general_directed_graph -i:"../data/Bitcoin/graphs_njp/au_graph_multidigraph.txt" -ior:"../data/Bitcoin/graphs_njp/au_graph_node_rank.txt" -p:0 -choice:4


Bitcoin full graph:
./process_temporal_graph -i:"../data/Bitcoin/edge_time.txt" -choice:0
./vertex_ordering_general_directed_graph -i:"WCC.txt" -ior:"predicted_rank.txt" -choice:6 -nothreads:4
./vertex_ordering_general_directed_graph -i:"../data/Bitcoin/edge_time_induced_rank_100000.txt" -ior:"../data/Bitcoin/predicted_rank_eqless_100000.txt" -choice:6 -nothreads:4
./process_temporal_graph -i:"../data/Bitcoin/WCC.txt" -inl:"../data/Bitcoin/predicted_rank_eqless_100000.txt" -choice:2
./process_temporal_graph -i:"../data/Bitcoin/edge_time_induced_rank_100000.txt" -choice:3

Mathoverflow graph:
./vertex_ordering_general_directed_graph -i:"../data/mathoverflow/sx-mathoverflow_digraph.txt" -ior:"../data/mathoverflow/sx-mathoverflow_node_rank.txt" -choice:6

Mathoverflow graph:
./vertex_ordering_general_directed_graph -i:"../data/SMS-A/SD01_digraph.txt" -ior:"../data/SMS-A/SD01_node_rank.txt" -choice:6

FBWall graph:
./vertex_ordering_general_directed_graph -i:"../data/FBWall/facebook-wall_digraph.txt" -ior:"../data/FBWall/facebook-wall_node_rank.txt" -choice:6

CollegeMsg graph:
./vertex_ordering_general_directed_graph -i:"../data/CollegeMsg/CollegeMsg_digraph.txt" -ior:"../data/CollegeMsg/CollegeMsg_node_rank.txt" -choice:6

ArXiv graph:
./vertex_ordering_general_directed_graph -i:"../data/ArXiv-HEP-TH_from_KDD03/hep-th-citations_weakly_connected.txt" -ior:"../data/ArXiv-HEP-TH_from_KDD03/hep-th-citations_orig_rank.csv" -choice:6

Dynamic Simple Wiki:
./vertex_ordering_general_directed_graph -i:"../data/dynamic-simplewiki/dynamic-simplewiki_cleaned.txt" -ior:"../data/dynamic-simplewiki/dynamic-simplewiki_data_cleaned.csv" -choice:6

Protein network:
./vertex_ordering_general_directed_graph -i:"../data/BioGrid_cerevisiae/graph_cerevisiae_cleaned.txt" -ior:"../data/BioGrid_cerevisiae/PH_rank.txt" -choice:6

./vertex_ordering_general_directed_graph -i:"../data/BioGrid_cerevisiae/graph_cerevisiae_cleaned_wapinski07.txt" -ior:"../data/BioGrid_cerevisiae/PH_rank_cerevisiae_wapinski07.txt" -choice:6

./vertex_ordering_general_directed_graph -i:"../data/BioGrid_cerevisiae/graph_cerevisiae_PPODv4_PTHR7-OrthoMCL_wagner1.0.txt" -ior:"../data/BioGrid_cerevisiae/PH_rank_cerevisiae_PPODv4_PTHR7-OrthoMCL_wagner1.0.txt" -choice:6

