#include "vertex_ordering.hpp"

TUInt64UInt64VH find_parallel_rank_PNGraph(PNGraph G)
{
	uint64_t no_nodes = G->GetNodes();

	PNGraph G_n = TNGraph::New();;
	*G_n = *G;

	TUInt64UInt64VH rank_new;
	TUInt64UInt64VH rank_new_1;

	TUInt64V MnDegV;
	uint64_t MnDeg;
	uint64_t depth_dag;
	uint64_t i = 0;

	while(((uint64_t)G_n->GetNodes()) >0)
	{
		MnDeg = 100*no_nodes;

		MnDegV.Clr();
		for(PNGraph::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++)
		{
			if (MnDeg > (uint64_t)NI.GetInDeg())
			{
				MnDegV.Clr(); MnDeg = (uint64_t)NI.GetInDeg();
			}
			if (MnDeg == (uint64_t)NI.GetInDeg())
			{
				MnDegV.Add((uint64_t)NI.GetId());
			}
		}
		rank_new.AddDat(i,MnDegV);

		for (uint64_t k = 0; k < (uint64_t)MnDegV.Len(); k++)
		{
			G_n->DelNode(MnDegV[k].Val);
		}
		i++;
	}
	depth_dag = i;
	for(uint64_t k =0; k<depth_dag; k++)
	{
		rank_new_1.AddDat((depth_dag-1-k),rank_new.GetDat(k));
	}
	return rank_new_1;
}

TUInt64UInt64VH find_parallel_rank_PUNGraph(PUNGraph G)
{
	uint64_t no_nodes = G->GetNodes();

	PUNGraph G_n = TUNGraph::New();;
	*G_n = *G;

	TUInt64UInt64VH rank_new;
	TUInt64UInt64VH rank_new_1;

	TUInt64V MnDegV;
	uint64_t MnDeg;
	uint64_t depth_dag;
	uint64_t i = 0;

	while(((uint64_t)G_n->GetNodes()) >0)
	{
		MnDeg = 100*no_nodes;

		MnDegV.Clr();
		for(PUNGraph::TObj::TNodeI NI = G_n->BegNI(); NI < G_n->EndNI(); NI++)
		{
			if (MnDeg > (uint64_t)NI.GetDeg())
			{
				MnDegV.Clr(); MnDeg = (uint64_t)NI.GetDeg();
			}
			if (MnDeg == (uint64_t)NI.GetDeg())
			{
				MnDegV.Add((uint64_t)NI.GetId());
			}
		}
		rank_new.AddDat(i,MnDegV);

		for (uint64_t k = 0; k < (uint64_t)MnDegV.Len(); k++)
		{
			G_n->DelNode(MnDegV[k].Val);
		}
		i++;
	}
	depth_dag = i;
	for(uint64_t k =0; k<depth_dag; k++)
	{
		rank_new_1.AddDat((depth_dag-1-k),rank_new.GetDat(k));
	}
	return rank_new_1;
}


void calculate_distance_metric_copying(const int& TimeN, const float& pr_beta, const int& m, const int& no_runs)
{
	// PNGraph G;
	PNEANet G;

	TUInt64UInt64H rank_seq_uni;
	TUInt64UInt64VH rank_parallel;
	TUInt64UInt64VH rank_parallel_new;
	TUInt64V rank_parallel_temp;
	TUInt64UInt64H rank_orig;
	float eta_seq = 0;
	float eta_par = 0;
	float eta_par_presn = 0;
	uint64_t G_no_nodes;
	float eta_t_par = 0;//redeclared for Kendaul-Tau distance

	float depth_parallel = 0;

	TUInt64V bb;
	uint64_t rank_parallel_len;
	float normln_temp;
	uint64_t no_pairs;
	float eta_t_par_presn = 0;
	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;

		// G = TSnap::GenCopyModel(TimeN, pr_beta);
		// G = GenCopyModel_jks(TimeN, pr_beta, m);
		G = GenCopyModel_undirected_jks(TimeN, pr_beta, m);

		G_no_nodes = (uint64_t)G->GetNodes();
		for(uint64_t i =0; i < G_no_nodes; i++) {
			rank_orig.AddDat(i) = i;
		}
		rank_parallel = find_parallel_rank(G);
		eta_t_par = FindThetaVH_kendaul_tau(rank_parallel, rank_orig,G_no_nodes,0);
		normln_temp = (G_no_nodes* (G_no_nodes-1)/(float)2);
		no_pairs = no_pairs_DAG(rank_parallel);
		eta_t_par_presn = eta_t_par * normln_temp/(float)no_pairs;

		eta_par = eta_par + eta_t_par;
		eta_par_presn = eta_par_presn + eta_t_par_presn;

		rank_parallel_len = rank_parallel.Len();
		depth_parallel = depth_parallel + rank_parallel_len;
	}

	eta_seq = eta_seq/no_runs;
	eta_par = eta_par/no_runs;
	eta_par_presn = eta_par_presn/no_runs;
	depth_parallel = depth_parallel/no_runs;

	TStr FuncStr = TStr::Fmt("%s:graph", __func__);
	PrintGStats(FuncStr.CStr(), G);
	std::cout<<" rho_par: "<<eta_par<< ", depth_bins: "<< depth_parallel <<std::endl;
	std::cout<<" theta_par: "<<eta_par_presn<< std::endl;
}

void calculate_distance_metric_powerlaw(const int& TimeN, const int& no_runs)
{
	PUNGraph G;

	TUInt64UInt64H rank_seq_uni;
	TUInt64UInt64VH rank_parallel;
	TUInt64UInt64VH rank_parallel_new;
	TUInt64V rank_parallel_temp;
	TUInt64UInt64H rank_orig;
	float eta_seq = 0;
	float eta_par = 0;
	float eta_par_presn = 0;
	uint64_t G_no_nodes;
	float eta_t_seq = 0;//redeclared for Kendaul-Tau distance
	float eta_t_par = 0;//redeclared for Kendaul-Tau distance

	float depth_parallel = 0;

	TUInt64V bb;
	uint64_t rank_parallel_len;
	float normln_temp;
	uint64_t no_pairs;
	float eta_t_par_presn = 0;
	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;

		G = TSnap::GenRndPowerLaw(TimeN, 3);

		G_no_nodes = (uint64_t)G->GetNodes();
		for(uint64_t i =0; i < G_no_nodes; i++) {
			rank_orig.AddDat(i) = i;
		}

		rank_parallel = find_parallel_rank_PUNGraph(G);

		eta_t_par = FindThetaVH_kendaul_tau(rank_parallel, rank_orig,G_no_nodes,0);
		normln_temp = (G_no_nodes* (G_no_nodes-1)/(float)2);
		no_pairs = no_pairs_DAG(rank_parallel);
		eta_t_par_presn = eta_t_par * normln_temp/(float)no_pairs;

		eta_seq = eta_seq + eta_t_seq;
		eta_par = eta_par + eta_t_par;
		eta_par_presn = eta_par_presn + eta_t_par_presn;

		rank_parallel_len = rank_parallel.Len();
		depth_parallel = depth_parallel + rank_parallel_len;

	}

	eta_seq = eta_seq/no_runs;
	eta_par = eta_par/no_runs;
	eta_par_presn = eta_par_presn/no_runs;
	depth_parallel = depth_parallel/no_runs;

	TStr FuncStr = TStr::Fmt("%s:graph", __func__);
	PrintGStats(FuncStr.CStr(), G);
	std::cout<<" rho_par: "<<eta_par<< ", depth_bins: "<< depth_parallel <<std::endl;
	std::cout<<" theta_par: "<<eta_par_presn<< std::endl;
}

int main(int argc, char* argv[]){
	Env = TEnv(argc, argv, TNotify::StdNotify);
	Env.PrepArgs(TStr::Fmt("cpm. build: %s, %s. Time: %s", __TIME__, __DATE__, TExeTm::GetCurTm()));
	TExeTm ExeTm;

	const int TimeN = Env.GetIfArgPrefixInt("-timen:", 100, "time_steps");
	const float pr_beta=  Env.GetIfArgPrefixFlt("-prbeta:", 1, "add node: copy edge or uniformly select end node");
	const int no_runs = Env.GetIfArgPrefixInt("-noruns:", 5000, "No. of runs");
	const int m = Env.GetIfArgPrefixInt("-m:", 5, "m");
	const int choice = Env.GetIfArgPrefixInt("-choice:", 0, "0:copying, 1:configuration_model with power law");

	printf("\n");

	switch (choice) {
	case 0:
		calculate_distance_metric_copying(TimeN, pr_beta, m, no_runs);
		break;
	case 1:
		calculate_distance_metric_powerlaw(TimeN, no_runs);
		break;
	default:
		printf("Can not determine what to execute!\n");
	}

	printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());
	return 0;
}
