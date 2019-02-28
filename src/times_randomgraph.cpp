#include "times.hpp"

void pr_success_node_0_in_bin_0(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs)
{
	PNEANet G;
	TUInt64UInt64VH rank_parallel;
	TUInt64V bin_0;
	uint64_t theta_par_pref = 0;
	uint64_t bin_0_len = 0;
	const bool self_loops_allowed = 0;
	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;

		G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);

		rank_parallel = find_parallel_rank(G);
		//Oldest node bin
		bin_0 = rank_parallel.GetDat(0);

		if(bin_0.IsIn(0))
		{
			theta_par_pref += 1;
			bin_0_len += bin_0.Len();
			// print_TUInt64V(bin_0,"bin_0");
		}
	}
	double prob_success = (double)theta_par_pref/no_runs;
	double avg_bin_0_len = (double)bin_0_len/theta_par_pref;

	printf("Prob of success: %0.3f; Avg size of bin 0: %0.3f\n", prob_success,avg_bin_0_len);
}

void perfect_correct_pairs_DAG(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs)
{
	PNEANet G;

	// TUInt64UInt64H rank_seq_uni;
	// TUInt64UInt64VH rank_parallel;
	// TUInt64UInt64VH rank_parallel_new;
	// TUInt64V rank_parallel_temp;
	//
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
	// for(int ii = 0; ii<no_runs; ii++)
	// {
	//      std::cout<<" run index: "<<ii<<std::endl;
	//
	//      G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
	//
	//      G_no_nodes = (uint64_t)G->GetNodes();
	//      for(uint64_t i =0; i < G_no_nodes; i++) {
	//              rank_orig.AddDat(i) = i;
	//      }
	//
	//      std::tie(G_DAG,rank_parallel) = find_parallel_rank_DAG(G);
	//      // write_dot(G,"original_G.dot.dat","G");
	//      // write_dot(G_DAG,"DAG_G.dot.dat","DAG(G)");
	//
	//      eta = count_perfect_correct_pairs_DAG(G_DAG, rank_parallel, rank_orig, G_no_nodes);
	//      eta_cert = eta_cert + eta.first;
	//      eta_ded = eta_ded + eta.second;
	//
	//      rank_parallel_len = rank_parallel.Len();
	//      depth_parallel = depth_parallel + rank_parallel_len;
	// }
	// eta_cert = eta_cert/no_runs;
	// eta_ded = eta_ded/no_runs;
	// depth_parallel = depth_parallel/no_runs;
	//
	// TStr FuncStr = TStr::Fmt("%s:graph", __func__);
	// PrintGStats(FuncStr.CStr(), G);
	// std::cout<<" eta_cert: "<<eta_cert<<", eta_ded: "<<eta_ded<<std::endl;
	// std::cout<<" depth_bins: "<< depth_parallel <<std::endl;
}

void calculate_rho_theta(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs, const float& p, const int& distance_metric)
{
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
	float eta_t_seq = 0;                                                                                                                                                                                                                                                                                            //redeclared for Kendaul-Tau distance
	float eta_t_par = 0;                                                                                                                                                                                                                                                                                            //redeclared for Kendaul-Tau distance

	float depth_parallel = 0;


	TUInt64V bb;
	uint64_t rank_parallel_len;
	const bool self_loops_allowed = 1;
	float normln_temp;
	uint64_t no_pairs;
	float eta_t_par_presn = 0;
	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;

		G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);

		G_no_nodes = (uint64_t)G->GetNodes();
		for(uint64_t i =0; i < G_no_nodes; i++) {
			rank_orig.AddDat(i) = i;
		}

		rank_seq_uni = find_seq_rank_unif_Hout(G);
		rank_parallel = find_parallel_rank(G);

		if(distance_metric == 1) {
			// eta_t_seq = FindThetaH_kendaul_tau(rank_seq_uni, rank_orig,G_no_nodes,p);
			eta_t_par = FindThetaVH_kendaul_tau(rank_parallel, rank_orig,G_no_nodes,p);
			normln_temp = (G_no_nodes* (G_no_nodes-1)/(float)2);
			no_pairs = no_pairs_DAG(rank_parallel);
			eta_t_par_presn = eta_t_par * normln_temp/(float)no_pairs;
		}
		else if(distance_metric == 0) {
			eta_t_seq = FindThetaH(rank_seq_uni, rank_orig,G_no_nodes);
			eta_t_par = FindThetaVH(rank_parallel, rank_orig,G_no_nodes);
		}
		else if(distance_metric == 2) {
		    				//TEST MODE
			float ttt_1,ttt_2;
			float tt = (G_no_nodes* (G_no_nodes-1)/(float)2);
			ttt_1 = FindThetaVH_kendaul_tau(rank_parallel, rank_orig,G_no_nodes,p);
			ttt_2 = FindThetaVH(rank_parallel, rank_orig,G_no_nodes);
			uint64_t no_pairs = no_pairs_DAG(rank_parallel);
			printf("ttt_1: %f, ttt_2: %f, no_pairs: %d, tt: %f\n", ttt_1,ttt_2,no_pairs,tt);
			float ttt_3 = ttt_1*tt/(float)no_pairs;
			printf("rho: %f, theta: %f, theta_Abram: %f\n", ttt_1,ttt_2,ttt_3);
		}
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
	if(distance_metric == 1)
	{
		// std::cout<<" rho_seq: "<<eta_seq<<std::endl;
		std::cout<<" rho_par: "<<eta_par<< ", depth_bins: "<< depth_parallel <<std::endl;
		std::cout<<" theta_par: "<<eta_par_presn<< std::endl;
	}
	else
	{
		std::cout<<" theta_seq: "<<eta_seq<<std::endl;
		std::cout<<" theta_par: "<<eta_par<< ", depth_bins: "<< depth_parallel <<std::endl;
	}
	// printf("%d\n",G->GetEdges());
}

void calculate_rho_theta_with_extra_pairs(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs, const float& p, const int& distance_metric)
{
	PNEANet G;

	TUInt64UInt64VH rank_parallel;
	TUInt64UInt64H rank_node;
	TUInt64UInt64H rank_orig;
	float rho_temp = 0;
	float theta_temp = 0;
	float rho = 0;
	float theta = 0;
	uint64_t G_no_nodes;
	float depth_parallel = 0;
	float normln_temp;
	uint64_t nopairsDAG;

	float rho_extra_temp = 0;
	float theta_extra_temp = 0;
	float rho_extra = 0;
	float theta_extra = 0;
	uint64_t nopairsDAGextra;

	uint64_t rank_parallel_len;
	const bool self_loops_allowed = 1;

	const char* outfile_theta_delta = "theta_delta.txt";
	std::ofstream ost1 {outfile_theta_delta};

	//For Perfect pair calculation
	PNEANet G_DAG;
	TUInt64UInt64VH rank_parallel_1;
	double rho_perf_temp = 0, theta_perf_temp = 0, temp_aaa, temp_bbb, rho_perf = 0, theta_perf = 0;

	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;

		G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);

		G_no_nodes = (uint64_t)G->GetNodes();
		for(uint64_t i = 0; i < G_no_nodes; i++) {
			rank_orig.AddDat(i) = i;
		}
		//Structure of rank_parallel = {rank:list of nodes},rank_node = {node:rank}
		std::tie(rank_parallel,rank_node) = find_parallel_rank_returnnoderank(G);

		rho_temp = FindThetaVH_kendaul_tau(rank_parallel,rank_orig,G_no_nodes,p);
		normln_temp = (G_no_nodes* (G_no_nodes-1)/(float)2);
		nopairsDAG = no_pairs_DAG(rank_parallel);
		rho += rho_temp;
		theta_temp = rho_temp * normln_temp/(float)nopairsDAG;
		theta += theta_temp;

		std::tie(rho_extra_temp,nopairsDAGextra) = Count_extra_pairs(G, G_no_nodes, rank_orig, rank_parallel, rank_node);
		rho_extra += rho_temp + rho_extra_temp;
		theta_extra_temp = (rho_temp + rho_extra_temp) * normln_temp/(float)(nopairsDAG+nopairsDAGextra);
		theta_extra += theta_extra_temp;

		// cout<<"entering perf-precision algorithm"<<endl;
		// std::tie(G_DAG,rank_parallel_1) = find_parallel_rank_DAG(G);
		// std::tie(rho_perf_temp,theta_perf_temp, temp_aaa, temp_bbb) = count_perfect_correct_pairs_DAG(G_DAG, rank_parallel, rank_orig, G_no_nodes);
		// rho_perf += rho_perf_temp;
		// theta_perf += theta_perf_temp;
		// cout<<"leaving perf-precision algorithm"<<endl;

		rank_parallel_len = rank_parallel.Len();
		depth_parallel = depth_parallel + rank_parallel_len;
		// ost1 << theta_extra_temp << "\t" << ((rho_temp + rho_extra_temp)/theta_extra_temp)<< "\n";
		ost1 << theta_perf_temp << "\t" << (rho_perf_temp/theta_perf_temp)<< "\n";
	}
	ost1.close();
	rho = rho/no_runs;
	theta = theta/no_runs;
	rho_extra = rho_extra/no_runs;
	theta_extra = theta_extra/no_runs;
	depth_parallel = depth_parallel/no_runs;
	rho_perf /= no_runs;
	theta_perf /= no_runs;

	TStr FuncStr = TStr::Fmt("%s:graph", __func__);
	PrintGStats(FuncStr.CStr(), G);
	std::cout<<" rho_peel: "<<rho<< ", theta_peel: "<<theta<<", depth_bins: "<< depth_parallel <<", density_peel: "<<(rho/theta)<< std::endl;
	std::cout<<" rho_peel+: "<<rho_extra<< " theta_peel+: "<<theta_extra<<", density_peel+: "<<(rho_extra/theta_extra)<<std::endl;
	std::cout<<" rho_perf "<<rho_perf<< " theta_perf: "<<theta_perf<<", density_perf: "<<(rho_perf/theta_perf)<<std::endl;
}

// Improved Recall and Theta calculations. Contains:
//1. rho and theta of perfect pairs
//2. rho and theta of guessed pairs, which is perfect pairs + guess of order between nodes of different bin+ same bin pairs
//3. rho and theta of Peeling algorithm
void calculate_rho_theta_improvement(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs, const float& p, const int& distance_metric)
{
	PNEANet G, G_DAG;

	TUInt64UInt64VH rank_parallel;
	TUInt64UInt64H rank_node;
	TUInt64UInt64H rank_orig;
	uint64_t G_no_nodes;
	float depth_parallel = 0;
	float rho_perf_tt,theta_perf_tt,rho_guess_tt,theta_guess_tt,rho_peel_tt,theta_peel_tt,density_guess_tt;
	float rho_perf = 0,theta_perf = 0,rho_guess = 0,theta_guess = 0,rho_peel = 0,theta_peel = 0;
	uint64_t rank_parallel_len;
	const bool self_loops_allowed = 1;

	const char *FName = "PA_graph.dat";
	const char *Desc = "PA Graph generated via Jithin's technique";

	const char* outfile_theta_delta = "theta_delta.txt";
	std::ofstream ost1 {outfile_theta_delta};

	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;

		G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
		// TSnap::SaveEdgeList(G, FName, Desc);

		// Load edge list
		// G = TSnap::LoadEdgeList<PNEANet>(FName);


		G_no_nodes = (uint64_t)G->GetNodes();
		for(uint64_t i = 0; i < G_no_nodes; i++) {
			rank_orig.AddDat(i) = i;
		}
		//Structure of rank_parallel = {rank:list of nodes},rank_node = {node:rank}
		std::tie(rank_parallel,rank_node) = find_parallel_rank_returnnoderank(G);
		G_DAG = create_DAG(G,rank_parallel);
		// write_dot(G,"original_G.dot.dat","G");
		// write_dot(G_DAG,"DAG_G.dot.dat","DAG(G)");
		std::tie(rho_perf_tt,theta_perf_tt,rho_guess_tt,theta_guess_tt,rho_peel_tt,theta_peel_tt) = FindRhoThetaPerfectGuess(G, G_DAG, rank_parallel, rank_orig, rank_node, G_no_nodes);


		rho_perf += rho_perf_tt;
		theta_perf += theta_perf_tt;

		rho_guess += rho_guess_tt;
		theta_guess += theta_guess_tt;
		density_guess_tt = rho_guess_tt/theta_guess_tt;

		rho_peel += rho_peel_tt;
		theta_peel += theta_peel_tt;

		ost1 << theta_guess_tt << "\t" << density_guess_tt<< "\n";
	}
	ost1.close();

	rho_perf /= no_runs;
	theta_perf /= no_runs;
	rho_guess /= no_runs;
	theta_guess /= no_runs;
	rho_peel /= no_runs;
	theta_peel /= no_runs;
	// depth_parallel = depth_parallel/no_runs;

	TStr FuncStr = TStr::Fmt("%s:graph", __func__);
	PrintGStats(FuncStr.CStr(), G);
	printf("rho(Rperf): %f, theta(Rperf): %f \n",rho_perf,theta_perf);
	printf("rho(Rguess): %f, theta(Rguess): %f \n",rho_guess,theta_guess);
	printf("rho(Rpeel): %f, theta(Rpeel): %f \n",rho_peel,theta_peel);
}

std::pair<double, double> max_avg_tree_size_DAG(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs)
{
	PNEANet G;

	TUInt64UInt64H rank_seq_uni;
	TUInt64UInt64VH rank_parallel;
	TUInt64UInt64VH rank_parallel_new;
	TUInt64V rank_parallel_temp;


	TUInt64UInt64H rank_orig;
	uint64_t G_no_nodes;
	TUInt64V bb;
	const bool self_loops_allowed = 1;

	double avg_tr = 0;
	double max_tr = 0;
	double avg_tt,max_tt;
	PNEANet G_DAG;
	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;

		G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);

		G_no_nodes = (uint64_t)G->GetNodes();
		for(uint64_t i =0; i < G_no_nodes; i++) {
			rank_orig.AddDat(i) = i;
		}

		std::tie(G_DAG,rank_parallel) = find_parallel_rank_DAG(G);
		write_dot(G,"original_G.dot.dat","G");
		write_dot(G_DAG,"DAG_G.dot.dat","DAG(G)");

		std::tie(avg_tt,max_tt) = count_tree_size_DAG(G_DAG, rank_parallel, G_no_nodes);

		avg_tr = avg_tr + avg_tt;
		max_tr = max_tr + max_tt;
	}
	avg_tr = avg_tr/no_runs;
	max_tr = max_tr/no_runs;

	// TStr FuncStr = TStr::Fmt("%s:graph", __func__);
	// PrintGStats(FuncStr.CStr(), G);
	// std::cout<<" avg_tr: "<<avg_tr<<", max_tr: "<<max_tr<<std::endl;
	return std::make_pair(avg_tr,max_tr);
}

std::pair<double, double> max_avg_tree_height_DAG(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs)
{
	PNEANet G;
	PNEANet G_DAG;
	TUInt64UInt64VH rank_parallel;
	uint64_t G_no_nodes;
	const bool self_loops_allowed = 1;

	double avg_tr = 0;
	double max_tr = 0;
	double avg_tt,max_tt;
	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;
		G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
		G_no_nodes = (uint64_t)G->GetNodes();
		std::tie(G_DAG,rank_parallel) = find_parallel_rank_DAG(G);
		// write_dot(G,"original_G.dot.dat","G");
		// write_dot(G_DAG,"DAG_G.dot.dat","DAG(G)");
		std::tie(avg_tt,max_tt) = count_tree_height_DAG(G_DAG, rank_parallel, G_no_nodes);
		avg_tr = avg_tr + avg_tt;
		max_tr = max_tr + max_tt;
	}
	avg_tr = avg_tr/no_runs;
	max_tr = max_tr/no_runs;
	return std::make_pair(avg_tr,max_tr);
}

double study_levels(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs)
{
	PNEANet G;
	TUInt64UInt64VH rank_parallel;
	double depth_parallel = 0;
	uint64_t rank_parallel_len;
	const bool self_loops_allowed = 1;
	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;

		G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
		rank_parallel = find_parallel_rank(G);
		rank_parallel_len = rank_parallel.Len();
		depth_parallel = depth_parallel + rank_parallel_len;

	}
	depth_parallel = depth_parallel/no_runs;
	// std::cout<<"------------------------------------------"<<std::endl;
	// std::cout<<"depth_bins: "<< depth_parallel <<std::endl;
	// TStr FuncStr = TStr::Fmt("%s:graph", __func__);
	// PrintGStats(FuncStr.CStr(), G);
	return depth_parallel;
}

void find_mean_variance_indegree_slbin(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs, const float& p, const int& distance_metric)
{
	PNEANet G, G_DAG;

	TUInt64UInt64VH rank_parallel;
	TUInt64UInt64H rank_node;
	TUInt64UInt64H rank_orig;
	uint64_t G_no_nodes;
	uint64_t rank_parallel_len;
	const bool self_loops_allowed = 1;
	TUInt64V bin_sl;
	const char *FName = "PA_graph.dat";
	const char *Desc = "PA Graph generated via Jithin's technique";
	PNEANet::TObj::TNodeI NodeI;
	int indeg;
	std::vector<double> indeg_sl_bin;
	double sum,mean,sq_sum,stdev;
	double mean_avg = 0, stdev_avg = 0;
	uint node_u;
	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<" run index: "<<ii<<std::endl;

		G = GenPrefAttachGeneral(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma);
		// Load edge list
		// G = TSnap::LoadEdgeList<PNEANet>(FName);


		G_no_nodes = (uint64_t)G->GetNodes();
		for(uint64_t i = 0; i < G_no_nodes; i++) {
			rank_orig.AddDat(i) = i;
		}
		rank_parallel = find_parallel_rank(G);
		int rank_parallel_len = rank_parallel.Len();

		int upp = (rank_parallel_len-2);
		int low = (rank_parallel_len-2);
		int len_upp_low =(upp-low+1);
		for(int64_t i = upp; i >= low; --i) {
			bin_sl = TUInt64V(rank_parallel.GetDat(i));
			int bin_sl_len = bin_sl.Len();
			for(int u_i = 0; u_i < bin_sl_len; ++u_i) {
				node_u = bin_sl[u_i].Val;
				NodeI = G->GetNI(node_u);
				indeg = NodeI.GetInDeg()/len_upp_low;
				indeg_sl_bin.push_back((float)indeg);
			}
		}

		sum = std::accumulate(indeg_sl_bin.begin(), indeg_sl_bin.end(), 0.0);
		mean = sum / indeg_sl_bin.size();
		sq_sum = std::inner_product(indeg_sl_bin.begin(), indeg_sl_bin.end(), indeg_sl_bin.begin(), 0.0);
		stdev = std::sqrt(sq_sum / indeg_sl_bin.size() - mean * mean);
		indeg_sl_bin.clear();

		mean_avg += mean;
		stdev_avg += stdev;
		cout<<mean_avg<<","<<stdev_avg<<endl;
	}

	mean_avg /= no_runs;
	stdev_avg /= no_runs;

	// depth_parallel = depth_parallel/no_runs;

	TStr FuncStr = TStr::Fmt("%s:graph", __func__);
	PrintGStats(FuncStr.CStr(), G);
	printf("mean_avg: %f, stdev_avg: %f \n",mean_avg,stdev_avg);
}

void estimate_prob_lemma1_sequential(const int& TimeN, const float& pr_alpha, const int& vec_p_1, const int& vec_p_2, const float& pr_beta, const int& vec_q_1, const int& vec_q_2, const float& pr_delta, const float& pr_gamma,const int& no_runs, const int& u, const int& v){
	PNEANet G, G_DAG;
	TUInt64UInt64VH rank_parallel;
	TUInt64UInt64H rank_seq_uni;
	TUInt64V bin_sl;
	const char *FName = "PA_graph.dat";

	// uint64_t u = 8;
	// uint64_t v = 6;

	// const bool self_loops_allowed = 1;
	// G = GenPrefAttachGeneral_undirected(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, self_loops_allowed);
	// const char *Desc = "PA Graph generated via Jithin's technique";
	// TSnap::SaveEdgeList(G, FName, Desc);
	// // //----------------------------------
	// // G_no_nodes = (uint64_t)G->GetNodes();
	// // TUInt64UInt64H rank_orig;
	// // for(uint64_t i = 0; i < G_no_nodes; i++) {
	// //  rank_orig.AddDat(i) = i;
	// // }
	// TUInt64UInt64H rank_node;
	// std::tie(rank_parallel,rank_node) = find_parallel_rank_returnnoderank(G);
	// G_DAG = create_DAG(G,rank_parallel);
	// // float rho_perf_tt,theta_perf_tt,rho_guess_tt,theta_guess_tt,rho_peel_tt,theta_peel_tt;
	// // std::tie(rho_perf_tt,theta_perf_tt,rho_guess_tt,theta_guess_tt,rho_peel_tt,theta_peel_tt) = FindRhoThetaPerfectGuess(G_DAG, rank_parallel, rank_orig, rank_node, G_no_nodes);
	// //----------------------------------
	// write_dot(G,"original_G.dot","G");
	// write_dot(G_DAG,"DAG_G.dot","DAG(G)");

	// Load edge list
	G = TSnap::LoadEdgeList<PNEANet>(FName);

	double emp_average = 0;
	for(int ii = 0; ii<no_runs; ii++)
	{
		std::cout<<"run of Sequential technique: "<<ii<<std::endl;
		rank_seq_uni = find_seq_rank_unif_Hout_node_rank(G);
		if(rank_seq_uni.GetDat(u)<rank_seq_uni.GetDat(v))
			emp_average += 1.0;
	}
	emp_average /= no_runs;

	printf("Pr(node:%d < node:%d) = %f\n",u,v,emp_average);
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
	const int vec_q_1 = Env.GetIfArgPrefixInt("-vecq1:", 1, "edge b/w existing nodes: lb for uniform dist of m");
	const int vec_q_2 = Env.GetIfArgPrefixInt("-vecq2:", 1, "edge b/w existing nodes: ub for uniform dist of m");
	const float pr_delta=  Env.GetIfArgPrefixFlt("-prdelta:", 1, "edge b/w existing nodes: source node, preferentially or uniformly");
	const float pr_gamma=  Env.GetIfArgPrefixFlt("-prgamma:", 1, "edge b/w existing nodes: terminal node, preferentially or uniformly");
	const int no_runs = Env.GetIfArgPrefixInt("-noruns:", 5000, "No. of runs");
	const float p =  Env.GetIfArgPrefixFlt("-p:", 0, "Parameter for partial Kendaul-Tau distance");
	const int distance_metric = Env.GetIfArgPrefixInt("-dm:", 1, "Distance metric; 1:precision and recall, 0:Bin based precision (old), 2:Test mode");
	const int choice = Env.GetIfArgPrefixInt("-choice:", 2, "0:pr_success_node_0_in_bin_0, 1:count_certain_correct_pairs_DAG, 2:calculate_distance_metric, 3:tree size, 4: tree height test, 5: study levels, 6: calculate_distance_metric_with_non-probability-1_pairs");
	const bool rgplot = Env.GetIfArgPrefixInt("-rgplot:", 0, "0:standard mode, 1:write results for various set of parameters to a file");

	const int node_u = Env.GetIfArgPrefixInt("-nodeu:", 0, "nodeu");
	const int node_v = Env.GetIfArgPrefixInt("-nodev:", 1, "nodev");

	printf("\n \n");

	double avg_tr,max_tr;

	switch (choice) {
		case 0:
		pr_success_node_0_in_bin_0(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs);
		break;

	case 1:                                                                                                                                                                                                                                                                                             // Compute Recall and Precision of R_{perf} and R_{peel}
	perfect_correct_pairs_DAG(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs);
	break;

	case 2:                                                                                                                                                                                                                                                                                             // Compute Recall and Precision of R_{peel}
	calculate_rho_theta(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs, p, distance_metric);
	break;

	case 3:
	if(rgplot == 0)
	{
		std::tie(avg_tr, max_tr) = max_avg_tree_size_DAG(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs);
		std::cout<<" avg_tr: "<<avg_tr<<", max_tr: "<<max_tr<<std::endl;
	}
	else {
		std::vector<int> TimeN_v{100, 500, 1000, 2000, 3000, 4000, 5000};
		std::vector<int> vec_p_1_v{5, 10, 15, 20, 25};
		const char* outfile = "results_rg_max_avg_tree.csv";
		std::ofstream ost {outfile};
			// if (!ost) error("can't open output file ",outfile);
		ost << "n" << "," << "m" << "," << "avg_tr" << ","<< "max_tr" << "\n";
		for(auto &i : vec_p_1_v) {
			for(auto &j: TimeN_v) {
				std::tie(avg_tr, max_tr) = max_avg_tree_size_DAG(j, pr_alpha, i, i, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs);
				ost << j << "," << i << "," << avg_tr << ","<< max_tr << "\n";
			}
		}
		ost.close();
	}
	break;
	case 4:
	if(rgplot == 0)
	{
		std::tie(avg_tr, max_tr) = max_avg_tree_height_DAG(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs);
		std::cout<<" avg_tr: "<<avg_tr<<", max_tr: "<<max_tr<<std::endl;
	}
	break;
	case 5:                                                                                                                                                                                                                                                                                             // study number of levels
	double depth_parallel;
	if(rgplot == 0)
	{
		depth_parallel = study_levels( TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma,no_runs);
		std::cout<<"depth_bins: "<< depth_parallel <<std::endl;
	}
	else
	{
		std::vector<int> TimeN_v{1000,2000,5000,10000,25000,50000,75000,100000,500000,1000000};
		std::vector<int> vec_p_1_v{5, 25, 50, 100};
		const char* outfile = "results_level_size.csv";
		std::ofstream ost {outfile};
			// if (!ost) error("can't open output file ",outfile);
		ost << "n" << "," << "m" << "," << "depth_parallel" << "\n";
		for(auto &i : vec_p_1_v)
		{
			for(auto &j: TimeN_v)
			{
				std::cout<<"(n,m): "<<j<<","<<i<<std::endl;
				depth_parallel = study_levels(j, pr_alpha, i, i, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs);
				ost << j << "," << i << "," << depth_parallel <<"\n";
			}
		}
		ost.close();
	}
	break;
	case 6:                                                                                                                                                                                                                                                                                             //Compute Recall and Precision using Peeling method. AND predict order in pairs across bins (apart from Probability 1 pairs). AND predict order in pairs inside bins.
	calculate_rho_theta_with_extra_pairs(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs, p, distance_metric);
	break;
	case 7:                                                                                                                                                                                                                                                                                             //Combines theta and rho calculations of perfect pairs, guessed pairs and peeling algorithm
	calculate_rho_theta_improvement(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs, p, distance_metric);
	break;
	case 8:                                                                                                                                                                                                                                                                                             //Miscellaneous tests
		// find_mean_variance_indegree_slbin(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs, p, distance_metric);
		//Estimate Probability(\pi^{-1} v < \pi^{-1} w) given by Lemma 1 using Sequential technique
	estimate_prob_lemma1_sequential(TimeN, pr_alpha, vec_p_1, vec_p_2, pr_beta, vec_q_1, vec_q_2, pr_delta, pr_gamma, no_runs, node_u, node_v);
	break;
	default:
	printf("Can not determine what to execute!\n");
}

printf("\nrun time: %s (%s)\n", ExeTm.GetTmStr(), TExeTm::GetCurTm());
return 0;
}
