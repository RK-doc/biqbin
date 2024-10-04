//============================================================================
// Name        : StableSetBaB.cpp
// Author      : Timotej Hrga
// Version     : 
// Copyright   :
// Description : Parallel version of Branch and Bound algorith for Stable set
//  			  problem using MPI based on serial code of Melanie Siebenhofer.
//============================================================================

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <list>
#include <limits>	//std::numeric::limits
#include <string>
#include <time.h>	//* clock_t, clock, CLOCKS_PER_SEC */
#include <signal.h>	//sigaction
#include <unistd.h> //alarm
#include <cmath>	//std::pow
//#include <iomanip>  //std::fixed, std::setprecision

#include "Subproblem.h"
#include "BoundComputationTrivial.h"
#include "BoundComputationTheta.h"
#include "BoundComputationExactSubgraph.h"
#include "BoundComputationBundle.h"


#include <mpi.h>

#define THETA
//#define DEBUG

#define SECONDS 7200

#define SEND_TO_WORKERS 1
#define GET_CHILDREN 2
#define FATHOM 3
#define DONE 4
#define NEW_VALUE 5


//global variables
int g_lowerBound;
LowerBoundComputation* g_lbndComp;
UpperBoundComputation* g_ubndComp;
Params params;
int counterBoundComputationFailed = 0;
int counterMosekSolutionNearOptimal = 0;

// FUNCTIONS
// read parameters
void read_params(){
	//program
	params.DFS = 1;
	params.minSizeGraphSizeForFullEnum = 3;
	params.subproblemWithVertexInSetFirst = 1;
	params.strategyBranchVar = 2; //1: random, 2: min, 3:max, 4:closest to mean, 5: closest to 0.5, 6: median
	params.lowerBoundForEveryNth = 1;
	params.lowerBoundInRoot = 3;
	params.maxExecutionTime = 60*60*1; //1 hour
	params.numberExactSubgraphs[3] = 500;
	params.numberExactSubgraphs[4] = 500;
	params.numberExactSubgraphs[5] = 500;
	params.numberExactSubgraphs[6] = 0;
	params.numberExactSubgraphs[7] = 0;
	
	//bundle
	params.mL = 0.075;
	params.epsStop = 0.0001;
	params.Lmax = 30;
	params.muSSFactor = 0.99;
	params.muNSFactor = 1.02;
	params.muMin =1e-10;
	params.mr = 0.5;

	//output
	params.printlevel = 1;
	params.outputEdgesSubproblem = 0;
	params.outputVerticesSubproblem = 0;
	params.outputVertexLocal = 0;
	params.outputVertexGlobal = 0;
	params.logfile = 0;
}

// read input
void read_input(char** argv, int& n, int& m, std::vector<Vertex>& vertices, std::vector<Edge>& edges){

	std::ifstream input(argv[1]);

	// read number of vertices
	input >> n;

	// read number of edges
	input >> m;

	//create list of vertices
	for(int i = 1; i <= n; i++){
		vertices.emplace_back(i,i);
	}

	//read edges and store in list, ignore weights
	int v1, v2, weight;
	for(int i = 0; i<m; i++)
	{
		input >> v1;
		input >> v2;
		input >> weight;
	
		//add edge to list
		if(v1 < v2)
		{
			edges.emplace_back(v1,v2);
		}
		else
		{
			edges.emplace_back(v2,v1);
		}
	}
	input.close();
}

// print info solution
void print_solution(int n, int sol[], int val, int counterBaBNodes, int counterPrune, double timeTotal, int maxWorkers)
{
	std::cout << "====================================================================" << std::endl;
	std::cout << "SOLUTION: "  << std::endl;

	std::cout << "optimal solution: (";
	for(int i = 0; i < n-1; i++){
		std::cout << sol[i] << ", ";
	}
	std::cout << sol[n-1] << ")" << std::endl;

	std::cout << "value: "<< val << std::endl;

	std::cout << "b&b nodes: " << counterBaBNodes << std::endl;
	std::cout << "pruned by bound: " << counterPrune << std::endl;
	std::cout << "total time: " << timeTotal << std::endl;
	std::cout << "max workers: " << maxWorkers << std::endl;
	std::cout << "====================================================================" << std::endl;
}

// print info solution when timer expired
void print_solution_expired(int n, int sol[], int val,int counterBaBNodes, double best_upper_bound, int left_subproblems)
{

	std::cout << "====================================================================" << std::endl;
	std::cout << "Time limit " << SECONDS <<  " exceeded.\n";
	std::cout << "APPROXIMATE SOLUTION"  << std::endl;

	std::cout << "(";
	for(int i = 0; i < n-1; i++){
		std::cout << sol[i] << ", ";
	}
	std::cout << sol[n-1] << ")" << std::endl;

	std::cout << "value: "<< val << std::endl;
	std::cout << "highest upper bound: " << best_upper_bound << std::endl;

	std::cout << "b&b nodes: " << counterBaBNodes << std::endl;
	std::cout << "left subproblems in queue: " << left_subproblems << std::endl;
	std::cout << "====================================================================" << std::endl;
}

void outputFile(int n, int sol[], int val, int counterBaBNodes, int counterPrune, double timeTotal, int maxWorkers)
{
	std::ofstream file_output("solution.txt");

	file_output << "====================================================================" << std::endl;
	file_output << "OPTIMAL SOLUTION FOUND."  << std::endl;

	file_output << "value: " << val << std::endl;

	file_output << "stable set: (";
	int numberOfPrintedElements = 0;

	for(int i = 0; i < n; i++){
		if (sol[i] != 0){
			numberOfPrintedElements++;
			if (numberOfPrintedElements < val)
				file_output << i+1 << ", ";
			else
				file_output << i+1;
		}	
	}
	file_output << ")" << std::endl;

	file_output << "====================================================================" << std::endl;
	file_output << "ADDITIONAL INFORMATION" << std::endl;

	file_output << "B&B nodes: " << counterBaBNodes << std::endl;
	file_output << "pruned by bound: " << counterPrune << std::endl;
	file_output << "total time: " << timeTotal << std::endl;
	file_output << "max workers: " << maxWorkers << std::endl;
	file_output << "====================================================================" << std::endl;

	file_output.close();
}

void outputFile_approx(int n, int sol[], int val,int counterBaBNodes, double best_upper_bound, int left_subproblems)
{
	std::ofstream file_output("solution.txt");

	file_output << "====================================================================" << std::endl;
	file_output << "Time limit " << SECONDS <<  " exceeded.\n";
	file_output << "APPROXIMATE SOLUTION"  << std::endl;
	
	file_output << "value: " << val << std::endl;

	file_output << "stable set: (";
	int numberOfPrintedElements = 0;

	for(int i = 0; i < n; i++){
		if (sol[i] != 0){
			numberOfPrintedElements++;
			if (numberOfPrintedElements < val)
				file_output << i+1 << ", ";
			else
				file_output << i+1;
		}	
	}
	file_output << ")" << std::endl;

	file_output << "highest upper bound: " << best_upper_bound << std::endl;
	file_output << "left subproblems in queue: " << left_subproblems << std::endl; 

	file_output << "====================================================================" << std::endl;

	file_output.close();
}


// initialization
inline static void initalize_bound_computation(){
#if defined THETA
	g_lbndComp = new BoundComputationTheta();
	g_ubndComp = new BoundComputationTheta();
	//std::cout << "using theta" << std::endl;
#elif defined EXACTSUB
	g_lbndComp = new BoundComputationExactSubgraph();
	g_ubndComp = new BoundComputationExactSubgraph();
	std::cout << "using exact subgraph" << std::endl;
#elif defined BUNDLE
	g_lbndComp = new BoundComputationBundle();
	g_ubndComp = new BoundComputationBundle();
	std::cout << "using bundle method" << std::endl;
#else
	g_lbndComp = new BoundComputationTrivial();
	g_ubndComp = new BoundComputationTrivial();
	//std::cout << "using trivial heuristic" << std::endl;
#endif
}


int findFreeWorker(int *busyWorkers, int numbWorkers)
{
	for(int i = 1; i <= numbWorkers-1; i++)
	{
		if(busyWorkers[i] == 0)
		{
			return i;
		}
	}
}


// MAIN FUNCTION
int main(int argc, char **argv)
{
	double time_start_MAIN = MPI_Wtime();
	double time_begin = MPI_Wtime();

	// Everybody sees original graph

	// number of processes = master + slaves
	int numbWorkers;
	// rank of each process
	int rank;

	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &numbWorkers);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	// tag of message
	int tag;

	//variables for input information
	int m_n, m_m;
	std::vector<Edge> m_edges;
	std::vector<Vertex> m_vertices;

	//solution vector
	int *sol;

	// received solution vector and lower bound
	int *recv_sol;
	int bound_lb;

	//double eps = std::numeric_limits<double>::epsilon();
	double eps = OPTEPS; //defined in makefile
	
	//read params
	read_params();

	//initialize object for bound computation
	initalize_bound_computation(); 

	
	//compute lower bound
	g_lowerBound = 0;

	//helper variable for current subproblem
	Subproblem p;
	

	if (rank == 0)
	{
		//only master reads the file and sends info to slaves
		//read input
		read_input(argv, m_n, m_m, m_vertices, m_edges);

		// extract vectors

		// vertices (Original graph: local numbers are equal to global numbers)
		std::vector<int> vertices_numbers;

		for (std::vector<Vertex>::iterator it = m_vertices.begin() ; it != m_vertices.end(); ++it)
		{
			vertices_numbers.push_back(it->gNumber);
		}

		// edges
		std::vector<int> edges_left;
		std::vector<int> edges_right;

		for (std::vector<Edge>::iterator it = m_edges.begin() ; it != m_edges.end(); ++it)
		{
			edges_left.push_back(it->vertex1);
			edges_right.push_back(it->vertex2);
		}
/*
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&m_n,1,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&m_m,1,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(vertices_numbers.data(),m_n,MPI_INT,0,MPI_COMM_WORLD);
		
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(edges_left.data(),m_m,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(edges_right.data(),m_m,MPI_INT,0,MPI_COMM_WORLD);
*/
		for(int i = 1; i <= numbWorkers-1; i++)
		{
			MPI_Send(&m_n, 1, MPI_INT, i, 1, MPI_COMM_WORLD);
			MPI_Send(&m_m, 1, MPI_INT, i, 2, MPI_COMM_WORLD);
			MPI_Send(vertices_numbers.data(),m_n,MPI_INT,i,3,MPI_COMM_WORLD);
			MPI_Send(edges_left.data(),m_m,MPI_INT,i,4,MPI_COMM_WORLD);
			MPI_Send(edges_right.data(),m_m,MPI_INT,i,5,MPI_COMM_WORLD);

		}

		
		//initialize variable(s) depending on number of vertices/edges
		sol = new int[m_n];
		recv_sol = new int[m_n];

		// best upper bound
		double best_upper_bound = 0;

		//create graph
		Graph g = Graph(m_vertices, m_edges);
		
		//#############
		g_lbndComp->getLowerBound(&g, g_lowerBound, sol);

		//create vector of open subproblems and insert first subproblem
		std::list<Subproblem> openProblems;
		openProblems.push_back(Subproblem(g));

		int *busyWorkers;
		busyWorkers = new int[numbWorkers];
	
		//master is busy
		busyWorkers[0] = 1;

		for (int i = 1; i<numbWorkers; ++i)
		{
			busyWorkers[i] = 0;
		}

		int firstTime = 0; // send original problem to first slave
		int numbFreeWorkers = numbWorkers-1;
		int iterations;	// number of problems master has to send
		int freeWorker;	// rank of free worker
		int source;	// rank of source
		int active;
		int over;	// tag to FINISH

		int n_bbNodes = 1;
		int m_counterPrune = 0;
		int maxWorkers = 0;

		//std::cout << "Ok before main while loop --> master" << std::endl;

		while(1)
		{
#ifdef DEBUG
			std::cout << "Waiting" << std::endl;
#endif
			
			
			if(firstTime == 0)
			{
				tag = SEND_TO_WORKERS;
				firstTime = 1;
			}
			else
			{
				double aaa = MPI_Wtime();
				MPI_Recv(&source, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&tag, 1, MPI_INT, source, 1, MPI_COMM_WORLD, &status);
#ifdef DEBUG
				std::cout << "Initial receive time: " << MPI_Wtime() - aaa << std::endl; 	
#endif	
			}
		

			if(tag == GET_CHILDREN)
			{	
#ifdef DEBUG
				std::cout << "GET_CHILDREN" << std::endl;
#endif
				// receive subproblem
				// receive branchign variable

				double time_start = MPI_Wtime();

				int branching_variable;
				MPI_Recv(&branching_variable, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);

				int n_vertices_subgraph;
				MPI_Recv(&n_vertices_subgraph, 1, MPI_INT, source, 3, MPI_COMM_WORLD, &status);
				
				std::vector<int> globalNumbersVertices;
				globalNumbersVertices.resize(n_vertices_subgraph);
				MPI_Recv(globalNumbersVertices.data(), n_vertices_subgraph, MPI_INT, source, 4, MPI_COMM_WORLD, &status);
				
				double bound_ub;
				int depth;
				int const_value;
				int n_stable_set;

				MPI_Recv(&bound_ub, 1, MPI_DOUBLE, source, 5, MPI_COMM_WORLD, &status);
				

				MPI_Recv(&depth, 1, MPI_INT, source, 6, MPI_COMM_WORLD, &status);
				MPI_Recv(&const_value, 1, MPI_INT, source, 7, MPI_COMM_WORLD, &status);
				MPI_Recv(&n_stable_set, 1, MPI_INT, source, 8, MPI_COMM_WORLD, &status);

				std::vector<int> stable_set;
				stable_set.resize(n_stable_set);
				MPI_Recv(stable_set.data(), n_stable_set, MPI_INT, source, 9, MPI_COMM_WORLD, &status);
	
				std::set<int> stable_set_converted(stable_set.begin(), stable_set.end());

				// decode subproblem
				p = Subproblem(globalNumbersVertices, m_edges, depth, bound_ub, const_value, stable_set_converted);

				// branch
				openProblems.push_back(Subproblem(p,branching_variable,1));
				openProblems.push_front(Subproblem(p,branching_variable,0));
				
				double time_end = MPI_Wtime();
#ifdef DEBUG
				std::cout << "Time: " << time_end - time_start << std::endl;
#endif
				busyWorkers[source] = 0;
				numbFreeWorkers++;

				n_bbNodes += 2;
				
				// before sending problems to others wait for all the messages
				int flag;

				MPI_Iprobe(MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &flag, &status);
		
				if (flag == 0)
				{
					tag = SEND_TO_WORKERS;
#ifdef DEBUG
					std::cout << "Sending to workers. No pending message " << std::endl;
#endif
				}
				else 
				{
					tag = -1;
#ifdef DEBUG
					std::cout << "Pending message." << std::endl;
#endif
				}


				if ((MPI_Wtime() - time_begin) > SECONDS)
				{
					over = 1;
					for(int i = 1; i <= numbWorkers-1; i++)
					{
						MPI_Send(&over, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
					}

					// delete all subproblems in queue which have lower upper bound than g_lowerBound and
					// output the highest upper_bound from the rest --> calculate gap
	
					int queue_length = openProblems.size();

					// create iterator for openProblems
					std::list<Subproblem>::iterator it;
				
					for (it = openProblems.begin(); it != openProblems.end(); it++){
						best_upper_bound = (it->upperBound() > best_upper_bound) ? it->upperBound() : best_upper_bound;				
					}

					print_solution_expired(m_n, sol, g_lowerBound, n_bbNodes, best_upper_bound, queue_length);
					outputFile_approx(m_n, sol, g_lowerBound, n_bbNodes, best_upper_bound, queue_length);

					break;

				}

			}	
			if(tag == FATHOM)
			{	
#ifdef DEBUG
				std::cout << "FATHOM" << std::endl;
#endif
				busyWorkers[source] = 0;
				numbFreeWorkers++;
				m_counterPrune++;
				
				tag = SEND_TO_WORKERS;
			}
			if(tag == NEW_VALUE)
			{
#ifdef DEBUG
				std::cout << "NEW_VALUE" << std::endl;
#endif
				// receive lower bound and vector
				MPI_Recv(&bound_lb, 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
				MPI_Recv(recv_sol, m_n, MPI_INT, source, 3, MPI_COMM_WORLD, &status);
				if(bound_lb > g_lowerBound)
				{
					g_lowerBound = bound_lb;
					for(int i = 0; i < m_n; i++)
					{
						sol[i]= recv_sol[i];
					}
				}

				MPI_Send(&g_lowerBound, 1, MPI_INT, source, 0, MPI_COMM_WORLD);
			}

			if(tag == SEND_TO_WORKERS)
			{
#ifdef DEBUG
				std::cout << "SEND_TO_WORKERS" << std::endl;
#endif
				over = 0;
				if(openProblems.empty() && (numbFreeWorkers == numbWorkers-1))
				{
					double time_end = MPI_Wtime();

					print_solution(m_n, sol, g_lowerBound, n_bbNodes, m_counterPrune, time_end-time_start_MAIN, maxWorkers);
					outputFile(m_n, sol, g_lowerBound, n_bbNodes, m_counterPrune, time_end-time_start_MAIN, maxWorkers);

					tag = DONE;
				}
				else
				{
					int itemCount = openProblems.size();
					iterations = (itemCount < numbFreeWorkers) ? itemCount : numbFreeWorkers;
					active = 1;
#ifdef DEBUG
					std::cout << "Iterations: " << iterations << std::endl;
#endif
					
					for(int k = 1; k <= iterations; k++)
					{
						// take element from list
						p = openProblems.back();
						openProblems.pop_back();

						// sned problem only if upper bound is greater than lower bound
						if (p.upperBound() > g_lowerBound)
						{
							// find free worker and set it to busy
							freeWorker = findFreeWorker(busyWorkers, numbWorkers);
#ifdef DEBUG
							std::cout << "Freeworker: " << freeWorker << std::endl;
#endif
							busyWorkers[freeWorker] = 1;
							numbFreeWorkers--;

							// need only global numbers of vertices
							m_vertices = p.graph().vertices();
							std::vector<int> globalNumbersVertices;

							for (std::vector<Vertex>::iterator it = m_vertices.begin() ; it != m_vertices.end(); ++it)
							{
								globalNumbersVertices.push_back(it->gNumber);
							}

							int n_vertices_subgraph = globalNumbersVertices.size();
							//std::cout << "rank " << rank << "number of subgraph vertices " << n_vertices_subgraph << std::endl;
							int depth = p.depth();
							double upper_bound = p.upperBound();
							int const_value = p.get_constValue();	
							//std::cout << "rank " << rank << "ok " << i << std::endl;
	
							// convert set stable set to vector
							std::set<int> temp_stable_set = p.stableSet();
							std::vector<int> stable_set(temp_stable_set.begin(), temp_stable_set.end());
							//stable_set.assign(temp_stable_set.begin(), temp_stable_set.end());

							//stable_set.assign(p.stableSet().begin(), p.stableSet().end());
							//std::vector<int> stable_set(p.stableSet().begin p.stableSet().end);

							int n_stable_set = stable_set.size();
						
							//std::cout << "rank " << rank << "ok_next" << std::endl;

							//std::cout << "rank " << rank << "before sending problem to other" << std::endl;

							// send subproblem to worker
							MPI_Send(&over, 1, MPI_INT, freeWorker, 0, MPI_COMM_WORLD);
							MPI_Send(&active, 1, MPI_INT, freeWorker, 1, MPI_COMM_WORLD);
							MPI_Send(&g_lowerBound, 1, MPI_INT, freeWorker, 2, MPI_COMM_WORLD);
							MPI_Send(&n_vertices_subgraph, 1, MPI_INT, freeWorker, 3, MPI_COMM_WORLD);
							MPI_Send(globalNumbersVertices.data(), n_vertices_subgraph, MPI_INT, freeWorker, 4, MPI_COMM_WORLD);
							MPI_Send(&upper_bound, 1, MPI_DOUBLE, freeWorker, 5, MPI_COMM_WORLD);
							MPI_Send(&depth, 1, MPI_INT, freeWorker, 6, MPI_COMM_WORLD);
							MPI_Send(&const_value, 1, MPI_INT, freeWorker, 7, MPI_COMM_WORLD);
							MPI_Send(&n_stable_set, 1, MPI_INT, freeWorker, 8, MPI_COMM_WORLD);
							MPI_Send(stable_set.data(), n_stable_set, MPI_INT, freeWorker, 9, MPI_COMM_WORLD);
						}						
					}
					maxWorkers = (maxWorkers >=  numbWorkers - 1 - numbFreeWorkers) ? maxWorkers :  numbWorkers - 1 - numbFreeWorkers;
#ifdef DEBUG
					std::cout << "Number of free workers: " << numbFreeWorkers << std::endl;
#endif
				}
			}

			if(tag == DONE)
			{
				over = 1;
				for(int i = 1; i <= numbWorkers-1; i++)
				{
					MPI_Send(&over, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				}
				
				break;
			}

			
		}	
		delete [] busyWorkers;
	}

	else
	{
		int over; // tag for FINISH
		int active;
/*
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&m_n,1,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&m_m,1,MPI_INT,0,MPI_COMM_WORLD);
*/
		// receive graph info
		// vertices
		std::vector<int> vertices_numbers;
		
		// edges
		std::vector<int> edges_left;
		std::vector<int> edges_right;
/*
		// resize
		vertices_numbers.resize(m_n);
		edges_left.resize(m_m);
		edges_right.resize(m_m);

		MPI_Barrier(MPI_COMM_WORLD);
		
		MPI_Bcast(vertices_numbers.data(),m_n,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Bcast(edges_left.data(),m_m,MPI_INT,0,MPI_COMM_WORLD);

		MPI_Barrier(MPI_COMM_WORLD);

		MPI_Bcast(edges_right.data(),m_m,MPI_INT,0,MPI_COMM_WORLD);
	*/	

		MPI_Recv(&m_n, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&m_m, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);

		// resize
		vertices_numbers.resize(m_n);
		edges_left.resize(m_m);
		edges_right.resize(m_m);

		MPI_Recv(vertices_numbers.data(),m_n,MPI_INT,0,3,MPI_COMM_WORLD,&status);
		MPI_Recv(edges_left.data(),m_m,MPI_INT,0,4,MPI_COMM_WORLD,&status);
		MPI_Recv(edges_right.data(),m_m,MPI_INT,0,5,MPI_COMM_WORLD,&status);


		for(int i = 0; i < m_n; ++i)
		{
			m_vertices.emplace_back(vertices_numbers.at(i),vertices_numbers.at(i));
		}

		for (int i = 0; i < m_m; ++i)
		{
			m_edges.emplace_back(edges_left.at(i), edges_right.at(i));
		}

		//initialize variable(s) depending on number of vertices/edges
		sol = new int[m_n];
		recv_sol = new int[m_n];

		//create graph
		Graph g = Graph(m_vertices, m_edges);

		//std::cout << "Slave ready to go." << std::endl;

		while(1)
		{
			
			MPI_Recv(&over, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
			if(over == 1)
			{
				break;
			}
			else
			{
				MPI_Recv(&active, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
				if(active == 1)
				{
					

					// receive subproblem with data (upper bound) + best lower bound (bound_lb)
					// dont need solution vector
					
					double time_start = MPI_Wtime();

					MPI_Recv(&bound_lb, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);

					int n_vertices_subgraph;
					MPI_Recv(&n_vertices_subgraph, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
					
					std::vector<int> globalNumbersVertices;
					globalNumbersVertices.resize(n_vertices_subgraph);
					MPI_Recv(globalNumbersVertices.data(), n_vertices_subgraph, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);
					
					double bound_ub;
					int depth;
					int const_value;
					int n_stable_set;

					MPI_Recv(&bound_ub, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &status);
					MPI_Recv(&depth, 1, MPI_INT, 0, 6, MPI_COMM_WORLD, &status);
					MPI_Recv(&const_value, 1, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
					MPI_Recv(&n_stable_set, 1, MPI_INT, 0, 8, MPI_COMM_WORLD, &status);

					std::vector<int> stable_set;
					stable_set.resize(n_stable_set);
					MPI_Recv(stable_set.data(), n_stable_set, MPI_INT, 0, 9, MPI_COMM_WORLD, &status);
					//std::cout << "End receiving" << std::endl;
					
					double time_end = MPI_Wtime();
					//std::cout << "Time receiving: " << time_end-time_start << " rank " << rank << std::endl;

					if(bound_ub < bound_lb + 1 - eps)
					{
						tag = FATHOM;
						MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
						MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
						//MPI_Send(&bound_lb, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
						
						
					}

					else
					{
						// decode the subproblem
						// transform stable set vector to set when constructing subproblem
						std::set<int> stable_set_converted(stable_set.begin(), stable_set.end());

						time_start = MPI_Wtime();
						double time_start_2 = MPI_Wtime();
						p = Subproblem(globalNumbersVertices, m_edges, depth, bound_ub, const_value, stable_set_converted);
						time_end = MPI_Wtime();
						//std::cout << "Time constructing problem: " << time_end-time_start << std::endl;
/*
						std::cout << "constructed subgraph vertices\n";
						std::vector<Vertex> temp_vertices = p.graph().vertices();
						for (auto itv : temp_vertices)
						{
							std::cout << itv.gNumber << " and " << itv.lNumber << std::endl;
						}
			
				
						std::cout << "constructed subgraph edges\n";
						// print edges
						std::vector<Edge> temp_edges = p.graph().edges();
						for (auto it : temp_edges)
						{
							std::cout << it.vertex1 << " and " << it.vertex2 << std::endl;
						}
*/
						// calculate upper and lower bound
						
						time_start = MPI_Wtime();
						p.computeUpperBound();
						time_end = MPI_Wtime();
						//std::cout << "Time computing upper bound: " << time_end-time_start << std::endl;

						time_start = MPI_Wtime();
						p.computeLowerBound(sol, m_n);
						time_end = MPI_Wtime();
						//std::cout << "Time computing lower bound: " << time_end-time_start << std::endl;

						// check if calculated lower bound is better than received one
						//if (g_lowerBound > bound_lb + 1 - eps)
						//	bound_lb = g_lowerBound;

						if (rank == 2)
						{
#ifdef DEBUG
							std::cout << "##############################################Cas: " << MPI_Wtime()-time_start_2 << std::endl;
#endif
						}


						//if((p.upperBound() > bound_lb) && (g_lowerBound > bound_lb + 1 - eps))
						if (g_lowerBound > bound_lb + 1 - eps)
						{
							tag = NEW_VALUE;
							MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
							MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
							MPI_Send(&g_lowerBound, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
							MPI_Send(sol, m_n, MPI_INT, 0, 3, MPI_COMM_WORLD);

							MPI_Recv(&bound_lb, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
						}

						if (p.upperBound() < bound_lb + 1 - eps)
						{
							tag = FATHOM;
							MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
							MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
						}
						else
						{
							// send problem back for branching
							tag = GET_CHILDREN;
							int branching_variable = p.nextBranchVar();
							double upper_bound = p.upperBound();

							MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
							MPI_Send(&tag, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
							MPI_Send(&branching_variable, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
							MPI_Send(&n_vertices_subgraph, 1, MPI_INT, 0, 3, MPI_COMM_WORLD);
							MPI_Send(globalNumbersVertices.data(), n_vertices_subgraph, MPI_INT, 0, 4, MPI_COMM_WORLD);
							MPI_Send(&upper_bound, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
							MPI_Send(&depth, 1, MPI_INT, 0, 6, MPI_COMM_WORLD);
							MPI_Send(&const_value, 1, MPI_INT, 0, 7, MPI_COMM_WORLD);
							MPI_Send(&n_stable_set, 1, MPI_INT, 0, 8, MPI_COMM_WORLD);
							MPI_Send(stable_set.data(), n_stable_set, MPI_INT,0, 9, MPI_COMM_WORLD);
							
						}
						
					}
				   
				}
			}
		}
	}



	//delete variables on heap
	delete [] sol;
	delete [] recv_sol;
	delete g_lbndComp;
	delete g_ubndComp;
	

	MPI_Finalize();
	return 0;
}


