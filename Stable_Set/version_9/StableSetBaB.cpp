// Parallel Stable Set solver based on code from Melanie Siebenhofer
// Uses MPI
// Written by Timotej Hrga
// timotej.hrga@gmail.com


/* MASTER (load coordinator): 
- at beginning broadcasts original graph to the workers --> when workers are encoding subproblem send to them, they need to know original graph
- branch original problem and send subproblem to each worker --> that way everybody starts working immedeately
- has bool array of free/busy workers
- has queue of subproblems just for use at the beggining. Later he does not store any subproblem
- manages best found optimal value and solution and distributes them when neccessary
- receives/sends info specifying which workers are free
*/

/* WORKER:
- has its own local queue of subproblems
- after branching ask master if any other worker is free and if better lower bound has been found --> send him the subproblem
- otherwise subproblem remains in workers queue
*/


// NOTE: global lower bound is exchanged when message NEW_VALUE is sent of when worker asked master for free slaves

#include <iostream>
#include <fstream>
#include <queue>
#include <algorithm>
#include <mpi.h>

#include "Subproblem.h"
#include "ReadAndOutput.h"


// MESSAGES
enum Message{
    SEND_FREEWORKERS,  	// some worker asked if there is any free worker --> send corresponding rank
    IDLE,             	// worker is free, his local queue of subproblems is empty
    NEW_VALUE         	// better lower bound found
};

// TAGS in MPI messages
enum Tags{
    OVER,                       // when sending info to finish
    MESSAGE,                    // for type of message
    TOKEN_FREEWORKER,           // when receiving/sending rank of free worker
    NUM_FREE_WORKERS,
    GLOBAL_LOWER_BOUND,         // when best lower bound is received/send
    RECEIVED_SOLUTION,          // best global solution vector
    NUM_VERTICES_SUBGRAPH,      // tags when transmiting subproblem to another worker
    GLOBAL_SUBVERTICES,
    UPPER_BOUND,
    CONST_VALUE,
    NUM_STABLE_SET,
    STABLE_SET
};



// global variables
int g_lowerBound = 1;                // global lower bound
Params params;                       // read parameters from file
BoundComputationTheta* g_bndComp;


// class to compare two subproblems
// needed for priority queue of open subproblems
// sort subproblems in priority queue by decreasing upper bound.
// The one with the HIGHEST (worst) upper bound is considered FIRST
class subproblemComp{
public:
    bool operator() (const Subproblem& lhs, const Subproblem& rhs){
        return lhs.upperBound() < rhs.upperBound();
    }
};


// when initial branching is done to distribute subproblems to the workers
// we branch on vertex with lowest/highest degree (params.
inline int get_vertex(Graph &g){
	std::vector<int> degrees(g.n(),0);	// vector of degrees of vertices
	
	for (std::vector<Edge>::iterator it = g.edges().begin(); it != g.edges().end(); ++it){
		++(degrees[it->vertex1-1]);
		++(degrees[it->vertex2-1]);
	}
	
	if (params.initialBranchingStrategy) // select the vertex with maximal degree	
		return ( std::max_element(degrees.begin(),degrees.end()) - degrees.begin() ) + 1;
	else
		return ( std::min_element(degrees.begin(),degrees.end()) - degrees.begin() ) + 1;
}



int main(int argc, char *argv[]) {
  
    /*******************************************************
     *********** BRANCH & BOUND: PARALLEL ALGORITHM ********
     ******************************************************/
    
    // Everybody sees original graph: master reads and distributes to workers via broadcast
    
    // number of processes = master + slaves
    int numbWorkers;
    
    // rank of each process: from 0 to numWorkers-1
    int rank;
    
    // MPI Start: start parallel environment
    MPI_Init(NULL,NULL);
    
    // get number of proccesses and corresponding ranks
    MPI_Comm_size(MPI_COMM_WORLD, &numbWorkers);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Status status;
    
    if (argc < 5){
    	if (rank == 0)
    		std::cout << "Not enough arguments to the program." << std::endl;
    		
    	MPI_Finalize(); 
    	return EXIT_SUCCESS;	
    }
    
    // start clock
    double m_timeTotal = MPI_Wtime();
    
    // time constants for temp output during the algorithm
	int TIME_INFO = std::stoi(argv[3]); 		// output first info after MINUTES minutes
	int TIME_TEMP = std::stoi(argv[4]); 			// output further info every hour
    int add_limit = TIME_TEMP;
    
    // tag of message
    Message message;
    
    // variables for input information
    int m_n, m_m;
    std::vector<Edge> m_edges;
    std::vector<Vertex> m_vertices;
    
    // solution vector
    // for MASTER it holds global solution
    // for SLAVE it holds  solution of the suproblem that is send to the master
    int *sol;
    
    // received lower bound --> to compare with the global bound
    int bound_lb = 0;
    
    // threshhold for comparing upper and lower bound
    double eps = 1e-5; 
    
    // initialize object for bound computation
    g_bndComp = new BoundComputationTheta();

     // read params
    read_params(params);
    
    // helper varibles to extract, send and receive original graph
    std::vector<int> vertices_numbers;
    std::vector<int> edges_left;
    std::vector<int> edges_right;
    
    // helper variable for current subproblem
    Subproblem p;
        
    // priority queue containing subproblems --> each of the workers has its local
    std::priority_queue<Subproblem, std::vector<Subproblem>, subproblemComp> openProblemsSorted;
    
    // tag to FINISH set to FALSE
    int over = 0; 

    /******************** MASTER PROCESS ********************/
    if (rank == 0){

        // only master counts the nodes in the tree 
        int counterBaBNodes = 1;
         
        // only master reads the file and sends the graph to slaves
        // read input
        int read_flag = read_input(argv, m_n, m_m, m_vertices, m_edges);
        
        // check if reading succeeded
        if (read_flag == -1){	// error
        	
			// send -1 to workers --> they will know that they have to finish
	    	over = -1;
	    	MPI_Bcast(&over, 1, MPI_INT, 0, MPI_COMM_WORLD);
	    	goto FINISH;
        }
        
    
        // create graph
        Graph g { m_vertices, m_edges };
        
        // compute density of the graph
    	double density = 2 * static_cast<double>(m_m) / (m_n*(m_n-1));
        
        if(m_n <= params.minSizeGraphSizeForFullEnum){  // FULL ENUMERATION 
        
	        std::vector<int> stableSet = g.getMaximumStableSet(m_n, g_lowerBound);
	       	g_lowerBound = stableSet.size();
	    
	        int *solution = new int[m_n](); 	// only to transform stableSet to array and print
	        
	        for (std::vector<int>::iterator it = stableSet.begin(); it != stableSet.end(); ++it)
	        	solution[*it-1] = 1;
	        
	    	print_solution(m_n, m_m, solution, g_lowerBound, counterBaBNodes, MPI_Wtime() - m_timeTotal, density, 1);
	    	outputFile(argv, m_n, m_m, solution, g_lowerBound, counterBaBNodes, MPI_Wtime() - m_timeTotal, density, 1);
	    	
	    	delete [] solution;
	    	
	    	// send -1 to workers --> they will know that they have to finish
	    	over = -1;
	    	MPI_Bcast(&over, 1, MPI_INT, 0, MPI_COMM_WORLD);	
	    	
	    	goto FINISH;    	
	    }
        
        
        // initialize solution vector
        sol = new int[m_n];
        
        // run heuristic to compute lower bound       
        g_bndComp->getLowerBound(&g, g_lowerBound, sol);           
    
        // send lower bound to all processes
        MPI_Bcast(&g_lowerBound, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        // prepare data for sending original graph to workers
        // vertices (Original graph: local numbers are equal to global numbers)
        for (std::vector<Vertex>::iterator it = m_vertices.begin(); it != m_vertices.end(); ++it)
        {
            vertices_numbers.push_back(it->gNumber);
        }

        for (std::vector<Edge>::iterator it = m_edges.begin(); it != m_edges.end(); ++it)
        {
            edges_left.push_back(it->vertex1);
            edges_right.push_back(it->vertex2);
        }

        // broadcast the graph
        MPI_Bcast(&m_n,1,MPI_INT,0,MPI_COMM_WORLD);                         // number of vertices
        MPI_Bcast(&m_m,1,MPI_INT,0,MPI_COMM_WORLD);                         // number of edges
        MPI_Bcast(vertices_numbers.data(),m_n,MPI_INT,0,MPI_COMM_WORLD);    // vertices 
        MPI_Bcast(edges_left.data(),m_m,MPI_INT,0,MPI_COMM_WORLD);          // edges: left and right vertex
        MPI_Bcast(edges_right.data(),m_m,MPI_INT,0,MPI_COMM_WORLD);
        
        // solution received solution from the worker
        int *recv_sol = new int[m_n];
        
        
        // put original problem into priority queue and branch as long as the number of MPI processors
        openProblemsSorted.emplace(g);
        
        // branch on vertex with highest/lowest degree --> specified in params
      	int nextBranchVar;
        for (int i = 1; i < numbWorkers-1; ++i){	// branch numbWorkers-2 times to get numbWorkers-1 subproblems
        
            //get next subproblem from list
            p = openProblemsSorted.top();
            openProblemsSorted.pop();
            
            // branch on vertex with the highest degree
            nextBranchVar = get_vertex(p.graph());
        
            openProblemsSorted.emplace(p,nextBranchVar,0);
            openProblemsSorted.emplace(p,nextBranchVar,1);
            
            counterBaBNodes += 2;
         }   
            
         // send subpoblems to the workers  
         for (int i = 1; i < numbWorkers; ++i){   
            p = openProblemsSorted.top();
            openProblemsSorted.pop();                  

            // extract data and send graph to the corresponding worker
            m_vertices = p.graph().vertices();
            std::vector<int> globalNumbersVertices;

            for (std::vector<Vertex>::iterator it = m_vertices.begin() ; it != m_vertices.end(); ++it)
            {
                globalNumbersVertices.push_back(it->gNumber);
            }

            int n_vertices_subgraph = globalNumbersVertices.size();
            double upper_bound = p.upperBound();
            int const_value = p.getConstValue();

            // convert set stable set to vector
            std::set<int> temp_stable_set = p.stableSet();
            std::vector<int> stable_set(temp_stable_set.begin(), temp_stable_set.end());

            int n_stable_set = stable_set.size();

            // send subproblem to worker
            MPI_Send(&n_vertices_subgraph, 1, MPI_INT, i, NUM_VERTICES_SUBGRAPH, MPI_COMM_WORLD);
            MPI_Send(globalNumbersVertices.data(), n_vertices_subgraph, MPI_INT, i, GLOBAL_SUBVERTICES, MPI_COMM_WORLD);
            MPI_Send(&upper_bound, 1, MPI_DOUBLE, i, UPPER_BOUND, MPI_COMM_WORLD);
            MPI_Send(&const_value, 1, MPI_INT, i, CONST_VALUE, MPI_COMM_WORLD);
            MPI_Send(&n_stable_set, 1, MPI_INT, i, NUM_STABLE_SET, MPI_COMM_WORLD);
            MPI_Send(stable_set.data(), n_stable_set, MPI_INT, i, STABLE_SET, MPI_COMM_WORLD);  
        }
        
        
        // array of busy workers: 0 = free, 1 = busy
        bool *busyWorkers = new bool[numbWorkers];
        
        for (int i = 0; i < numbWorkers; ++i){ // all workers + master are busy
            busyWorkers[i] = 1;
        }
        
        int numbFreeWorkers = 0;
        int source;
       
		bool print_approximate_sol = 0;			// is approximate solution already printed
		bool print_initial_info = 0;
	 
        /************* MAIN LOOP for master **************/
        do{

            // Check if time limit reached
            if ((MPI_Wtime() - m_timeTotal) > params.maxExecutionTime){
           		if (!print_approximate_sol){
           			// print approximate solution
	            	print_solution(m_n, m_m, sol, g_lowerBound, counterBaBNodes, MPI_Wtime() - m_timeTotal, density, 0);
	            	outputFile(argv, m_n, m_m, sol, g_lowerBound, counterBaBNodes, MPI_Wtime() - m_timeTotal, density, 0);		
	                
	                print_approximate_sol = 1;
                }      	                 
            }
            
            // Check if time for outputing info after start
            if (!print_initial_info){
	            if ((MPI_Wtime() - m_timeTotal) > TIME_INFO){
	            	outputFile_after_start(argv, m_n, m_m, sol, g_lowerBound,counterBaBNodes, MPI_Wtime() - m_timeTotal, density, 1);
	            	print_initial_info = 1;                 
	            }
            }
            
            // Check if time for outputing info during the algorithm
            if ((MPI_Wtime() - m_timeTotal) > TIME_TEMP){
            	outputFile_after_start(argv, m_n, m_m, sol, g_lowerBound,counterBaBNodes, MPI_Wtime() - m_timeTotal, density, 0);
            	TIME_TEMP += add_limit;                 
            }
            
            // wait for messages
            // extract source from status.MPI_SOURCE !!!
            MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, MESSAGE, MPI_COMM_WORLD, &status);
            source = status.MPI_SOURCE;

            if (message == IDLE){                       
                busyWorkers[source] = 0;
                numbFreeWorkers++;          
            }

            else if (message == SEND_FREEWORKERS){
                
                int workers_request; 	// get number of requested workers             
                MPI_Recv(&workers_request, 1, MPI_INT, source, TOKEN_FREEWORKER, MPI_COMM_WORLD, &status);
						    
				// compute number of freeworkers
				int num_workers_available = (workers_request < numbFreeWorkers) ? workers_request : numbFreeWorkers;
				
				std::vector<int> available_workers;
				
				for(int i = 1, j = num_workers_available; (i < numbWorkers) && j; ++i)    // master has rank 0 and is not considered
				{
    				if(busyWorkers[i] == 0){ // is free
        				available_workers.push_back(i);
        				--j;
        				busyWorkers[i] = 1; // set to busy
                    	numbFreeWorkers--;
    				}
				}
	    
	    		// worker branched subproblem in lcoal queue --> add 2 bab nodes
	    		counterBaBNodes += 2;	
						    
			    MPI_Send(&num_workers_available, 1, MPI_INT, source, NUM_FREE_WORKERS, MPI_COMM_WORLD);			     
			    MPI_Send(available_workers.data(), num_workers_available, MPI_INT, source, TOKEN_FREEWORKER, MPI_COMM_WORLD);
                MPI_Send(&g_lowerBound, 1, MPI_INT, source, GLOBAL_LOWER_BOUND, MPI_COMM_WORLD);       
            }

            else if(message == NEW_VALUE)
            {
                // receive lower bound and vector
                MPI_Recv(&bound_lb, 1, MPI_INT, source, GLOBAL_LOWER_BOUND, MPI_COMM_WORLD, &status);
                MPI_Recv(recv_sol, m_n, MPI_INT, source, RECEIVED_SOLUTION, MPI_COMM_WORLD, &status);
                
                if(bound_lb > g_lowerBound)	// update if better lower bound found
                {
                    g_lowerBound = bound_lb;
                    for(int i = 0; i < m_n; i++)
                    {
                        sol[i]= recv_sol[i];
                    }
                }
                
                // send update information to the sender -> other will receive it when asked for freeWorker
                MPI_Send(&g_lowerBound, 1, MPI_INT, source, GLOBAL_LOWER_BOUND, MPI_COMM_WORLD);   
            }

        } while (numbFreeWorkers != numbWorkers-1);


		if (!print_approximate_sol){ // print optimal solution
	        	print_solution(m_n, m_m, sol, g_lowerBound, counterBaBNodes, MPI_Wtime() - m_timeTotal, density, 1);
	        	outputFile(argv, m_n, m_m, sol, g_lowerBound, counterBaBNodes, MPI_Wtime() - m_timeTotal, density, 1);
   		}
        	
        //and send over messages to the workers
        over = 1;
        for(int i = 1; i < numbWorkers; ++i)
        {
            MPI_Send(&over, 1, MPI_INT, i, OVER, MPI_COMM_WORLD);
        }

        delete [] busyWorkers;
        delete [] recv_sol;
    }
    
    /******************** WORKER PROCESS ********************/
    else
    {
        // receive lower bound
        MPI_Bcast(&g_lowerBound, 1, MPI_INT, 0, MPI_COMM_WORLD);
        
        if (g_lowerBound == -1)		// master did full enumeration --> finish
        	goto FINISH;
        
        // receive original graph info
        MPI_Bcast(&m_n,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&m_m,1,MPI_INT,0,MPI_COMM_WORLD);

        // resize before receiving
        vertices_numbers.resize(m_n);
        edges_left.resize(m_m);
        edges_right.resize(m_m);

        MPI_Bcast(vertices_numbers.data(),m_n,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(edges_left.data(),m_m,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(edges_right.data(),m_m,MPI_INT,0,MPI_COMM_WORLD);
    
        // construct m_vertices and m_edges
        for(int i = 0; i < m_n; ++i){
            m_vertices.emplace_back(vertices_numbers[i],vertices_numbers[i]);
        }

        for (int i = 0; i < m_m; ++i){
            m_edges.emplace_back(edges_left[i], edges_right[i]);
        }
        
        // solution vector
        sol = new int[m_n];   
   
        // create graph
        Graph g { m_vertices, m_edges };
        
        int flag;
        
        /************* MAIN LOOP for worker **************/
        do {			
            // RECEIVE original problem from master or subproblem from other worker          
            int n_vertices_subgraph;
            MPI_Recv(&n_vertices_subgraph, 1, MPI_INT, MPI_ANY_SOURCE, NUM_VERTICES_SUBGRAPH, MPI_COMM_WORLD, &status);

            std::vector<int> globalNumbersVertices;
            globalNumbersVertices.resize(n_vertices_subgraph);
            MPI_Recv(globalNumbersVertices.data(), n_vertices_subgraph, MPI_INT, MPI_ANY_SOURCE, GLOBAL_SUBVERTICES, MPI_COMM_WORLD, &status);

            double bound_ub;
            int const_value;
            int n_stable_set;		

            MPI_Recv(&bound_ub, 1, MPI_DOUBLE, MPI_ANY_SOURCE, UPPER_BOUND, MPI_COMM_WORLD, &status);
            MPI_Recv(&const_value, 1, MPI_INT, MPI_ANY_SOURCE, CONST_VALUE, MPI_COMM_WORLD, &status);
            MPI_Recv(&n_stable_set, 1, MPI_INT, MPI_ANY_SOURCE, NUM_STABLE_SET, MPI_COMM_WORLD, &status);

            std::vector<int> stable_set;
            stable_set.resize(n_stable_set);
            MPI_Recv(stable_set.data(), n_stable_set, MPI_INT, MPI_ANY_SOURCE, STABLE_SET, MPI_COMM_WORLD, &status);
        
        	// transform stable set vector to set when constructing subproblem
            std::set<int> stable_set_converted(stable_set.begin(), stable_set.end());

            // START LOCAL QUEUE
            // decode the subproblem and put it into queue
            openProblemsSorted.emplace(globalNumbersVertices, m_edges, bound_ub, const_value, stable_set_converted);

			
            while( !openProblemsSorted.empty() ){			
				
				// check if time limit reached
				if ((MPI_Wtime() - m_timeTotal) > params.maxExecutionTime){
					goto TIME_LIMIT_WORKER;
				}
				
				
                // get next subproblem from list
                p = openProblemsSorted.top();
                openProblemsSorted.pop();

                if(p.upperBound() - g_lowerBound + eps >= 1){
                    // we are not able to prune by inhereted upper bound
                    
                    if( p.graph().n() > params.minSizeGraphSizeForFullEnum ){
                        
                        // compute upper and lower bound          
                        p.computeUpperBound();

                        // save "old" lower bound
                        bound_lb = g_lowerBound;
                                
                        // compute lower bound
                        p.computeLowerBound(sol, m_n);
                               
                        // check if calculated lower bound is better than received one
                        // update info with master
                        if (g_lowerBound > bound_lb){
                            message = NEW_VALUE;
                            
                            MPI_Send(&message, 1, MPI_INT, 0, MESSAGE, MPI_COMM_WORLD);
                            MPI_Send(&g_lowerBound, 1, MPI_INT, 0, GLOBAL_LOWER_BOUND, MPI_COMM_WORLD);
                            MPI_Send(sol, m_n, MPI_INT, 0, RECEIVED_SOLUTION, MPI_COMM_WORLD);
                            
                            MPI_Recv(&g_lowerBound, 1, MPI_INT, 0, GLOBAL_LOWER_BOUND, MPI_COMM_WORLD, &status);
                        }
                        
                        // compare: check if we can prune -> if not, branch
                        if(p.upperBound() - g_lowerBound + eps >= 1){
                            
                            // branch --> add new subproblems into list
                            openProblemsSorted.emplace(p,p.nextBranchVar(),0);
                            openProblemsSorted.emplace(p,p.nextBranchVar(),1);

							// leave 1 problem for this worker and the rest is distributed
			      			int workers_request = openProblemsSorted.size() - 1;
			      				 
			      			// check if other subproblems can be send to free workers --> ask master
			      			message = SEND_FREEWORKERS;
			      			
			      			int num_free_workers;
			      			std::vector<int> free_workers;
			      			
						    MPI_Send(&message, 1, MPI_INT, 0, MESSAGE, MPI_COMM_WORLD);
						    MPI_Send(&workers_request, 1, MPI_INT, 0, TOKEN_FREEWORKER, MPI_COMM_WORLD);
						    
						    MPI_Recv(&num_free_workers, 1, MPI_INT, 0, NUM_FREE_WORKERS, MPI_COMM_WORLD, &status);
						    
						    free_workers.resize(num_free_workers);
						    
						    MPI_Recv(free_workers.data(), num_free_workers, MPI_INT, 0, TOKEN_FREEWORKER, MPI_COMM_WORLD, &status);
						    MPI_Recv(&g_lowerBound, 1, MPI_INT, 0, GLOBAL_LOWER_BOUND, MPI_COMM_WORLD, &status);

					      	if ( num_free_workers != 0) {// free workers found
						  
							  	for (int i = 0; i < num_free_workers; ++i){
							  		// get next subproblem from list
            					  	p = openProblemsSorted.top();
            					  	openProblemsSorted.pop();

								  	// extract data and send graph to the corresponding worker
									m_vertices = p.graph().vertices();
								  	std::vector<int> globalNumbersVertices;

								  	for (std::vector<Vertex>::iterator it = m_vertices.begin() ; it != m_vertices.end(); ++it){
								    	globalNumbersVertices.push_back(it->gNumber);
								  	}

								 	int n_vertices_subgraph = globalNumbersVertices.size();
								  	double upper_bound = p.upperBound();
								  	int const_value = p.getConstValue();

								  	// convert set stable set to vector
								  	std::set<int> temp_stable_set = p.stableSet();
								  	std::vector<int> stable_set(temp_stable_set.begin(), temp_stable_set.end());

								  	int n_stable_set = stable_set.size();

								  	// send subproblem to free worker
									MPI_Send(&over, 1, MPI_INT, free_workers[i], OVER, MPI_COMM_WORLD);
								  	MPI_Send(&n_vertices_subgraph, 1, MPI_INT, free_workers[i], NUM_VERTICES_SUBGRAPH, MPI_COMM_WORLD);
								  	MPI_Send(globalNumbersVertices.data(), n_vertices_subgraph, MPI_INT, free_workers[i], GLOBAL_SUBVERTICES, MPI_COMM_WORLD);
								 	MPI_Send(&upper_bound, 1, MPI_DOUBLE, free_workers[i], UPPER_BOUND, MPI_COMM_WORLD);
								  	MPI_Send(&const_value, 1, MPI_INT, free_workers[i], CONST_VALUE, MPI_COMM_WORLD);
								  	MPI_Send(&n_stable_set, 1, MPI_INT, free_workers[i], NUM_STABLE_SET, MPI_COMM_WORLD);
								  	MPI_Send(stable_set.data(), n_stable_set, MPI_INT, free_workers[i], STABLE_SET, MPI_COMM_WORLD);
								}
					    	}   			      
			   		    }     
                    }
                    else{
                        //do a full enumeration
                        p.doFullEnumeration(sol, m_n);
                    }   
                }
            }

            // no more problems in local queue --> send info to master
            TIME_LIMIT_WORKER:		// tag for time limit reached --> all workers become idle
            message = IDLE;

            MPI_Send(&message, 1, MPI_INT, 0, MESSAGE, MPI_COMM_WORLD);
        
            // wait for info: stop (from master) or receive new subproblem from other worker
            MPI_Recv(&over, 1, MPI_INT, MPI_ANY_SOURCE, OVER, MPI_COMM_WORLD, &status);

        } while (over != 1);
    }
 
    //delete variables on heap
    delete [] sol;
    
    FINISH:				// tag for goto --> if initial graph was small we did full enumeration
    delete g_bndComp;

    // MPI finish
    MPI_Finalize();
    
    return EXIT_SUCCESS;
}
