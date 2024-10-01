/*
 * BoundComputationExactSubgraph.cpp
 *
 *  Created on: 06.09.2017
 *      Author: mesiebenhofe
 */

#include "BoundComputationExactSubgraph.h"

#include <random>
//#include <algorithm>

using namespace mosek::fusion;
using namespace monty;

BoundComputationExactSubgraph::BoundComputationExactSubgraph() {
	m_rootNode = true;

}

BoundComputationExactSubgraph::~BoundComputationExactSubgraph() {
}


int BoundComputationExactSubgraph::getUpperBound(Graph* g, double& upperBound, int& branchVar, int constValueUpperBound){
	m_vectorSet = false;
	int n = g->n();

	if(n==0){
		branchVar = 0;
		return 0;
	}

	//create model and add constraints
	//===================================================================================

	Model::t M = new Model("sdp"); auto _M = finally([&]() { M->dispose(); });

	//var X ... (n+1)x(n+1) matrix
	//---------------------------------------
	Variable::t X = M->variable(Domain::inPSDCone(n+1));

	//generate matrix C for objective function
	//----------------------------------------
	//C = diag(0*I_{1}, I_{n})  ... (n+1)x(n+1)
	auto rows = new_array_ptr<int,1>(n);
	auto cols = new_array_ptr<int,1>(n);
	auto values = new_array_ptr<double,1>(n);
	for(int i = 0; i<n; i++){
		(*rows)[i] = i+1;
		(*cols)[i] = i+1;
		(*values)[i] = 1;
	}
	auto C = Matrix::sparse(n+1,n+1,rows,cols,values);

	//objective function:
	//-------------------------------------------
	M->objective(ObjectiveSense::Maximize, Expr::dot(C, X));

	//add constraints:
	//--------------------------------------------
	//constraint 1: X_{1,1} = 1
	//..........................
	M->constraint(X->index(0,0)->asExpr(), Domain::equalsTo(1.0));

	//constraints type 1: X_{i,i} = x_{i}
	//............................
	for(int i = 1; i<n+1; i++){
		//X_{0,i} = X_{i,0} = X_{i,i} <==>
		//-0.5 X_{0,i} - 0.5 X_{i,0} + X_{i,i} = 0
		rows = new_array_ptr<int,1>({0,i,i});
		cols = new_array_ptr<int,1>({i,0,i});
		values = new_array_ptr<double,1>({-0.5,-0.5,1});
		auto A = Matrix::sparse(n+1,n+1,rows,cols,values);
		M->constraint(Expr::dot(A,X),Domain::equalsTo(0.0));
	}

	//constraints type 2: X_{i,j} = 0 } \forall (i,j) \in E <==>
	//0.5 X_{i,j} + 0.5 X_{j,i} = 0
	//...............................
	std::vector<Edge> edges = g->edges();
	for (std::vector<Edge>::iterator it = edges.begin() ; it != edges.end(); it++ ){
		rows = new_array_ptr<int,1>({it->vertex1,it->vertex2});
		cols = new_array_ptr<int,1>({it->vertex2, it->vertex1});
		values = new_array_ptr<double,1>({0.5,0.5});
		auto A = Matrix::sparse(n+1,n+1,rows,cols,values);
		M->constraint(Expr::dot(A,X),Domain::equalsTo(0.0));
	}


	//add exact subgraph constraints
	//=================================================================================================================
	//create adjacency
	//--------------------------
	int adjacency[n][n];
	for(int i = 0; i<n; i++){
		for(int j = 0; j<n; j++){
			adjacency[i][j] = 0;
		}
	}
	for(std::vector<Edge>::iterator ite = edges.begin() ; ite != edges.end(); ++ite){
		adjacency[ite->vertex1-1][ite->vertex2-1] = 1;
		adjacency[ite->vertex2-1][ite->vertex1-1] = 1;
	}
	//----------------------------

	//random device for std::shuffle
	std::random_device rd;
	std::mt19937 g1(rd());

	for (std::map<int,int>::iterator it= params.numberExactSubgraphs.begin(); it!=params.numberExactSubgraphs.end(); ++it){
		int sizeOfSubgraph = it->first;

#ifdef DEBUG
		//print size of subgraphs
		if(params.printlevel>3){
			std::cout << "size of subgraphs: " << sizeOfSubgraph << std::endl;
		}
#endif

		if(sizeOfSubgraph < n){
			std::vector<int>indicesGraph;
			for(int i = 0; i < n; i++){
				indicesGraph.push_back(i+1);
			}

#ifdef DEBUG
			//print number of subgraphs for this subgraph-size
			if(params.printlevel>3){
				std::cout << "add " << it->second << " subgraph constraints" << std::endl;
			}
#endif
			//add subgraph constraints it->second times
			for(unsigned iterSubgraphs = 0; iterSubgraphs < it->second; iterSubgraphs++){

				//shuffle indices
				std::shuffle(indicesGraph.begin(), indicesGraph.end(), g1);

				//subgraph = first sizeSubgraph entries in indices
				std::vector<int>Subgraph (indicesGraph.begin(),indicesGraph.begin() + sizeOfSubgraph) ;

#ifdef DEBUG
				if(params.printlevel>3){
					std::cout << "Subgraph: ";
					for(auto s: Subgraph){
						std::cout << s << " ";
					}
					std::cout << std::endl;
				}
#endif

				std::vector<std::vector<int>> stableSets = getStableSets(Subgraph, &adjacency[0][0], n);
				int numberOfStableSets = stableSets.size();

				//add variable lambda of size numberOfStableSets
				Variable::t lambda = M->variable(Domain::greaterThan(0,numberOfStableSets));

				//constraint \sum{lambda_{i}} = 1
				M->constraint(Expr::sum(lambda),Domain::equalsTo(1));

				//constraint = 0 matrix
				auto constraint = Expr::constTerm(Matrix::sparse(sizeOfSubgraph,sizeOfSubgraph));

				//build constraint
				int indexOfStableSet = 0;

				//generate stable set matrices and build expression for constraint
				for(std::vector<std::vector<int> >::iterator it1 = stableSets.begin() ; it1 != stableSets.end(); ++it1){
					std::vector<int> stableSet = *it1;
					//sparse matrix with sizeof(stableSet)^2 entries
					int numberMatrixEntries = stableSet.size()*stableSet.size();
					auto rows = new_array_ptr<int,1>(numberMatrixEntries);
					auto cols = new_array_ptr<int,1>(numberMatrixEntries);
					auto values = new_array_ptr<double,1>(numberMatrixEntries);

					int i = 0;
					for(std::vector<int>::iterator it2 = stableSet.begin() ; it2 != stableSet.end(); ++it2){
						//S_{j,j} = 1 \forall j \in stableSet
						(*rows)[i] = *it2 - 1;
						(*cols)[i] = *it2 - 1;
						(*values)[i] = 1;
						i++;
						//S_{i,j} = 1 \forall i,j \in stableSet, i \neq j
						for(std::vector<int>::iterator it3 = it2 + 1; it3 != stableSet.end(); ++it3){
							(*rows)[i] = *it2 - 1;
							(*cols)[i] = *it3 - 1;
							(*values)[i] = 1;
							i++;
							(*rows)[i] = *it3 - 1;
							(*cols)[i] = *it2 - 1;
							(*values)[i] = 1;
							i++;
						}
					}
					auto S = Matrix::sparse(sizeOfSubgraph,sizeOfSubgraph,rows,cols,values);

					//add to constraint
					auto subexpr = Expr::mul(lambda->index(indexOfStableSet),S); //avoid nesting
					constraint = Expr::add(constraint, subexpr);

					indexOfStableSet++;
				}

				//generate X_{subgraph}  (pick and reshape)
				//store required indices in arrays
				auto rows = new_array_ptr<int,1>(sizeOfSubgraph*sizeOfSubgraph);
				auto cols = new_array_ptr<int,1>(sizeOfSubgraph*sizeOfSubgraph);
				int iter = 0;
				for(std::vector<int >::iterator it1 = Subgraph.begin() ; it1 != Subgraph.end(); ++it1){
					for(std::vector<int >::iterator it2 = Subgraph.begin() ; it2 != Subgraph.end(); ++it2){
							(*rows)[iter] = *it1;
							(*cols)[iter] = *it2;
							iter++;
					}
				}

				auto subexpr = Var::reshape(X->pick(rows,cols),sizeOfSubgraph,sizeOfSubgraph)->asExpr();
				subexpr = Expr::neg(subexpr);
				constraint = Expr::add(constraint, subexpr);

				//add constraint
				M->constraint(constraint,Domain::equalsTo(Matrix::sparse(sizeOfSubgraph,sizeOfSubgraph)));


			}
		}
	}
	//================================================ added exact subgraph constraints ================================

#ifdef DEBUG
	if(params.printlevel > 3){
		std::cout << "added exact subgraph constraints" << std::endl;
	}
#endif

	//solve SDP
	//=======================================================================
	try{
		M->solve();
	}
	catch(...){
		counterBoundComputationFailed++;
		std::cerr << "Problem solving SDP" << std::endl;
		return -1;
	}
	//=============== solved SDP ===========================================


	//process mosek output
	//======================================================================
	//get solution status
	//------------------------
	//Model.getPrimalSolutionStatus()
	if(M->getPrimalSolutionStatus() == SolutionStatus::NearOptimal || M->getPrimalSolutionStatus() == SolutionStatus::Optimal){
		if(M->getPrimalSolutionStatus() == SolutionStatus::NearOptimal){
			counterMosekSolutionNearOptimal++;
		}

#ifdef DEBUG
		if(params.printlevel > 3){
			std::cout << "Problem (near) optimal" << std::endl;
		}
#endif
		// set variables
		//-------------------------------
		upperBound = M->primalObjValue();
		setVector(X, n);
		setBranchVar(branchVar, n);
		//setVectorAndBranchVar(X, n, branchVar);

		return 0;

	} else {
		counterBoundComputationFailed++;
		return -1;
	}
	//================== processed mosek output ==============================
}


std::vector<std::vector<int>> BoundComputationExactSubgraph::getStableSets(std::vector<int>& Subgraph, int* adjacency, int n){

	std::vector<std::vector<int>> stableSets;

	//-> check for all sets of vertices if stable
	//generate set
	int sizeOfSubgraph = Subgraph.size();
	std::vector<int> set;
	int mask = 1; //bitmask
	for(int i = 0; i<std::pow(2,sizeOfSubgraph); i++){
		//get new set
		set.clear();
		for(int j = 0; j<sizeOfSubgraph; j++){
			if( (i & (mask << j)) > 0  ){
				set.push_back(sizeOfSubgraph-j);
			}
		}

#ifdef DEBUG
		//print generated set
		if(params.printlevel>6){
			std::cout << "   Set: ";
			for(auto s: set){
				std::cout << Subgraph[s-1] << " ";
			}
		}
#endif
		//check if set stable
		bool stable = true;
		for(std::vector<int>::iterator it1 = set.begin() ; it1 != set.end() && stable; ++it1){
			for(std::vector<int>::iterator it2 = it1 + 1 ; it2 != set.end(); ++it2){
				//if(adjacency[Subgraph[*it1-1]-1][Subgraph[*it2-1]-1] == 1){
				if(adjacency[(Subgraph[*it1-1]-1) * n + (Subgraph[*it2-1]-1)] == 1){
					stable = false;
					break;
				}
			}
		}
		if(stable){
			stableSets.push_back(std::vector<int>(set));
#ifdef DEBUG
			if(params.printlevel>6){
				std::cout << " is stable" << std::endl;
			}
#endif
		}
#ifdef DEBUG
		else{
			if(params.printlevel>6){
				std::cout << " is not stable" << std::endl;
			}
		}
#endif
	}

	return stableSets;

}

