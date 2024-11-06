/*
 * BoundComputationTheta.cpp
 *
 *  Created on: 24.08.2017
 *      Author: mesiebenhofe
 */

#include "BoundComputationTheta.h"

#include <stdlib.h> //rand(), abs()
#include <algorithm> //std::nth_elemnt, std::find

#include <time.h>


using namespace mosek::fusion;
using namespace monty;

BoundComputationTheta::BoundComputationTheta() {
	///first computation: node in branch and bound tree is the root
	///--> compute lower bound (param in paramfile) times
	m_rootNode = true;
	m_vectorSet = false;
	std::srand(time(NULL)); //set seed
}

BoundComputationTheta::~BoundComputationTheta() {

}

int BoundComputationTheta::getUpperBound(Graph* g, double& upperBound, int& branchVar, int constValueUpperBound){
#ifdef DEBUG
	clock_t m_time;
	m_time = clock();
#endif

	m_vectorSet = false;
	int n = g->n();

	///case: empty graph
	if(n==0){
		branchVar = 0;
		return 0;
	}

	///create Mosek model
	Model::t M = new Model("sdp"); auto _M = finally([&]() { M->dispose(); });

	///var X ... (n+1)x(n+1) matrix, psd
	Variable::t X = M->variable(Domain::inPSDCone(n+1));

	///generate matrix C for objective function
	///C = diag(0*I_{1}, I_{n})  ... (n+1)x(n+1)
	///C ... sparse matrix with n entries
	auto rows = new_array_ptr<int,1>(n);
	auto cols = new_array_ptr<int,1>(n);
	auto values = new_array_ptr<double,1>(n);
	for(int i = 0; i<n; i++){
		(*rows)[i] = i+1;
		(*cols)[i] = i+1;
		(*values)[i] = 1;
	}
	auto C = Matrix::sparse(n+1,n+1,rows,cols,values);

	///objective function:
	M->objective(ObjectiveSense::Maximize, Expr::dot(C, X));

	///add constraints:
	///constraint 1: X_{1,1} = 1
	M->constraint(X->index(0,0)->asExpr(), Domain::equalsTo(1.0));

	///constraints type 1: X_{i,i} = x_{i}
	for(int i = 1; i<n+1; i++){
		///X_{0,i} = X_{i,0} = X_{i,i} <==>
		///-0.5 X_{0,i} - 0.5 X_{i,0} + X_{i,i} = 0 <==>
		/// <A,X> = 0,  a_{0,i} = -0.5, a_{i,0} = -0.5, a_{i,i} = 1
		rows = new_array_ptr<int,1>({0,i,i});
		cols = new_array_ptr<int,1>({i,0,i});
		values = new_array_ptr<double,1>({-0.5,-0.5,1});
		auto A = Matrix::sparse(n+1,n+1,rows,cols,values);
		M->constraint(Expr::dot(A,X),Domain::equalsTo(0.0));
	}

	///constraints type 2: X_{i,j} = 0  \forall (i,j) \in E <==>
	///0.5 X_{i,j} + 0.5 X_{j,i} = 0
	std::vector<Edge> edges = g->edges();
	for (std::vector<Edge>::iterator it = edges.begin() ; it != edges.end(); it++ ){
		rows = new_array_ptr<int,1>({it->vertex1,it->vertex2});
		cols = new_array_ptr<int,1>({it->vertex2, it->vertex1});
		values = new_array_ptr<double,1>({0.5,0.5});
		auto A = Matrix::sparse(n+1,n+1,rows,cols,values);
		M->constraint(Expr::dot(A,X),Domain::equalsTo(0.0));
	}

#ifdef DEBUG
	m_time = clock() - m_time;
	if(params.printlevel>4){
		std::cout << "time to create input for sdp: " << ((float)m_time)/CLOCKS_PER_SEC << std::endl;
	}
#endif


#ifdef DEBUG
	m_time = clock();
#endif
	//std::cout << "before solve" << std::endl;
	///solve SDP
	
	M->setSolverParam("intpntMultiThread", "off");
	//M->setSolverParam("numThreads", 2);
	M->solve();
	//std::cout << "after solve" << std::endl;
#ifdef DEBUG
	m_time = clock() - m_time;
	std::cout << "MOSEK SDP solver: " << m_time << std::endl;

	if(params.printlevel>4){
		std::cout << "time to solve sdp: " << ((float)m_time)/CLOCKS_PER_SEC << std::endl;
	}

	m_time = clock();
#endif

	///check solution status
	///if successful: set upper bound, vector and branching variable
	///returns 0 if (near) optimal, else -1

	if(M->getPrimalSolutionStatus() == SolutionStatus::NearOptimal || M->getPrimalSolutionStatus() == SolutionStatus::Optimal){
		if(M->getPrimalSolutionStatus() == SolutionStatus::NearOptimal){
			counterMosekSolutionNearOptimal++;
		}

#ifdef DEBUG
		if(params.printlevel>4){
			if(M->getPrimalSolutionStatus() == SolutionStatus::Optimal){
				std::cout << "Problem optimal" << std::endl;
			}
			else if(M->getPrimalSolutionStatus() == SolutionStatus::NearOptimal){
				std::cout << "Problem near optimal" << std::endl;
			}
		}
#endif

		///set upper bound
		upperBound = M->primalObjValue();
		setVector(X, n);
		setBranchVar(branchVar, n);
		//setVectorAndBranchVar(X, n, branchVar);
		return 0;
	}
	else{
		///case: not (near) optimal
		counterBoundComputationFailed++;
		return -1;
	}
}

int BoundComputationTheta::getLowerBound(Graph* g, int& lowerBound, int* sol){
	std::vector<int> stableSet;

	int n = g->n();

	if(n==0){
		return 0;
	}



	if(!m_vectorSet){
		///if vector not set -> set random
		m_vector.clear();
		for(int i = 0; i<n; i++){
			m_vector.push_back(std::rand()/RAND_MAX); // random number \in [0,1]
		}
	}

	getStableSetBasedOnVector(g, stableSet);

	///if root node --> perform heuristic again with random vector
	if(m_rootNode){
		std::vector<int> newStableSet;

		for(int i = 0; i<params.lowerBoundInRoot; i++){
			///set vector randomly
			m_vector.clear();
			for(int i = 0; i<n; i++){
				m_vector.push_back(std::rand()/RAND_MAX); // random number \in [0,1]
			}

			newStableSet.clear();
			getStableSetBasedOnVector(g,newStableSet);


			///better lower bound found? yes: change stableSet to new, bigger stable set
			if(newStableSet.size() > stableSet.size()){
				stableSet.swap(newStableSet);
			}
		}

		m_vectorSet = false;
		m_rootNode = false;
	} //end if(rootNode)

	///update lower bound
	if(lowerBound < (int) stableSet.size()){
		lowerBound = stableSet.size();
		///update solution vector
		for(int i = 0; i<n; i++){
			sol[i] = 0;
		}
		for(unsigned i = 0; i < stableSet.size(); i++){
			sol[stableSet[i] - 1] = 1;
		}
	}

	return 0;
}

///probiert: Model::t als objekt, referenz, zeiger übergeben - nicht funktioniert

void BoundComputationTheta::setVector(mosek::fusion::Variable::t X, int n){
	///set vector (if no errors occurred)
	m_vector.clear();
	for(int i = 0; i<n; i++){
		m_vector.push_back((*(X->index(0,i+1)->level()))[0]);
	}
	m_vectorSet = true;
}

void BoundComputationTheta::setBranchVar(int& branchVar, int n){
	//choosing branching variable
	//---------------------------

	//2) position of min in vector
	if(params.strategyBranchVar == 2){
		int index = 0;
		int i = 0;
		int min = m_vector[0];
		for(std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); it++){
			if(*it < min){
				min = *it;
				index = i;
			}
			i++;
		}
		branchVar = index + 1;
	}
	//3) position of max in vector
	else if(params.strategyBranchVar == 3){
		int index = 0;
		int i = 0;
		int max = m_vector[0];
		for(std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); it++){
			if(*it > max){
				max = *it;
				index = i;
			}
			i++;
		}
		branchVar = index + 1;
	}

	//4) closest to mean
	else if(params.strategyBranchVar == 4){
		float mean = 0;
		for(std::vector<float>::iterator it = m_vector.begin(); it!=m_vector.end(); it++){
			mean += *it;
		}
		mean = mean/n;

		int index = 0;
		int i = 0;
		int min = abs(m_vector[0] - mean);

		for (std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); it++){
			if(abs(*it - mean) < min){
				index = i;
			}
			i++;
		}
		branchVar = index + 1;
	}



	//5) closest to 0.5
	else if(params.strategyBranchVar == 5){
		///branch by fractional
		int index = 0;
		int i = 0;
		int min = abs(m_vector.at(0)-0.5);
		for (std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); it++){
			if(abs(*it-0.5) < min){
				min = abs(*it-0.5);
				index = i;
			}
			i++;
		}
		branchVar = index + 1;
	}


	//6) median
	else if(params.strategyBranchVar == 6){
		std::vector<float> vecCopy = m_vector;
		std::nth_element(vecCopy.begin(), vecCopy.begin()+vecCopy.size()/2, vecCopy.end()); //rearranges vector -> use copy
		float median = *(vecCopy.begin() + vecCopy.size()/2);
		int index = std::find(m_vector.begin(), m_vector.end(), median) - m_vector.begin(); //get position of median

		branchVar = index + 1;

	}

	//1) random
	else {
		///set random branching variable
		branchVar = std::rand() % n + 1;  ///< random branchvar (between 1 and n)
	}

}

/*
void BoundComputationTheta::setVectorAndBranchVar(Variable::t X, int n, int& branchVar){
	///set vector (if no errors occurred)
	m_vector.clear();
	for(int i = 0; i<n; i++){
		m_vector.push_back((*(X->index(0,i+1)->level()))[0]);
	}
	m_vectorSet = true;

	//choosing branching variable
	//---------------------------

	//2) position of min in vector
	if(params.strategyBranchVar == 2){
		int index = 0;
		int i = 0;
		int min = m_vector[0];
		for(std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); it++){
			if(*it < min){
				min = *it;
				index = i;
			}
			i++;
		}
		branchVar = index + 1;
	}
	//3) position of max in vector
	else if(params.strategyBranchVar == 3){
		int index = 0;
		int i = 0;
		int max = m_vector[0];
		for(std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); it++){
			if(*it > max){
				max = *it;
				index = i;
			}
			i++;
		}
		branchVar = index + 1;
	}

	//4) closest to mean
	else if(params.strategyBranchVar == 4){
		float mean = 0;
		for(std::vector<float>::iterator it = m_vector.begin(); it!=m_vector.end(); it++){
			mean += *it;
		}
		mean = mean/n;

		int index = 0;
		int i = 0;
		int min = abs(m_vector[0] - mean);

		for (std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); it++){
			if(abs(*it - mean) < min){
				index = i;
			}
			i++;
		}
		branchVar = index + 1;
	}



	//5) closest to 0.5
	else if(params.strategyBranchVar == 5){
		///branch by fractional
		int index = 0;
		int i = 0;
		int min = abs(m_vector.at(0)-0.5);
		for (std::vector<float>::iterator it = m_vector.begin() ; it != m_vector.end(); it++){
			if(abs(*it-0.5) < min){
				min = abs(*it-0.5);
				index = i;
			}
			i++;
		}
		branchVar = index + 1;
	}


	//6) median
	else if(params.strategyBranchVar == 6){
		std::vector<float> vecCopy = m_vector;
		std::nth_element(vecCopy.begin(), vecCopy.begin()+vecCopy.size()/2, vecCopy.end()); //rearranges vector -> use copy
		float median = *(vecCopy.begin() + vecCopy.size()/2);
		int index = std::find(m_vector.begin(), m_vector.end(), median) - m_vector.begin(); //get position of median

		branchVar = index + 1;

	}

	//1) random
	else {
		///set random branching variable
		branchVar = std::rand() % n + 1;  ///< random branchvar (between 1 and n)
	}

}
*/

void BoundComputationTheta::getStableSetBasedOnVector(Graph* g, std::vector<int>& stableSet){
	std::vector<Vertex> vertices = g->vertices();
	std::vector<Edge> edges = g->edges();

	stableSet.clear();
	std::vector<int>adjacentVertices;

	/**
	while(knotenmenge!=0){
		suche unter allen knoten in liste den mit größtem eintrag
		dieser Knoten kommt ins stable set
		lösche alle Knoten, die zu diesem adjazent sind (und deren Kanten)
	*/

	while(!vertices.empty()){
		///find vertex with greatest entry in vector
		int index = 0;
		for(int i = 1; i<(int)vertices.size(); i++){
			if(m_vector[vertices.at(i).lNumber - 1] > m_vector[vertices.at(index).lNumber-1]){
				index = i;
			}
		}
		///delete vertex from list
		int vertexNumber = vertices.at(index).lNumber;
		vertices.erase(vertices.begin()+index);

		///put vertex in stableSet
		stableSet.push_back(vertexNumber);

		///delete adjacent vertices and edges
		adjacentVertices.clear();
		///delete edges to adjacent vertices
		for (std::vector<Edge>::iterator it = edges.begin() ; it != edges.end(); ){
		    if(it->vertex1==vertexNumber){
		    	adjacentVertices.push_back(it->vertex2);
		    	it = edges.erase(it);
		    } else if(it->vertex2==vertexNumber){
		    	adjacentVertices.push_back(it->vertex1);
		    	it = edges.erase(it);
		    }
		    else {
		    	it++;
		    }
		}

		/// \forall adjacent vertices:
		/// 		delete edges
		///		delete vertex
		for (std::vector<int>::iterator itv = adjacentVertices.begin() ; itv != adjacentVertices.end(); itv++){
		    ///delete edges
			for(std::vector<Edge>::iterator ite = edges.begin(); ite != edges.end(); ){
		    	if(ite->vertex1 == *itv || ite->vertex2 == *itv){
		    		ite = edges.erase(ite);
		    	}
		    	else{
		    		ite++;
		    	}
		    }
			///delete vertex
		    for(std::vector<Vertex>::iterator itv2 = vertices.begin(); itv2 != vertices.end(); ){
		    	if(itv2->lNumber==*itv){
		    		itv2 = vertices.erase(itv2);
		    		break;
		    	} else {
		    		itv2++;
		    	}
		    }
		}
	}
}


