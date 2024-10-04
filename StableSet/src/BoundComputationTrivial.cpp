/*
 * BoundComputationTrivial.cpp
 *
 *  Created on: 10.08.2017
 *      Author: mesiebenhofe
 */

#include "BoundComputationTrivial.h"
#include <stdlib.h> //rand()
#include <vector>

BoundComputationTrivial::BoundComputationTrivial() {
	std::srand(time(NULL)); //set seed
}

BoundComputationTrivial::~BoundComputationTrivial() {

}

int BoundComputationTrivial::getUpperBound(Graph* g, double& upperBound, int& branchVar, int constValueUpperBound){
	//returns -1 on exit failure
	//upperBound = #vertices
	int n = g->n();
	upperBound = n;
	if(n==0){
		branchVar = 0;
		return 0;
	}
	branchVar = std::rand() % n + 1;  //random branchvar (between 1 and n)

	return 0;
}

int BoundComputationTrivial::getLowerBound(Graph* g, int& lowerBound, int* sol){
	//returns -1 on exit failure
	std::vector<int> stableSet;


	for (std::vector<Vertex>::iterator itv = g->vertices().begin() ; itv != g->vertices().end(); ++itv){
		if(isStableSet(g, &stableSet, itv->lNumber)){
			stableSet.push_back(itv->lNumber);
		}
	}


	//check if this stable set is bigger than the current lower bound
	//if true --> set lower bound to size of the new stable set and update the
	//solution vector according to the stable set we found
	if(lowerBound < (int) stableSet.size()){
		lowerBound = stableSet.size();
#ifdef DEBUG
		if(params.printlevel>1){
			std::cout << "found better bound \n";
		}
#endif
		//update local solution vector
		int n = g->n();
		for(int i = 0; i<n; i++){
			sol[i] = 0;
		}
		for (std::vector<int>::iterator itv = stableSet.begin() ; itv != stableSet.end(); ++itv){
			sol[*itv - 1] = 1;
		}

	}

	return 0;
}


//returns true if still stable with vertex in set
bool BoundComputationTrivial::isStableSet(Graph* g, std::vector<int>* Set, int vertex){
	//for all vertices in set
	for (std::vector<int>::iterator itv = Set->begin() ; itv != Set->end(); ++itv){
		//check if there is an edge with this vertex and the new vertex
		for (std::vector<Edge>::iterator ite = g->edges().begin() ; ite != g->edges().end(); ++ite){
			if((ite->vertex1 == vertex && ite->vertex2 == *itv) || (ite->vertex1 == *itv && ite->vertex2 == vertex)){
				return false;
			}
		}
	}
	return true;
}
