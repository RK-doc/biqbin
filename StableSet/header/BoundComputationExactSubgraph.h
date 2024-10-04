/*
 * BoundComputationExactSubgraph.h
 *
 *  Created on: 06.09.2017
 *      Author: mesiebenhofe
 */

#ifndef BOUNDCOMPUTATIONEXACTSUBGRAPH_H_
#define BOUNDCOMPUTATIONEXACTSUBGRAPH_H_

#include "BoundComputationTheta.h"

extern Params params;
extern int counterMosekSolutionNearOptimal;

class BoundComputationExactSubgraph: public BoundComputationTheta {
public:
	BoundComputationExactSubgraph();
	virtual int getUpperBound(Graph* g, double& upperBound, int& branchVar, int constValueUpperBound); //returns -1 on exit failure
	//int getLowerBound(Graph* g, int& lowerBound, int* sol); //returns -1 on exit failure

	virtual ~BoundComputationExactSubgraph();

protected:
	virtual std::vector<std::vector<int>> getStableSets(std::vector<int>& Subgraph, int* adjacency, int n);
};

#endif /* BOUNDCOMPUTATIONEXACTSUBGRAPH_H_ */
