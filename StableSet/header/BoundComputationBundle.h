/*
 * BoundComputationBundle.h
 *
 *  Created on: 20.09.2017
 *      Author: mesiebenhofe
 */

#ifndef BOUNDCOMPUTATIONBUNDLE_H_
#define BOUNDCOMPUTATIONBUNDLE_H_


#include "BoundComputationExactSubgraph.h"
#include "bdl_dgr.h"

extern int g_lowerBound;

using namespace mosek::fusion;
using namespace monty;

class BoundComputationBundle: public BoundComputationExactSubgraph {
public:
	BoundComputationBundle();

	virtual int getUpperBound(Graph* graph, double& upperBound, int& branchVar, int constValueUpperBound); //returns -1 on exit failure

	virtual ~BoundComputationBundle();

protected:
	std::vector<Matrix::t> getAllStableSetMatrices(std::vector<std::vector<int> >& stableSets, int sizeOfSubgraph);
	void getbAndyNeeded(int& b, std::vector<int>& yNeeded, int* adjacency,  std::vector<int>& Subgraph, int n);
	std::vector<std::vector<int>> getMatrixD(int b, int sizeSubgraph, std::vector<int>& yNeeded, std::vector<Matrix::t>& stableSetMatrices);
	void setVector(std::vector<double> X, int n);
	void searchForGoodSubgraphs(std::vector<double>& X, std::vector<double>& B, std::vector<int>& bestFoundSubgraph);
};

#endif /* BOUNDCOMPUTATIONBUNDLE_H_ */
