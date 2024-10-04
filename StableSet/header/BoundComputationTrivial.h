/*
 * BoundComputationTrivial.h
 *
 *  Created on: 10.08.2017
 *      Author: mesiebenhofe
 */

#ifndef BOUNDCOMPUTATIONTRIVIAL_H_
#define BOUNDCOMPUTATIONTRIVIAL_H_

#include "LowerBoundComputation.h"
#include "UpperBoundComputation.h"

extern Params params;


class BoundComputationTrivial: public LowerBoundComputation,
		public UpperBoundComputation {
public:

	int getUpperBound(Graph* g, double& upperBound, int& branchVar, int constValueUpperBound); //returns -1 on exit failure
	int getLowerBound(Graph* g, int& lowerBound, int* sol); //returns -1 on exit failure

	BoundComputationTrivial();
	virtual ~BoundComputationTrivial();

private:
	bool isStableSet(Graph* g, std::vector<int>* Set, int vertex);

};

#endif /* BOUNDCOMPUTATIONTRIVIAL_H_ */
