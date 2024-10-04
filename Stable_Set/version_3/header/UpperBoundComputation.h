/*
 * UpperBoundComputation.h
 *
 *  Created on: 10.08.2017
 *      Author: mesiebenhofe
 */

#ifndef UPPERBOUNDCOMPUTATION_H_
#define UPPERBOUNDCOMPUTATION_H_

#include "Graph.h"

class UpperBoundComputation {
public:
	virtual int getUpperBound(Graph* g, double& upperBound, int& branchVar, int constValueUpperBound)=0; //returns -1 on exit failure
	UpperBoundComputation();
	virtual ~UpperBoundComputation();
};

#endif /* UPPERBOUNDCOMPUTATION_H_ */
