/*
 * LowerBoundComputation.h
 *
 *  Created on: 10.08.2017
 *      Author: mesiebenhofe
 */

#ifndef LOWERBOUNDCOMPUTATION_H_
#define LOWERBOUNDCOMPUTATION_H_

#include "Graph.h"

class LowerBoundComputation {
public:
	virtual int getLowerBound(Graph* g, int& lowerBound, int* sol)=0; //returns -1 on exit failure
	LowerBoundComputation();
	virtual ~LowerBoundComputation();
};

#endif /* LOWERBOUNDCOMPUTATION_H_ */
