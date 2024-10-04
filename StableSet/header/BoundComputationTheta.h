/*
 * BoundComputationTheta.h
 *
 *  Created on: 24.08.2017
 *      Author: mesiebenhofe
 */

#ifndef BOUNDCOMPUTATIONTHETA_H_
#define BOUNDCOMPUTATIONTHETA_H_

#include "UpperBoundComputation.h"
#include "LowerBoundComputation.h"
#include <vector>
#include "fusion.h"

extern Params params;
extern int counterMosekSolutionNearOptimal;
extern int counterBoundComputationFailed;

class BoundComputationTheta: public UpperBoundComputation,
		public LowerBoundComputation {
public:
	virtual int getUpperBound(Graph* g, double& upperBound, int& branchVar, int constValueUpperBound); //returns -1 on exit failure
	virtual int getLowerBound(Graph* g, int& lowerBound, int* sol); //returns -1 on exit failure
	BoundComputationTheta();
	virtual ~BoundComputationTheta();

protected:
	std::vector<float> m_vector;
	bool m_vectorSet;
	bool m_rootNode;

	//virtual void setVectorAndBranchVar(mosek::fusion::Variable::t X, int n, int& branchVar);
	virtual void setVector(mosek::fusion::Variable::t X, int n);
	virtual void setBranchVar(int& branchVar, int n);

	virtual void getStableSetBasedOnVector(Graph* g, std::vector<int>& stableSet);


};

#endif /* BOUNDCOMPUTATIONTHETA_H_ */
