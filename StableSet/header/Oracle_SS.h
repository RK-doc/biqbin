/*
 * Oracle_SS.h
 *
 *  Created on: 26.09.2017
 *      Author: mesiebenhofe
 */

#ifndef ORACLE_SS_H_
#define ORACLE_SS_H_


#include "Graph.h"

#include <vector>
#include "fusion.h"

using namespace mosek::fusion;
using namespace monty;

extern int counterMosekSolutionNearOptimal;
extern int counterBoundComputationFailed;

class Oracle_SS {
public:
	Oracle_SS(Graph* graph, std::vector<std::vector<int>> exactSubgraphs, std::vector<int> tI,
			 std::vector<int> bI, std::vector<std::vector<int> > yINeeded, std::vector<std::vector<std::vector<int>>> DI);
	virtual ~Oracle_SS();
	int oracle(std::vector<double> y, double& fmax, double& h, std::vector<double>& g, std::vector<double>& XMatrix); //returns fmax, h and g at y


protected:
	int n;
	std::vector<Edge> edges;
	std::vector<Matrix::t> staticConstraints; ///< input for theta function which remains the same in each iteration


	std::vector<std::vector<int>> m_exactSubgraphs;
	std::vector<int> m_tI;
	std::vector<int> m_bI;
	std::vector<std::vector<int> > m_yINeeded;
	std::vector<std::vector<std::vector<int>>> m_DI;

	std::shared_ptr<monty::ndarray<int,1>> rowsSubX;
	std::shared_ptr<monty::ndarray<int,1>> colsSubX;

};

#endif /* ORACLE_SS_H_ */
