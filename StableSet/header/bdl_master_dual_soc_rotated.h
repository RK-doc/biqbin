/*
 * bdl_master_dual_soc_rotated.h
 *
 *  Created on: 04.10.2017
 *      Author: mesiebenhofe
 */

#ifndef BDL_MASTER_DUAL_SOC_ROTATED_H_
#define BDL_MASTER_DUAL_SOC_ROTATED_H_

#include "fusion.h"
#include "Structures.h"

using namespace mosek::fusion;
using namespace monty;

extern int counterMosekSolutionNearOptimal;
extern int counterBoundComputationFailed;
extern Params params;

class bdl_master_dual_soc_rotated {
public:
	bdl_master_dual_soc_rotated(std::vector<int>& tI, std::vector<int>& bI,
			std::vector<std::vector<std::vector<int>>>& DI);
	virtual ~bdl_master_dual_soc_rotated();
	int get_new_trial_point(std::vector<double>& yCenter, std::vector<std::vector<double>>& G,
							std::vector<double>& e, double tau, std::vector<double>& yTrial, std::vector<double>& alpha,
							double& timeQP);

protected:
	Matrix::t constConstraintMatrix;			//const submatrix of A
	std::shared_ptr<monty::ndarray<double,1>> b;   //right side of constraint Ax = b

	std::vector<int> m_bI;
	std::vector<int> m_tI;
	std::vector<std::vector<std::vector<int>>> m_DI;
	int m_q;
	int m_sum_bI;
	int m_sum_tI;

};

#endif /* BDL_MASTER_DUAL_SOC_ROTATED_H_ */
