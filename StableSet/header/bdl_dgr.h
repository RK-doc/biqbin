/*
 * bdl_dgr.h
 *
 *  Created on: 06.10.2017
 *      Author: mesiebenhofe
 */

#ifndef BDL_DGR_H_
#define BDL_DGR_H_

#include <string>
#include "Graph.h"
#include "Oracle_SS.h"
#include "bdl_master_dual_soc_rotated.h"

#include <algorithm>    // std::min_element, std::max_element


extern Params params;
//extern int g_lowerBound;

struct bundleElement{
	std::vector<double> y;
	double h;
	std::vector<double> g;
	std::vector<double> X;
	bundleElement(std::vector<double> yNew, double hNew, std::vector<double> gNew, std::vector<double> XNew){
		y = yNew;
		h = hNew;
		g = gNew;
		X = XNew;
	}
};

//struct for inputDetails?


class bdl_dgr {
public:
	bdl_dgr();
	virtual ~bdl_dgr();

	static int bundle_for_exact_subgraph_constraints(Oracle_SS oracle, bdl_master_dual_soc_rotated soc_rotated,
			std::vector<double>& yStart, double muStart, std::string MCorSS, double lowerBoundToStop,
			double& fOpt, std::vector<double>& yOpt, std::vector<double>& gOpt,
			std::vector<double>& XMatrix);
};

#endif /* BDL_DGR_H_ */
