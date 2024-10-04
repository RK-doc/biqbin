/*
 * bdl_dgr.cpp
 *
 *  Created on: 06.10.2017
 *      Author: mesiebenhofe
 */

#include "bdl_dgr.h"
#include <time.h>	//* clock_t, clock, CLOCKS_PER_SEC */
#include <iomanip> //std::fixed, std::setprecision
#include <cmath> //std::sqrt
#include <limits> //std::numeric_limits<double>::max()
#include <algorithm> //std::max_element, std::max


struct outputData {
  int iteration;
  bool nullstep;
  double fOpt;
  double fTrial;
  double normGlam;
  double delta;
  double mu;
  int r;
  double timeOracle;
  double timeQP;

} ;

bdl_dgr::bdl_dgr() {

}

bdl_dgr::~bdl_dgr() {

}

int bdl_dgr::bundle_for_exact_subgraph_constraints(Oracle_SS oracle, bdl_master_dual_soc_rotated soc_rotated,
		std::vector<double>& yStart, double muStart, std::string MCorSS, double lowerBoundToStop,
		double& fOpt, std::vector<double>& yOpt, std::vector<double>& gOpt, std::vector<double>& XMatrix){
//TODO: option eine zahl zu übergeben (g_lowerbound - const)

#ifdef DEBUG
	std::cout << std::setprecision(6) << std::fixed;
	if(params.printlevel > 3){
		std::cout << "----------------------------------" << std::endl;
		std::cout << "start bundle iteration: " << std::endl;
	}

	if(params.printlevel > 6){
		std::cout << "yStart = (";
		for(auto i: yStart){
			std::cout << i << ", ";
		}
		std::cout << ")" << std::endl;
	}


#endif

	outputData output;

	bool ss;
	if(MCorSS.compare("MC") == 0){
		ss = false;
	}
	else{
		ss = true;
	}


	if(ss){
		//create oracle --> besser als übergabeparameter?
	} else{
		std::cerr << "MC not implemented" << std::endl;
		return -1;
	}


	//create bundle (vector of bundle elements)
	std::vector<bundleElement> bundle;


	//bundle iteration
	//start
	double mu = muStart;
	double mL = params.mL;
	double epsStop = params.epsStop;
	clock_t time;

	//variables for kiwiel:
	//set up variables for kiwiel:
	int ikmu = 0;
	double epskv = std::numeric_limits<double>::max(); //infinity
	//other:
	double oldmu;
	double maxTerm;
	double mukint;
	std::vector<double> dk;
	std::vector<double> pk;
	double alphaTildapk;
	double alphaxkykPlus1;


#ifdef DEBUG
	if(params.printlevel > 3){
		std::cout << "mu = " << mu << std::endl;
		std::cout << "mL = " << mL << std::endl;
		std::cout << "epsStop = " << epsStop << std::endl;
	}
#endif

	//choose penalty in objective

	std::vector<double> yTrial = yStart;
	//call oracle for yTrial
	double fmax;
	double h;
	std::vector<double> gTrial;
	std::vector<double> alpha;
	std::vector<double> X;
	//--> call oracle (rückgabe -1 --> return -1)
	int status = oracle.oracle(yTrial, fmax, h, gTrial, X);
	if(status == -1){
		return -1;
	}



#ifdef DEBUG
	if(params.printlevel > 2){
		std::cout << "fStart = " << fmax + h  << std::endl;
		std::cout << "muStart = " << muStart << std::endl;
		std::cout << "it \t \t fOpt \t \t fTrial \t norm(dg) \t delta \t \t mu \t r \t tQP \t tOracle" << std::endl;
	}
	if(params.printlevel > 4){
		std::cout << "fmax = " << fmax << std::endl;
		std::cout << "h = " << h << std::endl;

		std::cout << "g = (";
		for(auto i: gTrial){
			std::cout << i << ", ";
		}
		std::cout << ")" <<std::endl;

	}
#endif

	std::vector<double> yCenter = yTrial;
	std::vector<double> gCenter = gTrial;

	bundle.emplace_back(yCenter, h, gCenter, X);


	//iteration
	std::vector<double> e;
	bool nullStep;

	double fTrial = h + fmax;
	double fCenter = fTrial;
	//double fCenter = h + fmax;

#ifdef DEBUG
	if(params.printlevel > 3){
		std::cout << "fCenter = fTrial = " << fCenter << std::endl;
	}
#endif

	for(int l = 0; l < params.Lmax; l++){
		output.iteration = l+1;
		//--> call get new trial point (rückgabe -1 --> return -1)
		//===========================================================
		//--> calculate e
		e.clear(); //TODO: set size of e to bundle.size()?
		//e_i = hCenter - h_i - <g_i,yCenter-y_i>
		double hCenter = bundle.back().h; //entry of center is the last entry in bundle
		for(uint i = 0; i < bundle.size(); i++){
			double tempValue = hCenter;
			tempValue -= bundle[i].h;       // h_i
			for(uint j = 0; j < bundle[i].g.size(); j++){
				tempValue -= bundle[i].g[j] * (yCenter[j] - bundle[i].y[j]);
			}
			e.push_back(tempValue);
		}

#ifdef DEBUG
		if(params.printlevel > 3){
			std::cout << "------------------------" << std::endl;
			std::cout << "Bundle iteration " << l+1 << std::endl;
		}
		if(params.printlevel > 6){
			std::cout << "e = (";
			for(auto i: e){
				std::cout <<  i << ", ";
			}
			std::cout << ")" << std::endl;
		}
#endif

		//create G  //TODO: nicht jedes mal neu erstellen?
		std::vector<std::vector<double>> G;
		for(uint i = 0; i < bundle.size(); i++){
			G.push_back(bundle[i].g);
		}

		//--> get new trial point
		//==============================================================
		output.mu = mu;
		output.r = bundle.size();
		time = clock();
		status = soc_rotated.get_new_trial_point(yCenter, G, e, 1/mu, yTrial, alpha, output.timeQP);
		time = clock() - time;
		if(status == -1){
			return -1;
		}

#ifdef DEBUG
		if(params.printlevel > 3){
			std::cout << "mu = " << mu << std::endl;
			std::cout << "r = " << bundle.size() << std::endl;
			std::cout << "time to solve SOCR = " << ((float)time)/CLOCKS_PER_SEC << std::endl;
		}

#endif
		//--> call oracle
		//==============================================================
		time = clock();
		status = oracle.oracle(yTrial, fmax, h, gTrial, X); //--> -1? return -1
		time = clock() - time;
		if(status == -1){
			return -1;
		}
		output.timeOracle = ((double)time)/CLOCKS_PER_SEC;




		//calculate fhat
		double fHat = fmax;
		//fhat += max{h_i + <g_i,y-y_i>: (y_i,h_i,g_i) \in Bundle}
		//fhat = optimale Lösung aus soc
		std::vector<double> temp;
		double tempValue;
		for(uint i = 0; i < bundle.size(); i++){
			tempValue = bundle[i].h;
			for(uint j = 0; j < bundle[i].g.size(); j++){
				tempValue += bundle[i].g[j] * (yTrial[j] - bundle[i].y[j]);
			}
			temp.push_back(tempValue);
		}

		fHat += *(std::max_element(temp.begin(),temp.end()));


		fTrial = h + fmax;
		double delta = fCenter - fHat; //f(yCenter_l) - (w+sum(vI))

		output.fTrial = fTrial;
		output.delta = delta;

		//calculate norm(glam) = abs(mu) * sqrt(\sum_{i} (yCenter_i - yTrial_i)^2)
		output.normGlam = 0;
		for(uint i = 0; i < yCenter.size(); i++){
			output.normGlam += (yCenter[i]-yTrial[i])*(yCenter[i]-yTrial[i]);
		}
		output.normGlam = std::sqrt(output.normGlam);
		output.normGlam *= std::abs(mu);




#ifdef DEBUG
		if(params.printlevel > 6){
			std::cout << "fmax = " << fmax << std::endl;
			std::cout << "fHat-fmax = " << fHat - fmax << std::endl;
			std::cout << "fhat = " << fHat <<std::endl;
			std::cout << "fCenter = " << fCenter << std::endl;
		}
		if(params.printlevel > 3){
			std::cout << "fTrial = " << fTrial << std::endl;
			std::cout << "fCenter = " << fCenter << std::endl;
			std::cout << "delta = " << delta << std::endl;
			std::cout << "time to solve oracle = " << ((float)time)/CLOCKS_PER_SEC << std::endl;
		}

#endif

		//nullstep oder serious step?
		//=======================================
		nullStep = fCenter - fTrial < mL * delta;


		//update mu:
		//-----------------
		//set variables:
		oldmu = mu;

		std::vector<double>DIMaxAttained;
		//TODO: use this computation of glam for normglam?
		for(int i = 0; i < yCenter.size(); i++){
			DIMaxAttained.push_back(mu*(yCenter[i]-yTrial[i]));
		}

		for(int j = 0; j < bundle.size(); j++){
			for(int i = 0; i < yCenter.size(); i++){
				DIMaxAttained[i] -= G[j][i]*alpha[j];
			}
		}

		//update serious step:
		//-------------------
		if(!nullStep){
			if((fTrial <= fCenter - params.mr*delta) && ikmu > 0){
				mu = std::max({mukint, 0.1*mu, params.muMin});
				//std::cout << "case 1: mu = " << mu << std::endl;
			}
			else if(ikmu > 3){
				mu = std::max(0.5*mu, params.muMin);
				//std::cout << "case 2: mu = " << mu << std::endl;
			}
			epskv = std::numeric_limits<double>::max(); //infinity
			if(oldmu == mu){
				ikmu = std::max(ikmu+1, 1);
			}
			else{
				ikmu = 1;
			}
			//std::cout << "ikmu = " << ikmu << std::endl;


			//old version: //TODO: param in paramfile to switch between the two versions
			//mu *= params.muSSFactor;


			//update center:
			//--------------
			yCenter = yTrial;
			gCenter = gTrial;
			fCenter = fTrial;
		}

		//update null step:
		//-------------------------
		else{
			//variables:
			//maxTerm = (fCenter - hCenter) - (DIMaxAttained'*yCenter)
			maxTerm = fCenter - hCenter;
			for(int i = 0; i < yCenter.size(); i++){
				maxTerm -= DIMaxAttained[i]*yCenter[i];
			}

			mukint = 2*mu*(1-(fTrial - fCenter + maxTerm)/(-std::abs(maxTerm - delta)));

			dk.clear();
			pk.clear();
			alphaTildapk = fCenter - fHat;
			alphaxkykPlus1 = hCenter - h;
			for(int i = 0; i < yCenter.size(); i++){
				dk.push_back(yTrial[i]-yCenter[i]);
				pk.push_back(-mu*dk[i]);
				//alphaTildapk = fCenter - fHat + pk'*dk
				alphaTildapk += pk[i]*dk[i];
				//alphaxkykPlus1 = hCenter - hTrial + gTrial'*dk
				alphaxkykPlus1 += gTrial[i]*dk[i];
			}
			//std::cout << "-DIMaxAtt'yCenter = " << maxTerm - (fCenter - hCenter) << std::endl;
			//std::cout << "fTrial - fCenter + maxTerm = " <<  fTrial - fCenter + maxTerm << std::endl;
			//std::cout << "maxterm - delta = " << maxTerm << "-" << delta << "= "<< (maxTerm - delta) <<std::endl;
			//std::cout << "-abs(maxterm-delta)= " << -std::abs(maxTerm-delta) << std::endl;
			//std::cout << "alphaTildapk = " << alphaTildapk << ", alphaxkykPlus1 = " << alphaxkykPlus1 << std::endl;

			double normpk = 0;
			for(int i = 0; i < pk.size(); i++){
				normpk += std::pow(pk[i],2);
			}
			normpk = std::sqrt(normpk);
			epskv = std::min(epskv, normpk + alphaTildapk);
			//std::cout << "epskv = " << epskv << std::endl;

			//set mu:
			//--------
			if(alphaxkykPlus1 > epskv && ikmu < -3){
				mu = std::min(mukint, 10*mu);
				//std::cout << "case 1: mu = " << mu << std::endl;
			}
			if(oldmu == mu){
				ikmu = std::min(ikmu-1, -1);
			}
			else{
				ikmu = -1;
			}
			//std::cout << "ikmu = " << ikmu << std::endl;



			//old version:
			//------------------
			//mu *= params.muNSFactor;
		}


		output.nullstep = nullStep;
		output.fOpt = fCenter;


#ifdef DEBUG
		if(params.printlevel > 3){
			std::cout << "nullstep = " << nullStep << std::endl;
		}
		if(params.printlevel > 2){
			std::cout << output.iteration << "\t";
			if(output.nullstep){
				std::cout << "N";
			} else {
				std::cout << "S";
			}
			std::cout << "\t" << output.fOpt << "\t" << output.fTrial << "\t" << output.normGlam << "\t" << output.delta <<
					"\t" << output.mu << "\t" << output.r << "\t"  << output.timeQP << "\t" << output.timeOracle << std::endl;
		}
#endif

		//stop bundle iteration?
		//===========================
		if(delta <= epsStop || fCenter - lowerBoundToStop + OPTEPS < 1){ //TODO: passt das so?
			//return fOpt, yOpt, gOpt
			yOpt = yCenter;
			fOpt = fCenter;
			gOpt = gCenter;
			//calculate X* = \sum_{i=0}^r \alpha_i X_i
			for(uint i = 0; i < bundle[0].X.size(); i++){
				XMatrix.push_back(0);
			}
			int j;
			for(uint i = 0; i < bundle.size(); i++){
				j = 0;
				for(std::vector<double>::iterator it = bundle[i].X.begin(); it != bundle[i].X.end(); it++){
					XMatrix[j] += (*it) * alpha[i];
					j++;
				}
			}

#ifdef DEBUG
			if(params.printlevel > 2){
				if(delta <= epsStop){
					std::cout << "delta = " << delta << " <= " << epsStop << std::endl;
				}
				else{
					std::cout << fCenter << "+ eps < 1 + " << lowerBoundToStop << std::endl;
				}
			}
			if(params.printlevel > 4){
				std::cout << "yOpt = (";
				for(auto i: yCenter){
					std::cout << i << ", ";
				}
				std::cout << ")" << std::endl;

				std::cout << "fOpt = " << fOpt << std::endl;

				std::cout << "gOpt = (";
				for(auto i: gOpt){
					std::cout << i << ", ";
				}
				std::cout << ")" << std::endl;
			}
#endif
			return 0;
		}

		//modify bundle:
		//====================================
		//in letzter Iteration nicht
		//delete elements from bundle:
		//delete all elements j from bundle (except last ~ center) where \alpha_j < 0.01
		if(l < params.Lmax - 1){ //do not delete from bundle in last iteration
			std::vector<bundleElement>::iterator bundleIterator = bundle.begin();
			for(std::vector<double>::iterator it = alpha.begin(); it != alpha.end() - 1; it++){
				if(*it < 0.01){
					bundleIterator = bundle.erase(bundleIterator);
				} else {
					bundleIterator++;
				}
			}

			if(nullStep){
				//current center needs to be last element in list (won't be deleted)
				bundleElement centerElement = bundle.back();
				bundle.pop_back();
				bundle.emplace_back(yTrial, h, gTrial, X);
				bundle.push_back(centerElement);
			}
			else{
				//new bundle elemente is center -> last element in list
				bundle.emplace_back(yTrial, h, gTrial, X);
			}
		}

	}

	//end of LMax iterations:
	//===================================
	yOpt = yCenter;
	fOpt = fCenter;
	gOpt = gCenter;

	//calculate X* = \sum_{i=0}^r \alpha_i X_i
	for(uint i = 0; i < bundle[0].X.size(); i++){
		XMatrix.push_back(0);
	}
	int j;
	for(uint i = 0; i < bundle.size(); i++){
		j = 0;
		for(std::vector<double>::iterator it = bundle[i].X.begin(); it != bundle[i].X.end(); it++){
			XMatrix[j] += *it * alpha[i];
			j++;
		}
	}

#ifdef DEBUG
	if(params.printlevel > 4){
		std::cout << "yOpt = (";
		for(auto i: yCenter){
			std::cout << i << ", ";
		}
		std::cout << ")" << std::endl;

		std::cout << "fOpt = " << fOpt << std::endl;

		std::cout << "gOpt = (";
		for(auto i: gOpt){
			std::cout << i << ", ";
		}
		std::cout << ")" << std::endl;
	}
#endif
	return 0;

}
