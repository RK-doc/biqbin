#ifndef HEURISTIC_H_
#define HEURISTIC_H_

//using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <vector>

#include "Structures.h"

extern Params params;

typedef struct {
  double *s;
  double *y;
  double rho;
  double a;
} l_bfgs;

typedef struct {
  int clique;
  int rank;
  int numresets;
  int maxiter;
  int maxtime;
  int M;
  int savesoln;
  int solvemethod;
  int randstart;
} paramsH;



class Heuristic {
public:
	Heuristic();
	virtual ~Heuristic();
	 std::vector<int> getBestStableSetFound(int n, std::vector<Edge>& edges, int par_rank, int par_numReset, int maxSeconds);

protected:


	 int     normalize(double*);
	 double  function(double*, double*, double*, double*);
	 int     gradient(double*, double*, double*, double*);
	 int     projgradient(double*);
	 double  dot(double*, double*, int);
	 int     addvec(double*, double*, double, double*, int);
	 double  infnorm(double*, int);
	 int     direction(void);
	 int     updateLBFGS(double);
	 int     stableset(double*, double*);
	 double  f_inexact_lsearch(double*, double*, double*, double*, double,
	         double*, double, double*, double*, double*, double*);
	 double  f_bracketing(double*, double*, double*, double*, double,
	         double, double, double, double, double, double, double*,
	         double, double*, double*, double*, double*);
	 double  f_sectioning(double*, double*, double*, double*, double*,
	         double*, double*, double*, double, double, double, double, double,
	         double, double, double, double, double, double, double, double*, int);
	 int     f_interpolation_quadratic(double, double, double, double, double,
	         double, double, double*);
	 double  min(double, double);
	 double  max(double, double);
	 clock_t GetTime(void);
	 double  TimeInSec(clock_t, clock_t);
	 double  reset(double*, double*, double*, double*);
	 double  resetfull(double*, double*, double*, double*);
	 int     checkmaximal(void);


	//variables
	 int    n, m, oldest, iter, outeriter, tempint, maxcliquesize, rank;
	 int    *edgesi, *edgesj, stablesetsize, *beststableset, beststablesize;
	 int    *set1, *set2, *set3, *set4;
	 int    size1, size2, size3, size4;
	 double tempval, dirgradprod;
	 double *lambda, sigma, *grad, *projgrad, *dir;
	 double  etabar, omegabar, tau, gammabar, alpha_omega;
	 double  beta_omega, alpha_eta, beta_eta, alpha_star, beta_star;
	 double  eta_star, omega_star, alpha_k, eta_k, omega_k;
	//FILE   *solnfile;
	 paramsH *par;
	 l_bfgs *vecs;
};



#endif /* HEURISTIC_H_ */
