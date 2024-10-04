/*
function [permOpt, costOpt ] = qap_simul2_c(H,X);
        input: H .. kxk SYMMETRIC matrix
               X .. nxn SYMMETRIC matrix (with n > k)
        output: permOpt ... a Permutation of 1, ..., n
                costOpt ... the cost of <X(permOpt,permOpt), H> =
                                       sum_{i,j = 1,...,n} h_{i,j}x_{perm(i),perm(j)}
                   in such a way, that the cost is minimised (over all permutations) with simulated annealing
        optional output: (if OUTPUT_COST_OF_ACCEPTED_SUBMATRIX is defined)
                         c ... a Matrix with the cost of all accepted permutations
        call with optional output: >> [a,b,c] = qap_simul2_c(H,X); figure; plot(c)
*/
/*ATTENTION: H and X have to be symmetric! */


#include <math.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mex.h"

#define H_IN prhs[0]
#define X_IN prhs[1]
#define PERM_OUT plhs[0]
#define COST_OUT plhs[1]


/*------- Parameters to change -------------*/

#define MITER_KOEFFICIENT_OF_N 1  /*inner iterations*/
#define TRIALS 1                 /*number of trials (restarts)   20 */
#define FT 0.7                    /*factor <1 to reduce temperature    0.6 */
#define FITER 1.15                 /* factor to increase inner trials    1.1 */  

/*Uncomment if optional output is required*/
/*#define OUTPUT_COST_OF_ACCEPTED_SUBMATRIX*/

/*------ End of Parameters to change ------ */


#ifdef OUTPUT_COST_OF_ACCEPTED_SUBMATRIX
#define COST_OF_ACCEPTED_SUBMATRIX_OUT plhs[2]
#endif


void qap_simul(
   int k,
   double **HtwoD, /* array with data from ARRAY_IN */
   int n,
   double **XtwoD,
   double optCost[],     /* vector with data for VECTOR_OUT */
   double optPerm[]      /* vector with data from VECTOR_OUT */
#ifdef OUTPUT_COST_OF_ACCEPTED_SUBMATRIX
   , double costOfAcceptedSubmatrix[]
#endif   
        )
      
{
    
    
  /* Definition of the parameters of simulated annealing */     
  int miter = round(MITER_KOEFFICIENT_OF_N*n);            /*inner iterations*/
  double trials = TRIALS;        /*number of trials (restarts) */
  double ft = FT;          /*factor <1 to reduce temperature*/
  double fiter = FITER;       /* factor to increase inner trials*/
  
  
  /*Definition of needed Variables*/
  int nr_trial, nr_perm, nr_index, col, row, nr_it_t;
  int randomNumber, smallerIndices;
  
  /*Initialisaton of random number generator*/
  /*ATTENTION: srand() is not called because qap_simul2_c is called several
   * times within a second. time() changes only every second, hence srand 
   * is always initalised with the same value and hence produces always the
   * same random numbers and hence always the same output for the optimal 
   * permutation.
   * But MATLAB does not close MEX-Files automatically, hence if srand is
   * not initialized after one run, then it remains in the state of the
   * last call and the random numbers produced are all different */
  /*srand(time(NULL));*/
  
  
  /*Calculating initial temperature t = ( sum(sum(abs(H))) * sum( sum( abs(X))))  / n / (n-1);*/
  double totalSumX = 0;
  for (col=0; col < n; col++)
  {
      for (row=col+1; row < n; row++)
      {
          totalSumX += fabs(XtwoD[row][col]);
      }
  }
  totalSumX *= 2;
  for (col=0; col < n; col++)
  {
      totalSumX += fabs(XtwoD[col][col]);
  }
  /*mexPrintf("beginning totalSumX: %f\n", totalSumX);
  mexEvalString("drawnow;");*/
  double totalSumH = 0;
  for (col=0; col < k; col++)
  {
      for (row=col+1; row < k; row++)
      {
          totalSumH += fabs(HtwoD[row][col]);
      }
  }
  totalSumH *= 2;
  for (col=0; col < k; col++)
  {
      totalSumH += fabs(HtwoD[col][col]);
  }
  /*mexPrintf("beginning totalSumH: %f\n", totalSumH);
  mexEvalString("drawnow;");*/
  double t = totalSumX/((double)n)*totalSumH/((double)(n-1));
  /*mexPrintf("beginning t: %f\n", t);
  mexEvalString("drawnow;");*/
  
  double bestFoundCost = 10000000000;
  int bestFoundPerm[n];
  
 
  
  /*Iteration for one starting permutation*/
  for(nr_trial = 1; nr_trial <= trials; nr_trial++)
  {
#ifdef OUTPUT_COST_OF_ACCEPTED_SUBMATRIX
     int indexOfCostOfAcceptedSubmatrix = 500*(nr_trial - 1);
#endif
      
      
      /*Generate a random pemutation of size n*/
      int perm[n];
      for(nr_perm = 1; nr_perm <= n; nr_perm++)
      {
          perm[nr_perm-1] = nr_perm;
      }
      for(nr_perm = n; nr_perm > 0; nr_perm--)
      {
          randomNumber = (rand() % nr_perm) + 1;
          int temp = perm[nr_perm - 1];
          perm[nr_perm - 1] = perm[randomNumber - 1];
          perm[randomNumber - 1] = temp;
      }
      /*mexPrintf("Permutation: [");
      for(nr_perm = 0; nr_perm < n; nr_perm++)
      {
          mexPrintf("%d, ",(int) perm[nr_perm]);
      }
      mexPrintf("]\n");
      mexEvalString("drawnow;");
      */
      
      /*Initalisating current values*/
      double t1 = t;
      int m1 = miter;
      /*mexPrintf("t1 =%f, m1 = %d\n", t1, m1);
      mexEvalString("drawnow;");*/
      
      double curBestSol = 0;
      for (col=0; col < k; col++)
      {
          for (row=col+1; row < k; row++)
          {
              curBestSol += HtwoD[row][col]*XtwoD[perm[row]-1][perm[col]-1];
          }
      }
      curBestSol = curBestSol*2;
      for (col=0; col < k; col++)
      {
          curBestSol += HtwoD[col][col]*XtwoD[perm[col]-1][perm[col]-1];
      }
      /*mexPrintf("curBestSol=%f\n", curBestSol);
      mexEvalString("drawnow;");*/
      
      
      

      int fertig = 1;
      
      
      /*Iterations over temperature t1*/
      while(fertig > 0) /*while things change*/
      {
          fertig = 0;
          
          
          /* m1 Iterations at constant temperature*/
          for(nr_it_t = 0; nr_it_t < m1; nr_it_t++)
          {
              
              
              /*Generate two random variables i1 in (1,...,k) and i2 (1,...,n) with i1 <= i2*/
              int i1 = (rand() % k) + 1;
              int i2 = (rand() % n) + 1;
              if(i2 < i1)
              {
                  int temp = i1;
                  i1 = i2;
                  i2 = temp;
              }
              /* mexPrintf("i1=%d, i2=%d\n", i1, i2); */
              
              
              
              /*Compute the new objective function value*/
              double tempSol = 0;
              double changeOfSol = 0;
              for (col=0; col < k; col++)
              {
                  changeOfSol += HtwoD[i1-1][col]*(XtwoD[perm[i2-1]-1][perm[col]-1] - XtwoD[perm[i1-1]-1][perm[col]-1]);
              }
              if(i2 <= k)
              {
                  for (col=0; col < k; col++)
                  {
                      changeOfSol += HtwoD[i2-1][col]*(XtwoD[perm[i1-1]-1][perm[col]-1] - XtwoD[perm[i2-1]-1][perm[col]-1]);
                  }
              }
              changeOfSol *= 2;
              double tempValue1 = (XtwoD[perm[i1-1]-1][perm[i1-1]-1] - 2*XtwoD[perm[i1-1]-1][perm[i2-1]-1] + XtwoD[perm[i2-1]-1][perm[i2-1]-1]);
              double tempValue2 = HtwoD[i1-1][i1-1];
              if(i2 <= k)
              {
                  tempValue2 += -2*HtwoD[i2-1][i1-1] + HtwoD[i2-1][i2-1];
              }
              /*changeOfSol -= HtwoD[i1-1][i1-1]*(XtwoD[perm[i2-1]-1][perm[i1-1]-1] - XtwoD[perm[i1-1]-1][perm[i1-1]-1]);
              if(i2 <= k)
              {
                  changeOfSol -= (HtwoD[i2-1][i2-1]*(XtwoD[perm[i1-1]-1][perm[i2-1]-1] - XtwoD[perm[i2-1]-1][perm[i2-1]-1]) +
                                  2*HtwoD[i2-1][i1-1]*(XtwoD[perm[i1-1]-1][perm[i1-1]-1] - XtwoD[perm[i2-1]-1][perm[i1-1]-1]));
              }*/
              changeOfSol += tempValue1*tempValue2;
              tempSol = curBestSol + changeOfSol;
              /*mexPrintf("changeOfSol=%f\n", changeOfSol);
              mexPrintf("tempSol=%f\n", tempSol);*/
              
              
              /* Determine whether the swap is accepted */
              int accept;
              if(changeOfSol > 0) /*If solution is worse accept with certain probability*/
              {
                  double dt1 = changeOfSol/t1;
                  if(dt1>5)
                  {
                      accept = 0;
                  }
                  else
                  {
                      double prob = exp(-dt1);
                      if( (rand()/(RAND_MAX + 1.)) < prob)
                      {
                          accept = 1;
                      }
                      else
                      {
                          accept = 0;
                      }
                  }
              }
              else /* If solution is better accept */
              {
                  accept = 1;
              }
              
              
              /* Do the uptdate if the swap is accepted */
              if(accept == 1)
              {
                  /*mexPrintf("acctepted with changeOfSol=%f\n", changeOfSol);*/
                  if(fabs(changeOfSol) > 0.0001)
                  {
                      fertig = 1;
                  }
                  int temp = perm[i1-1];
                  perm[i1-1] = perm[i2-1];
                  perm[i2-1] = temp;
                  curBestSol = tempSol;
                  
                  
#ifdef OUTPUT_COST_OF_ACCEPTED_SUBMATRIX
                  costOfAcceptedSubmatrix[indexOfCostOfAcceptedSubmatrix] = curBestSol;
                  indexOfCostOfAcceptedSubmatrix++;
#endif
                  
                  
                  /* If the just found solution is the best one of all, store it */
                  if(curBestSol < bestFoundCost)
                  {
                      bestFoundCost = curBestSol;
                      for(nr_perm = 1; nr_perm <= n; nr_perm++)
                      {
                          bestFoundPerm[nr_perm-1] = perm[nr_perm-1];
                          /*mexPrintf("bestFoundPerm[nr_perm-1]=%d\n", bestFoundPerm[nr_perm-1]);*/
                      }
                      /*mexPrintf("With cost %f\n",bestFoundCost);
                      mexEvalString("drawnow;");*/
                  }
                  
              }

              
              
          } /* End: m1 Iterations at constant temperature*/
          
          t1 *= ft;
          m1 = round(m1*fiter);
          /*mexPrintf("t1 =%f, m1 = %d\n", t1, m1);
          mexEvalString("drawnow;");*/
          
      } /*End: Iterations over temperature t1*/
      
  } /*End: Iteration for one starting permutation*/
  
  
  optCost[0] = bestFoundCost;
  for(nr_perm = 1; nr_perm <= n; nr_perm++)
  {
      optPerm[nr_perm-1] = bestFoundPerm[nr_perm-1];
  }
       
  
}


void mexFunction( int nlhs, mxArray *plhs[],
		          int nrhs, const mxArray *prhs[])
{
    
    /* checking input/output arguments*/
    if (nrhs != 2)
    {
       mexErrMsgTxt("2 Input Arguments required:\n H ... k x k small symmetric matrix\n X ... n x n big symmetric matrix\n");
    }
#ifdef OUTPUT_COST_OF_ACCEPTED_SUBMATRIX
    if (nlhs != 3)
    {
       mexErrMsgTxt("3 Output Arguments required:\n perm ... best found permutation\n cost ... cost of this permutation, i.e. <H,X(perm,perm)>\n cost_history ... the costs of all accepted permuatations\n");
    }    
#else
    if (nlhs != 2)
    {
       mexErrMsgTxt("2 Output Arguments required:\n perm ... best found permutation\n cost ... cost of this permutation, i.e. <H,X(perm,perm)>\n");
    }    
#endif        


    
    double *in1, *out0, *out1, *mw, *position, *d;
	int n1;
    
   
    /* Setting the variables*/
    int k = mxGetM(H_IN); /*k = number of rows = number of cols of H */
    int n = mxGetM(X_IN); /*n = number of rows = number of cols of X */
    int row, col, i; /*loop indices */
    double **XtwoD;
    double **HtwoD;
    

    /* Allocate memory for X and H */
    XtwoD = (double **) malloc(n*sizeof(double *));
    for (i=0; i<n; i++)
    {
        XtwoD[i] = (double *) malloc(n*sizeof(double));
    }
    
    
    HtwoD = (double **) malloc(k*sizeof(double *));
    for (i=0; i<k; i++)
    {
        HtwoD[i] = (double *) malloc(k*sizeof(double));
    }
    
    /* Allocate memory for return arguments */ 
    PERM_OUT = mxCreateDoubleMatrix(n, 1, mxREAL);
    COST_OUT = mxCreateDoubleMatrix(1, 1, mxREAL);
    
    
#ifdef OUTPUT_COST_OF_ACCEPTED_SUBMATRIX
    COST_OF_ACCEPTED_SUBMATRIX_OUT = mxCreateDoubleMatrix(500, TRIALS, mxREAL);
#endif
    
    
    /* Convert X_IN and H_IN to a 2x2 C array */
    /*MATLAB stores a two-dimensional matrix in memory as a one-
    dimensional array.  If the matrix is size MxN, then the
    first M elements of the one-dimensional array correspond to
    the first column of the matrix, and the next M elements
    correspond to the second column, etc. The following loop
    converts from MATLAB format to C format: */
      
    for (col=0; col < n; col++)
    {
        for (row=0; row < n; row++)
        {
            XtwoD[row][col] = (mxGetPr(X_IN))[row+col*n];
        }
    }
    
    for (col=0; col < k; col++)
    {
        for (row=0; row < k; row++)
        {
            HtwoD[row][col] = (mxGetPr(H_IN))[row+col*k];
        }
    }
    

    qap_simul(k, HtwoD, n, XtwoD, mxGetPr(COST_OUT), mxGetPr(PERM_OUT)
#ifdef OUTPUT_COST_OF_ACCEPTED_SUBMATRIX
    , mxGetPr(COST_OF_ACCEPTED_SUBMATRIX_OUT)
#endif
    );
    
    
   /* deallocate the arrays */
   for (i=0; i<k; i++)
   {
       free(HtwoD[i]);
   }
   free(HtwoD);
   
   for (i=0; i<n; i++)
   {
       free(XtwoD[i]);
   }
   free(XtwoD);
   
}













