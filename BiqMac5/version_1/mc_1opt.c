/* 
 * improves given cut.
 *
 * input:  L, n, cutvector x
 * output: cost, xnew
 *
 * NOTE: x is of type double, though it consists only of 
 * elements {+1,-1}. 
 *
 * translation from mc_1opt.m (f. rendl)
 * nov 2004, a.wiegele
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"


int mc_1opt(L,n,px,pcost)
    double *L;
    int n;
    /* input/output */
    double **px;     
    /* output */
    double *pcost;
{   
  /* vars for call to blas- and lapack-routines */
  double alpha,beta;
  int inc=1;
  char up='u';
  
  /* and some other var's */
  int i,j;
  double *Lx, *diagL;
  double best, cost;
  int bestarg;
  double *pLx,*pL;
  
  /**********************************
   *                                *         
   * initialization                 *
   *                                *         
   **********************************/
  
  /* allocate space for Lx and diagL */
  Lx=(double *)malloc((n+1)*sizeof(double));
  diagL=(double *)malloc((n+1)*sizeof(double));
  if ((Lx == NULL)||(diagL == NULL))
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };

  alpha=1.0;beta=0.0;               
  dspmv_(&up,&n,&alpha,L,*px+1,&inc,&beta,Lx+1,&inc); /* Lx=L*x */
  
  cost=ddot_(&n,*px+1,&inc,Lx+1,&inc); /* cost = x'*Lx */
  diag(L,&diagL,n); /* diagL=diag(L) */
  
  deltamax(diagL,*px,Lx,n,&best,&bestarg);
  
  while (best > EPSILON)
    {
      if ((*px)[bestarg] > 0)
	{
	  j=bestarg*(bestarg-1)/2;
	  /* scan column until pivot */
	  for (i=1,pLx=Lx,pL=L+(j-1);i<=bestarg;i++)
	    {
	      pLx++; pL++;
	      (*pLx)-=2*(*pL);
	    };
	  /* scan row from pivot on */
	  for (i=bestarg;i<n;i++)
	    {
	      pLx++; pL+=i;
	      (*pLx)-=2*(*pL);
	    };
	  
	  /*
	    for (i=1;i<=n;i++)
	    if (i <= bestarg)
	    Lx[i]-=2*L[ijtokp(i,bestarg,n)];
	    else
	    Lx[i]-=2*L[ijtokp(bestarg,i,n)];
	  */
	  (*px)[bestarg]=-1.0;
	}
      else
	{
	  j=bestarg*(bestarg-1)/2;
	  /* scan column until pivot */
	  for (i=1,pLx=Lx,pL=L+(j-1);i<=bestarg;i++)
	    {
	      pLx++; pL++;
	      (*pLx)+=2*(*pL);
	    };
	  /* scan row from pivot on */
	  for (i=bestarg;i<n;i++)
	    {
	      pLx++; pL+=i;
	      (*pLx)+=2*(*pL);
	    };
	  
	  /*
	  for (i=1;i<=n;i++)
	    if (i<= bestarg)
	      Lx[i]+=2*L[ijtokp(i,bestarg,n)];
	    else
	      Lx[i]+=2*L[ijtokp(bestarg,i,n)];
	  */
	  (*px)[bestarg]=1.0;
	};
      cost+=4*best;
      deltamax(diagL,*px,Lx,n,&best,&bestarg);
    }
  
  *pcost=cost;
  
  /* free and return */
  free(diagL);free(Lx);
  return(0);
}

/*
 *
 */
void deltamax(diagL,x,Lx,n,best,bestarg)
     double *diagL,*x,*Lx;
     int n;
     double *best;
     int *bestarg;
{
  int i;
  double *px,*pdiagL,*pLx;
  
  *bestarg=1;
  *best=( x[1] > 0 ? diagL[1]-Lx[1] : diagL[1]+Lx[1] );
  
  for (i=2,px=x+2,pdiagL=diagL+2,pLx=Lx+2;i<=n;i++,px++,pdiagL++,pLx++)
    if (*px > 0.0)    /* x[i] = 1 */
      {
	if  (*best < *pdiagL-*pLx) 
	  {
	    *best=*pdiagL-*pLx; 
	    *bestarg=i; 
	  };
      }
    else        /* x[i] = -1 */
      {
	if (*best < *pdiagL+*pLx)
	  {
	    *best=*pdiagL+*pLx;
	    *bestarg=i;
	  };
      };
}
