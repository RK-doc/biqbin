/*
 *
 * solve lambda eta problem for fixed beta, G, gamma
 * [d, lam, eta, t] = lam_eta( beta, G, gamma, t);  
 * note: eta never used in parent function!
 *
 * translation from lam_eta.m (f. rendl)
 * dec 2004, a. wiegele
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"


int lam_eta(k,m,beta,G,gamma,pt,pd,plam)
    /* input */
    int k;        /* number of subgrads */
    int m;        /* number of constr */
    double *beta;
    double *G;
    double *gamma;
    /* input/output */
    double *pt;
    /* output */
    double **pd;
    double **plam;
{
  /* vars for call to blas- and lapack-routines */
  int inc=1;
  double alpha,betab; 
  char tr='T',nt='N',up='U'; 
  
  /* some other var's */
  int ret,i,j;
  int done,lam_cnt;
  double *Q,*Qdmy;
  double *c;
  double *tmp,*eta,*d_curr;
  double dir_pre,dir_curr;
  double *p;
  
  
  /* note: allocation of memory for d and lam is done in bdl_mc2! */
  /*    *plam=(double *)malloc((k+1)*sizeof(double));
        *pd=(double *)malloc((m+1)*sizeof(double));   */
  
  /**********************************
   *                                *
   * initialization                 *
   *                                *
   **********************************/
  
  /* allocate space for matrix Qdmy and Q */
  Qdmy=(double *)malloc(k*k*sizeof(double));
  Q=(double *)malloc(k*(k+1)/2*sizeof(double));
  if ((Qdmy == NULL)||(Q == NULL))
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };
  
  /* Q=t*G'*G */
  betab=0.0;
  dsyrk_(&up,&tr,&k,&m,pt,G,&m,&betab,Qdmy,&k);
  /* transform Qdmy to Q in packed format */
  for (i=1,p=Q;i<=k;i++)
    for (j=1;j<=i;j++,p++)
      *p=Qdmy[ijtok(j,i,k)];  // Q[ijtokp(j,i,k)]=Qdmy[ijtok(j,i,k)];


  eta=(double *)calloc((m+1),sizeof(double));
  c=(double *)malloc((k+1)*sizeof(double));
  tmp=(double *)malloc((m+1)*sizeof(double));
  d_curr=(double *)malloc((m+1)*sizeof(double));
  if ((eta == NULL)||(c == NULL)||(tmp == NULL)||(d_curr == NULL))
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };

  /******************
   *                *
   * main loop      *
   *                *
   ******************/
  
  done=1;
  dir_pre=0;
  lam_cnt=0;

  while (done == 1)
    {
      lam_cnt++;
      /* c=-t*G'*eta + beta; */
      dcopy_(&k,beta+1,&inc,c+1,&inc);  /* c=beta; */
      alpha=-*pt;betab=1.0;
      dgemv_(&tr,&m,&k,&alpha,G,&m,eta+1,&inc,&betab,c+1,&inc);
      
      ret=solvelambda(Q,c,k,plam,1); /* 1 means no ouput in solvelambda */
      
      /* tmp = G*lam; */
      alpha=1.0;betab=0.0;
      dgemv_(&nt,&m,&k,&alpha,G,&m,*plam+1,&inc,&betab,tmp+1,&inc);

      dcopy_(&m,tmp+1,&inc,eta+1,&inc); /* eta=tmp; */
      alpha=-1.0/(*pt);
      daxpy_(&m,&alpha,gamma+1,&inc,eta+1,&inc);  /* eta=-1/(*pt)*gamma+tmp */ 

      for (i=1,p=eta+1;i<=m;i++,p++)
	if ((*p)< 0.0)
	  (*p)=0.0;

      dcopy_(&m,eta+1,&inc,d_curr+1,&inc);        /* d_curr=eta; */
      alpha=-1.0;
      daxpy_(&m,&alpha,tmp+1,&inc,d_curr+1,&inc); /* d_curr=-tmp+eta; */
      dscal_(&m,pt,d_curr+1,&inc);                /* d_curr=(*pt)*d_curr; */

      dir_curr=dnrm2_(&m,d_curr+1,&inc);

      if (Abs(dir_pre-dir_curr)/(1+dir_curr) < EPSILON)
	done=0;
      if (lam_cnt > 50)
	{
	  done=0;
	  //printf(" *warning* max iter in lam_eta %8.4f\n",Abs(dir_pre-dir_curr));
	  (*pt)*=.95;
	};
      dir_pre=dir_curr;
    };
  dcopy_(&m,d_curr+1,&inc,*pd+1,&inc);
  
  /* free memory and return */
  free(d_curr);free(tmp);
  free(c);free(eta);
  free(Q);free(Qdmy);      
  return(0);  
}
