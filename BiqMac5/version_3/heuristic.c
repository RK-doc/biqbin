/*
 *
 * translation from heuristic.m (f. rendl)
 * dec 2004, a.wiegele
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"

int heuristic(L,X,n,pxh,pfh,con,eq,rhs,bg,sp,co)
    /* input */
    double *L;         /* Laplace mat */
    double *X;         /* primal sol */
    int n;             /* dim */
    /* input/output */
    double **pxh;      /* entry: cut vector, exit: (improved) cut */
    /* output */
    double *pfh;       /* value of cut xh */
    /* input */
    int con;           /* number of additional constraints of the problem */
    int *eq;
    double *rhs;
    int *bg;
    int *sp;
    double *co;
{

  /* vars for call to blas- and lapack-routines */
  double alpha,beta;
  int inc=1;
  int nn, info;
  char up='U';

  /* some other var's */
  int i,j,ret,done;
  int its,seed;
  double *Xs, *v;
  double f0, *cut;
  double *p;
  int valid=0;
 
  /**********************************
   *                                *
   * initialization                 *
   *                                *
   *                                *
   **********************************/


  seed = 10000;
  srand(seed);

  its=(2*n>500 ? 2*n : 500);
  nn=n*(n+1)/2;

  /* allocate space for Xs */
  Xs=(double *)malloc(nn*sizeof(double));
  v=(double *)malloc(n*n*sizeof(double));
  cut=(double *)malloc((n+1)*sizeof(double));
  if ((Xs == NULL)||(v == NULL)||(cut == NULL))
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };
  
  /* fh = xh'*L*xh */
  alpha=1.0; beta=0.0;
  dspmv_(&up,&n,&alpha,L,*pxh+1,&inc,&beta,cut+1,&inc);
  *pfh=ddot_(&n,*pxh+1,&inc,cut+1,&inc);

  
  /************************
   *                      *
   * main loop            *
   *                      *
   ************************/  
  
  for (done=0;done<2;)
    {
      done+=1;
      
      /* Xs = alpha*X + (1-alpha)*xh*xh' */
      for (i=0,p=Xs;i<nn;i++,p++)
	(*p)=0.0;
      alpha=.3+.6*( (double)rand()/(double)(RAND_MAX) );
      daxpy_(&nn,&alpha,X,&inc,Xs,&inc);   /* Xs = alpha*X */
      alpha=1-alpha;
      dspr_(&up,&n,&alpha,*pxh+1,&inc,Xs); /* Xs = Xs + (1-alpha)*xh*xh' */
      
      /* Xs = chol(Xs) */
      dpptrf_(&up,&n,Xs,&info);
      if (info < 0)
	{
	  printf("\n heuristic.c: problems in call to dpptrf! %d\n",info);
	  exit(10);
	};
      if (info > 0)
	{
	  printf("\n heuristic.c: problems computing cholesky! %d\n",info);
	  exit(10);
	};
      /* have to move cholesky factor Xs to matrix v (nxn) for the call to mc_cut */
      
      /* ...weil mc_cut hat als input eine (k x n) matrix v mit Xs =v'*v, 
	 die nicht notwendigerweise upper tri ist */
      for (i=0,p=v;i<n*n;i++,p++)
	*p=0.0;
      for (i=1,p=Xs;i<=n;i++)
	for (j=1;j<=i;j++,p++)
	  {
	    v[ijtok(j,i,n)]=*p;
	    /*       v[ijtok(i,j,n)]=Xs[ijtokp(i,j,n)];
		     if (j != i)
		     v[ijtok(j,i,n)]=0.0;*/
	  };

      ret=mc_cut(L,n,v,n,its,con,eq,rhs,bg,sp,co,&cut,&f0);

      if (ret > 0) 
	{  /* in case of constraints, the cut is valid. 
	      in case of no constraints, ret=1 anyway */

	  if (f0 > *pfh + 1.0e-10)
	    {
	      *pfh=f0;
	      /* printf("   better cut:         %8.4f\n",*pfh); */
	      dcopy_(&n,cut+1,&inc,*pxh+1,&inc);
	      valid=1;
	      done=0; /* do another round */
	    };
	}
      else
	{ /* the extra constraints are not valid for this cut.
	     do one more try with a different matrix Xs */
	  /* Xs = .5*X + .5*cut*cut' */
          //done=0;
	};
    };
  
  /* free memory */
  free(cut);free(v);free(Xs);
  return(valid);
}
