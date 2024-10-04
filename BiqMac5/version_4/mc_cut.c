/*
 *
 * translation from mc_cut.m (f. rendl)
 * dec 2004, a.wiegele
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"

int mc_cut(L,n,v,k,trials,con,eq,rhs,bg,sp,co,pcut,pcost)
    /* input */
    double *L;
    int n;
    double *v;
    int k;
    int trials;
    int con;    /* number of additional constraints of the problem */
    int *eq;
    double *rhs;
    int *bg;
    int *sp;
    double *co;

    /* output */
    double **pcut; 
    double *pcost;

{
  /* vars for call to blas- and lapack-routines */
  double alpha,beta;
  int inc=1;
  char tr='t',up='U';
  
  /* some other var's */
  int i,j,ret;
  double feas, *x, *y, *r;
  int seed;
  double *p;
  int valid=0;
  
  /**********************************
   *                                *
   * initialization                 *
   *                                *
   *                                *
   **********************************/
  
  seed = 10000;    /* choose a seed value */
  srand(seed);     /* initialize random number generator */
  
  /* allocate space for x, y, *pcut,r  */
  x=(double *)malloc((n+1)*sizeof(double));
  y=(double *)calloc((n+1),sizeof(double));
  r=(double *)malloc((k+1)*sizeof(double));
  if ((x == NULL)||(y == NULL)||(r == NULL))
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };

  *pcost=-1.0e10;
  
  /************************
   *                      *
   * main loop            *
   *                      *
   ************************/  
  
  for (i=1;i<=trials;i++)
    { 
      for (j=1,p=r+1;j<=k;j++,p++)     /* r= rand(k,1)-.5 */
	*p = ( (double)rand() / (double)(RAND_MAX) )-0.5;
      /* y=v'*r */
      alpha=1.0;beta=0.0;

	  // glej: http://www.netlib.org/lapack/explore-3.1.1-html/dgemv.f.html
	  //	y := alpha*A'*x + beta*y
	  // v našem: y+1 := alpha*v'*(r+1) + beta*(y+1)
	  // v je dimenzije k*n

      dgemv_(&tr,&k,&n,&alpha,v,&n,r+1,&inc,&beta,y+1,&inc);
      
      /* y = sign(y) */
      for (j=1,p=y+1;j<=n;j++,p++)
	*p=(*p<0.0 ? -1.0 : 1.0 );

      dcopy_(&n,y+1,&inc,x+1,&inc);  /* x=y; */

      if (con > 0)
	{
	  /* check if the found solution is feasible 
	     and if it is, compute the cut-value */
	  
	  double *Lx;
	  Lx=(double *)calloc((n+1),sizeof(double));
	  if (Lx == NULL)
	    {
	      printf("Storage allocation failed!\n");
	      exit(10);
	    };
	  
	  alpha=1.0;beta=0.0;

	  dspmv_(&up,&n,&alpha,L,x+1,&inc,&beta,Lx+1,&inc); /* Lx=L*x */
	  feas=ddot_(&n,x+1,&inc,Lx+1,&inc); /* cost = x'*Lx */

	  if (feas > *pcost) {
	    if (valid_cut(n,x,con,eq,rhs,bg,sp,co))
	      {
		*pcost=feas;
		dcopy_(&n,x+1,&inc,*pcut+1,&inc);  /* cut=x */
		valid=1;
	      }
	    else {
	      int i1, k1;
	      double vr, vf;
	      for (j=1,p=x+1;j<=n;j++,p++) {
		for(i1=0;i1<con;i1++) {
		  vr = rhs[i1];
		  vf = 0.0;
		  for(k1=bg[i1];k1<bg[i1+1];k1++) {
		    if(l2i1(sp[k1])+1 == j) vf += co[k1]*x[l2i2(sp[k1])+1];
		    else if(l2i2(sp[k1])+1 == j) vf += co[k1]*x[l2i1(sp[k1])+1];
		    else vr -= co[k1]*x[l2i1(sp[k1])+1]*x[l2i2(sp[k1])+1];
		  }
		  if(!Zero(vf) && vr*eq[i1] >= 0) *p += vr/vf;
		}
		*p=(*p<0.0 ? -1.0 : 1.0 );
	      }

	      alpha=1.0;beta=0.0;               
	      dspmv_(&up,&n,&alpha,L,x+1,&inc,&beta,Lx+1,&inc); /* Lx=L*x */
	      feas=ddot_(&n,x+1,&inc,Lx+1,&inc); /* cost = x'*Lx */
	      
	      if (feas > *pcost && valid_cut(n,x,con,eq,rhs,bg,sp,co))
		{
		  *pcost=feas;
		  dcopy_(&n,x+1,&inc,*pcut+1,&inc);  /* cut=x */
		  valid=1;
		}
	    }
	  }
	  free(Lx);

	}
      else // if (con==0)  
	{
	  /* NOTE: mc_1opt does not work yet for handling
	     constraints and so it happens, that a feasible solution gets
	     'destroyed' by mc_1opt */
	  ret=mc_1opt(L,n,&x,&feas);     /* improve cut */
      
	  if (feas > *pcost)
	    {
	      *pcost=feas;
	      dcopy_(&n,x+1,&inc,*pcut+1,&inc);  /* cut=x */
	    };
	};
    }; /* end of for loop */
  

  /* free and return */
  free(r);free(y);free(x);
  
  if (con > 0)
    return(valid);
  else 
    return(1);
}
