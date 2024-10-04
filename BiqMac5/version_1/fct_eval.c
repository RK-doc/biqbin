/* 
 * 
 * evaluates max-cut function with A(X) <= b and X in Elliptope
 * f = <L -A'(gamma),X>  + b'*gamma    function value
 * g = b - A(X)                        subgradient
 * call:  [f, x, g] = fct_eval( gamma, L, b, T); 
 *
 * translation from fct_eval.m (f. rendl)
 * dec 2004, a.wiegele
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"


int fct_eval(n,m,gamma,mu,con,eq,rhs,bg,sp,co,L,T,pf,pX,pg,ph)
    /* input */
    int n,m;
    double *gamma;
    double *mu;
    /* constraints */
    int con;
    int *eq;
    double *rhs;
    int *bg;
    int *sp;
    double *co;
    /* */
    double *L;
    int *T;
    /* output */
    double *pf;//zgornja meja!! (ub oz. bound_ub)
    double **pX;
    double **pg;
    double **ph;
{
    /* vars for call to blas- and lapack-routines */
    int inc=1;
    double alpha;

    /* and some other vars... */
    int nn,i;
    double *L0;
    double *Atg;
    double *y;
    double *b;
    double phi;
    /* vars for the additional constraints */
    double *Btm;
    double *d;

    nn=n*(n+1)/2;

    y=(double *)malloc((n+1)*sizeof(double));
    if (y == NULL)
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };

    L0=(double *)malloc(nn*sizeof(double));
    if (L0 == NULL)
      {
		printf("Storage allocation failed!\n");
		exit(10);
      };

    dcopy_(&nn,L,&inc,L0,&inc); /* L0=L */

    if (con>0)
      { 
		Btm=(double *)malloc(nn*sizeof(double));
		if (Btm == NULL)
		  {
			printf("Storage allocation failed!\n");
			exit(10);
		  };

		op_bt(n,con,eq,bg,sp,co,mu,&Btm);
	
		alpha=-1.0;
		daxpy_(&nn,&alpha,Btm,&inc,L0,&inc); /* L0=L-Bt(mu) */

		d=(double *)malloc((con+1)*sizeof(double));
		if (d == NULL)
		  {
				printf("Storage allocation failed!\n");
				exit(10);
		  };

		for (i=1;i<=con;i++)
		  if (eq[i-1]<1)
			d[i]=rhs[i-1];
		  else
			d[i]=-rhs[i-1];
		free(Btm);
      };


    if (m>0)
        {
            Atg=(double *)malloc(nn*sizeof(double));
            b=(double *)malloc((m+1)*sizeof(double));
            if ((Atg == NULL)||(b == NULL))
                {
                    printf("Storage allocation failed!\n");
                    exit(10);
                };
            for (i=1;i<=m;i++)           /* b=-ones(m,1) */    
              b[i]=-1.0;

            mc_at(n,m,MAX_INEQ,T,gamma,&Atg);     /* A^T(gamma) */

            alpha=-1.0;
            daxpy_(&nn,&alpha,Atg,&inc,L0,&inc); /* L0=L0-Atg */
			mc_psdpk(L0,n,1,pX,&y,&phi);
	    
			*pf=phi + ddot_(&m,b+1,&inc,gamma+1,&inc);  /* function value */

            /* g=b-ax */
            mc_a(n,m,MAX_INEQ,T,*pX,pg);                                      
            alpha=-1.0;
            daxpy_(&m,&alpha,b+1,&inc,*pg+1,&inc);   /* g=-b+g; */
            dscal_(&m,&alpha,*pg+1,&inc);            /* g=-g;   */

            free(b);free(Atg);
        }
    else /* m=0 */
      {
	    mc_psdpk(L0,n,1,pX,&y,&phi);

	    *pf=phi;    /* function value e'*y */
        };

    if (con>0)
      {
		/* h=d-bx */
		op_b(n,con,eq,bg,sp,co,*pX,ph);
		alpha=-1.0;
		daxpy_(&con,&alpha,d+1,&inc,*ph+1,&inc);
		dscal_(&con,&alpha,*ph+1,&inc);

		(*pf)+=ddot_(&con,d+1,&inc,mu+1,&inc);

		free(d);
      };

    free(L0);free(y);

    return(0);
}
