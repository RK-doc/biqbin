/* 
 * 
 * input: L (Laplace matrix)
 *        initial set of constraints: defined by (b,T,hash)
 *        initial guess for gamma  
 * perform (at most) itmax iterations for given data
 * call: [ fopt,gamma,xi,t,x] = bdl_mc2( L, b, T, hash, gamma, itmax,t,xi,Xin);
 *
 * translation from bdl_mc2.m (f. rendl)
 * jan 2005, a.wiegele
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"


int bdl_mc2(n,m,L,T,hash,con,eq,rhs,bg,sp,co,best,pgamma,pmu,itmax,pt,pxi,pX,pfopt,px,pkin,printlevel,gap_relax)
     int n;            /* dim(L) */
     int m;            /* length(b), i.e. number of constraints */
     double *L;
     double best;
     int *T;
     unsigned int *hash;
     /* constraints */
     int con;
     int *eq;
     double *rhs;
     int *bg;
     int *sp;
     double *co;
     /* */
     double **pgamma;
     double **pmu;
     int itmax;
     double *pt;
     double **pxi;
     double **pX;    
     double *pfopt;
     double **px;   /* sym mat of dim n */
     int *pkin;       /* size(Xin,2) */
     int printlevel;
     double gap_relax;
{
    /* vars for call to blas- and lapack-routines */
    int inc=1;
    double alpha,betab;
    char tr='T';

    /* and some other vars */
    double f,*g;
    double *F,*G;
    double *axi;
    double *b;

    int ret,i,cnt,j;
    int k,nn;

    double *beta;
    double *dgamma,*lam,*gamma_test;
    double ftest,*xtest,*gtest;
    double *gv;

    int arglmax;
    double lmax;

    /* for additional constraints */
    double *h,*mu_test,*htest,*H,*d,*dmu;
    double *bxi,*hv;

    /* internal data:
       gamma ... current point
       X = (x1,...) ...collection of interesting trial points
       G = (g1,...) ...corresponding subgrads (gi= b-A(xi))
       F = (f1,...) ...function values for x1,...
       [H = (h1,...) ...subgrads part II (hi = d-B(xi))]
    */

    /* allocate memory */
    nn=n*(n+1)/2;
    F=(double *)malloc((MAX_BDL+1)*sizeof(double));
    G=(double *)malloc(m*MAX_BDL*sizeof(double)); /* m by k matrix */
    if ((G == NULL)||(F == NULL))
        {
            printf("Storage allocation failedGF!\n");
            exit(10);
        };

    beta=(double *)malloc((MAX_BDL+1)*sizeof(double));
    if (beta == NULL)
        {
            printf("Storage allocation failedbeta!\n");
            exit(10);
        };

    if (con > 0)
      {
	H=(double *)malloc(con*MAX_BDL*sizeof(double));
	d=(double *)malloc((con+1)*sizeof(double));
	if ((H == NULL)||(d == NULL))
	  {
            printf("Storage allocation failedHd!\n");
            exit(10);
	  };
	for (i=1;i<=con;i++)
	  if (eq[i-1]<1)
	    d[i]=rhs[i-1];
	  else
	    d[i]=-rhs[i-1];
      };

    /*******************
     *                 *
     * initialization  *
     *                 *
     *******************/

    /* start: initial function and subgrad evaluation at gamma */
    g=(double *)malloc((m+1)*sizeof(double));
    b=(double *)malloc((m+1)*sizeof(double));
    if ((g == NULL)||(b == NULL))
        {
            printf("Storage allocation failedgb!\n");
            exit(10);
        };

    for (i=1;i<=m;i++)     /* b=-ones(m,1); */
        b[i]=-1.0;

    if (con > 0)
      {
	h=(double *)malloc((con+1)*sizeof(double));
	if (h == NULL)
	  {
            printf("Storage allocation failedh!\n");
            exit(10);
	  };
      };

    /*ret=fct_eval(n,m,*pgamma,L,T,&f,px,&g);*/
    ret=fct_eval(n,m,*pgamma,*pmu,con,eq,rhs,bg,sp,co,L,T,&f,px,&g,&h);

    /*
     * initial F,G,X (first kin elements are given by
     * input pX, last two elements are associated with xi
     * and the last element with px from the last fct_eval 
     */

    *pfopt=f;     /* initial function */
    F[(*pkin)+1]=traceprod(L,*pxi,n);
    F[(*pkin)+2]=traceprod(L,*px,n);

    axi=(double *)malloc((m+1)*sizeof(double));
    if (axi == NULL)
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };

    ret=mc_a(n,m,MAX_INEQ,T,*pxi,&axi);            /* A(xi) */
    dcopy_(&m,b+1,&inc,G+(*pkin)*m,&inc);          /* G(:,kin+1)=b; */
    alpha=-1.0;
    daxpy_(&m,&alpha,axi+1,&inc,G+(*pkin)*m,&inc); /* G(:,kin+1)=b-axi */

    dcopy_(&m,g+1,&inc,G+((*pkin)+1)*m,&inc);      /* G(:,kin+2)=g */

    if (con>0)
      {
	bxi=(double *)malloc((con+1)*sizeof(double)); 
	if (bxi == NULL)
	  {
            printf("Storage allocation failed!\n");
            exit(10);
	  };

	ret=op_b(n,con,eq,bg,sp,co,*pxi,&bxi);             /* B(xi) */

	dcopy_(&con,d+1,&inc,H+(*pkin)*con,&inc);          /* H(:,kin+1)=d; */
	alpha=-1.0;
	daxpy_(&con,&alpha,bxi+1,&inc,H+(*pkin)*con,&inc); /* H(:,kin+1)=d-bxi; */

	dcopy_(&con,h+1,&inc,H+((*pkin)+1)*con,&inc);      /* H(:,kin+2)=h; */

      };

    /* (*pX) = [(*pX) xi x(:)] */
    dcopy_(&nn,*pxi,&inc,(*pX)+(*pkin)*nn,&inc);    
    dcopy_(&nn,*px,&inc,(*pX)+((*pkin)+1)*nn,&inc);

    for (i=1;i<=*pkin;i++)
        {
            F[i]=traceprod(L,*pX+(i-1)*nn,n);
            ret=mc_a(n,m,MAX_INEQ,T,*pX+(i-1)*nn,&axi);  /* A(xi) */

            dcopy_(&m,b+1,&inc,G+(i-1)*m,&inc);          /* G(:,i)=b; */
            alpha=-1.0;
            daxpy_(&m,&alpha,axi+1,&inc,G+(i-1)*m,&inc); /* G(:,i)=b-axi */

        };

    if (con > 0)
      {
	for (i=1;i<=*pkin;i++)
	  {
	    ret=op_b(n,con,eq,bg,sp,co,*pX+(i-1)*nn,&bxi);   /* B(xi) */
	    dcopy_(&con,d+1,&inc,H+(i-1)*con,&inc);          /* H(:,kin+1)=d; */
	    alpha=-1.0;
	    daxpy_(&con,&alpha,bxi+1,&inc,H+(i-1)*con,&inc); /* H(:,kin+1)=d-bxi; */
	  };
      };

    k=(*pkin)+2;  /* current bundle size */
    if (k > MAX_BDL)
        {
            printf("Bundle size too large! adjust MAX_BDL\n");
            exit(10);
        };

    if (printlevel > 4)
        {
            printf(" %3.0d  %12.3f %12.3f      %7.5f %3.0d\n",0,*pfopt,*pfopt,*pt,k);
        };

    /****************
     *              *
     * main loop    *
     *              *
     ****************/

    /* allocate memory for dgamma, lam, gamma_test */
    dgamma=(double *)malloc((m+1)*sizeof(double));
    lam=(double *)malloc((MAX_BDL+1)*sizeof(double));
    gamma_test=(double *)malloc((m+1)*sizeof(double));
    if ((dgamma == NULL)||(lam == NULL)||(gamma_test == NULL))
        {
            printf("Storage allocation failed! dgamma,lam,gamma_test\n");
            exit(10);
        };
    xtest=(double *)malloc(nn*sizeof(double));
    gtest=(double *)malloc((m+1)*sizeof(double));
    gv=(double *)malloc((m+1)*sizeof(double));
    if ((xtest == NULL)||(gtest == NULL)||(gv == NULL))
        {
            printf("Storage allocation failed! xtest,gtest,gv\n");
            exit(10);
        };

    if (con > 0)
      {
	htest=(double *)malloc((con+1)*sizeof(double));
	mu_test=(double *)malloc((con+1)*sizeof(double));
	dmu=(double *)malloc((con+1)*sizeof(double));
	hv=(double *)malloc((con+1)*sizeof(double));
	if ((htest == NULL)||(mu_test == NULL)||(dmu == NULL)||(hv == NULL))
	  {
            printf("Storage allocation failed! htest,mu_test\n");
            exit(10);
	  };
      };

    for (cnt=1;cnt<=itmax;cnt++)
        {
            /* 
             * compute lam, eta and dgamma 
             */

            /* beta = -F - G'*gamma; */
            dcopy_(&k,F+1,&inc,beta+1,&inc);  /* beta = F */
            alpha=-1.0;betab=-1.0;            /* beta = -G'*gamma-beta */
	    if (m > 0)
	      dgemv_(&tr,&m,&k,&alpha,G,&m,*pgamma+1,&inc,&betab,beta+1,&inc); /* NOTE: G is m by k, but memory m by MAX_BDL! */
	    else
	      dscal_(&k,&alpha,beta+1,&inc);

	    if (con > 0)
	      {
		//for (i=1;i<=con;i++)
		//printf("mu=%lf\n",(*pmu)[i]);
		/* beta = -G'*gamma-H'*mu-F */
		alpha=-1.0;betab=1.0;
		dgemv_(&tr,&con,&k,&alpha,H,&con,*pmu+1,&inc,&betab,beta+1,&inc);
	      };
	    //for (i=1;i<=(k<20 ? k : 20);i++)
	    //printf("beta=%lf\n",beta[i]);

	    if (con == 0) /* con == 0 */
	      ret=lam_eta(k,m,beta,G,*pgamma,pt,&dgamma,&lam);
	    else /* con > 0 */
	      {
		double *GH,*gamma_mu,*dgamma_mu;
		GH=(double *)malloc((con+m)*MAX_BDL*sizeof(double));
		gamma_mu=(double *)malloc((con+m+1)*sizeof(double));
		dgamma_mu=(double *)malloc((con+m+1)*sizeof(double));
		if ((GH == NULL)||(gamma_mu == NULL)||(dgamma_mu == NULL))
		  {
		    printf("Storage allocation failed!\n");
		    exit(10);
		  };
		/* GH = [G;H] ...(m+con) by k */
		/*alpha=m*k;
		betab=con*k;
		dcopy_(&alpha,G,&inc,GH,&inc);
		dcopy_(&betab,H,&inc,GH+m*k,&inc);*/
		for (i=1;i<=k;i++)
		  {
		    dcopy_(&m,G+(i-1)*m,&inc,GH+(i-1)*(m+con),&inc);
		    dcopy_(&con,H+(i-1)*con,&inc,GH+(i-1)*(m+con)+m,&inc);
		  };
		/* gamma_mu =  [gamma;mu]; */
		dcopy_(&m,*pgamma+1,&inc,gamma_mu+1,&inc);
		dcopy_(&con,*pmu+1,&inc,gamma_mu+m+1,&inc);

		ret=lam_eta(k,m+con,beta,GH,gamma_mu,pt,&dgamma_mu,&lam);
		//for(i=1;i<=(k<20 ? k : 20);i++)
		//printf("gamma_mu=%lf, dgamma_nu=%lf \n",gamma_mu[i],dgamma_mu[i]);
		dcopy_(&m,dgamma_mu+1,&inc,dgamma+1,&inc);
		dcopy_(&con,dgamma_mu+m+1,&inc,dmu+1,&inc);

		free(dgamma_mu);free(gamma_mu);free(GH);
	      };
	    //for (i=1;i<=20;i++)
	    //printf("beta=%lf gamma=%lf  dgamma=%lf\n",beta[i],(*pgamma)[i],dgamma[i]);
            /*
             * make a step in the direction of -dgamma:
             * gamma_test=gamma+dgamma;
             */
            dcopy_(&m,*pgamma+1,&inc,gamma_test+1,&inc);  /* gamma_test=gamma */
            alpha=1.0;                                    /* gamma_test=dgamma+gamma_test */
            daxpy_(&m,&alpha,dgamma+1,&inc,gamma_test+1,&inc);

	    if (con > 0)
	      {
		/* make a step in the direction of -dmu:
		 * mu_test=mu+dmu;
		 */
		dcopy_(&con,*pmu+1,&inc,mu_test+1,&inc);       /* mu_test=mu; */
		alpha=1.0;                                     /* mu_test=dmu+mu_test */
		daxpy_(&con,&alpha,dmu+1,&inc,mu_test+1,&inc);
	      };

            /*
             * compute function at gamma_test 
             */
            
            /*ret=fct_eval(n,m,gamma_test,L,T,&ftest,&xtest,&gtest);*/
	    ret=fct_eval(n,m,gamma_test,mu_test,con,eq,rhs,bg,sp,co,L,T,&ftest,&xtest,&gtest,&htest);
 
	    /* check if the bound is below the current best cut value */
	    if (ftest < best+gap_relax)
	      {
		*pfopt=ftest;
		if (printlevel > 4)
		  printf("...bound below cut+%f => stop bundle and return to master program \n",gap_relax);
		break;
	      }

            /* xi = .1*xtest(:) + .9*xi; */
            alpha=0.9;
            dscal_(&nn,&alpha,*pxi,&inc);  /* xi=.9*xi */
            alpha=0.1;                       /* xi=.1*xtest+xi */
            daxpy_(&nn,&alpha,xtest,&inc,*pxi,&inc); 

            ret=mc_a(n,m,MAX_INEQ,T,*pxi,&axi);   /* A(xi) */

            /* gv = b-axi;    subgradient at xi */
            dcopy_(&m,b+1,&inc,gv+1,&inc);          /* gv=b; */
            alpha=-1.0;
            daxpy_(&m,&alpha,axi+1,&inc,gv+1,&inc); /* gv=b-axi */

	    if (con > 0)
	      {
		/* hv = d-bxi; */
		dcopy_(&con,d+1,&inc,hv+1,&inc);
		alpha=-1.0;
		daxpy_(&con,&alpha,bxi+1,&inc,hv+1,&inc);
	      };

            /*
             * update f and gamma, if progress is made 
             */

            /* lmax = max(lam); */
            dmax_val(k,lam,&arglmax,&lmax);

            if (ftest < *pfopt) 
                {
                    /* serious step */

                    j=1;
                    for (i=1;i<=k;i++)
                        {
                            if (lam[i] > 0.01*lmax)
                                {
                                    if (i > j)
                                        {
                                            dcopy_(&m,G+(i-1)*m,&inc,G+(j-1)*m,&inc);
                                            dcopy_(&nn,*pX+(i-1)*nn,&inc,*pX+(j-1)*nn,&inc);
                                            F[j]=F[i];
					    if (con > 0)
					      dcopy_(&con,H+(i-1)*con,&inc,H+(j-1)*con,&inc);
                                        }
                                    j++;
                                };
                        }

		    if (j+1 > MAX_BDL)
		      {
			printf("Bundle size too large! adjust MAX_BDL\n");
			exit(10);
		      };

                    dcopy_(&m,gv+1,&inc,G+(j-1)*m,&inc);
                    dcopy_(&m,gtest+1,&inc,G+j*m,&inc);

                    dcopy_(&nn,*pxi,&inc,*pX+(j-1)*nn,&inc);
                    dcopy_(&nn,xtest,&inc,*pX+j*nn,&inc);

                    F[j]=traceprod(L,*pxi,n);
                    F[j+1]=traceprod(L,xtest,n);

		    if (con > 0)
		      {
			dcopy_(&con,hv+1,&inc,H+(j-1)*con,&inc);
			dcopy_(&con,htest+1,&inc,H+j*con,&inc);
			dcopy_(&con,mu_test+1,&inc,*pmu+1,&inc);
		      };

                    dcopy_(&m,gamma_test+1,&inc,*pgamma+1,&inc);
                    *pfopt=ftest;
                    dcopy_(&nn,xtest,&inc,*px,&inc);   /* x=xtest */
                    (*pt)*=1.01;
                    k=j+1;
                }
            else
                {
                    /* null step */
 
                    j=1;
                    for (i=1;i<=k-1;i++)
                        {
                            if (lam[i] > 0.01*lmax)
                                {
                                    if (i > j)
                                        {
                                            dcopy_(&m,G+(i-1)*m,&inc,G+(j-1)*m,&inc);
                                            dcopy_(&nn,*pX+(i-1)*nn,&inc,*pX+(j-1)*nn,&inc);
                                            F[j]=F[i];
					    if (con > 0)
					      dcopy_(&con,H+(i-1)*con,&inc,H+(j-1)*con,&inc);
                                        }
                                    j++;
                                };
                        };

		    if (j+2 > MAX_BDL)
		      {
			printf("Bundle size too large! adjust MAX_BDL\n");
			exit(10);
		      };

                    /* on last position the current last element
                       (old: k, new: j+2) */
                    dcopy_(&m,G+(k-1)*m,&inc,G+(j+1)*m,&inc);
                    dcopy_(&nn,*pX+(k-1)*nn,&inc,*pX+(j+1)*nn,&inc);
                    F[j+2]=F[k];
		    if (con > 0)
		      dcopy_(&con,H+(k-1)*con,&inc,H+(j+1)*con,&inc);

                    /* and in between gtest and gv, xtest and xv */
                    dcopy_(&m,gtest+1,&inc,G+(j-1)*m,&inc);
                    dcopy_(&m,gv+1,&inc,G+j*m,&inc);

                    dcopy_(&nn,xtest,&inc,*pX+(j-1)*nn,&inc);
                    dcopy_(&nn,*pxi,&inc,*pX+j*nn,&inc);

                    F[j]=traceprod(L,xtest,n);
                    F[j+1]=traceprod(L,*pxi,n);

		    if (con > 0)
		      {
			dcopy_(&con,htest+1,&inc,H+(j-1)*con,&inc);
			dcopy_(&con,hv+1,&inc,H+j*con,&inc);
		      };

                    (*pt)/=1.02;
                    k=j+2;
                };
            
            /* display some info */
            if (printlevel > 4)
                {
                    printf(" %3.0d  %12.3f %12.3f      %7.5f %3.0d\n",cnt,*pfopt,ftest,*pt,k);
                };
        };

    *pkin=k;

    if (printlevel > 5)
        {
            printf(" output from bdl_mc2 before returning to bdl_main:\n");
            printf(" fopt(=f)=%8.4f t=%8.4f bdlsize=%d\n",*pfopt,*pt,*pkin);
            printf(" gamma:\n");
            for (i=1;i<=(m<20 ? m : 20);i++)
                printf(" %8.4f ",(*pgamma)[i]);
            printf("\n");
        };
    
    /* free memory and return */
    if (con > 0)
      {
	free(bxi);
	free(hv);free(dmu);free(mu_test);free(htest);
      };
    free(gv);free(gtest);free(xtest);
    free(gamma_test);free(lam);free(dgamma);
    free(axi);
    free(b);free(g);free(beta);
    if (con > 0)
      {
	free(h);
	free(d);free(H);
      };
    free(G);free(F);
    return(0);
}
