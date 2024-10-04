/*
 * bundle method: main routine
 *
 * used as a subroutine for G. Rinaldi's bab_for_klu.c 
 *
 * see  I. Fischer, G. Gruber, F. Rendl, R. Sotirov
 * "Computational experience with a bundle approach for 
 * semidefinite cutting plane relaxations of Max-Cut and
 * Equipartition", In: Mathematical Programming,
 * Serie B, vol.105 (2006), pp. 451-469.
 * translation from bdl_main.m (f. rendl)
 * oct 2005, a. wiegele
 *
 */

/* note: 
for easier indexing, i allocate for vectors of length k 
(k+1)*sizeof(type) memory and ignore the 0-element;
so v[i] is the i-th element. (that's why in the blas/lapack routines 
i have the arguments v+1.)
for sym. matrices of dim k, i allocate the space they need, i.e. k*(k+1)/2 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "declarations.h"

int bdl_main(P,L,n,forec,con,eq,rhs,bg,sp,co,pbest,pf,pcut,primalX,pm,T,hash,gamma)
    /* input */
    T_param P;
    double *L;
    int n;
    int forec;
    /* constraints */
    int con;
    int *eq;
    double *rhs;
    int *bg;
    int *sp;
    double *co;
    /* output */
    double *pbest;
    double *pf;//zgornja meja!! (ub oz. bound_ub)
    double **pcut;
    double **primalX;
    /* (in/)output */
    int *pm;
    int T[MAX_INEQ*4];
    unsigned int *hash;
    double *gamma;

{    
    /* vars for call to blas/lapack routines */
    int inc=1;
    double alpha;

    /*  */
    int ret,i,j,cont;

    /* some counters */
    int main_loop,iterates,iterdone;
    /* some numbers */
    int nn,in_it,m3,m2=0;
    /*  */
    double f;
    double *X;

    double f0,f00,*x0,*g0;

    double *x;
    double t;
    double *cut,*cut1,best,best1;

    int bdlsize;   /* number of cols in matrix X ... nn-by-bdlsize */

    double gamma_mean,*gamma_new,gamma_max;
    int arggmax;
    double *pp;

    /* variables needed for additional constraints */
    double *mu;
    double *h0;

	//izpis Laplaceove matrike
    if (P.printlev > 5)
        {
            printf("\n\n matrix L: \n");
            for (j=1;j<=n;j++)
                {
                    for (i=1;i<=n;i++)
                        {
                            if (i>j)
                                printf("%4.2f ",L[ijtokp(j,i,n)]);
                            else
                                printf("%4.2f ",L[ijtokp(i,j,n)]);
                        }
                    printf("\n");
                };  
        };

    /* memory allocation */

    cut=(double *)malloc((n+1)*sizeof(double));
    cut1=(double *)malloc((n+1)*sizeof(double));
    if ((cut == NULL)||(cut1 == NULL))
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };
    /*    T=(int *)malloc(MAX_INEQ*4*sizeof(int)); 
    hash=(unsigned int *)malloc((MAX_INEQ+1)*sizeof(unsigned int));
    gamma=(double *)malloc((MAX_INEQ+1)*sizeof(double));   */
    gamma_new=(double *)malloc((MAX_INEQ+1)*sizeof(double));
    if (gamma_new == NULL)
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };
    g0=(double *)malloc((MAX_INEQ+1)*sizeof(double));
    x0=(double *)malloc(n*(n+1)/2*sizeof(double));

    if ((x0 == NULL)||(g0 == NULL))
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };


    if (con > 0)
      {
		mu=(double *)calloc((con+1),sizeof(double));
		h0=(double *)calloc((con+1),sizeof(double));
		if ((mu == NULL)||(h0 == NULL))
		  {
				printf("Storage allocation failed!\n");
				exit(10);
		  };
      };

    /**********************************
     *                                *
     * initial function evaluation    *
     *                                *
     **********************************/

    nn=n*(n+1)/2;

    if (*pm==0)
        {
             /* initial fct value f0 and initial matrix x0
               (x0 is output of mc_psdpk!) */
            if (P.printlev > 4)
                printf(" ...initial function evaluation\n");
			ret=fct_eval(n,*pm,gamma,mu,con,eq,rhs,bg,sp,co,L,T,&f0,&x0,&g0,&h0);
                    /* note: because of the fact that m=0, gamma,T are not 
                       used in fct_eval. also, output g0 remains unchanged. */

			if (!forec) /* i want to skip the heuristic */
                {
				for (i=1,pp=cut+1;i<=n;i++,pp++)
				(*pp)=1.0; //cut[i]=1.0;

				if (P.printlev > 4)     
						printf(" ...calling rounding heuristic:\n");
				if(!heuristic(L,x0,n,&cut,&best,con,eq,rhs,bg,sp,co)||(best < *pbest)) {
					dcopy_(&n,*pcut,&inc,cut+1,&inc);
							best=*pbest;         /* use cut given by the masterprogram if it's 
												  * better than the one given by the heuristic
												  */
					}
                }
			else 
               {
				dcopy_(&n,*pcut,&inc,cut+1,&inc);
				best=*pbest;
               };
            
            /* check if the bound is already below the best cut. in this case you can return to 
             * bound.c and the node will be fathomed.
             */
  
            if (f0 < best+P.gap_relax)
                {
                    if (P.printlev > 4)
                        printf(" initial bound below best cut!\n");
                    
                    /* copy values to output parameters, free and return */
                    
                    *pbest=best;
                    dcopy_(&n,cut+1,&inc,*pcut,&inc);
                    *pf=f0;
                    dcopy_(&nn,x0,&inc,*primalX,&inc);
                    
				if (con > 0)
				  {
					free(h0);free(mu);
				  };
						free(x0); 
						free(g0);   
						free(gamma_new);
						free(cut1);free(cut);
                    
						return(1);
                }
            

            if (P.printlev > 5)
                {
                    printf(" matrix x0 (after first call to heuristic): \n");
                    for (j=1;j<=n;j++)
                        {
                            for (i=1;i<=n;i++)
                                {
                                    if (i>j)
                                        printf("%lf ",x0[ijtokp(j,i,n)]);
                                    else
                                        printf("%lf ",x0[ijtokp(i,j,n)]);
                                }
                            printf("\n");
                        };  
                };
            
            if (P.printlev > 4)  
                printf(" ...finding initial set of violated ineq.\n");
            /* *pm=0 if called from bound main */
            ret=separation(P,x0,n,pm,T,&hash,&gamma);             
	    

            if (*pm > MAX_INEQ)
                {
                    printf(" number of max. violated inequalities m=%d exceeded. Please adjust MAX_INEQ!\n",*pm);
                    exit(10);
                };
            

            if (P.printlev > 2)
                {
                    printf("           initial bound: %10.2f\n",f0);
                    printf(" current best cut, ratio: %10.2f %6.3f\n",best,(f0/best-1)*100);
                    printf("    initial inequalities: %8.0d\n",*pm);
                };
        }
    else     /*  *pm > 0, i.e. start bundle with given gamma,T,hash */
        {
            f0=*pf;
            best=*pbest;
            dcopy_(&n,*pcut,&inc,cut+1,&inc);      /* cut=given cut */
            dcopy_(&nn,*primalX,&inc,x0,&inc);     /* x0 =given X */
        };

    if (*pm > 0)
      t=0.3*(f0-best)/ddot_(pm,gamma+1,&inc,gamma+1,&inc); 
    else
      t=0.3*(f0-best);
    if (P.printlev > 5)
      printf("initial t=%lf\n",t);
    alpha=t*0.1;
    dscal_(pm,&alpha,gamma+1,&inc); /* gamma=0.1*t*gamma */

    /******************
     *                *
     * main loop      *
     *                *
     ******************/

    x=(double *)malloc(nn*sizeof(double));
    X=(double *)malloc(nn*MAX_BDL*sizeof(double));
    if ((x == NULL)||(X == NULL))
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };

    f=f0;
    f00=f0;
    bdlsize=0;
    dcopy_(&nn,x0,&inc,x,&inc);  /* x=x0 */

    if (P.printlev > 3)
        printf("        f      cut     ratio  t_left t_now it:\n");

    P.min_outer_it=( P.min_outer_it<3 ? 3 : P.min_outer_it);
    iterates=( forec==0 ? P.min_outer_it : P.outer_it);

    in_it=P.inner_it;

    for (cont=1,iterdone=0;cont; )
        {
	  for (main_loop=1;(main_loop<=iterates);main_loop++,iterdone++)
                {
                    t*=0.95;     /* reduce t */
                    
                    /* remember last function value (for deciding if one should stop 
                       earlier or maybe do some extra rounds) */
                    f00=f0;
                    f0=f;            

		    if ((*pm == 0)&&(valid_cut(n,cut,con,eq,rhs,bg,sp,co))) 
		      /* no more violated triangles found (may happen in forecast) */
		      {
                            
                            if (P.printlev > 2)
                                printf("   %8.3f %8.3f %7.3f %6.0d %6.0d %3.0d\n\n",f,best,(f/best-1)*100,m2,*pm,iterdone+1);
                            
                            /* copy values to output parameters, free and return */
                            
                            *pbest=best;
                            dcopy_(&n,cut+1,&inc,*pcut,&inc);
                            *pf=f;
                            dcopy_(&nn,x0,&inc,*primalX,&inc);
                            
                            free(X);free(x);
			    if (con > 0)
			      {
				free(h0);free(mu);
			      };
			    free(x0); 
                            free(g0);   
                            free(gamma_new);
                            free(cut1);free(cut);
                            
                            return(0);

		      }

                    /* start bundle with given constraints */
                    ret=bdl_mc2(n,*pm,L,T,hash,con,eq,rhs,bg,sp,co,best,&gamma,&mu,in_it,&t,&x,&X,&f,&x0,&bdlsize,P.printlev,P.gap_relax);
                    
                    /* check if the bound is below current best cut+(1-error) (in case of integer weights) 
		       or if the bound is below current best cut+error (in case of rational weights) */
                    if (f < best+P.gap_relax)
                        {

                            if (P.printlev > 2)
                                printf("   %8.3f %8.3f %7.3f %6.0d        %3.0d\n\n",f,best,(f/best-1)*100,*pm,iterdone+1);

                            /* copy values to output parameters, free and return */
                            
                            *pbest=best;
                            dcopy_(&n,cut+1,&inc,*pcut,&inc);
                            *pf=f;
                            dcopy_(&nn,x0,&inc,*primalX,&inc);
                            
                            free(X);free(x);
							if (con > 0)
							  {
								free(h0);free(mu);
							  };
							free(x0); 
										free(g0);   
										free(gamma_new);
										free(cut1);free(cut);
                            
							return(2); /* return and fathom node */
                        };

                    /* x0= .1*x0 + .9*x       .1 from volume paper */
                    alpha=0.1;
                    dscal_(&nn,&alpha,x0,&inc);
                    alpha=0.9;
                    daxpy_(&nn,&alpha,x,&inc,x0,&inc);

                    /* purge inactive constraints */
                    gamma_mean=0.0;
                    for (i=1;i<=*pm;i++)
                        gamma_mean+=gamma[i];
                    gamma_mean/=(*pm);
                    
                    j=1;
                    for (i=1;i<=*pm;i++)
                        {
						if (gamma[i] > 0.01*gamma_mean)
                                {
                                    if (i > j)
                                        {
                                            hash[j]=hash[i];
                                            gamma[j]=gamma[i];
                                            T[ijtok(j,1,MAX_INEQ)]=T[ijtok(i,1,MAX_INEQ)];
                                            T[ijtok(j,2,MAX_INEQ)]=T[ijtok(i,2,MAX_INEQ)];
                                            T[ijtok(j,3,MAX_INEQ)]=T[ijtok(i,3,MAX_INEQ)];
                                            T[ijtok(j,4,MAX_INEQ)]=T[ijtok(i,4,MAX_INEQ)];
                                        };
                                    j++;
                                };
					};
                    *pm=j-1;   /* new number of constraints */
                    m2=*pm;
                    

                    /* generate new constraints */
		    ret=separation(P,x0,n,pm,T,&hash,&gamma_new);             
		    if (*pm > MAX_INEQ)
                        {
                            printf(" number of max. violated inequalities m=%d exceeded. Please adjust MAX_INEQ!\n",*pm);
                            exit(10);
                        };
                    
                    if (P.printlev > 3)
                        printf("   %8.3f %8.3f %7.3f %6.0d %6.0d %3.0d\n",f,best,(f/best-1)*100,m2,*pm,iterdone+1);

                    /* gamma = [gamma; gamma_new*t*.02]; */
                    m3=*pm-m2;
                    if (m3 < 0)
                        {
                            printf(" ERROR: after separation you 'lost' some inequalities...\n");
                            exit(10);
                        };

                    alpha=0.02*t;
                    dscal_(&m3,&alpha,gamma_new+1,&inc);
                    dcopy_(&m3,gamma_new+1,&inc,gamma+m2+1,&inc);

		    if ((*pm == 0)&&(con == 0))  /* no more violated triangles found (may happen in forecast) */
		      {
                            if (P.printlev == 2)
                                printf("   %8.3f %8.3f %7.3f %6.0d %6.0d %3.0d\n\n",f,best,(f/best-1)*100,m2,*pm,iterdone+1);
                            
                            /* copy values to output parameters, free and return */
                            
                            *pbest=best;
                            dcopy_(&n,cut+1,&inc,*pcut,&inc);
                            *pf=f;
                            dcopy_(&nn,x0,&inc,*primalX,&inc);
                            
                            free(X);free(x);
			    if (con > 0)
			      {
				free(h0);free(mu);
			      };
			    free(x0); 
                            free(g0);   
                            free(gamma_new);
                            free(cut1);free(cut);
                            
                            return(0);

		      }


		    in_it++;  // increase in_it by one
                    if (in_it>P.max_inner_it) // but do never more than max_inner_it iterations
		      in_it=P.max_inner_it;

                    if (P.printlev > 5)
                        {
                            printf(" end of iteration in bdl_main:\n");
                            printf(" gamma:\n");
                            for (i=1;i<=(*pm<20 ? *pm : 20);i++)
                                printf(" %8.4f\n ",gamma[i]);
                            printf("\n");
                        };
                };

            switch (cont)
                {
                case 1:     /* after a minimum of outer iterations, check if it's worth doing 
                             * the other outer iterations as well...  */
                    if (forec)        /* if we're doing forecast on branching decisions, 
                                       * just do the demanded iterations and stop. */
                            cont=0;
                    else
                        {
                            iterates=P.outer_it-iterdone;
                            if (f-(iterates+1)*(f0-f) > best+P.gap_relax)
                                cont=0;
                            else
                                cont++;
                        };
                    break;
                case 2:   /* do some extra rounds? */
                    iterates=P.outer_it+P.extra_it-iterdone;
                    if (f-(iterates+1)*(f0-f) > best+P.gap_relax)
                        cont=0;
                    else
                        cont++;
                    break;
                default:  /* have already done enough iterations */
                    cont=0;
                    break;
                };

        }; /* end outer loop */
            
    if (!forec)
        {

            /* 
             * a final computation of the bounds
             */
            in_it=( n>500 ? in_it : in_it+3 );
            t*=1.5;
            
            ret=bdl_mc2(n,*pm,L,T,hash,con,eq,rhs,bg,sp,co,best,&gamma,&mu,in_it,&t,&x,&X,&f,&x0,&bdlsize,P.printlev,P.gap_relax);
            
            /* x0= .1*x0 + .9*x       .1 from volume paper */
            alpha=0.1;
            dscal_(&nn,&alpha,x0,&inc);
            alpha=0.9;
            daxpy_(&nn,&alpha,x,&inc,x0,&inc);
            
            /* 
             * final purge 
             */
            dmax_val(*pm,gamma,&arggmax,&gamma_max); 
            
            j=1;
            for (i=1;i<=*pm;i++)
                {
                    if (gamma[i] > 0.0001*gamma_max)
                        {
                            if (i > j)
                                {
                                    hash[j]=hash[i];
                                    gamma[j]=gamma[i];
                                    T[ijtok(j,1,MAX_INEQ)]=T[ijtok(i,1,MAX_INEQ)];
                                    T[ijtok(j,2,MAX_INEQ)]=T[ijtok(i,2,MAX_INEQ)];
                                    T[ijtok(j,3,MAX_INEQ)]=T[ijtok(i,3,MAX_INEQ)];
                                    T[ijtok(j,4,MAX_INEQ)]=T[ijtok(i,4,MAX_INEQ)];
                                };
                            j++;
                        };
                };
            *pm=j-1;   /* new number of constraints */
            m2=*pm;    /* needed only for printing some output */
            
            if (P.printlev > 4)
                printf("\n ...calling rounding heuristic again:\n");

            dcopy_(&n,cut+1,&inc,cut1+1,&inc);
            if (heuristic(L,x0,n,&cut1,&best1,con,eq,rhs,bg,sp,co)&&(best1 > best))
                {
                    best=best1;
                    dcopy_(&n,cut1+1,&inc,cut+1,&inc);
                };
    
            if (P.printlev > 2)
                {
                    printf(" ...info from last call to mc2:\n");
                    printf("   %8.3f %8.3f %7.3f %6.0d        %3.0d\n",f,best,(f/best-1)*100,*pm,iterdone+1);
                    printf("\n");
                };
        };

    /* copy values to output parameters, free and return */

    *pbest=best;
    dcopy_(&n,cut+1,&inc,*pcut,&inc);
    *pf=f;
    dcopy_(&nn,x0,&inc,*primalX,&inc);

    free(X);free(x);
    if (con > 0)
      {
	free(h0);free(mu);
      };
    free(x0); 
    free(g0);   
    free(gamma_new);
    free(cut1);free(cut);
    return(0);
}
