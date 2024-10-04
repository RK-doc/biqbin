/*
 * december 2004, a. wiegele
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "declarations.h"

int compare_long(const void *a, const void *b)
{
    return(((*(unsigned int *)a)<(*(unsigned int *)b)) ? -1: ((*(unsigned int *)a)==(*(unsigned int *)b))? 0: 1);
}


int separation(P,X,n,pm,pT,phash,pgamma_new)
    /* input */
    T_param P;
    double *X;
    int n;
    int *pm;

    /* in/output */
    int pT[MAX_INEQ*4];
    unsigned int **phash;

    /* output */
    double **pgamma_new;

{
    /* var's for call to blas/lapack */
    int inc=1;
    double alpha;

    /* some other var's */
    int i,j,ret;
    int new_ineq,nn;
    double *X1;    
    unsigned int *h_tmp;
    double *bn;
    double *pp;

    /* output of tri_sep2 */
    int *Tn;
    unsigned int *hashn;
    double *g_new;

    int printlevel;

    nn=n*(n+1)/2;
    new_ineq=n*5;
    new_ineq=( new_ineq< P.heap_size ? P.heap_size : new_ineq );

    /*
     * find new violated constraints in X
     */
    /* transform X to 0-1 model */
    X1=(double *)malloc(nn*sizeof(double));
    if (X1 == NULL)
      {
          printf("Storage allocation failed!\n");
          exit(10);
      };
    for (i=0,pp=X1;i<nn;i++,pp++)
        (*pp)=0.5; // X1[i]=0.5;
    alpha=-0.5;              
    daxpy_(&nn,&alpha,X,&inc,X1,&inc); /* X1=(ones(n)-X)/2 */

    /* copy hash to h_tmp and sort h_tmp */
    h_tmp=(unsigned int *)malloc(((*pm)+1)*sizeof(unsigned int));
    if (h_tmp == NULL)
        {
          printf("Storage allocation failed!\n");
          exit(10);
        };    
    for (i=1;i<=*pm;i++)
        h_tmp[i]=(*phash)[i];
    qsort(h_tmp+1,*pm,sizeof(unsigned int),compare_long);

	// increase print level to output
    printlevel=-10;
    if (printlevel > 1)
        {
            printf(" separation: befor calling tri_sep2 new_ineq=%d m=%d\n",new_ineq,*pm);
            for (i=1;i<=*pm;i++)
                {
                    for (j=1;j<=4;j++)
                        printf(" %d ",pT[ijtok(i,j,MAX_INEQ)]);//(*pT)[ijtok(i,j,MAX_INEQ)]);
                    printf("   %d \n",h_tmp[i]);
                };
        };

    tri_sep2(P,X1,n,&new_ineq,h_tmp,*pm,&Tn,&hashn,&g_new); /* mem allocated in tri_sep2 */
    /* in new_ineq is the length of hashn stored */

    printlevel=-10;
    if (printlevel > 1)
        {
            printf(" separation: after calling tri_sep2 new_ineq=%d m=%d\n",new_ineq,*pm);
            for (i=1;i<=new_ineq;i++)
                {
		                     for (j=1;j<=4;j++)
		      printf(" %d ",Tn[ijtok(i,j,new_ineq)]);
                    printf("   %lf \n",g_new[i]);
                };
        };


    /* if (new_ineq==0)
       printf(" separation: nothing new violated.\n"); */
    if (new_ineq>0) //else /* new_ineq > 0 */
        {
            /* right hand side is -1 */
            bn=(double *)malloc((new_ineq+1)*sizeof(double));
            if (bn == NULL)
                {
                    printf("Storage allocation failed!\n");
                    exit(10);
                };

            for (i=1;i<=new_ineq;i++)
                bn[i]=-1.0;    

            ret=mc_a(n,new_ineq,new_ineq,Tn,X,pgamma_new);

            alpha=-1.0;    /* gamma_new=ax-bn */
            daxpy_(&new_ineq,&alpha,bn+1,&inc,*pgamma_new+1,&inc);

            /*
             * merge new with old constraints
             */
	    //printf("in separation:\n");

            for (i=1;i<=new_ineq;i++)
                {
		  //printf(" %lf\n",(*pgamma_new)[i]);
                    (*phash)[i+(*pm)]=hashn[i];
                    pT[ijtok(i+(*pm),1,MAX_INEQ)]=Tn[ijtok(i,1,new_ineq)];
                    pT[ijtok(i+(*pm),2,MAX_INEQ)]=Tn[ijtok(i,2,new_ineq)];
                    pT[ijtok(i+(*pm),3,MAX_INEQ)]=Tn[ijtok(i,3,new_ineq)];
                    pT[ijtok(i+(*pm),4,MAX_INEQ)]=Tn[ijtok(i,4,new_ineq)];
                };
            (*pm)+=new_ineq;

            free(bn);
        };

    /* free and return (also free memory, allocated in tri_sep2!) */
    free(g_new);free(hashn);free(Tn);
    free(h_tmp);free(X1);
    return(0);
}

