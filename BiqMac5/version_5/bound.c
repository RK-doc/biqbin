/*
 * bounding routine for bab_for_klu.c
 * find bound and best cut, decide on which edge to branch 
 *
 * feb 2005, a. wiegele
 * oct 2005, a. wiegele (added branching rule 7, 8)
 * april 2006, a. wiegele (added branching rule 9, 10)
 * may 2006, a. wiegele (removed branching rules 7, 8)
 * dec 2006, a. wiegele (added branching rule 7, 'opposite' of 2)
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"


void bound(P,n,cost,con,eq_tmp,rhs_tmp,bg_tmp,sp_tmp,co_tmp,bestvalue,lb,ub,cut,left,right,fix)
     /* input paramters */
     T_param P;
     int n;
     double *cost; /* cost matrix, stored lower triangle row-wise */
     int con; /* number of constraints */
     int *eq_tmp; /* for i=0,...,con-1 
		 eq[i]=-1: constraint i is "<="
		 eq[i]=0:  constraint i is "=="
		 eq[i]=1:  constraint i is ">=" */
     double *rhs_tmp; /* for i=0,...,con-1 
		     rhs[i] is the right hand side of constraint i */
     int *bg_tmp; /* for i=0,...,con (!) 
		 bg[i] is the index of the first index of a non-zero in constraint i */
     int *sp_tmp; /* for j=0,...,bg[con]-1
		 sp[i] is the support of the j-th non-zero (lower triangle, 
		 row-wise, including main diagonal */
     double *co_tmp; /* for j=0,...,bg[con]-1 
		    co[i] is the coefficient of the j-th non-zero */
     double bestvalue; /* the value of the best cut known a priori or
                        * already computed in the outside world
                        */
     /* output parameters */
     double *lb,*ub;  // the l.b./u.b. computed by the procedure
     double *cut;   // The node +/-1 vector defining the best cut
     int *left,*right; /* The left/right end node of the branching variable
                        * 0 <= *left <= n - 1, 0 <= *right <= n - 1
                        */
     int *fix;      /* 
                     *  a) If you want the previous behavior, don't do anything
                     *  b) If you want to generate subproblem Join(i,j) only, set
                     *     *bound_fix = 0;
                     *  c) If you want to generate subproblem Separate(i,j) only, set
                     *     *bound_fix = 1.
                     */

{
    double *L;  /* laplace matrix, stored upper triangle colwise */
    double *primalX;  /* output from bdl_main */

    int i,j,ret;
    double *dist,dmy;

    /* to make a forecast of the branching decisions we need the constraints: */
    double *gamma;
    int T[MAX_INEQ*4];
    unsigned int *hash;
    int m;    /* length of gamma, hash, T */

   /* constraints */
    int equations; 
    int *eq;
    double *rhs;
    int *bg;
    int *sp;
    double *co;
    int inc=1;
    int k,l1,l2,supp;

    if (P.printlev > 5)
      {
	printf(" n= %d\n",n);
	// izpis matrike cen
	for (i=1;i<=n;i++)
	  {
	    for (j=1;j<=n;j++)
	      {
		if (i==j)
		  printf(" %d ",0);
		else
		  if (i>j)
		    printf(" %8.4f ",cost[ijtokr(j,i,n)]);
		  else
		    printf(" %8.4f ",cost[ijtokr(i,j,n)]);
	      };
	    printf("\n");
	  };
      };

    /*
     * make laplace matrix and store it as 
     * packed matrix (upper triangle, columnwise) 
     */
    L=(double *)malloc(n*(n+1)/2*sizeof(double));
    primalX=(double *)malloc(n*(n+1)/2*sizeof(double));
    dist=(double *)malloc(n*sizeof(double));
    if ((L == NULL)||(primalX == NULL)||(dist == NULL))
      {
	printf("Storage allocation failed!\n");
	exit(10);
      };
	// ustvari Laplaceovo matriko (shranjeno kot packed matrix)
    for (i=1;i<=n;i++)
      L[ijtokp(i,i,n)]=0.0;

    for (i=1;i<=n;i++)
      {
		for (j=i+1;j<=n;j++)
		  {
			dmy=cost[ijtokr(i,j,n)]/4.0;
			L[ijtokp(i,j,n)]=-dmy;
			L[ijtokp(i,i,n)]+=dmy;
			L[ijtokp(j,j,n)]+=dmy;
		  };
      };

    if (P.printlev > 5)
      {
		if (con > 0)
		  {
			printf("number of equations=%d\n",con);
			for (i=0;i<con;i++)
			  {
			for(j=bg_tmp[i];j<bg_tmp[i+1];j++)
			  printf(" %lf*x[%d] ",co_tmp[j],sp_tmp[j]);
			if (eq_tmp[i]<0)
			  printf(" <= ");
			else
			  if (eq_tmp[i]==0)
				printf(" == ");
			  else
				printf(" >= ");
			printf("%lf\n",rhs_tmp[i]);
			  };
		  };
	  }; 

    
    for (i=0,equations=0,supp=0;i<con;i++)
      if (eq_tmp[i] == 0)
		{
		  equations++;
		  supp+=bg_tmp[i+1]-bg_tmp[i];
		};

	if (equations > 0)
		{
		eq=(int *)malloc((con+equations)*sizeof(int));
		rhs=(double *)malloc((con+equations)*sizeof(double));
		bg=(int *)malloc((con+equations+1)*sizeof(int));
		sp=(int *)malloc((bg_tmp[con]+supp)*sizeof(int));
		co=(double *)malloc((bg_tmp[con]+supp)*sizeof(double));

		scopy_(&con,eq_tmp,&inc,eq,&inc);
		dcopy_(&con,rhs_tmp,&inc,rhs,&inc);
		i=con+1;
		scopy_(&i,bg_tmp,&inc,bg,&inc);
			i=bg[con];
		scopy_(&i,sp_tmp,&inc,sp,&inc);
		dcopy_(&i,co_tmp,&inc,co,&inc);

		for (i=0,j=0;i<con;i++)
			if (eq_tmp[i] == 0)
			{
				eq[i]=1;
				eq[con+j]=-1;
				rhs[con+j]=rhs_tmp[i];
				k=bg_tmp[i+1]-bg_tmp[i]; /* number coeff in constr i */
				bg[con+j+1]=bg[con+j]+k;
				l1=bg[i];
				l2=bg[con+j];
				scopy_(&k,sp+l1,&inc,sp+l2,&inc);
				dcopy_(&k,co+l1,&inc,co+l2,&inc);
				j++;
			};
		con+=equations;

		if (P.printlev > 5)
			{
			for (i=0;i<con;i++)
				{
			for(j=bg[i];j<bg[i+1];j++)
				printf(" %lf*x[%d] ",co[j],sp[j]);
			if (eq[i]<0)
				printf(" <= ");
			else
				if (eq[i]==0)
				printf(" == ");
				else
				printf(" >= ");
			printf("%lf\n",rhs[i]);
				};
			};
		}
	else
		{
		eq=eq_tmp;
		rhs=rhs_tmp;
		bg=bg_tmp;
		sp=sp_tmp;
		co=co_tmp;
		};



    /* allocate mem for gamma and hash */
    gamma=(double *)malloc((MAX_INEQ+1)*sizeof(double));
    hash=(unsigned int *)malloc((MAX_INEQ+1)*sizeof(unsigned int));
    if ((hash == NULL)||(gamma == NULL))
      {
		printf("Storage allocation failed!\n");
		exit(10);
      };

    /* call bdl_main */
    *lb=bestvalue;
    m=0;  /* and therefore T, hash, gamma only output. if m>0, then T,hash,gamma input as well */

    ret=bdl_main(P,L,n,0,con,eq,rhs,bg,sp,co,lb,ub,&cut,&primalX,&m,T,hash,gamma);
        
    /* value of ret:  1...initial bound+error was below best+1 (node will be fathomed)
     *                2...bound+error below best+1 (node will be fathomed)
     *                0...normal exit
     */
 
	// izpis ustvarjene "matrike" primalX
    if (P.printlev > 5)
      {
	printf("\n\n primal X: \n");
	for (j=1;j<=(n>20 ? 20 : n);j++)
	  {
	    for (i=1;i<=(n>20 ? 20 : n);i++)
	      {
			if (i>j)
			  printf("%6.4f ",primalX[ijtokp(j,i,n)]);
			else
			  printf("%6.4f ",primalX[ijtokp(i,j,n)]);
	      }
			printf("\n");
		};  
		printf("\n");
      };

    /* 
     * determine branching variable if   
     *   *ub >= *lb + (1-error)   (integers)
     *   *ub >= *lb + error       (reals)
     * (otherwise the node will be fathomed)
     */ 

    if (1)//(*ub >= *lb+P.gap_relax)
      {

	switch (P.branch_rule)
	  {
	  case 2:
	    {
	      for (i=1;i<=n;i++)
		{
		  dist[i-1]=0.0;
		  for (j=1;j<=n;j++)
		    {
		      if (i<j)
			dmy=1-Abs(primalX[ijtokp(i,j,n)]);
		      else
			dmy=1-Abs(primalX[ijtokp(j,i,n)]);
		      dmy*=dmy;
		      dist[i-1]+=dmy;
		    };
		};
	      
	      if (P.printlev > 5)
		{
		  printf(" distances in primal X\n");
		  for (i=0;i<( n > 20 ? 20 : n);i++)
		    printf(" %8.4f\n",dist[i]);
		};
	      
	      min_two(n,dist,left,right);
	      
	      break;
	    }
	
	  case 3:
	    {
	      *left=1;
	      *right=2;
	      dmy=Abs(primalX[1]); // = primalX[ijtokp(1,2,n)];
	      for (j=3;j<=n;j++)
			for (i=1;i<j;i++)
			  {
				if (Abs(primalX[ijtokp(i,j,n)])<dmy)
				  {
					*left=i;
					*right=j;
					dmy=Abs(primalX[ijtokp(i,j,n)]);
				  };
			};
	      if (P.printlev > 5)
		printf(" R3: the value of X on the edge to branch is              %8.4e\n",dmy);
	      (*left)--;
	      (*right)--;   // node-index from 0..n-1 !!
	      break;	    
	    }
	
	  case 6:
	    {
	      // find 'first' edge
	      i=1;j=2;
	      while (Abs(L[ijtokp(i,j,n)])<EPSILON)
		{ 
		  i++;
		  if (i==j)
		    {
		      j++;
		      i=1;
		    }
		  if ((i==n)&&(j==n))
		    {
		      printf("Graph contains no edges?!\n");
		      exit(10);
		    };
		};
	      
	      *left=i;
	      *right=j;
	      dmy=Abs(primalX[ijtokp(i,j,n)]); 
	      
	      for (j=*right;j<=n;j++)
		for (i=*left+1;i<j;i++)
		  {
		    if ((Abs(primalX[ijtokp(i,j,n)])<dmy)&&(Abs(L[ijtokp(i,j,n)])>0.0))
		      {
			*left=i;
			*right=j;
			dmy=Abs(primalX[ijtokp(i,j,n)]);
		      };
		  };
	    
	      if (P.printlev > 5)
		printf(" R6: the value of X on the edge to branch is              %8.4e\n",dmy);
	      (*left)--;
	      (*right)--;   // node-index from 0..n-1 !!
	      break;
	    }
	
	  case 7:
	    { /* new 19/12/06: 'opposite' of R2 */
	      for (i=1;i<=n;i++)
		{
		  dist[i-1]=0.0;
		  for (j=1;j<=n;j++)
		    {
		      if (i<j)
			dmy=primalX[ijtokp(i,j,n)];
		      else
			dmy=primalX[ijtokp(j,i,n)];
		      dmy*=dmy;
		      dist[i-1]+=dmy;
		    };
		};
	      
	      if (P.printlev > 1)
		{
		  printf(" distances in primal X\n");
		  for (i=0;i<( n > 20 ? 20 : n);i++)
		    printf(" %8.4f\n",dist[i]);
		};
	      
	      min_two(n,dist,left,right);

	      if (*left<*right)
		printf( "branching X[i,j]=%lf, dist[i]=%lf, dist[j]=%lf\n",primalX[ijtokp((*left+1),(*right+1),n)],dist[*left],dist[*right]);
	      else
		printf( "branching X[i,j]=%lf, dist[i]=%lf, dist[j]=%lf\n",primalX[ijtokp((*left+1),(*right+1),n)],dist[*left],dist[*right]);
	      
	      break;
	    }
	
	  case 9:
	    {
	      int lt,rt;

	      /* choose branching edge according to rule 2 */
	      for (i=1;i<=n;i++)
		{
		  dist[i-1]=0.0;
		  for (j=1;j<=n;j++)
		    {
		      if (i<j)
			dmy=1-Abs(primalX[ijtokp(i,j,n)]);
		      else
			dmy=1-Abs(primalX[ijtokp(j,i,n)]);
		      dmy*=dmy;
		      dist[i-1]+=dmy;
		    };
		};
	      
	      min_two(n,dist,left,right);

	      lt=(*left < *right ? *left+1 : *right+1);
	      rt=(*left < *right ? *right+1 : *left+1);	      
	      if (Abs(primalX[ijtokp(lt,rt,n)])<0.6)
		{
		  /* this is not a "good" branching decision. better choose according to rule 3 */

		  *left=1;
		  *right=2;
		  dmy=Abs(primalX[1]); // = primalX[ijtokp(1,2,n)];
		  for (j=3;j<=n;j++)
		    for (i=1;i<j;i++)
		      {
			if (Abs(primalX[ijtokp(i,j,n)])<dmy)
			  {
			    *left=i;
			    *right=j;
			    dmy=Abs(primalX[ijtokp(i,j,n)]);
			  };
		      };

		  (*left)--;
		  (*right)--;   // node-index from 0..n-1 !!
		};

	      break;
	    }
	  case 10:
	    {
	      int lt,rt;

	      /* check, if there is an obvious branching decision (rule 2) */
	      for (i=1;i<=n;i++)
		{
		  dist[i-1]=0.0;
		  for (j=1;j<=n;j++)
		    {
		      if (i<j)
			dmy=1-Abs(primalX[ijtokp(i,j,n)]);
		      else
			dmy=1-Abs(primalX[ijtokp(j,i,n)]);
		      dmy*=dmy;
		      dist[i-1]+=dmy;
		    };
		};
	      
	      min_two(n,dist,left,right);

	      lt=(*left < *right ? *left+1 : *right+1);
	      rt=(*left < *right ? *right+1 : *left+1);	      
	      if (Abs(primalX[ijtokp(lt,rt,n)])<0.8)
		{
		  /* there seem to be no obvious branching decision, so let's choose upon rule R3 */

		  *left=1;
		  *right=2;
		  dmy=Abs(primalX[1]); // = primalX[ijtokp(1,2,n)];
		  for (j=3;j<=n;j++)
		    for (i=1;i<j;i++)
		      {
			if (Abs(primalX[ijtokp(i,j,n)])<dmy)
			  {
			    *left=i;
			    *right=j;
			    dmy=Abs(primalX[ijtokp(i,j,n)]);
			  };
		      };

		  (*left)--;
		  (*right)--;   // node-index from 0..n-1 !!
		};


	      break;
	    }
	  };
      }
    else
      {   /* node will be fathomed and no branching decision needed. but 
	     let's give some numbers to *left and *right. */
	*left=0;*right=0; 
      };
    
    if (P.printlev > 5)
      printf(" left node = %d, right node = %d\n",*left,*right);
    
    /* free and return */
    free(hash);
    free(gamma);
    
    free(dist);
    free(primalX);
    free(L);
    if (equations > 0)
      {
	free(eq);free(rhs);free(bg);free(sp);free(co);
      };
}

void min_two(k,x,argmin1,argmin2)     
     int k;  // length of x
     double *x;
     int *argmin1;
     int *argmin2;
{ 
    int i;
    double min1,min2;

    if (k<2)
        {
            printf(" primal Matrix has dimension < 2! \n");
            exit(5);
        };

    *argmin1=(x[0] < x[1] ? 0 : 1 );
    min1=x[*argmin1];                // the smaller of the first 2 elements
    *argmin2=(x[0] < x[1] ? 1 : 0 );
    min2=x[*argmin2];                // the bigger of the first 2 elements

    for (i=2;i<k;i++)
        { 
            if (x[i]<min1) 
                {
                    *argmin2=*argmin1;
                    *argmin1=i;
                    min2=min1;
                    min1=x[i];
                }
            else
                {
                    if (x[i]<min2)
                        {
                            *argmin2=i;
                            min2=x[i];
                        };
                };

        };
}
