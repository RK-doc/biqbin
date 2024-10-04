/* 
 * solver for 
 * (P) min <lam,c> + 1/2 <lam,Q lam>  (D) max <b,y> - 1/2 <lam,Q lam>
 *     s.t. <e,lam> = 1                   s.t. z = c + Q lam - y e
 *          lam >= 0                           z >= 0, y free
 *
 * input:  Q,c,k
 * output: lam,z,y
 *
 * translation from solve_lambda.m (oct 06,2000, i.fischer & g.gruber)
 * nov 2004, a.wiegele
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"
#define MAX_ITER 30


int solvelambda(Q,c,k,plam,silent)
    /* input data */
    double *Q;
    double *c;
    int k;
    /* lam (the pointer on lam resp.) 
       is the output of the program. */
    double **plam;    /* primal solution */
    /* printlevel */
    int silent;
{
  /* vars for call to blas- and lapack-routines */
  int inc=1;
  double alpha,beta; 
  char up='U'; 
  int info,*ipiv;
  double *ap;
  
  /* and some other var's */
  double *z;        /* y,z... dual solution */
  double y;
  int i,j,done;
  double *tmp,*Qlam,lamQlam;
  int argmintmp;
  double mintmp;
  double mu=0.0;
  double A_lam,*AT_y;
  double P_cost,D_cost,pd_gap;
  
  double *M,*rhs;
  double *dw,*dlam,dy,*dz; /* dw is sol of dspsv, dlam is pointer to dw, dy is dw[k+1] */
  int kinc,nn;             /* kinc=k+1; nn=(k+1)*(k+2)/2; */
  double alpha_p,alpha_d;  /* step length */
  double *p,*pz,*pM,*pQ;
  
  /**********************************
   *                                *
   * initialization:                *
   *                                *
   * starting triplet (infeasible). *
   * primal/dual cost, gap          *
   **********************************/
  
  /* allocate space for lam and z */
  /* *plam=(double *)malloc((k+1)*sizeof(double)); ... done in bdl_mc2!!! */
  z=(double *)malloc((k+1)*sizeof(double));
  if (z == NULL)
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };
  
  /* lam = (1/k)*e */
  for (i=1,p=(*plam)+1;i<=k;i++,p++)
    (*p)=1.0/k; // (*plam)[i]=1.0/k;

  /* allocate space for tmp and Qlam */
  Qlam = (double *)calloc((k+1),sizeof(double));
  tmp = (double *)malloc((k+1)*sizeof(double));
  if ((Qlam == NULL)||(tmp == NULL))
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };

  /* tmp = Q*lam + c */
  alpha=1.0; beta=0.0;
  dspmv_(&up,&k,&alpha,Q,*plam+1,&inc,&beta,Qlam+1,&inc); /* Qlam = Q*lam; */
  dcopy_(&k,c+1,&inc,tmp+1,&inc);                         /* tmp=c */
  daxpy_(&k,&alpha,Qlam+1,&inc,tmp+1,&inc);               /* tmp=Qlam+tmp */
  
  /* mintmp = min(tmp) */
  dmin_val(k,tmp,&argmintmp,&mintmp); 

  if (mintmp > 1)
    { 
      /* y=0; z=tmp; */
      y=0;
      dcopy_(&k,tmp+1,&inc,z+1,&inc);
    }
  else
    {   /* y=mintmp-1; z=tmp-y; */
      y=mintmp-1;
      dcopy_(&k,tmp+1,&inc,z+1,&inc);
      for (i=1,p=z+1;i<=k;i++,p++)
	(*p)-=y; // z[i]-=y;
    };
  
  /* mu=(z'*lam)/k; mu=mu/2; */
  mu=ddot_(&k,z+1,&inc,*plam+1,&inc);    
  mu/=(2*k);
  
  /* A_lam=e'*lam=1; AT_y=y*e; */
  A_lam=1.0;    
  AT_y=(double *)malloc((k+1)*sizeof(double));
  if (AT_y == NULL)
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };
  for (i=1,p=AT_y+1;i<=k;i++,p++)
    (*p)=y; // AT_y[i]=y;  
  
  /* D_cost = y - 1/2*lam'*Q*lam; P_cost = lam'*c + 1/2*lam'*Q*lam; */
  lamQlam=ddot_(&k,*plam+1,&inc,Qlam+1,&inc);
  D_cost=y-0.5*lamQlam;
  P_cost=ddot_(&k,*plam+1,&inc,c+1,&inc)+0.5*lamQlam;
  pd_gap=P_cost-D_cost;
  
  /*****************************
   *                           *
   * main loop                 *
   *                           *
   *****************************/

  /* allocate space for M, rhs and dw, dz */
  M=(double *)malloc((k+1)*(k+2)/2*sizeof(double));
  rhs=(double *)malloc((k+2)*sizeof(double));
  dw=(double *)malloc((k+2)*sizeof(double));
  dz=(double *)malloc((k+1)*sizeof(double));
  if ((M == NULL)||(rhs == NULL)||(dw == NULL)||(dz == NULL))
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };
  
  /* allocate space for int array ipiv which is used in 
   * the call to dspsv. also, allocate space for ap. we will
   * copy M to ap, because in dspsv on output the matrix is changed, 
   * but we need M in each iteration...
   */
  ipiv=(int *)malloc((k+2)*sizeof(int));
  ap=(double *)malloc((k+1)*(k+2)/2*sizeof(double));
  if ((ipiv == NULL)||(ap == NULL))
    {
      printf("Storage allocation failed!\n");
      exit(10);
    };
  
  /* M=[-Q e; e' 0] 
   * this is a symmetric matrix, so we have to add only 
   * vector e (i.e. last col of M) to the upper triangle of Q
   *
   * in each iteration only the main diag of M is changed.
   */
  
  nn=k*(k+1)/2;
  dcopy_(&nn,Q,&inc,M,&inc);  /* M=Q; */
  alpha=-1.0;
  dscal_(&nn,&alpha,M,&inc);  /* M=-M; */
  
  for (j=k*(k+1)/2;j<(k+1)*(k+2)/2-1;j++)
    M[j]=1.0;
  M[(k+1)*(k+2)/2-1]=0.0;

  if (silent == 0)
    {
      printf("\n iter alpha_p alpha_d  log(mu)  P_cost   D_cost\n");
      printf("                       %5.4f  %5.4f  %5.4f\n",log10(mu),P_cost,D_cost);
    };

  for (i=1,done=1,kinc=k+1;i<=MAX_ITER && done>0;i++)
    {
      /* M(i,i)=-Q(i,i)-z(i)/lam(i);  */
      for (j=1,pM=M,pQ=Q,pz=z+1,p=(*plam)+1;j<=k;j++,p++,pz++)
	{
	  /* M[ijtokp(j,j,k)]=-Q[ijtokp(j,j,k)]-z[j]/(*plam)[j]; */
	  (*pM)=-(*pQ)-(*pz)/(*p);
	  pM+=(j+1);
	  pQ+=(j+1);
	};

      /* rhs=[ c-AT_y+Q*lam-mu./lam; 1-A_lam ];  */
      for (j=1,p=(*plam)+1,pz=rhs+1;j<=k;j++,p++,pz++)
	(*pz)=-mu/(*p); //rhs[j]=-mu/(*plam)[j];
      alpha=1.0;
      daxpy_(&k,&alpha,tmp+1,&inc,rhs+1,&inc);
      alpha=-1.0;
      daxpy_(&k,&alpha,AT_y+1,&inc,rhs+1,&inc);
      rhs[k+1]=1.0-A_lam;
      
      /* copy rhs to dw and M to ap */
      dcopy_(&kinc,rhs+1,&inc,dw+1,&inc);
      nn=kinc*(kinc+1)/2;
      dcopy_(&nn,M,&inc,ap,&inc);
      
      /* solve the linear equation */
      dspsv_(&up,&kinc,&inc,ap,ipiv+1,dw+1,&kinc,&info);
      if (info != 0)
	{
	    fprintf(stderr,"solvelambda: problems in solving a system of linear equations! %d\n",info);
	  //fwrite(&k,sizeof(int),1,stderr);
	  //fwrite(Q,sizeof(double),k*(k+1)/2,stderr);
	  //fwrite(c,sizeof(double),k+1,stderr); /* note that c starts with c[1] */
	  exit(10);
	};
      dlam=dw;        /* dlam=dw(1:k); */
      dy=dw[kinc];    /* dy=dw(k+1);   */
      
      /*   dz=1./lam .* (mu*e-lam.*z - z.*dlam);  */
      for (j=1;j<=k;j++)
	dz[j]=1/(*plam)[j]*(mu-(*plam)[j]*z[j]-z[j]*dlam[j]);
      
      /*
       * find steplengths alpha_p and alpha_d s.t.   
       * lam+alpha_p*dlam >= 0 and z+alpha_d*dz >= 0 
       */
      alpha_p=steplength(k,dlam,*plam);
      alpha_d=steplength(k,dz,z);
      
      if (alpha_p > 0)
	if (0.99/alpha_p < 1)
	  alpha_p=0.99/alpha_p;
	else
	  alpha_p=1;
      else
	alpha_p=1;
      if (alpha_d > 0)
	if (0.99/alpha_d < 1)
	  alpha_d=0.99/alpha_d;
	else
	  alpha_d=1;
      else
	alpha_d=1;
      
      /*
       * update 
       */
      
      daxpy_(&k,&alpha_p,dlam+1,&inc,*plam+1,&inc); /* lam=lam+alpha_p*dlam */
      y=y+alpha_d*dy;                               /* y  =y  +alpha_d*dy   */
      daxpy_(&k,&alpha_d,dz+1,&inc,z+1,&inc);       /* z  =z  +alpha_d*dz   */
      
      A_lam=0.0;
      for (j=1,p=(*plam)+1,pz=AT_y+1;j<=k;j++,p++,pz++)      
	{
	  A_lam+=(*p);                 /* A_lam=e'*lam; */
	  (*pz)=y; //AT_y[j]=y;        /* AT_y=y*e;     */
	}
      
      /* tmp = Q*lam + c */
      alpha=1.0; beta=0.0;
      dspmv_(&up,&k,&alpha,Q,*plam+1,&inc,&beta,Qlam+1,&inc); /* Qlam = Q*lam; */
      dcopy_(&k,c+1,&inc,tmp+1,&inc);                         /* tmp=c */
      daxpy_(&k,&alpha,Qlam+1,&inc,tmp+1,&inc);               /* tmp=Qlam+tmp */
      
      /* mu=(z'*lam)/k; */
      mu=(ddot_(&k,z+1,&inc,*plam+1,&inc))/k;    
      mu*=0.4; /* reduce mu */
      
      if (alpha_p+alpha_d > 1.8)
	mu*=0.2;
      
      /* D_cost = y - 1/2*lam'*Q*lam; P_cost = lam'*c + 1/2*lam'*Q*lam; */
      lamQlam=ddot_(&k,*plam+1,&inc,Qlam+1,&inc);
      D_cost=y-0.5*lamQlam;
      P_cost=ddot_(&k,*plam+1,&inc,c+1,&inc)+0.5*lamQlam;
      pd_gap=P_cost-D_cost;
      
      if ((Abs(pd_gap) < EPSILON6)||(Abs(pd_gap) < Abs(P_cost)*EPSILON6))
	done=0;
      
      if (silent == 0)
	{
	  /* display current iteration */
	  printf("   %2d  %5.4f  %5.4f  %5.4f  %5.4f  %5.4f\n",i,alpha_p,alpha_d,log10(mu),P_cost,D_cost);
	};
    };
  
  if (silent == 0)
    {
      printf("\n    *** solution ***\n");
      for (j=1;j<=k;j++)
	printf(" lam%d = %lf z = %lf \n",j,(*plam)[j],z[j]);
      printf("\n y = %lf\n",y); 
    };
  
  
  
  /* free memory and return   */
  free(ap);free(ipiv); 
  free(dz);free(dw);free(rhs);free(M);
  free(AT_y);free(tmp);free(Qlam);
  free(z);
  return(0);
}

/*
 * compute max(-dx./x);
 */
double steplength(k,dx,x)
     int k;
     double *dx;
     double *x;
{ 
    int i,argmax;
    double *dmy, max;

    dmy=(double *)malloc((k+1)*sizeof(double));
    if (dmy == NULL)
        {
            printf("Storage allocation failed!\n");
            exit(10);
        };
    for (i=1;i<=k;i++)
        dmy[i]=-dx[i]/x[i];

    dmax_val(k,dmy,&argmax,&max);
    free(dmy);
    return(max);
}

void dmax_val(k,x,argmax,max)
     int k;
     double *x;
     int *argmax;
     double *max;
{ 
    int i;
    *argmax=1;
    *max=x[1];
    for (i=2;i<=k;i++)
        {
            if (x[i]>(*max))
                {
                    *argmax=i;
                    *max=x[i];
                };
        };
}

void dmin_val(k,x,argmin,min)
     int k;
     double *x;
     int *argmin;
     double *min;
{ 
    int i;
    *argmin=1;
    *min=x[1];
    for (i=2;i<=k;i++)
        {
            if (x[i]<(*min))
                {
                    *argmin=i;
                    *min=x[i];
                };
        };
}
