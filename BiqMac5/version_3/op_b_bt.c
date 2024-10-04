/* 
 * 
 * may 2006, a. wiegele
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"

/*
 * op_b(n,m,T,X,ax);
 * input: n ... number of nodes
 *        con,eq,bg,sp,co ... defining the operator B
 *        X ... n*(n+1)/2 vector (corresponds to matrix X=(x_ij)
 * output: bx: B(X)
 */

int op_b(n,con,eq,bg,sp,co,X,pbx)
     int n,con;
     int *eq;
     int *bg;
     int *sp;
     double *co;
     double *X;
     double **pbx;
{
  int i,j;	
  double *p,*pco;
  int *s;

  for (i=0,p=(*pbx)+1,pco=co,s=sp;i<con;i++,p++)
    {
      (*p)=0.0;
      if (eq[i]<1)
	for (j=0;j<(bg[i+1]-bg[i]);j++,pco++,s++)
	  (*p)+=(*pco)*X[*s];
      else
	for (j=0;j<(bg[i+1]-bg[i]);j++,pco++,s++)
	  (*p)-=(*pco)*X[*s];

    };
  return(0);
}

/*
 *  op_bt(n,con,eq,bg,sp,co,m,nu,pBtn);
 *  input: 
 *         con,eq,bg,sp,co ... defining the operator B
 *         n ... number of nodes
 *         nu ... con vector (corresponds to dual variables)
 *  output: Btn ... B^t(nu)  (matrix as n*(n+1)/2 vector)
 *  
 */

int op_bt(n,con,eq,bg,sp,co,nu,pBtn)
     int n,con;
     int *eq;
     int *bg;
     int *sp;
     double *co;
     double *nu;
     double **pBtn;
{
  int i,j,nn;	
  double *p;
  int *s;
  
  nn=n*(n+1)/2;
  
  for (i=0,p=(*pBtn);i<nn;i++,p++) (*p)=0.0;
  for (i=0,p=co,s=sp;i<con; i++)
    {
      if (eq[i]<1)
	for (j=0;j<(bg[i+1]-bg[i]);j++,p++,s++)
	  {
	    //	    printf("add 0.5*%lf*%lf\n",nu[i+1],*p);
	    (*pBtn)[(*s)] += 0.5*(nu[i+1])*(*p);
	  }
      else
	for (j=0;j<(bg[i+1]-bg[i]);j++,p++,s++)
	  {
	    //      printf("subtract 0.5*%lf*%lf\n",nu[i+1],(*p));
	    (*pBtn)[(*s)] -= 0.5*(nu[i+1])*(*p);
	  };

    };
  return(0);  
}

int valid_cut(n,cut,con,eq,rhs,bg,sp,co)  /* NOTE: cut is in cut[1]..cut[n] !! */
    int n;
    double *cut;
    int con;
    int *eq;
    double *rhs;
    int *bg;
    int *sp;
    double *co;
{
  int i,valid,j;
  double lhs;
  
  for(i=0,valid=1;valid && i<con;i++) 
    {
      lhs = 0.0;
      for(j=bg[i];j<bg[i+1];j++)
	lhs+=co[j]*cut[l2i1(sp[j])+1]*cut[l2i2(sp[j])+1];
      if((eq[i] >= 0 && Positive(rhs[i]-lhs))
	 || (eq[i] <= 0 && Negative(rhs[i]-lhs))) 
	valid = 0; 
    };

  return(valid);
};
