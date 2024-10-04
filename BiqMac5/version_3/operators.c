/* 
  This file contains the operators diag, Diag and <A,B>
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"

void Diag(pX,y,n)
     double **pX;
     double *y;
     int n;
{
  int i,j;
  double *z,*v;
  
  for (i=1,z=(*pX),v=y+1;i<=n;i++,v++)
    {
      for (j=1;j<i;j++,z++)
	*z=0.0;
      *z=*v;
      z++;
    };
}

void diag(X,py,n)
     double *X;
     double **py;
     int n;
{
  int i;
  double *z,*v;
  
  for (i=1,z=X,v=(*py)+1;i<=n;i++,v++,z+=i)
    *v=*z;
}


/* traceproduct of matrices stored in packed format */
double traceprod(A,B,n)
     double *A,*B;
     int n;
{
  int nn,inc=1,i;
  double ret;
  double *za,*zb;

  nn=n*(n+1)/2;
  ret=2*ddot_(&nn,A,&inc,B,&inc);
  for (i=1,za=A,zb=B;i<=n;i++,za+=i,zb+=i)
    ret-=(*za)*(*zb);
  return(ret);
}
