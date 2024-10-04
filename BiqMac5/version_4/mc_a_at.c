/* 
 * 
 * this is from mex-files mc_a and mc_at (f. rendl)
 * changes: matrix X is now stored in packed format!
 *
 * dec 2004, a.wiegele
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "declarations.h"

/*
 * mc_a(n,m,T,X,ax);
 * input: n ... number of nodes, m ... number of triangles
 *        m_mem ... number of allocated space for T (i.e. m <= m_mem);
 *        (necessary, because sometimes m_mem rows are allocated, but only m in use)
 *        T ... m x 4 matrix describing triangle support+type
 *        X ... n*(n+1)/2 vector (corresponds to matrix X=(x_ij)
 * output: ax: A(X)
 */

int mc_a(n,m,m_mem,T,X,pax)
     int n,m,m_mem,*T;
     double *X;
     double **pax;
{
  double xij,xik,xjk;
  int i,j,k,r,typ;	
  double *p;

  for (r=1,p=(*pax)+1;r<=m;r++,p++)
    {
      i=T[ijtok(r,1,m_mem)];
      j=T[ijtok(r,2,m_mem)];
      k=T[ijtok(r,3,m_mem)];
      typ=T[ijtok(r,4,m_mem)];
      xij=(i<=j ? X[ijtokp(i,j,n)] : X[ijtokp(j,i,n)]);
      xik=(i<=k ? X[ijtokp(i,k,n)] : X[ijtokp(k,i,n)]);
      xjk=(j<=k ? X[ijtokp(j,k,n)] : X[ijtokp(k,j,n)]);
      
      (*p)=-(X[ijtokp(i,i,n)]+X[ijtokp(j,j,n)]+X[ijtokp(k,k,n)]);  /* note: -(xii+xjj+xkk) = -3 */
      
      switch (typ)
	{
	case 0:
	  (*p)-= 2.0*(  xij +xik +xjk);     
	  break;
	case 1:
	  (*p)-= 2.0*(  xij -xik -xjk);
	  break;
	case 2:
	  (*p)-= 2.0*( -xij +xik -xjk);
	  break;
	case 3:
	  (*p)-= 2.0*( -xij -xik +xjk);
	  break;
	}
    };
  return(0);
}

/*
 *  mc_at(n,m,m_memT,gamma,pAtg);
 *  input: n ... number of nodes, m ... number of triangles
 *         m_mem ... number of allocated space for T (i.e. m <= m_mem);
 *         (necessary, because sometimes m_mem rows are allocated, but only m in use)
 *         T ... m x 4 matrix describing triangle support+type
 *         gamma ... m vector (corresponds to dual variables
 *  output: Atg ... A^t(gamma)  (matrix as n*n vector)
 *  
 */

int mc_at(n,m,m_mem,T,gamma,pAtg)
     int n,m,m_mem,*T;
     double *gamma;
     double **pAtg;
{
  double g;
  int i,j,k,r,typ,nn;	
  double *p;
  
  nn=n*(n+1)/2;
  
  for (r=0,p=(*pAtg);r<nn;r++,p++) (*p)=0.0;
  for (r=1; r<=m; r++)
    {
      i=T[ijtok(r,1,m_mem)];
      j=T[ijtok(r,2,m_mem)];
      k=T[ijtok(r,3,m_mem)];
      typ=T[ijtok(r,4,m_mem)];
      
      /* note: matrix A corrsponding to a triangle looks (for example) like 
       *             -1 -1  1
       *             -1 -1  1
       *              1  1 -1     (...typ 2)
       */
      
      g=gamma[r];
      (*pAtg)[ijtokp(i,i,n)] -= g;
      (*pAtg)[ijtokp(j,j,n)] -= g;
      (*pAtg)[ijtokp(k,k,n)] -= g;
      
      switch (typ)
	{
	case 0:
	  {
	    (*pAtg)[(i<=j ? ijtokp(i,j,n) : ijtokp(j,i,n))] -= g;
	    (*pAtg)[(i<=k ? ijtokp(i,k,n) : ijtokp(k,i,n))] -= g;
	    (*pAtg)[(j<=k ? ijtokp(j,k,n) : ijtokp(k,j,n))] -= g;
	    break;
	  }
	case 1:
	  {
	    (*pAtg)[(i<=j ? ijtokp(i,j,n) : ijtokp(j,i,n))] -= g;
	    (*pAtg)[(i<=k ? ijtokp(i,k,n) : ijtokp(k,i,n))] += g;
	    (*pAtg)[(j<=k ? ijtokp(j,k,n) : ijtokp(k,j,n))] += g;
	    break;
	  }
	case 2:
	  {
	    (*pAtg)[(i<=j ? ijtokp(i,j,n) : ijtokp(j,i,n))] += g;
	    (*pAtg)[(i<=k ? ijtokp(i,k,n) : ijtokp(k,i,n))] -= g;
	    (*pAtg)[(j<=k ? ijtokp(j,k,n) : ijtokp(k,j,n))] += g;
	    break;
	  }
	case 3:
	  {
	    (*pAtg)[(i<=j ? ijtokp(i,j,n) : ijtokp(j,i,n))] += g;
	    (*pAtg)[(i<=k ? ijtokp(i,k,n) : ijtokp(k,i,n))] += g;
	    (*pAtg)[(j<=k ? ijtokp(j,k,n) : ijtokp(k,j,n))] -= g;
	    break;
	  }
	};
    };
  return(0);  
}
