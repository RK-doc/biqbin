/*
function [T, gamma] = trianglesep( y, n, maxtri);
        input: matrix y (n,n) matrix as long rowvector
               maxtri .. triangles requested
                  
        output: [T,gamma] data structure of new triangles
		gamma: violation
       last modified: jan 2018
*/

#include <math.h>
#include <stdio.h>
#include "mex.h"
/*#include "heapsort.h"*/

#define MAX_INEQ 20000

#define COMPARE(x,y) (((x)<(y)) ? -1: ((x)==(y))? 0: 1)


/*
 * function prototypes
 */

int binsearch(long list[], int num, int left, int right);


/*
 * subroutines
 */

void heapify_n(i,n,ind,d)
int i,n,*ind;
double *d;
{
 int k,h,j;

 ind--; i++;
 k=i;
 while(2*k<=n){
     if (2*k==n){
         if (d[ind[k]]>d[ind[2*k]]){
             h=ind[k];
             ind[k]=ind[2*k];
             ind[2*k]=h;
         }
         break;
     }
     if ((d[ind[k]]>d[ind[2*k]])||
         (d[ind[k]]>d[ind[2*k+1]])){
         if (d[ind[2*k]]<d[ind[2*k+1]]) j=2*k;
         else j=2*k+1;
         h=ind[k];
         ind[k]=ind[j];
         ind[j]=h;
         k=j;
     }
     else break;
 }
}


int binsearch(long list[], int num, int left, int right)
{
/* durchsuche list[0]<=list[1]<=...<=list[n-1] nach num 
gib +1 zurueck falls gefunden, sonst gib 0 zurueck*/
int middle;
while (left <= right) {
  middle=floor((left+right)/2);
  switch(COMPARE(list[middle],num)) 
  { case -1: left = middle+1;break;
    case 0 : return 1;right=0;
    case 1 : right=middle-1;
  }
}
return 0;
}

 

void mexFunction( int nlhs, mxArray *plhs[],
		  int nrhs, const mxArray *prhs[])
{
  double *y, *n, *maxtri, *T, *gamma;
  int l, vmax, ftest=0;
	
  /* call with 3 inputs and 3 outputs */
  if (nlhs != 2 || nrhs != 3) {
      mexErrMsgTxt("input-output arguments inconsistent.");
  }

  y = mxGetPr(prhs[0]);
  n = mxGetPr(prhs[1]);
  maxtri = mxGetPr(prhs[2]);
  
  vmax = maxtri[0]; l = n[0];
  if (vmax>MAX_INEQ) {vmax = MAX_INEQ;}
  if (mxGetM(prhs[0]) != l*l) {
      mexErrMsgTxt("y and n inconsistent.");
  }

  {
    double threshold,x, yij, yik, yjk, yi,yj,yk,yy,violl[MAX_INEQ];
    int ind[MAX_INEQ], i, j, k, ii, cntr;
    int T_tmp[MAX_INEQ*4], b_tmp[MAX_INEQ];
    int flag,i1,typ;
    
    for (i=0; i<vmax; i++) {
	violl[i] = -1.0;
	ind[i] = i;
    }

    cntr = 0;
    threshold = .001;   /* violation should be larger than threshold */

    for (i=1; i<=l-2; i++) 
    { for (j=i+1; j<=l-1; j++)
      { for (k=j+1; k<=l; k++)
	{ yij = y[(i-1)*l+j-1];
	  yik = y[(i-1)*l+k-1];
	  yjk = y[(j-1)*l+k-1];

/* Typ 0 */
      x = yij + yik + yjk -2.0;
	  if (x>threshold) 
	   { if (x>violl[ind[0]])
	     { ii = ind[0];
	       violl[ii] = x;
	       T_tmp[ii] = i; T_tmp[vmax+ii] = j; 
	       T_tmp[2*vmax+ii] = k; T_tmp[3*vmax+ii] = 0;
	       heapify_n(0,vmax,ind,violl);
	       cntr++;
	     }
	   }
/* Typ 1 */
	  x = yij - yik - yjk;
	  typ=1;

	  if (x>threshold)
	   { if (x>violl[ind[0]])
	     { ii = ind[0];
	       violl[ii] = x;
	       T_tmp[ii] = i; T_tmp[vmax+ii] = j; 
	       T_tmp[2*vmax+ii] = k; T_tmp[3*vmax+ii] = typ;
	       heapify_n(0,vmax,ind,violl);
	       cntr++;
	     }
	   }
/* Typ 2 */
	  x =  - yij + yik - yjk;
	  typ=2;
	 
	  if (x>threshold)
	   { if (x>violl[ind[0]])
	     { ii = ind[0];
	       violl[ii] = x;
	       T_tmp[ii] = i; T_tmp[vmax+ii] = j; 
	       T_tmp[2*vmax+ii] = k; T_tmp[3*vmax+ii] = typ;
	       heapify_n(0,vmax,ind,violl);
	       cntr++;
	     }
	   }
/* Typ 3 */
	  x =  - yij - yik + yjk;
	  typ=3;
	  
	  if (x>threshold)
	   { if (x>violl[ind[0]])
	     { ii = ind[0];
	       violl[ii] = x;
	       T_tmp[ii] = i; T_tmp[vmax+ii] = j; 
	       T_tmp[2*vmax+ii] = k; T_tmp[3*vmax+ii] = typ;
	       heapify_n(0,vmax,ind,violl);
	       cntr++;
	     }
	   }

	}
      }
    }

  ende:
    if (cntr > vmax) cntr= vmax;
    plhs[0] = mxCreateDoubleMatrix(cntr*4, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(cntr, 1, mxREAL);

    T = mxGetPr(plhs[0]);
    gamma = mxGetPr(plhs[1]);
    
    i = 0;
    for (k=0; k<vmax; k++) {
	j = ind[k];
	if (violl[j]>-1.) {
	    gamma[i] = violl[j];
	    T[i] = (double) T_tmp[j];
	    T[cntr+i] = (double) T_tmp[vmax+j];
	    T[2*cntr+i] = (double) T_tmp[vmax*2+j];
	    T[3*cntr+i] = (int) T_tmp[vmax*3+j];
	    i++;}
    }


  }
}
