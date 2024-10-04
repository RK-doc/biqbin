/*
 * function [T, b, hash, gamma] = tri_sep2(X,n,vmax,h_tmp,m);
 *      input: matrix X (n,n) in sym packed format
 *             vmax ... number of desired violated triangles
 *                (vmax=0) -> just checking if there are violated triangles
 *             h_tmp...sorted old hash-array
 *      output: [T,b,hash] data structure of new triangles
 *	gamma: violation
 *
 * adapted from the mex-file tri_sep2.c (helmberg, rendl)
 * dec 2004, a. wiegele
 *
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "declarations.h"
#include "gb_flip.h"

/*
 * function prototypes
 */
int binsearch(unsigned int list[], unsigned int num, int left, int right);


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
  while(2*k<=n)
    {
      if (2*k==n)
		{
		  if (d[ind[k]]>d[ind[2*k]])
			{
			  h=ind[k];
			  ind[k]=ind[2*k];
			  ind[2*k]=h;
			}
		  break;
		}
      if ((d[ind[k]]>d[ind[2*k]])||(d[ind[k]]>d[ind[2*k+1]]))
		{
		  if (d[ind[2*k]]<d[ind[2*k+1]]) 
			j=2*k;
		  else 
			j=2*k+1;
		  h=ind[k];
		  ind[k]=ind[j];
		  ind[j]=h;
		  k=j;
		}
      else break;
    }
}

void heapsort(n,ind,d)
int n,*ind;
double *d;
{
 int i,h;

/* //build heap  */
// build_heap(n,ind,d);

/* //extract greatest and rebuild heap  */
 for(i=n-1;i>=1;i--){
     h=ind[i];ind[i]=ind[0];ind[0]=h;
     heapify_n(0,i,ind,d);
 }
}


int binsearch(unsigned int list[], unsigned int num, int left, int right)
{
/* durchsuche list[0]<=list[1]<=...<=list[n-1] nach num 
gib +1 zurueck falls gefunden, sonst gib 0 zurueck*/
  unsigned int middle;
  while (left <= right) 
    {
      middle=floor((left+right)/2);
      switch(COMPARE(list[middle],num)) 
	{ 
	case -1: left = middle+1;break;
	case 0 : return 1;right=0;
	case 1 : right=middle-1;
	}
    }
  return 0;
}

 
int tri_sep2(P,X,n,pvmax,h_tmp,m,pT,phash,pgamma)
    /* input */
    T_param P;
    double *X; /* sym packed matrix */
    int n;     /* dim of X */
    int *pvmax;  /* number of desired violated triangles, on output length hash array */
    unsigned int *h_tmp; /* hash array (sorted) */  
    int m;     /* length of hash array */
    /* output */
    int **pT;
    unsigned int **phash;
    double **pgamma;
{
    /* and some other var's */
    double x,yij,yik,yjk,violl[MAX_INEQ];
    double yi,yj,yk,yy;
    int ind[MAX_INEQ],i,j,k,ii,cntr;
    int T_tmp[MAX_INEQ*4], b_tmp[MAX_INEQ];
	
    unsigned int h_new[MAX_INEQ];

    int flag,typ,ftest=0;
    unsigned int hash_curr;
    unsigned int *h_tmpC;

    int pT_size,T_tmp_size;

    flag=1;
    if ((*pvmax)<1)
      {
		ftest=1;
		(*pvmax)=1;
      };
    if ((*pvmax)>MAX_INEQ)
      (*pvmax)=MAX_INEQ;
    T_tmp_size=*pvmax;
    
    /* dynamic memory allocation */
    h_tmpC=(unsigned int *)malloc(m*sizeof(unsigned int));
    
    for (j=0;j<m;j++)                            /* move h_tmp to h_tmpC */
      h_tmpC[j]=(unsigned int)h_tmp[j+1];   /* (long-array) */
    
    for (i=0;i<(*pvmax);i++) 
      {
		violl[i]=-1.0;
		ind[i]=i;
      }

    cntr=0;
    
    for (i=1;i<=n-2;i++) 
      { 
		for (j=i+1;j<=n-1;j++)
		  { 
			for (k=j+1;k<=n;k++)
			  { 
				yij = X[ijtokp(i,j,n)];
				yik = X[ijtokp(i,k,n)];
				yjk = X[ijtokp(j,k,n)];
				yi = X[ijtokp(i,i,n)];
				yj = X[ijtokp(j,j,n)];
				yk = X[ijtokp(k,k,n)];
				yy = yi + yj + yk;

				/* Typ 0 */
				x =  yy*0.5 + yij + yik + yjk - 2.0;
				hash_curr = (i*500+j)*500+k;
		
				if (m<1) 
				  flag=0;
				else 
				  {       
		    
					if ((hash_curr<h_tmpC[0])||(hash_curr>h_tmpC[m-1])) 
					  flag=0;
					else 
					  { 
			
					if (x>P.viol_tol)
					  flag=binsearch( h_tmpC, hash_curr, 0, m-1);
					else
					  flag=0; 
					  };
				  };
		
		
				if ((x>P.viol_tol) && (flag<1))
				  if (x>violl[ind[0]])
					{ 
					  ii = ind[0];
					  violl[ii] = x;
					  T_tmp[ijtok(ii+1,1,T_tmp_size)] = i; 
					  T_tmp[ijtok(ii+1,2,T_tmp_size)] = j; 
					  T_tmp[ijtok(ii+1,3,T_tmp_size)] = k; 
					  T_tmp[ijtok(ii+1,4,T_tmp_size)] = 0;
					  h_new[ii] = hash_curr;
					  b_tmp[ii] = 2;
					  heapify_n(0,T_tmp_size,ind,violl);
					  cntr++;
					};
		
				/* Typ 1 */
				x = yy*0.5 + yij - yik - yjk;
				typ=1;
				hash_curr = ((typ*500+i)*500+j)*500+k;
				if (m<1) 
				  flag=0;
				else 
				  {
					if ((hash_curr<h_tmpC[0])||(hash_curr>h_tmpC[m-1]))
					  flag=0;
					else 
					  {
					if (x>P.viol_tol)
					  flag=binsearch( h_tmpC, hash_curr, 0, m-1);
					else
					  flag=0; 
					  };
				  };
		
				if ((x>P.viol_tol) && (flag<1))
				  if (x>violl[ind[0]])
					{ 
					  ii = ind[0];
					  violl[ii] = x;
					  T_tmp[ijtok(ii+1,1,T_tmp_size)] = i; 
					  T_tmp[ijtok(ii+1,2,T_tmp_size)] = j; 
					  T_tmp[ijtok(ii+1,3,T_tmp_size)] = k; 
					  T_tmp[ijtok(ii+1,4,T_tmp_size)] = typ;
					  h_new[ii] = hash_curr;
					  b_tmp[ii] = 0;
					  heapify_n(0,T_tmp_size,ind,violl);
					  cntr++;
					};
		
				/* Typ 2 */
				x = yy*0.5 - yij + yik - yjk;
				typ=2;
				hash_curr =((typ*500+i)*500+j)*500+k;
				if (m<1) 
				  flag=0;
				else 
				  {
					if ((hash_curr<h_tmpC[0])||(hash_curr>h_tmpC[m-1])) 
					  flag=0;
					else
					  {
					if (x>P.viol_tol)
					  flag=binsearch( h_tmpC, hash_curr, 0, m-1);
					else
					  flag=0; 
					  };
				  };
		
				if ((x>P.viol_tol) && (flag<1))
				  if (x>violl[ind[0]])
					{ 
					  ii = ind[0];
					  violl[ii] = x;
					  T_tmp[ijtok(ii+1,1,T_tmp_size)] = i; 
					  T_tmp[ijtok(ii+1,2,T_tmp_size)] = j; 
					  T_tmp[ijtok(ii+1,3,T_tmp_size)] = k; 
					  T_tmp[ijtok(ii+1,4,T_tmp_size)] = typ;
					  h_new[ii] = hash_curr;
					  b_tmp[ii] = 0;
					  heapify_n(0,T_tmp_size,ind,violl);
					  cntr++;
					};
		
				/* Typ 3 */
				x = yy*0.5 - yij - yik + yjk;
				typ=3;
				hash_curr =  ((typ*500+i)*500+j)*500+k;
				if (m<1) 
				  flag=0;
				else 
				  {
					if ((hash_curr<h_tmpC[0])||(hash_curr>h_tmpC[m-1]))
					  flag=0;
					else 
					  {
					if (x>P.viol_tol)
					  flag=binsearch( h_tmpC, hash_curr, 0, m-1);
					else
					  flag=0;
					  };
				  };
		
				if ((x>P.viol_tol) && (flag<1))
				  if (x>violl[ind[0]])
					{ 
					  ii = ind[0];
					  violl[ii] = x;
					  T_tmp[ijtok(ii+1,1,T_tmp_size)] = i; 
					  T_tmp[ijtok(ii+1,2,T_tmp_size)] = j; 
					  T_tmp[ijtok(ii+1,3,T_tmp_size)] = k; 
					  T_tmp[ijtok(ii+1,4,T_tmp_size)] = typ;
					  h_new[ii] = hash_curr;
					  b_tmp[ii] = 0;
					  heapify_n(0,T_tmp_size,ind,violl);
					  cntr++;
					};
		
				if (ftest>0 && cntr >0) goto ende;

				/*		for (dd=0;dd<cntr+2;dd++)
				  {
					printf("violl(%d)=%6.4f  ind(%d)=%d\n",ind[dd],violl[ind[dd]],dd,ind[dd]);
					}; printf("*\n"); */
		
			  };
		  };
      };
    
  ende:
    //printf("found additional %d violated triangles\n",cntr);
    if (cntr > *pvmax) 
      cntr=*pvmax;

    /* free */
    free(h_tmpC);
    

    pT_size=cntr;
    if (cntr > P.max_gen_tr)
      pT_size=P.max_gen_tr;
    /* allocate memory for return values */
    *pT=(int *)malloc((pT_size+1)*4*sizeof(int));
    *phash=(unsigned int *)malloc((pT_size+1)*sizeof(unsigned int));
    *pgamma=(double *)malloc((pT_size+1)*sizeof(double));
    if ((*pT == NULL)||(*phash == NULL)||(*pgamma == NULL))
      {
	printf("Storage allocation failed!\n");
	exit(10);
      };

    
    if (cntr > P.max_gen_tr) //(*pvmax > P.max_gen_tr)
      {
	/* we found more violated triangles than we want to add,
	   so we will sort the heap and then take the first mandatory_tr and 
	   randomly choose from the remaining (cntr-mandatory_tr) 
	   if we didn't find enough triangles to fill up the heap, the least violation is in 
	   violl[ind[cntr-1]] */

	heapsort(*pvmax,ind,violl);

	if (P.mandatory_tr < P.max_gen_tr)
	  {
	    for (k=P.mandatory_tr;k<P.max_gen_tr;k++)
	      {
		int h;
		j=gb_unif_rand(cntr-k); //j=gb_unif_rand(*pvmax-k);
		h=ind[k+j];
		ind[k+j]=ind[k];
		ind[k]=h;
		//printf("%d %d\n",k,k+j);
	      };

	  };
	*pvmax=P.max_gen_tr;
	cntr=P.max_gen_tr;
      };

    i=0;
    for (k=0;k<(*pvmax);k++)
      {
	j = ind[k];
	if (violl[j]>-1.) 
	  {    
	    (*pgamma)[i+1] = violl[j];
	    //printf(" %8.4f \n",violl[j]);
	    (*phash)[i+1] = h_new[j];
	    (*pT)[ijtok(i+1,1,pT_size)] = T_tmp[ijtok(j+1,1,T_tmp_size)];  
	    (*pT)[ijtok(i+1,2,pT_size)] = T_tmp[ijtok(j+1,2,T_tmp_size)];
	    (*pT)[ijtok(i+1,3,pT_size)] = T_tmp[ijtok(j+1,3,T_tmp_size)];
	    (*pT)[ijtok(i+1,4,pT_size)] = T_tmp[ijtok(j+1,4,T_tmp_size)];
	    i++;
	  };
      }


    *pvmax=i; /* length of hash array */

    return(0);
}
