/*
                       heapsort.cpp
*/

#include "heapsort.h"

/*
//--------------------------------------------------------------
//                          heapsort
//--------------------------------------------------------------
//sorts the list of doubles non decreasingly by rearranging the
//indices in the array ind. n is the number of elements
*/

void heapify(i,n,ind,d)
int i,n,*ind;
double *d;
{
 int k,h,j;

 ind--; i++;
 k=i;
 while(2*k<=n){
     if (2*k==n){
         if (d[ind[k]]<d[ind[2*k]]){
             h=ind[k];
             ind[k]=ind[2*k];
             ind[2*k]=h;
         }
         break;
     }
     if ((d[ind[k]]<d[ind[2*k]])||
         (d[ind[k]]<d[ind[2*k+1]])){
         if (d[ind[2*k]]>d[ind[2*k+1]]) j=2*k;
         else j=2*k+1;
         h=ind[k];
         ind[k]=ind[j];
         ind[j]=h;
         k=j;
     }
     else break;
 }
}

void build_heap(n,ind,d)
int n,*ind;
double *d;
{
 int i;
 for(i=n/2-1;i>=0;i--) heapify(i,n,ind,d);
}

void heapsort(n,ind,d)
int n,*ind;
double *d;
{
 int i,h;

/* //build heap  */
 build_heap(n,ind,d);

/* //extract greatest and rebuild heap  */
 for(i=n-1;i>=1;i--){
     h=ind[i];ind[i]=ind[0];ind[0]=h;
     heapify(0,i,ind,d);
 }
}
