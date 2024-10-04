@* Main Function for Branch and Bound Skeleton
   for the independent set problem

@c
  @<include files@>@/
@;
int 
main()
{@/ 
  @<variables for the main@>@/
  @<read the data and call bab_sdp@>@/
  @<output the result@>@/
  return 0;
@;
}

int primal(int n,double *bound_cut,double *cost,void *user_info) { 
  return 1;
};

int repair_cut(int n,double *bound_cut,int con,
               int *eq,double *rhs,int *bg,int *sp,double *co) {
  return 0;
};

@ @<include files@>=
#include <stdio.h>
#include <stdlib.h>

int bab_sdp(int,int,int*,int*,double*,
	    int,int*,double*,int*,double*,int*,
	    int*,int,int,int*,void*);

@ @<variables for the main@>=
int    n,
       m;
int    *_h,
       *_t,
       *start_cut,
       *best_cut;
double *_w;
 
@ @<read the data and call bab_sdp@>=
{ int i;
  scanf("%d %d",&n,&m);
  _h = (int*) malloc (sizeof(int) * (m+n));
  _t = (int*) malloc (sizeof(int) * (m+n));
  _w = (double*) malloc (sizeof(double) * (m+n));
  for (i=0; i<n; i++) {
    _h[i] = 1;
    _t[i] = i+2;
    _w[i] = 2.0;
  }
  for (i=0; i<m; i++) {
    scanf("%d %d",&_h[n+i],&_t[n+i]);
    _h[n+i]++;
    _t[n+i]++;
    _w[_h[n+i]-2] -= 1.0;
    _w[_t[n+i]-2] -= 1.0;
    _w[n+i] = 1.0;
  }
  start_cut = (int*) malloc (sizeof(int) * (1+n));
  for (i=0; i<=n; i++) start_cut[i]= 0;
  best_cut = (int*) malloc (sizeof(int) * (1+n));
  bab_sdp(1+n,n+m,_h,_t,_w,
          0,0,0,0,0,0,
	  best_cut,0,0,start_cut,0);
  free(start_cut);
}

@ @<output the result@>=
{
  int i;
  for(i=0; i<m; i++)
    if((best_cut[_h[n+i]-1]^best_cut[0])
       &&(best_cut[_t[n+i]-1]^best_cut[0]))
      best_cut[_h[n+i]-1] = !best_cut[_h[n+i]-1];
  free(_h);
  free(_t);
  free(_w);
  for (i=0; i<n; i++) if(best_cut[i+1]^best_cut[0]) printf(" %d",i+1);
  printf("\n");
  free(best_cut);
}

@*

@* Index.


