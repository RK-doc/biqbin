@* Main Function for Branch and Bound Skeleton.

@c

  @<include files@>@/
@;
int 
main()
{@/ 
  @<variables for the main@>@/
  @<read the data@>@/
  @<call bab_sdp@>@/
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
       mm;
int    *_h,
       *_t,
       *start_cut,
       *best_cut;
double *_w;  
int    con,
       *eq,
       *bg,
       *sp;
double *rhs,
       *co;

 
@ @<read the data@>=
{ int i,j,k,nz,dmy;
  scanf("%d %d %d %d",&n,&mm,&con,&nz);
  _h = (int*) malloc (sizeof(int) * mm);
  _t = (int*) malloc (sizeof(int) * mm);
  _w = (double*) malloc (sizeof(double) * mm);
  for (i=0; i<mm; i++) {
    scanf("%d %d %lf",&_h[i],&_t[i],&_w[i]);
    if(_h[i] < 1 || _h[i] > n) {
      printf(" INPUT ERROR: node %d of edge %d is out of range.\n",_h[i],i);
      exit(11);
    }
    if(_t[i] < 1 || _t[i] > n) {
      printf(" INPUT ERROR: node %d of edge %d is out of range.\n",_t[i],i);
      exit(11);
    }
  }
  eq = (int*) malloc (sizeof(int) * con);
  bg = (int*) malloc (sizeof(int) * (con+1));
  rhs = (double*) malloc (sizeof(double) * con);
  sp = (int*) malloc (sizeof(int) * nz);
  co = (double*) malloc (sizeof(double) * nz);
  bg[0]=0;
  for (i=0; i<con; i++) {
    scanf("%d %d %lf",&eq[i],&k,&rhs[i]); /* =,<=,>=flag, number of coeff in constr. i, rhs */
    bg[i+1]=bg[i]+k;
    for (j=0; j<k; j++) {
      scanf("%d %lf",&dmy,&co[bg[i]+j]);
      if((dmy>mm)||(dmy<0)) {
	printf(" INPUT ERROR: edge number out of range.\n");
	exit(11);
      }
      sp[bg[i]+j]=dmy;
    }
  }
  if (bg[con]!=nz) {
    printf(" INPUT ERROR: inconsistency in number of nonzero coefficients in the constraints.\n");
    exit(11);
  }
}

@ @<call bab_sdp@>=
{ int i;
  start_cut = (int*) malloc (sizeof(int) * n);
  for(i=0; i<n; i++) 
    start_cut[i] = -1;
  best_cut = (int*) malloc (sizeof(int) * n);
  bab_sdp(n,mm,_h,_t,_w, /* the graph */ 
	  con,eq,rhs,bg,co,sp, /* the constraints */
	  best_cut,      /* the returned max-cut */
	  -1,            /* the print level */
	  1,             /* check integrality */
	  start_cut,     /* the initial feasible solution */
          0);            /* the user information */ 
  free(start_cut);
}

@ @<output the result@>=
{
  int i;
  double best_value = 0.0;
  for(i=0; i<mm; i++)
    if(best_cut[_h[i]-1] != best_cut[_t[i]-1]) 
      best_value += _w[i];
  free(_h);
  free(_t);
  free(_w);
  printf("cut value: %lf\n",best_value);
  printf("one side of the cut:\n");
  for (i=0; i<n; i++)
    if (best_cut[i] == -1) printf(" %d",i + 1);
  printf("\n");
  free(best_cut);
}

@*

@* Index.


