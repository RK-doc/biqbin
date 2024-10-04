@* Main Function for Branch and Bound Skeleton
   for quadratic knapsack problem

@c

  @<include files@>@/
@;

@<primal heuristics@>@/

int main()
{@/ 
  double constant = 0;
  @<variables for the main@>@/
  @<read the data and call bab_sdp@>@/
  @<output the result@>@/
  return 0;
@;
}

@ @<include files@>=
#include <stdio.h>
#include <stdlib.h>

#define EPS 1e-07

int bab_sdp(int,int,int*,int*,double*,
	    int,int*,double*,int*,double*,int*,
	    int*,int,int,int*,void*);

@ @<primal heuristics@>=
int primal(int n2,double *bound_cut, double **cost, void *user_info) {
  int i, bpos;
  double viol, best;
  double *orig;
  orig = (double*) user_info;
  viol = -orig[n2-1];
  for(i = 0; i < n2-1; i++) if(bound_cut[i+1]*bound_cut[0] < 0) viol += orig[i];
  best = 1.0;
  while(viol > EPS && best > EPS) {
    best = -1.0;
    for(i = 0; i < n2-1; i++) if(-orig[i]*bound_cut[i+1]*bound_cut[0] > best) {
      // better: lhs divided by increase of obj function
      best = -orig[i]*bound_cut[i+1]*bound_cut[0];
      bpos = i+1;
    }
    if(best > EPS) {
      bound_cut[bpos] = -bound_cut[bpos];
      viol -= best;
    }
  }
  return 1;
};

int repair_cut(int n,double *bound_cut,int con,
               int *eq,double *rhs,int *bg,int *sp,double *co) {
  double viol, best;
  int i, bpos;
  if(rhs[0] < 0) return 0;
  viol = -rhs[0];
  for(i = 0; i < bg[1]; i++) if(bound_cut[sp[i]+1]*bound_cut[0] < 0) viol += co[i];
  best = 1.0;
  while(viol > EPS && best > EPS) {
    best = -1.0;
    for(i = 0; i < bg[1]; i++) if(-co[i]*bound_cut[sp[i]+1]*bound_cut[0] > best) {
      // better: lhs divided by increase of obj function
      best = -co[i]*bound_cut[sp[i]+1]*bound_cut[0];
      bpos = sp[i]+1;
    }
    if(best > EPS) {
      bound_cut[bpos] = -bound_cut[bpos];
      viol -= best;
    }
  }
  return 1;
};

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
  int *sense,*bg,*sp;
  double *rhs,*co;
  double lhs;
  double w;
  double *orig;
  scanf("%d %d",&n,&m);
  _h = (int*) malloc (sizeof(int) * (m+n));
  _t = (int*) malloc (sizeof(int) * (m+n));
  _w = (double*) malloc (sizeof(double) * (m+n));
  for (i=0; i<n; i++) {
    _h[i] = 1;
    _t[i] = i+2;
    _w[i] = 0;
  }
  for (i=0; i<m; i++) {
    scanf("%d %d %lf",&_h[n+i],&_t[n+i],&w);
    _h[n+i]++;
    _t[n+i]++;
    _w[_h[n+i]-2] += w;
    _w[_t[n+i]-2] += w;
    _w[n+i] = -w;
    constant -= w;
  }
  start_cut = (int*) malloc (sizeof(int) * (1+n));
  sense = (int*) malloc (sizeof(int));
  bg = (int*) malloc (sizeof(int) * 2);
  co = (double*) malloc (sizeof(double) * n);
  sp = (int*) malloc (sizeof(int) * n);
  rhs = (double*) malloc (sizeof(double));
  orig = (double*) malloc (sizeof(double) * (n+1));
  sense[0] = -1;
  bg[0] = 0;
  start_cut[0] = 0;
  lhs = 0;
  for(i=0; i<n; i++) {
    scanf("%lf",&co[i]);
    orig[i] = co[i];
    sp[i] = i+1;
    start_cut[i+1] = (co[i] < 0);
    lhs += co[i]*start_cut[i+1];
    co[i] = -co[i];
  }
  scanf("%lf",&rhs[0]);
  orig[n] = rhs[0];
  if(rhs[0] < lhs-1e-6) {
    printf("problem is infeasible!\n");
    exit(0);
  }
  rhs[0] *= 2.0;
  for(i=0; i<n; i++) rhs[0] += co[i];
  bg[1] = n;
  best_cut = (int*) malloc (sizeof(int) * (1+n));
  bab_sdp(1+n,n+m,_h,_t,_w,
          1,sense,rhs,bg,co,sp,
	  best_cut,0,0,start_cut,orig);
  free(start_cut);
  free(sense);
  free(bg);
  free(co);
  free(sp);
  free(rhs);
  free(orig);
}

@ @<output the result@>=
{
  int i;
  double best_value = constant;
  for(i=0; i<n+m; i++)
    if(best_cut[_h[i]+1]^best_cut[_t[i]+1]) best_value -= _w[i];
    else best_value += _w[i];
  free(_h);
  free(_t);
  free(_w);
  printf("%lf\n",-best_value/4);
  for (i=0; i<n; i++)
    printf(" x(%d) = %d\n",i+1,best_cut[i+1]^best_cut[0]);
  free(best_cut);
}

@*

@* Index.


