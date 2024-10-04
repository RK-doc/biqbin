@* Main Function for Branch and Bound Skeleton.

@c

  @<include files@>@/
@;
int 
main(int argc, char *argv[])
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
#include <declarations.h>

int bab_sdp(int,int,
	    int,int*,double*,int*,double*,int*,
	    int*,int,int,void*,int,char**);

@ Macro definition

@ @<variables for the main@>=
int    n,
       mm;
int    *_h,
       *_t,
       *best_cut;
double *_w;

 
@ @<read the data@>=
{ 
}

@ @<call bab_sdp@>=
{
  bab_sdp(n,mm, /* the graph */ 
	  0,0,0,0,0,0, /* the constraints */
	  best_cut,      /* the returned max-cut */
	  -1,             /* the print level */
	  1,             /* check integrality */
      0,            /* the user information */ 
	  argc,
	  argv);
}

@ @<output the result@>=
{
  
}

@*

@* Index.


