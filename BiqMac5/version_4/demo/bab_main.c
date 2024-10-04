/*1:*/
#line 3 "./demo/bab_main.w"


/*2:*/
#line 27 "./demo/bab_main.w"

#include <stdio.h> 
#include <stdlib.h> 
#include <declarations.h> 

int bab_sdp(int,int,
int,int*,double*,int*,double*,int*,
int*,int,int,void*,int,char**);

/*:2*/
#line 5 "./demo/bab_main.w"


int
main(int argc,char*argv[])
{
/*4:*/
#line 38 "./demo/bab_main.w"

int n,
mm;
int*_h,
*_t,
*best_cut;
double*_w;


/*:4*/
#line 10 "./demo/bab_main.w"

/*5:*/
#line 47 "./demo/bab_main.w"

{
}

/*:5*/
#line 11 "./demo/bab_main.w"

/*6:*/
#line 51 "./demo/bab_main.w"

{
bab_sdp(n,mm,
0,0,0,0,0,0,
best_cut,
-1,
1,
0,
argc,
argv);
}

/*:6*/
#line 12 "./demo/bab_main.w"

/*7:*/
#line 63 "./demo/bab_main.w"

{

}

/*:7*/
#line 13 "./demo/bab_main.w"

return 0;

}

int primal(int n,double*bound_cut,double*cost,void*user_info){
return 1;
};

int repair_cut(int n,double*bound_cut,int con,
int*eq,double*rhs,int*bg,int*sp,double*co){
return 0;
};

/*:1*/
