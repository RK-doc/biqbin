// definicije spremenljivk in funkcij, v komentarjih zapisano kje se uporabijo
#define EPSILON 1.0e-5   /* bound, mc_2opt, lam_eta */
#define EPSILON6 1.0e-6  /* solvelambda, mc_psdpk */
#define EPSILON8 1.0e-8  /* mc_cut */
#define MAX_INEQ 50000   /* bound, bdl_main, bdl_mc2, fct_eval, tri_sep2, separation */
#define MAX_BDL 400      /* bdl_main, bdl_mc2 */
#define COMPARE(x,y) (((x)<(y)) ? -1: ((x)==(y))? 0: 1)
#define ijtokp(iii,jjj,lda) ((iii+jjj*(jjj-1)/2)-1) 
#define ijtokr(iii,jjj,lda) ((jjj-1)*(jjj-2)/2+iii-1)
#define ijtok(iiii,jjjj,lda) ((jjjj-1)*lda+iiii-1)
#define Abs(i) ((i)>0 ? (i):-(i))
#define TAIL_PERC 0.05
#define MAX_TAIL 5

#define EPS 1e-8

#define e2l(i,j) (i>=j?((i+1)*(i))/2+(j):((j+1)*(j))/2+i)  

#define l2i2(i) ((int) (-0.5+sqrt(2*(double) (i)+0.25)+EPS)) 

#define l2i1(i) (i-e2l(0,l2i2(i))) 

#define Round(a) (a>=0?(int) (a+0.5) :(int) (a-0.5) ) 

#define Positive(a) (a> EPS) 

#define Negative(a) (a<-EPS) 

#define Zero(a) (!Positive(a) &&!Negative(a))

#define Integer(a) (Zero((a)-Round(a)))

typedef struct {
  int    max_bb_nodes;
  double gap_relax;
  int    inner_it;
  int    max_inner_it;
  int    outer_it;
  int    extra_it;
  int    min_outer_it;
  double viol_tol;
  int    heap_size;
  int    mandatory_tr;
  int    max_gen_tr;
  int    printlev;
  int    branch_rule;
               } T_param; 


/* the comments give the progs where the routines are called */

/* our own routines */
int valid_cut(int n,double *cut,int con,int *eq,double *rhs,
	      int *bg,int *sp,double *co);
int op_b(int n,int con,int *eq,int *bg,int *sp,double *co,double *X,
	 double **pbx);
int op_bt(int n,int con,int *eq,int *bg,int *sp,double *co,double *nu,
	  double **pBtn);
void min_two(int k,double *x,int *argmin1,int *argmin2); /* bound */

void diag(double *X,double **py,int n);      /* mc_1opt, mc_psdpk */

void Diag(double **pX,double *y,int n);      /* mc_psdpk */

double traceprod(double *A,double *B,int n); /* mc_psdpk, fct_eval, bdl_mc2 */

void deltamax(double *diagL,double *x,double *Lx,
              int n,double *best,int *bestarg);            /* mc_1opt */

int mc_1opt(double *L,int n,double **px,double *pcost);    /* mc_cut */

int mc_cut(double *L,int n,double *v,int k,int trials,
	   int con,int *eq,double *rhs,int *bg,int *sp,double *co,
           double **pcut,double *pcost);    /* heuristic */

int heuristic(double *L,double *X,int n,
              double **pxh,double *pfh,
	      int con,int *eq,double *rhs,int *bg,int *sp,double *co); 
                                                           /* bdl_main */

void my_dgsmm(int *pn,double *pA,double *pB);              /* mc_psdpk */

int mc_at(int n,int m,int m_mem,int *T,double *gamma,
          double **pAtg);                             /* fct_eval */

int mc_a(int n,int m,int m_mem,int *T,double *X,double **pax);  
                                                    /* fct_eval, bdl_mc2 */

int fct_eval(int n,int m,double *gamma,double *mu,
	     int con,int *eq,double *rhs,int *bg,int *sp,double *co,
	     double *L,int *T,
             double *pf,double **pX,double **pg,double **ph);   
                                                    /* bdl_main, bdl_mc2 */

int mc_psdpk(double *L,int n,int silent,double **pX,
             double **py,double *pphi);                 /* fct_eval */

int solvelambda(double *Q,double *c,int k,double **plam,
                int silent);                               /* lam_eta */

int tri_sep2(T_param P,double *X,int n,int *pvmax,unsigned int *h_tmp,
             int m,int **pT,unsigned int **phash,
             double **pgamma);                             /* separation */

int separation(T_param P,double *X,int n,int *pm,int pT[MAX_INEQ*4],unsigned int **phash,
               double **pgamma_new);                      /* bdl_main */

int lam_eta(int k,int m,double *beta,double *G,double *gamma,
            double *pt,double **pd,double **plam);        /* bdl_mc2 */

double steplength(int k,double *dx,double *x);            /* solvelambda */

void dmax_val(int k,double *x,int *argmax,double *max);   /* steplength, bdl_mc2 */

void dmin_val(int k,double *x,int *argmin,double *min);   /* solvelambda */


int bdl_mc2(int n,int m,double *L,int *T,unsigned int *hash,
	    int con,int *eq,double *rhs,int *bg,int *sp,double *co,
	    double best,
            double **pgamma,double **pmu,int itmax,double *pt,double **pxi,
            double **pX,double *pfopt,double **px,int *pkin,
            int printlevel,double gap_relax);             /* bdl_main */

int bdl_main(T_param P,double *L,int n,int forec,
	     int con,int *eq,double *rhs,int *bg,int *sp,double *co,
             double *pbest,double *pf,double **pcut,double **primalX,
             int *pm,int T[MAX_INEQ*4],unsigned int *hash,double *gamma);
                                                          /* bound */

/* level 1 blas */
void dscal_();      /* fct_eval, bdl_main, bdl_mc2 */
void dcopy_();      /* mc_cut, heuristic, mc_psdpk, fct_eval, solvelambda, lam_eta, bdl_mc2 */
void scopy_();      /* bound */
void daxpy_();      /* heuristic, mc_psdpk, fct_eval, solvelambda, separation, bdl_mc2 */
double ddot_();     /* mc_1opt, heuristic, fct_eval, solvelambda, bdl_main */
double dnrm2_();    /* lam_eta */

/* level 2 blas */
void dgemv_();      /* mc_cut, lam_eta, bdl_mc2 */
void dspmv_();      /* mc_1opt, heuristic, mc_psdpk, solvelambda */
void dspr_();       /* heuristic */

/* level 3 blas */
void dsymm_();      /* my_dgsymm */
void dsyrk_();      /* lam_eta */

/* lapack */
void dpptrf_();     /* heuristic, mc_psdpk */
void dpptri_();     /* mc_psdpk */
void dspsv_();      /* mc_psdpk, solvelambda */


// OPIS BLAS:
// The BLAS (Basic Linear Algebra Subprograms) are routines that provide standard building
// blocks for performing basic vector and matrix operations. The Level 1 BLAS perform scalar,
// vector and vector-vector operations, the Level 2 BLAS perform matrix-vector operations,
// and the Level 3 BLAS perform matrix-matrix operations. Because the BLAS are efficient,
// portable, and widely available, they are commonly used in the development of high quality
// linear algebra software, LAPACK for example

// klièemo funkcijo iz BLAS
// glej: https://stackoverflow.com/questions/14470799/calling-ddot-function-in-blas-library
// v komentarjih so funkcije, ki uporabijo dane BLAS funkcije
