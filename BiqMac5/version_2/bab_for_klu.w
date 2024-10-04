@* A Branch and Bound Skeleton.

@c

  @<utility functions@>@/
  @<include files@>@/
  @<data structures@>@/
@;

void bound(T_param,int,double*,
	   int,int*,double*,int*,int*,double*,
	   double,double*,double*,double*,int*,int*,int*);

int primal(int,double*,double**,void*);

int 
bab_sdp(int n,int mm,int *head,int *tail,double *weight,
        int num_con,int *sense,double *rhs,int *bg,double *co,int *sp,
        int *best_cut,int prlev,int check_integrality,void *user_info)
{@/ 
  @<variables for bab_sdp@>@/
  @<read the data@>@/
  @<read parameters@>@/
  @<initializations@>@/
  @<main loop@>@/
  @<output the result@>@/
  @<free@>@/
  return 0;
@;
}

@ @<include files@>=
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <gb_flip.h>
#include <declarations.h>

@ @<data structures@>=

typedef struct tree_node {
  int name;          /* node name */
  int flip;          /* -1 if this node represents the flipped contraction; +1 otherwise */
  int level;         /* node level (0 at root) */
  int n_sons;        /* number of sons */
  int log_shrink;    /* number of logical shrinkings in the node */
  struct tree_node *father; /* pointer to the father node */
  double  correction;/* switching correction (0 at root) */

  double  ub;        /* upper bound of this node */
  int left,          /* left  end node of edge to be contracted */
      right;         /* right end node of edge to be contracted */
               } Tree_node;

typedef struct active_node{
  Tree_node   *node;
  struct active_node *pred;
  struct active_node *succ;
               } Active_node;

typedef struct {
  int      n;
  int     *alive;
  int     *ref;
  int     *flip;
  int     *link;
  double **m;
               } Matrix;

typedef struct {
  int   left;
  int   right;
  int   flip;
               } Stack;

@ Macro definitions

@d BIQMAC 0

@d alloc_mat(A,n)
{
  int i;
  A.n = n;
  A.alive = (int *) malloc (sizeof(int) * n);
  A.ref = (int*) malloc (sizeof(int) *n);
  A.flip = (int*) malloc (sizeof(int) *n);
  A.link = (int*) malloc (sizeof(int) *n);
  for (i=0; i<n; i++) {
    A.alive[i]=1;
    A.flip[i]=1;
  }
  A.m = (double **) malloc (sizeof(double) * n);
  A.m[1] = (double *) calloc ((n*(n-1)/2), sizeof(double));
  for (i=2; i<n; i++)
    A.m[i] = A.m[i-1] + i - 1;
}

@d free_mat(A)
{
  free(A.m[1]);
  free(A.m);
  free(A.ref);
  free(A.flip);
  free(A.link);
  free(A.alive);
}

@d shrink(A,i,j,f)
{
  int small, large, l;
  small = i < j ? i : j;
  large = i < j ? j : i;

  star_value[small] = f * (star_value[small] - A.m[large][small]) + 
                          (star_value[large] - A.m[large][small]); 
  if (f == -1) {
    for (l=0; l<small; l++)
      if (A.alive[l] && (l != large))
        star_value[l] -= (2 * A.m[small][l]);
    for (l=small+1; l<n; l++)
      if (A.alive[l] && (l != large))
        star_value[l] -= (2 * A.m[l][small]);
  }
  for (l=0; l<small; l++) 
    A.m[small][l] = f * A.m[small][l] + A.m[large][l];
  for (l=small+1; l<large; l++) 
    A.m[l][small] = f * A.m[l][small] + A.m[large][l];
  for (l=large+1; l<A.n; l++) 
    A.m[l][small] = f * A.m[l][small] + A.m[l][large];
  A.alive[large] = 0;
  A.ref[large] = small;
  A.flip[small] *= f;
  A.link[large] = A.flip[small];
}

@d shrink_by_log_imp(A,i,j,f)
{
  if(f == -1)
    cur_node->correction += star_value[i];
  if(P.printlev >= 1) {
    printf("shrinking [%d,%d] by logical implications, ",f*(i+1),j+1);
    printf("correction becomes %f\n",cur_node->correction);
  }
  shrink(A,i,j,f);
  log_stack[log_n].left = i;
  log_stack[log_n].right = j;
  log_stack[log_n].flip = f;
  log_n++;
  cur_node->log_shrink++;
}

@d unshrink(A,i,j,f)
{
  int small, large, l;
  small = i < j ? i : j;
  large = i < j ? j : i;
  for (l=0; l<small; l++) 
    A.m[small][l] = f * (A.m[small][l] - A.m[large][l]);
  for (l=small+1; l<large; l++) 
    A.m[l][small] = f * (A.m[l][small] - A.m[large][l]);
  for (l=large+1; l<A.n; l++) 
    A.m[l][small] = f * (A.m[l][small] - A.m[l][large]);

  if (f == -1) {
    for (l=0; l<small; l++)
      if (A.alive[l] && (l != large))
        star_value[l] += (2 * A.m[small][l]);
    for (l=small+1; l<n; l++)
      if (A.alive[l] && (l != large))
        star_value[l] += (2 * A.m[l][small]);
  }
  star_value[small] = f * (star_value[small] - star_value[large] 
                           + A.m[large][small]) + A.m[large][small];

  A.alive[large] = 1;
  A.ref[large] = -1;
  A.flip[small] *= f;
  A.link[large] = 0;
}

@d node_info(p)
{
  if (p) {
    if (p->father) {
      printf(" #%d son of %d [%d,%d] <ub %4.1f> ^corr. %4.1f^ {gap %5.2f} time: %.2f",p->name,
             (p->father)->name,p->flip*((p->father)->left+1),
             (p->father)->right+1, p->ub,p->correction, 
             upper_bound-best_value, (cputime()-starttime)/100.);
/*
      printf(" #%d,%d (%d) [%d,%d] <%4.1f> ^%4.1f^ {%5.2f}",p->name,p->n_sons,
             (p->father)->name,p->flip*((p->father)->left+1),
             (p->father)->right+1, p->ub,p->correction, 
             upper_bound-best_value);
*/
    }
    else {
/*      printf(" Root *%d <%4.1f>",p->n_sons,p->ub); */
      printf(" Root <ub %4.1f>",p->ub);
    }
  }
  else {
    printf(" No such node!! ");
  }
}

@d Getparam(a,b)
  if (!(strcmp (st, #a))) fscanf (pf, #b, &(P.a))

@d con_info(i)
{
  int j;
  for(j= tmp_bg[i];j<tmp_bg[i+1];j++)
    printf(" + (%f)*(%d,%d)",tmp_co[j],l2i1(tmp_sp[j]),l2i2(tmp_sp[j]));
  if(tmp_sense[i] == -1) printf(" <= ");
  else if(tmp_sense[i] == 1) printf(" >= ");
  else printf(" == ");
  printf("%f\n",tmp_rhs[i]);
}

@d compute_orig_cost(v,c)
{
  int i, j;
  v = 0.0;
  for (i=0; i<n; i++) {
    for (j=0; j<i; j++) {
      v += orig_cost.m[i][j]*(double)(c[i]*c[j] < 0);
    }
  }
}

@ @<variables for bab_sdp@>=

  int          m,
               n_supernodes,
               n_bbnodes;
  int         *original_name;
  double      *star_value;
  Tree_node   *cur_node,
              *prev_node,
              *kill_node;
  Active_node *active,
              *last_active;
  Stack       *stack,
              *log_stack;
  int          log_n;
  Matrix       cost,
               tmp_cost,
               orig_cost;
  double       best_value,
               upper_bound,
               bound_ub,
               bound_lb,
               prim_lb;
  int          bound_left,
               bound_right,
               bound_fix;
  double      *bound_cut;
  double      *prim_cut;
  T_param      P;
  int          starttime;
  int          tmp_con;
  int         *tmp_sense;
  double      *tmp_rhs;
  int         *tmp_bg;
  double      *tmp_co;
  int         *tmp_sp;
  int         *new_pos;
  int          p_res;

@ @<read the data@>=
{ int i, j, k, large, small;
  double c;
  int mult_edges=0;
  int loops=0;
  m = (n * (n - 1)) / 2;
  alloc_mat(cost,n);
  alloc_mat(orig_cost,n);
  P.gap_relax = 1.-EPS;
  for (i=0; i<mm; i++) {
    j = head[i]-1;
    k = tail[i]-1;
    c = weight[i];
    if (j!=k) {
      if(check_integrality && !Integer(c)) 
        P.gap_relax = EPS;
      large = (j > k ? j : k);
      small = (j > k ? k : j);
      if (!Zero(cost.m[large][small]))
	mult_edges++;
      cost.m[large][small] += c; 
      orig_cost.m[large][small] += c;
    }
    else	     
      loops++;
  }
  if (loops > 0)
      printf("Warning: the graph contains loops. Loops have been ignored.\n");
  if (mult_edges > 0)
      printf("Warning: multiple edges! Edge weights have been summed up.\n");
#ifdef BIQMAC
  printf("The input graph has %d nodes and %d edges, ",n,mm-mult_edges-loops);
  if (P.gap_relax > EPS+EPS)
      printf("edge weights are integers.\n\n");
  else
      printf("edge weights are reals.\n\n");
#endif
  if(num_con) {
    for (i=0; i<bg[num_con]; i++) {
      sp[i]--;
      if(head[sp[i]] < tail[sp[i]]) sp[i] = (tail[sp[i]]*(tail[sp[i]]-1))/2+head[sp[i]]-1;
      else sp[i] = (head[sp[i]]*(head[sp[i]]-1))/2+tail[sp[i]]-1;
    }
  }
}

@ @<read parameters@>=
{ FILE *pf;
  int eof;
  char st[128];

  P.max_bb_nodes   = 1000000;
  P.inner_it       = 3;
  P.max_inner_it   = 10;
  P.outer_it       = 5;
  P.extra_it       = 0;
  P.min_outer_it   = 5;
  P.viol_tol       = 1.0e-3;
  P.heap_size      = 1000;
  P.mandatory_tr   = 1000;
  P.max_gen_tr     = 1000;
  P.branch_rule    = 3;
  P.printlev       = 1;
  pf = fopen ("param", "r");
  if (pf) {
    do {
      eof = (fscanf (pf, "%s", st) == EOF);
      Getparam (max_bb_nodes,    %d);
      Getparam (gap_relax,       %lf);
      Getparam (inner_it,        %d);
      Getparam (max_inner_it,    %d);
      Getparam (outer_it,        %d);
      Getparam (extra_it,        %d);
      Getparam (min_outer_it,    %d);
      Getparam (viol_tol,        %lf);
      Getparam (heap_size,       %d);
      Getparam (mandatory_tr,    %d);
      Getparam (max_gen_tr,      %d);
      Getparam (printlev,        %d);
      Getparam (branch_rule,     %d);
    } while (strcmp (st, "end.") && !eof);
    fclose (pf);
  }
  if(prlev >= 0) P.printlev = prlev;
  if (P.mandatory_tr > P.heap_size ||
      P.mandatory_tr > P.max_gen_tr ||
      P.max_gen_tr   > P.heap_size    ) {
    printf ("parameter inconsistency: mandatory_tr < max_gen_tr < heap_size\n"); 
    exit(1);
  }
  if(P.printlev >= 2) {
    printf ("max_bb_nodes =    %d\n", P.max_bb_nodes);
    printf ("gap_relax =       %lf\n",P.gap_relax);
    printf ("inner_it =        %d\n", P.inner_it);
    printf ("max_inner_it =    %d\n", P.max_inner_it);
    printf ("outer_it =        %d\n", P.outer_it);
    printf ("extra_it =        %d\n", P.extra_it);
    printf ("min_outer_it =    %d\n", P.min_outer_it);
    printf ("viol_tol =        %f\n", P.viol_tol);
    printf ("heap_size =       %d\n", P.heap_size);
    printf ("mandatory_tr =    %d\n", P.mandatory_tr);
    printf ("max_gen_tr =      %d\n", P.max_gen_tr);
    printf ("printlev =        %d\n", P.printlev);
    printf ("branch_rule =     %d\n", P.branch_rule);
    printf ("\n\n");
  }
}

@ @<initializations@>=
{ int i, j;
  double cij;
 
  starttime = cputime();

  original_name = (int *) malloc (sizeof(int) * n);
  star_value = (double *) malloc (sizeof(double) * n);
  bound_cut = (double *) malloc (sizeof(double) * n);
  prim_cut = (double *) malloc (sizeof(double) * n);
  alloc_mat(tmp_cost,n);

  if(num_con) {
    tmp_sense = (int *) malloc (sizeof(int) * num_con);
    tmp_rhs = (double *) malloc (sizeof(double) * num_con);
    tmp_bg = (int *) malloc (sizeof(int) * (num_con + 1));
    tmp_co = (double *) malloc (sizeof(double) * bg[num_con]);
    tmp_sp = (int *) malloc (sizeof(int) * bg[num_con]);
  }
  new_pos = (int *) malloc (sizeof(int) * n);

  n_bbnodes = 0;
  best_value = 0.0;
  upper_bound = 0.0;
  
  for (i=0; i<n; i++) {
    star_value[i] = 0.0;
    for (j=0; j<i; j++) {
      cij = cost.m[i][j];
      star_value[i] += cij;
      star_value[j] += cij;
      if (cij > 0) upper_bound += cij;
    }
  }
  compute_orig_cost(best_value,best_cut);
  for (i=0; i<n; i++) prim_cut[i] = (double)best_cut[i];
  p_res = primal(n,prim_cut,orig_cost.m,user_info);
  if (p_res > 0 && valid_cut(n,prim_cut-1,num_con,sense,rhs,bg,sp,co)) {
    compute_orig_cost(prim_lb,prim_cut);
    if(prim_lb > best_value) {
      best_value = prim_lb;
      for (i=0; i<n; i++) best_cut[i] = Round(prim_cut[i]);
    }
  }
  if(P.printlev >= 2) {
    printf("initialize best cut -> %f\n",best_value);
    for(i= 0;i<n;i++) printf(" %d",best_cut[i]);
    printf("\n");
  }
  bound_ub = upper_bound;
  gb_init_rand(1);
  
  cur_node = (Tree_node *) malloc (sizeof(Tree_node));
  cur_node->name = n_bbnodes++;
  cur_node->flip = 1;
  cur_node->level = 0;
  cur_node->n_sons = 0;
  cur_node->log_shrink = 0;
  cur_node->father = NULL;
  cur_node->correction = 0.0;
  cur_node->ub = upper_bound;

  prev_node = cur_node;

  active = (Active_node *) malloc (sizeof(Active_node));
  active->node = cur_node;
  active->pred = NULL;
  active->succ = NULL;
  last_active = active;
  stack = (Stack *) malloc (sizeof(Stack) * n);
  log_stack = (Stack *) malloc (sizeof(Stack) * n);
  log_n = 0;
}

@ @<main loop@>=
{
  int fathomed;
  Tree_node *common_ancestor;
  while(p_res < 2 && active) {
    @<remove from active list@>@;
    @<find common ancestor@>@;
    @<undo shrinkings@>@;
    @<do shrinkings@>@;
    @<prepare constraints@>@;
    @<prepare data for bounding procedure@>@;
    @<bound@>@;
    @<apply primal heuristics@>@;
    @<update current node information@>@;
    if (!fathomed) {
#ifdef BIQMAC
	if ((cputime()-starttime)>3*3600*100)
	    {
		printf("Time limit of three hours exceeded.\n\n");
		printf("Summary of current status:\n");
		printf("Upper bound: %.2f \n",upper_bound);
		break;
	    };
#endif 
      @<branch@>@;
      if (n_bbnodes >= P.max_bb_nodes) {
      printf("maximum number of nodes in the branch and bound tree reached.\n");
      break;}
    }
    else {
    fathom:
      @<fathom current node@>@;
    }
  }
}

@ @<remove from active list@>=
{
  Active_node *tmp;
  upper_bound = active->node->ub;
  if(P.printlev >= 2) {
    printf("Processing node");
    node_info(active->node);
    printf("\n");
  }
  tmp = active;
  cur_node = active->node;
  active = active->succ;
  if (active) active->pred = NULL;
  free((char *) tmp);
  if (!active) last_active = NULL;
  if(upper_bound < best_value + P.gap_relax) goto fathom;
}

@ @<find common ancestor@>=
{
  Tree_node *p_n, *c_n;
  p_n=prev_node;
  c_n=cur_node;
  while (c_n->level > p_n->level)
    c_n = c_n->father;
  while (p_n->level > c_n->level)
    p_n = p_n->father;
  while (p_n != c_n) {
    c_n = c_n->father;
    p_n = p_n->father;
  }
  common_ancestor = p_n;
}

@ @<undo shrinkings@>=
{
  Tree_node *p_n;
  for (p_n=prev_node; p_n != common_ancestor; p_n = p_n->father) {
    while(p_n->log_shrink > 0) {
      log_n--;
      if(log_stack[log_n].flip == -1)
        p_n->correction -= star_value[log_stack[log_n].left];
      if(P.printlev >= 1) { 
	printf("unshrinking [%d,%d] by logical implications, ",
	       log_stack[log_n].flip*(log_stack[log_n].left+1),
	       log_stack[log_n].right+1);
        printf("correction becomes %f\n",p_n->correction);
      }
      unshrink(cost,log_stack[log_n].left,log_stack[log_n].right,
	       log_stack[log_n].flip);
      p_n->log_shrink--;
    }
    unshrink(cost,(p_n->father)->left,(p_n->father)->right,p_n->flip);
  }
}
@    @<do shrinkings@>=
{
  Tree_node *c_n;
  int n_positions, i;
  n_positions = cur_node->level - common_ancestor->level;
  i = n_positions;
  for (c_n = cur_node; c_n != common_ancestor; c_n = c_n->father) {
    stack[--i].left = (c_n->father)->left;
    stack[i].right = (c_n->father)->right;
    stack[i].flip = c_n->flip;
  }
  for (i=0; i< n_positions; i++)
    shrink(cost,stack[i].left,stack[i].right,stack[i].flip);
}

@ @<prepare constraints@>=
{ int i, j, lo1, lo2, l1, l2, pos;
  start:
  pos = 0;
  for(i= 0;i<n;i++) if(cost.alive[i]) new_pos[i] = pos++;
  tmp_con = 0;
  if(num_con) {
    pos = 0;
    for(i= 0;i<num_con;i++) {
      tmp_sense[tmp_con] = sense[i];
      tmp_rhs[tmp_con] = rhs[i];
      tmp_bg[tmp_con] = pos;
      for(j= bg[i];j<bg[i+1];j++) if(!Zero(co[j])) {
	l1= lo1= l2i1(sp[j]);
	l2= lo2= l2i2(sp[j]);
	while(l1 != l2 && (!cost.alive[l1] || !cost.alive[l2])) {
	  if(!cost.alive[l2] && l1==cost.ref[l2]) { 
	    tmp_rhs[tmp_con] -= cost.link[lo2]*co[j];
	    l2= l1;
	  }
	  else if(!cost.alive[l1] && l2==cost.ref[l1]) {
	    tmp_rhs[tmp_con] -= cost.link[lo1]*co[j];
	    l1= l2;
	  }
	  else if((l1 < l2 || cost.alive[l1]) && !cost.alive[l2])
	    l2= cost.ref[l2];
	  else
	    l1= cost.ref[l1];
	}
	if(l1 != l2) {
	  tmp_sp[pos] = e2l(l1,l2);
	  tmp_co[pos] = cost.flip[lo1]*cost.flip[lo2]*co[j];
	  pos++;
	}
      }
      if(pos > tmp_bg[tmp_con]) {
	int all_eq = 1;
	int sp1 = tmp_sp[tmp_bg[tmp_con]];
	double co1 = tmp_co[tmp_bg[tmp_con]];
	for(j= tmp_bg[tmp_con]+1; j < pos; j++) {
	  co1 += tmp_co[j];
	  if(tmp_sp[j] != sp1) all_eq = 0;
	}
	if(all_eq) {
	  if(Zero(co1)) {
	    pos = tmp_bg[tmp_con];
	    if((tmp_sense[tmp_con] >= 0 && Positive(tmp_rhs[tmp_con]))
	       || (tmp_sense[tmp_con] <= 0 && Negative(tmp_rhs[tmp_con]))) {
              if(P.printlev >= 1)
                printf("Constraints are infeasible.\n");
	      goto fathom;
            }
	  }
	  else {
	    pos = tmp_bg[tmp_con];
	    tmp_rhs[tmp_con] /= co1;
	    if(co1 < 0) tmp_sense[tmp_con] = -tmp_sense[tmp_con];
	    if(tmp_rhs[tmp_con] > 1.+EPS) {
	      if(tmp_sense[tmp_con] >= 0) {
                if(P.printlev >= 1)
                  printf("Constraints are infeasible.\n");
                goto fathom;
              }
	    }
	    else if(tmp_rhs[tmp_con] > 1.-EPS) {
	      if(tmp_sense[tmp_con] >= 0) {
		shrink_by_log_imp(cost,l2i1(sp1),l2i2(sp1),1);
		goto start;
	      }
	    }
	    else if(tmp_rhs[tmp_con] > -1.+EPS) {
	      if(tmp_sense[tmp_con] > 0) {
		shrink_by_log_imp(cost,l2i1(sp1),l2i2(sp1),1);
		goto start;
	      }
	      else if(tmp_sense[tmp_con] < 0) {
		shrink_by_log_imp(cost,l2i1(sp1),l2i2(sp1),-1);
		goto start;
	      }
	      else {
                if(P.printlev >= 1)
                  printf("Constraints are infeasible.\n");
		goto fathom;
	      }
	    }
	    else if(tmp_rhs[tmp_con] > -1.-EPS) {
	      if(tmp_sense[tmp_con] <= 0) {
		shrink_by_log_imp(cost,l2i1(sp1),l2i2(sp1),-1);
		goto start;
	      }
	    }
	    else {
	      if(tmp_sense[tmp_con] <= 0) {
                if(P.printlev >= 1)
                  printf("Constraints are infeasible.\n");
                goto fathom;
              }
	    }
	  }
	}
      }
      if(pos == tmp_bg[tmp_con]) {
	if((tmp_sense[tmp_con] <= 0 && Negative(tmp_rhs[tmp_con]))
	   || (tmp_sense[tmp_con] >= 0 && Positive(tmp_rhs[tmp_con]))) {
          if(P.printlev >= 1)
            printf("Constraints are infeasible.\n");
          goto fathom;
        }
      }
      else tmp_con++;
    }
    for(i= 0; i<pos; i++)
      tmp_sp[i]= e2l(new_pos[l2i1(tmp_sp[i])],new_pos[l2i2(tmp_sp[i])]); 
    tmp_bg[tmp_con] = pos;
  } 
}
  
@ @<prepare data for bounding procedure@>=
{
   int i, j, l, k;
   n_supernodes = 1;
   original_name[0] = 0;
   for (i=1, l=1; i<n; i++)
     if (cost.alive[i]) {
       n_supernodes++;
       original_name[l] = i;
       for (j=0, k=0; j<i; j++)
         if (cost.alive[j]) {
           tmp_cost.m[l][k++] = cost.m[i][j];
         }
       l++;
     }
   for (i=0; i<n_supernodes; i++) 
     bound_cut[i] = (double) (best_cut[original_name[i]]*cost.flip[original_name[i]]);
   if(P.printlev >= 5) {
     printf("subproblem:\n");
     printf(" n = %d\n",n_supernodes);
     for(i= 0;i<n;i++) if(cost.alive[i])
     printf(" %d -> %d\n",i+1,new_pos[i]);
     k= 0;
     for(j= 1;j<n_supernodes;j++) for(i= 0;i<j;i++) {
       if(!Zero(tmp_cost.m[1][k]))
       printf(" (%d,%d): %f\n",i,j,tmp_cost.m[1][k]);
       k++;
     }
     printf("e = %d\n",tmp_con);
     for(i= 0;i<tmp_con;i++) con_info(i);
     printf("bound cut:");
     for(i= 0;i<n_supernodes;i++) printf(" %.0f",bound_cut[i]);
     printf("\n");
   }
}

@ @<bound@>=
{
  int tmp, i;
  bound_fix = 2; // if bound does not touch it no fixing will be done
  if(n_supernodes == 1) {
    bound_lb = bound_ub = 0.0;
    bound_cut[0] = -1;
    bound_left = bound_right = 0;
  }
  else bound(P,n_supernodes,tmp_cost.m[1],
	     tmp_con,tmp_sense,tmp_rhs,tmp_bg,tmp_sp,tmp_co,
	     best_value-cur_node->correction,
	     &bound_lb,&bound_ub,bound_cut,
	     &bound_left,&bound_right,&bound_fix);
  bound_left  = original_name[bound_left];
  bound_right = original_name[bound_right];
  if (bound_left > bound_right) {
    tmp = bound_left;
    bound_left  = bound_right;
    bound_right = tmp;
  }
  bound_ub += cur_node->correction;
  bound_lb += cur_node->correction;
  for(i=n_supernodes-1; i >= 0; i--)
    bound_cut[original_name[i]] = bound_cut[i]*cost.flip[original_name[i]];
  for(i=0; i<n; i++) if(!cost.alive[i])
    bound_cut[i] = bound_cut[cost.ref[i]]*cost.link[i];
}

@ @<apply primal heuristics@>=
{
  int i;
  for (i=0; i<n; i++) prim_cut[i] = bound_cut[i];
  p_res = primal(n,prim_cut,orig_cost.m,user_info);
  if (p_res > 0 && valid_cut(n,prim_cut-1,num_con,sense,rhs,bg,sp,co)) {
    compute_orig_cost(prim_lb,prim_cut);
    if(prim_lb > bound_lb) {
      bound_lb = prim_lb;
      for (i=0; i<n; i++) bound_cut[i] = prim_cut[i];
    }
  }
}

@ @<update current node information@>=
   {
     cur_node->ub = bound_ub;
     if((bound_ub > best_value) && (bound_lb > best_value + EPS)) {
       best_value = bound_lb;
       @<update best cut@>@;
       @<prune tree@>@;
     }
     fathomed = (p_res == 2) || (bound_ub < best_value + P.gap_relax);
   }

@  @<update best cut@>=
{
  int i;
  for(i=0; i<n; i++) {
    if(P.printlev >= 3) printf(" %f",bound_cut[i]);
    best_cut[i] = Round(bound_cut[i]);
  }
  if(P.printlev >= 3) printf("\n");
  if(P.printlev >= 2) {
    printf("update best cut -> %f\n",best_value);
    for(i= 0;i<n;i++) printf(" %d",best_cut[i]);
    printf("\n");
  }    
}

@ @<prune tree@>=
{ 
  Active_node *current;
  for (current=last_active; current && (current->node)->ub<=best_value; 
       current=current->pred) {
    kill_node=current->node;
    @<kill tree-node@>@;
    last_active = current->pred;
    if (last_active) {
      last_active->succ = NULL;
    }
    else {
      active = NULL;
    }
  }
}

@ @<branch@>=
{
  Tree_node *tmp_t;
  Active_node *tmp_a[2], *current, *tmp_first, *tmp_last;
  double ub;
  int son, flip;

  tmp_a[0] = tmp_a[1] = (Active_node *) NULL;
  cur_node->left = bound_left;
  cur_node->right = bound_right;
  cur_node->n_sons = 0;
  ub = cur_node->ub;
  for (son = 0; son < 2; son++) {
    if (son == 0) { // setting the edge to value 1
      flip = -1;
      if (bound_fix == 0) continue;
    } else {        // setting the edge to value 0
      flip = 1;
      if (bound_fix == 1) continue;
    }
    cur_node->n_sons++;
    tmp_t = (Tree_node *) malloc (sizeof(Tree_node));
    tmp_t->name = n_bbnodes++;
    tmp_t->flip = flip;
    tmp_t->level = cur_node->level + 1;
    tmp_t->n_sons = 0;
    tmp_t->log_shrink = 0;
    tmp_t->father = cur_node;
    tmp_t->ub = ub;
    tmp_t->correction = cur_node->correction;
    if (flip == -1)
      tmp_t->correction += star_value[bound_left]; 

    tmp_a[son] = (Active_node *) malloc (sizeof(Active_node));
    tmp_a[son]->node = tmp_t;
  }
  if (tmp_a[0]) {
    tmp_a[0]->pred = NULL;
    tmp_a[0]->succ = tmp_a[1];
    tmp_first = tmp_a[0];
  } else {
    tmp_first = tmp_a[1];
  }
  if (tmp_a[1]) {
    tmp_a[1]->pred = tmp_a[0];
    tmp_a[1]->succ = NULL;
    tmp_last = tmp_a[1];
  } else {
    tmp_last = tmp_a[0];
  }
  if (!active) {
    active = tmp_first;
  }
  else {
    current = active;
    while ((current) && ((current->node)->ub > ub)) {
      tmp_first->pred = current;
      current = current->succ;
    }
    tmp_last->succ = current;
    if (current) {
      if (!(current->pred)) {
        active = tmp_first;
      }
      else {
        (current->pred)->succ = tmp_first;
      }
      current->pred=tmp_last;
    }
    else {
      (tmp_first->pred)->succ=tmp_first;
    }
  }
  if (!tmp_last->succ) last_active = tmp_last;
  prev_node = cur_node;
  if(P.printlev >= 2) {
    printf("Adding nodes:");
    if (tmp_a[0]) {
      node_info(tmp_a[0]->node);
      printf("\n             ");
    }
    if (tmp_a[1]) {
      node_info(tmp_a[1]->node);
      printf("\n");
    }
  }
}

@ @<fathom current node@>=
{
 Tree_node *tmp;
 cur_node->n_sons = 0;
 while (cur_node->father && !cur_node->n_sons) {
   if(P.printlev >= 2) {
     printf("Fathoming node");
     node_info(cur_node);
     printf("\n");
   }
   tmp = cur_node->father;
   tmp->n_sons--;
   while(cur_node->log_shrink > 0) {
     log_n--;
     if(P.printlev >= 1) {
       printf("unshrinking [%d,%d] by logical implications\n",
               log_stack[log_n].flip*(log_stack[log_n].left+1),
	       log_stack[log_n].right+1);
     }
     unshrink(cost,log_stack[log_n].left,log_stack[log_n].right,
	      log_stack[log_n].flip);
     cur_node->log_shrink--;
   }
   unshrink(cost,tmp->left,tmp->right,cur_node->flip);
   free((char *) cur_node);
   cur_node = tmp;
 }
 prev_node = cur_node;
}

@ @<kill tree-node@>=
{
 Tree_node *tmp;
 kill_node->n_sons = 0;
 while (kill_node && !kill_node->n_sons) {
   if(P.printlev >= 2) {
     printf("Killing node");
     node_info(kill_node);
     printf("\n");
   }
   tmp = kill_node->father;
   tmp->n_sons--;
   free((char *) kill_node);
   kill_node = tmp;
 }
}

@ @<print matrix@>=
{
  int i, j;
  printf("\n");
  for (i=0; i<n; i++) if (cost.alive[i]) {
    printf("%4d | ",i + 1);
    for (j=0; j < i; j++) if (cost.alive[j])
      printf("%4d ", (int) cost.m[i][j]);
    printf("\n");
  }
  printf("       ");
  for (j=0; j<n; j++) if (cost.alive[j])
    printf("---- ");
  printf("\n       ");
  for (j=0; j<n; j++) if (cost.alive[j])
    printf("%4d ",j + 1);
  printf("\n\n");
}

@ @<free@>=
{
  free(cur_node);
  free(original_name);
  free(star_value);
  free(bound_cut);
  free(prim_cut);
  free(new_pos);
  free(stack);
  free(log_stack);
  if(num_con) {
    free(tmp_sense);
    free(tmp_rhs);
    free(tmp_bg);
    free(tmp_sp);
    free(tmp_co);
  }
  free_mat(cost);
  free_mat(tmp_cost);
  free_mat(orig_cost);
}

@ @<output the result@>=
{
  if(P.printlev >= 1) { 
    printf("Total number of branch-and-bound nodes: %d\n",n_bbnodes);
    printf("Total time: %.2f sec\n",(cputime() - starttime)/100.);
  }
}

@*

@<utility functions@>=

#include <sys/types.h>
#include <sys/times.h>

int cputime()
{ struct tms now;
  int t;
  times(&now);
  t= (int) now.tms_utime;
  return t;
}
@* Index.


