#define SECONDS 50000

#define BIQMAC 0 \

//https://stackoverflow.com/questions/31589392/from-c-malloc-matrix-to-java-array
#define alloc_mat(A,n)  \
{ \
	int i; \
	A.n= n; \
	A.alive= (int*) malloc(sizeof(int) *n) ; \
	A.ref= (int*) malloc(sizeof(int) *n) ; \
	A.flip= (int*) malloc(sizeof(int) *n) ; \
	A.link= (int*) malloc(sizeof(int) *n) ; \
	for(i= 0;i<n;i++) { \
		A.alive[i]= 1; \
		A.flip[i]= 1; \
	} \
	A.m= (double**) malloc(sizeof(double) *n) ; \
	A.m[1]= (double*) calloc((n*(n-1) /2) ,sizeof(double) ) ; \
	for(i= 2;i<n;i++)  \
	A.m[i]= A.m[i-1]+i-1; \
} \

#define free_mat(A)  \
{ \
	free(A.m[1]) ; \
	free(A.m) ; \
	free(A.ref) ; \
	free(A.flip) ; \
	free(A.link) ; \
	free(A.alive) ; \
} \

#define shrink(A,i,j,f)  \
{ \
	int small,large,l; \
	small= i<j?i:j; \
	large= i<j?j:i; \
	 \
	star_value[small]= f*(star_value[small]-A.m[large][small]) + \
	(star_value[large]-A.m[large][small]) ; \
	if(f==-1) { \
		for(l= 0;l<small;l++)  \
			if(A.alive[l]&&(l!=large) )  \
				star_value[l]-= (2*A.m[small][l]) ; \
		for(l= small+1;l<n;l++)  \
			if(A.alive[l]&&(l!=large) )  \
				star_value[l]-= (2*A.m[l][small]) ; \
	} \
	for(l= 0;l<small;l++)  \
		A.m[small][l]= f*A.m[small][l]+A.m[large][l]; \
	for(l= small+1;l<large;l++)  \
		A.m[l][small]= f*A.m[l][small]+A.m[large][l]; \
	for(l= large+1;l<A.n;l++)  \
		A.m[l][small]= f*A.m[l][small]+A.m[l][large]; \
	A.alive[large]= 0; \
	A.ref[large]= small; \
	A.flip[small]*= f; \
	A.link[large]= A.flip[small]; \
} \

#define unshrink(A,i,j,f)  \
{ \
	int small,large,l; \
	small= i<j?i:j; \
	large= i<j?j:i; \
	for(l= 0;l<small;l++)  \
		A.m[small][l]= f*(A.m[small][l]-A.m[large][l]) ; \
	for(l= small+1;l<large;l++)  \
		A.m[l][small]= f*(A.m[l][small]-A.m[large][l]) ; \
	for(l= large+1;l<A.n;l++)  \
		A.m[l][small]= f*(A.m[l][small]-A.m[l][large]) ; \
	 \
	if(f==-1) { \
		for(l= 0;l<small;l++)  \
			if(A.alive[l]&&(l!=large) )  \
				star_value[l]+= (2*A.m[small][l]) ; \
		for(l= small+1;l<n;l++)  \
			if(A.alive[l]&&(l!=large) )  \
			star_value[l]+= (2*A.m[l][small]) ; \
	} \
	star_value[small]= f*(star_value[small]-star_value[large] \
	+A.m[large][small]) +A.m[large][small]; \
	 \
	A.alive[large]= 1; \
	A.ref[large]= -1; \
	A.flip[small]*= f; \
	A.link[large]= 0; \
} \

#define node_info(p,gap)  \
{ \
	if(p) { \
		if(p->father) { \
			printf("# son of [%d,%d] <ub %4.1f> ^corr. %4.1f^ {gap %5.2f} time: %.2f\n",\
			p->flip*((p->father)->left+1) , \
			(p->father)->right+1,p->ub,p->correction, \
			gap,(double)(time(NULL)-start)) ; \
		 \
		 \
		 \
		 \
		 \
		 \
		} \
		else{ \
		 \
			printf(" Root <ub %4.1f>",p->ub) ; \
		} \
	} \
	else{ \
	printf(" No such node!! ") ; \
	} \
} \

#define Getparam(a,b)  \
	if(!(strcmp(st,#a) ) ) fscanf(pf,#b,&(P.a) )  \

#define con_info(i)  \
{ \
	int j; \
	for(j= tmp_bg[i];j<tmp_bg[i+1];j++)  \
		printf(" + (%f)*(%d,%d)",tmp_co[j],l2i1(tmp_sp[j]) ,l2i2(tmp_sp[j]) ) ; \
	if(tmp_sense[i]==-1) printf(" <= ") ; \
	else if(tmp_sense[i]==1) printf(" >= ") ; \
	else printf(" == ") ; \
	printf("%f\n",tmp_rhs[i]) ; \
} \

#define compute_orig_cost(v,c)  \
{ \
	int i,j; \
	v= 0.0; \
	for(i= 0;i<n;i++) { \
		for(j= 0;j<i;j++) { \
			v+= orig_cost.m[i][j]*(double) (c[i]*c[j]<0) ; \
		} \
	} \
} \

#define SEND_TO_WORKERS 1
#define GET_CHILDREN 2
#define FATHOM 3
#define DONE 4
#define NEW_VALUE 5
#define TIME_LIMIT 6

#include <sys/types.h> 
#include <sys/times.h>

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 
#include <string.h> 
#include <gb_flip.h> 
#include <declarations.h> 
#include <stdbool.h>
#include <mpi.h>


typedef struct tree_node{
	int *name;
	int flip;
	int level;
	int n_sons;
	struct tree_node *father;
	struct tree_node *left_ch;
	struct tree_node *right_ch;
	double correction;
	double ub;
	int left,
	right;
}Tree_node;


typedef struct{
	int n;
	int*alive;
	int*ref;
	int*flip;
	int*link;
	double**m;
}Matrix;


void bound(T_param,int,double*,
int,int*,double*,int*,int*,double*,
double,double*,double*,double*,int*,int*,int*);

int primal(int,double*,double**,void*);

//dynamic queue
struct queue{
	int *node_name;
	struct queue *next;
	double ub;
};
struct queue *front = NULL;
struct queue *rear = NULL;
int itemCount = 0;

bool isEmpty()
{
	if(front == NULL && rear == NULL)
		return true;
	else
		return false;
}

struct queue *dequeue()
{
	if(isEmpty())
		return NULL;
	else
	{
		struct queue *y = front;
		front = front->next;
		if(front == NULL)
			rear = NULL;
		y->next = NULL;
		itemCount--;

		return y;
	}
}

void enqueue(int *y, int n, double ub)
{
	struct queue *x = malloc(sizeof(struct queue));
	x->node_name = (int*) malloc(n*sizeof(int));
	x->ub = ub;
	int i;
	for(i = 0; i < n; i++)
		x->node_name[i] = y[i];

	if(isEmpty())
	{
		front = x;
		rear = x;
	}
	else
	{
		struct queue *y = front;
		if(y->next == NULL) //one element in queue
		{
			if(ub > y->ub)
			{
				x->next = y;
				front = x;
			}
			else
			{
				x->next = NULL;
				y->next = x;
				rear = x;
			}
		}
		else //two or more elements in queue
		{
			if(ub > y->ub)
			{
				x->next = y;
				front = x;
			}
			else
			{
				struct queue *z = front->next;
				while(z != NULL && (z->ub > ub))
				{
					z = z->next;
					y = y->next;
				}
				if(z == NULL) //all elements in queue have biger ub than node
				{
					x->next = NULL;
					rear->next = x;
					rear = x;
				}
				else //we insert new element between y and z!
				{
					x->next = z;
					y->next = x;
				}
			}
		}
	}
	itemCount++;
}

int numbWorkers; //number of processes (master + workers)
int rank; //process rank
int tag;

void freeToBusy(int *busyWorkers, int i)
{
	busyWorkers[i] = 1;
}
void busyToFree(int *busyWorkers, int i)
{
	busyWorkers[i] = 0;
}
int findFreeWorker(int *busyWorkers)
{
	int i;
	for(i = 1; i <= numbWorkers-1; i++)
	{
		if(busyWorkers[i] == 0)
		{
			return i;
		}
	}
}

int
bab_sdp(int n,int mm,
int num_con,int*sense,double*rhs,int*bg,double*co,int*sp,
int*best_cut,int prlev,int check_integrality,void*user_info,int argc, char **argv)
{

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &numbWorkers);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	int w;
	for(w = 1; w < argc; w++) //for each argument of function do this
	{
		front = NULL;
		rear = NULL;
		itemCount = 0;

		int m,
		n_supernodes,
		n_bbnodes;
		int*original_name;
		double*star_value;
		Tree_node*cur_node;
		Matrix cost,
		tmp_cost,
		orig_cost;
		double best_value,
		upper_bound,
		bound_ub,
		bound_lb,
		prim_lb;
		int bound_left,
		bound_right,
		bound_fix;
		double*bound_cut;
		double*prim_cut;
		T_param P;
		int starttime;
		int tmp_con;
		int*tmp_sense;
		double*tmp_rhs;
		int*tmp_bg;
		double*tmp_co;
		int*tmp_sp;
		int*new_pos;
		int p_res;

		int i,j,k,l,large,small;
		double c;
		int mult_edges= 0;
		int loops= 0;
		P.gap_relax= 1.-EPS;

		int *head;
		int *tail;
		double *weight;

		if(rank == 0)
			printf("graph: %s\n", argv[w]);
		FILE *mojaDatoteka = fopen(argv[w], "r");

		if(mojaDatoteka)
		{
			fscanf(mojaDatoteka,"%d %d", &n, &mm);
			head = (int*)malloc(sizeof(int)*mm);
			tail = (int*)malloc(sizeof(int)*mm);
			weight = (double*)malloc(sizeof(double)*mm);
			best_cut = (int*)calloc(n, sizeof(int));
			for(i = 0; i < mm; i++)
			{
				fscanf(mojaDatoteka, "%d %d %lf", &head[i], &tail[i], &weight[i]);
				//if(rank == 0)
					//printf("%d %d %lf\n", head[i], tail[i], weight[i]);
			}
		}
		else
		{
			printf("NAPAKA!\n");
		}
		
		m= (n*(n-1))/2;
		alloc_mat(cost,n);
		alloc_mat(orig_cost,n);

		for(i = 0; i < n; i++){
			cost.ref[i] = 0;
			cost.link[i] = 0;
		}

		if(rank == 0)
		{
			for(i = 0; i < mm; i++){
				j = head[i] - 1;
				k = tail[i] - 1;
				c = weight[i];
				if(j != k){
					if(check_integrality && !Integer(c))
						P.gap_relax = EPS;
					large = (j > k ? j : k);
					small = (j > k ? k : j);
					if(!Zero(cost.m[large][small]))
						mult_edges++;
					cost.m[large][small] += c;
					orig_cost.m[large][small] += c;
				}
				else
					loops++;
			}
			if(loops > 0)
				printf("Warning: the graph contains loops. Loops have been ignored.\n");
			if(mult_edges > 0)
				printf("Warning: multiple edges! Edge weights have been summed up.\n");
			#ifdef BIQMAC
				printf("The input graph has %d nodes and %d edges, ", n, mm-mult_edges-loops);
				if(P.gap_relax > EPS+EPS)
				printf("edge weights are integers.\n\n");
				else
				printf("edge weights are reals.\n\n");
			#endif
		
			free(head);
			free(tail);
			free(weight);
		}


		FILE*pf;
		int eof;
		char st[128];

		P.max_bb_nodes= 1000000;
		P.inner_it= 3;
		P.max_inner_it= 10;
		P.outer_it= 5;
		P.extra_it= 0;
		P.min_outer_it= 5;
		P.viol_tol= 1.0e-3;
		P.heap_size= 1000;
		P.mandatory_tr= 1000;
		P.max_gen_tr= 1000;
		P.branch_rule= 3;
		P.printlev= 1;
		pf= fopen("param","r");
		if(pf){
			do{
				eof= (fscanf(pf,"%s",st)==EOF);
				Getparam(max_bb_nodes,%d);
				Getparam(gap_relax,%lf);
				Getparam(inner_it,%d);
				Getparam(max_inner_it,%d);
				Getparam(outer_it,%d);
				Getparam(extra_it,%d);
				Getparam(min_outer_it,%d);
				Getparam(viol_tol,%lf);
				Getparam(heap_size,%d);
				Getparam(mandatory_tr,%d);
				Getparam(max_gen_tr,%d);
				Getparam(printlev,%d);
				Getparam(branch_rule,%d);
			}while(strcmp(st,"end.")&&!eof);
			fclose(pf);
		}
		if(prlev>=0)
			P.printlev= prlev;
		if(P.mandatory_tr> P.heap_size||
		P.mandatory_tr> P.max_gen_tr||
		P.max_gen_tr> P.heap_size){
			printf("parameter inconsistency: mandatory_tr < max_gen_tr < heap_size\n");
			exit(1);
		}
		if(P.printlev>=2){
			printf("max_bb_nodes =    %d\n",P.max_bb_nodes);
			printf("gap_relax =       %lf\n",P.gap_relax);
			printf("inner_it =        %d\n",P.inner_it);
			printf("max_inner_it =    %d\n",P.max_inner_it);
			printf("outer_it =        %d\n",P.outer_it);
			printf("extra_it =        %d\n",P.extra_it);
			printf("min_outer_it =    %d\n",P.min_outer_it);
			printf("viol_tol =        %f\n",P.viol_tol);
			printf("heap_size =       %d\n",P.heap_size);
			printf("mandatory_tr =    %d\n",P.mandatory_tr);
			printf("max_gen_tr =      %d\n",P.max_gen_tr);
			printf("printlev =        %d\n",P.printlev);
			printf("branch_rule =     %d\n",P.branch_rule);
			printf("\n\n");
		}

		double cij;

		//wallclock
		time_t start = time(NULL);
		//maximum number of busy workers
		int maxWorkers = 0;

		original_name= (int*)malloc(sizeof(int)*n);
		star_value= (double*)malloc(sizeof(double)*n);
		double *orig_star_value = (double*)malloc(sizeof(double)*n);
		bound_cut= (double*)calloc(n ,sizeof(double));
		prim_cut= (double*)malloc(sizeof(double)*n);
		alloc_mat(tmp_cost,n);
		new_pos= (int*)malloc(sizeof(int)*n);
		int *busyWorkers;
		int *node_name;

		if(rank == 0)
		{
			busyWorkers = (int*) calloc(numbWorkers, sizeof(int));
			//master is busy
			busyWorkers[0] = 1;
			n_bbnodes = 0;
			best_value = 0.0;
			upper_bound = 0.0;
			
			for(i = 0; i < n; i++){
				star_value[i] = 0.0;
				for(j = 0; j < i; j++){
					cij = cost.m[i][j];
					star_value[i] += cij;
					star_value[j] += cij;
					orig_star_value[i] += cij;
					orig_star_value[j] += cij;
					if(cij > 0)
						upper_bound += cij;
				}
			}
			compute_orig_cost(best_value, best_cut);
			for(i = 0; i < n; i++)
				prim_cut[i]= (double)best_cut[i];
			p_res= primal(n, prim_cut,orig_cost.m, user_info);
			if(p_res> 0 && valid_cut(n, prim_cut-1, num_con, sense, rhs, bg, sp, co)){
				compute_orig_cost(prim_lb, prim_cut);
				if(prim_lb > best_value){
					best_value = prim_lb;
					for(i = 0; i < n; i++)
						best_cut[i]= Round(prim_cut[i]);
				}
			}
			bound_ub = upper_bound;
			//gb_init_rand(1);

			//we create root node
			Tree_node *root = (Tree_node*)malloc(sizeof(Tree_node));
			root->name = (int*)calloc(n, sizeof(int));
			root->flip = 1;
			root->level = 0;
			root->n_sons = 0;
			root->father = NULL;
			root->left_ch = NULL;
			root->right_ch = NULL;
			root->correction = 0.0;
			root->ub= upper_bound;

			n_bbnodes++;

			enqueue(root->name,n,root->ub);

			int counter = 0;
			int numbFreeWorkers = numbWorkers-1;
			int iterations;
			int freeWorker;
			int source;
			struct queue *q;
			int active;
			int over;
			int maxQueueSize = itemCount;
			node_name = (int*)malloc(sizeof(int)*n);

			while(1)
			{
				if(counter == 0)
				{
					tag = SEND_TO_WORKERS;
					counter = 1;
				}
				else
				{
					MPI_Recv(&source, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(&tag, 1, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
				}
				if(tag == GET_CHILDREN)
				{
					//printf("\nTag for children \n");
					double upper_b;
					MPI_Recv(&upper_b, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(node_name, n, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(star_value, n, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(&bound_left, 1, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
					//printf("bound_left = %d\n", bound_left);
					MPI_Recv(&bound_right, 1, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
					//printf("bound_right = %d\n", bound_right);

					double potencial_bound_lb;
					MPI_Recv(&potencial_bound_lb, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
					if(potencial_bound_lb > bound_lb)
						bound_lb = potencial_bound_lb;

					cur_node = root;
					i = 0;
					while(node_name[i] != 0)
					{
						if(node_name[i] == 1){
							cur_node = cur_node->left_ch;
						}
						if(node_name[i] == -1){
							cur_node = cur_node->right_ch;
						}
						i++;
					}

					cur_node->ub = upper_b;
					cur_node->left = bound_left;
					cur_node->right = bound_right;
					cur_node->n_sons += 2;

					Tree_node *ch1 = (Tree_node*)malloc(sizeof(Tree_node));
					ch1->name = (int*)malloc(sizeof(int)*n);

					for(i = 0; i < n; i++)
					{
						ch1->name[i] = node_name[i];
					}
					ch1->name[cur_node->level] = 1;

					ch1->flip = 1;
					ch1->level= cur_node->level+1;
					ch1->n_sons = 0;
					ch1->father = cur_node;
					ch1->left_ch = NULL;
					ch1->right_ch = NULL;
					ch1->correction = cur_node->correction;
					ch1->ub = cur_node->ub;

					cur_node->left_ch = ch1;
					
					enqueue(ch1->name,n,ch1->ub);

					Tree_node *ch2 = (Tree_node*)malloc(sizeof(Tree_node));
					ch2->name = (int*)malloc(sizeof(int)*n);

					for(i = 0; i < n; i++)
					{
						ch2->name[i] = node_name[i];
					}
					ch2->name[cur_node->level] = -1;

					ch2->flip = -1;
					ch2->level = cur_node->level+1;
					ch2->n_sons = 0;
					ch2->father = cur_node;
					ch2->correction = cur_node->correction;
					ch2->correction += star_value[cur_node->left];
					ch2->ub= cur_node->ub;

					cur_node->right_ch = ch2;

					enqueue(ch2->name,n,ch2->ub);

					busyToFree(busyWorkers, source);
					numbFreeWorkers++;

					if(itemCount > maxQueueSize)
						maxQueueSize = itemCount;

					n_bbnodes += 2;
					tag = SEND_TO_WORKERS;
				}
				if(tag == FATHOM)
				{
					//printf("\nTag for fathom \n");
					MPI_Recv(node_name, n, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
					double potencial_bound_lb;
					MPI_Recv(&potencial_bound_lb, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
					if(potencial_bound_lb > bound_lb)
						bound_lb = potencial_bound_lb;

					cur_node = root;

					i = 0;
					while(node_name[i] != 0)
					{
						if(node_name[i] == 1){
							cur_node = cur_node->left_ch;
						}
						if(node_name[i] == -1){
							cur_node = cur_node->right_ch;
						}
						i++;
					}

					Tree_node*tmp;
					cur_node->n_sons= 0;
					while(cur_node->father&&!cur_node->n_sons){
						//printf("Fathoming node");
						//node_info(cur_node,cur_node->ub-best_value);
						tmp= cur_node->father;
						tmp->n_sons--;
						free(cur_node->name);
						free(cur_node);
						cur_node= tmp;
					}

					//we free worker
					busyToFree(busyWorkers, source);
					numbFreeWorkers++;

					tag = SEND_TO_WORKERS;
				}
				if(tag == NEW_VALUE)
				{
					//printf("\nTag for new value \n");
					MPI_Recv(&bound_lb, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
					MPI_Recv(bound_cut, n, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
					if(bound_lb > best_value)
					{
						best_value = bound_lb;
						for(i = 0; i < n; i++)
						{
							best_cut[i]= Round(bound_cut[i]);
						}
					}

					MPI_Send(&best_value, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD);
				}
				if(tag == SEND_TO_WORKERS)
				{
					//printf("\nTag for sending to workers \n");
					over = 0;
					if(itemCount == 0 && numbFreeWorkers == numbWorkers-1)
					{
						printf("BEST VALUE = %lf\n", best_value);
						printf("BEST CUT: ");
						for(i = 0; i < n; i++)
						{
							printf("%d ", best_cut[i]);
						}
						printf("\n");
						printf("ONE SIDE OF THE CUT: ");
						for(i = 0; i < n; i++)
						{
							if(best_cut[i] == 1)
								printf("%d ", i+1);
						}
						printf("\n");
						tag = DONE;
					}
					else
					{
						if(maxWorkers < (numbWorkers-numbFreeWorkers))
							maxWorkers = (numbWorkers-numbFreeWorkers);
						if(maxQueueSize < itemCount)
							maxQueueSize = itemCount;
						iterations = itemCount < numbFreeWorkers ? itemCount : numbFreeWorkers;
						active = 1;
						for(k = 1; k <= iterations; k++)
						{
							//we take element from queue
							q = dequeue();
							freeWorker = findFreeWorker(busyWorkers);
							freeToBusy(busyWorkers, freeWorker);
							numbFreeWorkers--;
							cur_node = root;
							
							//we overwrite cost
							for(i = 1; i < n; i++)
							{
								for(j = 0; j < i; j++)
								{
									cost.m[i][j] = orig_cost.m[i][j];
								}
							}
							for(i = 0; i < n; i++)
							{
								cost.alive[i] = 1;
								cost.flip[i] = 1;
								cost.ref[i] = 0;
								cost.link[i] = 0;
								star_value[i] = orig_star_value[i];
							}
							i = 0;
							while(q->node_name[i] != 0)
							{
								if(q->node_name[i] == 1){
									shrink(cost, cur_node->left, cur_node->right, cur_node->left_ch->flip);
									cur_node = cur_node->left_ch;
								}
								if(q->node_name[i] == -1){
									shrink(cost, cur_node->left, cur_node->right, cur_node->right_ch->flip)
									cur_node = cur_node->right_ch;
								}
								i++;
							}

							MPI_Send(&over, 1, MPI_INT, freeWorker, 0, MPI_COMM_WORLD);
							MPI_Send(&active, 1, MPI_INT, freeWorker, 1, MPI_COMM_WORLD);
							MPI_Send(cur_node->name, n, MPI_INT, freeWorker, 2, MPI_COMM_WORLD);
							MPI_Send(&(cur_node->flip), 1, MPI_INT, freeWorker, 3, MPI_COMM_WORLD);
							MPI_Send(&(cur_node->level), 1, MPI_INT, freeWorker, 4, MPI_COMM_WORLD);
							MPI_Send(&(cur_node->n_sons), 1, MPI_INT, freeWorker, 5, MPI_COMM_WORLD);
							MPI_Send(&(cur_node->correction), 1, MPI_DOUBLE, freeWorker, 6, MPI_COMM_WORLD);
							MPI_Send(&(cur_node->ub), 1, MPI_DOUBLE, freeWorker, 7, MPI_COMM_WORLD);
							MPI_Send(star_value, n, MPI_DOUBLE, freeWorker, 8, MPI_COMM_WORLD);
							MPI_Send(cost.alive, n, MPI_INT, freeWorker, 9, MPI_COMM_WORLD);
							MPI_Send(cost.flip, n, MPI_INT, freeWorker, 10, MPI_COMM_WORLD);
							MPI_Send(cost.ref, n, MPI_INT, freeWorker, 11, MPI_COMM_WORLD);
							MPI_Send(cost.link, n, MPI_INT, freeWorker, 12, MPI_COMM_WORLD);
							MPI_Send(&(cost.m[1][0]), (n*(n-1)/2), MPI_DOUBLE, freeWorker, 13, MPI_COMM_WORLD);
							MPI_Send(&best_value, 1, MPI_DOUBLE, freeWorker, 14, MPI_COMM_WORLD);
							MPI_Send(best_cut, n, MPI_INT, freeWorker, 15, MPI_COMM_WORLD);
							MPI_Send(&bound_lb, 1, MPI_DOUBLE, freeWorker, 16, MPI_COMM_WORLD);
							MPI_Send(&bound_ub, 1, MPI_DOUBLE, freeWorker, 17, MPI_COMM_WORLD);
						}
					}
				}
				if(tag == DONE)
				{
					over = 1;
					for(i = 1; i <= numbWorkers-1; i++)
					{
						MPI_Send(&over, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
					}
					printf("B&B nodes = %d\n", n_bbnodes);
					printf("maxQueueSize = %d\n", maxQueueSize);
					printf("maxWorkers = %d\n", maxWorkers);
					time_t end = time(NULL);
					printf("TIME: %.2f sec\n\n\n", (double)(end - start));
					break;
				}
				if(tag == TIME_LIMIT)
				{
					over = 1;
					for(i = 1; i <= numbWorkers-1; i++)
					{
						MPI_Send(&over, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
					}

					struct queue *p;
					Tree_node* tmp1;
					printf("Time limit %d sec exceeded.\n", SECONDS);
					printf("Fathoming nodes...\n");
					while(itemCount != 0)
					{
						p = dequeue();
						tmp1 = root;

						i = 0;
						while(p->node_name[i] != 0)
						{
							if(p->node_name[i] == 1){
								tmp1 = tmp1->left_ch;
							}
							if(p->node_name[i] == -1){
								tmp1 = tmp1->right_ch;
							}
							i++;
						}

						Tree_node*tmp;
						tmp1->n_sons= 0;
						while(tmp1->father&&!tmp1->n_sons){
							//node_info(tmp1,tmp1->ub-best_value);
							tmp= tmp1->father;
							tmp->n_sons--;
							free(tmp1->name);
							free(tmp1);
							tmp1= tmp;
						}
						free(p->node_name);
						free(p);
					}
					printf("best find value = %lf\n", best_value);
					printf("B&B nodes = %d\n", n_bbnodes);
					printf("maxQueueSize = %d\n", maxQueueSize);
					printf("maxWorkers = %d\n", maxWorkers);
					printf("TIME: %.2f sec\n\n\n", (double)(time(NULL) - start));
					break;
				}
			}

		}
		else
		{
			int over;
			int active;
			cur_node = (Tree_node*)malloc(sizeof(Tree_node));
			cur_node->name = (int*)malloc(sizeof(int)*n);

			while(1)
			{
				MPI_Recv(&over, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
				if(over == 1)
				{
					break;
				}
				else
				{
					MPI_Recv(&active, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
					if(active == 1)
					{
						MPI_Recv(cur_node->name, n, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
						MPI_Recv(&(cur_node->flip), 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
						MPI_Recv(&(cur_node->level), 1, MPI_INT, 0, 4, MPI_COMM_WORLD, &status);
						MPI_Recv(&(cur_node->n_sons), 1, MPI_INT, 0, 5, MPI_COMM_WORLD, &status);
						MPI_Recv(&(cur_node->correction), 1, MPI_DOUBLE, 0, 6, MPI_COMM_WORLD, &status);
						MPI_Recv(&(cur_node->ub), 1, MPI_DOUBLE, 0, 7, MPI_COMM_WORLD, &status);
						MPI_Recv(star_value, n, MPI_DOUBLE, 0, 8, MPI_COMM_WORLD, &status);
						MPI_Recv(cost.alive, n, MPI_INT, 0, 9, MPI_COMM_WORLD, &status);
						MPI_Recv(cost.flip, n, MPI_INT, 0, 10, MPI_COMM_WORLD, &status);
						MPI_Recv(cost.ref, n, MPI_INT, 0, 11, MPI_COMM_WORLD, &status);
						MPI_Recv(cost.link, n, MPI_INT, 0, 12, MPI_COMM_WORLD, &status);
						MPI_Recv(&(cost.m[1][0]), (n*(n-1)/2), MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, &status);
						MPI_Recv(&best_value, 1, MPI_DOUBLE, 0, 14, MPI_COMM_WORLD, &status);
						MPI_Recv(best_cut, n, MPI_INT, 0, 15, MPI_COMM_WORLD, &status);
						MPI_Recv(&bound_lb, 1, MPI_DOUBLE, 0, 16, MPI_COMM_WORLD, &status);
						MPI_Recv(&bound_ub, 1, MPI_DOUBLE, 0, 17, MPI_COMM_WORLD, &status);

						upper_bound = cur_node->ub;
						//printf("upper_bound = %lf\n", upper_bound);
						//printf("best_value + gap = %lf\n", best_value + P.gap_relax);
						if(upper_bound < best_value + P.gap_relax)
						{
							tag = FATHOM;
							MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
							MPI_Send(&tag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
							MPI_Send(cur_node->name, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
							MPI_Send(&bound_lb, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
						}
						else
						{
							int pos = 0;
							for(i = 0; i < n; i++)
							{
								if(cost.alive[i])
								{
									new_pos[i] = pos++;
								}
							}
							n_supernodes = 1;
							original_name[0] = 0;
							for(i = 1, l = 1; i < n; i++)
							{
								if(cost.alive[i])
								{
									n_supernodes++;
									original_name[l] = i;
									for(j = 0, k = 0; j < i; j++)
									{
										if(cost.alive[j])
											tmp_cost.m[l][k++] = cost.m[i][j];
									}
									l++;
								}
							}
							for(i = 0; i < n_supernodes; i++)
							{
								bound_cut[i]= (double)(best_cut[original_name[i]]*cost.flip[original_name[i]]);
							}
							bound_fix = 2;
							if(n_supernodes==1)
							{
								bound_lb = bound_ub= 0.0;
								bound_cut[0]= -1;
								bound_left = bound_right= 0;
							}
							else
							{
								tmp_con = 0;
								bound_ub = cur_node->ub;
								bound(P,n_supernodes,tmp_cost.m[1],
									tmp_con,tmp_sense,tmp_rhs,tmp_bg,tmp_sp,tmp_co,
									best_value-cur_node->correction,
									&bound_lb,&bound_ub,bound_cut,
									&bound_left,&bound_right,&bound_fix);
							}
							bound_left = original_name[bound_left];
							bound_right= original_name[bound_right];

							if(bound_left> bound_right){
								int tmp= bound_left;
								bound_left= bound_right;
								bound_right= tmp;
							}

							bound_ub+= cur_node->correction;
							bound_lb+= cur_node->correction;

							for(i= n_supernodes-1;i>=0;i--){
								bound_cut[original_name[i]]= bound_cut[i]*cost.flip[original_name[i]];
							}
							for(i= 0;i<n;i++)
								if(!cost.alive[i]){
									bound_cut[i]= bound_cut[cost.ref[i]]*cost.link[i];
								}

							for(i = 0; i < n; i++)
								prim_cut[i] = bound_cut[i];
							p_res = primal(n,prim_cut,orig_cost.m,user_info);
							if(p_res> 0 && valid_cut(n,prim_cut-1,num_con,sense,rhs,bg,sp,co)){
								compute_orig_cost(prim_lb, prim_cut);
								if(prim_lb > bound_lb){
									bound_lb = prim_lb;
									for(i = 0; i < n; i++)
										bound_cut[i]= prim_cut[i];
								}
							}

							cur_node->ub = bound_ub;
							if((bound_ub> best_value)&&(bound_lb> best_value+EPS)){
								tag = NEW_VALUE;
								MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
								MPI_Send(&tag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

								MPI_Send(&bound_lb, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
								MPI_Send(bound_cut, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

								MPI_Recv(&best_value, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
					
							}
							//printf("bound_ub = %lf\n", bound_ub);
							//printf("best_value+P.gap_relax = %lf\n", best_value+P.gap_relax);
							int fathomed = (p_res== 2)||(bound_ub < best_value+P.gap_relax);
							if(!fathomed){
								if((double)(time(NULL)-start) > SECONDS)
								{									
									tag = TIME_LIMIT;
									MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
									MPI_Send(&tag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
								}
								else
								{		
									tag = GET_CHILDREN;
									MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
									MPI_Send(&tag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
									MPI_Send(&bound_ub, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
									MPI_Send(cur_node->name, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
									MPI_Send(star_value, n, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
									MPI_Send(&bound_left, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
									MPI_Send(&bound_right, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
									MPI_Send(&bound_lb, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
								}					
							}
							else
							{
								tag = FATHOM;
								MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
								MPI_Send(&tag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
								MPI_Send(cur_node->name, n, MPI_INT, 0, 0, MPI_COMM_WORLD);
								MPI_Send(&bound_lb, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
							}
						}
					}
				}
			}
		}
		
		if(rank != 0)
		{
			free(cur_node->name);
			free(cur_node);
		}
		if(rank == 0)
		{
			free(busyWorkers);
			free(node_name);
		}
		free(original_name);
		free(star_value);
		free(orig_star_value);
		free(bound_cut);
		free(prim_cut);
		free(new_pos);
		//if(num_con){
			//free(tmp_sense);
			//free(tmp_rhs);
			//free(tmp_bg);
			//free(tmp_sp);
			//free(tmp_co);
		//}
		free_mat(cost);
		free_mat(tmp_cost);
		free_mat(orig_cost);
	}

	MPI_Finalize();

	return 0;
}
