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

#define SEND_TO_WORKERS 1
#define GET_CHILDREN 2
#define FATHOM 3
#define DONE 4
#define NEW_VALUE 5
#define TIME_LIMIT 6

// time limit
#define SECONDS 85000

// uncomment to print additional information about the graph
//#define BIQMAC \

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

#define alloc_Tree_node(node,n)	\
{ \
	alloc_mat(node->cost,n); \
	node->star_value = (double*) malloc(sizeof(double)*n); \
} \

#define free_Tree_node(node)  \
{ \
	free_mat(node->cost); \
	free(node->star_value); \
} \

#define shrink(A,i,j,f,star_value) \
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

#define unshrink(A,i,j,f,star_value) \
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





typedef struct matrix{
	int n;
	int*alive;
	int*ref;
	int*flip;
	int*link;
	double**m;
}Matrix;

typedef struct tree_node{
	int flip;
	Matrix cost;
	double *star_value;
	double correction;
	double ub;
	int left,
	right;
}Tree_node;

void bound(T_param,int,double*,int,int*,double*,int*,int*,double*,double,double*,double*,double*,int*,int*,int*);
int primal(int,double*,double**,void*);

//dynamic queue:
struct queue{
	Tree_node *node;
	struct queue *next;
};

struct queue *front = NULL;
struct queue *rear = NULL;
int itemCount = 0; //number of items in queue

bool isEmpty()
{
	if((front == NULL && rear == NULL)||(front == rear->next))
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
		struct queue *x = front;
		front = front->next;
		if(front == NULL)
			rear = NULL;
		itemCount--;
		x->next = NULL;

		return x;
	}
}

void enqueue(Tree_node *node)
{
	struct queue *x = malloc(sizeof(struct queue));
	x->node = node;
	x->next = NULL;
	if(rear == NULL)
	{
		front = x;
		rear = x;
	}
	else
	{
		if(front == rear->next)//when one element in queue and we delete it
		{
			rear->next = x;
			front = x;
			rear = rear->next;
		}
		else
		{
			rear->next = x;
			rear = rear->next;
		}
	}
	itemCount++;
}
//end of queue

int numbWorkers; //number of nodes (processes)
int rank; //rank of node (process)
int tag; //for saving the tag number

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

bool compare(int *a, int *b, int n)
{
	int i;
	for(i = 0; i < n; i++)
	{
		if(a[i] != b[i])
			return false;
	}
	return true; 
	
}


// print solution: optimal or approximate
void print_solution(int n, int m, int sol[], double val, int counterBaBNodes, double timeTotal, double density, int optimal) {
	
	printf("{\n");

	printf("\t\"GraphMD\": {\n");
        printf("\t\t\"ObjectTypeString\": \"Graph\",\n");
        printf("\t\t\"Vertices\": %d,\n", n);
        printf("\t\t\"Edges\": %d,\n", m);
        printf("\t\t\"Density\": %.3lf\n", density);

        printf("\t},\n");


	printf("\t\"ExecutionMD\": {\n");
	printf("\t\t\"ExecutionTime\": %.3lf,\n", timeTotal);
	
	if (optimal)
		printf("\t\t\"SolutionType\": \"Optimal\",\n");
	else
		printf("\t\t\"SolutionType\": \"Approximate\",\n");
	
	printf("\t\t\"Solution\": %lf,\n", val);
	
	int index = 1;
	printf("\t\t\"OneSideOfTheCut\": \"(");
	
	int i, j;
	
	int size_of_cut = 0;			// get size of one side of a cut
	
	for (j = 0; j < n; ++j){
		if ( sol[j] == 1 ){
			++size_of_cut;
		}	
	}
	
	for(i = 0; i < n; i++){
        if ( (sol[i] == 1) ){
        	if (index != size_of_cut){
        		printf("%d, ", i+1);
            	++index;
            }
            else
            	printf("%d)\",\n", i+1);
        }    
    }
   
   	printf("\t\t\"BabNodes\": %d\n", counterBaBNodes);
    	printf("\t}\n");
	
	printf("}\n");
	
}

void outputFile(char *argv[], int n, int m, int sol[], double val, int counterBaBNodes, double timeTotal, double density, int optimal){
	
	char fileOutput[200];
	
	if (optimal)
		snprintf(fileOutput, sizeof(fileOutput), "%s%s", argv[2], "_solution_optimal.txt");
	else
		snprintf(fileOutput, sizeof(fileOutput), "%s%s", argv[2], "_solution_approximate.txt");
	
	FILE *file = fopen(fileOutput, "w");
	
	fprintf(file, "{\n");


	fprintf(file, "\t\"GraphMD\": {\n");
        fprintf(file, "\t\t\"ObjectTypeString\": \"Graph\",\n");
        fprintf(file,"\t\t\"Vertices\": %d,\n", n);
        fprintf(file, "\t\t\"Edges\": %d,\n", m);
        fprintf(file, "\t\t\"Density\": %.3lf\n", density);

        fprintf(file, "\t},\n");

	fprintf(file, "\t\"ExecutionMD\": {\n");
	fprintf(file, "\t\t\"ExecutionTime\": %.3lf,\n", timeTotal);
	
	if (optimal)
		fprintf(file, "\t\t\"SolutionType\": \"Optimal\",\n");
	else
		fprintf(file, "\t\t\"SolutionType\": \"Approximate\",\n");
	
	fprintf(file, "\t\t\"Solution\": %lf,\n", val);
	
	int index = 1;
	fprintf(file, "\t\t\"OneSideOfTheCut\": \"(");
	
	int i, j;
	
	int size_of_cut = 0;			// get size of one side of a cut
	
	for (j = 0; j < n; ++j){
		if ( sol[j] == 1 ){
			++size_of_cut;
		}	
	}
	
	for(i = 0; i < n; i++){
        if ( (sol[i] == 1) ){
        	if (index != size_of_cut){
        		fprintf(file, "%d, ", i+1);
            	++index;
            }
            else
            	fprintf(file, "%d)\",\n", i+1);
        }    
    }
   
   	fprintf(file, "\t\t\"BabNodes\": %d\n", counterBaBNodes);
    	fprintf(file, "\t}\n");
	fprintf(file, "}\n");
	
	fclose(file);
}

void outputFile_after_start(char *argv[], int n, int m, int sol[], double val,  int counterBaBNodes, double timeTotal, double density, int begin){

	char fileOutput[200];	

	if (begin){
	
		snprintf(fileOutput, sizeof(fileOutput), "%s%s", argv[2], "_start_info.txt");
		
		FILE *file = fopen(fileOutput, "w");
		
		fprintf(file, "{\n");

		fprintf(file, "\t\"GraphMD\": {\n");
                fprintf(file, "\t\t\"ObjectTypeString\": \"Graph\",\n");
                fprintf(file,"\t\t\"Vertices\": %d,\n", n);
                fprintf(file, "\t\t\"Edges\": %d,\n", m);
                fprintf(file, "\t\t\"Density\": %.3lf\n", density);

                fprintf(file, "\t},\n");

		fprintf(file, "\t\"ExecutionMD\": {\n");
		fprintf(file, "\t\t\"ExecutionTime\": %.3lf,\n", timeTotal);
		
		fprintf(file, "\t\t\"SolutionType\": \"Approximate\",\n");
		
		fprintf(file, "\t\t\"Solution\": %lf,\n", val);
		
		int index = 1;
		fprintf(file, "\t\t\"OneSideOfTheCut\": \"(");
		
		int i, j;
		
		int size_of_cut = 0;			// get size of one side of a cut
		
		for (j = 0; j < n; ++j){
			if ( sol[j] == 1 ){
				++size_of_cut;
			}	
		}
		
		for(i = 0; i < n; i++){
	        if ( (sol[i] == 1) ){
	        	if (index != size_of_cut){
	        		fprintf(file, "%d, ", i+1);
	            	++index;
	            }
	            else
	            	fprintf(file, "%d)\",\n", i+1);
	        }    
	    }
	   
	   	fprintf(file, "\t\t\"BabNodes\": %d\n", counterBaBNodes);
	    	fprintf(file, "\t}\n");

		fprintf(file, "}\n");
		
		fclose(file);

	}
	else {
		static int count = 1;
		
		snprintf(fileOutput, sizeof(fileOutput), "%s%s%d%s", argv[2], "_temporary_solution_", count, ".txt");
		
		FILE *file = fopen(fileOutput, "w");
		
		fprintf(file, "{\n");

                fprintf(file, "\t\"GraphMD\": {\n");
                fprintf(file, "\t\t\"ObjectTypeString\": \"Graph\",\n");
                fprintf(file,"\t\t\"Vertices\": %d,\n", n);
                fprintf(file, "\t\t\"Edges\": %d,\n", m);
                fprintf(file, "\t\t\"Density\": %.3lf\n", density);

                fprintf(file, "\t},\n");

                fprintf(file, "\t\"ExecutionMD\": {\n");
                fprintf(file, "\t\t\"ExecutionTime\": %.3lf,\n", timeTotal);

                fprintf(file, "\t\t\"SolutionType\": \"Approximate\",\n");

                fprintf(file, "\t\t\"Solution\": %lf,\n", val);

                int index = 1;
                fprintf(file, "\t\t\"OneSideOfTheCut\": \"(");

                int i, j;

                int size_of_cut = 0;                    // get size of one side of a cut

                for (j = 0; j < n; ++j){
                        if ( sol[j] == 1 ){
                                ++size_of_cut;
                        }
                }

                for(i = 0; i < n; i++){
                if ( (sol[i] == 1) ){
                        if (index != size_of_cut){
                                fprintf(file, "%d, ", i+1);
                        ++index;
                    }
                    else
                        fprintf(file, "%d)\",\n", i+1);
                }
            }

                fprintf(file, "\t\t\"BabNodes\": %d\n", counterBaBNodes);
                fprintf(file, "\t}\n");

                fprintf(file, "}\n");

		fclose(file);
		count++;
	}
}

/*********** MAIN ROUTINE *********/
int bab_sdp(int n,int mm,int num_con,int*sense,double*rhs,int*bg,double*co,int*sp,int*best_cut,int prlev,int check_integrality,void*user_info,int argc, char **argv){

	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &numbWorkers);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;

	if (argc < 5){
    	if (rank == 0)
    		printf("Not enough arguments to the program.\n");
    		
    	MPI_Finalize(); 
    	return EXIT_SUCCESS;	
    }

	front = NULL;
	rear = NULL;
	itemCount = 0;

	int m,
	n_supernodes,
	n_bbnodes;
	int*original_name;
	Matrix tmp_cost,
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

	#ifdef BIQMAC
	if(rank == 0)
		printf("graph: %s\n", argv[1]);
	#endif
		
	FILE *myFile = fopen(argv[1], "r");

	if(myFile)
	{
		fscanf(myFile,"%d %d", &n, &mm);
		head = (int*)malloc(sizeof(int)*mm);
		tail = (int*)malloc(sizeof(int)*mm);
		weight = (double*)malloc(sizeof(double)*mm);
		best_cut = (int*)calloc(n, sizeof(int));
		for(i = 0; i < mm; i++)
		{
			fscanf(myFile, "%d %d %lf", &head[i], &tail[i], &weight[i]);
		}
	}
	else
	{
		if (rank == 0)
			printf("Problem opening the file!\n");
			
		MPI_Finalize();
		return EXIT_FAILURE;
	}
	
	m= (n*(n-1))/2;

	alloc_mat(orig_cost,n);

	//we prepare root node
	Tree_node *root = (Tree_node*)malloc(sizeof(Tree_node));

	if(rank == 0)
	{
		//we allocate root node and fill it with some data
		alloc_Tree_node(root,n);
		root->flip = 1;
		root->correction = 0.0;

		for(i = 0; i < mm; i++){
			j = head[i] - 1;
			k = tail[i] - 1;
			c = weight[i];
			if(j != k){
				if(check_integrality && !Integer(c))
					P.gap_relax = EPS;
				large = (j > k ? j : k);
				small = (j > k ? k : j);
				if(!Zero(root->cost.m[large][small]))
					mult_edges++;
				root->cost.m[large][small] += c;
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
	}

	free(head);
	free(tail);
	free(weight);

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

	// wallclock
	double start = MPI_Wtime();
	
	// time constants for temp output during the algorithm
	int TIME_INFO = strtol(argv[3], NULL, 10); 		// output first info after MINUTES minutes
	int TIME_TEMP = strtol(argv[4], NULL, 10);			// output further info every hour
    int add_limit = TIME_TEMP;
	
	//max busy workers in algorithm
	int maxWorkers = 0;

	original_name= (int*)malloc(sizeof(int)*n);
	bound_cut= (double*)calloc(n ,sizeof(double));
	prim_cut= (double*)malloc(sizeof(double)*n);
	alloc_mat(tmp_cost,n);
	Tree_node *cur_node = (Tree_node*)malloc(sizeof(Tree_node));
	alloc_Tree_node(cur_node,n);
	new_pos= (int*)malloc(sizeof(int)*n);
	int *busyWorkers;

	// density of the graph
	double density = (double)2*mm / (n*(n-1));

	if(rank == 0)
	{
		busyWorkers = (int*) calloc(numbWorkers, sizeof(int));
		//master is busy
		busyWorkers[0] = 1;
		n_bbnodes = 0;
		best_value = 0.0;
		upper_bound = 0.0;
		
		for(i = 0; i < n; i++){
			root->star_value[i] = 0.0;
			for(j = 0; j < i; j++){
				cij = root->cost.m[i][j];
				root->star_value[i] += cij;
				root->star_value[j] += cij;
				if(cij > 0)
					upper_bound += cij;
			}
		}
		root->ub= upper_bound;

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
		n_bbnodes++;

		//we add root to queue
		enqueue(root);

		int firstTime = 0;
		int numbFreeWorkers = numbWorkers-1;
		int iterations;
		int freeWorker;
		int source;
		int active;
		int over;
		int maxQueueSize = itemCount;

		int print_initial_info = 0;
		while(1)
		{
			if(firstTime == 0)
			{
				tag = SEND_TO_WORKERS;
				firstTime = 1;
			}
			else
			{
				// Check if time for outputing info after start
	            if (!print_initial_info){
		            if ((MPI_Wtime() - start) > TIME_INFO){
		            	outputFile_after_start(argv, n, mm, best_cut, best_value, n_bbnodes, MPI_Wtime() - start, density, 1);
		            	print_initial_info = 1;                 
		            }
	            }
	            
	            // Check if time for outputing info during the algorithm
	            if ((MPI_Wtime() - start) > TIME_TEMP){
	            	outputFile_after_start(argv, n, mm, best_cut, best_value, n_bbnodes, MPI_Wtime() - start, density, 0);
	            	TIME_TEMP += add_limit;                 
	            }
			
				MPI_Recv(&source, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(&tag, 1, MPI_INT, source, 0, MPI_COMM_WORLD, &status);
			}
			if(tag == GET_CHILDREN)
			{
				//printf("\nTag for getting children\n");
				MPI_Recv(&(cur_node->flip), 1, MPI_INT, source, 2, MPI_COMM_WORLD, &status);
				MPI_Recv(&(cur_node->correction), 1, MPI_DOUBLE, source, 3, MPI_COMM_WORLD, &status);
				MPI_Recv(&(cur_node->ub), 1, MPI_DOUBLE, source, 4, MPI_COMM_WORLD, &status);
				MPI_Recv(cur_node->star_value, n, MPI_DOUBLE, source, 5, MPI_COMM_WORLD, &status);
				MPI_Recv(cur_node->cost.alive, n, MPI_INT, source, 6, MPI_COMM_WORLD, &status);
				MPI_Recv(cur_node->cost.flip, n, MPI_INT, source, 7, MPI_COMM_WORLD, &status);
				MPI_Recv(cur_node->cost.ref, n, MPI_INT, source, 8, MPI_COMM_WORLD, &status);
				MPI_Recv(cur_node->cost.link, n, MPI_INT, source, 9, MPI_COMM_WORLD, &status);
				MPI_Recv(&(cur_node->cost.m[1][0]), (n*(n-1)/2), MPI_DOUBLE, source, 10, MPI_COMM_WORLD, &status);
				MPI_Recv(&(cur_node->left), 1, MPI_INT, source, 11, MPI_COMM_WORLD, &status);
				//printf("bound_left = %d\n", cur_node->left);
				MPI_Recv(&(cur_node->right), 1, MPI_INT, source, 12, MPI_COMM_WORLD, &status);
				//printf("bound_right = %d\n", cur_node->right);

				double potencial_bound_lb;
				MPI_Recv(&potencial_bound_lb, 1, MPI_DOUBLE, source, 13, MPI_COMM_WORLD, &status);
				//intf("potencial_bound_lb = %lf\n", potencial_bound_lb);
				//intf("bound_lb = %lf\n", bound_lb);
				if(potencial_bound_lb > bound_lb)
					bound_lb = potencial_bound_lb;


				Tree_node *child1 = (Tree_node*)malloc(sizeof(Tree_node));
				alloc_Tree_node(child1, n);

				child1->flip = 1;
				child1->correction = cur_node->correction;
				child1->ub = cur_node->ub;
				shrink(cur_node->cost, cur_node->left, cur_node->right, 1, cur_node->star_value);
				for(i = 0; i < n; i++)
				{
					child1->cost.alive[i] = cur_node->cost.alive[i];
					child1->cost.ref[i] = cur_node->cost.ref[i];
					child1->cost.link[i] = cur_node->cost.link[i];
					child1->cost.flip[i] = cur_node->cost.flip[i];
					child1->star_value[i] = cur_node->star_value[i];
				}
				for(i = 1; i < n; i++)
				{
					for(j = 0; j < i; j++)
					{
						child1->cost.m[i][j] = cur_node->cost.m[i][j];
					}
				}
		
				enqueue(child1);

				unshrink(cur_node->cost, cur_node->left, cur_node->right, 1, cur_node->star_value);			

				Tree_node *child2 = (Tree_node*)malloc(sizeof(Tree_node));
				alloc_Tree_node(child2, n);

				child2->flip = -1;
				child2->correction = cur_node->correction;
				child2->correction += cur_node->star_value[cur_node->left];

				child2->ub = cur_node->ub;
				shrink(cur_node->cost, cur_node->left, cur_node->right, -1, cur_node->star_value);
				for(i = 0; i < n; i++)
				{
					child2->cost.alive[i] = cur_node->cost.alive[i];
					child2->cost.ref[i] = cur_node->cost.ref[i];
					child2->cost.link[i] = cur_node->cost.link[i];
					child2->cost.flip[i] = cur_node->cost.flip[i];
					child2->star_value[i] = cur_node->star_value[i];
				}
				for(i = 1; i < n; i++)
				{
					for(j = 0; j < i; j++)
					{
						child2->cost.m[i][j] = cur_node->cost.m[i][j];
					}
				}

				enqueue(child2);

				busyToFree(busyWorkers, source);
				numbFreeWorkers++;

				if(itemCount > maxQueueSize)
					maxQueueSize = itemCount;

				n_bbnodes += 2;
				tag = SEND_TO_WORKERS;
			}
			if(tag == FATHOM)
			{
				//printf("\nZahtevek za fathom \n");
				double potencial_bound_lb;
				MPI_Recv(&potencial_bound_lb, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
				if(potencial_bound_lb > bound_lb)
					bound_lb = potencial_bound_lb;

				busyToFree(busyWorkers, source);
				numbFreeWorkers++;

				tag = SEND_TO_WORKERS;
			}
			if(tag == NEW_VALUE)
			{
				//printf("\nTag for new value \n");
				MPI_Recv(&bound_lb, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
				MPI_Recv(bound_cut, n, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
				//printf("bound_lb = %lf\n",bound_lb);
				//printf("best_value = %lf\n", best_value);
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
					// print optimal solution
					print_solution(n, mm, best_cut, best_value, n_bbnodes, MPI_Wtime() - start, density,1);
					outputFile(argv, n, mm, best_cut, best_value, n_bbnodes, MPI_Wtime() - start, density, 1);
				
					tag = DONE;
				}
				else
				{
					if(maxWorkers < (numbWorkers-numbFreeWorkers))
						maxWorkers = (numbWorkers-numbFreeWorkers);
					iterations = itemCount < numbFreeWorkers ? itemCount : numbFreeWorkers;
					active = 1;
					for(k = 1; k <= iterations; k++)
					{
						struct queue *elt;
						Tree_node *cur_node1;
						elt = dequeue();
						cur_node1 = elt->node;

						freeWorker = findFreeWorker(busyWorkers);
						freeToBusy(busyWorkers, freeWorker);
						numbFreeWorkers--;

						MPI_Send(&over, 1, MPI_INT, freeWorker, 0, MPI_COMM_WORLD);
						MPI_Send(&active, 1, MPI_INT, freeWorker, 1, MPI_COMM_WORLD);
						MPI_Send(&(cur_node1->flip), 1, MPI_INT, freeWorker, 2, MPI_COMM_WORLD);
						MPI_Send(&(cur_node1->correction), 1, MPI_DOUBLE, freeWorker, 3, MPI_COMM_WORLD);
						MPI_Send(&(cur_node1->ub), 1, MPI_DOUBLE, freeWorker, 4, MPI_COMM_WORLD);
						MPI_Send(cur_node1->star_value, n, MPI_DOUBLE, freeWorker, 5, MPI_COMM_WORLD);
						MPI_Send(cur_node1->cost.alive, n, MPI_INT, freeWorker, 6, MPI_COMM_WORLD);
						MPI_Send(cur_node1->cost.flip, n, MPI_INT, freeWorker, 7, MPI_COMM_WORLD);
						MPI_Send(cur_node1->cost.ref, n, MPI_INT, freeWorker, 8, MPI_COMM_WORLD);
						MPI_Send(cur_node1->cost.link, n, MPI_INT, freeWorker, 9, MPI_COMM_WORLD);
						MPI_Send(&(cur_node1->cost.m[1][0]), (n*(n-1)/2), MPI_DOUBLE, freeWorker, 10, MPI_COMM_WORLD);
						MPI_Send(&best_value, 1, MPI_DOUBLE, freeWorker, 11, MPI_COMM_WORLD);
						MPI_Send(best_cut, n, MPI_INT, freeWorker, 12, MPI_COMM_WORLD);
						MPI_Send(&bound_lb, 1, MPI_DOUBLE, freeWorker, 13, MPI_COMM_WORLD);
						MPI_Send(&bound_ub, 1, MPI_DOUBLE, freeWorker, 14, MPI_COMM_WORLD);

						//we free allocated element of queue
						free_Tree_node(cur_node1);
						free(elt);
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
				
				break;
			}
			if(tag == TIME_LIMIT)
			{
				over = 1;
				for(i = 1; i <= numbWorkers-1; i++)
				{
					MPI_Send(&over, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
				}

				while(itemCount != 0)
				{
					struct queue *elt;
					Tree_node *cur_node1;
					elt = dequeue();
					cur_node1 = elt->node;

					free_Tree_node(cur_node1);
					free(elt);
				}
				
				// print approximate solution
				print_solution(n, mm, best_cut, best_value, n_bbnodes, MPI_Wtime() - start, density,0);
				outputFile(argv, n, mm, best_cut, best_value, n_bbnodes, MPI_Wtime() - start, density, 0);			
			
				break;
			}
		}

	}
	else
	{
		int over;
		int active;

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
					MPI_Recv(&(cur_node->flip), 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
					MPI_Recv(&(cur_node->correction), 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD, &status);
					MPI_Recv(&(cur_node->ub), 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD, &status);
					MPI_Recv(cur_node->star_value, n, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD, &status);
					MPI_Recv(cur_node->cost.alive, n, MPI_INT, 0, 6, MPI_COMM_WORLD, &status);
					MPI_Recv(cur_node->cost.flip, n, MPI_INT, 0, 7, MPI_COMM_WORLD, &status);
					MPI_Recv(cur_node->cost.ref, n, MPI_INT, 0, 8, MPI_COMM_WORLD, &status);
					MPI_Recv(cur_node->cost.link, n, MPI_INT, 0, 9, MPI_COMM_WORLD, &status);
					MPI_Recv(&(cur_node->cost.m[1][0]), (n*(n-1)/2), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD, &status);
					MPI_Recv(&best_value, 1, MPI_DOUBLE, 0, 11, MPI_COMM_WORLD, &status);
					MPI_Recv(best_cut, n, MPI_INT, 0, 12, MPI_COMM_WORLD, &status);
					MPI_Recv(&bound_lb, 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD, &status);
					MPI_Recv(&bound_ub, 1, MPI_DOUBLE, 0, 14, MPI_COMM_WORLD, &status);
					upper_bound = cur_node->ub;
					if(upper_bound < best_value + P.gap_relax)
					{
						tag = FATHOM;
						MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
						MPI_Send(&tag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
						MPI_Send(&bound_lb, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
					}
					else
					{
						int pos = 0;
						for(i = 0; i < n; i++)
						{
							if(cur_node->cost.alive[i])
							{
								new_pos[i] = pos++;
							}
						}
						n_supernodes = 1;
						original_name[0] = 0;
						for(i = 1, l = 1; i < n; i++)
						{
							if(cur_node->cost.alive[i])
							{
								n_supernodes++;
								original_name[l] = i;
								for(j = 0, k = 0; j < i; j++)
								{
									if(cur_node->cost.alive[j])
										tmp_cost.m[l][k++] = cur_node->cost.m[i][j];
								}
								l++;
							}
						}
						for(i = 0; i < n_supernodes; i++)
						{
							bound_cut[i]= (double)(best_cut[original_name[i]]*cur_node->cost.flip[original_name[i]]);
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
							bound_cut[original_name[i]]= bound_cut[i]*cur_node->cost.flip[original_name[i]];
						}
						for(i= 0;i<n;i++)
							if(!cur_node->cost.alive[i]){
								bound_cut[i]= bound_cut[cur_node->cost.ref[i]]*cur_node->cost.link[i];
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
						int fathomed = (p_res== 2)||(bound_ub < best_value+P.gap_relax);
						if(!fathomed){
							if(MPI_Wtime()-start > SECONDS)//seconds
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
								MPI_Send(&(cur_node->flip), 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
								MPI_Send(&(cur_node->correction), 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
								MPI_Send(&(cur_node->ub), 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
								MPI_Send(cur_node->star_value, n, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
								MPI_Send(cur_node->cost.alive, n, MPI_INT, 0, 6, MPI_COMM_WORLD);
								MPI_Send(cur_node->cost.flip, n, MPI_INT, 0, 7, MPI_COMM_WORLD);
								MPI_Send(cur_node->cost.ref, n, MPI_INT, 0, 8, MPI_COMM_WORLD);
								MPI_Send(cur_node->cost.link, n, MPI_INT, 0, 9, MPI_COMM_WORLD);
								MPI_Send(&(cur_node->cost.m[1][0]), (n*(n-1)/2), MPI_DOUBLE, 0, 10, MPI_COMM_WORLD);
								MPI_Send(&bound_left, 1, MPI_INT, 0, 11, MPI_COMM_WORLD);
								MPI_Send(&bound_right, 1, MPI_INT, 0, 12, MPI_COMM_WORLD);
								MPI_Send(&bound_lb, 1, MPI_DOUBLE, 0, 13, MPI_COMM_WORLD);
							}			
						}
						else
						{
							tag = FATHOM;
							MPI_Send(&rank, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
							MPI_Send(&tag, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
							MPI_Send(&bound_lb, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
						}
					}
				}
			}
		}
	}
	
	free(original_name);
	free(bound_cut);
	free(prim_cut);
	free(new_pos);
	free_Tree_node(cur_node);
	free(cur_node);
	free(best_cut);
	if(rank == 0)
		free(busyWorkers);
	//if(num_con){
		//free(tmp_sense);
		//free(tmp_rhs);
		//free(tmp_bg);
		//free(tmp_sp);
		//free(tmp_co);
	//}
	free_mat(tmp_cost);
	free_mat(orig_cost);

	MPI_Finalize();

	return 0;
}

