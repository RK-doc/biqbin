#########################################
## Make file for Stable Set problem ##
## module load openmpi
## module load intel-mkl
## export OMP_NUM_THREADS=2
#########################################
CC				= mpicc

OBJECTS			= bab_for_klu.o \
				bound.o \
				bdl_main.o \
				bdl_mc2.o \
				separation.o \
				tri_sep2.o \
				heuristic.o \
				mc_cut.o \
				mc_1opt.o \
				fct_eval.o \
				mc_a_at.o \
				mc_psdpk.o \
				operators.o \
				lam_eta.o \
				solvelambda.o \
				gb_flip.o \
				op_b_bt.o

LIBS			= -mkl -L. -lbabsdp
CFLAGS			= -O3 -fopenmp -DNDEBUG -I.

.c.o:
	$(CC) $(CFLAGS) -c $*.c -o $*.o

all: main

main: lib demo/bab_main.o
	$(CC) -o bab demo/bab_main.o $(LIBS)

lib: $(OBJECTS)
	ar cr libbabsdp.a $(OBJECTS)

clean:
	rm -f *.o libbabsdp.a bab demo/bab_main.o
