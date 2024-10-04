#########################################
## Make file for Stable Set problem ##
#########################################
CC				= mpic++
PROG            = StableSetBaB

OBJS            = StableSetBaB.o \
				BoundComputationTheta.o \
				Subproblem.o \
				bpm.o \
				ipm_Lovasz_theta.o \
				ReadAndOutput.o \
				Graph.o \
				Heuristic.o

# headers
ARMADILLO_INC	= ./armadillo-9.850.1/include
INCLUDES        = -I$(ARMADILLO_INC)

LIBS            = -DARMA_USE_OPENMP -DARMA_DONT_USE_WRAPPER -DARMA_NO_DEBUG -DARMA_OPTIMISE_SOLVE_SYMPD -lopenblas -lm
CFLAGS        	= -O3 -std=c++11 -fopenmp -DNDEBUG 

.PRECIOUS: %.o

all: $(PROG)

%.o: %.cpp
	$(CC) -c -o $@ $^ $(INCLUDES) $(LIBS) $(CFLAGS)

$(PROG): $(OBJS)
	$(CC) $(CFLAGS) -o $(PROG) $(OBJS) $(INCLUDES) $(LIBS) 

clean:
	rm -f *.o core
