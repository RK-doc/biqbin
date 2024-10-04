

CFLAGS = -g -O3 -Wall -I.

#LDFLAGS= -static

CC=mpicc

o = .o

OBJECTS=bab_for_klu$o bound$o bdl_main$o bdl_mc2$o \
	separation$o tri_sep2$o heuristic$o mc_cut$o mc_1opt$o \
	fct_eval$o mc_a_at$o mc_psdpk$o \
	operators$o lam_eta$o solvelambda$o gb_flip$o op_b_bt$o

#.c.o:
#	$(CC) $(CFLAGS) -I. -c $*.c -o $*.o

all: main

main: lib demo/bab_main$o
	$(CC) $(LDFLAGS) -o bab demo/bab_main.o -L. -L/usr/lib/atlas -L/usr/lib/sse2 -L/usr/lib/atlas/sse2 \
	-lbabsdp -llapack -lf77blas -lcblas -latlas -lm -lgfortran

con:  lib demo/bab_con$o
	$(CC) $(LDFLAGS) -o con demo/bab_con.o -L. -L/usr/lib/sse2 -L/usr/lib/atlas/sse2 \
	-lbabsdp -llapack -lf77blas -lcblas -latlas -lg2c -lm -lgfortran

qkn:  lib demo/bab_qkn$o
	$(CC) $(LDFLAGS) -o qkn demo/bab_qkn.o -L. -L/usr/lib/sse2 -L/usr/lib/atlas/sse2 \
	-lbabsdp -llapack -lf77blas -lcblas -latlas -lg2c -lm -lgfortran

lib: $(OBJECTS)
	ar cr libbabsdp.a $(OBJECTS)

clean:
	rm -f *.o libbabsdp.a bab
