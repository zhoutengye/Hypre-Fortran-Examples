FC = gnu
F90 = mpif90


LINKOPTS  = $(COPTS)
LIBS      = -L$(HYPRE_DIR)/lib -lHYPRE -lm
FOPTS     = -g
FINCLUDES = -I$(HYPRE_DIR)/include
FFLAGS    = $(FOPTS) $(FINCLUDES)
LFLAGS    = $(LINKOPTS) $(LIBS) -lstdc++

ifeq ($(FC),$(filter $(FC), intel))
	HYPRE_DIR = /home/yeldon/hypre/intel
	COPTS     = -g -fPIC
else ifeq ($(FC),$(filter $(FC), gnu))
	HYPRE_DIR = /home/yeldon/hypre/gnu
	COPTS     = -g -Wall -fPIC
endif

.SUFFIXES: .f90 .o

.f90.o:
	$(F90) $(FFLAGS) -c $<

.o.exec:
	$(F90) -o $@ $^ $(LFLAGS)

ALLPROGS = ex1 ex2 ex3

all: $(ALLPROGS)

ex1: ex1.o
	$(F90) -o $@ $^ $(LFLAGS)

ex2: ex2.o
	$(F90) -o $@ $^ $(LFLAGS)

ex3: ex3.o
	$(F90) -o $@ $^ $(LFLAGS)

clean:
	rm -f *.o
	rm -f $(ALLPROGS)
