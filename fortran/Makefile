
.PHONY: all clean force

FORT=gfortran

FORTOPTS= -Wall -Wextra -fbounds-check -funroll-loops -finline -ffast-math -march=native -mtune=native -O2

LDLIBS= -lsatlas -L/opt/hdf5/lib -lhdf5_fortran
#LDLIBS= -lblas # openBLAS
INC= -I/opt/hdf5/include
LDADDPATH= -Wl,-rpath=/opt/hdf5/lib

all: compla.o test

compla.o: compla.f08
	$(FORT) -c compla.f08 $(FORTOPTS) $(LDLIBS) $(INC) $(LDADDPATH)

test: test.f08 compla.o force
	$(FORT) -c test.f08 $(FORTOPTS) $(INC)
	$(FORT) -o test $(FORTOPTS) $(LDLIBS) $(LDADDPATH) test.o compla.o
	./test

# phony target to force make to build targets that depend on force
force:

clean:
	rm -f test save_stuff.h5
	rm -f *.o *.mod *.out *.a