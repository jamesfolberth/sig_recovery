
.PHONY: all clean force

FORT=gfortran

FORTOPTS= -Wall -Wextra -fbounds-check -funroll-loops -finline -ffast-math -march=native -mtune=native -O2

LDLIBS= -lsatlas -L/opt/hdf5/lib -lhdf5_fortran
#LDLIBS= -lblas # openBLAS
INC= -I/opt/hdf5/include
LDADDPATH= -Wl,-rpath=/opt/hdf5/lib

all: compla.o tropp_fig2_comp

compla.o: compla.f08
	$(FORT) -c compla.f08 $(FORTOPTS) $(LDLIBS) $(INC) $(LDADDPATH)

omp.o: omp.f08
	$(FORT) -c omp.f08 $(FORTOPTS) $(LDLIBS) $(INC) $(LDADDPATH)

test: test.f08 compla.o force
	$(FORT) -c test.f08 $(FORTOPTS) $(INC)
	$(FORT) -o test $(FORTOPTS) $(LDLIBS) $(LDADDPATH) test.o compla.o
	./test

tropp_fig1_comp: tropp_fig1_comp.f08 compla.o omp.o force
	$(FORT) -c tropp_fig1_comp.f08 $(FORTOPS) $(INC)
	$(FORT) -o tropp_fig1_comp $(FORTOPTS) $(LDLIBS) $(LDADDPATH) tropp_fig1_comp.o compla.o omp.o
	./tropp_fig1_comp

tropp_fig2_comp: tropp_fig2_comp.f08 compla.o omp.o force
	$(FORT) -c tropp_fig2_comp.f08 $(FORTOPS) $(INC)
	$(FORT) -o tropp_fig2_comp $(FORTOPTS) $(LDLIBS) $(LDADDPATH) tropp_fig2_comp.o compla.o omp.o
	./tropp_fig2_comp

# phony target to force make to build targets that depend on force
force:

clean:
	rm -f test tropp_fig1_comp
	rm -f *.o *.mod *.out *.a *.h5
