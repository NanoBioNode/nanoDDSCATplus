# Makefile for pip. Uses the default gfortran, which should be present
# on any Unix. Other compilers and optimization flags may also be used
# - should be adjusted below.

CF      = ifort
FFLAGS  = -fast
EXTRA_FLAGS = -openmp
FSOURCE = FEM-Geo-Wr.f90 geometry.f90 ivread_wr.f90 limits_cor3.f90 \
	  point-inside-polyhedron.f90
PROG    = pip
BIN     = ../../bin
SETUP   = . /etc/environ.sh; use -e -r intel-fc-14.0.0

FOBJECTS := $(FSOURCE:.f90=.o)


.PHONY: all clean install distclean

all: $(PROG)

$(PROG): $(FOBJECTS)
	$(SETUP); $(CF) -o $@ $(FOBJECTS) $(EXTRA_FLAGS)

$(FOBJECTS): %.o: %.f90
	$(SETUP); $(CF) -c $(FFLAGS) $(EXTRA_FLAGS) $<

install: all
	cp -f $(PROG) $(BIN)

clean:
	rm -f *.o $(PROG)

distclean: clean
	rm -f $(BIN)/$(PROG)

