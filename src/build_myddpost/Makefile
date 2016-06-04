# Makefile for myDDPOSTPROCESS:
# A modified version of the Makefile for DDSCAT

# define the following:
PRECISION	= dp
OPENMP		=
FC		= gfortran
FFLAGS	= -O2
LFLAGS	=
PROG    = myddpostprocess
BIN     = ../../bin
#

#%.o: %.f90 ddprecision.mod ddcommon_1.mod cgmodule.mod
#	$(FC) -c $(FFLAGS) $(OPENMP) $< -o $@

myDDPOSTPROCESS.o: myDDPOSTPROCESS.f90 ddprecision.mod readnf_bcom.mod \
	readnf_ecom.mod
	$(FC) -c $(FFLAGS) $(OPENMP) myDDPOSTPROCESS.f90 \
	-o myDDPOSTPROCESS.o
	
readnf.o: readnf.f90 ddprecision.mod readnf_bcom.mod readnf_ecom.mod
	$(FC) -c $(FFLAGS) $(OPENMP) readnf.f90 \
	-o readnf.o

# dependencies for myDDPOSTPROCESS
myOBJS3 =	myDDPOSTPROCESS.o\
	readnf_bcom.o\
	readnf_ecom.o\
	readnf.o

all:	$(PROG)

myddpostprocess: ddprecision.mod readnf_bcom.mod readnf_ecom.mod $(myOBJS3)
	$(FC) -o myddpostprocess \
	$(myOBJS3) $(LFLAGS) 

#--------------- modules ---------------------------------------------

ddprecision.mod: ddprecision.f90
	cpp -P -traditional-cpp -D$(PRECISION) ddprecision.f90 \
	ddprecision_cpp.f90
	$(FC) -c $(FFLAGS) ddprecision_cpp.f90 -o ddprecision.o
	rm ddprecision_cpp.f90

readnf_bcom.mod: readnf_bcom.f90
	$(FC) -c readnf_bcom.f90 -o readnf_bcom.o

readnf_ecom.mod: readnf_ecom.f90
	$(FC) -c readnf_ecom.f90 -o readnf_ecom.o

#---------------------------------------------------------------------

install: all
	cp -f $(PROG) $(BIN)

clean:	
	rm -f *.o make.out*  *.mod $(PROG)

distclean: clean
	rm -f $(BIN)/$(PROG)

veryclean: clean
	rm -f $(PROG)
