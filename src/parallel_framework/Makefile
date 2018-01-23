FFTWINC = /home3/feathern/FFTW3/include
FFTWLIB = /home3/feathern/FFTW3/lib
F90  = /usr/local/mpich2-1.0.7/bin/mpif90
F90FLAGS = -ip -ipo -r8 -O3 -fpp -FR -fp-model strict -I$(FFTWINC)
LIBFLAGS = -L/home3/feathern/FFTW3/lib/ -lfftw3 -lm -lc 
F90FLAGS = -r8 -O0 -CB -g -traceback -fpp -FR 
##LIBFLAGS = -lmpi -lmpio
PROG =	pfw

.SUFFIXES: .o .F90 .f .F
OBJ = MPI_Base.o All_to_All.o MPI_LAYER.o Load_Balance.o Parallel_Framework.o \


$(PROG) : $(OBJ)

.F90.o :
	$(F90) $(F90FLAGS) -c $<
clean : 
	rm -f *.o *.mod

