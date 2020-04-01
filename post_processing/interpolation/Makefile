# For Mac OS 10.9 "Mavericks" or linux gfortran
.SUFFIXES: .o .F90

MPIF90 = gfortran
PROG = interp3d

BASEOPTS = -cpp -fopenmp -ffree-form -ffree-line-length-none -O3

OPTSF90_MAC = $(BASEOPTS) -march=native -mtune=native -Wa,-q

OPTSF90_FST = $(BASEOPTS) -march=native -mtune=native

OPTSF90_ERR = -cpp -fopenmp -ffree-form -ffree-line-length-none -O0 -fbacktrace -Wall -fcheck=all

ifeq (${MAKECMDGOALS},clean)
     DOPTS = 3
else
     DOPTS = 1
endif

OBJSMPI = Main.F90 Interpolation.o Read_CMD.o

ifeq (${DOPTS},0)
     OPTSF90 = $(OPTSF90_ERR)
else ifeq (${DOPTS},1)
     OPTSF90 = $(OPTSF90_FST)
else ifeq (${DOPTS},2)
     OPTSF90 = $(OPTSF90_MAC)
endif

main: $(if $(strip $(DOPTS)), $(OBJSMPI), $(error Specify the DOPTS flag to select the desired optimization level)) 
	$(MPIF90) $(OPTSF90)  $(OBJSMPI) -o $(PROG)

clean:
	touch temp.o temp.mod
	rm *.o *.mod $(PROG)

.F90.o:
	$(MPIF90) $(OPTSF90) -c $< 

