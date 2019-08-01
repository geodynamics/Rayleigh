# CU/NCAR Janus
DOPTS = 1

.SUFFIXES: .o .f90

MPIF90 = ifort

BASEOPTS = -openmp -real-size 32 -O3 -ipo -align all -align rec16byte -ftz -vec-threshold5 -auto -opt-subscript-in-range -ansi-alias -fno-fnalias -fno-alias -opt-prefetch -no-prec-div -unroll

COPTS  = -openmp -O3 -m64 -ipo -xSSE3

OPTSF90AMD = -m64 -xSSE3 $(BASEOPTS)

OPTSF90WES = -xHost $(BASEOPTS)

OPTSF90_ERR = -fpp -openmp -O0 -traceback -g -check all -debug all -warn all


ifeq (${DOPTS},0)
	OPTSF90=$(OPTSF90_ERR)
	LIBS = -lstdc++
	EXEC="interpolator.dbg"
	OBJSMPI = interp.o
else ifeq (${DOPTS},1)
	OPTSF90=$(OPTSF90WES)
	LIBS = -lstdc++
	EXEC="interpolator.wes"
	OBJSMPI = interp.o
else ifeq (${DOPTS},2)
	OPTSF90=$(OPTSF90AMD)
	EXEC="interpolator.amd"
	LIBS = -lstdc++
	OBJSMPI = interp.o amd_optimization.o
else
	OPTSF90=$(OPTSF90_ERR)
	LIBS = -lstdc++
	EXEC="interpolator.err"
	OBJSMPI = interp.o
endif

main:	driver.f90 $(OBJSMPI)
	$(MPIF90) $(OPTSF90) -o $(EXEC) driver.f90 $(OBJSMPI) $(LIBS)

clean:
	touch temp.o temp.mod
	rm *.o *.mod

.f90.o:
	$(MPIF90) $(OPTSF90) -c $< 

.cpp.o:
	icpc $(COPTS) -c $<

