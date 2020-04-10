PROG = interp3d

BASEOPTS_GNU = -cpp -fopenmp -ffree-form -ffree-line-length-none -O3
OPTSF90_GNU = $(BASEOPTS_GNU) -march=native -mtune=native
F90_GNU = gfortran

BASEOPTS_INTEL1 = -qopenmp -real-size 32 -O3 -ipo -align all -align rec16byte -ftz -vec-threshold5 -auto 
BASEOPTS_INTEL2 = -qopt-subscript-in-range -ansi-alias -fno-fnalias -fno-alias -qopt-prefetch -no-prec-div -unroll
OPTSF90_INTEL = -xHost $(BASEOPTS_INTEL1) $(BASEOPTS_INTEL2)
F90_INTEL = ifort

interp3d.gnu: 
	$(F90_GNU) $(OPTSF90_GNU) -c Interpolation.F90
	$(F90_GNU) $(OPTSF90_GNU) -c Read_CMD.F90
	$(F90_GNU) $(OPTSF90_GNU)  Main.F90 Interpolation.o Read_CMD.o -o $(PROG)

interp3d.intel: 
	$(F90_INTEL) $(OPTSF90_INTEL) -c Interpolation.F90
	$(F90_INTEL) $(OPTSF90_INTEL) -c Read_CMD.F90
	$(F90_INTEL) $(OPTSF90_INTEL)  Main.F90 Interpolation.o Read_CMD.o -o $(PROG)

clean:
	touch temp.o temp.mod
	rm *.o *.mod $(PROG)


