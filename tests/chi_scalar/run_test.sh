#!/usr/bin/env bash

cd tests/chi_scalar

#cd bench
#mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
#cd ..

# generate initial conditions
../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_t_init \
   -e 'import numpy as np; x = 2*radius - rmin - rmax; rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) - x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

cd T
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

cd chi
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

# after both versions have run, we test the output for errors
PYTHONPATH=../../post_processing:../../pre_processing:$PYTHONPATH python test_output.py

