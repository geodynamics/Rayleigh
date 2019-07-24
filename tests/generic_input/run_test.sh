#!/usr/bin/env bash

cd tests/generic_input
cd base
mpirun -np 4 --oversubscribe ../../../bin/rayleigh.dbg
cd ..
cd script
../../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_t_init \
   -e 'import numpy as np; x = 2*radius - rmin - rmax; rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) - x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'
mpirun -np 4 --oversubscribe ../../../bin/rayleigh.dbg
cd ..
PYTHONPATH=../../post_processing:$PYTHONPATH python test_output.py

