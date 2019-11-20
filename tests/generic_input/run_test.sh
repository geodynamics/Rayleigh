#!/usr/bin/env bash

cd tests/generic_input

# first we run the "base" case, this just sets up and runs Rayleigh for one time-step with hard-coded initial conditions 
# for the Christensen et al., 2001 benchmark case 1
cd base
mpirun -np 4 ../../../bin/rayleigh.dbg
../../../post_processing/convert_full3d_to_vtu.py
cd ..

# next we run the same thing but using generic input files to set the initial conditions
cd script
# we generate the generic input files two different ways so that both are tested and demonstrated...
# first we generate the temperature initial condition using rayleigh_spectral_input.py in "script" mode:
../../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_t_init \
   -e 'import numpy as np; x = 2*radius - rmin - rmax; rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) - x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'
# then we use a custom python script using rayleigh_spectral_input.py as a module to write the magnetic initial conditions
PYTHONPATH=../../../pre_processing:$PYTHONPATH python generate_magnetic_input.py
# finally we run Rayleigh
mpirun -np 4 ../../../bin/rayleigh.dbg
../../../post_processing/convert_full3d_to_vtu.py
cd ..

# onto testing bcs... again, first the base case using hard-coded Rayleigh options
cd bcs_base
../../../pre_processing/rayleigh_spectral_input.py -m 0 0 0 0.0+0.j -o zero_init_vol
mpirun -np 4 ../../../bin/rayleigh.dbg
cd ..

# then again onto the same thing but using generic input files to set the boundary conditions
cd bcs_script
../../../pre_processing/rayleigh_spectral_input.py -m 0 0 0 0.0+0.j -o zero_init_vol
../../../pre_processing/rayleigh_spectral_input.py -m 1 0 5.5225155517445783 -m 1 1 1.3771166096309384+0.j -o ctop_init_bc
../../../pre_processing/rayleigh_spectral_input.py -m 1 0 15.778615862127371 -m 1 1 3.9346188846598249+0.j -o cbottom_init_bc
../../../pre_processing/rayleigh_spectral_input.py -e '5.0' -o five_init_bc
../../../pre_processing/rayleigh_spectral_input.py -e '-5.0' -o mfive_init_bc
mpirun -np 4 ../../../bin/rayleigh.dbg
cd ..

# after both versions have run, we test the output for errors
PYTHONPATH=../../post_processing:../../pre_processing:$PYTHONPATH python test_output.py

