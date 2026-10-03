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

cp -r T/Checkpoints T.check/.

cd T.check
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

cd chi
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

cp -r chi/Checkpoints chi.check/.

cd chi.check
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

# chi.magnetic repeats the chi setup with magnetism = .true. (magnetic field
# held at zero) to test magnetism + active/passive scalar. 
# It runs 19 iterations, forcing a checkpoint at
# iteration 17 (in addition to the usual final-iteration one at 19); its own
# iteration-18/19 output is the reference for the restart tests
# below. chi.magnetic.check restarts from the iteration-17 checkpoint and
# runs one more step to iteration 18, which test_output.py compares directly
# against chi.magnetic's iteration-18 output.
cd chi.magnetic
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

cp -r chi.magnetic/Checkpoints chi.magnetic.check/.

cd chi.magnetic.check
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

# chi.restart_magnetism_on restarts from chi's checkpoint (magnetism=.false.)
# with magnetism turned on (B held at zero), and runs 2 steps (to iteration
# 19, not 1: iteration 18 just re-reports the checkpoint's own state).
# test_output.py compares its iteration-19 result against
# chi.magnetic's iteration-19 result: since B == 0 on both sides, turning
# magnetism on only at the restart vs. having it on for the whole run should
# be indistinguishable.
cp -r chi/Checkpoints chi.restart_magnetism_on/.

cd chi.restart_magnetism_on
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

# chi.magnetic.restart_magnetism_off is the reverse: restarts from
# chi.magnetic's checkpoint (magnetism=.true., B == 0, there) with magnetism
# turned off. Its C/A files are never read (magnetism=.false.). 
# Also runs 2 steps to iteration 19, compared the same way against
# chi.magnetic's iteration-19 result.
cp -r chi.magnetic/Checkpoints chi.magnetic.restart_magnetism_off/.

cd chi.magnetic.restart_magnetism_off
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

# T.restart_no_scalars restarts from T's checkpoint (n_active_scalars=2,
# n_passive_scalars=2 there) into a run with no scalar fields. 
# Compared against T.check's iteration-18 result.
cp -r T/Checkpoints T.restart_no_scalars/.

cd T.restart_no_scalars
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

# after all versions have run, we test the output for errors
PYTHONPATH=../../post_processing:../../pre_processing:$PYTHONPATH python3 test_output.py

