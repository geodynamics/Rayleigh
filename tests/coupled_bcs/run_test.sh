#!/usr/bin/env bash

cd tests/coupled_bcs

# we run a const case, where the bcs on the field that everything is coupled to are just constant
cd const
PYTHONPATH=../../../pre_processing:$PYTHONPATH python3 generate_input.py
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

# after running, we test the output for errors
PYTHONPATH=../../post_processing:../../pre_processing:$PYTHONPATH python3 test_output.py

