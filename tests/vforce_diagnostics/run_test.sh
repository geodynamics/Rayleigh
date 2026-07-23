#!/usr/bin/env bash

cd tests/vforce_diagnostics

mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../bin/rayleigh.dbg

# after the run, we test the output for errors
PYTHONPATH=../../post_processing:../../pre_processing:$PYTHONPATH python3 test_output.py
