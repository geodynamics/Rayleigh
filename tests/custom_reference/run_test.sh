#!/usr/bin/env bash

cd tests/custom_reference

# Run a test where the full reference state is specified
cd basic
cp ../../../post_processing/reference_tools.py .
python3 gen_case0.py

mpirun -np 4 ../../../bin/rayleigh.dbg

