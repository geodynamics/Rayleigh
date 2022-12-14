#!/usr/bin/env bash

cd tests/custom_reference

# Run a test where the full reference state is specified
cd basic
cp ../../../post_processing/reference_tools.py .
python3 gen_case0.py
mpirun -np 4 ../../../bin/rayleigh.dbg


# Augment the usual case 0 with heating and higher Ra
cd ..
cd augment
cp ../../../post_processing/reference_tools.py .
python3 augment.py
mpirun -np 4 ../../../bin/rayleigh.dbg
