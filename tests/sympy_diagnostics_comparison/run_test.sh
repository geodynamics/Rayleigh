#!/usr/bin/env bash

cd tests/sympy_diagnostics_comparison

# no-slip case: single-iteration no-slip velocity/pressure/temperature IC plus
# a magnetic (c_init/a_init) IC, compared against an independent sympy
# reconstruction of the raw diagnostic output (see no_slip/compare.py)
cd no_slip
echo "================="
echo "running no_slip"
echo "================="
python3 generate_input.py
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

# free-slip case: same idea, no magnetism (see free_slip/compare.py) but the
# velocity field gets a toroidal component to test the curl of the viscous
# force
cd free_slip
echo "================="
echo "running free_slip"
echo "================="
python3 generate_input.py
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

# after both scenarios have run, we test the output for errors
(cd no_slip && python3 compare.py) && (cd free_slip && python3 compare.py)
