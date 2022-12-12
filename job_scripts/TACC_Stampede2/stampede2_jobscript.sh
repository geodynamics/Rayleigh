#!/bin/bash
#SBATCH -J geodynamo-test  # Job name
#SBATCH -o log.o%j         # Name of stdout output file
#SBATCH -e error.e%j       # Name of stderr error file
#SBATCH -p skx-dev         # Queue (partition) name; skx-dev for testing; skx-normal for production.
#SBATCH -N 1               # Total # of nodes. 1-4 for skx-dev. >=4 for skx-normal
#SBATCH --ntasks-per-node 48
#SBATCH -t 00:10:00        # Run time (hh:mm:ss) max 2h on skx-dev, max 48h on skx-normal
#SBATCH --mail-user=
#SBATCH --mail-type=none   # Send no email

module list

# Launch MPI code...

export OMP_NUM_THREADS=1
# Replace -n X with correct number of MPI ranks, or remove to use all ranks
# requested on the nodes above.

ibrun -n 48 ./rayleigh.opt

# Compilation Notes
# At the time this documentation was written, the default-loaded module 
# stack works without problems (Sep 2022).
#
# To configure and compile the code, the following commands should suffice:
#
# FC=mpifc CC=mpicc ./configure  (choose 'AVX512')
# make -j
# make install
