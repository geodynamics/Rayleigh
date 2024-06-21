#!/bin/bash
#SBATCH -J geodynamo-test  # Job name
#SBATCH -o log.o%j         # Name of stdout output file
#SBATCH -e error.e%j       # Name of stderr error file
#SBATCH -p skx-dev         # Queue (partition) name; skx-dev for testing; skx, icx, or spr for production.
#SBATCH -N 1               # Total # of nodes. 1-4 for skx-dev. >=4 for skx-normal
#SBATCH --ntasks-per-node 48 # 48 for skx-dev and skx, 80 for icx, 112 for spr
#SBATCH -t 00:10:00        # Run time (hh:mm:ss) max 2h on skx-dev, max 48h on everything else
#SBATCH --mail-user=
#SBATCH --mail-type=none   # Send no email

module list

# Launch MPI code...

export OMP_NUM_THREADS=1

# Replace this with the path and name of the Rayleigh binary
export RAYLEIGH_BINARY=./rayleigh

# Replace -n X with correct number of MPI ranks, or remove to use all ranks
# requested on the nodes above. This assumes you submit
# the job from the directory of your main_input file.
ibrun -n 48 ${RAYLEIGH_BINARY}
