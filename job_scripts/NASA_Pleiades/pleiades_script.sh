#PBS -S /bin/bash
#PBS -q long                # 'long' queue allows jobs up to 5 days (use 'normal' for 8 hours or less).
#PBS -j oe
#PBS -W group_list=s1234    # account number
#PBS -m e
#PBS -N nrho5_ff_5          # job name
#PBS -l walltime=120:00:00  # run for 120 hours
#PBS -l select=74:ncpus=28:mpiprocs=28:model=bro  # Request 74 Broadwell nodes.  28 MPI ranks per node.

module purge
module load comp-intel
module load mpi-hpe
export OMP_NUM_THREADS=1

# Link the Rayleigh executables from your home directory (modify the path below accordingly)
ln -s /home4/nfeather/Ra_share/Rayleigh/bin/* .
sleep 10

# Run the code
# Make sure that the number of cores is <= select*mpiprocs
mpiexec -np 2048 ./rayleigh.avx2 -nprow 64 -npcol 32


# Compilation Notes
#
# To configure and build Rayleigh, the following commands should suffice:
#
# module purge
# module load comp-intel
# module load mpi-hpe
# ./configure --FC=mpif90 --CC=icc  (select 'ALL')
# make -j
# make install

# Additional notes
# (1) Pleiades is a hetergenous cluster, composed of many (primarily Intel) processor types.
#     We suggest selecting the 'ALL' option when configuring Rayleigh to ensure that a unique
#     executable is created for each of the possible vectorization options.3
#
# (2) This script is set up to run Rayleigh on 2048 Pleiades Broadwell cores.
#     Each Broadwell node has 28 cores, so 74 nodes (2072 cores) must be requested
#     in order to have access to 2048 cores.  This is typical of most Pleiades node
#     configurations as most node types do not possess a core-count that is a power of 2.
#
# (3) To change the node type, change (i) the PBS 'select' line and the Rayleigh 
#     executable correspondingly.  Some alternative examples are:
#
#   a) Sandybridge Nodes (128 nodes, 2048 cores):  
#                   PBS -l select=128:ncpus=16:mpiprocs=16:model=san
#                   mpiexec -np 2048 ./rayleigh.avx -nprow 64 -npcol 32
#   b) Ivybridge Nodes (103 nodes, 2060 cores):
#                   PBS -l select=103:ncpus=20:mpiprocs=20:model=ivy
#                   mpiexec -np 2048 ./rayleigh.avx -nprow 64 -npcol 32
#   c) Haswell Nodes (86 nodes, 2064 cores):
#                   PBS -l select=86:ncpus=24:mpiprocs=24:model=has
#                   mpiexec -np 2048 ./rayleigh.avx2 -nprow 64 -npcol 32

