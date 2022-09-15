Steps to use this container on Stampede2:

module load tacc-singularity
singularity pull docker://gassmoeller/rayleigh-tacc
ibrun -n X singularity run rayleigh-tacc_latest.sif /opt/Rayleigh/bin/rayleigh.opt

Timing of the native Rayleigh installation (Intel compiler, MKL) on the skx-dev queue and $WORK file system with 4 cores:

//////////////////////////////////////////////
   Measured Timings for Process 0  (seconds)

 Elapsed time:        13.0618
  Column time:         0.7478
     Row time:         6.1813
Legendre time:         1.2880
     FFT time:         0.8845
   Solve time:         0.4607
    rlma time:         0.2792
    rlmb time:         0.1068
  pspace time:         1.2199
  psolve time:         1.5894
    dphi time:         0.1590
captured time:        12.9166

     iter/sec:        30.5471
//////////////////////////////////////////////
TACC:  Shutdown complete. Exiting. 


Timing of the singularity container (standard LAPACK/BLAS, Intel compiler) on the skx-dev queue and $SCRATCH file system with 4 cores:

//////////////////////////////////////////////
   Measured Timings for Process 0  (seconds)

 Elapsed time:        21.9945
  Column time:         3.8781
     Row time:         7.4885
Legendre time:         3.8041
     FFT time:         2.2408
   Solve time:         0.6536
    rlma time:         0.2052
    rlmb time:         0.0768
  pspace time:         1.1225
  psolve time:         2.1921
    dphi time:         0.1550
captured time:        21.8168

     iter/sec:        18.1409
//////////////////////////////////////////////

Difference is likely that I have no access to MKL in the docker container (since it has to be based on Ubuntu 18.04).