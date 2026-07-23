#!/usr/bin/env bash

set -x

cd tests/chi_custom_reference

# generate initial conditions
# ../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o const_chi_init \
#    -e '1.0'

# ../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o diff_t_init \
#    -e 'rmax*rmin/radius - rmin'

../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_t_init \
   -e 'import numpy as np; x = 2*radius - rmin - rmax; rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) - x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

# ../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_dtdr_init \
#    -e 'import numpy as np; x = 2*radius - rmin - rmax; xp = 2; -rmax*rmin*(radius**-2) + 210*0.1*(- 6*x*xp + 12*(x**3)*xp - 6*(x**5)*xp)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_dtdr_rmin_init \
   -e 'import numpy as np; x = rmin - rmax; xp = 2; -rmax/rmin + 210*0.1*(- 6*x*xp + 12*(x**3)*xp - 6*(x**5)*xp)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_dtdr_rmax_init \
   -e 'import numpy as np; x = rmax - rmin; xp = 2; -rmin/rmax + 210*0.1*(- 6*x*xp + 12*(x**3)*xp - 6*(x**5)*xp)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_tp_init \
   -e 'import numpy as np; x = 2*radius - rmin - rmax; 210*0.1*(1 - 3*x*x + 3*(x**4) - x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

# ../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_dtpdr_init \
#    -e 'import numpy as np; x = 2*radius - rmin - rmax; xp = 2; 210*0.1*(- 6*x*xp + 12*(x**3)*xp - 6*(x**5)*xp)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_dtpdr_rmin_init \
   -e 'import numpy as np; x = rmin - rmax; xp = 2; 210*0.1*(- 6*x*xp + 12*(x**3)*xp - 6*(x**5)*xp)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

../../pre_processing/rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 64 -nr 48 -o bench_dtpdr_rmax_init \
   -e 'import numpy as np; x = rmax - rmin; xp = 2; 210*0.1*(- 6*x*xp + 12*(x**3)*xp - 6*(x**5)*xp)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'
   

cd T.full
PYTHONPATH=../../../post_processing:../../../pre_processing:$PYTHONPATH python3 gen_ref_state.py
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

cd T.split
PYTHONPATH=../../../post_processing:../../../pre_processing:$PYTHONPATH python3 gen_ref_state.py
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

cd chi.full
PYTHONPATH=../../../post_processing:../../../pre_processing:$PYTHONPATH python3 gen_ref_state.py
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..

cd chi.split
PYTHONPATH=../../../post_processing:../../../pre_processing:$PYTHONPATH python3 gen_ref_state.py
mpirun -np 4 $RAYLEIGH_TEST_MPI_PARAMS ../../../bin/rayleigh.dbg
cd ..


# after both versions have run, we test the output for errors
PYTHONPATH=../../post_processing:../../pre_processing:$PYTHONPATH python3 test_output.py

