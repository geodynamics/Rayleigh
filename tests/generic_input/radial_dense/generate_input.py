### Use rayleigh_spectral_input.py to generate generic input.
### Christensen Benchmark case 1.

from rayleigh_spectral_input import *

rmin = 0.5; rmax = 1.0

si = SpectralInput(n_theta=64,n_r=48)
si.transform_from_rtp_function(lambda radius: 1.0 - (rmax/radius)*(radius - rmin)/(rmax - rmin), rmin=rmin, rmax=rmax)
si.write('radial_t_init')
