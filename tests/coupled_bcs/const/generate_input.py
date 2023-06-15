### Use rayleigh_spectral_input.py to generate generic input.

from rayleigh_spectral_input import *

rmin = 0.5; rmax = 1.0

si = SpectralInput(n_theta=1,n_r=48)
si.transform_from_rtp_function(lambda radius: 2.0 - (rmax/radius)*(radius - rmin)/(rmax - rmin), rmin=rmin, rmax=rmax)
si.write('radial_prof_init')

si = SpectralInput()
si.add_mode(0.0+0.j, 0, 0, 0)
si.write('constant_init')

