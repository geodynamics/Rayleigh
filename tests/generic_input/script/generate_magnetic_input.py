### Use rayleigh_spectral_input.py to generate generic input.
### Christensen Benchmark case 1.

from rayleigh_spectral_input import *

ar = 0.35; sd=1.0
rmin, rmax = radial_extents(aspect_ratio=ar, shell_depth=sd)

### Poloidal magnetic potential BC.
si_BC = SpectralInput(n_theta=64,n_r=48)
def func_BC(theta,phi,radius):
    Br = 5./8*(8*rmax-6*radius-2*rmin**4/radius**3)*np.cos(theta)
    r2Br = radius**2*Br
    return r2Br
si_BC.transform_from_rtp_function(func_BC, aspect_ratio=ar, shell_depth=sd)
# si_BC now contains r**2*B_r, because this is only related to the horizontal gradient of C, 
# we need to scale by l*(l+1) to get the correct potential
for ell in range(1,si_BC.coeffs.shape[1]):
    si_BC.coeffs[:,ell,:] = si_BC.coeffs[:,ell,:]/(ell*(ell+1))
si_BC.write('bench_c_init')

### Toroidal magnetic potential BA.
si_BA = SpectralInput(n_theta=64,n_r=48)
def func_BA(theta,phi,radius):
    P20 = 0.5*(3*np.cos(theta)**2-1)
    rdotcurlB = 5*np.sin(np.pi*(radius-rmin))/radius*P20*4
    return radius**2*rdotcurlB
si_BA.transform_from_rtp_function(func_BA, aspect_ratio=ar, shell_depth=sd)
# si_BA now contains r**2*B_r, because this is only related to the horizontal gradient of A,
# we need to scale by l*(l+1) to get the correct potential
for ell in range(1,si_BA.coeffs.shape[1]):
    si_BA.coeffs[:,ell,:] = si_BA.coeffs[:,ell,:]/(ell*(ell+1))
si_BA.write('bench_a_init')
