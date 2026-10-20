"""Build w_init (velocity), p_init (pressure), and t_init (temperature) for
a single-mode assess smooth free-slip poloidal field. Wl is rescaled to
Rayleigh's W convention; the pressure Pl needs no rescaling.

t_init is just the spherical harmonic Y(l,m) (radial part = 1):
Rayleigh's own Buoyancy_Coeff(r) = (Ra/Pr)*(r/Rp)**gravity_power already
supplies the r**k radial dependence of assess's density perturbation
(delta_rho ~ r**k * Y(l,m) / Rp**k) once gravity_power is set to k in
main_input -- no separate radial profile is needed in t_init.
"""

import sys
import os
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'common'))

from assess_fields import SmoothFreeSlipPoloidal
from rayleigh_utils import write_mode_file, write_modes_file
from sympy_utils import r

rmin, rmax = 0.5, 1.0
l, m, k = 8, 4, 9
rho = 1.0  # Boussinesq, reference_type=1

# The velocity also has an axisymmetric (m=0) mode so that the mean (m=0) and
# fluctuating parts of the viscous force curl are both nontrivial.
l0, m0 = 2, 0

sol = SmoothFreeSlipPoloidal(l=l, m=m, k=k, Rp=rmax, Rm=rmin, nu=1.0, g=1.0)
sol0 = SmoothFreeSlipPoloidal(l=l0, m=m0, k=k, Rp=rmax, Rm=rmin, nu=1.0, g=1.0)

Wl = -rho * r * sol.Pl
Wl0 = -rho * r * sol0.Pl

write_modes_file('w_init', [(Wl, l, m), (Wl0, l0, m0)], rmin, rmax, n_r=48)
print(f"wrote w_init (l={l}, m={m}, k={k}) + (l={l0}, m={m0}, k={k})")

# Toroidal velocity, so that the radial curl of the viscous force is nonzero:
# a non-axisymmetric (lz, mz) mode and an axisymmetric (lz0, 0) mode, sharing a
# polynomial radial profile (exactly representable in Chebyshev) that satisfies
# the free-slip condition d(Z/r^2)/dr = 0 at both boundaries.
lz, mz = 3, 2
lz0, mz0 = 4, 0
Zl = rho * r**2 * (1 + 2*r**3 - 3*(rmin + rmax)*r**2 + 6*rmin*rmax*r)
write_modes_file('z_init', [(Zl, lz, mz), (Zl/2, lz0, mz0)], rmin, rmax, n_r=48)
print(f"wrote z_init (l={lz}, m={mz}) + (l={lz0}, m={mz0})")

write_mode_file('p_init', sol.Pl_pressure, l, m, rmin, rmax, n_r=48)
print(f"wrote p_init (l={l}, m={m}, k={k})")

write_mode_file('t_init', sp.Integer(1), l, m, rmin, rmax, n_r=48)
print(f"wrote t_init (l={l}, m={m})")
