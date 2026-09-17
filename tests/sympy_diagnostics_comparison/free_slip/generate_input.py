"""Build w_init (velocity), p_init (pressure), and t_init (temperature) for
a single-mode assess smooth free-slip poloidal field. Wl is rescaled to
Rayleigh's W convention; the pressure Pl needs no rescaling.

t_init is just the bare spherical harmonic Y(l,m) (radial part = 1):
Rayleigh's own Buoyancy_Coeff(r) = (Ra/Pr)*(r/Rp)**gravity_power already
supplies the r**k radial dependence of assess's density perturbation
(delta_rho ~ r**k * Y(l,m) / Rp**k) once gravity_power is set to k in
main_input -- so no separate radial profile is needed in t_init
itself.
"""

import sys
import os
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'common'))

from assess_fields import SmoothFreeSlipPoloidal
from rayleigh_utils import write_mode_file
from sympy_utils import r

rmin, rmax = 0.5, 1.0
l, m, k = 8, 4, 9
rho = 1.0  # Boussinesq, reference_type=1

sol = SmoothFreeSlipPoloidal(l=l, m=m, k=k, Rp=rmax, Rm=rmin, nu=1.0, g=1.0)

Wl = -rho * r * sol.Pl

write_mode_file('w_init', Wl, l, m, rmin, rmax, n_r=48)
print(f"wrote w_init (l={l}, m={m}, k={k})")

write_mode_file('p_init', sol.Pl_pressure, l, m, rmin, rmax, n_r=48)
print(f"wrote p_init (l={l}, m={m}, k={k})")

write_mode_file('t_init', sp.Integer(1), l, m, rmin, rmax, n_r=48)
print(f"wrote t_init (l={l}, m={m})")
