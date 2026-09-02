"""Build w_init (velocity), p_init (pressure), t_init (temperature), and
c_init/a_init (magnetic poloidal/toroidal potentials) for a single-mode
assess smooth zero-slip poloidal field. Wl is rescaled to Rayleigh's W
convention; the pressure Pl needs no rescaling.

t_init is just the bare spherical harmonic Y(l,m) (radial part = 1):
Rayleigh's own Buoyancy_Coeff(r) = (Ra/Pr)*(r/Rp)**gravity_power already
supplies the r**k radial dependence of assess's density perturbation
(delta_rho ~ r**k * Y(l,m) / Rp**k) once gravity_power is set to k in
main_input -- so no separate radial profile is needed in t_init
itself.

c_init/a_init use the exact Christensen et al. 2001 benchmark 1 insulating
field (l_C,m_C=1,0 poloidal; l_A,m_A=2,0 toroidal).
"""

import sys
import os
import sympy as sp
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'common'))

from assess_fields import SmoothZeroSlipPoloidal
from rayleigh_utils import write_mode_file
from magnetic_field import benchmark_insulating_radial
from sympy_utils import r

rmin, rmax = 0.5, 1.0
l, m, k = 8, 4, 9
rho = 1.0  # Boussinesq, reference_type=1
l_C, m_C = 1, 0   # poloidal magnetic potential mode -- matches magnetic_init_type=1
l_A, m_A = 2, 0   # toroidal magnetic potential mode -- matches magnetic_init_type=1

sol = SmoothZeroSlipPoloidal(l=l, m=m, k=k, Rp=rmax, Rm=rmin, nu=1.0, g=1.0)

Wl = -rho * r * sol.Pl

write_mode_file('w_init', Wl, l, m, rmin, rmax, n_r=48)
print(f"wrote w_init (l={l}, m={m}, k={k})")

write_mode_file('p_init', sol.Pl_pressure, l, m, rmin, rmax, n_r=48)
print(f"wrote p_init (l={l}, m={m}, k={k})")

write_mode_file('t_init', sp.Integer(1), l, m, rmin, rmax, n_r=48)
print(f"wrote t_init (l={l}, m={m})")

rfunc1, rfunc2 = benchmark_insulating_radial(rmin, rmax)

write_mode_file('c_init', rfunc1, l_C, m_C, rmin, rmax, n_r=48)
print(f"wrote c_init (l={l_C}, m={m_C})")

write_mode_file('a_init', rfunc2, l_A, m_A, rmin, rmax, n_r=48)
print(f"wrote a_init (l={l_A}, m={m_A})")
