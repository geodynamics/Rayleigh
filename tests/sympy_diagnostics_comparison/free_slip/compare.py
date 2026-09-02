"""Compare Rayleigh's output (full3d and Point_Probes) against sympy reconstructions of the same quantities.
"""

import argparse
import sys
import os
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'common'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'post_processing'))

from assess_fields import SmoothFreeSlipPoloidal, Y
from velocity_field import velocity_field_quantities, velocity_from_W
from velocity_field_codes import QUANTITY_CODES as VELOCITY_CODES
from momentum_forces import v_grad_v_force, coriolis_force, viscous_force, pressure_force, buoyancy_force
from momentum_force_codes import QUANTITY_CODES as MOMENTUM_CODES
from curl_momentum_forces import curl_v_grad_v_force, curl_buoyancy_force, curl_coriolis_force, curl_pressure_force, curl_viscous_force
from curl_momentum_force_codes import QUANTITY_CODES as CURL_MOMENTUM_CODES
from compare_utils import compare_full3d_and_probes
from sympy_utils import r, theta, phi

parser = argparse.ArgumentParser()
parser.add_argument('--vtu', action='store_true',
                     help='also write a Paraview-loadable VTU (rayleigh_/analytic_/diff_ per quantity, '
                          'see common/vtu_utils.py) from the exact same quantity_codes/numeric used below')
args = parser.parse_args()

rmin, rmax = 0.5, 1.0
l, m, k = 8, 4, 9
rho = 1.0
Ekman_Number = 1.0
coriolis_coeff = 2.0 / Ekman_Number   # ref%Coriolis_Coeff, reference_type=1
pfactor = 1.0 / Ekman_Number          # ref%dpdr_w_term/ref%density, rotation=.true.
Rayleigh_Number, Prandtl_Number = 1.0, 1.0
gravity_power = k                     # matches assess's r**k density perturbation
buoyancy_coeff = (Rayleigh_Number / Prandtl_Number) * (r / rmax)**gravity_power  # ref%Buoyancy_Coeff(r)

sol = SmoothFreeSlipPoloidal(l=l, m=m, k=k, Rp=rmax, Rm=rmin, nu=1.0, g=1.0)

Wl = -rho * r * sol.Pl
quantities = velocity_field_quantities(Wl, l, m, rho)
numeric = {name: sp.lambdify((r, theta, phi), expr, 'numpy') for name, expr in quantities.items()}

vr, vt, vp = velocity_from_W(Wl, l, m, rho)
numeric['v_grad_v_r'], numeric['v_grad_v_theta'], numeric['v_grad_v_phi'] = v_grad_v_force(vr, vt, vp, rho)
numeric['Coriolis_Force_r'], numeric['Coriolis_Force_theta'], numeric['Coriolis_Force_phi'] = \
    coriolis_force(vr, vt, vp, coriolis_coeff, rho)
numeric['viscous_Force_r'], numeric['viscous_Force_theta'], numeric['viscous_Force_phi'] = \
    viscous_force(vr, vt, vp, mu_visc=1.0)
numeric['pressure_Force_r'], numeric['pressure_Force_theta'], numeric['pressure_Force_phi'] = \
    pressure_force(sol.P, pfactor)

Theta = Y(l, m)  # t_init radial part = 1
numeric['buoyancy_force'] = buoyancy_force(Theta, buoyancy_coeff)
numeric['curl_v_grad_v_r'], numeric['curl_v_grad_v_theta'], numeric['curl_v_grad_v_phi'] = \
    curl_v_grad_v_force(vr, vt, vp, rho)
numeric['curl_buoyancy_force_theta'], numeric['curl_buoyancy_force_phi'] = \
    curl_buoyancy_force(Theta, buoyancy_coeff)
numeric['curl_coriolis_force_r'], numeric['curl_coriolis_force_theta'], numeric['curl_coriolis_force_phi'] = \
    curl_coriolis_force(vr, vt, vp, coriolis_coeff, rho)
numeric['curl_pressure_force_theta'], numeric['curl_pressure_force_phi'] = \
    curl_pressure_force(sol.P, pfactor)
numeric['curl_viscous_force_r'], numeric['curl_viscous_force_theta'], numeric['curl_viscous_force_phi'] = \
    curl_viscous_force(vr, vt, vp, mu_visc=1.0)

quantity_codes = dict(VELOCITY_CODES)
quantity_codes.update(MOMENTUM_CODES)
quantity_codes.update(CURL_MOMENTUM_CODES)

ok = compare_full3d_and_probes(quantity_codes, numeric)

if args.vtu:
    from vtu_utils import write_comparison_vtu
    write_comparison_vtu(quantity_codes, numeric, out_path='Spherical_3D/compare_velocity.vtu')

if not ok:
    print("\nERROR: one or more quantities exceeded their tolerance.")
    sys.exit(1)
print("\nPASS")
sys.exit(0)
