"""Compare Rayleigh's output (full3d and Point_Probes) against sympy reconstructions of the same quantities.
"""

import argparse
import sys
import os
import sympy as sp

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'common'))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'post_processing'))

from assess_fields import SmoothFreeSlipPoloidal, Y
from velocity_field import velocity_quantities_from_v, velocity_from_W, velocity_from_Z
from velocity_field_codes import QUANTITY_CODES as VELOCITY_CODES
from momentum_forces import v_grad_v_force, coriolis_force, viscous_force, pressure_force, buoyancy_force
from momentum_force_codes import QUANTITY_CODES as MOMENTUM_CODES
from curl_momentum_forces import curl_v_grad_v_force, curl_buoyancy_force, curl_coriolis_force, curl_pressure_force, curl_viscous_force
from curl_momentum_force_codes import QUANTITY_CODES as CURL_MOMENTUM_CODES, MEAN_FLUCTUATING_CODES, \
    CURL_VISCOUS_DERIVED_CODES
from compare_utils import compare_full3d_and_probes
from sympy_utils import r, theta, phi

parser = argparse.ArgumentParser()
parser.add_argument('--vtu', action='store_true',
                     help='also write a Paraview-loadable VTU (rayleigh_/analytic_/diff_ per quantity, '
                          'see common/vtu_utils.py) from the exact same quantity_codes/numeric used below')
args = parser.parse_args()

print("=================")
print("testing free_slip")
print("=================")

rmin, rmax = 0.5, 1.0
l, m, k = 8, 4, 9
rho = 1.0
Ekman_Number = 1.0
coriolis_coeff = 2.0 / Ekman_Number   # ref%Coriolis_Coeff, reference_type=1
pfactor = 1.0 / Ekman_Number          # ref%dpdr_w_term/ref%density, rotation=.true.
Rayleigh_Number, Prandtl_Number = 1.0, 1.0
gravity_power = k                     # matches assess's r**k density perturbation
buoyancy_coeff = (Rayleigh_Number / Prandtl_Number) * (r / rmax)**gravity_power  # ref%Buoyancy_Coeff(r)

l0, m0 = 2, 0   # axisymmetric velocity mode (see generate_input.py)

sol = SmoothFreeSlipPoloidal(l=l, m=m, k=k, Rp=rmax, Rm=rmin, nu=1.0, g=1.0)
sol0 = SmoothFreeSlipPoloidal(l=l0, m=m0, k=k, Rp=rmax, Rm=rmin, nu=1.0, g=1.0)

Wl = -rho * r * sol.Pl
Wl0 = -rho * r * sol0.Pl

# toroidal modes (see generate_input.py)
lz, mz = 3, 2
lz0, mz0 = 4, 0
Zl = rho * r**2 * (1 + 2*r**3 - 3*(rmin + rmax)*r**2 + 6*rmin*rmax*r)

def vsum(*vs):
    return tuple(sum(c) for c in zip(*vs))

# fluctuating (m != 0) and mean (m = 0) parts of the velocity
vr1, vt1, vp1 = vsum(velocity_from_W(Wl, l, m, rho), velocity_from_Z(Zl, lz, mz, rho))
vr0, vt0, vp0 = vsum(velocity_from_W(Wl0, l0, m0, rho), velocity_from_Z(Zl/2, lz0, mz0, rho))
vr, vt, vp = vr1 + vr0, vt1 + vt0, vp1 + vp0

quantities = velocity_quantities_from_v(vr, vt, vp)
numeric = {name: sp.lambdify((r, theta, phi), expr, 'numpy') for name, expr in quantities.items()}

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
cvf = curl_viscous_force(vr, vt, vp, mu_visc=1.0)
numeric['curl_viscous_force_r'], numeric['curl_viscous_force_theta'], numeric['curl_viscous_force_phi'] = cvf
for comp, f in zip(('r', 'theta', 'phi'), cvf):
    numeric[f'curl_viscous_force_{comp}_squared'] = (lambda g: lambda *x: g(*x)**2)(f)
numeric['curl_viscous_force_abs'] = lambda *x: (cvf[0](*x)**2 + cvf[1](*x)**2 + cvf[2](*x)**2)**0.5
numeric['curl_viscous_pforce_r'], numeric['curl_viscous_pforce_theta'], numeric['curl_viscous_pforce_phi'] = \
    curl_viscous_force(vr1, vt1, vp1, mu_visc=1.0)
numeric['curl_viscous_mforce_r'], numeric['curl_viscous_mforce_theta'], numeric['curl_viscous_mforce_phi'] = \
    curl_viscous_force(vr0, vt0, vp0, mu_visc=1.0)

quantity_codes = dict(VELOCITY_CODES)
quantity_codes.update(MOMENTUM_CODES)
quantity_codes.update(CURL_MOMENTUM_CODES)
quantity_codes.update(MEAN_FLUCTUATING_CODES)
quantity_codes.update(CURL_VISCOUS_DERIVED_CODES)

ok = compare_full3d_and_probes(quantity_codes, numeric)

if args.vtu:
    from vtu_utils import write_comparison_vtu
    write_comparison_vtu(quantity_codes, numeric, out_path='Spherical_3D/compare_velocity.vtu')

if not ok:
    print("\nERROR: one or more quantities exceeded their tolerance.")
    sys.exit(1)
print("\nPASS")
sys.exit(0)
