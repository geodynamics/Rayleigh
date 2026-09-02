"""Tolerance model for the quantity comparisons.

Every quantity matches Rayleigh's output to floating-point precision
(~1e-13 relative) *except* the second-theta-derivative-related velocity and
magnetic quantities, and the viscous force, which lose accuracy near the poles.

POLE_LIMITED_TOLERANCES gives each affected quantity its own
tol(theta) = floor + scale / sin(theta)**2. floor covers each quantity's own
equatorial error; scale covers its pole-amplified growth. Both are 
calibrated per quantity against the observed errors on free_slip and no_slip 
(48x64 grid, l=8,m=4,k=9), with a several-fold safety margin.
"""

import numpy as np

TIGHT_TOL = 1.0e-10

# name -> (floor, scale) for tol(theta) = floor + scale / sin(theta)**2
POLE_LIMITED_TOLERANCES = {
    'dv_theta_d2r':  (3.0e-5, 3.0e-5),
    'dv_phi_d2r':    (2.5e-5, 3.5e-5),
    'dv_r_d2t':      (1.5e-6, 4.0e-6),
    'dv_theta_d2t':  (1.5e-4, 1.7e-4),
    'dv_phi_d2t':    (6.0e-6, 2.7e-5),
    'dv_r_d2rt':     (5.0e-6, 5.0e-6),
    'dv_theta_d2rt': (4.0e-5, 7.0e-5),
    'dv_phi_d2rt':   (6.0e-5, 6.0e-5),
    'dv_theta_d2rp': (4.0e-6, 4.5e-6),
    'dv_phi_d2rp':   (3.0e-6, 5.0e-6),

    'viscous_Force_r':     (1.5e-6, 6.0e-6),
    'viscous_Force_theta': (1.5e-4, 1.6e-4),
    'viscous_Force_phi':   (2.0e-5, 3.5e-5),

    'curl_viscous_force_r':     (3.0e-3, 3.0e-3),
    'curl_viscous_force_theta': (7.0e-4, 7.0e-4),
    'curl_viscous_force_phi':   (2.5e-3, 2.5e-3),

    'curl_v_grad_v_r':     (1.0e-8, 1.0e-8),
    'curl_v_grad_v_theta': (8.0e-9, 8.0e-9),
    'curl_v_grad_v_phi':   (1.5e-8, 1.5e-8),

    'curl_j_cross_b_r':     (5.0e-11, 5.0e-11),
    'curl_j_cross_b_theta': (2.5e-3, 3.5e-3),
    'curl_j_cross_b_phi':   (2.0e-3, 2.0e-3),

    'db_r_d2t':      (1.0e-3, 2.0e-3),
    'db_theta_d2r':  (3.0e-3, 3.0e-3),
    'db_phi_d2r':    (2.5e-2, 2.0e-2),
    'db_r_d2rt':     (7.0e-4, 7.0e-4),

    'db_theta_dr':   (3.5e-8, 3.5e-8),
    'db_r_d2r':      (6.0e-8, 6.0e-8),
    'db_theta_d2t':  (3.0e-11, 3.0e-11),
    'db_phi_d2t':    (5.0e-12, 5.0e-11),
    'db_theta_d2rt': (1.5e-8, 1.5e-8),
    'db_phi_d2rt':   (6.0e-11, 6.0e-11),
}


def tolerance(name, theta):
    """Elementwise error tolerance for quantity `name` at colatitude(s)
    `theta`."""
    if name not in POLE_LIMITED_TOLERANCES:
        return TIGHT_TOL
    floor, scale = POLE_LIMITED_TOLERANCES[name]
    return floor + scale / np.sin(theta)**2
