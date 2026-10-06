"""Tolerance model for the quantity comparisons.

Tolerances are relative to the size of each quantity:

    tol = rel_tol(name) * max|analytic| + ABS_FLOOR

where max|analytic| is taken over the whole full3d grid (and reused for the
point probes), rel_tol is DEFAULT_REL_TOL unless the quantity appears in
REL_TOLERANCES, and ABS_FLOOR only matters for quantities that are
analytically (close to) zero, whose errors are roundoff-level.  Being relative,
the tolerances don't need recalibrating when the test fields' amplitudes
change.

Each REL_TOLERANCES entry is ~5x the largest relative error observed on
free_slip and no_slip (48x64 grid; Christensen benchmark magnetic field;
free_slip velocity with poloidal (l=8,m=4) + (l=2,m=0) and toroidal
(l=3,m=2) + (l=4,m=0) modes).

One small theta dependence remains: a few quantities have extra error in the
one or two grid rows nearest the poles, at every radius rather than just at the
boundaries.  It is small enough to fit under a flat tolerance, so no
theta-dependent term is used:
  - db_theta_d2t/db_phi_d2t (relative ~1e-9 at the poles, ~1e-11 elsewhere):
    ~1e-12 radial-truncation errors in the first-derivative inputs are divided
    by sin(theta) (or sin^2) in the d2t formulas.
  - db_r_d2t (similar size): Lap_1 multiplies the transform's flat,
    roundoff-level spectral noise (~2e-13 relative) by l(l+1), up to ~1700 at
    l_max, and for this m=0 field the high-l harmonics peak at the poles.
  - dv_phi_d2t and curl_v_grad_v_r (~1e-10 relative locally): the same kind of
    excess from free_slip's axisymmetric toroidal velocity.  The
    non-axisymmetric velocity modes vanish at the poles like sin^m(theta) and
    don't show it.
"""

# Most quantities match to <~2e-11 relative.
DEFAULT_REL_TOL = 1.0e-10

# Absolute floor, for quantities that are analytically ~zero.
ABS_FLOOR = 1.0e-12

# name -> relative tolerance (observed max relative error in comments)
REL_TOLERANCES = {
    # Second radial derivatives of the velocity, and the viscous force built
    # from them.
    'dv_theta_d2r':        1.5e-8,   # 2.6e-9
    'dv_phi_d2r':          1.5e-8,   # 2.6e-9
    'dv_r_d2r':            7.0e-10,  # 1.3e-10
    'dv_r_d2t':            1.5e-10,  # 2.3e-11
    'viscous_Force_r':     3.0e-10,  # 5.7e-11
    'viscous_Force_theta': 2.0e-9,   # 4.1e-10
    'viscous_Force_phi':   2.0e-9,   # 4.1e-10

    # The curl of the viscous force needs third radial derivatives.
    'curl_viscous_force_r':      1.0e-9,  # 1.5e-10
    'curl_viscous_force_theta':  1.0e-6,  # 1.8e-7
    'curl_viscous_force_phi':    1.0e-6,  # 1.7e-7
    'curl_viscous_pforce_r':     1.0e-9,  # 1.9e-10
    'curl_viscous_pforce_theta': 1.0e-6,  # 1.8e-7
    'curl_viscous_pforce_phi':   1.0e-6,  # 1.7e-7
    'curl_viscous_mforce_r':     1.0e-9,  # 1.4e-10
    'curl_viscous_mforce_theta': 1.0e-6,  # 1.8e-7
    'curl_viscous_mforce_phi':   1.0e-6,  # 1.0e-7
    'curl_viscous_force_r_squared':     1.5e-9,  # 3.0e-10
    'curl_viscous_force_theta_squared': 1.5e-6,  # 2.5e-7
    'curl_viscous_force_phi_squared':   1.5e-6,  # 2.4e-7
    'curl_viscous_force_abs':           5.0e-8,  # 9.1e-9

    # Magnetic field radial derivatives: the benchmark field's radial profile
    # is not resolved to machine precision by the Chebyshev expansion.
    'db_theta_dr':   5.0e-9,   # 8.7e-10
    'db_r_d2r':      4.0e-9,   # 7.4e-10
    'db_theta_d2r':  1.5e-6,   # 3.0e-7
    'db_phi_d2r':    1.0e-9,   # 1.6e-10
    'db_r_d2rt':     2.0e-10,  # 3.0e-11
    'db_theta_d2rt': 5.0e-9,   # 8.8e-10

    # Second theta derivatives of B (see the note on the poles above).
    'db_r_d2t':      5.0e-9,   # 1.1e-9
    'db_theta_d2t':  3.0e-10,  # 5.7e-11
    'db_phi_d2t':    2.0e-10,  # 3.0e-11

    # Lorentz force and its curl inherit curl(B)'s radial-derivative floor.
    'j_cross_b_r':          2.0e-9,  # 3.4e-10
    'j_cross_b_theta':      1.0e-9,  # 1.5e-10
    'curl_j_cross_b_theta': 1.5e-9,  # 2.5e-10
    'curl_j_cross_b_phi':   2.0e-7,  # 4.5e-8
}


def tolerance(name, scale):
    """Absolute error tolerance for quantity `name`, whose analytic values have
    maximum magnitude `scale`."""
    return REL_TOLERANCES.get(name, DEFAULT_REL_TOL) * scale + ABS_FLOOR
