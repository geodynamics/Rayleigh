"""Shared sympy symbols and operators for the analytic (r, theta, phi) side
of the diagnostics comparisons in this test folder.
"""

import sympy as sp

r, theta, phi = sp.symbols('r theta phi', positive=True)


def curl(Fr, Ft, Fp):
    """Standard spherical curl of a vector field with physical (orthonormal)
    components (Fr, Ft, Fp), each a sympy expression in (r, theta, phi).
    theta = co-latitude, phi = longitude.
    """
    curl_r = (1 / (r * sp.sin(theta))) * (sp.diff(Fp * sp.sin(theta), theta) - sp.diff(Ft, phi))
    curl_t = (1 / r) * (sp.diff(Fr, phi) / sp.sin(theta) - sp.diff(r * Fp, r))
    curl_p = (1 / r) * (sp.diff(r * Ft, r) - sp.diff(Fr, theta))
    return curl_r, curl_t, curl_p


def cross(Ar, At, Ap, Br, Bt, Bp):
    """Cross product A x B of two vector fields with physical (orthonormal)
    components, each a sympy expression in (r, theta, phi)."""
    Cr = At * Bp - Ap * Bt
    Ct = Ap * Br - Ar * Bp
    Cp = Ar * Bt - At * Br
    return Cr, Ct, Cp
