"""Rayleigh's velocity poloidal-toroidal decomposition.

Convention:

    rho(r) * v = curl(curl(W(r,theta,phi) * rhat)) + curl(Z(r,theta,phi) * rhat)

For a single mode W(r,theta,phi) = Wl(r) * Y_lm(theta,phi), Z=0 (purely
poloidal), `velocity_from_W` below computes v_r, v_theta, v_phi as sympy
expressions in (r, theta, phi) directly from this relation (matching, for a
single mode, the known reduced form v_r = l(l+1)/(rho*r^2) * Wl(r) * Y_lm 
-- derived here by curling twice, so it generalises to whatever
W is passed in).
"""

import sympy as sp

from sympy_utils import r, theta, phi, curl
from assess_fields import Y


def velocity_from_W(Wl, l, m, rho=1):
    """Given the radial part Wl(r) of a purely poloidal, single (l,m)-mode
    potential (a sympy expression in r), return (v_r, v_theta, v_phi) as
    sympy expressions in (r, theta, phi).
    """
    W = Wl * Y(l, m)
    curlWr, curlWt, curlWp = curl(W, 0, 0)  # first curl of W*rhat
    Pr, Pt, Pp = curl(curlWr, curlWt, curlWp)  # = rho*v
    return sp.simplify(Pr / rho), sp.simplify(Pt / rho), sp.simplify(Pp / rho)


def velocity_field_quantities(Wl, l, m, rho=1):
    """v_r, v_theta, v_phi and their r/theta/phi
    derivatives, as sympy expressions in (r, theta, phi).
    """
    v = dict(zip(('r', 'theta', 'phi'), velocity_from_W(Wl, l, m, rho)))

    out = {'v_r': v['r'], 'v_theta': v['theta'], 'v_phi': v['phi']}
    for comp, f in v.items():
        out[f'dv_{comp}_dr'] = sp.diff(f, r)
        out[f'dv_{comp}_dt'] = sp.diff(f, theta)
        out[f'dv_{comp}_dp'] = sp.diff(f, phi)
        out[f'dv_{comp}_dtr'] = sp.diff(f, theta) / r
        out[f'dv_{comp}_dprs'] = sp.diff(f, phi) / (r * sp.sin(theta))
        out[f'dv_{comp}_d2r'] = sp.diff(f, r, 2)
        out[f'dv_{comp}_d2t'] = sp.diff(f, theta, 2)
        out[f'dv_{comp}_d2p'] = sp.diff(f, phi, 2)
        out[f'dv_{comp}_d2rt'] = sp.diff(f, r, theta)
        out[f'dv_{comp}_d2rp'] = sp.diff(f, r, phi)
        out[f'dv_{comp}_d2tp'] = sp.diff(f, theta, phi)
    return out
