"""Sympy reconstructions of the curl-momentum-equation forces in
src/Diagnostics/curl_momentum_equation_codes.F, built from the
sympy_utils.curl() operator.
"""

import sympy as sp

from sympy_utils import r, theta, phi, curl, cross
from momentum_forces import ADotGradB


def curl_v_grad_v_force(vr, vt, vp, rho=1):
    """curl_v_grad_v_r, curl_v_grad_v_theta, curl_v_grad_v_phi as
    numpy-callables f(r,theta,phi). Built from momentum_forces.py's
    ADotGradB(v,v,v) (the (v.grad)v advection term), scaled by rho and
    curled. momentum_forces.py's v_grad_v_force's ell=0 subtraction (on its
    r-component only) is *not* needed here: that correction is a function of
    r alone, which has identically zero curl in every component (r, theta,
    and phi), so this is curl(rho*ADotGradB(v,v,v)), not curl of 
    v_grad_v_r/theta/phi.
    """
    Cr, Ct, Cp = ADotGradB(vr, vt, vp, vr, vt, vp)
    Fr, Ft, Fp = rho * Cr, rho * Ct, rho * Cp
    curl_r, curl_t, curl_p = curl(Fr, Ft, Fp)
    return (sp.lambdify((r, theta, phi), curl_r, 'numpy'),
            sp.lambdify((r, theta, phi), curl_t, 'numpy'),
            sp.lambdify((r, theta, phi), curl_p, 'numpy'))


def curl_buoyancy_force(Theta, buoyancy_coeff):
    """curl_buoyancy_force_theta, curl_buoyancy_force_phi as numpy-callables
    f(r,theta,phi). `Theta` is the temperature field (a sympy
    expression in r,theta,phi); `buoyancy_coeff` is a sympy expression of r). 
    curl_buoyancy_force_r (always 0) is discarded.
    """
    F = buoyancy_coeff * Theta
    _, curl_t, curl_p = curl(F, 0, 0)
    return sp.lambdify((r, theta, phi), curl_t, 'numpy'), sp.lambdify((r, theta, phi), curl_p, 'numpy')


def curl_pressure_force(P, pfactor, dlnrho=0):
    """curl_pressure_force_theta, curl_pressure_force_phi as numpy-callables
    f(r,theta,phi). curl_pressure_force_r (always 0) is discarded.
    """
    Fr = -pfactor * sp.diff(P, r) + pfactor * dlnrho * P
    Ft = -pfactor * sp.diff(P, theta) / r
    Fp = -pfactor * sp.diff(P, phi) / (r * sp.sin(theta))
    _, curl_t, curl_p = curl(Fr, Ft, Fp)
    return sp.lambdify((r, theta, phi), curl_t, 'numpy'), sp.lambdify((r, theta, phi), curl_p, 'numpy')


def curl_viscous_force(vr, vt, vp, mu_visc=1):
    """curl_viscous_force_r, curl_viscous_force_theta, curl_viscous_force_phi
    as numpy-callables f(r,theta,phi): curl(Fr,Ft,Fp) where (Fr,Ft,Fp) is 
    the full viscous_Force vector from momentum_forces.py's viscous_force() 
    (mu_visc * -curl(curl(v)), valid for a constant density/viscosity Boussinesq 
    reference state.
    """
    curlvr, curlvt, curlvp = curl(vr, vt, vp)
    lap_r, lap_t, lap_p = curl(curlvr, curlvt, curlvp)
    Fr, Ft, Fp = mu_visc * -lap_r, mu_visc * -lap_t, mu_visc * -lap_p
    curl_r, curl_t, curl_p = curl(Fr, Ft, Fp)
    return (sp.lambdify((r, theta, phi), curl_r, 'numpy'),
            sp.lambdify((r, theta, phi), curl_t, 'numpy'),
            sp.lambdify((r, theta, phi), curl_p, 'numpy'))


def curl_coriolis_force(vr, vt, vp, coriolis_coeff, rho=1):
    """curl_coriolis_force_r, curl_coriolis_force_theta,
    curl_coriolis_force_phi as numpy-callables f(r,theta,phi), for a
    constant-density (dlnrho=0) reference state.
    """
    C = coriolis_coeff
    Fr = rho * C * sp.sin(theta) * vp
    Ft = rho * C * sp.cos(theta) * vp
    Fp = -rho * C * (sp.cos(theta) * vt + sp.sin(theta) * vr)
    curl_r, curl_t, curl_p = curl(Fr, Ft, Fp)
    return (sp.lambdify((r, theta, phi), curl_r, 'numpy'),
            sp.lambdify((r, theta, phi), curl_t, 'numpy'),
            sp.lambdify((r, theta, phi), curl_p, 'numpy'))


def curl_j_cross_b_force(Br, Bt, Bp, lorentz_coeff):
    """curl_j_cross_b_r, curl_j_cross_b_theta, curl_j_cross_b_phi as
    numpy-callables f(r,theta,phi).
    """
    Jr, Jt, Jp = curl(Br, Bt, Bp)
    Fr, Ft, Fp = tuple(lorentz_coeff * c for c in cross(Jr, Jt, Jp, Br, Bt, Bp))
    curl_r, curl_t, curl_p = curl(Fr, Ft, Fp)
    return (sp.lambdify((r, theta, phi), curl_r, 'numpy'),
            sp.lambdify((r, theta, phi), curl_t, 'numpy'),
            sp.lambdify((r, theta, phi), curl_p, 'numpy'))
