"""Sympy reconstructions of the full-field momentum-equation forces.

Assumed Boussinesq (reference_type=1) with constant density (rho=1, dlnrho=0)
and constant viscosity (nu_type=1 default -> dmudr=0):

  v_grad_v  = rho * (v . grad) v                         (ADotGradB(v, v))
  Coriolis  = Coriolis_Coeff * rho * (-zhat x v)
  viscous   = mu_visc * (vector Laplacian of v)
            = mu_visc * (-curl(curl(v)))                 (exact for div(v)=0)
  pressure  = -pfactor * grad(P)                          (dlnrho=0 drops the
                                                             density-gradient term)
  buoyancy  = buoyancy_coeff(r) * Theta                    (purely radial)
  j_cross_b = lorentz_coeff * (curl(B) x B)                 (requires magnetism=.true.
                                                             and a B-field reconstruction,
                                                             see magnetic_field.py)

The *_r components of v_grad_v, Coriolis, and pressure (and the buoyancy
force, which is purely radial) are defined in Rayleigh with their ell=0 
(spherical-mean) part subtracted off. For a single
non-axisymmetric (m != 0) mode, that ell=0 part is zero for any
quantity *linear* in v, P, or Theta (integrating cos(m*phi)/sin(m*phi) over
a full period is exactly zero for m != 0, regardless of any theta-dependent
weighting) -- true for Coriolis_Force_r, viscous_Force_r, pressure_Force_r,
and buoyancy_force, so they are returned as plain sympy expressions.
v_grad_v is *quadratic* in v, so a pure m mode can self-interact to produce
m=0 (axisymmetric) content, which v_grad_v_r has subtracted -- computed via symbolic
phi-integration (kills every m!=0 term exactly) followed by quadrature over cos(theta).
j_cross_b is likewise quadratic (in B), so its r-component gets the same
ell=0 subtraction, whether or not B itself happens to be axisymmetric.
That correction has no closed form, so *every* function below
returns plain numpy-callables f(r, theta, phi) (each
already lambdified internally) rather than sympy expressions, for a
consistent interface.
"""

import numpy as np
import sympy as sp

from sympy_utils import r, theta, phi, curl, cross


def ADotGradB(Ar, At, Ap, Br, Bt, Bp):
    """(A . grad) B in spherical coordinates."""
    dBrdr, dBrdt, dBrdp = sp.diff(Br, r), sp.diff(Br, theta), sp.diff(Br, phi)
    dBtdr, dBtdt, dBtdp = sp.diff(Bt, r), sp.diff(Bt, theta), sp.diff(Bt, phi)
    dBpdr, dBpdt, dBpdp = sp.diff(Bp, r), sp.diff(Bp, theta), sp.diff(Bp, phi)
    oor = 1 / r
    csc = 1 / sp.sin(theta)
    cot = sp.cos(theta) / sp.sin(theta)
    Cr = Ar * dBrdr + oor * (At * (dBrdt - Bt) + Ap * (csc * dBrdp - Bp))
    Ct = Ar * dBtdr + oor * (At * (dBtdt + Br) + Ap * (csc * dBtdp - cot * Bp))
    Cp = Ar * dBpdr + oor * (At * dBpdt + Ap * (csc * dBpdp + Br + cot * Bt))
    return Cr, Ct, Cp


def _lambdify3(expr):
    return sp.lambdify((r, theta, phi), expr, 'numpy')


def _remove_ell0(expr, n_theta=64):
    """Return a numpy-callable f(r, theta, phi) = expr - <expr>_ell0(r),
    where <.>_ell0(r) is the spherical average at each radius. If expr's 
    phi-average is identically zero, this reduces to a plain
    lambdify with no quadrature overhead.
    """
    phi_avg = sp.integrate(expr, (phi, 0, 2 * sp.pi)) / (2 * sp.pi)
    f_full = _lambdify3(expr)
    if phi_avg == 0:
        return f_full

    f_avg = sp.lambdify((r, theta), phi_avg, 'numpy')
    costheta, weights = np.polynomial.legendre.leggauss(n_theta)
    thetas = np.arccos(costheta)

    def ell0_of_r(rval):
        # ell0 depends on r alone; a full3d grid repeats each r value many
        # times (once per theta,phi point), so make unique before the quadrature
        # rather than repeating.
        shape = rval.shape
        unique_r, inverse = np.unique(rval, return_inverse=True)
        out_unique = np.empty(unique_r.size)
        for i, rv in enumerate(unique_r):
            out_unique[i] = 0.5 * np.sum(f_avg(np.full_like(thetas, rv), thetas) * weights)
        return out_unique[inverse].reshape(shape)

    def numeric(rv, tv, pv):
        rv, tv, pv = np.broadcast_arrays(np.asarray(rv, dtype='float64'), tv, pv)
        return f_full(rv, tv, pv) - ell0_of_r(rv)

    return numeric


def v_grad_v_force(vr, vt, vp, rho=1, n_theta=64):
    """v_grad_v_r, v_grad_v_theta, v_grad_v_phi as numpy-callables f(r,theta,phi)."""
    Fr, Ft, Fp = ADotGradB(vr, vt, vp, vr, vt, vp)
    Fr, Ft, Fp = rho * Fr, rho * Ft, rho * Fp
    return _remove_ell0(Fr, n_theta), _lambdify3(Ft), _lambdify3(Fp)


def coriolis_force(vr, vt, vp, coriolis_coeff, rho=1):
    """Coriolis_Force_r, Coriolis_Force_theta, Coriolis_Force_phi as
    numpy-callables f(r,theta,phi). (The r-component's ell=0
    subtraction is not necessary for a single m!=0 mode so is not applied.)"""
    C = coriolis_coeff
    Fr = rho * C * sp.sin(theta) * vp
    Ft = rho * C * sp.cos(theta) * vp
    Fp = -rho * C * (sp.cos(theta) * vt + sp.sin(theta) * vr)
    return _lambdify3(Fr), _lambdify3(Ft), _lambdify3(Fp)


def viscous_force(vr, vt, vp, mu_visc=1):
    """viscous_Force_r, viscous_Force_theta, viscous_Force_phi as
    numpy-callables f(r,theta,phi). For constant density and constant
    viscosity (dmudr=0, dlnrho=0). Since v is assumed
    divergence-free (v = curl(curl(W rhat))/rho, and div(curl(...))
    is identically zero), the vector Laplacian identity
    grad(div v) - curl(curl v) reduces to -curl(curl(v)) exactly. 
    (The r-component's ell=0 subtraction is not necessary for a single m!=0 mode.)
    """
    curlvr, curlvt, curlvp = curl(vr, vt, vp)
    lap_r, lap_t, lap_p = curl(curlvr, curlvt, curlvp)
    Fr, Ft, Fp = mu_visc * -lap_r, mu_visc * -lap_t, mu_visc * -lap_p
    return _lambdify3(Fr), _lambdify3(Ft), _lambdify3(Fp)


def pressure_force(P, pfactor):
    """pressure_Force_r, pressure_Force_theta, pressure_Force_phi as
    numpy-callables f(r,theta,phi) (with dlnrho=0, dropping the density-gradient
    term in the r-component). (The r-component's ell=0 subtraction is
    not necessary for a single m!=0 mode.)
    """
    Fr = -pfactor * sp.diff(P, r)
    Ft = -pfactor * sp.diff(P, theta) / r
    Fp = -pfactor * sp.diff(P, phi) / (r * sp.sin(theta))
    return _lambdify3(Fr), _lambdify3(Ft), _lambdify3(Fp)


def j_cross_b_force(Br, Bt, Bp, lorentz_coeff, n_theta=64):
    """j_cross_b_r, j_cross_b_theta, j_cross_b_phi as numpy-callables
    f(r,theta,phi): lorentz_coeff * (curl(B) x B). j_cross_b is quadratic
    in B (like v_grad_v is quadratic in v), so its r-component gets the
    same ell=0 (spherical-mean) subtraction as v_grad_v_force's.
    """
    Jr, Jt, Jp = curl(Br, Bt, Bp)
    Fr, Ft, Fp = tuple(lorentz_coeff * c for c in cross(Jr, Jt, Jp, Br, Bt, Bp))
    return _remove_ell0(Fr, n_theta), _lambdify3(Ft), _lambdify3(Fp)


def buoyancy_force(Theta, buoyancy_coeff):
    """buoyancy_force as a numpy-callable f(r,theta,phi): c2*f2(r)*Theta. 
    `Theta` is the temperature field (a sympy expression in r,theta,phi); 
    `buoyancy_coeff` is ref%Buoyancy_Coeff(r) (a sympy expression in r alone). 
    (The ell=0 subtraction is not necessary for a single m!=0 mode.)
    """
    F = buoyancy_coeff * Theta
    return _lambdify3(F)
