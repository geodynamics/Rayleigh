"""Rayleigh's magnetic poloidal-toroidal decomposition and the exact radial
profiles used by magnetic_init_type=1.

Convention:

    B = curl(curl(C(r,theta,phi) * rhat)) + curl(A(r,theta,phi) * rhat)

-- the same poloidal/toroidal decomposition as velocity's
rho*v = curl(curl(W rhat)) + curl(Z rhat) (see velocity_field.py), but for
B there is no density factor, and (unlike our purely-poloidal velocity test
fields, which had Z=0) the benchmark field is both poloidal (C)
*and* toroidal (A).
"""

import sympy as sp

from sympy_utils import r, theta, phi, curl
from assess_fields import Y


def benchmark_insulating_radial(rmin, rmax):
    """The radial profiles rfunc1 (l=1,m=0, poloidal C) and rfunc2
    (l=2,m=0, toroidal A) used by magnetic_init_type=1
    (Benchmark_Insulating_Init), transcribed from
    Initial_Conditions.F90's rfunc1/rfunc2. These satisfy Rayleigh's
    default potential-field boundary conditions. Returns (rfunc1, rfunc2) as 
    radial-only sympy expressions.
    """
    r_outer, r_inner = rmax, rmin

    nrm1 = sp.Rational(5, 8) * sp.sqrt(sp.pi / 3)
    nrm2 = sp.sqrt(sp.pi / 5) * sp.Rational(4, 3) * 5

    rfunc1 = nrm1 * r**2 * (8 * r_outer - 6 * r - 2 * r_inner**4 / r**3)
    rfunc2 = nrm2 * r * sp.sin(sp.pi * (r - r_inner))
    return rfunc1, rfunc2


def benchmark_insulating_CA(rmin, rmax):
    """C (l=1,m=0, poloidal) and A (l=2,m=0, toroidal) fields used
    by magnetic_init_type=1 (Benchmark_Insulating_Init). Returns (C, A) as full 
    sympy expressions in (r, theta, phi) -- radial profile times the 
    appropriate Y_lm.
    """
    rfunc1, rfunc2 = benchmark_insulating_radial(rmin, rmax)
    C = rfunc1 * Y(1, 0)
    A = rfunc2 * Y(2, 0)
    return C, A


def generic_radial_profile(l, rmin, rmax, amp=1.0):
    """A smooth, arbitrary radial profile on [rmin, rmax] for use as either
    C_l(r) or A_l(r), for any l/m independent of the velocity solution's mode
    (unlike benchmark_insulating_CA, which is pinned to l=1/l=2).
    """
    return amp * (r - rmin) * (rmax - r) * r**l


def poloidal_potential_profile(l, rmin, rmax, amp=1.0):
    """C_l(r) radial profile for the poloidal magnetic potential, satisfying
    Rayleigh's default magnetic boundary conditions -- matching
    a potential (insulating) field at both boundaries:

        dC/dr + (l/rmax)*C = 0        at r = rmax (outer)
        dC/dr - ((l+1)/rmin)*C = 0    at r = rmin (inner)

    Constructed as an arbitrary smooth bump (same one used by
    generic_radial_profile) corrected by a combination of r**l and
    r**(-l-1) -- the two homogeneous solutions of the boundary operators
    above -- chosen (by solving the resulting 2x2 linear system) to cancel
    the bump's own boundary residuals exactly, the same way
    assess_fields.py's solutions are built to satisfy velocity's free-slip/
    no-slip boundary conditions rather than being arbitrary.
    """
    base = generic_radial_profile(l, rmin, rmax, amp)
    alpha, beta = sp.symbols('alpha beta')
    Cl = base + alpha * r**l + beta * r**(-l - 1)

    outer_bc = sp.diff(Cl, r).subs(r, rmax) + sp.nsimplify(l) / rmax * Cl.subs(r, rmax)
    inner_bc = sp.diff(Cl, r).subs(r, rmin) - sp.nsimplify(l + 1) / rmin * Cl.subs(r, rmin)

    sol = sp.solve([sp.nsimplify(outer_bc), sp.nsimplify(inner_bc)], [alpha, beta])
    return sp.simplify(Cl.subs(sol))


def magnetic_field_from_CA(C, A):
    """Given the poloidal potential C and toroidal potential A 
    sympy expressions in r,theta,phi), return (B_r, B_theta, B_phi) 
    as sympy expressions in (r, theta, phi).
    """
    curlCr, curlCt, curlCp = curl(C, 0, 0)   # first curl of C*rhat
    PBr, PBt, PBp = curl(curlCr, curlCt, curlCp)  # poloidal part of B
    TBr, TBt, TBp = curl(A, 0, 0)  # toroidal part of B (single curl of A*rhat)
    Br, Bt, Bp = PBr + TBr, PBt + TBt, PBp + TBp
    return sp.simplify(Br), sp.simplify(Bt), sp.simplify(Bp)


def magnetic_field_quantities(C, A):
    """b_r, b_theta, b_phi and their full-field r/theta/phi
    derivatives, as sympy expressions in (r, theta, phi). Derivatives
    are obtained by differentiating the b_r/b_theta/b_phi expressions from
    `magnetic_field_from_CA` directly with sympy, matching the same
    r/theta/phi (and the 1/r, 1/(r*sin(theta)) scaled "tr"/"prs") derivative
    conventions used for velocity.
    """
    b = dict(zip(('r', 'theta', 'phi'), magnetic_field_from_CA(C, A)))

    out = {'b_r': b['r'], 'b_theta': b['theta'], 'b_phi': b['phi']}
    for comp, f in b.items():
        out[f'db_{comp}_dr'] = sp.diff(f, r)
        out[f'db_{comp}_dt'] = sp.diff(f, theta)
        out[f'db_{comp}_dp'] = sp.diff(f, phi)
        out[f'db_{comp}_dtr'] = sp.diff(f, theta) / r
        out[f'db_{comp}_dprs'] = sp.diff(f, phi) / (r * sp.sin(theta))
        out[f'db_{comp}_d2r'] = sp.diff(f, r, 2)
        out[f'db_{comp}_d2t'] = sp.diff(f, theta, 2)
        out[f'db_{comp}_d2p'] = sp.diff(f, phi, 2)
        out[f'db_{comp}_d2rt'] = sp.diff(f, r, theta)
        out[f'db_{comp}_d2rp'] = sp.diff(f, r, phi)
        out[f'db_{comp}_d2tp'] = sp.diff(f, theta, phi)
    return out
