"""Closed-form spherical-shell Stokes solutions, transcribed from the
`assess` package (https://github.com/stephankramer/assess),
specifically `assess/spherical.py` and `assess/smooth.py`, expressed as
sympy expressions in (r, theta, phi).

assess defines its solutions through a poloidal potential
    P(r, theta, phi) = Pl(r) * Y_lm(theta, phi)
with velocity
    u = curl(r_vec x grad(P))
(r_vec = r * rhat, the position vector).

Rayleigh's velocity poloidal potential W uses a different convention
(rho*v = curl(curl(W * rhat))) but the two are related for a single (l, m) 
mode by:

    W_l(r) = -rho(r) * r * Pl(r)

This is derived in common/derive_poloidal.py.

Only the "smooth" (r^k forcing) family, free-slip and zero-slip (no-slip),
is transcribed here.
"""

import sympy as sp

from sympy_utils import r, theta, phi


def Y(l, m):
    """Real-valued spherical harmonic Y_lm(theta, phi) as a sympy expression,
    matching the normalization sqrt((2l+1)/(4 pi) (l-m)!/(l+m)!) P_l^m(cos
    theta) cos(m phi) used by both assess and Rayleigh's generic spectral
    input (see `compute_plms` in pre_processing/rayleigh_spectral_input.py).

    Built from sympy's real spherical harmonic Znm, which carries an extra
    sqrt(2) normalization for m>0 relative to our convention (Znm is
    (Y_lm + conj(Y_lm))/sqrt(2), whereas our Y is plain Re(Y_lm)).
    """
    z = sp.Znm(l, m, theta, phi).expand(func=True)
    if m != 0:
        z = z / sp.sqrt(2)
    return sp.simplify(z)


def dYdtheta(l, m):
    return sp.diff(Y(l, m), theta)


def dYdphi(l, m):
    return sp.diff(Y(l, m), phi)


def coefficients_sphere_smooth_fs(Rp, Rm, k, l, g, nu):
    """A, B, C, D, E coefficients of the smooth (r^k forcing), free-slip
    spherical poloidal solution. Transcribed from
    assess/smooth.py::coefficients_sphere_smooth_fs."""
    alpha = Rm / Rp
    A = 0.5 * Rp**(-l + 3) * (alpha**(k + 3) - alpha**(-l + 1)) * g / ((alpha**l - alpha**(-l + 1)) * (k + l + 2) * (k - l + 3) * (2 * l + 1) * nu)
    B = 0.5 * Rp**(l + 4) * (alpha**(k + 4) - alpha**(l + 3)) * g / ((alpha**(l + 3) - 1 / alpha**l) * (k + l + 4) * (k - l + 1) * (2 * l + 1) * nu)
    C = -0.5 * Rp**(-l + 1) * (alpha**(k + 4) - 1 / alpha**l) * g / ((alpha**(l + 3) - 1 / alpha**l) * (k + l + 4) * (k - l + 1) * (2 * l + 1) * nu)
    D = -0.5 * Rp**(l + 2) * (alpha**(k + 3) - alpha**l) * g / ((alpha**l - alpha**(-l + 1)) * (k + l + 2) * (k - l + 3) * (2 * l + 1) * nu)
    E = g / (Rp**k * (k + l + 4) * (k + l + 2) * (k - l + 3) * (k - l + 1) * nu)
    return A, B, C, D, E


def coefficients_sphere_smooth_ns(Rp, Rm, k, l, g, nu):
    """A, B, C, D, E coefficients of the smooth (r^k forcing), zero-slip
    (no-slip) spherical poloidal solution. Transcribed from
    assess/smooth.py::coefficients_sphere_smooth_ns."""
    alpha = Rm / Rp
    Gamma = ((alpha**(l + 1) + alpha**(l - 3)) * (2 * l + 1)**2 - 2 * alpha**(l - 1) * (2 * l + 3) * (2 * l - 1) - 4 * alpha**(3 * l) - 4 * alpha**(-l - 2)) * (k + l + 4) * (k + l + 2) * (k - l + 3) * (k - l + 1)
    A = ((alpha**(k + 2) + alpha**(l - 1)) * (k + l + 2) * (2 * l + 3) - (alpha**k + alpha**(l + 1)) * (k + l + 4) * (2 * l + 1) - 2 * (alpha**(k + 2 * l + 3) + alpha**(-l - 2)) * (k - l + 1)) * Rp**(-l + 3) * g / (Gamma * nu)
    B = ((alpha**(k + 2 * l + 1) + alpha**(l + 1)) * (k - l + 3) * (2 * l + 1) - (alpha**(k + 2 * l + 3) + alpha**(l - 1)) * (k - l + 1) * (2 * l - 1) - 2 * (alpha**(k + 2) + alpha**(3 * l)) * (k + l + 2)) * Rp**(l + 4) * g / (Gamma * nu)
    C = -((alpha**(k + 2) + alpha**(l - 3)) * (k + l + 2) * (2 * l + 1) - (alpha**k + alpha**(l - 1)) * (k + l + 4) * (2 * l - 1) - 2 * (alpha**(k + 2 * l + 1) + alpha**(-l - 2)) * (k - l + 3)) * Rp**(-l + 1) * g / (Gamma * nu)
    D = -((alpha**(k + 2 * l + 1) + alpha**(l - 1)) * (k - l + 3) * (2 * l + 3) - (alpha**(k + 2 * l + 3) + alpha**(l - 3)) * (k - l + 1) * (2 * l + 1) - 2 * (alpha**k + alpha**(3 * l)) * (k + l + 4)) * Rp**(l + 2) * g / (Gamma * nu)
    E = g / (Rp**k * (k + l + 4) * (k + l + 2) * (k - l + 3) * (k - l + 1) * nu)
    return A, B, C, D, E


class _SmoothPoloidalBase:
    """Shared radial part Pl(r) (and sympy-differentiated derivatives) and
    velocity reconstruction for assess's smooth spherical poloidal
    potential, for a single (l, m, k) mode. The functional form of Pl is the
    same for the free-slip and zero-slip families; only the A,B,C,D,E
    coefficients (computed by the subclass) differ.

    Also builds assess's pressure solution P(r,theta,phi) = (G r^l + H
    r^(-l-1) + K r^(k+1)) Y_lm, transcribed from
    assess/spherical.py::SphericalStokesSolutionSmooth. Unlike Pl (which
    needs rescaling into Rayleigh's W convention) assess's
    P can be used as Rayleigh's P directly.

    Raises ValueError for (k, l) combinations assess is not valid for.
    """

    _coefficients = None  # set by subclass: staticmethod(Rp, Rm, k, l, g, nu) -> (A,B,C,D,E)

    def __init__(self, l, m, k, Rp, Rm, nu=1.0, g=1.0):
        if (k + 1) * (k + 2) == l * (l + 1) or (k + 3) * (k + 4) == l * (l + 1):
            raise ValueError(f"Smooth solution not implemented for k={k}, l={l}")
        self.l, self.m, self.k = l, m, k
        self.Rp, self.Rm, self.nu, self.g = Rp, Rm, nu, g
        A, B, C, D, E = self._coefficients(Rp, Rm, k, l, g, nu)
        self.Pl = A * r**l + B * r**(-l - 1) + C * r**(l + 2) + D * r**(-l + 1) + E * r**(k + 3)
        self.dPldr = sp.diff(self.Pl, r)

        # Pressure solution coefficients (assess/spherical.py::SphericalStokesSolutionSmooth.__init__).
        self.G = -2 * nu * (l + 1) * (2 * l + 3) * C
        self.H = -2 * nu * l * (2 * l - 1) * D
        self.K = -g * (k + 2) / ((k + 1) * (k + 2) - l * (l + 1)) / Rp**k
        self.Pl_pressure = self.G * r**l + self.H * r**(-l - 1) + self.K * r**(k + 1)  # radial part only
        self.P = self.Pl_pressure * Y(l, m)  # full (r,theta,phi) field

    def u_r(self):
        """assess's own velocity - not usable because Rayleigh's poloidal convention differs from assess's."""
        return -self.l * (self.l + 1) * self.Pl * Y(self.l, self.m) / r

    def u_theta(self):
        return -(self.Pl / r + self.dPldr) * dYdtheta(self.l, self.m)

    def u_phi(self):
        return -(self.Pl / r + self.dPldr) / sp.sin(theta) * dYdphi(self.l, self.m)


class SmoothFreeSlipPoloidal(_SmoothPoloidalBase):
    """assess's smooth, free-slip (stress-free) spherical poloidal
    solution -- matches Rayleigh's `no_slip_boundaries = .false.` default."""
    _coefficients = staticmethod(coefficients_sphere_smooth_fs)


class SmoothZeroSlipPoloidal(_SmoothPoloidalBase):
    """assess's smooth, zero-slip (no-slip) spherical poloidal solution --
    matches Rayleigh's `no_slip_boundaries = .true.`."""
    _coefficients = staticmethod(coefficients_sphere_smooth_ns)
