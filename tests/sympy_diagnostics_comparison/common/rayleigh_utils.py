"""Shared helper for writing Rayleigh generic spectral input files
(w_init_file, c_init_file, ...) from an analytic radial function.
"""

import sys
import os
import numpy as np
import sympy as sp

from sympy_utils import r

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'pre_processing'))
from rayleigh_spectral_input import SpectralInput, compute_gamma, compute_tns  # noqa: E402


def write_mode_file(filename, radial_expr, l, m, rmin, rmax, n_r):
    """Chebyshev-decompose an analytic radial function (a sympy expression
    in r) on [rmin, rmax] into a sparse (n, l, m)
    generic spectral input file, writing only the single (l, m) column of
    Chebyshev modes actually needed.

    Reuses the exact Chebyshev collocation points/basis (`compute_gamma`,
    `compute_tns`) and forward-transform normalization from
    pre_processing/rayleigh_spectral_input.py.
    """
    return write_modes_file(filename, [(radial_expr, l, m)], rmin, rmax, n_r)


def write_modes_file(filename, modes, rmin, rmax, n_r):
    """As write_mode_file, for a superposition of modes given as a list of
    (radial_expr, l, m) tuples (distinct (l, m) pairs)."""
    n_max = n_r - 1
    gamma = compute_gamma(n_r)
    rx = np.cos(gamma)
    radius = (rx - rx[-1]) * (rmax - rmin) / (rx[0] - rx[-1]) + rmin
    tns = compute_tns(n_max, gamma)  # shape (n_max+1, n_r)

    si = SpectralInput()  # no lm_max/n_theta supplied -> sparse (n,l,m) list storage
    for radial_expr, l, m in modes:
        radial_num = sp.lambdify(r, radial_expr, 'numpy')
        vals = np.broadcast_to(np.asarray(radial_num(radius), dtype='float64'), radius.shape)
        # same forward Chebyshev transform as SpectralInput.transform_from_rtp_data,
        # specialised to a single (l, m) column (no angular quadrature).
        coeffs_n = (2.0 / n_r) * np.sum(vals[np.newaxis, :] * tns, axis=1)
        for n in range(n_max + 1):
            si.add_mode(coeffs_n[n], n=n, l=l, m=m)
    si.write(filename)
    return si
