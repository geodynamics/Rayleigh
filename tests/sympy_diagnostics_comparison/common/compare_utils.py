"""Shared full3d + Point_Probes comparison, used by compare.py.
"""

import numpy as np

from tolerances import tolerance
from rayleigh_diagnostics import Spherical_3D_multi, Point_Probes


def compare_full3d_and_probes(quantity_codes, numeric, full3d_path='Spherical_3D/', probe_path='Point_Probes/'):
    """quantity_codes: {code: name}. numeric: {name: callable(r, theta, phi)}.
    Prints a full3d and point-probe report and returns True iff everything
    is within tolerance (see common/tolerances.py).
    """
    error = False

    # -------------------------------------------------------------
    # Full 3D grid comparison
    # -------------------------------------------------------------
    F = Spherical_3D_multi('00000001', path=full3d_path)
    radius = F.rs       # (nr,)
    thetas = F.thetas   # (ntheta,)
    phis = np.linspace(0, 2 * np.pi, F.nphi, endpoint=False)

    PHI, THETA, RADIUS = np.meshgrid(phis, thetas, radius, indexing='ij')  # each (nphi, ntheta, nr)

    print("--- full3d ---")
    for code, name in sorted(quantity_codes.items()):
        rayleigh_val = F.vals[f'{code:05d}']  # (nphi, ntheta, nr)
        analytic = np.broadcast_to(numeric[name](RADIUS, THETA, PHI), rayleigh_val.shape)
        tol = tolerance(name, THETA)

        diff = np.abs(rayleigh_val - analytic)
        maxratio = (diff / tol).max()
        maxdiff = diff.max()
        maxval = np.abs(analytic).max()
        print(f"{name} ({code}): max|Rayleigh - analytic| = {maxdiff:.3e}  (max|analytic| = {maxval:.3e},"
              f" max diff/tol = {maxratio:.3f})")
        if np.any(diff > tol):
            print(f"ERROR: {name} mismatch exceeds tolerance!")
            error = True

    # -------------------------------------------------------------
    # Point probe comparison: same quantities, at fixed (r, phi) and three
    # colatitudes -- one at the equator, two near the poles.
    # -------------------------------------------------------------
    P = Point_Probes('00000001', path=probe_path)
    probe_radius = P.radius[0]
    probe_phi = P.phi[0]
    probe_thetas = np.arccos(P.costheta)

    print("\n--- point probes ---")
    for ti, th in enumerate(probe_thetas):
        location = "equator" if abs(np.sin(th) - 1.0) < 0.1 else "near-pole"
        print(f"probe theta={th:.4f} (sin={np.sin(th):.4f}, {location}):")
        for code, name in sorted(quantity_codes.items()):
            rayleigh_val = P.vals[0, ti, 0, P.lut[code], 0]
            analytic = numeric[name](probe_radius, th, probe_phi)
            tol = tolerance(name, th)
            diff = abs(rayleigh_val - analytic)
            status = "OK" if diff <= tol else "FAIL"
            if status == "FAIL":
                error = True
            print(f"  {name} ({code}): |Rayleigh - analytic| = {diff:.3e}  (tol = {tol:.3e})  {status}")

    return not error
