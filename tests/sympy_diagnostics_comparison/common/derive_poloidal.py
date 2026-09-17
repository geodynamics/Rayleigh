"""Symbolic derivation of Rayleigh's poloidal velocity relation

    rho*v = curl(curl(W(r,theta,phi) * rhat))

used by velocity_field.py. Tests, for a single mode
W = Wl(r)*Y_lm(theta,phi), that the radial component reduces to the known
Rayleigh relation v_r = l(l+1)/(rho*r^2) * Wl(r) * Y_lm (matching
`SBUFFA(IDX2,vr) = l_l_plus1(...)*SBUFFA(IDX2,vr)*Over_RhoRSQ(r)` in
src/Physics/Sphere_Hybrid_Space.F90), and reads off the corresponding
v_theta, v_phi formulas. Not part of the automated comparison.
"""

import sympy as sp

from sympy_utils import r, theta, phi, curl

Wl = sp.Function('Wl')(r)
Y = sp.Function('Y')(theta, phi)
W = Wl * Y

# A = W * rhat (unit vector): physical components (W, 0, 0)
curlWr, curlWt, curlWp = curl(W, 0, 0)
Pr, Pt, Pp = curl(curlWr, curlWt, curlWp)
Pr, Pt, Pp = sp.simplify(Pr), sp.simplify(Pt), sp.simplify(Pp)

print("P = curl(curl(W * rhat)) = rho * v, for W = Wl(r)*Y(theta,phi):")
print(" Pr =", Pr)
print(" Pt =", Pt)
print(" Pp =", Pp)

horiz_lap_Y = (1 / sp.sin(theta)) * sp.diff(sp.sin(theta) * sp.diff(Y, theta), theta) + (1 / sp.sin(theta)**2) * sp.diff(Y, phi, 2)
print()
print("Pr matches -1/r^2 * (horizontal Laplacian of Y) * Wl(r):",
      sp.simplify(Pr - (-1 / r**2) * horiz_lap_Y * Wl) == 0)
print("(horizontal Laplacian of Y_lm = -l(l+1) Y_lm, giving Pr = l(l+1)/r^2 * Wl(r) * Y_lm,")
print(" i.e. v_r = l(l+1)/(rho*r^2) * Wl(r) * Y_lm -- matches Rayleigh's Over_RhoRSQ(r) scaling.)")
print()
print("Pt = (1/r) * dWl/dr * dY/dtheta  ->  v_theta = 1/(rho*r) * dWl/dr * dY_lm/dtheta")
print("Pp = 1/(r*sin(theta)) * dWl/dr * dY/dphi  ->  v_phi = 1/(rho*r*sin(theta)) * dWl/dr * dY_lm/dphi")
