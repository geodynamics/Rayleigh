.. raw:: latex

   \clearpage

.. _sec:installation:

Overview
==========

Rayleigh solves the MHD equations, in a rotating frame, within spherical shells,
using the anelastic or Boussinesq approximations.
Derivatives in Rayleigh are calculated using a spectral transform scheme.
Spherical harmonics are used as basis functions in the horizontal direction.
Chebyshev polynomials are employed in radius.
Time-stepping is accomplished used the semi-implicit Crank-Nicolson method
for the linear terms, and the Adams-Bashforth method for the nonlinear terms.
Both methods are second-order in time.

This document serves as a guide to installation and running Rayleigh.
Rayleigh's diagnostics package is discussed in the companion document
Rayleigh/doc/Diagnostics\_Plotting.{html,pdf}

Rayleigh is written by Nicholas Featherstone, with
National-Science-Foundation support through the Geodynamo Working Group
of the Computational Infrastructure for Geodynamics (CIG, PI: Louise Kellogg).

he CIG Geodynamo Working Group Members are:
Jonathon Aurnou, Benjamin Brown, Bruce Buffett, Nicholas Featherstone,
Gary Glatzmaier, Eric Heien, Moritz Heimpel, Lorraine Hwang, Louise Kellogg,
Hiroaki Matsui, Peter Olson, Krista Soderlund, Sabine Stanley, Rakesh Yadav.

Rayleigh's implementation of the pseudo-spectral algorithm and its
parallel design would not have been possible without earlier work by
Gary Glatzmaier and Thomas Clune, described in:

#. Glatzmaier, G.A., 1984, *J. Comp. Phys.*, 55, 461

#. Clune, T.C., Elliott, J.R., Miesch, M.S., and Toomre, J.,1999, *Parallel Comp.*, **25**, 361
