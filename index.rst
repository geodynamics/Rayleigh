*Rayleigh*:  MHD in Spherical Geometry  
======================================

Rayleigh solves the magnetohydrodynamic (MHD) equations, in a rotating frame, within spherical shells,
using the anelastic or Boussinesq approximations.
Derivatives in Rayleigh are calculated using a spectral transform scheme.
Spherical harmonics are used as basis functions in the horizontal direction.
Chebyshev polynomials are employed in radius.
Time-stepping is accomplished using the semi-implicit Crank-Nicolson method
for the linear terms, and the Adams-Bashforth method for the nonlinear terms.
Both methods are second-order in time.

This documentation is structured into the following sections:

.. toctree::
   :maxdepth: 1

   doc/source/User_Guide/index.rst
   doc/source/citing_rayleigh
   doc/source/accessing_and_sharing_data
   doc/source/research_enabled_by_rayleigh
   doc/source/quick_reference
   doc/source/getting_help

