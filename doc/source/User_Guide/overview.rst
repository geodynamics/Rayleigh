.. raw:: latex

   \clearpage

.. _Overview:

Overview
==========

Rayleigh solves the MHD equations, in a rotating frame, within spherical shells,
using the anelastic or Boussinesq approximations.
Derivatives in Rayleigh are calculated using a spectral transform scheme.
Spherical harmonics are used as basis functions in the horizontal direction.
Chebyshev polynomials are employed in radius.
Time-stepping is accomplished using the semi-implicit Crank-Nicolson method
for the linear terms, and the Adams-Bashforth method for the nonlinear terms.
Both methods are second-order in time.

This document serves as a guide to installation and running Rayleigh.
Rayleigh's diagnostics package is discussed in the companion document :ref:`DValues2`.

Referencing
-----------
Rayleigh's implementation of the pseudo-spectral algorithm and its
parallel design would not have been possible without earlier work by
Gary Glatzmaier and Thomas Clune described in: :cite:`GLATZMAIER1984461`,
:cite:`Glatzmaier1999`

  Glatzmaier, G.A., 1984,Numerical simulations of stellar convective dynamos. I. the model and method,
  *J. Comp. Phys.*, 55(3), 461-484. ISSN 0021-9991, doi:10.1016/0021-9991(84)90033-0.

  Clune, T.C., Elliott, J.R., Miesch, M.S.,Toomre, J., and Glatzmaier, G.A., 1999,
  Computational aspects of a code to study rotating turbulent convection in
  spherical shells, Parallel Comp., 25, 361-380.


Rayleigh has been written by Nicholas Featherstone with contributions by others.
We ask that you cite the appropriate references if you publish results that were obtained in some
part using Rayleigh.  To cite other versions of the code, please see: https://geodynamics.org/cig/abc

Please cite the code as:

  Featherstone, N. (2018), Rayleigh Version 1.0.0, Computational Infrastructure for Geodynamics,
  DOI: 10.5281/zenodo.1236565

.. code-block::

  @Software{nicholas_featherstone_2018_1236565,
	author = "Featherstone, N.",
	title="Rayleigh 1.0.0",
	year="2018",
	organization="",
	optkeywords="Rayleigh",
	doi="http://doi.org/10.5281/zenodo.1236565",
	opturl="https://doi.org/10.5281/zenodo.1236565"}

Please also cite the following references:

  Featherstone, N.A.; Hindman, B.W. (2016), The spectral
  amplitude of stellar convection and its scaling in the
  high-rayleigh-number regime, *The Astrophysical Journal*, 818 (1) ,
  32, DOI: 10.3847/0004-637X/818/1/32

  Matsui, H. et al., 2016, Performance benchmarks for
  a next generation numerical dynamo model, *Geochem., Geophys., Geosys.*, 17,1586
  DOI: 10.1002/2015GC006159

.. code-block::

  @Article{,
  author = "Featherstone, N.A. and Hindman, B.W.",
  title="The Spectral Amplitude Of Stellar Convection And Its Scaling In The High-Rayleigh-Number Regime",
  year="2016",
  journal="The Astrophysical Journal",
  volume="818",
  number="1",
  pages="32",
  optkeywords="Rayleigh",
  issn="1538-4357",
  doi="http://doi.org/10.3847/0004-637X/818/1/32",
  opturl="http://stacks.iop.org/0004-637X/818/i=1/a=32?key=
  crossref.a90f82507dd0eeb7a6e7562d1e4b0210"}

  @Article{Matsui_etal_2016,
  author = "Matsui, H. and Heien, E. and Aubert, J. and Aurnou, J.M. and Avery, M. and Brown, B. and Buffett, B.A. and Busse, F. and Christensen, U.R. and Davies, C.J. and Featherstone, N. and Gastine, T. and Glatzmaier, G.A. and Gubbins, D. and Guermond, J.-L. and Hayashi, Y.-Y. and Hollerbach, R. and Hwang, L.J. and Jackson, A. and Jones, C.A. and Jiang, W. and Kellogg, L.H. and Kuang, W. and Landeau, M. and Marti, P.H. and Olson, P. and Ribeiro, A. and Sasaki, Y. and Schaeffer, N. and Simitev, R.D. and Sheyko, A. and Silva, L. and Stanley, S. and Takahashi, F. and Takehiro, S.-ichi and Wicht, J. and Willis, A.P.",
  title="Performance benchmarks for a next generation numerical dynamo model",
  year="2016",
  journal="Geochemistry, Geophysics, Geosystems",
  volume="17",
  number="5",
  pages="1586-1607",
  optkeywords="Calypso",
  issn="1525-2027",
  doi="http://doi.org/10.1002/2015GC006159",
  opturl="http://doi.wiley.com/10.1002/2015GC006159"
  }


Acknowledging
-------------
Rayleigh is written by Nicholas Featherstone with
National-Science-Foundation support through the Geodynamo Working Group
of the Computational Infrastructure for Geodynamics (CIG, geodynamics.org).

The CIG Geodynamo Working Group Members are:
Jonathon Aurnou, Benjamin Brown, Bruce Buffett, Nicholas Featherstone,
Gary Glatzmaier, Eric Heien, Moritz Heimpel, Lorraine Hwang, Louise Kellogg,
Hiroaki Matsui, Peter Olson, Krista Soderlund, Sabine Stanley, Rakesh Yadav.

Please acknowledge CIG as follows:

.. note::

  Rayleigh is hosted and receives support from the Computational
  Infrastructure for Geodynamics (CIG) which is supported by the
  National Science Foundation awards NSF-0949446 and NSF-1550901.
