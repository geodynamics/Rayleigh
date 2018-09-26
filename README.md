
# Rayleigh Version 0.9.1 #
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1236565.svg)](https://doi.org/10.5281/zenodo.1236565)



Rayleigh is a 3-D convection code designed for the study of dynamo behavior in spherical geometry.  It evolves the incompressible and anelastic MHD equations in spherical geometry using a pseudo-spectral approach.  Rayleigh employs spherical harmonics in the horizontal direction and Chebyshev polynomials in the radial direction.  The code has undergone extensive accuracy testing using the [Christensen et al. (2001)](http://adsabs.harvard.edu/abs/2001PEPI..128...25C) Boussinesq benchmarks and the [Jones et al. (2011)](http://adsabs.harvard.edu/abs/2011Icar..216..120J) anelastic benchmarks.   Rayleigh has been developed with NSF support through the Computational Infrastructure for Geodynamics ([CIG](https://geodynamics.org/cig/news/newsletters/may-2016/)).


Contributing to Rayleigh
------------------------

Rayleigh is a community project that lives by the participation of its members
-- i.e., including you! It is our goal to build an inclusive and participatory
community so we are happy that you are interested in participating! We have
collected a set of guidelines and advice on how to get involved in the
community and keep them in the [CONTRIBUTING.md](CONTRIBUTING.md) file in
Rayleigh's repository.


Parallelization
-----------------------------
The pseudo-spectral nature of Rayleigh means that its parallelization necessarily relies heavily on *global communication* patterns.  That said, Rayleigh's parallelization is based around a 2-D domain decomposition and large-message-size all-to-alls.  These features allow the code to overcome many of the obstacles that traditionally limit the scalability of spectral methods.   The end result is a pseudo-spectral code optimized for petascale machines.  Rayleigh's pure-MPI mode has demonstrated highly efficient strong scaling on  131,000 cores of the Mira Blue Gene/Q supercomputer for problems with approximately 2048^3 grid points (2048 spherical harmonics).  Performance numbers from Mira are shown below.  A summary of Rayleigh's performance and how it compares against other popular dynamo codes (albeit at at smaller process counts) may be found in the recent performance benchmark results of [Matsui et al. (2016)](http://onlinelibrary.wiley.com/doi/10.1002/2015GC006159/full).
 
 Getting Started
 ----------------
The following documents form the initial Rayleigh documentation.  Note that the html version of the Diagnostics_Plotting document, which is automatically generated from a Jupyter notebook, is more cleanly formatted than the pdf version.

| Document | Description |
|----------|-------------|
| INSTALL | in-depth installation instructions |
| doc/Rayleigh_User_Guide.pdf |  Overview of compiling and running Rayleigh |
| doc/Diagnostic_Plotting.{html,pdf} | Guide to generating and plotting Rayleigh output (text version) |
| etc/analysis/Diagnostic_Plotting.ipynb | Same as above (interactive Jupyter notebook) |
| doc/Rayleigh_Output_Variables.pdf |  Tables of Rayleigh output menu codes |

Authors
--------
Rayleigh was written by Nicholas Featherstone with NSF support through CIG.  Please see the *ACKNOWLEDGEMENT* file for citation information.

License
-------
Raleigh v0.9.1 is released under the [GPL v3 license.](https://www.gnu.org/licenses/gpl-3.0.en.html)
