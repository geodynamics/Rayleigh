![rotator_and_cutaway.jpg](https://bitbucket.org/repo/Rp975y/images/1513682443-rotator_and_cutaway.jpg)

# What is Rayleigh? #

Rayleigh is a 3-D convection code designed for the study of planetary and stellar dynamos.  In particular, it evolves the incompressible and anelastic MHD equations in spherical geometry using a pseudo-spectral approach.  Rayleigh employs spherical harmonics in the horizontal direction and Chebyshev polynomials in the radial direction.  The code has undergone extensive accuracy testing using the [Christensen et al. (2001)](http://adsabs.harvard.edu/abs/2001PEPI..128...25C) Boussinesq benchmarks and the [Jones et al. (2011)](http://adsabs.harvard.edu/abs/2011Icar..216..120J) anelastic benchmarks.   Rayleigh has been developed with NSF support through the Computational Infrastructure for Geodynamics ([CIG](https://geodynamics.org/cig/news/newsletters/may-2016/)).

# What Makes Rayleigh Unique? #
The pseudo-spectral nature of Rayleigh means that its parallelization necessarily relies heavily on **global communication** patterns.  That said, Rayleigh's parallelization is based around a 2-D domain decomposition and large-message-size all-to-alls.  These features allow the code to overcome many of the obstacles that traditionally limit the scalability of spectral codes.   The end result is a pseudo-spectral code optimized for petascale machines.  Rayleigh's pure-MPI mode has demonstrated highly efficient strong scaling on  131,000 cores of the Mira Blue Gene/Q supercomputer for problems with approximately 2048^3 grid points (2048 spherical harmonics).  Performance numbers from Mira are shown below.  A summary of Rayleigh's performance and how it compares against other popular dynamo codes (albeit at at smaller process counts) may be found in the recent performance benchmark results of [Matsui et al. (2016)](http://onlinelibrary.wiley.com/doi/10.1002/2015GC006159/full).
 
![mira_performance.jpg](https://bitbucket.org/repo/Rp975y/images/3897197863-mira_performance.jpg)


# When Will Rayleigh Be Released?#
Rayleigh is scheduled for release on June 19, 2016 at the [CIG All Hands Meeting](https://geodynamics.org/cig/events/calendar/2016-cig-all-hands-meeting/)  In the meantime, please pardon the mess.  The code and (and presently incomplete) documentation are in a state of flux.  See the Wiki and the documentation directory (in the source) for tips.  Some videos corresponding to the images above can be found [here](http://www.youtube.com/user/feathern24).