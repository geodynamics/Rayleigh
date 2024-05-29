.. _citing_rayleigh:

Citing Rayleigh
===============

We ask that you cite the appropriate references if you publish results that were obtained in some
part using Rayleigh.  Receiving citations for Rayleigh is important to demonstrate the relevance of our work to our funding agencies and is a matter of fairness to all the developers that have donated their effort and time to make Rayleigh what it is today.

Please cite the code as:
  Featherstone, Nicholas A., Edelmann, Philipp V. F., Gassmoeller, Rene, Matilsky, Loren I., & Wilson, Cian R. (2024). Rayleigh Version 1.2.0 (v1.2.0). Zenodo. https://doi.org/10.5281/zenodo.6522806

To cite other versions of the code, please see: https://geodynamics.org/resources/rayleigh/howtocite

.. code-block::

  @Software{featherstone_et_al_2024,
	author = "{Featherstone}, N.~A. and {Edelmann}, P.~V.~F. and {Gassmoeller}, R. and {Matilsky}, L.~I. and {Wilson}, C.~R.",
	title="Rayleigh 1.2.0",
	year="2024",
	organization="",
	optkeywords="Rayleigh",
	doi="http://doi.org/10.5281/zenodo.6522806",
	opturl="https://doi.org/10.5281/zenodo.6522806"}

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

Rayleigh's development is supported by the
National Science Foundation through the Dynamo Working Group of the
Computational Infrastructure for Geodynamics (CIG,
https://geodynamics.org/groups/dynamo).

Please acknowledge CIG support in your work as follows:

  .. note::

    Rayleigh is hosted and receives support from the Computational
    Infrastructure for Geodynamics (CIG) which is supported by the
    National Science Foundation awards NSF-0949446, NSF-1550901 and NSF-2149126.

Publishing
----------

Open research statements are now a common requirement when publishing research. 
These support reuse, validation, and citation and often take the form of *Data availability, Data access, 
Code availability, Open Research*, and *Software availability* statements.
We recommend depositing input files that allow your published research to be reproduced and 
output model data in support of your research outcomes and figures. 
In addition, consider depositing model files that may be reused by others.

Remember to cite software and data in your text as well as in your Data Availability or similar statement.

Files should be deposited in an approved repository.

Additional information on `Publishing <https://geodynamics.org/software/software-bp/software-publishing` is available on the CIG website.

Data
~~~~

**Input parameters**

* Main_input, *_input 
* Data files for custom: 
  * profiles, 
  * boundary conditions, 
  * generic initial conditions, 
  * reference states (coefficients)
* Basic simulation information e.g. grid, job

**Model output**

Data products/checkpoints for the cases used in your publication.

Repository
~~~~~~~~~~

The Rayleigh Simulation Library (RSL), a repository for accessing published Rayleigh datasets has been 
established using the Open Science Framework (OSF) at the University of Colorado Boulder. 
For more information on this repository and preparing your datasets see the RSL home page:
https://osf.io/j275z/

Template
~~~~~~~~

  We use Rayleigh version number (Featherstone et al., XXXX; Featherstone and Hindman, 2016, Matsui et al., 2016) 
  which is available for download through its software landing page https://geodynamics.org/resources/rayleigh 
  or from Zenodo<insert PID>. Model data necessary to reproduce these results including <insert a description> 
  can be downloaded from Zenodo <insert PID> (Authors, YYYY).

  Featherstone, N.A.; Hindman, B.W. (2016), The spectral amplitude of stellar convection and its scaling in the high-rayleigh-number regime, The Astrophysical Journal, 818 (1) , 32, DOI: 10.3847/0004-637X/818/1/32

  Matsui, H. et al., 2016, Performance benchmarks for a next generation numerical dynamo model, Geochem., Geophys., Geosys., 17,1586 DOI: 10.1002/2015GC006159

  Authors (ZZZZ), ….


Where XXXX refers to the appropriate year of the software version cited and Authors (ZZZZ) is the citation to the data.

`IOP <https://publishingsupport.iopscience.iop.org/iop-publishing-standard-data-policy/>` (The Astrophysical Journal) recommends the following form:
  The data that support the findings of this study are openly available at the following 
  URL/DOI: [insert web link or DOI to the data].

See above or https://geodynamics.org/resources/rayleigh/howtocite for the citation to the version used.

Published examples
~~~~~~~~~~~~~~~~~~

https://doi.org/10.5281/zenodo.7117668 


Acknowledging
-------------

.. include:: ../../../AUTHORS

Rayleigh's implementation of the pseudo-spectral algorithm and its
parallel design would not have been possible without earlier work by
Gary Glatzmaier and Thomas Clune described in: :cite:`GLATZMAIER1984461`,
:cite:`Glatzmaier1999`

  Glatzmaier, G.A., 1984, Numerical simulations of stellar convective dynamos. I. the model and method,
  *J. Comp. Phys.*, 55(3), 461-484. ISSN 0021-9991, doi:10.1016/0021-9991(84)90033-0.

  Clune, T.C., Elliott, J.R., Miesch, M.S.,Toomre, J., and Glatzmaier, G.A., 1999,
  Computational aspects of a code to study rotating turbulent convection in
  spherical shells, Parallel Comp., 25, 361-380.
