.. raw:: latex

   \clearpage

.. _model_setup:

Setting Up a Model
==================

This section details the basics of running a custom model
in Rayleigh. For a commplete list of all Rayleigh input 
parameters, see :ref:`namelists`.

.. _setup_prep:
Preparation
-----------

Each simulation run using Rayleigh should have its own directory. The
code is run from within that directory, and any output is stored in
various subdirectories created by Rayleigh at run time. Wherever you
create your simulation directory, ensure that you have sufficient space
to store the output.

**Do not run Rayleigh from within the source code directory.
Do not cross the beams: no running two models from within the same
directory.**

After you create your run directory, you will want to copy (cp) or soft
link (ln -s ) the executable from Rayleigh/bin to your run directory.
Soft-linking is recommended; if you recompile the code, the executable
remains up-to-date. If running on an IBM machine, copy the script named
Rayleigh/etc/make_dirs to your run directory and execute the script.
This will create the directory structure expected by Rayleigh for its
outputs. This step is unnecessary when compiling with the Intel, GNU,
AOCC, or Cray compilers.

Next, you must create a main_input file. This file contains the
information that describes how your simulation is run. Rayleigh always
looks for a file named main_input in the directory that it is launched
from. Copy one of the sample input files from the
Rayleigh/input_examples/ into your run directory, and rename it to
main_input. The file named *benchmark_diagnostics_input* can be used to
generate output for the diagnostics plotting tutorial (see
§\ :ref:`diagnostics`).

Finally, Rayleigh has some OpenMP-related logic that is still in
development. We do not support Rayleigh’s OpenMP mode at this time, but
on some systems, it can be important to explicitly disable OpenMP in
order to avoid tripping any OpenMP flags used by external libraries,
such as Intel’s MKL. Please be sure and run the following command before
executing Rayleigh. This command should be precede *each* call to
Rayleigh.

::

   export OMP_NUM_THREADS=1 (bash)
   setenv OMP_NUM_THREADS 1 (c-shell)

.. _setup_grid:

Grid Setup
----------

By default, Rayleigh employs a single-domain Chebyshev decomposition in radius and a spherical-harmonic decomposition in the :math:`\theta-\phi` directions. Additionally, multiple Chebyshev domains or a finite-difference scheme may be alternatively employed in radius.  We focus on the default mode first.

Standard Grid Specification
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The number of radial grid points is denoted by
:math:`N_r`, and the number of :math:`\theta` grid points by
:math:`N_\theta`. The number of grid points in the :math:`\phi`
direction is always :math:`N_\phi=2\times N_\theta`. :math:`N_r` and
:math:`N_\theta` may each be defined in the problemsize_namelist of
main_input:

::

   &problemsize_namelist
    n_r = 48
    n_theta = 96
   /

:math:`N_r` and :math:`N_\theta` may also be specified at the command
line (overriding the values in main_input) via:

::

   mpiexec -np 8 ./rayleigh.opt -nr 48 -ntheta 96

If desired, the number of spherical harmonic degrees :math:`N_\ell` or the maximal spherical harmonic degree
:math:`\ell_\mathrm{max}\equiv N_\ell-1` may be specified in lieu of
:math:`N_\theta`.  The example above may equivalently be written as

::

   &problemsize_namelist
    n_r = 48
    l_max = 63
   /

or

::

   &problemsize_namelist
    n_r = 48
    n_l = 64
   /

The radial domain bounds are determined by the namelist variables
``rmin`` (the lower radial boundary) and ``rmax`` (the upper
radial boundary):

::

   &problemsize_namelist
    rmin = 1.0
    rmax = 2.0
   /

Alternatively, the user may specify the shell depth (``rmax-rmin``)
and aspect ratio (``rmin/rmax``) in lieu of ``rmin`` and
``rmax``. The preceding example may then be written as:

::

   &problemsize_namelist
    aspect_ratio = 0.5
    shell_depth = 1.0
   /

Note that the interpretation of ``rmin`` and ``rmax`` depends on
whether your simulation is dimensional or nondimensional. We discuss
these alternative formulations in §\ :ref:`physics_math`

Using Multiple Chebyshev Domains in Radius
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is possible to run Rayleigh with multiple, stacked domains in the
radial direction. Each of these is discretized using their own set of
Chebyshev polynomials. The boundaries and number of polynomials can be
set for each domain indiviadually, which makes it possible to control
the radial resolution at different radii.

To use this feature the problem size has to be specified using
``domain_bounds`` and ``ncheby`` instead of ``rmin``, ``rmax``, and
``n_r``. ``ncheby`` takes a comma-separated list of the number of radial
points to use in each domain. ``domain_bounds`` takes a comma-separated
list of the radii of the domain boundaries, starting with the smallest
radius. It has one element more than the number of domains. This is an
example of two radial domains, one covering the radii 1 to 2 with 16
radial points, the other the radii 2 to 4 with 64 radial points.

::

   &problemsize_namelist
    domain_bounds = 1.0, 2.0, 4.0
    ncheby = 16, 64
   /

Radial values in the diagnostic output will be repeated at the inner
domain boundaries. Most quantities are forced to be continuous at these
points.


Employing a Finite-Difference Approach in Radius
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Rayleigh's default behavior is to employ a Chebyshev collocation scheme in radius.   If desired, a finite-difference method can be applied instead.  This mode is activated by setting the value of ``chebyshev`` to .false. in the ``numerical_controls_namelist``.  At present, Rayleigh's finite-difference scheme employs a five-point stencil with 4th-order accuracy in the interior points.  Boundary derivatives are taken with second-order accuracy.   By default, a uniform radial grid is assumed.  Consider the following example:

::

   &problemsize_namelist
    rmin = 1.0
    rmax = 2.0
    n_r = 4
   /
   &numerical_controls_namelist
    chebyshev=.false.
   /

This results in the uniform grid:

::

  radius = 1.000 , 1.333 , 1.667 , 2.000
      dr = 0.333 , 0.333 , 0.333

An example input file using a uniform radial grid and a finite-difference scheme is provided in ``input_examples\main_input_mhd_jones_FD``.  If desired, a nonuniform grid can also be generated.  There are two ways to do this: via *main_input* and via a grid-description file. 


Using Main_Input to Specify a Nonuniform Grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. _nonuniform grids:

The first method of specifying a nonuniform grid is to join together a series of uniformly-gridded subdomains with different grid spacings. This is accomplished using the ``dr_weight`` and ``nr_count`` parameters.  ``nr_count`` indicates the number of gridpoints within each subregion, and ``dr_weight`` indicates the relative size of the grid spacing within each region. Consider the following example:

::

   &problemsize_namelist
    rmin = 1.0
    rmax = 2.0
    n_r = 4
    dr_weights = 0.1,0.3,0.2
    nr_count = 2,4,2
   /

This example defines a nonuniform grid ranging from 1.0 to 2.0 with 8 gridpoints (Rayleigh will reset the value of ``n_r`` to be the total of nr_count).  The grid spacing within the first 2-point region will be 1/3 of that in the second, 4-point region.  Similarly, the grid-spacing in the third, 2-point region will be 2/3 that of the second region and twice that of the first region.  The resulting radial grid and spacing is:

::

  radius = 1.000 , 1.059 , 1.235 , 1.412 , 1.588 , 1.765 , 1.882 , 2.   
      dr = 0.059 , 0.176 , 0.176 , 0.176 , 0.176 , 0.118 , 0.118

Note that for n_r points, there are n_r-1 spaces between gridpoints.  Rayleigh's convention is to apply dr_weights(1) nr_count(1)-1 times. As a result, specifying a symmetric ``nr_count`` will lead to assymetry in the grid spacing.  We can adjust this by adding one to nr_count(1) and subtracting one from nr_count(2) so that we have:
::

   &problemsize_namelist
    rmin = 1.0
    rmax = 2.0
    n_r = 4
    dr_weights = 0.1,0.3,0.2
    nr_count = 3,3,2
   /

This results in the symmetric grid:

::

  radius = 1.000 , 1.067 , 1.133 , 1.333 , 1.533 , 1.733 , 1.867 , 2.000
      dr = 0.067 , 0.067 , 0.200 , 0.200 , 0.200 , 0.133 , 0.133   

Be sure to leave the ``nr_count`` and ``dr_weights`` parameters unset in ``main_input`` if you wish to use a uniform grid in radius. 

.. _grid file section:

Using a Grid-Description File to Specify a Nonuniform Grid
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**NOTE:** The functionality described below is currently incompatible with Rayleigh's ensemble mode.  

An arbitrary radial grid may also be generated using Python and then stored to a file that is read when Rayleigh initializes.  To do so, import the *reference_tools* module and define a custom grid as illustrated by the code snippet below.
::
   import numpy  # Import necessary modules
   import reference_tools as rt 

   ri = 0.5   # Inner radius
   ro = 1.5   # Outer radius
   nr = 128   # Number of radial points
   radius = numpy.linspace(ri,ro,nr)  # The radial grid

   my_grid = rt.radial_grid(radius)   # Instantiate the grid object

   my_grid.write('grid_layout_128.dat')   # Store contents to file

Note that we could have generated the grid in either ascending or descending order.  The *write* method accounts for the grid-ordering before storing its contents to the file.  Now that we have created a grid-description file ('grid_layout_128.dat' in this example), we indicate the relevant filename in ``main_input`` using the ``radial_grid_file`` parameter:

::

   &problemsize_namelist
    n_theta = 32
    radial_grid_file = 'grid_layout_128.dat'
   /
   &numerical_controls_namelist
    chebyshev=.false.
   /
There are two important points to be aware of:

1.  When ``radial_grid_file`` is specified, all information concerning the grid structure is derived from that file.  The values of ``rmin``, ``rmax``, ``N_R`` etc. are completely ignored.  For this reason, we strongly suggest indicating ``N_R`` in the grid file's name. 
2. In the event that ``radial_grid_file``, ``nr_count`` and ``dr_weights`` are simultaneously specified, the grid-description file takes precedence. 

There is one exception to point 1 above because there may be instances where the same grid structure is useful for problems with different values of *rmin* and *rmax*.  If desired, the grid stored in *radial_grid_file* can be rescaled to a new *rmin* and *rmax* by setting the *rescale_radial_grid* keyword to true:
::

   &problemsize_namelist
    n_theta = 32
    rmin = 1.0
    rmax = 2.0
    rescale_radial_grid = .true.
    radial_grid_file = 'grid_layout_128.dat'
   /
   &numerical_controls_namelist
    chebyshev=.false.
   /

The example above will generate a grid identical to that stored in the grid file, but rescaled to run from r=1 to r=2, rather than r=0.5 to r=1.5 as specified in the original Python code.



.. _numerical_controls:

Numerical Controls 
------------------

The Numerical_Controls namelist was added to facilitate fine-control over some aspects of Rayleigh's parallelization and is documented in :ref:`namelists`.  Two numerical_controls parameters worth mentioning that are particularly important for setting up a new model are the ``chebyshev`` and ``bandsolve`` keywords.   

The value of ``chebyshev`` is set to ``.true.`` by default.  When set to ``.false.``, a finite-difference scheme will be employed in radius rather than a Chebyshev collocation scheme.

The value of the ``bandsolve`` keyword is also set to ``.false.`` by default.  When set to ``.true.``, the otherwise dense matrices used in the implicit timestepping scheme will be recast in banded or block-banded form for the finite-difference and Chebyshev schemes respectively.  This can save memory and may offer performance gains.   Note that this mode has no effect for models run in Chebyshev mode with only 1 or 2 Chebyshev domains in radius.  A minimum of three Chebyshev domains is required before any memory savings is possible.


.. _physics_controls:

Physics Controls
----------------

Many physical effects can be turned on or off in Rayleigh. The details of 
what physics you want to include will depend on the type of model you want 
to run. Be careful, however, that if you are adapting an input file from 
the benchmark described in :ref:`benchmark` that you set 
:code:`benchmark_mode` to 0 or omit it entirely, as this will override 
other input flags in favor of running the specified benchmark.

A number of logical variables can be used to turn certain physics on
(value = .true.) or off ( value = .false.). These variables are
described in Table table_logicals_, with default
values indicated in brackets.

  .. _table_logicals:

.. centered:: **Table. Logicals.**

Variables in the Physical_Controls_Namelist
that may be specified to control run behavior (defaults indicated in
brackets)

   +-----------------------------------+-----------------------------------+
   | Variable [Default value]          | Description                       |
   +===================================+===================================+
   | magnetism [.false.]               | Turn magnetism on or off          |
   +-----------------------------------+-----------------------------------+
   | rotation [.false.]                | Turn rotation on or off (pressure |
   |                                   | is not scaled by E when off)      |
   +-----------------------------------+-----------------------------------+
   | lorentz_forces [.true.]           | Turn Lorentz forces on or off     |
   |                                   | (magnetism must be .true.)        |
   +-----------------------------------+-----------------------------------+
   | viscous_heating [.true.]          | Turn viscous heating on or off    |
   |                                   | (inactive in Boussinesq mode)     |
   +-----------------------------------+-----------------------------------+
   | ohmic_heating [.true.]            | Turn ohmic heating off or on      |
   |                                   | (inactive in Boussinesq mode)     |
   +-----------------------------------+-----------------------------------+

.. _initial_conditions:

Initial Conditions
------------------

A Rayleigh simulation may be initialized with a random thermal and/or
magnetic field, or it may be restarted from an existing checkpoint file
(see §\ :ref:`checkpointing` for a detailed
discussion of checkpointing). This behavior is controlled through the
**initial_conditions_namelist** and the **init_type** and
**magnetic_init_type** variables. The init_type variable controls the
behavior of the velocity and thermal fields at initialization time.
Available options are:

-  init_type=-1 ; read velocity and thermal fields from a checkpoint
   file

-  init_type=1 ; Christensen et al. (2001) case 0 benchmark init (
   {:math:`\ell=4,m=4`} temperature mode)

-  init_type=6 ; Jones et al. (2011) steady anelastic benchmark (
   {:math:`\ell=19,m=19`} entropy mode)

-  init_type=7 ; random temperature or entropy perturbation

-  init_type=8 ; user generated temperature or entropy perturbation
   (see Generic Initial Conditions below)

When initializing a random thermal field, all spherical harmonic modes
are independently initialized with a random amplitude whose maximum
possible value is determined by the namelist variable **temp_amp**. The
mathematical form of of this random initialization is given by

.. _eq_init:

.. math::

   T(r,\theta,\phi) = \sum_\ell \sum_m  c_\ell^m f(r)g(\ell)\mathrm{Y}_\ell^m(\theta,\phi),

where the :math:`c_\ell^m`\ ’s are (complex) random amplitudes,
distributed normally within the range [-temp_amp, temp_amp]. The radial
amplitude :math:`f(r)` is designed to taper off to zero at the
boundaries and is given by

.. math:: f(r) = \frac{1}{2}\left[1-\mathrm{cos}\left( 2\pi\frac{r-rmin}{rmax-rmin} \right)   \right].

The amplitude function :math:`g(\ell)` concentrates power in the
central band of spherical harmonic modes used in the simulation. It is
given by

.. math:: g(\ell) = \mathrm{exp}\left[  - 9\left( \frac{ 2\,\ell-\ell_\mathrm{max} }{ \ell_\mathrm{max} }  \right)^2 \right],

which is itself, admittedly, a bit random.

When initializing using a random thermal perturbation, it is important
to consider whether it makes sense to separately initialize the
spherically-symmetric component of the thermal field with a profile that
is in conductive balance. This is almost certainly the case when running
with fixed temperature conditions. The logical namelist variable
**conductive_profile** can be used for this purpose. It’s default value
is .false. (off), and its value is ignored completely when restarting
from a checkpoint. To initialize a simulation with a random temperature
field superimposed on a spherically-symmetric, conductive background
state, something similar to the following should appear in your
main_input file:

::

   &initial_conditions_namelist
   init_type=7
   temp_amp = 1.0d-4
   conductive_profile=.true.
   /
   
Alternatively, you may wish to specify an ell=0 initial thermal profile
that is neither random nor conductive.  To create your own profile, follow the example found in
Rayleigh/examples/custom_thermal_profile/custom_thermal_profile.ipynb.   Then, use the following combination
of input parameters in main_input:

::

   &initial_conditions_namelist
   init_type=7
   temp_amp = 1.0d-4
   custom_thermal_file = 'my_custom_profile.dat' 
   /

This will use the radial profile stored in my_custom_profile.dat for the ell=0 component of entropy/temperature
Random values will be used to initialize all other modes.

Magnetic-field initialization follows a similar pattern. Available
values for magnetic_input type are:

-  magnetic_init_type = -1 ; read magnetic field from a checkpoint file

-  magnetic_init_type = 1 ; Christensen et al. (2001) case 0 benchmark
   init

-  magnetic_init_type = 7 ; randomized vector potential

-  magnetic_init_type=8 ; user generated magnetic potential fields
   (see Generic Initial Conditions below)

For the randomized magnetic field, both the poloidal and toroidal
vector-potential functions are given a random power distribution
described by Equation eq_init_. Each mode’s random
amplitude is then determined by namelist variable **mag_amp**. This
variable should be interpreted as an approximate magnetic field strength
(it’s value is rescaled appropriately for the poloidal and toroidal
vector potentials, which are differentiated to yield the magnetic
field).

When initializing all fields from scratch, a main_input file should
contain something similar to:

::

   &initial_conditions_namelist
   init_type=7
   temp_amp = 1.0d-4
   conductive_profile=.true.  ! Not always necessary (problem dependent) ...
   magnetic_init_type=7
   mag_amp = 1.0d-1
   /


.. _sec:generic_ic:

Generic Initial Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~

The user can input any initial conditions from data files generated by
a python routine "rayleigh_spectral_input.py", which can be called as
a script or imported as a python class.

The available generic initial conditions options are

::

   &initial_conditions_namelist
   init_type=8
   T_init_file = '<filename>'  !! Temperature
   W_init_file = '<filename>'  !! Poloidal velocity potential
   Z_init_file = '<filename>'  !! Toroidal velocity potential
   P_init_file = '<filename>'  !! `Pressure` potential
   magneic_init_type=8
   C_init_file = '<filename>'  !! Poloidal magnetic potential
   A_init_file = '<filename>'  !! Toroidal magnetic potential

   /

where `T_init_file` is a user generated initial temperature field and
<filename> is the name of the file generated by the python script.  If
`T_init_file` is not specified the initial field will be zero by
default.  The same for the other fields.  Fields T, W, Z, and P are
only initialized from the file if `init_type=8`.  Fields C and A are
only initialized from file if `magnetic_init_type=8`.

To generate a generic initial condition input file, for example, if a user wanted to specify a single mode in that input file then they could just run the script:

::

   rayleigh_spectral_input.py -m 0 0 0 1.+0.j -o example


to specify (n,l,m) = (0,0,0) to have a coefficient 1.+0.j and output it to the file example.

This could also be done using the python as a module. In a python
shell this would look like:

::

   from rayleigh_spectral_input import *
   si = SpectralInput()
   si.add_mode(1., n=0, l=0, m=0)
   si.write('example')


For a more complicated example, e.g. the hydrodynamic benchmark from
Christensen et al. 2001, the user can specify functions of theta, phi
and radius that the python will convert to spectral:

::

   rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 96 -nr 64 -o example \
    -e 'import numpy as np; x = 2*radius - rmin - rmax;
    rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) -
    x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

in "script" mode.

Alternatively, in "module" mode in a python shell:

::

   from rayleigh_spectral_input import *
   si = SpectralInput(n_theta=96, n_r=64)
   rmin, rmax = radial_extents(aspect_ratio=0.35, shell_depth=1.0)
   def func(theta, phi, radius):
      x = 2*radius - rmin - rmax
      return rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) - x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)
   si.transform_from_rtp_function(func, aspect_ratio=0.35, shell_depth=1.0)
   si.write('example')


The above commands will generate a file called `example` which can be
called by

::

   &initial_conditions_namelist
   init_type=8
   T_init_file = 'example'

Note that these two examples will have produced different data formats - the first one sparse (listing only the mode specified) and the second one dense (listing all modes).

For more examples including magnetic potentials see `tests/generic_input`.

Custom Reference States
---------------------

If desired, the constant and nonconstant equation coefficients enumerated :ref:`here <equations_solved>` may be completely or partially specified by the user.  This allows the user to specify diffusivity profiles, background states, or nondimensionalizations that are not supplied by Rayleigh.  Two use cases are supported:  

1.  One of Rayleigh's predefined reference states may be used, but with some coefficients supplied instead through an auxilliary coefficients file.  For each predefined reference type, the coefficients file may be used to override volumetric heating (:math:`\,c_{10} + f_6\,`) and the transport coefficients :math:`\nu,\kappa,\eta` (:math:`\,c_5 + f_3,\,\, c_6 + f_5,\,\, c_7 + f_7\,` ).  For Boussinesq runs, the buoyancy term may also be modified (:math:`\,c_2 + f_2\,`). 

2.  Nonconstant coefficents may be completely specified through the coefficients file.  In this mode, activated by setting reference_type=4, the user must fully specify nonconstant coefficients :math:`f_1-f_7`.  


In either case, constant coefficients may be defined within the coefficients file or, read in from main_input, or some combination of the two.  Moreover, the radial variation of transport coefficients, as specified by nu_type, kappa_type, and eta_type flags is respected.  We elaborate on this behavior below.

Creating a Coefficients File
~~~~~~~~~~~~~~~~~~~~~~~~~~

The first step in modifying Rayleigh's equation coefficients is to generate an equation coefficients file.
This file will be used alongside options defined in main_input to determine which combination of coefficients are overridden.
In order create your coefficients file, you will need to create an instance of the equation_coefficients class, provided
in post_processing/reference_tools.py.  Constant and nonconstant coefficients may then be set through  set_constant and set_function methods respectively.

The equation_coefficient class is instantiated by passing a radial grid to its init method.  This grid can be cast in ascending or descending order, but it should generally possess a much finer mesh than what you plan to use in Rayleigh.
Nonconstant coefficients specified in the coefficients file will be interpolated onto the Rayleigh grid at input time.   

The file structure created through the class's write method contains a record of those functions and contants that have been set. Rayleigh uses this information at runtime along with main_input to to perform consistency checks and to determine the values ultimately assigned to each constant coefficient. 

The sample code below defines a file with sufficient information to alter the viscous, heating, and buoyancy functions of a Rayleigh-provided reference state.  This information would be insufficient for use with reference_type=4, but several example notebooks handling that scenario are provided below. 

.. code-block::

    import numpy
    from reference_tools import equation_coefficients

    #Define a name for your equation coefficients file
    ofile = 'my_coeffs.dat'

    # Define the radial grid.  We suggest using a uniform,
    # but finer radial mesh than what you plan for Rayleigh.
    # Rayleigh's radial domain bounds should match or fall 
    # within the domain bounds used for this radial grid.
    nr = 2048                     # number of radial points
    ri = 0.5                      # Inner radius
    ro = 1.0                      # Outer radius [aspect ratio = 0.5]
    radius=numpy.linspace(ri,ro,nr, dtype='float64')

    #Instantiate an equation_coefficients object
    eqc = equation_coefficients(radius)

    # Set the buoyancy, heating, and viscosity functions
    # These particular choices may be questionable!
    buoy = radius
    nu   = radius**2
    heat = radius**3
    eqc.set_function(buoy , 2)  # set function 2
    eqc.set_function(nu   , 3)  # set function 3
    eqc.set_function(heat , 6)  # set function 6

    # Set the corresponding constants
    cbuoy = 10.0
    cnu   = 20.0
    cheat = 30.0
    eqc.set_constant(cbuoy , 2)   # set constant 2
    eqc.set_constant(cnu   , 5)   # set constant 5
    eqc.set_constant(cheat , 10)  # set constant 10

    #Generate the coefficients file
    eqc.write(ofile)
    
     

Constant Coefficients: Runtime Control
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

While constant coefficients may be specified via the coefficients file, many of these coefficients represent simulation
"control knobs" that the user may wish to modify at run-time.   For instance, the user may want to frequently use a particular profile for viscous diffusion (:math:`f_3`), but would like to vary its amplitude (:math:`c_5`) between simulations without generating a new coefficients file.  Rayleigh provides the opportunity to override 
all constant coefficients, or a subset of them, through the main_input file.

Consider the example below.

::

   &Reference_Namelist
    ...
    custom_reference_file='mycoeffs.dat'
    override_constants=T
    ra_constants( 2) = 1.0
    ra_constants( 5) = 10.0
    ra_constants(10) = 14.0
    ...
   /

In this example, values of constant coefficients :math:`c_2,\,c_5,\,c_{10}` will be determined entirely via the main_input file and assigned the values
of 1.0, 10.0, and 14.0 respectively.  Values specified in mycoeffs.dat will be ignored completely.   

This behavior is dictated by the override_constants flag, which instructs Rayleigh to ignore ALL constant coefficients specified in the coefficients files.   If a coefficient is not specified in
main_input, its value will be set to Rayleigh's internal default value of 0.   Consider the following example

::

   &Reference_Namelist
    ...
    custom_reference_file='mycoeffs.dat'
    override_constants=T
    ra_constants( 2) = 1.0
    ra_constants(10) = 14.0
    ...
   /

The resulting values of :math:`c_2,\,c_5,\,c_{10}` will be 1.0, 0.0, and 14.0 respectively.  The constant :math:`c_5` will not be set to 20.0 (the value specified in the coefficients file).


To specify a subset of constants, use the override_constant flag for each constant you wish to override, as shown below.

::

   &Reference_Namelist
    ...
    custom_reference_file='mycoeffs.dat'
    override_constant( 2) = T
    override_constant(10) = T
    ra_constants( 2) = 1.0
    ra_constants(10) = 14.0
    ...
   /

In this case, the values of constants :math:`c_2` and :math:`c_{10}` will be taken the main_input file.  The value of :math:`c_5` will be taken from the coefficients file.   If a constant's override flag is set, but its value is not specified in main_input, the default value of zero will be used. 

Augmenting a Rayleigh-Provided Reference State
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When augumenting one of Rayleigh's internal reference-state types, set the with_custom_reference flag (Reference_Namelist) to true in main_input.   In addition, assign a list of values to with_custom_constants and with_custom_functions.  As an example, to modify the heating and buoyancy profiles using entirely information provided through the equation coefficients file, main_input would contain the following
::

   &Reference_Namelist
    ...
    reference_type=1
    custom_reference_file='mycoeffs.dat'
    with_custom_reference=T
    with_custom_constants=2,10
    with_custom_functions=2,6
    ...
   /

These flags can be used in tandem with the override flags to specify values via main_input.  For example, the following input combination would set a value of :math:`c_2` of 13.0
::

   &Reference_Namelist
    ...
    reference_type=1
    custom_reference_file='mycoeffs.dat'
    with_custom_reference=T
    with_custom_constants=2,10
    with_custom_functions=2,6
    override_constant(2)=T
    ra_constants(2) = 13.0
    ...
   /



Specifing an Entire Custom Reference State
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To specify a full set of custom equation coefficients, set reference_type to 4.  Constant coefficients
may be overridden, if desired, and as described above.  Note that you must fully specify nonconstant coefficients :math:`f_1-f_7`.
If desired, you may also specify their logarithmic derivatives on the fine mesh (see the anelastic notebooks below).  This is optional,
however, as Rayleigh will compute those funtions if not provided.

::

   &Reference_Namelist
    ...
    reference_type = 4
    custom_reference_file='mycoeffs.dat'
    override_constant( 2) = T
    override_constant(10) = T
    ra_constants( 2) = 1.0
    ra_constants(10) = 14.0
    ...
   /



Behavior of Transport Coefficients
...........................................................
Transport coefficients may also be specified as desired, but nu_type, kappa_type, and eta_type still behave as described :ref:`below <transport>`.
If you wish to specify a custom diffusivity profile, set the corresponding type to 3.  In that case, the corresponding nonconstant coefficient MUST be set in the equation coefficients file.  Moreover, if reference_type=4, these corresponding constant must be set in either the coefficients file or in main_input (regardless of the diffusion type specified).  

For diffusion types 2 and 3, if the reference_type is not 4, the value of {nu,kappa,eta}_top normally used by that reference_type will be invoked if the corresponding constant coefficient is not set.

A Note on Volumetric Heating
......................
Finally, if specifying a custom form for the volumetric heating, please ensure that heating_type is set to a positive, nonzero value in the reference_namelist.  Otherwise, reference heating will be deactivated.  Any Rayleigh-initialization of the heating function that takes place initially will be overridden by the with_custom_reference or reference_type=4 flags. 


Custom Reference State Examples
~~~~~~~~~~~~~~~~

The notebooks below provide several examples of how to generate a custom-equation-coefficient file.
These notebooks are located in the examples/custom_reference_states subdirectory of the main Rayleigh directory.
Each notebook has an accompanying main_input file, also located in this directory.

.. toctree::
  :titlesonly:

  ../../../examples/custom_reference_states/Boussinesq_CZ.ipynb  


  ../../../examples/custom_reference_states/Boussinesq_Dynamo_Viscous.ipynb  


  ../../../examples/custom_reference_states/Anelastic_Dim_CZ.ipynb  


  ../../../examples/custom_reference_states/Anelastic_NonDim_CZ.ipynb  


  ../../../examples/custom_reference_states/MESA-input-1Msun-ZAMS.ipynb


.. _boundary_conditions:

Boundary Conditions & Internal Heating
--------------------------------------

Boundary conditions are controlled through the
**Boundary_Conditions_Namelist**. All Rayleigh simulations are run with
impenetrable boundaries. These boundaries may be either no-slip or
stress-free (default). If you want to employ no-slip conditions at both
boundaries, set **no_slip_boundaries = .true.**. If you wish to set
no-slip conditions at only one boundary, set **no_slip_top=.true.** or
**no_slip_bottom=.true.** in the Boundary_Conditions_Namelist.

By default, magnetic boundary conditions are set to match to a potential field at
each boundary.

By default, the thermal anomoly :math:`\Theta` is set to a fixed value
at each boundary. The upper and lower boundary-values are specified by
setting **T_top** and **T_bottom** respectively in the
Boundary_Conditions_Namelist. Their defaults values are 1 and 0
respectively.

Alternatively, you may specify a constant value of
:math:`\partial\Theta/\partial r` at each boundary. This is accomplished
by setting the variables **fix_dTdr_top** and **fix_dTdr_bottom**.
Values of the gradient may be enforced by setting the namelist variables
**dTdr_top** and **dTdr_bottom**. Both default to a value of zero.

Generic Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Boundary conditions for temperature, :math:`T`, and the magnetic poloidal potential, :math:`C`,
may also be set using generic spectral input.  As with initial conditions (see :ref:`sec:generic_ic`)
this is done by generating a generic input file using the **rayleigh_spectral_input.py** script.
The file output from this script can then be applied using:

-  **fix_Tvar_top** and **T_top_file** (overrides **T_top** for a constant boundary condition)
   to set a fixed upper :math:`T` boundary condition
-  **fix_dTdr_top** and **dTdr_top_file** (overrides **dTdr_top**) to set a fixed upper :math:`T` gradient boundary condition
-  **fix_Tvar_bottom** and **T_bottom_file** (overrides **T_bottom**) to set a fixed lower :math:`T` boundary condition
-  **fix_dTdr_bottom** and **dTdr_bottom_file** (overrides **dTdr_bottom**) to set a fixed lower :math:`T` gradient boundary
   condition
-  **fix_poloidal_top** and **C_top_file** (overrides **impose_dipole_field**) to set a fixed upper :math:`C` boundary condition
-  **fix_poloidal_bottom** and **C_bottom_file** (overrides **impose_dipole_field**) to set a fixed lower :math:`C` boundary condition

For example, to set a :math:`C` boundary condition on both boundaries with modes (l,m) = (1,0) and (1,1) set to pre-calculated
values run:

::
   
   rayleigh_spectral_input.py -m 1 0 2.973662220170157 -m 1 1 0.5243368809294343+0.j -o ctop_init_bc
   rayleigh_spectral_input.py -m 1 0 8.496177771914736 -m 1 1 1.4981053740840984+0.j -o cbottom_init_bc

which will generate generic spectral input files `ctop_init_bc` and `cbottom_init_bc`.  Set these to be used as the boundary
conditions in `main_input` using:

::

   &Boundary_Conditions_Namelist
   fix_poloidalfield_top = .true.
   fix_poloidalfield_bottom = .true.
   C_top_file = 'ctop_init_bc'
   C_bottom_file = 'cbottom_init_bc'
   /

This can be seen being applied in `tests/generic_input`.

Internal Heating
~~~~~~~~~~~~~~~~

The internal heating function :math:`Q(r)` is activated and described by
two variables in the **Reference_Namelist**. These are **Luminosity**
and **heating_type**. Note that these values are part of the
**Reference_Namelist** and not the **Boundary_Conditions** namelist.
Three heating types (0,1, and 4) are fully supported at this time.
Heating type zero corresponds to no heating. This is the default.

**Heating_type=1:** This heating type is given by :

.. math::

   \label{eq:heating}
   %\frac{\partial \Theta}{\partial t}=\gamma\left( 1 -\frac{\hat{\rho}(r_\mathrm{max})\,\hat{T}(r_\mathrm{max})  }{\hat{\rho}(r)\, \hat{T}(r)} \right),
   Q(r)= \gamma\,\hat{\rho}(r)\, \hat{T}(r)

where :math:`\gamma` is a normalization constant defined such that

.. _eq_lum:

.. math::


   %4\pi r_o^2 \hat{\rho}\hat{T}\kappa(r)\frac{\partial \Theta}{\partial r}=\mathrm{Luminosity}
   4\pi\int_{r=r_\mathrm{min}}^{r=r_\mathrm{max}} Q(r)\,  r^2 dr = \mathrm{Luminosity}.

This heating profile is particularly useful for emulating radiative
heating in a stellar convection zone.

**Heating_type=4:** This heating type corresponds a heating that is
variable in radius, but constant in *energy density*. Namely

.. math:: \hat{\rho}\hat{T}\frac{\partial \Theta}{\partial t}=\gamma.

The constant :math:`\gamma` in this case is also set by enforcing
Equation eq_lum_.

**Note:** If internal heating is used in combination with **fix_dTdr_top**, then the value of :math:`\partial\Theta/\partial r` 
at the upper boundary is set by Rayleigh.  Any value for **dTdr_top** specified in main_input is ignored.  This is done to ensure consistency with the internal 
heating and any flux passing through the lower boundary due
to the use of a fixed-flux condition.  To override this behavior, set **adjust_dTdr_top** to .false. in the
**Boundary_Conditions** namelist.

.. _output_controls:

Output Controls
---------------

Rayleigh comes bundled with an in-situ diagnostics package that allows
the user to sample a simulation in a variety of ways, and at
user-specified intervals throughout a run. This package is comprised of
roughly 17,000 lines of code (about half of the Rayleigh code base). Here we will
focus on generating basic output, but we refer the user to the section 
:ref:`plotting` and :ref:`quantityCodes` for more information.

Rayleigh can compute the quantities listed in :ref:`quantityCodes` in a variety of
averages, slices, and spectra, which are collectively called data products. These 
are sorted by Rayleigh into directories as 
follows:

* G_Avgs: The quantity is averaged over the entire simulation volume.
* Shell_Avgs: The quantity is averaged over each spherical shell and output as a function of radius.
* AZ_Avgs: The quantity is averaged over longitude and output as a function of radius and latitude.
* Shell_Slices: The quantity at the specified radii is output as a function of latitude and longitude.
* Equatorial_Slices: The quantity at the specified latitudes is output as a function of radius and longitude.
* Meridional_Slices: The quantity at the specified longitudes is output as a function of radius and latitude.
* Spherical_3D: The qunatity over the entire domain. Careful -- these files can be quite large.
* Shell_Spectra: The quantity's spherical harmonic coefficents at the specified radii.
* Point_Probes: The quantity at a specified radius, latitude, and longtiude.
* SPH_Mode_Sampels: The quantity's spherical harmonic coefficents at the specified radii and :math:`\ell`.

In addition Rayleigh can output `Checkpoints`, which are the data required to restart 
Rayleigh and will be discussed in detail in :ref:`checkpointing`, and `Timings`, which 
contain information about the performance of the run.

Output in Rayleigh is controled through the `io_controls_namelist`. For each of the data 
products listed, the output is specified using the following pattern:

* \_values: The quantity codes desired (seperated by commas)
* \_frequency: The frequency in iterations those quantities will be output.
* \_nrec: Number of records to be combined into a single file.
* \_levels: Radial indicies at which the quantities will be output.
* \_indices: Latitudinal indicies at which the quantities will be output.
* \_ell: The spherical harmonic degree at which the quantities will be output.
* \_r, \_theta, \_phi: The radial, latitudinal, and longitudinal indicies at which the quantities will be output.

For example, if you wanted to output shell slice data for quantities 1, 2, 10, and 2711 at radial indicices 2 and 54 every 100 iterations and have 4 records per file, you would set

::

    shellslice_levels    = 2,54
    shellslice_values    = 1,2,10,2711
    shellslice_frequency = 100
    shellslice_nrec      = 4

Files output in this way will have the filename of their iteration.


.. _transport:

Transport coefficents
---------------------

Transport coefficients (viscosity, thermal diffusivity) are specified in the transport namelist.

   .. _table_anelastic_trans:


.. centered:: **Table. Anelastic Transport.**

Variables in the Transport_Namelist
that must be specified when running in dimensional anelastic mode. In
addition, **reference_type=2** must also be specified in the
Reference_Namelist.

   +-----------------------------------+-----------------------------------+
   | Variable [Default value]          | Description                       |
   +===================================+===================================+
   | nu_top [1.0]                      | kinematic viscosity at rmax,      |
   |                                   | :math:`\nu(rmax)`                 |
   +-----------------------------------+-----------------------------------+
   | nu_type [1]                       | determines whether :math:`\nu` is |
   |                                   | constant with radius (1) or       |
   |                                   | varies with density (2)           |
   +-----------------------------------+-----------------------------------+
   | nu_power [0.0]                    | exponent in :                     |
   |                                   | :math:`\nu(r) = \left( \frac{\rho |
   |                                   | (r)}{\rho(r=rmax)} \right)^       |
   |                                   | {nu\_power}`;                     |
   |                                   | use with nu_type=2                |
   +-----------------------------------+-----------------------------------+
   | kappa_top [1.0]                   | thermal diffusivity at rmax,      |
   |                                   | :math:`\kappa(rmax)`              |
   +-----------------------------------+-----------------------------------+
   | kappa_type [1]                    | determines whether :math:`\kappa` |
   |                                   | is constant with radius (1) or    |
   |                                   | varies with density (2)           |
   +-----------------------------------+-----------------------------------+
   | kappa_power [0.0]                 | exponent in :                     |
   |                                   | :math:`\kappa(r) = \left( \frac{\ |
   |                                   | rho(r)}{\rho(r=rmax)} \right)^    |
   |                                   | {kappa\_power}`;                  |
   |                                   | use with kappa_type=2             |
   +-----------------------------------+-----------------------------------+
   | eta_top [1.0]                     | magnetic diffusivity at rmax,     |
   |                                   | :math:`\eta(rmax)`                |
   +-----------------------------------+-----------------------------------+
   | eta_type [1]                      | determines whether :math:`\eta`   |
   |                                   | is constant with radius (1) or    |
   |                                   | varies with density (2)           |
   +-----------------------------------+-----------------------------------+
   | eta_power [0.0]                   | exponent in :                     |
   |                                   | :math:`\eta(r) = \left( \frac{    |
   |                                   | \rho(r)}{\rho(r=rmax)} \right)^   |
   |                                   | {eta\_power}`;                    |
   |                                   | use with eta_type=2               |
   +-----------------------------------+-----------------------------------+


.. _examples:

Examples from Recent Publications
---------------------------------

*A Solar-like Case*

This is the main_input file from Case 39 from:

`Hindman, Bradley W., Nicholas A. Featherstone, and Keith Julien. 2020. “Morphological 
Classification of the Convective Regimes in Rotating Stars.” The Astrophysical Journal 
898 (2): 120. https://doi.org/10.3847/1538-4357/ab9ec2.`



::
   
   &problemsize_namelist
   n_r = 64
   n_theta = 192
   nprow = 32
   npcol = 32
   rmin = 5.0d10
   rmax = 6.5860209d10
   /
   &numerical_controls_namelist
   /
   &physical_controls_namelist
   rotation  = .true.
   magnetism = .false.
   /
   &temporal_controls_namelist
   max_time_step = 1000.0d0
   max_iterations = 5000000
   checkpoint_interval = 50000
   quicksave_interval = 10000
   num_quicksaves = 4
   cflmin = 0.4d0
   cflmax = 0.6d0
   /
   &io_controls_namelist
   /
   &output_namelist
   !shellslice_levels    = 16,32,48,64,80,96,112
   !shellslice_values    = 1                                               ! Codes needed for standard output routines
   shellslice_levels    = 8,16,24,32,40,48,56,64,72,80,88,96,104,112,120
   shellslice_values    = 1,2,3,301,302,303,304,305,306,307,308,309,401,501,502,2701,2702,2703,2704,2705,2706,2707,2708,2709,2710,2711
   shellslice_frequency = 10000
   shellslice_nrec      = 1

   !shellspectra_values    = 1,2,3                                         ! Codes needed for standard output routines
   shellspectra_levels    = 16,32,48,64,80,96,112
   shellspectra_values    = 1,2,3,301,302,303,304,305,306,307,308,309,401,501,502,503,504,2701,2702,2703,2704,2705,2706,2707,2708,2709,2710,2711
   shellspectra_frequency = 10000
   shellspectra_nrec      = 1

   !azavg_values    = 1,2,3,201,202                                        ! Codes needed for standard output routines
   azavg_values    = 1,2,3,201,202,401,405,409,501,502,1433,1455,1470,1923,1935,1943,2701,2702,2703,2704,2705,2706,2707,2708,2709,2710,2711,2712,2713,2714,2715
   azavg_frequency = 1000
   azavg_nrec = 10

   !shellavg_values    = 1,2,3,501,502,1433,1455,1470,1923,1935            ! Codes needed for standard output routines
   shellavg_values    = 1,2,3,401,405,409,501,502,1433,1455,1470,1923,1935,2701,2702,2703,2704,2705,2706,2707,2708,2709,2710,2711,2712,2713,2714,2715
   shellavg_frequency = 100
   shellavg_nrec = 100

   !globalavg_values = 401,402,403,404,405,406,407,408,409,410,411,412      ! Codes needed for standard output routines
   globalavg_values = 401,402,403,404,405,406,407,408,409,410,411,412,413,417,421,2701,2702,2703,2704,2705,2706,2707
   globalavg_frequency = 100
   globalavg_nrec = 100

   !equatorial_values    = 1,3      					! Codes needed for standard output routines
   equatorial_values    = 1,2,3,4,5,6,201,203,301,302,303,304,305,306,307,308,309,401,501,502,503,504,2701,2702,2703,2704,2705,2706,2707,2708,2709,2710,2711
   equatorial_frequency = 10000
   equatorial_nrec      = 1

   full3d_values = 4
   full3d_frequency = 9000000
   /

   &Boundary_Conditions_Namelist
   no_slip_boundaries = .false.
   strict_L_Conservation = .false.
   dtdr_bottom = 0.0d0
   T_Top    = 0.0d0
   T_Bottom = 851225.7d0
   fix_tvar_top = .true.
   fix_tvar_bottom = .false.
   fix_dtdr_bottom = .true.
   /
   &Initial_Conditions_Namelist
   init_type = 7
   magnetic_init_type = -1
   mag_amp = 1.0d0
   temp_amp = 1.0d1
   temp_w = 0.01d4
   !restart_iter = 0	! restart from latest checkpoint of any flavor
   /
   &Test_Namelist
   /
   &Reference_Namelist
   reference_type = 2
   heating_type = 1
   luminosity = 3.846d33
   poly_n = 1.5d0
   poly_Nrho = 3.0d0
   poly_mass = 1.98891D33
   poly_rho_i = 0.18053428d0
   pressure_specific_heat = 3.5d8
   angular_velocity = 5.74d-6	! Sidereal period of 12.7 days (twice the sidereal Carrington rate)
   /
   &Transport_Namelist
   nu_top    = 4.d12
   kappa_top = 4.d12
   /

