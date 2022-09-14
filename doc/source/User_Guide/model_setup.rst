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

.. _numerical_controls:

Numerical Controls 
------------------

Rayleigh has several options that control aspects of the numerical method
used. For begining users these can generallly be left to default values.


.. _physics_controls:

Physics Controls
----------------

Many physical effects can be turned on or off in Rayleigh. The details of 
what physics you want to include will depend on the type of model you want 
to run. Be careful, however, that if you are adapting an input file from 
the benchmark described in :ref:`benchmark` that you set 
:code:`benchmark_mode` to 0 or omit it entirely, as this will override 
other input flags in favor of running the specified benchmark. For more
information on running benchmarks, see :ref:`cookbooks`.

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

Rayleigh can compute the quanties listed in :ref:`quantityCodes` in a variety of
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
* \_indices: Latitudinal indicies at which the quanties will be output.
* \_ell: The spherical harmonic degree at which the quantities will be output.
* \_r, \_theta, \_phi: The radial, latitudinal, and longitudinal indicies at which the quanties will be output.

For example, if you wanted to output shell slice data for quantities 1, 2, 10, and 2711 at radial indicices 2 and 54 every 100 iterations and have 4 records per file, you would set

::

    shellslice_levels    = 2,54
    shellslice_values    = 1,2,10,2711
    shellslice_frequency = 100
    shellslice_nrec      = 4

Files output in this way will have the filename of their iteration.

.. _examples:

Exampels from Recent Publications
---------------------------------

Need text here.

.. _surveys:

Ensemble Mode
-------------

Rayleigh can also be used to run multiple simulations under the umbrella
of a single executable. This functionality is particularly useful for
running parameter space studies, which often consist of mulitple,
similarly-sized simulations, in one shot. Moreover, as some queuing
systems favor large jobs over small jobs, an ensemble mode is useful for
advancing multiple small simulations through the queue in a reasonable
timeframe.

Running Rayleigh in ensemble mode is relatively straightforward. To
begin with, create a directory for each simulation as you normally
would, and place an appropriately modified main_input into each
directory. These directories should all reside within the same parent
directory. Within that parent directory, you should place a copy of the
Rayleigh executable (or a softlink). In addition, you should create a
text file named **run_list** that contains the name of each simulation
directory, one name per line. An ensemble job may then be executed by
calling Rayleigh with **nruns** command line flag as:

::

   user@machinename ~/runs/ $ mpiexec -np Y ./rayleigh.opt -nruns X

Here, Y is the total number of cores needed by all X simulations listed
in run_list.

**Example:** Suppose you wish to run three simulations at once from
within a parent directory named *ensemble* and that the simulation
directories are named run1, run2, and run3. When performing an *ls* from
within *ensemble*, you should see 5 items.

::

   user@machinename ~/runs/ $ cd ensemble
   user@machinename ~/runs/ensemble $ ls
   rayleigh.opt          run1          run2          run3          run_list

In this example, the contents of run_list should be the *local* names of
your ensemble run-directories, namely run1, run2, and run3.

::

   user@machinename ~runs/ensemble $ more run_list
   run1
   run2
   run3
             <--  place an empty line here

Note that some Fortran implementations will not read the last line in
run_list unless it ends in a newline character. Avoid unexpected crashes
by hitting "enter" following your final entry in run_list.

Before running Rayleigh, make sure you know how many cores each
simulation needs by examining the main_input files:

::

   user@machinename ~runs/ensemble $ head run1/main_input
   &problemsize_namelist
    n_r = 128
    n_theta = 192
    nprow = 16
    npcol = 16
   /

   user@machinename ~runs/ensemble $ head run2/main_input
   &problemsize_namelist
    n_r = 128
    n_theta = 384
    nprow = 32
    npcol = 16
   /

   user@machinename ~runs/ensemble $ head run3/main_input
   &problemsize_namelist
    n_r = 64
    n_theta = 192
    nprow = 16
    npcol = 16
   /

In this example, we need a total of 1024 cores (256+512+256) to execute
three simulations, and so the relevant call to Rayleigh would be:

::

   user@machinename ~/runs/ $ mpiexec -np 1024 ./rayleigh.opt -nruns 3

**Closing Notes:** When running in ensemble mode, it is *strongly
recommended* that you redirect standard output for each simulation to a
text file (see §\ :ref:`io`). Otherwise, all simulations
write to the same default (machine-dependent) log file, making it
difficult to read. Moreover, some machines such as NASA Pleiades will
terminate a run if the log file becomes too long. This is easy to do
when multiple simulations are writing to the same file.

Finally, The flags -nprow and -npcol **are ignored** when -nruns is
specified. The row and column configuration for all simulations needs to
be specified in their respective main_input files instead.

