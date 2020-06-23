.. raw:: latex

   \clearpage

.. _running:

Running the Code
================

Whenever you run a new simulation, a similar series of steps must be
performed. A summary of the typical Rayleigh work flow is:

#. Create a unique directory for storing simulation output

#. Create a main_input file

#. Copy or soft link the Rayleigh executable into the simulation
   directory

#. Modify main_input as desired

#. Run the code

#. Examine output and restart simulation as necessary

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
outputs. This step is unnecessary when compiling with the Intel or GNU
compilers.

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

Code Execution and Load-Balancing
---------------------------------

Rayleigh is parallelized using MPI and a 2-D domain decomposition. The
2-D domain decomposition means that we envision the MPI Ranks as being
distributed in rows and columns. The number of MPI ranks within a row is
*nprow* and the number of MPI ranks within a column is *npcol*. When
Rayleigh is run with N MPI ranks, the following constraint must be
satisfied:

.. math:: \mathrm{N} = \mathrm{npcol} \times \mathrm{nprow}   .

If this constraint is not satisfied , the code will print an error
message and exit. The values of *nprow* and *npcol* can be specified in
*main_input* or on the command line via the syntax:

::

   mpiexec -np 8 ./rayleigh.opt -nprow 4 -npcol 2

Load Balancing
~~~~~~~~~~~~~~

Rayleigh’s performance is sensitive to the values of *nprow* and
*npcol*, as well as the number of radial grid points :math:`N_r` and
latitudinal grid points :math:`N_\theta`. If you examine the main_input
file, you will see that it is divided into Fortran namelists. The first
namelist is the problemsize_namelist. Within this namelist, you will see
a place to specify nprow and npcol. Edit main_input so that nprow and
npcol agree with the N you intend to use (or use the command-line syntax
mentioned above). The dominate effect on parallel scalability is the
number of messages sent per iteration. For optimal message counts, nprow
and npcol should be as close to one another in value as possible.

#. N = nprow :math:`\times` npcol.

#. nprow and npcol should be equal or within a factor of two of one
   another.

The value of nprow determines how spherical harmonics are distributed
across processors. Spherical harmonics are distributed in
high-\ :math:`m`/low-:math:`m` pairs, where :math:`m` is the azimuthal
wavenumber. Each process is responsible for all :math:`\ell`-values
associated with those :math:`m`\ ’s contained in memory.

The value of npcol determines how radial levels are distributed across
processors. Radii are distributed uniformly across processes in
contiguous chunks. Each process is responsible for a range of radii
:math:`\Delta r`.

The number of spherical harmonic degrees :math:`N_\ell` is defined by

.. math:: N_\ell = \frac{2}{3}N_\theta

For optimal load-balancing, *nprow* should divide evenly into
:math:`N_r` and *npcol* should divide evenly into the number of
high-\ :math:`m`/low-:math:`m` pairs (i.e., :math:`N_\ell/2`). Both
*nprow* and *npcol* must be at least 2.

In summary,

#. :math:`nprow \ge 2`.

#. :math:`npcol \ge 2`.

#. :math:`n \times npcol = N_r` (for integer :math:`n`).

#. :math:`k \times nprow = \frac{1}{3}N_\theta` (for integer :math:`k`).

Specifying Resolution & Domain Bounds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As discussed, the number of radial grid points is denoted by
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
:math:`rmin` (the lower radial boundary) and :math:`rmax` (the upper
radial boundary):

::

   &problemsize_namelist
    rmin = 1.0
    rmax = 2.0
   /

Alternatively, the user may specify the shell depth (:math:`rmax-rmin`)
and aspect ratio (:math:`rmin/rmax`) in lieu of :math:`rmin` and
:math:`rmax`. The preceding example may then be written as:

::

   &problemsize_namelist
    aspect_ratio = 0.5
    shell_depth = 1.0
   /

Note that the interpretation of :math:`rmin` and :math:`rmax` depends on
whether your simulation is dimensional or nondimensional. We discuss
these alternative formulations in §\ :ref:`physics`

Controlling Run Length & Time Stepping
--------------------------------------

A simulation’s runtime and time-step size can be controlled using the
**temporal_controls** namelist. The length of time for which a
simulation runs before completing is controlled by the namelist variable
**max_time_minutes**. The maximum number of time steps that a simulation
will run for is determined by the value of the namelist
**max_iterations**. The simulation will complete when it has run for
*max_time_minutes minutes* or when it has run for *max_iterations time
steps* – whichever occurs first.

An orderly shutdown of Rayleigh can be manually triggered by creating a file
with the name set in **terminate_file** (i.e., running the command *touch
terminate* in the default setting). If the file is found, Rayleigh will stop
after the next time step and write a checkpoint file. The existence of
**terminate_file** is checked every **terminate_check_interval** iterations.
The check can be switched off completely by setting
**terminate_check_interval** to -1. Both of these options are set in the
**io_controls_namelist**. With the appropriate job script this feature can be
used to easily restart the code with new settings without losing the current
allocation in the queuing system. A **terminate_file** left over from
a previous run is automatically deleted when the code starts.

Time-step size in Rayleigh is controlled by the Courant-Friedrichs-Lewy
condition (CFL; as determined by the fluid velocity and Alfvén speed). A
safety factor of **cflmax** is applied to the maximum time step
determined by the CFL. Time-stepping is adaptive. An additional variable
**cflmin** is used to determine if the time step should be increased.

The user may also specify the maximum allowed time-step size through the
namelist variable **max_time_step**. The minimum allowable time-step
size is controlled through the variable **min_time_step**. If the CFL
condition is less than this value, the simulation will exit.

Let :math:`\Delta t` be the current time-step size, and let
:math:`t_\mathrm{CFL}` be the maximum time-step size as determined by
the CFL limit. The following logic is employed by Rayleigh when
calculating the time-step size:

-  IF { :math:`\Delta_t\ge \mathrm{cflmax}\times t_\mathrm{CFL}` } THEN
   { :math:`\Delta_t` is set to
   :math:`\mathrm{cflmax}\times t_\mathrm{CFL}` }.

-  IF { :math:`\Delta_t\le \mathrm{cflmin}\times t_\mathrm{CFL}` } THEN
   { :math:`\Delta_t` is set to
   :math:`\mathrm{cflmax}\times t_\mathrm{CFL}` }.

-  IF{ :math:`t_\mathrm {CFL}\ge \mathrm{max\_time\_step}` } THEN {
   :math:`\Delta_t` is set to max_time_step }

-  IF{ :math:`t_\mathrm {CFL}\le \mathrm{min\_time\_step}` } THEN {
   Rayleigh Exits }

The default values for these variables are:

::

   &temporal_controls_namelist
   max_iterations = 1000000
   max_time_minutes = 1d8
   cflmax = 0.6d0
   cflmin = 0.4d0
   max_time_step = 1.0d0
   min_time_step = 1.0d-13
   /
