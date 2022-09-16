.. raw:: latex

   \clearpage

.. _run_rayleigh:

Running Rayleigh
================

After setting up a custom `main_input` file, now it is time to run the new model.
This section focuses on the basics of running a Rayleigh model.

.. _load_balance:

Load-Balancing
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


.. _checkpointing:

Checkpointing
-------------

We refer to saved states in Rayleigh as **checkpoints**. A single
checkpoint consists of 13 files when magnetism is activated, and 9 files
when magnetism is turned off. A checkpoint written at time step *X*
contains all information needed to advance the system to time step
*X+1*. Checkpoint filenames end with a suffix indicating the contents of
the file (see Table table_checkpoints_). Each
checkpoint filename possess a prefix as well. Files belonging to the
same checkpoint share the same prefix. A checkpoint file collection,
written at time step 10,000 would look like that shown in Table
table_checkpoints_.

  .. _table_checkpoints:

.. centered:: **Table. Checkpoints.**

Example checkpoint file collection for a
time step 10,000. File contents are indicated.

   +-----------------------------------+-----------------------------------+
   | Filename                          | Contents                          |
   +===================================+===================================+
   | 00010000_W                        | Poloidal Stream function (at time |
   |                                   | step 10,000)                      |
   +-----------------------------------+-----------------------------------+
   | 00010000_Z                        | Toroidal Stream function          |
   +-----------------------------------+-----------------------------------+
   | 00010000_P                        | Pressure                          |
   +-----------------------------------+-----------------------------------+
   | 00010000_S                        | Entropy                           |
   +-----------------------------------+-----------------------------------+
   | 00010000_C                        | Poloidal Vector Potential         |
   +-----------------------------------+-----------------------------------+
   | 00010000_A                        | Toroidal Vector Potential         |
   +-----------------------------------+-----------------------------------+
   | 00010000_WAB                      | Adams-Bashforth (A-B) terms for   |
   |                                   | radial momentum (W) equation      |
   +-----------------------------------+-----------------------------------+
   | 00010000_ZAB                      | A-B terms for radial vorticity    |
   |                                   | (Z) equation                      |
   +-----------------------------------+-----------------------------------+
   | 00010000_PAB                      | A-B terms for horizontal          |
   |                                   | divergence of momentum (dWdr)     |
   |                                   | equation                          |
   +-----------------------------------+-----------------------------------+
   | 00010000_SAB                      | A-B terms for Entropy equation    |
   +-----------------------------------+-----------------------------------+
   | 00010000_CAB                      | A-B terms for C-equation          |
   +-----------------------------------+-----------------------------------+
   | 00010000_AAB                      | A-B terms for A-equation          |
   +-----------------------------------+-----------------------------------+
   | 00010000_grid_etc                 | grid and time-stepping info       |
   +-----------------------------------+-----------------------------------+

These files contain all information needed to advance the system state
from time step 10,000 to time step 10,001. Checkpoints come in two
flavors, denoted by two different prefix conventions: **standard
checkpoints** and **quicksaves**. This section discusses how to generate
and restart from both types of checkpoints.

Standard Checkpoints
~~~~~~~~~~~~~~~~~~~

Standard checkpoints are intended to serve as regularly spaced restart
points for a given run. These files begin with an 8-digit prefix
indicating the time step at which the checkpoint was created.

The frequency with which standard checkpoints are generated can be
controlled by modifying the **checkpoint_interval** variable in the
**temporal_controls_namelist**. For example, if you want to generate a
checkpoint once every 50,000 time steps, you would modify your
main_input file to read:

::

   &temporal_controls_namelist
    checkpoint_interval = 50000  ! Checkpoint every 50,000 time steps
   /

The default value of checkpoint_interval is 1,000,000, which is
typically much larger than what you will use in practice.

Restarting from a checkpoint is accomplished by first assigning a value
of -1 to the variables **init_type** and/or **magnetic_init_type** in
the **initial_conditions_namelist**. In addition, the time step from
which you wish to restart from should be specified using the
**restart_iter** variable in the initial_conditions_namelist. The
example below will restart both the magnetic and hydrodynamic variables
using the set of checkpoint files beginning with the prefix 00005000.

::

   &initial_conditions_namelist
    init_type = -1             !Restart magnetic and hydro variables from time step 5,000
    magnetic_init_type = -1
    restart_iter = 5000
   /

When both values are set to -1, hydrodynamic and magnetic variables are
read from the relevant checkpoint files. Alternatively, magnetic and
hydrodynamic variables may be initialized separately. This allows you to
add magnetism to an already equilibrated hydrodynamic case, for
instance. The example below will initialize the system with a random
magnetic field, but it will read the hydrodynamic information (W,Z,S,P)
from a checkpoint created at time step 5,000:

::

   &initial_conditions_namelist
    init_type = -1            ! Restart hydro from time step 5,000
    magnetic_init_type = 7    ! Add a random magnetic field
    restart_iter = 5000
   /

In addition to specifying the checkpoint time step manually, you can
tell Rayleigh to simply restart using the last checkpoint written by
assigning a value of zero to restart_iter:

::

   &initial_conditions_namelist
    init_type = -1
    magnetic_init_type = 7
    restart_iter = 0        ! Restart using the most recent checkpoint
   /

In this case, Rayleigh reads the **last_checkpoint** file (found within
the Checkpoints directory) to determine which checkpoint it reads.

.. _quicksaves:

Quicksaves
~~~~~~~~~~

In practice, Rayleigh checkpoints are used for two purposes: (1)
guarding against unexpected crashes and (2) supplementing a run’s record
with a series of restart points. While standard checkpoints may serve
both purposes, the frequency with which checkpoints are written for
purpose (1) is often much higher than that needed for purpose (2). This
means that a lot of data culling is performed at the end of a run or, if
disk space is a particularly limiting factor, during the run itself. For
this reason, Rayleigh has a **quicksave** checkpointing scheme in
addition to the standard scheme. Quicksaves can be written with
high-cadence, but require low storage due to the rotating reuse of
quicksave files.

The cadence with which quicksaves are written can be specified by
setting the **quicksave_interval** variable in the
**temporal_controls_namelist**. Alternatively, the elapsed wall time (in
minutes) that passes between quicksaves may be controlled by specifying
the **quicksave_minutes** variable. If both quicksave_interval and
quicksave_minutes are specified, quicksave_minutes takes precedence.

What distinguishes quicksaves from standard checkpoints is that only a
specified number of quicksaves exist on the disk at any given time. That
number is determined by the value of **num_quicksaves**. Quicksave files
begin with the prefix *quicksave_XX*, where XX is a 2-digit code,
ranging from 1 through num_quicksaves and indicates the quicksave
number. Consider the following example:

::

   &temporal_controls_namelist
    checkpoint_interval = 35000  ! Generate a standard checkpoint once every 35,000 time steps
    quicksave_interval = 10000   ! Generate a quicksave once every 10,000 time steps
    num_quicksaves = 2           ! Keep only two quicksaves on disk at a time
   /

At time step 10,000, a set of checkpoint files beginning with prefix
quicksave_01 will be generated. At time step 20,000, a set of checkpoint
files beginning with prefix quicksave_02 will be generated. Following
that, at time step 30,000, another checkpoint will be generated, *but it
will overwrite the existing quicksave_01 files*. At time step 40,000,
the quicksave_02 files will be overwritten, and so forth. Because the
**num_quicksaves** was set to 2, filenames with prefix quicksave_03 will
never be generated.

Note that checkpoints beginning with an 8-digit prefix (e.g., 00035000)
are still written to disk regularly and are not affected by the
quicksave checkpointing. On time steps where a quicksave and a standard
checkpoint would both be written, only the standard checkpoint is
written. Thus, at time step 70,000 in the example above, a standard
checkpoint would be written, and the files beginning with quicksave_01
would remain unaltered.

Restarting from quicksave_XX may be accomplished by specifying the value
of restart_iter to be -XX (i.e., the negative of the quicksave you wish
to restart from). The following example shows how to restart the
hydrodynamic variables from quicksave_02, while also initializing a
random magnetic field.

::

   &initial_conditions_namelist
    init_type = -1         ! Restart hydro variables from a checkpoint
    magnetic_init_type = 7 ! Initialize a random magnetic field
    restart_iter = -2      ! Restart from quicksave number 2
   /

Note that the file last_checkpoint contains the number of last
checkpoint written. This might be a quicksave or a standard checkpoint.
Specifying a value of zero for restart_iter thus works with quicksaves
and standard checkpoints alike.

Checkpoint Logs
~~~~~~~~~~~~~~~

When checkpoints are written, the number of the most recent checkpoint
is appended to a file named **checkpoint_log**, found in the Checkpoints
directory. The checkpoint log can be used to identify the time step
number of a quicksave file that otherwise has no identifying
information. While this information is also contained in the grid_etc
file, those are written in unformatted binary and cumbersome to access
from the terminal command line.

An entry in the log of "00050000 02" means that a checkpoint was written
at time step 50,000 to quicksave_02. An entry lacking a two-digit number
indicates that a standard checkpoint was written at that time step. The
most recent entry in the checkpoint log always comes at the end of the
file.

.. _control_time:

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

.. _log_file:

The Log File
------------

Section needs to be written.

.. _io_control:

I/O Control
-----------

Some aspects of Rayleigh's I/O can be controlled through variables found in the io_controls namelist.


I/O Format Controls
~~~~~~~~~~~~~~~~~~~


By default, integer output is reported with 8 digits and padded with leading zeros.  This includes integer iteration
numbers reported to stdout at each timestep and integer-number filenames created through diagnostics and checkpointing
output.  If desired, the number of digits may be controlled through the **integer_output_digits** variable.   When
reading in a Checkpoint created with a different number of digits, set the **integer_input_digits** variable to an appropriate
value.  

At several points in the code, floating-point output is sent to stdout.  This output is formatted using scienific notation, with
three digits to the right of the decimal place.     The number of digits after the decimal can be controlled through the
**decimal_places** variable.

As an example, the following combination of inputs
::
   &temporal_controls_namelist
   checkpoint_interval=10
   /
   &io_controls_namelist
   integer_output_digits=5
   integer_input_digits=3
   decimal_places=5
   /
   &initial_conditions_namelist
   init_type=-1
   restart_iter=10
   /

would restart from checkpoint files with the prefix formatted as:

::

   Checkpoints/010_grid_etc.

It would generate status line, shell_slice output, and checkpoints formatted as:

::

   Iteration:  00033   DeltaT:  1.00000E-04   Iter/sec:  2.68500E+00
   Shell_Slices/00020
   Checkpoints/00020_grid_etc.

**Developer's Note:**  The format codes generated through the values of these three variables are declared (with
descriptive comments) in Controls.F90.   For integer variables that may take on a negative value, additional format codes with one extra
digit (for the negative sign) are also provided.


I/O Redirection
~~~~~~~~~~~~~~~

Rayleigh writes all text output (e.g., error messages, iteration
counter, etc.) to stdout by default. Different computing centers handle
stdout in different ways, but typically one of two path is taken. On
some machines, a log file is created immediately and updated
continuously as the simulation runs. On other machines, stdout is
buffered on-node and written to disk only when the run has terminated.

There are situations where it can be advantageous to have a regularly
updated log file whose update frequency may be controlled. This feature
exists in Rayleigh and may be accessed by assigning values to
**stdout_flush_interval** and **stdout_file** in the io controls
namelist.

::

   &io_controls_namelist
   stdout_flush_interval = 1000
   stdout_file = 'routput'
   /

Set stdout_file to the name of a file that will contain Rayleigh’s text
output. In the example above, a file named *routput* will be appear in
the simulation directory and will be updated periodically throughout the
run. The variable stdout_flush_interval determines how many lines of
text are buffered before they are flushed to routput. Rayleigh prints
time-step information during each time step, and so setting this
variable to a relatively large number (e.g., 100+) prevents excessive
disk access from occurring throughout the run. In the example above, a
text buffer flush will occur once 1000 lines of text have been
accumulated.

Changes in the time-step size and self-termination of the run will also
force a text-buffer flush. Unexpected crashes and sudden termination by
the system job scheduler do not force a buffer flush. Note that the
default value of stdout_file is **‘nofile’**. If this value is
specified, output will directed to normal stdout.

To save on disk space for logs of very long runs, the number of status outputs
can be reduced by specifying **statusline_interval** in the
**io_controls_namelist**. This causes only every n-th status line to be
written.


.. _run_termination:

Run Termination
---------------

Section needs to be written.
