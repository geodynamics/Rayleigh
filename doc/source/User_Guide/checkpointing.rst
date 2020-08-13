.. raw:: latex

   \clearpage

.. _checkpointing:

Checkpointing
=============

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
--------------------

Standard checkpoints are intended to serve as regularly spaced restart
points for a given run. These files begin with an 8-digit prefix
indicating the time step at which the checkpoint was created.

Generating Standard Checkpoints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Restarts From Standard Checkpoints
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
----------

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

Generating Quicksaves
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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

Restarting from Quicksaves
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
-------------------

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
