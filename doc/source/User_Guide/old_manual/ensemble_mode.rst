.. raw:: latex

   \clearpage

Ensemble Mode
=============

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
