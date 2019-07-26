.. raw:: latex

   \clearpage

.. _cookbooks:

Cookbooks
=========

*This section has been written by Maria Weber.*

In this section, we explain and document example Rayleigh main input
files for a variety of problems. The example input files can be found in
directory **Rayleigh/input_examples/**. In some cases, we will also show
example diagnostic outputs. See
§\ :ref:`diagnostics` for more information on
generating Rayleigh diagnostic routines with Python.

Standard benchmarks that generate minimal output files are discussed in the next four
benchmarks:

* :ref:`cookbookCase0Minimal`
* :ref:`cookbookCase1Minimal`
* :ref:`cookbookHydroAnelastic`
* :ref:`cookbookMhdAnelastic`

.. _cookbookCase0Minimal:

Simple Boussinesq non-MHD benchmark: c2001_case0_minimal
--------------------------------------------------------

This is the standard benchmark test when running Rayleigh on a new
machine, as described in §\ :ref:`benchmarking`.
Christensen et al. (2001) :cite:`CHRISTENSEN200125` describes two Boussinesq tests that Rayleigh’s
results may be compared against. Case 0 in Christensen et al. (2001) :cite:`CHRISTENSEN200125`
solves for Boussinesq (non-dimensional) non-magnetic convection, and we
will discuss the input parameters necessary to set up this benchmark in
Rayleigh below. Rayleigh’s input parameters are grouped in so-called
namelists, which are subcategories of related input parameters that will
be read upon program start and assigned to Fortran variables with
identical names. Below are the first four Fortran namelists in the input
file **c2001_case0_minimal**.

::

   &problemsize_namelist
    n_r = 64
    n_theta = 96
    nprow = 16
    npcol = 32
   /
   &numerical_controls_namelist
   /
   &physical_controls_namelist
    benchmark_mode = 1
    benchmark_integration_interval = 100
    benchmark_report_interval = 5000
   /
   &temporal_controls_namelist
    max_iterations = 25000
    checkpoint_interval = 100000
    quicksave_interval = 10000
    num_quicksaves = 2
   /

In namelist ``problemsize_namelist``, the number of radial grid points
is denoted by ``n_r`` and the number of :math:`\theta` grid points by
``n_theta``. For optimal load-balancing, the number of MPI ranks
distributed within a row is denoted by ``nprow`` and within a column is
``npcol``. See §\ :ref:`running` for instructions on
appropriately defining these values.

| When running a benchmark, set ``benchmark_mode`` under
  ``physical_controls_namelist`` to the code number for the
  corresponding benchmark you want to run. When benchmark mode is
  active, custom inputs are overridden and reset to their benchmark
  appropriate value (see §\ :ref:`benchmarking`).
  Setting ``benchmark_mode = 1`` defines the appropriate Case 0
  Christensen et al. (2001) :cite:`CHRISTENSEN200125` initial conditions. A benchmark report is
  written every 5000 time steps by setting
  ``benchmark_report_interval = 5000``. The benchmark reports are text
  files found within directory **path_to_my_sim/Benchmark_Reports/** and
  numbered according to the appropriate time step. The
| ``benchmark_integration_interval`` variable sets the interval at which
  measurements are taken to calculate the values reported in the
  benchmark reports.

Within ``temporal_controls_namelist``, the maximum number of iterations
is set with ``max_interations``. Checkpoints are written at time step
intervals set by ``checkpoint_interval``. In this case, the checkpoint
interval is larger than the maximum number of iterations, so no
checkpoint will be written. The interval at which quicksaves are written
is set by variable ``quicksave_interval`` and the number of quicksaves
saved on disk at a time is set by ``num_quicksaves``. See
§\ :ref:`quicksaves` for more information on
quicksaves.

| Upon completion of this benchmark, verify that your installation is
  working correctly by comparing the file
| **path_to_my_sim/Benchmark_Reports/00025000** to Table Benchmark High in §\ :ref:`benchmarking`. All values should have a
  percent difference of less than 1.

.. _cookbookCase1Minimal:

Simple Boussinesq MHD benchmark: c2001_case1_minimal
----------------------------------------------------

The MHD Boussinesq benchmark with an insulating inner core of
Christensen et al. (2001) :cite:`CHRISTENSEN200125` is denoted as Case 1 and is specified with
input file **c2001_case1_minimal**. Only the namelists modified compared
to Case 0 (§\ :ref:`cookbookCase0Minimal` above) are shown
below.

::

   &physical_controls_namelist
    benchmark_mode = 2
    benchmark_integration_interval = 100
    benchmark_report_interval = 10000
   /
   &temporal_controls_namelist
    max_iterations = 150000
    checkpoint_interval = 100000
    quicksave_interval = 10000
    num_quicksaves = 2
   /

In this example, ``benchmark_mode = 2`` sets the benchmark-appropriate
values for Christensen et al. (2001) :cite:`CHRISTENSEN200125` Case 1. The variable
``benchmark_integration_interval`` remains the same as Case 0 above, but
the ``benchmark_report_interval`` has been increased in this MHD
problem. Here, ``max_iterations`` has also been increased compared to
Case 0 such that it is now larger than ``checkpoint_interval``. As such,
checkpoint files for time step 100000 will be written in directory
**path_to_my_sim/Checkpoints/00100000**. Upon completion of this
benchmark, verify that your installation is working correctly by looking
at file **path_to_my_sim/Benchmark_Reports/00150000**. All values should
have a percent difference of less than 1.

.. _cookbookHydroAnelastic:

Steady anelastic non-MHD benchmark: j2011_steady_hydro_minimal
--------------------------------------------------------------

Jones et al. (2011) describes a benchmark for an anelastic hydrodynamic
solution that is steady in a drifting frame. This benchmark is specified
for Rayleigh with input file **j2011_steady_hydro_minimal**. Below are
the relevant Fortran namelists.

::

   &problemsize_namelist
    n_r = 128
    n_theta = 192
    nprow = 32
    npcol = 16
   /
   &numerical_controls_namelist
   /
   &physical_controls_namelist
    benchmark_mode = 3
    benchmark_integration_interval = 100
    benchmark_report_interval = 10000
   /
   &temporal_controls_namelist
    max_iterations = 200000
    checkpoint_interval = 100000
    quicksave_interval = 10000
    num_quicksaves = 2
   /

Suggested problem size values are given in ``problemsize_namelist``,
along with variables for ``physical_controls_namelist`` and
``temporal_controls_namelist``. The variable ``benchmark_mode = 3``
designates appropriate input conditions for the Jones et al. (2011)
anelastic hydrodynamic benchmark. Upon completion of this benchmark,
verify that your installation is working correctly by looking at file
**path_to_my_sim/Benchmark_Reports/00200000**. All values should have a
percent difference of less than 1.

.. _cookbookMhdAnelastic:

Steady anelastic MHD benchmark: j2011_steady_mhd_minimal
--------------------------------------------------------

The anelastic MHD benchmark described in Jones et al. (2011) can be run
with main input file **j2011_steady_mhd_minimal**. The Fortran namelists
differing from the Jones et al. (2011) anelastic hydro benchmark
(§`:ref:cookbookHydroAnelastic` above) are shown here.

::

   &physical_controls_namelist
    benchmark_mode = 4
    benchmark_integration_interval = 100
    benchmark_report_interval = 10000
   /
   &temporal_controls_namelist
    max_iterations = 5000000
    checkpoint_interval = 100000
    quicksave_interval  = 25000
    num_quicksaves = 2
   /

You may wish to modify the problem size within ``problemsize_namelist``
(particularly ``nprow`` and ``npcol``), explained in more detail in
§\ :ref:`cookbookCase0Minimal`. The variable
``benchmark_mode = 4`` designates appropriate input conditions for the
Jones et al. (2011) anelastic MHD benchmark. Here, ``max_iterations``
has also been increased compared to the anelastic hydro benchmark of
Jones et al. (2011), as well as ``quicksave_interval``. Upon completion
of this benchmark, verify that your installation is working correctly by
looking at file **path_to_my_sim/Benchmark_Reports/05000000**. All
values should have a percent difference of less than 1.
