.. raw:: latex

   \clearpage

.. _getting_started:

Getting Started
=================================

Accessing Rayleigh
------------------

You can download the Rayleigh source code from `Rayleigh's GitHub respository <https://github.com/geodynamics/Rayleigh>`_ .

Setting up a Rayleigh Development Environment
---------------------------------------------

When running Rayleigh on HPC resources, always compile the software with the recommended compiler and link against
libraries optimized for the architecture you are running on. We provide example
instructions for some common systems at `Installation Instructions for HPC systems`_.

When developing Rayleigh or editing its documentation, however, such optimizations are rarely necessary.  Instead, it is sufficient for the code and documentation to compile.  For this purpose, we recommend setting up a `conda environment`_ or using our `Docker container`_.  Instructions for setting up an environment on Linux and Mac OS are provided below.

Conda Environment
~~~~~~~~~~~~~~~~~

First, if you don't have Conda, you should download and install the version appropriate for your architecture `here. <https://docs.conda.io/en/latest/miniconda.html>`_

Once you have Conda installed, create a Conda environment using the environment files we provide in Rayleigh's main directory.

.. code-block:: bash

    conda config --set channel_priority strict
    conda env create -f environment.yml
    conda activate radev

Because Rayleigh installs a number of different packages the first of these lines simplifies the
search for a compatible selection of package versions. The second line creates a new environment
named "radev" and installs all necessary packages.  The third line activates the environment.
This command will likely take a while (a few minutes).

If you want to undo the change to the channel priority setting afterwards, you can run

.. code-block:: bash

    conda config --set channel_priority flexible

after installing the environment (only necessary if you depend on the more flexible conda
solver).

Configuration and Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Building the documentation is the same on Linux and Mac.

.. code-block:: bash

    conda activate radev
    cd /path/to/Rayleigh
    make doc

Once the documetation builds, you can access it by opening Rayleigh/doc/build/html/index.html in your web browser.

Building the code is again the same on Linux and Mac. Execute the following:

.. code-block:: bash

    conda activate radev
    cd /path/to/Rayleigh
    ./configure -openblas --FC=mpifort
    make

At this point, you can run "make install," and run the code using mpirun as you normally would (keep the radev environment active when doing this).



Docker Container
~~~~~~~~~~~~~~~

Docker provides a standardized way to build, distribute and run containerized environments on Linux, macOS, and Windows. To get started you should install Docker on your system following the instructions from `here <https://www.docker.com/get-started>`_. On Linux you can likely also install it from a distribution package (e.g., ``docker-io`` on Debian/Ubuntu). `Podman <https://podman.io/>`_ is an alternative runtime that can run Docker containers and can be used as a drop-in replacement for Docker.
Singularity/Apptainer are other available alternatives that are more commonly used on HPC systems.

Launching the container
^^^^^^^^^^^^^^^^^^^^^^^
You can launch our pre-built container that is hosted on Docker Hub from a terminal.
This container is set up to get used to Rayleigh not to run productive models with it.

.. code-block:: bash
    docker run -it --rm geodynamics/rayleigh:latest bash

This command will create a terminal inside the container and drop you in a directory
that contains a pre-compiled version of Rayleigh. You can run input examples or
tests by executing `rayleigh.opt` or `rayleigh.dbg` and look at the output files, but
all files will be deleted when you `exit` the container.

.. note:: If you use Apptainer/Singularity instead of docker you can keep the model output files, because Apptainer by default mounts the current directory into the container. The command to run Rayleigh inside the container is ``mpirun -np X apptainer exec geodynamics/rayleigh:latest rayleigh.opt``` (assuming you have a Rayleigh input file in the current directory).

We also provide a container with a development environment for Rayleigh that allows you to change the code, build the
documentation and the code, and to keep model outputs.
The following command is for GNU/Linux and macOS users.

.. code-block:: bash

   ./docker-devel
   # This runs the following command:
   # docker run -it --rm -v $HOME:/work -e HOSTUID=$UID -e HOSTGID=$GROUPS -e HOSTUSER=$USER geodynamics/rayleigh-devel-jammy:latest

This will give you a shell inside the container and mount your home directory at ``/work``. You can clone, configure, build, and run the code and analyze the outputs using Python inside the container. Any changes below ``/work`` will be reflected in your home directory. Any other changes to the container will be deleted once you exit the shell.

.. note:: Your user has ``sudo`` rights within the container. This allows to install packages using the ``apt`` command or modify the system in any other way.

Windows users should run the script ``docker-devel.bat`` instead.

Configuration and Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. note:: All these commands are run inside the Docker container and assume you have a copy of Rayleigh at ``$HOME/path/to/Rayleigh`` (which corresponds to ``/work/path/to/Rayleigh`` inside the container).

Building the documentation

.. code-block:: bash

    cd /work/path/to/Rayleigh
    make doc

Building the code

.. code-block:: bash

    cd /work/path/to/Rayleigh
    ./configure --with-fftw=/usr
    make

Updating the container
^^^^^^^^^^^^^^^^^^^^^^
On the first launch of the container, your local Docker engine will automatically download our pre-built container from Docker Hub. Subsequent launches will just use this container and will not check for updates. You can download a newer version of the container using the following command.

.. code-block:: bash

    docker pull geodynamics/rayleigh-devel-jammy:latest

Building the container
^^^^^^^^^^^^^^^^^^^^^^
.. note:: This step purely optional. You only need to do this if you cannot pull the container from Docker Hub or you want to modify the Dockerfile.

To build the container you have to run this command from your host system (i.e., not from inside the container).

.. code-block:: bash

   cd docker
   docker build -t geodynamics/rayleigh-devel-jammy:latest rayleigh-devel-jammy

You can check the newly built container is there using this command.

.. code-block:: bash

    docker images

Spack Environment
~~~~~~~~~~~~~~~~~

`Spack <https://github.com/spack/spack>`_ can be used to create a development environment to build the code in a local directory. First set up Spack using the instructions in :ref:`spack-setup`

Afterwards create a new environment, activate it and set the status of the Rayleigh package to development. We select ``$PWD`` as the path, so run this command from the base directory of your git clone.

.. code-block:: bash

    spack env create rayleigh
    spack env activate rayleigh
    spack add rayleigh@master
    spack develop -p "$PWD" rayleigh@master

A subsequent ``spack install`` will install necessary dependencies and build Rayleigh in the selected directory.

.. _install_rayleigh:

Installing Rayleigh
-------------------


A detailed explanation of the installation process may be found in the
root directory of the code repository at:

https://github.com/geodynamics/Rayleigh/blob/main/INSTALL.

We provide an abbreviated version of those instructions here.

Third-Party Dependencies
~~~~~~~~~~~~~~~~~~~~~~~

In order to compile Rayleigh, you will need to have MPI (Message Passing
Interface) installed along with a Fortran 2003-compliant compiler.
Rayleigh has been successfully compiled with GNU, Intel, IBM, AOCC, and
Cray compilers (PGI has not been tested yet). Rayleigh’s configure script
provides native support for the Intel, GNU, AOCC, and Cray compilers. See
Rayleigh/INSTALL for an example of configuration using the IBM compiler.

Rayleigh depends on the following third party libraries:

#. BLAS (Basic Linear Algebra Subprograms)

#. LAPACK (Linear Algebra PACKage)

#. FFTW 3.x (Fastest Fourier Transform in the West)

You will need to install these libraries before compiling Rayleigh. If
you plan to run Rayleigh on Intel processors, we suggest installing
Intel’s Math Kernel Library (MKL) in lieu of installing these libraries
individually. The Math Kernel Library provides optimized versions of
BLAS, LAPACK, and FFTW. It has been tuned, by Intel, for optimal
performance on Intel processors. At the time of this writing, MKL is
provided free of charge. You may find it
`here <https://software.intel.com/en-us/mkl>`__.

Compilation
~~~~~~~~~~~

Rayleigh is compiled using the standard Linux installation scheme of
configure/make/make-install. From within the Rayleigh directory, run
these commands:

#. **./configure** – See Rayleigh/INSTALL or run ./configure --help to
   see relevant options.

#. **make** – This produces the code. You can run **make -j** to build several
   files in parallel and speed up the build this way.

#. **make install** – This places the Rayleigh executables in
   Rayleigh/bin. If you would like to place them in (say)
   /home/my_rayleigh/bin, run configure as: **./configure
   –prefix=/home/my_rayleigh**, i.e., the executables will be placed in the
   **$(prefix)/bin** directory.

For most builds, two executables will be created: rayleigh.opt and
rayleigh.dbg. Use them as follows:

#. When running production jobs, use **rayleigh.opt**.

#. If you encounter an unexpected crash and would like to report the
   error, rerun the job with **rayleigh.dbg**. This version of the code
   is compiled with debugging symbols. It will (usually) produce
   meaningful error messages in place of the gibberish that is output
   when rayleigh.opt crashes.

If *configure* detects the Intel compiler, you will be presented with a
number of choices for the vectorization option. If you select *all*,
rayleigh.opt will not be created. Instead, rayleigh.sse, rayleigh.avx,
etc. will be placed in Rayleigh/bin. This is useful if running on a
machine with heterogeneous node architectures (e.g., Pleiades). If you
are not running on such a machine, pick the appropriate vectorization
level, and rayleigh.opt will be compiled using that vectorization level.

The default behavior of the **make** command is to build both the
optimized, **rayleigh.opt**, and the debug versions, **rayleigh.dbg**. As
described above, if Intel is used and *all* is selected, every version will
be compiled. To build only a single version, the **target=<target>** option
may be used at the **make** stage, for example:

#. **make target=opt** - build only the optimized version, **rayleigh.opt**

#. **make target=dbg** - build only the debug version, **rayleigh.dbg**

#. **make target=avx** - build only the AVX version, **rayleigh.avx**

When building a single target, the final name of the executable can be changed
with the **output=<output>** option during the **make install** command. For example,
to build the optimized version and name the executable **a.out**:

#. **make target=opt** - only build the optimized version

#. **make target=opt output=a.out install** - install the optimized version as **a.out**

Inspection of the **$(prefix)/bin** directory (specified at configure time with the -prefix
option) will show a new file named **a.out**.

If both the optimized version and the debug version have already been built, they
can be renamed at install time as:

#. **make** - build both optimized and debug version (or all versions)

#. **make target=opt output=a.out.opt install** - install and rename the optimized version

#. **make target=dbg output=a.out.dbg install** - install and rename the debug version

The **output** option is only respected when a particular **target** is specified. Running
**make output=a.out install** will install all **rayleigh.*** executables, they will not
be renamed.

.. _spack-setup:

Alternative: Installation using Spack
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Spack is a package management tool designed to support multiple versions and
configurations of software on a wide variety of platforms and environments. It can be used to build Rayleigh with different compilers and a custom set of libraries for MPI, LAPACK, and FFTW. It can automatically build dependencies itself or use those provided by the HPC environment.

To set up Spack in your environment follow the instructions in the `documentation <https://spack.readthedocs.io/en/latest/getting_started.html>`_. Add local `compilers <https://spack.readthedocs.io/en/latest/getting_started.html#compiler-configuration>`_ and `packages <https://spack.readthedocs.io/en/latest/getting_started.html#system-packages>`_ as desired.

The next step has only to be performed once to add the Rayleigh package repository. Run this from the base directory of the Rayleigh repository.

.. code-block:: bash

    spack repo add spack-repo

Afterwards you can just install Rayleigh and its dependencies using:

.. code-block:: bash

    spack install rayleigh

Once the build succeeded the package can be loaded using the following command, which will make the ``rayleigh.opt`` and ``rayleigh.dbg`` executables available in the ``PATH`` and can be run to start simulations as usual.

.. code-block:: bash

    spack load rayleigh

There are many ways in which to modify the compiler and dependencies being used. They can be found in the `Spack documentation <https://spack.readthedocs.io/en/latest/index.html>`_.

As an example you can install Rayleigh using MKL for LAPACK and FFTW using:

.. code-block:: bash

    spack install rayleigh ^intel-mkl

To see the dependencies being installed you can use:

.. code-block:: bash

    spack spec rayleigh ^intel-mkl

.. _benchmark:


.. _hpc_installation_instructions:

Installation on HPC systems
-----------------------------------------

Given the amount of computational resources required to simulate convection in highly turbulent parameter regimes, many users will want to run Rayleigh in a HPC environment.  Here we provide instructions for compilation on two widely-used, national-scale supercomputing systems:  TACC Stampede2 and NASA Pleiades.   

Example jobscripts containing the necessary commands to compile and run Rayleigh on various systems may be found in *Rayleigh/job_scripts/*.

.. _stampede2:

TACC Stampede2
~~~~~~~~~~~~~~

Installing Rayleigh on NSF's Stampede 2 system is straightforward. At the time
this documentation is written (Sep 2022) the loaded default modules work out of
the box for Rayleigh. In case the modules change in the future here is a listed
for reference:

.. code-block:: bash

  1) intel/18.0.2      3) impi/18.0.2   5) autotools/1.1    7) cmake/3.20.2   9) TACC
  2) libfabric/1.7.0   4) git/2.24.1    6) python2/2.7.15   8) xalt/2.10.37

After cloning a Rayleigh repository, rayleigh can be configured and compiled as:

.. code-block:: bash

   FC=mpifc CC=mpicc ./configure  # select 'AVX512'
   make -j
   make install

We suggest choosing 'AVX512' at the configure menu.  This vectorization is supported by both the Skylake and Ice Lake nodes available on Stampede2.  An example jobscript for Stampede2 may be found in *Rayleigh/job_scripts/TACC_Stampede2*.

Using the Apptainer container system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We provide a precompiled container that provides an alternative way to use Rayleigh on the  TACC computing systems Stampede3 and Frontera.

To activate the container system and download the container:

.. code-block:: bash

   module load tacc-apptainer
   apptainer pull docker://gassmoeller/rayleigh:tacc-latest

This will create a file `rayleigh_tacc-latest.sif` that you can think of as a Rayleigh executable.
To run Rayleigh models using the downloaded container run (on a compute node):

.. code-block:: bash

   ibrun -n N_CORES apptainer run rayleigh_tacc-latest.sif rayleigh.opt

Make sure to replace N_CORES with the number of requested cores (or remove -n option to run with the
total number of cores requested). Also make sure to provide the correct path to `rayleigh_tacc-latest.sif`
if it is not in the current directory.

.. _pleiades:

NASA Pleiades
~~~~~~~~~~~~~

Installation on NASA's Pleiades cluster is similarly straightforward.  After cloning the repository, Rayleigh can be configured and compiled via the following commands:

.. code-block:: bash

   module purge
   module load comp-intel
   module load mpi-hpe
   ./configure --FC=mpif90 --CC=icc  # select 'ALL'
   make -j
   make install
   
We suggest using the default Intel and MPI compilers provided by Pleiades as in the example above.  As of December, 2022, this corresponded to the following version combination:

.. code-block:: bash

   1) comp-intel/2020.4.304   2) mpi-hpe/mpt.2.25

Note that Pleiades is a heterogeneous cluster, composed of many (primarily Intel) processor types. We suggest selecting the 'ALL' option when configuring Rayleigh to ensure that a unique executable is created for each of the possible vectorization options.  An example jobscript for Pleiades may be found in *Rayleigh/job_scripts/NASA_Pleiades*.



Verifying Your Installation
-------------------

Rayleigh has been programmed with internal testing suite so that its
results may be compared against benchmarks described in Christensen et al. (2001)
:cite:`CHRISTENSEN200125` and Jones et al. (2011)
:cite:`JONES2011120`

We recommend running a benchmark whenever running Rayleigh on a new
machine for the first time, or after recompiling the code. The
Christensen et al. (2001) :cite:`CHRISTENSEN200125` reference describes two Boussinesq tests that
Rayleigh’s results may be compared against. The Jones et al. (2011) :cite:`JONES2011120`
reference describes anelastic tests. Rayleigh has been tested
successfully against two benchmarks from each of these papers. Input
files for these different tests are enumerated in Table table_benchmark_
below. In addition to the
input files listed in Table table_benchmark_,
input examples appropriate for use as a template for new runs are
provided with the *\_input* suffix (as opposed to the *minimal* suffix.
These input files still have benchmark_mode active. Be sure to turn this
flag off if not running a benchmark.



**Important:** If you are not running a benchmark, but only wish to
modify an existing benchmark-input file, delete the line containing the
text “*benchmark_mode=X*.” When benchmark mode is active, custom inputs,
such as Rayleigh number, are overridden and reset to their
benchmark-appropriate values. For example, setting ``benchmark_mode = 1`` defines the appropriate Case 0
Christensen et al. (2001) :cite:`CHRISTENSEN200125` initial conditions. A benchmark report is
written every 5000 time steps by setting
``benchmark_report_interval = 5000``. The benchmark reports are text
files found within directory **path_to_my_sim/Benchmark_Reports/** and
numbered according to the appropriate time step. The
| ``benchmark_integration_interval`` variable sets the interval at which
measurements are taken to calculate the values reported in the
benchmark reports.

**We suggest using the c2001_case0_minimal input file for installation
verification**. Algorithmically, there is little difference between the
MHD, non-MHD, Boussinesq, and anelastic modes of Rayleigh. As a result,
when installing the code on a new machine, it is normally sufficient to
run the cheapest benchmark, case 0 from Christensen 2001 :cite:`CHRISTENSEN200125`.

To run this benchmark, create a directory from within which to run your
benchmark, and follow along with the commands below. Modify the
directory structure a each step as appropriate:

#. mkdir path_to_my_sim

#. cd path_to_my_sim

#. cp
   path_to_rayleigh/Rayleigh/input_examples/c2001_case0_minimal   main_input

#. cp path_to_rayleigh/Rayleigh/bin/rayleigh.opt   rayleigh.opt (or use
   *ln -s* in lieu of *cp*)

#. mpiexec -np **N** ./rayleigh.opt -nprow **X** -npcol **Y** -nr **R**
   -ntheta **T**

For the value **N**, select the number of cores you wish to run with.
For this short test, 32 cores is more than sufficient. Even with only
four cores, the lower-resolution test suggested below will only take
around half an hour. The values **X** and **Y** are integers that
describe the process grid. They should both be at least 2, and must
satisfy the expression

.. math:: N=X \times Y.

Some suggested combinations are {N,X,Y} = {32,4,8}, {16,4,4}, {8,2,4},
{4,2,2}. The values **R** and **T** denote the number of radial and
latitudinal collocation points respectively. Select either {R,T}={48,64}
or {R,T}={64,96}. The lower-resolution case takes about 3 minutes to run
on 32 Intel Haswell cores. The higher-resolution case takes about 12
minutes to run on 32 Intel Haswell cores.

Once your simulation has run, examine the file
path_to_my_sim/Benchmark_Reports/00025000. You should see output similar
to that presented in Tables table_benchmark_low_ or table_benchmark_high_ . Your numbers may differ
slightly, but all values should have a % difference of less than 1. If
this condition is satisfied, your installation is working correctly.

.. _table_benchmark_low:

.. centered:: **Benchmark Low**

Rayleigh benchmark report for Christensen
et al. (2001) :cite:`CHRISTENSEN200125` case 0 when run with nr=48 and ntheta=64. Run time was
approximately 3 minutes when run on 32 Intel Haswell cores.

Run command:

.. code-block::

 mpiexec -np 32 ./rayleigh.opt -nprow 4 -npcol 8 -nr 48 -ntheta 64

+-----------------+------------+------------+--------------+-----------+
| Observable      | Measured   | Suggested  | % Difference | Std. Dev. |
+=================+============+============+==============+===========+
| Kinetic Energy  | 58.347827  | 58.348000  | -0.000297    | 0.000000  |
+-----------------+------------+------------+--------------+-----------+
| Temperature     | 0.427416   | 0.428120   | -0.164525    | 0.000090  |
+-----------------+------------+------------+--------------+-----------+
| Vphi            | -10.118053 | -10.157100 | -0.384434    | 0.012386  |
+-----------------+------------+------------+--------------+-----------+
| Drift Frequency | 0.183272   | 0.182400   | 0.477962     | 0.007073  |
+-----------------+------------+------------+--------------+-----------+


.. _table_benchmark_high:


.. centered:: **Benchmark High**

Rayleigh benchmark report for Christensen
et al. (2001) :cite:`CHRISTENSEN200125` case 0 when run with nr=64 and ntheta=96. Run time was
approximately 12 minutes when run on 32 Intel Haswell cores.

Run command:

.. code-block::

  mpiexec -np 32 ./rayleigh.opt -nprow 4 -npcol 8 -nr 64 -ntheta 96

+-----------------+------------+------------+--------------+-----------+
| Observable      | Measured   | Suggested  | % Difference | Std. Dev. |
+=================+============+============+==============+===========+
| Kinetic Energy  | 58.347829  | 58.348000  | -0.000294    | 0.000000  |
+-----------------+------------+------------+--------------+-----------+
| Temperature     | 0.427786   | 0.428120   | -0.077927    | 0.000043  |
+-----------------+------------+------------+--------------+-----------+
| Vphi            | -10.140183 | -10.157100 | -0.166551    | 0.005891  |
+-----------------+------------+------------+--------------+-----------+
| Drift Frequency | 0.182276   | 0.182400   | -0.067994    | 0.004877  |
+-----------------+------------+------------+--------------+-----------+

.. _available_benchmarks:

Available Benchmarks
------------------



.. _table_benchmark:

.. centered:: **Benchmark**

Benchmark-input examples useful for verifying Rayleigh’s installation.
Those from Christensen et al. (2001) :cite:`CHRISTENSEN200125`
are Boussinesq. Those from Jones et al. (2011) :cite:`JONES2011120` are anelastic. Examples are found
in the directory: Rayleigh/input_examples/

+-----------------------+-----------------+--------------------------------+--------------------------------+
| Paper                 | Benchmark       | Input File                     | Specify in the main_input file |
+=======================+=================+================================+================================+
| Christensen et al.    | Case 0          | c2001_case0_minimal            | benchmark_mode = 1             |
+-----------------------+-----------------+--------------------------------+--------------------------------+
| Christensen et al.    | Case 1(MHD)     | c2001_case1_minimal            | benchmark_mode = 2             |
+-----------------------+-----------------+--------------------------------+--------------------------------+
| Jones et al. 2011     | Steady Hydro    | j2011_steady_hydro_minimal     | benchmark_mode = 3             |
+-----------------------+-----------------+--------------------------------+--------------------------------+
| Jones et al. 2011     | Steady MHD      | j2011_steady_MHD_minimal       | benchmark_mode = 4             |
+-----------------------+-----------------+--------------------------------+--------------------------------+
| Breuer et al. 2010    | Case 0          | b2010_case0_*T_input           |                                |
+-----------------------+-----------------+--------------------------------+--------------------------------+

Standard benchmarks that generate minimal output files are discussed in the next four
benchmarks:

* :ref:`cookbookCase0Minimal`
* :ref:`cookbookCase1Minimal`
* :ref:`cookbookHydroAnelastic`
* :ref:`cookbookMhdAnelastic`
* :ref:`cookbookDDCBreuer` 
.. _cookbookCase0Minimal:

Boussinesq non-MHD Benchmark: c2001_case0_minimal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is the standard benchmark test when running Rayleigh on a new
machine.
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


.. _cookbookCase1Minimal:

Boussinesq MHD Benchmark: c2001_case1_minimal
~~~~~~~~~~~~~~~~~~~~~~~

The MHD Boussinesq benchmark with an insulating inner core of
Christensen et al. (2001) :cite:`CHRISTENSEN200125` is denoted as Case 1 and is specified with
input file **c2001_case1_minimal**. Only the namelists modified compared
to Case 0 (\ :ref:`cookbookCase0Minimal` above) are shown
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



.. _cookbookHydroAnelastic:

Steady Anelastic non-MHD Benchmark: j2011_steady_hydro_minimal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

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



.. _cookbookMhdAnelastic:

Steady Anelastic MHD Benchmark: j2011_steady_mhd_minimal
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The anelastic MHD benchmark described in Jones et al. (2011) can be run
with main input file **j2011_steady_mhd_minimal**. The Fortran namelists
differing from the Jones et al. (2011) anelastic hydro benchmark
(§:ref:`cookbookHydroAnelastic` above) are shown here.

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

.. _cookbookDDCBreuer:

Steady Thermal-Chemical Boussinesq Convection Benchmark: b2010_case0_*T_input
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This is a Boussinesq convection benchmark described in Breuer et al. (2010) :cite:`Breuer2010`
in a dual buoyancy system that allows both thermal and chemical buoyancy sources. 
The case 0 contains three input lists that describes varying contributions of 
thermal vs chemical Rayleigh numbers whereas the total Rayleigh number stays the same. 
This benchmark is specified for Rayleigh with input file b2010_case0_*T_input. 
Below is an example for 80% thermal and 20% chemical convection scene for the 
relevant Fortran namelists:

::

   &problemsize_namelist
    n_r = 128
    n_theta = 192
    nprow = 32
    npcol = 16 
   &Reference_Namelist
    Ekman_Number = 1.0d-3
    Rayleigh_Number = 4.8d4
    Prandtl_Number = 3.0d-1
    chi_a_Rayleigh_Number = -1.2d5
    chi_a_Prandtl_Number = 3.0d0
