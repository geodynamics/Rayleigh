.. raw:: latex

   \clearpage

.. _installation:

Compiling and Installing Rayleigh
=================================

A detailed explanation of the installation process may be found in the
root directory of the code repository at:

  Rayleigh/INSTALL.

We provide an abbreviated version of those instructions here.

Third-Party Dependencies
------------------------

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
-----------

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

Verifying Your Installation
---------------------------

Rayleigh comes with a benchmarking mode that helps you verify that the
installation is performing correctly. If you are running Rayleigh for
the first time, or running on a new machine, follow along with the
example in §\ :ref:`benchmarking`, that you receive an accurate benchmark report before running a custom
model.

.. _spack-setup:

Alternative: Installation using Spack
---------------------------------------

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
