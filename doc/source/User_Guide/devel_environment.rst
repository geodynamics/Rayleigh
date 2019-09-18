Setting up a Rayleigh Development Environment
==============================================

When running Rayleigh on HPC resources, always compile the software with the recommended compiler and link against
libraries optimized for the architecture you are running on.

When developing Rayleigh or editing its documentation, however, such optimizations are rarely necessary.  Instead, it is sufficient for the code and documentation to compile.  For this purpose, we recommend setting up a conda environment.  Instructions for setting up an environment on Linux and Mac OS are provided below.   First, if you don't have Conda, you should download and install the version appropriate for your architecture `here. <https://docs.conda.io/en/latest/miniconda.html>`_

Once you have Conda installed, create a Conda environment named (say) radev

.. code-block:: bash

    conda create -n radev python=3
    conda activate radev

Once your environment is created and active, you are ready to install the packages required to compile Rayleigh and its documentation.  From here, the instructions for Linux and Mac differ slightly.

Package Setup:  Linux
-----------------------------

.. code-block:: bash

    conda activate radev  [ if you haven't done this already ]
    conda install -c conda-forge matplotlib \
        jupyter scipy sphinx sphinxcontrib-bibtex \
        nbsphinx pandoc recommonmark sphinx \
        gcc_linux-64 gfortran_linux-64 mkl fftw mpich

Package Setup:  Mac
-----------------------------

.. code-block:: bash

    conda activate radev  [ if you haven't done this already ]
    conda install -c conda-forge matplotlib \
        jupyter scipy sphinx sphinxcontrib-bibtex \
        nbsphinx pandoc recommonmark sphinx \
        clang_osx-64 gfortran_osx-64 mkl fftw mpich


MKL Setup: Linux and Mac
--------------------------
Once your packages are installed, you will most likely want to have the MKLROOT environment variable set whenever you activate your Conda environment.  To do this first, identify where Conda is located by running "which," and examining the output.

.. code-block:: bash

    which conda
    /custom/software/miniconda3/bin/conda   <<<  This is my output

The directory that we want to use for MKLROOT is located one level back, under the envs directory. In my case, I would set MKLROOT as

.. code-block:: bash

    export MKLROOT=/custom/software/miniconda3/envs/radev

Note that this is Bash syntax (use setenv if running c-shell).  Note that there should be no spaces on either side of the "=" sign.
If you stop here, you will have to do this every time you activate your development environment.   To have this happen automatically,
you only need to add two small scripts to radev/etc/conda/activate.d and radev/etc/conda/deactivate.d directories.   Scripts in these
directories are automatically executed when your conda environment is activated and deactivated, respectively.  

Change to your activate.d directory (for me, this was /custom/software/miniconda3/envs/radev/etc/conda/activate.d) and create a file named
activate_mkl.sh with the following three lines:

.. code-block:: bash

    #!/bin/bash
    export MKLSAVE=$MKLROOT
    export MKLROOT=/custom/software/miniconda3/envs/radev  [modify your path appropriately]

In the deactivate.d directory, create a file named deactivate_mkl.sh with the following two lines:

.. code-block:: bash

    #!/bin/bash
    export MKLROOT=$MKLSAVE

Now, try it out.

.. code-block:: bash

    conda deactivate
    echo $MKLROOT
    conda activate radev
    echo $MKLROOT

The MKLSAVE variable is used so that a separate MKL installation on your machine, if one exists,
is properly reset in your environment following deactivation.

Configuration and Compilation
-------------------------------
Building the documentation is the same on Linux and Mac.

.. code-block:: bash

    conda activate radev
    cd /path/to/Rayleigh
    make doc

Once the documetation builds, you can access it by opening Rayleigh/doc/build/html/index.html in your web browser.

Building the code different slightly on Linux and Mac.  For Linux, execute the following:

.. code-block:: bash

    conda activate radev
    cd /path/to/Rayleigh
    ./configure -conda-mkl --FC=mpifort
    make

For Mac, run:

.. code-block:: bash

    conda activate radev
    cd /path/to/Rayleigh
    ./configure -mac-mkl -conda-mkl --FC=mpifort
    make

At this point, you can run "make install," and run the code using mpirun as you normally would (keep the radev environment active when doing this).




