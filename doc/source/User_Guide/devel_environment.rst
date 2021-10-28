Setting up a Rayleigh Development Environment
==============================================

When running Rayleigh on HPC resources, always compile the software with the recommended compiler and link against
libraries optimized for the architecture you are running on.

When developing Rayleigh or editing its documentation, however, such optimizations are rarely necessary.  Instead, it is sufficient for the code and documentation to compile.  For this purpose, we recommend setting up a `conda environment`_ or using our `Docker container`_.  Instructions for setting up an environment on Linux and Mac OS are provided below.

Conda Environment
-----------------

First, if you don't have Conda, you should download and install the version appropriate for your architecture `here. <https://docs.conda.io/en/latest/miniconda.html>`_

Once you have Conda installed, create a Conda environment using the environment files we provide in Rayleigh's main directory.

.. code-block:: bash

    conda env create -f environment.yml
    conda activate radev

This command will likely take a while (a few minutes) and will install all necessary packages to compile Rayleigh.

MKL Setup: Linux and Mac
^^^^^^^^^^^^^^^^^^^^^^^^
Once your packages are installed, you will most likely want to have the ``MKLROOT`` environment variable set whenever you activate your Conda environment.  To do this we set ``MKLROOT`` to the location of the currently activated conda environment from the enviroment variable ``CONDA_PREFIX``.

.. code-block:: bash

    export MKLROOT="$CONDA_PREFIX"

Note that this is Bash syntax (use setenv if running c-shell).  Note that there should be no spaces on either side of the "=" sign.
If you stop here, you will have to do this every time you activate your development environment.   To have this happen automatically,
you only need to add two small scripts to radev/etc/conda/activate.d and radev/etc/conda/deactivate.d directories.   Scripts in these
directories are automatically executed when your conda environment is activated and deactivated, respectively.  

Change to your activate.d directory (for me, this was /custom/software/miniconda3/envs/radev/etc/conda/activate.d) and create a file named
activate_mkl.sh with the following three lines:

.. code-block:: bash

    #!/bin/bash
    export MKLSAVE="$MKLROOT"
    export MKLROOT="$CONDA_PREFIX"

In the deactivate.d directory, create a file named deactivate_mkl.sh with the following two lines:

.. code-block:: bash

    #!/bin/bash
    export MKLROOT="$MKLSAVE"

Now, try it out.

.. code-block:: bash

    conda deactivate
    echo $MKLROOT
    conda activate radev
    echo $MKLROOT

The MKLSAVE variable is used so that a separate MKL installation on your machine, if one exists,
is properly reset in your environment following deactivation.

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
    ./configure -conda-mkl --FC=mpifort
    make

At this point, you can run "make install," and run the code using mpirun as you normally would (keep the radev environment active when doing this).



Docker Container
----------------
Docker provides a standardized way to build, distribute and run containerized environments on Linux, macOS, and Windows. To get started you should install Docker on your system following the instructions from `here <https://www.docker.com/get-started>`_. On Linux you can likely also install it from a distribution package (e.g., ``docker-io`` on Debian/Ubuntu).

Launching the container
^^^^^^^^^^^^^^^^^^^^^^^
You can download our pre-built container from Docker Hub and launch it using the command from the main Rayleigh directory. The following command is for GNU/Linux and macOS users.

.. code-block:: bash

   ./docker-devel
   # This runs the following command:
   # docker run -it --rm -v $HOME:/work -e HOSTUID=$UID -e HOSTGID=$GROUPS -e HOSTUSER=$USER geodynamics/rayleigh-devel-bionic:latest

This will give you a shell inside the container and mount your home directory at ``/work``. You can clone, configure, build, and run the code and analyze the outputs using Python inside the container. Any changes below ``/work`` will be reflected in your home directory. Any other changes to the container will be deleted once you exit the shell.

.. note:: Your user has ``sudo`` rights within the container. This allows to install packages using the ``apt`` command or modify the system in any other way.

Windows users should run the script ``docker-devel.bat`` instead.

Configuration and Compilation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
.. note:: All these commands are run inside the Docker container and assume you have a copy of Rayleigh at ``$HOME/path/to/Rayleigh`` (which corresponds to ``/root/path/to/Rayleigh`` inside the container).

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

    docker pull geodynamics/rayleigh-devel-bionic:latest

Building the container
^^^^^^^^^^^^^^^^^^^^^^
.. note:: This step purely optional. You only need to do this if you cannot pull the container from Docker Hub or you want to modify the Dockerfile.

To build the container you have to run this command from your host system (i.e., not from inside the container).

.. code-block:: bash

   cd docker
   docker build -t geodynamics/rayleigh-devel-bionic:latest rayleigh-devel-bionic

You can check the newly built container is there using this command.

.. code-block:: bash

    docker images
