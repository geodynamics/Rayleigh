.. raw:: latex

   \clearpage

.. _troubleshooting:

Troubleshooting
===============

If you have questions that go beyond this manual, there are a number of
resources:

-   For questions on the source code of Rayleigh,
    portability, installation, new or existing features, etc., use the
    Rayleigh forum at
    <https://community.geodynamics.org/c/rayleigh>.

-   In case of more general questions about mantle convection, you can ask on
    the CIG mantle convection forum at
    <https://community.geodynamics.org/c/dynamo>.

-   If you have specific questions about Rayleigh
    that are not suitable for public and archived forums, you can contact the
    primary developers as listed at <https://github.com/geodynamics/Rayleigh/#more-information>.

.. _compile_error:

Compiling Errors
----------------

Need text here.

Using MKL in a Conda environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you use MKL inside a Conda environment (e.g. inside the ``radev`` environment set up by the "environment.yml" file inside Rayleigh),
you will most likely want to have the ``MKLROOT`` environment variable set whenever you activate your Conda environment.
To do this we set ``MKLROOT`` to the location of the currently activated conda environment from the enviroment variable ``CONDA_PREFIX``.

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

If you now run ``configure -conda-mkl`` Rayleigh should correctly pick up the conda MKL installation.

.. _seg_fault:

Segmentation Fault Crashes
--------------------------

Need test here.

.. _timestep_crash:

Timestep Crashes
----------------

Need text here.

.. _ringing:

Numerical Ringing
-----------------

Need text here.
