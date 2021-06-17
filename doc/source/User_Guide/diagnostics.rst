.. _diagnostics:

Diagnostic Outputs
==================

Rayleigh comes bundled with an in-situ diagnostics package that allows
the user to sample a simulation in a variety of ways, and at
user-specified intervals throughout a run. This package is comprised of
roughly 17,000 lines of code (about half of the Rayleigh code base), and
it is complex enough that we describe it in two other documents. We
refer the user to :

#. The diagnostics plotting manual, provided in two formats:

   -  Rayleigh/post_processing/Diagnostic_Plotting.ipynb (Jupyter Python
      notebook format; recommended for interactive use)

   -  Rayleigh/post_processing/Diagnostics_Plotting.html (recommended for optimal
      viewing; generated from the .ipynb file) `[html] <../../../post_processing/Diagnostic_Plotting.ipynb>`_


#. :ref:`DValues2` – This companion documentation
   provides the output menu system referred to in the main diagnostics
   documentation.

A number of stand-alone Python plotting examples may also be found in
the Rayleigh/post_processing/ directory.

The Lookup Table (LUT)
----------------------

Rayleigh has on the order of 1,000 possible diagnostic quantities available to the
user. As discussed in the examples above, the user specifies which diagnostic outputs
to compute by providing the appropriate quantity codes in the input file. Internally,
Rayleigh uses the quantity codes similarly to array indices. The purpose of the
lookup table is to map the quantity code to the correct position in the output data
array, you should never assume the quantities will be output in any particular order.
The user may have only requested two quantity codes, for example, 1 and 401.
The output data array will be of size 2 along the axis corresponding to the quantities.
The lookup table could map 401 to the first entry and 1 to the second entry.

The standard way to interact with the lookup table is to know the quantity code and
explicitly use it. Here we describe an alternative method. Each quantity code entry
(:ref:`quantityCodes`) has an equation, a code, and a name. There are some python
scripts in the post_processing directory that allow you to use the name, instead of
the code, when interacting with the lookup table:

    + lut.py
    + generate_mapping.py
    + lut_shortcuts.py

The lut.py file is the main user-interface and contains various utility routines,
including functions to convert between codes and names. The generate_mapping.py file
is responsible for generating the mapping between codes and names. The lut_shortcuts.py
allows users to define their own mapping, allowing a conversion from a user-defined name
to the desired quantity code.

The mapping has already been generated and is stored in the lut_mapping.py file. For
developers or anyone wanting to re-generate the mapping, use the generate_mapping.py file:

.. code-block:: bash

    python generate_mapping.py /path/to/Rayleigh

This will parse the Rayleigh directory tree and generate the standard mapping between
quantity codes and their associated names stored in the new file lut_mapping.py. Only
quantity codes that are defined within the Rayleigh source tree will be included.
Rayleigh does not need to be compiled before generating the mapping.

If a user has a custom directory where output diagnostics are defined, the above command
will not capture the custom diagnostic codes. To include custom quantities, the user
must generate the mapping themselvese with the generate_mapping.py file:

.. code-block:: bash

    python generate_mapping.py /path/to/Rayleigh/ --custom-dir=/path/to/custom/

Note that the Rayleigh directories are identical between the two calls, the only addition
is the custom-dir flag. This command will generate a new mapping stored in the file
lut_mapping_custom.py and will include all of the standard output quantities as well as
the custom diagnostics.

Without using this mapping technique, plotting something like the kinetic energy could
appear as:

.. code-block:: python

    from rayleigh_diagnostics import G_Avgs, build_file_list

    files = build_file_list(0, 10000000, path='G_Avgs')
    g = G_Avgs(filename=files[0], path='')

    ke_code = g.lut[401] # must use quantity code in lookup table

    ke = g.data[:, ke_code] # extract KE as a function of time

With the newly generated mapping, the above code could be rewritten as:

.. code-block:: python

    from rayleigh_diagnostics import G_Avgs, build_file_list

    from lut import lookup # <-- import helper function from main interface

    files = build_file_list(0, 10000000, path='G_Avgs')
    g = G_Avgs(filename=files[0], path='')

    ke_code = lookup('kinetic_energy') # use quantity *name* in lookup table

    ke = g.data[:, ke_code] # extract KE as a function of time, same as before

There is one drawback to using the quantity names: the naming scheme is somewhat
random and they can be quite long strings. This is where the lut_shortcuts.py
can be very useful. This allows users to define their own names to use in the mapping.
These are defined in the lut_shortcuts.py file and always take the form:

.. code-block:: python

    shortcuts['custom_name'] = 'rayleigh_name'

where custom_name is defined by the user, and rayleigh_name is the quantity name that
Rayleigh uses. The main dictionary must be named 'shortcuts'. With an entry like:

.. code-block:: python

    shortcuts['ke'] = 'kinetic_energy'

the above example for extracting the kinetic energy is even more simple:

.. code-block:: python

    from rayleigh_diagnostics import G_Avgs, build_file_list

    from lut import lookup # <-- import helper function from main interface

    files = build_file_list(0, 10000000, path='G_Avgs')
    g = G_Avgs(filename=files[0], path='')

    ke_code = lookup('ke') # user defined *name* in lookup table

    ke = g.data[:, ke_code] # extract KE as a function of time, same as before

