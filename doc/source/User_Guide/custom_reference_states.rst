Custom reference states
=======================




New words

If desired, the constant and nonconstant equation coefficients denoted  :ref:`here <equations_solved>` may be completely or partially specified by the user.  This allows the user to specify a background state and/or nondimensionalization note supplied by Rayleigh.  Three use cases are presently supported.  

1.  Nonconstant coefficents are completely specified through an auxilliary input file (described below).  Constant coefficients may be defined within this file or, read in from main_input, or some combination of the two.
2.  One of Rayleigh's predefined reference states may be used.  If desired, the buoyancy and/or volumetric heating coefficients may be overridden using the contents of the auxilliary input file.
3.  One of Rayleigh's predefined reference states may be used, but the transport coefficients (:math:`\nu,\kappa,\eta`) or, equivalently, (:math:`c_5, f_3, c_6, f_5, c_7, f_7`)  may be specified through a combination of values set in main_input and the auxilliary input file.

In all cases, the radial various of transport coefficients, as specified by nu_type, kappa_type, and eta_type flags is respected.

Generating an Equation Coefficients File
------------------------------------------
Words

Case 1:  Specifying a Full set of Equation Coefficients
---------------------------------------------------------
More words.


Case 2:  Augumenting a Rayleigh-provided Reference State
...........................................................
Even more words

Case 3:  Specifying Custom Transport-Coefficient Profiles
...........................................................


Constants are first read in from main_input.
Constants default to zero if not set in main_input.
Constants my be set in main_input by specifying:
::

   &Reference_Namelist
    ra_constants(1) = 1.0d0
   /

Constants set in custom reference file take precedence over main_input unless:
override_constants = T
or
override_constant(x) = T



::

   &problemsize_namelist
    n_r = 48
    n_theta = 96
   /


Example Notebooks
...................

The notebooks below provide several examples of how to generate a custom-equation-coefficient file.
These notebooks are located in the examples/custom_reference_states subdirectory of the main Rayleigh directory.
Each notebook has an accompanying main_input file, also located in this directory.

.. toctree::
  :titlesonly:

  ../../../examples/custom_reference_states/Boussinesq_CZ.ipynb  


  ../../../examples/custom_reference_states/Boussinesq_Dynamo_Viscous.ipynb  


  ../../../examples/custom_reference_states/Anelastic_Dim_CZ.ipynb  


  ../../../examples/custom_reference_states/Anelastic_NonDim_CZ.ipynb  


  ../../../examples/custom_reference_states/MESA-input-1Msun-ZAMS.ipynb
