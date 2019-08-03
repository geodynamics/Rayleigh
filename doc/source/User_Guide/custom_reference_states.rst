Custom Reference States
=======================

If desired, the constant and nonconstant equation coefficients enumerated :ref:`here <equations_solved>` may be completely or partially specified by the user.  This allows the user to specify custom diffusivity profiles, background state, or nondimensionalization that is not supplied by Rayleigh.  Two use cases are currently supported.  

1.  Nonconstant coefficents are completely specified through an auxilliary input file (described below).  Constant coefficients may be defined within this file or, read in from main_input, or some combination of the two.

2.  One of Rayleigh's predefined reference states may be used.  For each predefined reference type, the auxilliary input file may be used to override volumetric heating (:math:`\,c_{10} + f_6\,`) and the transport coefficients :math:`\nu,\kappa,\eta` (:math:`\,c_5 + f_3,\,\, c_6 + f_5,\,\, c_7 + f_7\,` ).  For Boussinesq runs only, the buoyancy term may also be modified (:math:`\,c_2 + f_2\,`). 

In each use case, the radial variation of transport coefficients, as specified by nu_type, kappa_type, and eta_type flags is respected.

Create an Equation Coefficients File
.........................................

The first step in modifying Rayleigh's equation coefficients is to generate an equation coefficients file.
This file will be used alongside options defined in main_input to determine which combination of coefficients are overridden.

The equation_coefficient class records which constant and nonconstant coefficients have been set.  Rayleigh uses this info to...

.. code-block::

    import numpy
    from reference_tools import equation_coefficients

    #Define a name for your equation coefficients file
    ofile = 'my_coeffs.dat'

    # Define the radial grid.  We suggest using a uninform,
    # but finer radial mesh than what you plan for Rayleigh.
    nr = 2048                     # number of radial points
    ri = 0.5                      # Inner radius
    ro = 1.0                      # Outer radius [aspect ratio = 0.5]
    radius=numpy.linspace(ri,ro,nr)

    #Instantiate an equation_coefficients object
    eqc = equation_coefficients(radius)

    # Set the buoyancy, heating, and viscosity functions
    # These particular choices may be questionable!
    buoy = radius
    nu   = radius**2
    heat = radius**3
    eqc.set_function(buoy , 2)  # set function 2
    eqc.set_function(nu   , 3)  # set function 3
    eqc.set_function(heat , 6)  # set function 6

    # Set the corresponding constants
    cbuoy = 10.0
    cnu   = 20.0
    cheat = 30.0
    eqc.set_constant(cbuoy , 2)   # set constant 2
    eqc.set_constant(cnu   , 5)   # set constant 3
    eqc.set_constant(cheat , 10)  # set constant 6

    #Generate the coefficients file
    eqc.write(ofile)
    
     

Heating type note!

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
