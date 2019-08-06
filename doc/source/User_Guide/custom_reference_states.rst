Custom Reference States
=======================

If desired, the constant and nonconstant equation coefficients enumerated :ref:`here <equations_solved>` may be completely or partially specified by the user.  This allows the user to specify diffusivity profiles, background states, or nondimensionalizations that are not supplied by Rayleigh.  Two use cases are supported:  

1.  One of Rayleigh's predefined reference states may be used, but with some coefficients supplied instead through an auxilliary coefficients file.  For each predefined reference type, the coefficients file may be used to override volumetric heating (:math:`\,c_{10} + f_6\,`) and the transport coefficients :math:`\nu,\kappa,\eta` (:math:`\,c_5 + f_3,\,\, c_6 + f_5,\,\, c_7 + f_7\,` ).  For Boussinesq runs, the buoyancy term may also be modified (:math:`\,c_2 + f_2\,`). 

2.  Nonconstant coefficents may be completely specified through the coefficients file.  In this mode, activated by setting reference_type=4, the user must fully specify nonconstant coefficients :math:`f_1-f_7`.  


In either case, constant coefficients may be defined within the coefficients file or, read in from main_input, or some combination of the two.  Moreover, the radial variation of transport coefficients, as specified by nu_type, kappa_type, and eta_type flags is respected.  We elaborate on this behavior below.

Creating a Coefficients File
.........................................

The first step in modifying Rayleigh's equation coefficients is to generate an equation coefficients file.
This file will be used alongside options defined in main_input to determine which combination of coefficients are overridden.
In order create your coefficients file, you will need to create an instance of the equation_coefficients class, provided
in post_processing/reference_tools.py.  Constant and nonconstant coefficients may then be set through  set_constant and set_function methods respectively.

The equation_coefficient class is instantiated by passing a radial grid to its init method.  This grid can be cast in ascending or descending order, but it should generally possess a much finer mesh than what you plan to use in Rayleigh.
Nonconstant coefficients specified in the coefficients file will be interpolated onto the Rayleigh grid at input time.   

The file structure created through the class's write method contains a record of those functions and contants that have been set. Rayleigh uses this information at runtime along with main_input to to perform consistency checks and to determine the values ultimately assigned to each constant coefficient. 

The sample code below defines a file with sufficient information to alter the viscous, heating, and buoyancy functions of a Rayleigh-provided reference state.  This information would be insufficient for use with reference_type=4, but several example notebooks handling that scenario are provided below. 

.. code-block::

    import numpy
    from reference_tools import equation_coefficients

    #Define a name for your equation coefficients file
    ofile = 'my_coeffs.dat'

    # Define the radial grid.  We suggest using a uniform,
    # but finer radial mesh than what you plan for Rayleigh.
    # Rayleigh's radial domain bounds should match or fall 
    # within the domain bounds used for this radial grid.
    nr = 2048                     # number of radial points
    ri = 0.5                      # Inner radius
    ro = 1.0                      # Outer radius [aspect ratio = 0.5]
    radius=numpy.linspace(ri,ro,nr, dtype='float64')

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
    eqc.set_constant(cnu   , 5)   # set constant 5
    eqc.set_constant(cheat , 10)  # set constant 10

    #Generate the coefficients file
    eqc.write(ofile)
    
     

Constant Coefficients: Runtime Control
---------------------------------------------
While constant coefficients may be specified via the coefficients file, many of these coefficients represent simulation
"control knobs" that the user may wish to modify at run-time.   For instance, the user may want to frequently use a particular profile for viscous diffusion (:math:`f_3`), but would like to vary its amplitude (:math:`c_5`) between simulations without generating a new coefficients file.  Rayleigh provides the opportunity to override 
all constant coefficients, or a subset of them, through the main_input file.

Consider the example below.

::

   &Reference_Namelist
    ...
    custom_reference_file='mycoeffs.dat'
    override_constants=T
    ra_constants( 2) = 1.0
    ra_constants( 5) = 10.0
    ra_constants(10) = 14.0
    ...
   /

In this example, values of constant coefficients :math:`c_2,\,c_5,\,c_{10}` will be determined entirely via the main_input file and assigned the values
of 1.0, 10.0, and 14.0 respectively.  Values specified in mycoeffs.dat will be ignored completely.   

This behavior is dictated by the override_constants flag, which instructs Rayleigh to ignore ALL constant coefficients specified in the coefficients files.   If a coefficient is not specified in
main_input, its value will be set to Rayleigh's internal default value of 0.   Consider the following example

::

   &Reference_Namelist
    ...
    custom_reference_file='mycoeffs.dat'
    override_constants=T
    ra_constants( 2) = 1.0
    ra_constants(10) = 14.0
    ...
   /

The resulting values of :math:`c_2,\,c_5,\,c_{10}` will be 1.0, 0.0, and 14.0 respectively.  The constant :math:`c_5` will not be set to 20.0 (the value specified in the coefficients file).


To specify a subset of constants, use the override_constant flag for each constant you wish to override, as shown below.

::

   &Reference_Namelist
    ...
    custom_reference_file='mycoeffs.dat'
    override_constant( 2) = T
    override_constant(10) = T
    ra_constants( 2) = 1.0
    ra_constants(10) = 14.0
    ...
   /

In this case, the values of constants :math:`c_2` and :math:`c_{10}` will be taken the main_input file.  The value of :math:`c_5` will be taken from the coefficients file.   If a constant's override flag is set, but its value is not specified in main_input, the default value of zero will be used. 

Augmenting a Rayleigh-Provided Reference State
---------------------------------------------------------
When augumenting one of Rayleigh's internal reference-state types, set the with_custom_reference, with_custom_constants, and with_custom_functions flag in the Reference Namelist.  As an example, to modify the heating and buoyancy profiles using entirely information provided through the equation coefficients file, main_input would contain the following
::

   &Reference_Namelist
    ...
    reference_type=1
    custom_reference_file='mycoeffs.dat'
    with_custom_reference=T
    with_custom_constants=2,10
    with_custom_functions=2,6
    ...
   /

These flags can be used in tandem with the override flags to specify values via main_input.  For example, the following input combination would set a value of :math:`c_2` of 13.0
::

   &Reference_Namelist
    ...
    reference_type=1
    custom_reference_file='mycoeffs.dat'
    with_custom_reference=T
    with_custom_constants=2,10
    with_custom_functions=2,6
    override_constant(2)=T
    ra_constants(2) = 13.0
    ...
   /



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
