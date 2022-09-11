
.. _namelists:

Main_Input Namelists
==============================================

This page provides a quick reference for all support main_input namelist variables.


Problemsize
-------------------------
This namelist is used to specify the grid.


.. include:: ../Namelist_Definitions/problemsize_namelist.txt 

Numerical Controls
-------------------------
This namelist provides access to Rayleigh's run-time optimization options.

.. include:: ../Namelist_Definitions/numerical_controls_namelist.txt 

Physical Controls
-------------------------
This namelist controls the physical effects used in a Rayleigh simulation.

.. include:: ../Namelist_Definitions/physical_controls_namelist.txt 

Temporal Controls
-------------------------
This namelist controls timing, time-stepping, and checkpointing in Rayleigh.

.. include:: ../Namelist_Definitions/temporal_controls_namelist.txt 

IO Controls
-------------------------
This namelist provides various options to control Rayleigh's input and output cadence and structure.

.. include:: ../Namelist_Definitions/io_controls_namelist.txt 

Output
-------------------------
This namelist is described in extensive detail in Rayleigh/post_processing/Diagnostic_Plotting.ipynb.  Please see that document for a discussion of these namelist variables and the general structure of Rayleigh's output.


Boundary Conditions
-------------------------
This namelist provides those options necessary to determine the boundary conditions employed in a Rayleigh model.

.. include:: ../Namelist_Definitions/boundary_conditions_namelist.txt 

Initial Conditions
-------------------------
All variables necessary to initialize velocity, temperature, pressure, and magnetic field are supplied here.

.. include:: ../Namelist_Definitions/initial_conditions_namelist.txt 


Reference
-------------------------
This namelist provides options to control the properties of Rayleigh's background state.

.. include:: ../Namelist_Definitions/reference_namelist.txt 

Transport
-------------------------
This namelist enables control of Rayleigh's diffusivities.

.. include:: ../Namelist_Definitions/transport_namelist.txt 
