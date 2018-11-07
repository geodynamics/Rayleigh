#######################################
# Set up how to compile the source files
########################################

# Add the source files
SET(src
    # First all Fortran files
    Data_Structures/Structures.F90
    Diagnostics/Diagnostics_ADotGradB.F90
    Diagnostics/Diagnostics_Angular_Momentum.F90
    Diagnostics/Diagnostics_Base.F90
    Diagnostics/Diagnostics_Current_Density.F90
    Diagnostics/Diagnostics_Custom.F90
    Diagnostics/Diagnostics_Energies.F90
    Diagnostics/Diagnostics_Induction.F90
    Diagnostics/Diagnostics_Inertial_Forces.F90
    Diagnostics/Diagnostics_Interface.F90
    Diagnostics/Diagnostics_KE_Flux.F90
    Diagnostics/Diagnostics_Linear_Forces.F90
    Diagnostics/Diagnostics_Lorentz_Forces.F90
    Diagnostics/Diagnostics_Magnetic_Field.F90
    Diagnostics/Diagnostics_Mean_Correction.F90
    Diagnostics/Diagnostics_Miscellaneous.F90
    Diagnostics/Diagnostics_Poynting_Flux.F90
    Diagnostics/Diagnostics_Second_Derivatives.F90
    Diagnostics/Diagnostics_Thermal_Energies.F90
    Diagnostics/Diagnostics_Thermal_Equation.F90
    Diagnostics/Diagnostics_Thermodynamic_Gradients.F90
    Diagnostics/Diagnostics_TurbKE_Budget.F90
    Diagnostics/Diagnostics_Velocity_Diffusion.F90
    Diagnostics/Diagnostics_Velocity_Field.F90
    Diagnostics/Diagnostics_Vorticity_Field.F90
    IO/Spherical_IO.F90
    Math_Layer/Chebyshev_Polynomials_Alt.F90
    Math_Layer/Chebyshev_Polynomials.F90
    Math_Layer/Finite_Difference.F90
    Math_Layer/Fourier_Derivatives.F90
    Math_Layer/Fourier_Transform.F90
    Math_Layer/Legendre_Polynomials.F90
    Math_Layer/Legendre_Transforms.F90
    Math_Layer/Linear_Solve.F90
    Math_Layer/Math_Constants.F90
    Math_Layer/Math_Utility.F90
    Math_Layer/Spectral_Derivatives.F90
    Math_Layer/Theta_Derivatives.F90
    Math_Layer/Timing.F90
    Parallel_Framework/All_to_All.F90
    Parallel_Framework/General_MPI.F90
    Parallel_Framework/ISendReceive.F90
    Parallel_Framework/Load_Balance.F90
    Parallel_Framework/MPI_Base.F90
    Parallel_Framework/MPI_LAYER.F90
    Parallel_Framework/Parallel_Framework.F90
    Parallel_Framework/SendReceive.F90
    Parallel_Framework/Spherical_Buffer.F90
    Physics/Benchmarking.F90
    Physics/BoundaryConditions.F90
    Physics/Checkpointing.F90
    Physics/ClockInfo.F90
    Physics/Controls.F90
    Physics/Fields.F90
    Physics/Initial_Conditions.F90
    Physics/Input.F90
    Physics/Main.F90
    Physics/ProblemSize.F90
    Physics/ReferenceState.F90
    Physics/Run_Parameters.F90
    Physics/Sphere_Driver.F90
    Physics/Sphere_Hybrid_Space.F90
    Physics/Sphere_Linear_Terms.F90
    Physics/Sphere_Physical_Space.F90
    Physics/Sphere_Spectral_Space.F90
    Physics/Timers.F90
    Physics/TransportCoefficients.F90
    Plugins/STABLE_plugin.F90
    Test_Suite/Test_Cheby.F90
    Test_Suite/Test_SHT.F90
    Test_Suite/TestSuite.F90
    Utility/BufferedOutput.F90
    Utility/General_IO.F90
    Utility/MakeDir.F90

    # Then also the C files
    Utility/cmkdir.c
    )

# Make sure both C and Fortran compilers know where to find include files
INCLUDE_DIRECTORIES(Include)

# Set compiler flags for the Fortran files
SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-line-length-140")
SET(CMAKE_Fortran_FLAGS_DEBUG "-fbounds-check -fbacktrace")



# Define the executable in terms of the source files, and ensure that
# it links against the correct libraries.
ADD_EXECUTABLE(${RayleighEXE} ${src})
TARGET_LINK_LIBRARIES(${RayleighEXE} blas lapack fftw3)
