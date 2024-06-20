# Changelog
All notable changes **following the 1.0.0 release** of *Rayleigh* will be documented in this file.  Changes in 1.0.0 and prior versions are summarized.  The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
with the addition of author(s), date of change and optionally the relevant issue. 

Add new entries a the bottom of the current list in the subheading. Item format: 
- Description. [Name; date; relevant github issue tag(s) and or pull requests]

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added

- The documentation for solving for active and passive scalar fields has been expanded in the "Under Development" section of the User Guide. \[Cian Wilson; 6-20-2024; [#541](https://github.com/geodynamics/Rayleigh/pull/541)\]

- A docker image for Stampede3 and Frontera (TACC) has been added to the repository.   \[Rene Gassmoeller; 6-20-2024; [#402](https://github.com/geodynamics/Rayleigh/pull/402)\]

- An image gallery has been added to the documentation.   \[Brandon Lazard; 6-20-2024; [#531](https://github.com/geodynamics/Rayleigh/pull/531)\]

- The streamfunction equations solved by Rayleigh have been added to the documentation.   \[Tami Rogers; 6-20-2024; [#533](https://github.com/geodynamics/Rayleigh/pull/533)\]

### Changed

- The logic in the configure script that determines the library-link flags for compilation on different systems has been simplified.   \[Rene Gassmoeller; 6-19-2024; [#534](https://github.com/geodynamics/Rayleigh/pull/534)\]

### Fixed

- Rayleigh no longer attempts to update the record count and close a diagnostics file (e.g., G_Avgs) that failed to open correctly.  \[Nick Featherstone; 6-18-2024; [#510](https://github.com/geodynamics/Rayleigh/pull/510) , [#523](https://github.com/geodynamics/Rayleigh/pull/523)  \]

- Fixed a bug that was causing documentation build to become progressively slower for subsequent builds due to recursive inclusion of the doc\build directory  \[Philipp Edelmann; 6-17-2024; [#527](https://github.com/geodynamics/Rayleigh/pull/527)\]

- Before an existing checkpoint file is overwritten, its file size is first set to zero.  This fixes a bug that caused Rayleigh to generate checkpoint files that were too large in some instances due to resolution changes \[Nick Featherstone; 6-18-2024; [#528](https://github.com/geodynamics/Rayleigh/pull/528)  \]

-  The docker container build script now accounts for the new requirement that container names are all lower case.   \[Rene Gassmoeller; 6-19-2024; [#530](https://github.com/geodynamics/Rayleigh/pull/530)\]

- The standard c++ libraries are no longer linked when compiling Rayleigh; \[Brandon Lazard; 6-19-2024; [#532](https://github.com/geodynamics/Rayleigh/pull/532)\]



## [1.2.0] - 5-29-2024

### Added

- Rayleigh's configure script now supports the MKL package provided by the Debian and Ubuntu package repositories.  This feature is accessed by invoking the -debian-mkl flag when running configure.  \[Philipp Edelmann; 9-14-2022; [#386](https://github.com/geodynamics/Rayleigh/pull/386)\]

- Rayleigh's equation set can now be extended to evolve multiple passive and/or active scalar fields.  \[Cian Wilson, Nick Featherstone, Loren Matilsky, Rafael Fuentes and Maria Camisassa; 11-22-2022; [#408](https://github.com/geodynamics/Rayleigh/pull/408); 12-2-2022 [#415](https://github.com/geodynamics/Rayleigh/pull/408)\; 1-5-2023 [#429](https://github.com/geodynamics/Rayleigh/pull/429) 1-30-2023 [#442](https://github.com/geodynamics/Rayleigh/pull/442)]

- A new Docker image of Rayleigh, based on Ubuntu Jammy, is now generated whenever a PR is merged into master or a new release is created.  \[Rene Gassmoeller and Philipp Edelmann; 9-14-2022 and 9-15-2022; [#385](https://github.com/geodynamics/Rayleigh/pull/385), [#387](https://github.com/geodynamics/Rayleigh/pull/387), [#391](https://github.com/geodynamics/Rayleigh/pull/391), [#392](https://github.com/geodynamics/Rayleigh/pull/392) ]

- Docker images and references based on Ubuntu Bionic (18.04) have been removed.   \[Rene Gassmoeller; 6-14-2023; [#463](https://github.com/geodynamics/Rayleigh/pull/463)\]

- Custom profiles for the background dSdr and dTdr can now be specified when the *with_custom_reference* flag is set to True.  \[Nick Featherstone; 1-5-2023; [#416](https://github.com/geodynamics/Rayleigh/pull/416)\]

- A new nondimensional polytropic reference state is now accessible by specifying reference_type=5. \[Loren Matilsky; 6-14-2023; [#450](https://github.com/geodynamics/Rayleigh/pull/450)\]

- The acknowledgement section of the documentation was updated to reflect the new CIG grant number (NSF-2149126) \[Lorraine Hwang; 5-1-2023; [#454](https://github.com/geodynamics/Rayleigh/pull/454)\]

- A new benchmark, based on the results of Breuer et al. 2010, which incorporates multiple contributions to the buoyancy from different chemical species, is now supported in Rayleigh.  It can be run using input files b2010_case0_0T_input, b2010_case0_80T_input and b2010_case0_100T input. \[Cian Wilson; 6-15-2023; [#467](https://github.com/geodynamics/Rayleigh/pull/467)\]

- A new Jupyter notebook that demonstrates how to perform spectral transforms and related operations on Rayeleigh output has been added to the code base.  It can be found in post_processing/spectral_utils_tutorial.ipynb.  \[Catherine Blume; 6-16-2023; [#470](https://github.com/geodynamics/Rayleigh/pull/470)\]

- Two new Jupyter notebooks demonstrating how ultraspherical polynomials might be used to build a sparse linear system for the implicit solve have been placed in the examples/ultraspherical_polynomials directory. \[Kyle Augustson; 6-20-2023; [#476](https://github.com/geodynamics/Rayleigh/pull/476)\]

- An optional finite-difference scheme in radius is now available.  Updated documentation is provided in the 'Setting Up A Model' section of the User Guide.  An example input file may be found in Rayleigh/input_examples/main_input_mhd_jones_FD. \[Rathish Ratnasingam, Philipp Edelmann, Nick Featherstone; 6-29-2023; [474](https://github.com/geodynamics/Rayleigh/pull/474)\]

- Two new initial conditions have been added to Rayleigh that allow existing thermal and magnetic fields to be augmented at restart with user-supplied profiles for entropy, and the toroidal and poloidal flux functions. \[Bhishek Manek; 12-5-2023; [488](https://github.com/geodynamics/Rayleigh/pull/488)\]

- A new constant, c_11, has been added to Rayleigh's equation set.  This constant multiplies the background entropy profile so that it may be rescaled at runtime, if desired.  \[Loren Matilsky; 12-5-2023; [489](https://github.com/geodynamics/Rayleigh/pull/489)\]

- Rayleigh now supports the MPI_f08 module (in addition to the older MPI module).  To use MPI_f08, call Rayleigh's configure script with the --mpi-f08 flag.  \[Philipp Edelmann; 2-29-2024; [497](https://github.com/geodynamics/Rayleigh/pull/497)\]
  
### Changed
- plot_Shell_Slices.py now generates a Mollweide plot. \[Nick Nelson; 9-12-2022; [#372](https://github.com/geodynamics/Rayleigh/pull/372)\]
  
- Several substantial changes were made to the Documentation in order to improve navigation and readability.  Primary contributors to this effort were Nick Nelson, Rene Gassmoeller, Nick Featherstone.  Changes were applied via a series of pull requests during Fall 2022:  [#378](https://github.com/geodynamics/Rayleigh/pull/378), [#379](https://github.com/geodynamics/Rayleigh/pull/379), [#380](https://github.com/geodynamics/Rayleigh/pull/380), [#381](https://github.com/geodynamics/Rayleigh/pull/381), [#388](https://github.com/geodynamics/Rayleigh/pull/388), [#389](https://github.com/geodynamics/Rayleigh/pull/389), [#390](https://github.com/geodynamics/Rayleigh/pull/390)\, [#394](https://github.com/geodynamics/Rayleigh/pull/394), [#395](https://github.com/geodynamics/Rayleigh/pull/395), [#396](https://github.com/geodynamics/Rayleigh/pull/396), [#397](https://github.com/geodynamics/Rayleigh/pull/397), [#398](https://github.com/geodynamics/Rayleigh/pull/398), [#399](https://github.com/geodynamics/Rayleigh/pull/399), [#400](https://github.com/geodynamics/Rayleigh/pull/400), [#401](https://github.com/geodynamics/Rayleigh/pull/401), [#404](https://github.com/geodynamics/Rayleigh/pull/404), [#405](https://github.com/geodynamics/Rayleigh/pull/405), [#406](https://github.com/geodynamics/Rayleigh/pull/406), [#407](https://github.com/geodynamics/Rayleigh/pull/407), [#411](https://github.com/geodynamics/Rayleigh/pull/411), [#412](https://github.com/geodynamics/Rayleigh/pull/412), [#414](https://github.com/geodynamics/Rayleigh/pull/414), [#417](https://github.com/geodynamics/Rayleigh/pull/417), [#418](https://github.com/geodynamics/Rayleigh/pull/418)
  
- The Readthedocs build environment is now managed using Mamba instead of Conda \[Nick Featherstone; 11-14-2022; [#413](https://github.com/geodynamics/Rayleigh/pull/413)\]

- The AZ_Avgs and Shell_Avgs classes in rayleigh_diagnostics.py now support time-averaging and concatenation of multiple files.  The time-averaged or concatenated data can be output to a new AZ_Avgs or Shell_Avgs file with the same format as that used by Rayleigh.  Run help(AZ_Avgs) or help(Shell_Avgs) to see examples of calling syntax.  \[Nick Featherstone; 5-7-2024; [#499](https://github.com/geodynamics/Rayleigh/pull/499), [#500](https://github.com/geodynamics/Rayleigh/pull/500)\]

- 
### Fixed
- plot_G_Avgs.py now correctly handles multiple files containing differing numbers of timesteps \[Philipp Edelmann; 9-13-2022; [#383](https://github.com/geodynamics/Rayleigh/pull/383)\]

- Several bugs were fixed that prevented Rayleigh constants and functions from being properly reflected in the equation_coefficients file in some instances. \[Loren Matilsky; 1-30-2023, [#443](https://github.com/geodynamics/Rayleigh/pull/443); 3-29-2023, [#447](https://github.com/geodynamics/Rayleigh/pull/447);\]

- The documentation now properly reflects the fact that the constant c_5 is appears in front of the viscous heating term. \[Loren Matilsky; 1-31-2023; [#444](https://github.com/geodynamics/Rayleigh/pull/444)\]

- Enthalpy flux is now output correctly in nondimensional models.  \[Loren Matilsky; 6-27-2023; [#480](https://github.com/geodynamics/Rayleigh/pull/480)\]

- A bug preventing Rayleigh from respecting the users choice of heating function in some instances has been corrected.  \[Loren Matilsky; 4-3-2023; [#492](https://github.com/geodynamics/Rayleigh/pull/492)\]

- The equation_coefficients class in rayleigh_diagnostics.py is now properly version-aware when reading in new or old equation_coefficients files \[Loren Matilsky; 5-7-2024, [#501](https://github.com/geodynamics/Rayleigh/pull/443); 3-29-2023, [#447](https://github.com/geodynamics/Rayleigh/pull/501);\]

## [1.1.0] - 4-29-2022
### Added
- Added new main_input class to rayleigh_diagnostics.py.   Usage examples are provided in Rayleigh/examples/main_input_class. \[Nick Featherstone; 1-11-2022; [#352](https://github.com/geodynamics/Rayleigh/pull/352)\]

- New namelist option m_balance_contiguous = True/False in the numerical_controls_namelist allows using faster transpose routines for 2a3a and 3b2b directions. When True, this option stores the high/low m values in a contiguous fashion, suitable for better vector operations within loops. The default (False) is to use the original high-low pairing. \[Ryan Orvedahl; 4-9-2022; [#363](https://github.com/geodynamics/Rayleigh/pull/363)\] 

- Add some unit tests and have them run within github actions. \[Ryan Orvedahl; 3-31-2022; [#361](https://github.com/geodynamics/Rayleigh/pull/361)\]

### Changed
- The rotation rate is now accounted for when computing the maximum allowable timestep.  Now, the timestep may not exceed CFLMax/(c1*4), where CFLMax is the CFL safety factor.  c1 is the 1st Rayleigh constant, effectively the inverse of the rotational timescale (c1 = 2 Omega for reference_type=2 and 2/Ek for reference_type=1). \[Nick Featherstone; 12-11-2021; [#348](https://github.com/geodynamics/Rayleigh/pull/348)\]

- When running with both internal heating and fix_dTdr_top=.true., the value of dTdr_top is now set internally by Rayleigh.  This is done to ensure consistency with the internal heating as well as any flux passing through the lower boundary (if fixed flux lower boundary conditions are applied).   To override this behavior, set adjust_dTdr to .false. in the Boundary_Conditions namelist. \[Nick Featherstone; 12-13-2021; [#349](https://github.com/geodynamics/Rayleigh/pull/349)\]
### Fixed
- Strict_L_conservation is now set to False when no_slip_boundaries or no_slip_top are true. \[Nick Featherstone; 4-29-2022; [#364](https://github.com/geodynamics/Rayleigh/pull/364)\]

- Fix a bug that incorrectly computed l_max when n_theta was provided. L_max was off by at most one. \[Ryan Orvedahl; 3-28-2022; [#360](https://github.com/geodynamics/Rayleigh/pull/360)\]

## [1.0.1] - 12-11-2021
### Fixed
- Radially varying magnetic diffusivity is now correctly evolved using the Crank-Nicholson scheme for the toroidal magnetic scalar potential *A*. \[Nick Featherstone; 12-03-2021; [#346](https://github.com/geodynamics/Rayleigh/pull/346)\]

## [1.0.0] - 11-12-2021
### Added
- Online, up-to-date-documentation is now available [here](https://rayleigh-documentation.readthedocs.io/en/latest/index.html).  The source files for this documentation are provided as part of the repository and are stored in the ‘doc’ directory.
- Customizability enhancements:  Users may specify custom boundary conditions and/or initial conditions as described [here](https://rayleigh-documentation.readthedocs.io/en/latest/doc/source/User_Guide/physics.html?highlight=generic#generic-boundary-conditions) and [here](https://rayleigh-documentation.readthedocs.io/en/latest/doc/source/User_Guide/physics.html?highlight=generic#generic-initial-conditions).  Users are also free to fully specify the set of constant and nonconstant coefficients that define the system of PDEs solved by Rayleigh, as described in the [documentation](https://rayleigh-documentation.readthedocs.io/en/latest/doc/source/User_Guide/custom_reference_states.html).

### Changed
- Portable checkpoint format:  The checkpoint format has changed such that each checkpoint resides in its own directory.  All files and information needed to restart a model, including a copy of main_input, are now stored within each numbered checkpoint directory.

## [0.9.1] - 04-28-2018
### Added
- The configure script now accepts the --with-custom and -devel flags (run ./configure --help for details)
 
### Fixed
- Several output-quantity codes were double-booked. The output quantity tables have been adjust accordingly.  This error produced issues with quantity codes in the following ranges:

    - 1900--2000 (kinetic energy equation codes)
    - 1400--1600 (thermal energy equation codes)
    - 1600--1800 (magnetic energy equation codes)


## [0.9.0] - 01-23-2018 - Initial Release


