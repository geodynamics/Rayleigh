# Changelog
All notable changes **following the 1.0.0 release** of *Rayleigh* will be documented in this file.  Changes in 1.0.0 and prior versions are summarized.  The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
with the addition of author(s), date of change and optionally the relevant issue. 

Add new entries a the bottom of the current list in the subheading. Item format: 
- Description. [Name; date; relevant github issue tag(s) and or pull requests]

This project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Added
- None yet

### Changed
- The rotation rate is now accounted for when computing the maximum allowable timestep.  The timestep may be no larger than c1/4*CFLMax, CFLMax is the CFL safety factor, and where c1 is the 1st Rayleigh constant and the inverse of the rotational timescale (2 Omega for reference_type=2 and 2/Ek for reference_type=1). \[Nick Featherstone; 12-11-2021; [#348](https://github.com/geodynamics/Rayleigh/pull/348)\]

### Fixed
- None yet

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


