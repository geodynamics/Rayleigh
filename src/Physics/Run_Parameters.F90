!
!  Copyright (C) 2018 by the authors of the RAYLEIGH code.
!
!  This file is part of RAYLEIGH.
!
!  RAYLEIGH is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  RAYLEIGH is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!

Module Run_Parameters

    Implicit None

    ! This included file is generated during the configure/make process
    ! It contains machine-specific information regarding
    ! the compiler, git repo, etc. that Rayleigh outputs at runtime.
    include "Run_Param_Header.F"

Contains

    Subroutine Write_Run_Parameters()

        Use Controls, Only : my_path, magnetism, jobinfo_file
        Use Controls, Only : numerical_controls_namelist, physical_controls_namelist, &
                           temporal_controls_namelist, io_controls_namelist
        Use Parallel_Framework, Only: pfi
        Use ProblemSize, Only : my_rank, problemsize_namelist, nprow, npcol, &
                              n_r, n_theta, l_max, rmin, rmax, &
                              ndomains, ncheby, dealias_by, domain_bounds
        Use Spherical_IO, Only : output_namelist
        Use BoundaryConditions, Only : boundary_conditions_namelist
        Use Initial_Conditions, Only : initial_conditions_namelist
        Use TestSuite, Only : test_namelist
        Use ReferenceState, Only : reference_namelist, reference_type, Rayleigh_Number, &
                                 Ekman_Number, Prandtl_Number, Magnetic_Prandtl_Number, &
                                 Modified_Rayleigh_Number, poly_n, poly_nrho, poly_mass, &
                                 Angular_Velocity, poly_rho_i, pressure_specific_heat
        Use TransportCoefficients, Only : transport_namelist

        Implicit None

        Character(len=16 ) :: date, time, zone_in
        Integer :: values(8)
        Integer :: io, ierr, i
        Logical :: file_exist
        Integer :: ntrim = 4

        Call Date_and_Time(date, time, zone_in, values) ! get current date/time

        io = 15
        999  Format (75('='))
        1001 Format (a,a)
        1002 Format (a,i6)
        1003 Format (a,i4.4,'-',i2.2,'-',i2.2)
        1004 Format (a,i2.2,':',i2.2,':',i2.2)
        1005 Format (a,ES11.4)

        If (my_rank .eq. 0) Then
            Inquire(file=Trim(my_path)//Trim(jobinfo_file), exist=file_exist)
            Open(unit=io, file=Trim(my_path)//Trim(jobinfo_file), form='formatted', &
               action='write', access='sequential', status='unknown', iostat=ierr)
            If (ierr .eq. 0) Then
                Write(io,*)
                Write(io,999)
                Write(io,*) "Job Information"
                Write(io,999)
                Write(io,1003) "Date = ", values(1), values(2), values(3)
                Write(io,1004) "Time = ", values(5), values(6), values(7)
                Write(io,1001) "Path = ", Trim(my_path)
                Write(io,*)
                Write(io,999)
                Write(io,*) "Compiler Information"
                Write(io,999)

                ! the string slicing is to remove the beginning "---short_name---"
                ! as well as the trailing "---Char(0)"
                Write(io,1001) "Compiler Location = ", &
                             Trim(FC_location(18:Len_Trim(FC_location)-ntrim))
                Write(io,1001) "Compiler Version  = ", &
                             Trim(FC_version(17:Len_Trim(FC_version)-ntrim))
                Write(io,*)
                Write(io,1001) "Compiler Flags    = ", &
                             Trim(build_fflags(19:Len_Trim(build_fflags)-ntrim))
                Write(io,1001) "Compiler Library  = ", &
                             Trim(build_lib(16:Len_Trim(build_lib)-ntrim))
                Write(io,*)
                Write(io,999)
                Write(io,*) "Build Information"
                Write(io,999)
                Write(io,1001) "Build Date       = ", &
                             Trim(build_date(17:Len_trim(build_date)-ntrim))
                Write(io,1001) "Build Machine    = ", &
                             Trim(build_machine(20:Len_trim(build_machine)-ntrim))
                Write(io,1001) "Build Directory  = ", &
                             Trim(build_dir(16:Len_trim(build_dir)-ntrim))
                Write(io,1001) "Custom Directory = ", &
                             Trim(build_custom_dir(17:Len_trim(build_custom_dir)-ntrim))
                Write(io,*)
                Write(io,1001) "Build Version = ", &
                             Trim(build_version(20:Len_trim(build_version)-ntrim))
                Write(io,1001) "Git Hash = ", &
                             Trim(build_git_commit(15:Len_trim(build_git_commit)-ntrim))
                Write(io,1001) "Git URL = ", &
                             Trim(build_git_url(14:Len_trim(build_git_url)-ntrim))
                Write(io,1001) "Git Branch = ", &
                             Trim(build_git_branch(17:Len_trim(build_git_branch)-ntrim))
                Write(io,*)

                Write(io,999)
                Write(io,*) "Preamble Information"
                Write(io,999)
                Write(io,*) "MPI"
                Write(io,1002) "    NCPU  : ", pfi%wcomm%np
                Write(io,1002) "    NPROW : ", nprow
                Write(io,1002) "    NPCOL : ", npcol
                Write(io,*)
                Write(io,*) "Grid"
                Write(io,1002) "    N_R                : ", n_r
                Write(io,1002) "    N_THETA            : ", n_theta
                Write(io,1002) "    Ell_MAX            : ", l_max
                Write(io,1005) "    R_MIN              : ", rmin
                Write(io,1005) "    R_MAX              : ", rmax
                Write(io,1002) "    Chebyshev Domains  : ", ndomains
                Write(io,*)
                Do i=1,ndomains
                    Write(io,1002) "      Domain ", i
                    Write(io,1002) "        Grid points           : ", ncheby(i)
                    If (dealias_by(i) .eq. -1) Then
                        Write(io,1002) "        Dealiased Polynomials : ", (ncheby(i)*2)/3
                    Else
                        Write(io,1002) "        Dealiased Polynomials : ", ncheby(i)-dealias_by(i)
                    Endif
                    Write(io,1005) "        Domain Lower Bound    : ", domain_bounds(i)
                    Write(io,1005) "        Domain Upper Bound    : ", domain_bounds(i+1)
                Enddo
                Write(io,*)
                Write(io,*) "Reference State"
                If (reference_type .eq. 1) Then
                    Write(io,1001) "    Reference type          : ", "Boussinesq (Non-dimensional)"
                    Write(io,1005) "    Rayleigh Number         : ", Rayleigh_Number
                    Write(io,1005) "    Ekman Number            : ", Ekman_Number
                    Write(io,1005) "    Prandtl Number          : ", Prandtl_Number
                    If (magnetism) Then
                        Write(io,1005) "    Magnetic Prandtl Number : ", Magnetic_Prandtl_Number
                    Endif
                Else If (reference_type .eq. 2) Then
                    Write(io,1001) "    Reference type           : ", "Polytrope (Non-dimensional)"
                    Write(io,1005) "    Modified Rayleigh Number : ", Rayleigh_Number
                    Write(io,1005) "    Ekman Number             : ", Ekman_Number
                    Write(io,1005) "    Prandtl Number           : ", Prandtl_Number
                    If (magnetism) Then
                        Write(io,1005) "    Magnetic Prandtl Number  : ", Magnetic_Prandtl_Number
                    Endif
                    Write(io,1005) "    Polytropic Index         : ", poly_n
                    Write(io,1005) "    Density Scaleheights     : ", poly_nrho
                Else If (reference_type .eq. 3) Then
                    Write(io,1001) "    Reference type                : ", "Polytrope (Dimensional)"
                    Write(io,1005) "    Angular Velocity (rad/s)      : ", Angular_Velocity
                    Write(io,1005) "    Inner-Radius Density (g/cm^3) : ", poly_rho_i
                    Write(io,1005) "    Interior Mass (g)             : ", poly_mass
                    Write(io,1005) "    Polytropic Index              : ", poly_n
                    Write(io,1005) "    Density Scaleheights          : ", poly_nrho
                    Write(io,1005) "    CP (erg g^-1 cm^-3 K^-1)      : ", pressure_specific_heat
                Else If (reference_type .eq. 4) Then
                    Write(io,1001) "    Reference type                : ", "Custom"
                Endif
                Write(io,*)

                Write(io,999)
                Write(io,*) "Namelist Information"
                Write(io,999)
                Write(io,nml=problemsize_namelist)
                Write(io,nml=numerical_controls_namelist)
                Write(io,nml=physical_controls_namelist)
                Write(io,nml=temporal_controls_namelist)
                Write(io,nml=io_controls_namelist)
                Write(io,nml=output_namelist)
                Write(io,nml=boundary_conditions_namelist)
                Write(io,nml=initial_conditions_namelist)
                Write(io,nml=test_namelist)
                Write(io,nml=reference_namelist)
                Write(io,nml=transport_namelist)

            Endif

            Close(unit=io)

        Endif

    End Subroutine Write_Run_Parameters

End Module Run_Parameters
