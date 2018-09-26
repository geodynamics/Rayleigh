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

!This module allows benchmarking measurements to be carried out in-situ
!Benchmark 1: Christensen et al. Model X (hydro)
!Benchmark 2: Christensen et al. Model Y (mhd)
!Benchmark 3: Jones et al. hydro + steady
!Benchmark 4: Jones et al. mhd + steady
!Benchmark 5: Jones et al. mhd + unsteady (In Development)

!  This module is intended for checking ACCURACY of Rayleigh's results.
!  This module has NOT been programmed with EFFICIENCY in mind (much).
Module Benchmarking
    Use ProblemSize
    Use Controls
    Use Spherical_IO
    Use Fields
    Use Legendre_Polynomials, Only : gl_weights
    Use ReferenceState
    Use TransportCoefficients
    Use Math_Constants
    Use BoundaryConditions
    Use Initial_Conditions
    Use TransportCoefficients
    Implicit None

    Integer, Private :: nobs, msymm
    Integer, Private :: max_numt, numt_ind, global_count, num_int
    Integer, Private :: report_interval = 90000000
    Integer, Private :: integration_interval =90000000
    Real*8 :: mag_factor
    Real*8, Allocatable :: time_series(:,:), time_saves(:), iter_saves(:), obs_series(:,:)
    Integer :: drift_sign, num_rep
    Integer :: r_one, r_two, theta_one, theta_two
    Integer :: strip_owners(1:4), btags(1:4), tvals(1:4), rvals(1:4)
    Integer :: num_strips = 4  !The number of strips we intend to average over
    Logical :: have_strip(1:4) , have_reference_strips = .false.
    Real*8, Allocatable :: strips(:,:) ! equatorial-mid-shell strips used for tracking pattern drift rate
    Real*8, Allocatable :: xnow(:), xlast(:), xref(:), drifts(:,:), report_sdev(:)
    Real*8, Allocatable :: observations(:), report_vals(:), suggested_vals(:)
    Real*8 :: volume_norm = 1.0d0
    Real*8 :: drift
    Real*8 :: drift_reference_time, previous_time

    Character*80, Allocatable :: report_names(:)
    Character*120 :: benchmark_name
Contains

    Subroutine Benchmark_Input_Reset()
        Implicit None
        Integer :: mode_remember, integration_remember, report_remember
        Integer :: init_remember, restart_remember, minit_remember
        ! This routine re-initializes input values to their benchmark values

        If (benchmark_mode .gt. 0) Then

            mode_remember = benchmark_mode  ! Keep track of a few things before restoring defaults
            integration_remember = benchmark_integration_interval
            report_remember = benchmark_report_interval
            init_remember = init_type
            restart_remember = restart_iter
            minit_remember = magnetic_init_type

            Call Restore_Transport_Defaults()
            Call Restore_InitialCondition_Defaults()
            Call Restore_BoundaryCondition_Defaults()
            Call Restore_Reference_Defaults()
            !Call Restore_Temporal_Defaults()
            Call Restore_Physics_Defaults()

            benchmark_mode = mode_remember
            benchmark_integration_interval = integration_remember
            benchmark_report_interval = report_remember
        Endif

        If (benchmark_mode .eq. 1) Then
            !Christsensen et al. Hydro (case 0)

            ! Domain Size
            shell_depth = 1.0d0
            aspect_ratio = 0.35d0

            !Temporal Controls
            rotation = .true.
            viscous_heating = .false.
            max_time_step = 1.0d-4
            alpha_implicit = 0.50001d0
            cflmin = 0.4d0
            cflmax = 0.6d0

            !Boundary Conditions
            no_slip_boundaries = .true.
            strict_L_Conservation = .false.
            dtdr_bottom = 0.0d0
            T_Top    = 0.0d0
            T_Bottom = 1.0d0
            fix_tvar_top = .true.
            fix_tvar_bottom = .true.
            fix_dtdr_bottom = .false.

            !Initial Conditions
            init_type = 1
            If (init_remember .eq. -1) Then
                 ! Allow for restarts
                 init_type = -1
                 restart_iter = restart_remember
            Endif

            !Reference_Namelist
            Ekman_Number = 1.0d-3
            Rayleigh_Number = 1.0d5
            Prandtl_Number = 1.0d0
            reference_type = 1
            heating_type = 0
            gravity_power = 1.0d0
            !dimensional_reference = .false.


        Endif

        If (benchmark_mode .eq. 2) Then
            ! Christensen et al. MHD (case 1)
            ! Domain Size
            shell_depth = 1.0d0
            aspect_ratio = 0.35d0

            !Physical Controls
            rotation = .true.
            magnetism = .true.
            viscous_heating = .false.
            ohmic_heating = .false.

            !Temporal Controls
            max_time_step = 1.0d-4
            alpha_implicit = 0.50001d0
            cflmin = 0.4d0
            cflmax = 0.6d0

            !Boundary Conditions
            no_slip_boundaries = .true.
            strict_L_Conservation = .false.
            dtdr_bottom = 0.0d0
            T_Top    = 0.0d0
            T_Bottom = 1.0d0
            fix_tvar_top = .true.
            fix_tvar_bottom = .true.
            fix_dtdr_bottom = .false.

            !Initial Conditions
            init_type = 1
            magnetic_init_type = 1
            If ( (init_remember .eq. -1) .or. (minit_remember .eq. -1) ) Then
                 ! Allow for restarts (assume hydro and mhd are both restarted)
                 init_type = -1
                 magnetic_init_type = -1
                 restart_iter = restart_remember
            Endif

            !Reference_Namelist
            Ekman_Number = 1.0d-3
            Rayleigh_Number = 1.0d5
            Prandtl_Number = 1.0d0
            Magnetic_Prandtl_Number = 5.0d0
            reference_type = 1
            heating_type = 0
            gravity_power = 1.0d0
            !dimensional_reference = .false.



        Endif

        If (benchmark_mode .eq. 3) Then
            ! Jones et al. Hydro
            ! Domain Size
            rmin = 2.45d9
            rmax = 7.0d9


            !Physical Controls
            rotation = .true.


            !Temporal Controls
            max_time_step = 30.0d0
            alpha_implicit = 0.50001d0
            cflmin = 0.4d0
            cflmax = 0.6d0

            !Boundary Conditions
            no_slip_boundaries = .false.
            strict_L_Conservation = .false.
            dtdr_bottom = 0.0d0
            T_Top    = 0.0d0
            T_Bottom = 851225.7d0
            fix_tvar_top = .true.
            fix_tvar_bottom = .true.
            fix_dtdr_bottom = .false.

            !Initial Conditions
            init_type = 6
            If (init_remember .eq. -1) Then
                 ! Allow for restarts
                 init_type = -1
                 restart_iter = restart_remember
            Endif

            !Reference_Namelist
            reference_type = 2
            heating_type = 0
            luminosity = 3.846d33
            poly_n = 2.0d0
            poly_Nrho = 5.0d0
            poly_mass = 1.9D30
            poly_rho_i = 1.1d0
            pressure_specific_heat = 1.0509d8
            !dimensional = .true.
            angular_velocity = 1.76d-4

            !Transport Namelist
            nu_top    = 3.64364d12
            kappa_top = 3.64364d12


        Endif


        If (benchmark_mode .eq. 4) Then
            ! Jones et al. MHD (steady)
            ! Domain Size
            rmin = 2.45d9
            rmax = 7.0d9


            !Physical Controls
            rotation = .true.
            magnetism = .true.

            !Temporal Controls
            max_time_step = 200.0d0
            alpha_implicit = 0.50001d0
            cflmin = 0.4d0
            cflmax = 0.6d0

            !Boundary Conditions
            no_slip_boundaries = .false.
            strict_L_Conservation = .true.
            dtdr_bottom = 0.0d0
            T_Top    = 0.0d0
            T_Bottom = 774268.3d0
            fix_tvar_top = .true.
            fix_tvar_bottom = .true.
            fix_dtdr_bottom = .false.

            !Initial Conditions
            init_type = 7
            magnetic_init_type = 7
            mag_amp = 1.0d0
            temp_amp = 1.0d1
            temp_w = 0.01d4
            conductive_profile=.true.

            If ( (init_remember .eq. -1) .or. (minit_remember .eq. -1) ) Then
                 ! Allow for restarts (assume hydro and mhd are both restarted)
                 init_type = -1
                 magnetic_init_type = -1
                 restart_iter = restart_remember
            Endif


            !Reference_Namelist
            reference_type = 2
            heating_type = 0
            luminosity = 3.846d33
            poly_n = 2.0d0
            poly_Nrho = 3.0d0
            poly_mass = 1.9D30
            poly_rho_i = 1.1d0
            pressure_specific_heat = 1.0509d8
            !dimensional = .true.
            angular_velocity = 1.76d-4

            !Transport Namelist
            nu_top    = 7.28728d12
            kappa_top = 7.28728d12
            eta_top = 1.457456d11

        Endif




    End Subroutine Benchmark_Input_Reset

    Subroutine Initialize_Benchmarking
        Implicit None

        Integer :: p,i
        Integer :: your_row_rank, your_col_rank
        Integer :: your_r_min, your_r_max
        Integer :: your_theta_min, your_theta_max
        Logical :: have_r_one, have_r_two, have_theta_one, have_theta_two


        global_count = 0
        numt_ind = 1
        If (magnetism) Then
            num_int = 6
            nobs = 4
        Else
            num_int = 3
            nobs = 3
        Endif

        If (benchmark_mode .eq. 1) Then
            benchmark_name = 'Christensen et al. 2001  (Non-MHD, Case 0)'

            drift_sign = 1 ! Prograde Drift
            integration_interval = 100
            report_interval = 10000
            max_numt = report_interval/integration_interval
            mag_factor = 1.0d0/(2*ekman_number*magnetic_prandtl_number)
            msymm = 4
            num_rep = 4
            volume_norm = 1
            Allocate(report_names(1:4))
            Allocate(report_vals(1:4))
            Allocate(report_sdev(1:4))
            Allocate(suggested_vals(1:4))
            report_names(1) = '  Kinetic Energy  : '
            report_names(2) = '  Temperature     : '
            report_names(3) = '  Vphi            : '
            report_names(4) = '  Drift Frequency : '

            !Suggested values from Jones et al. 2000
            suggested_vals(1) = 58.348d0
            suggested_vals(2) = 0.42812d0
            suggested_vals(3) = -10.1571d0
            suggested_vals(4) = 0.1824d0

            ! Ideally, we override namelist values with benchmark values here

        Endif

        If (benchmark_mode .eq. 2) Then
            benchmark_name = 'Christensen et al. 2001  (MHD, Case 1)'

            drift_sign = -1 ! Retrograde Drift
            integration_interval = 100 !100
            report_interval = 10000 ! 10000
            max_numt = report_interval/integration_interval
            mag_factor = 1.0d0/(2*ekman_number*magnetic_prandtl_number)
            msymm = 4
            num_rep = 6
            volume_norm = 1
            Allocate(report_names(1:6))
            Allocate(report_vals(1:6))
            Allocate(report_sdev(1:6))
            Allocate(suggested_vals(1:6))
            report_names(1) = '  Kinetic Energy  : '
            report_names(2) = '  Magnetic Energy : '
            report_names(3) = '  Temperature     : '
            report_names(4) = '  Vphi            : '
            report_names(5) = '  Btheta          : '
            report_names(6) = '  Drift Frequency : '

            !Suggested values from Jones et al. 2000
            suggested_vals(1) = 30.773d0
            suggested_vals(2) = 626.41d0
            suggested_vals(3) = 0.37338d0
            suggested_vals(4) = -7.6250d0
            suggested_vals(5) = -4.9289d0
            suggested_vals(6) = -3.1017d0
            ! Ideally, we override namelist values with benchmark values here

        Endif

        If (benchmark_mode .eq. 3) Then
            benchmark_name = 'Jones et al. 2001  (Hydrodynamic Case)'

            drift_sign = 1 ! Prograde Drift
            integration_interval = 100 !100
            report_interval = 10000 ! 10000
            max_numt = report_interval/integration_interval
            msymm = 19
            num_rep = 6
            volume_norm = (four_pi/3.0d0)*(radius(1)**3-radius(n_r)**3)

            Allocate(report_names(1:6))
            Allocate(report_vals(1:6))
            Allocate(report_sdev(1:6))
            Allocate(suggested_vals(1:6))
            report_names(1) = '  Kinetic Energy  : '
            report_names(2) = '  Zonal KE        : '
            report_names(3) = '  Meridional KE   : '
            report_names(4) = '  Entropy         : '
            report_names(5) = '  Vphi            : '
            report_names(6) = '  Drift Frequency : '

            !Suggested values from Jones et al. 2000
            suggested_vals(1) = 5.57028d35
            suggested_vals(2) = 6.38099d34
            suggested_vals(3) = 1.49825d32
            suggested_vals(4) = 7.9452d5
            suggested_vals(5) = 690.27
            suggested_vals(6) = 3.10512d-6
            ! Ideally, we override namelist values with benchmark values here

        Endif

        If (benchmark_mode .eq. 4) Then
            benchmark_name = 'Jones et al. 2001  (Steady Dynamo Case)'

            drift_sign = 1 ! Prograde Drift
            integration_interval = 100
            report_interval = 10000
            max_numt = report_interval/integration_interval
            msymm = 7
            num_rep = 10
            volume_norm = (four_pi/3.0d0)*(radius(1)**3-radius(n_r)**3)
            mag_factor = over_eight_pi
            Allocate(report_names(1:num_rep))
            Allocate(report_vals(1:num_rep))
            Allocate(report_sdev(1:num_rep))
            Allocate(suggested_vals(1:num_rep))
            report_names(1) = '  Kinetic Energy  : '
            report_names(2) = '  Zonal KE        : '
            report_names(3) = '  Meridional KE   : '
            report_names(4) = '  Magnetic Energy : '
            report_names(5) = '  Zonal ME        : '
            report_names(6) = '  Meridional ME   : '
            report_names(7) = '  Entropy         : '
            report_names(8) = '  Vphi            : '
            report_names(9) = '  Btheta          : '
            report_names(10) = '  Drift Frequency : '


            !Suggested values from Jones et al. 2000
            suggested_vals(1) = 8.03623d36
            suggested_vals(2) = 1.15318d36
            suggested_vals(3) = 1.01587d33
            suggested_vals(4) = 6.13333d36
            suggested_vals(5) = 4.62046d36
            suggested_vals(6) = 3.24927d35
            suggested_vals(7) = 6.0893d5
            suggested_vals(8) = -2942.2d0
            suggested_vals(9) = 272.92d0
            suggested_vals(10) = 4.30760d-6
            ! Ideally, we override namelist values with benchmark values here

        Endif

        If (my_rank .eq. 0) Then
            If (benchmark_mode .gt. 0) Then
                Call stdout%print(" ")
                Call stdout%print(" -- Benchmarking Mode is Activated.")
                Call stdout%print(" -- Selected Benchmark :  "//trim(benchmark_name))
                Call stdout%print(" ")
            Endif
        Endif
        If (benchmark_integration_interval .gt. 0) Then
            If (benchmark_report_interval .gt. 0) Then
                integration_interval = benchmark_integration_interval
                report_interval = benchmark_report_interval
                max_numt = report_interval/integration_interval
            Endif
        Endif

        !NOTE:  for dimensional anelastic runs, use mag_factor = over_eight_pi

        ! We need to bracket the equator and the radial center of the domain

        Do i = 1, 4
            btags(i)=100+i  !MPI Tags for each of the four strips
        Enddo

        !n_theta should always be even
        theta_one = n_theta/2
        theta_two = theta_one+1

        !n_r is always even in Chebyshev runs - assume it's even for now
        r_one = n_r/2
        r_two = r_one+1

        rvals(1:2) = r_one
        rvals(3:4) = r_two
        tvals(1) = theta_one
        tvals(2) = theta_two
        tvals(3) = theta_one
        tvals(4) = theta_two

        !///////////////////////////////////////////////////////////
        ! Figure out if I own one of the desired strips
        have_r_one = .false.
        have_r_two = .false.
        have_theta_two = .false.
        have_theta_one = .false.
        have_strip(1:4) = .false.


        If ((r_one .le. my_r%max) .and. (r_one .ge. my_r%min) ) have_r_one = .true.
        If ((r_two .le. my_r%max) .and. (r_two .ge. my_r%min) ) have_r_two = .true.
        If ((theta_one .le. my_theta%max) .and. (theta_one .ge. my_theta%min) ) have_theta_one = .true.
        If ((theta_two .le. my_theta%max) .and. (theta_two .ge. my_theta%min) ) have_theta_two = .true.


        If (have_r_one) Then
            if (have_theta_one) have_strip(1) = .true.
            if (have_theta_two) have_strip(2) = .true.
        Endif

        If (have_r_two) Then
            if (have_theta_one) have_strip(3) = .true.
            if (have_theta_two) have_strip(4) = .true.
        Endif



        If (my_rank .eq. 0) Then
            Allocate(strips(1:n_phi,1:nobs))
            Allocate(observations(1:nobs))
            Allocate(obs_series(1:max_numt,1:nobs))
            Allocate(time_series(1:max_numt,1:num_int))  ! Possible that only rank 0 needs these -- check
            Allocate(time_saves(1:max_numt))
            Allocate(iter_saves(1:max_numt))
            Allocate(drifts(1:max_numt,1:2))
            Allocate(xnow(1:msymm),xlast(1:msymm),xref(1:msymm))

            Do p = 0, (npcol*nprow)-1
                your_row_rank = mod(p,nprow)
                your_col_rank = p/nprow

                your_r_min = pfi%all_1p(your_col_rank)%min
                your_r_max = pfi%all_1p(your_col_rank)%max
                If ((r_one .le. your_r_max) .and. (r_one .ge. your_r_min) ) Then
                    your_theta_min = pfi%all_2p(your_row_rank)%min
                    your_theta_max = pfi%all_2p(your_row_rank)%max
                    If ((theta_one .le. your_theta_max) .and. &
                        &  (theta_one .ge. your_theta_min) ) Then

                        strip_owners(1) = p

                    Endif
                    If ((theta_two .le. your_theta_max) .and. &
                        & (theta_two .ge. your_theta_min) ) Then

                        strip_owners(2) = p

                    Endif
                Endif
                If ((r_two .le. your_r_max) .and. (r_two .ge. your_r_min) ) Then
                    your_theta_min = pfi%all_2p(your_row_rank)%min
                    your_theta_max = pfi%all_2p(your_row_rank)%max
                    If ((theta_one .le. your_theta_max) .and. &
                        &  (theta_one .ge. your_theta_min) ) Then

                        strip_owners(3) = p

                    Endif
                    If ((theta_two .le. your_theta_max) .and. &
                        & (theta_two .ge. your_theta_min) ) Then

                        strip_owners(4) = p

                    Endif
                Endif


            Enddo
        Endif

    End Subroutine Initialize_Benchmarking

    Subroutine Benchmark_Checkup(buffer,iteration, current_time)
        Implicit None
        Integer, Intent(In) :: iteration
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,:)
        Real*8, Intent(In) :: current_time
        Real*8 :: tmp, tmp2, tmp3, time_passed, over_n_phi
        Real*8 :: rel_diff, mean_value, sdev_value

        Integer :: i,p,t,r, funit, iter_start, iter_end
        Real*8, Allocatable :: ell0_values(:,:), volume_integrals(:), volume_sdev(:)
        Real*8, Allocatable :: qty(:,:,:,:), obs_sdev(:)
        Character*120 :: report_file

        Character*14 :: val_str, sdev_str, rel_str, sug_str, dt_str
        Character*8 :: fmtstr = '(F14.6)'
        Character*8 :: iter_string

        If (benchmark_mode .gt. 0) Then
        If (mod(iteration,integration_interval) .eq. 0) Then
            !Write(6,*)'Integrating!'
            !First we grab the volume-integrated quantities

            !We keep a time-series of volume-integrated quantities in memory
            over_n_phi = 1.0d0/dble(n_phi)


            Allocate(qty(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:num_int))

            qty(:,:,:,:) = 0.0d0


            !KE
            qty(1:n_phi,:,:,1) = buffer(1:n_phi,:,:,vphi)**2
            qty(1:n_phi,:,:,1) = qty(1:n_phi,:,:,1)+buffer(1:n_phi,:,:,vr)**2
            qty(1:n_phi,:,:,1) = qty(1:n_phi,:,:,1)+buffer(1:n_phi,:,:,vtheta)**2
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    Do p = 1, n_phi
                        qty(p,r,t,1) = qty(p,r,t,1)*ref%density(r)*0.5d0
                    Enddo
                Enddo
            Enddo


            !Zonal KE
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    ! compute mean v_phi here
                    tmp = 0.0d0
                    Do p = 1, n_phi
                        tmp = tmp+buffer(p,r,t,vphi)
                    Enddo
                    tmp = tmp*over_n_phi
                    tmp = 0.5d0*ref%density(r)*tmp**2
                    qty(:,r,t,2) = tmp
                Enddo
            Enddo

            !Meridional KE
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    ! compute mean v_phi here
                    tmp = 0.0d0
                    tmp2 = 0.0d0
                    Do p = 1, n_phi
                        tmp = tmp+buffer(p,r,t,vr)
                        tmp2 = tmp2+buffer(p,r,t,vtheta)
                    Enddo
                    tmp = tmp*over_n_phi
                    tmp2 = tmp2*over_n_phi
                    tmp = 0.5d0*ref%density(r)*(tmp**2+tmp2**2)
                    qty(:,r,t,3) = tmp
                Enddo
            Enddo

            If (magnetism) Then

                ! ME
                qty(1:n_phi,:,:,4) = buffer(1:n_phi,:,:,bphi)**2
                qty(1:n_phi,:,:,4) = qty(1:n_phi,:,:,4)+buffer(1:n_phi,:,:,br)**2
                qty(1:n_phi,:,:,4) = (qty(1:n_phi,:,:,4) &
                    & + buffer(1:n_phi,:,:,btheta)**2)


                ! Zonal ME
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        ! compute mean b_phi here
                        tmp = 0.0d0
                        Do p = 1, n_phi
                            tmp = tmp+buffer(p,r,t,bphi)
                        Enddo
                        tmp = tmp*over_n_phi
                        tmp = tmp**2
                        qty(:,r,t,5) = tmp
                    Enddo
                Enddo

                ! Meridional ME
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        tmp = 0.0d0
                        tmp2 = 0.0d0
                        Do p = 1, n_phi
                            tmp = tmp+buffer(p,r,t,br)
                            tmp2 = tmp2+buffer(p,r,t,btheta)
                        Enddo
                        tmp = tmp*over_n_phi
                        tmp2 = tmp2*over_n_phi
                        tmp = (tmp**2+tmp2**2)
                        qty(:,r,t,6) = tmp
                    Enddo
                Enddo

            Endif ! (magnetism)

            Allocate(ell0_values(my_r%min:my_r%max,1:num_int))
            Allocate(volume_integrals(1:num_int), volume_sdev(1:num_int))
            Allocate(obs_sdev(1:nobs))
            Call ComputeEll0(qty,ell0_values)   ! Requires communication across row
            DeAllocate(qty)

            Call Compute_Radial_Average(ell0_values,volume_integrals) !Requires communication across column
            DeAllocate(ell0_values)




            !////////////////////////////////////////////////////////
            Call Assemble_Strips(buffer)
            If (my_rank .eq. 0) Then
                !Now, add these into the time_series and update the counter
                If (numt_ind .gt. max_numt) numt_ind = 1


                Call Point_Observations(current_time)
                If (.not. have_reference_strips) Then
                    xref(:) = xnow(:)
                    xlast(:) = xnow(:)
                    drift_reference_time = current_Time
                    previous_time = current_time
                    have_reference_strips = .true.
                Endif

                Do i = 1, nobs
                    obs_series(numt_ind,i) = observations(i)
                Enddo

                time_series(numt_ind,1:num_int) = volume_integrals(1:num_int)
                time_saves(numt_ind) = current_time
                iter_saves(numt_ind) = iteration





                numt_ind = numt_ind+1
                global_count = MIN(global_count+1, max_numt)

                !The first drift measurement is always garbage
                !  - duplicate 2nd measurement for averaging purposes
                If (global_count .eq. 2) drifts(1,1:2) = drifts(2,1:2)

                If (Mod(iteration,report_interval) .eq. 0) Then
                    ! Generate a benchmark report

                    ! Re-task the volume_integrals array
                    Do i = 1, num_int
                        Call get_moments(time_series(1:global_count,i),mean_value,sdev_value)
                        volume_integrals(i) = mean_value ! SUM(time_series(1:global_count,i))/global_count
                        volume_sdev(i) = sdev_value
                    Enddo
                    Do i = 1, nobs
                        Call get_moments(obs_series(1:global_count,i),mean_value,sdev_value)
                        observations(i) = mean_value
                        obs_sdev(i) = sdev_value
                    Enddo

                    Write(dt_str,fmtstr)( MAXVAL(time_saves(1:global_count))- &
                        & MINVAL(time_saves(1:global_count)) )
                    iter_start = MINVAL(iter_saves(1:global_count))
                    iter_end   = MAXVAL(iter_saves(1:global_count))


                    Write(iter_string,'(i8.8)')iteration
                    funit = 88
                    report_file = Trim(my_path)//'Benchmark_Reports/'//TRIM(iter_string)
                    Open(unit = funit, file = report_file,action="write", status="REPLACE", FORM = 'FORMATTED')

            Write(funit,*)'///////////////////////////////////////////////////////////////////////////'
            Write(funit,*)'              RAYLEIGH ACCURACY BENCHMARK SUMMARY               '
            Write(funit,*)' '
            Write(funit,*)'  Benchmark:  '//TRIM(benchmark_name)
            Write(funit,*)' '
            Write(funit,*)'  Radial Resolution      N_R = ', N_R
            Write(funit,*)'  Angular Resolution N_theta = ', n_theta
            Write(funit,*)' '
            If (benchmark_mode .lt. 3) Then
                Write(funit,*)'  Averaging Interval (Viscous Diffusion Times) : ', dt_str
            Else
                Write(funit,*)'  Averaging Interval (seconds) : ', dt_str
            Endif
            Write(funit,*)' '
            Write(funit,*)'  Beginning Iteration : ', iter_start
            Write(funit,*)'  Ending Iteration    : ', iter_end
            Write(funit,*)'  Number of Samples   : ', global_count
            Write(funit,*)'----------------------------------------------------------------------------'
            Write(funit,*)'  Observable      |    Measured    | Suggested   | % Difference |  Std. Dev.'
            Write(funit,*)'----------------------------------------------------------------------------'

                    If (benchmark_mode .eq. 1) Then
                        report_vals(1) = volume_integrals(1)
                        report_vals(2) = observations(3)
                        report_vals(3) = observations(2)


                        report_sdev(1) = volume_sdev(1)
                        report_sdev(2) = obs_sdev(3)
                        report_sdev(3) = obs_sdev(2)

                        Call get_moments(drifts(1:global_count,2),mean_value,sdev_value)
                        report_vals(4) = mean_value
                        report_sdev(4) = sdev_value

                    Endif


                    If (benchmark_mode .eq. 2) Then
                        report_vals(1) = volume_integrals(1)
                        report_vals(2) = volume_integrals(4)*mag_factor
                        report_vals(3) = observations(3)
                        report_vals(4) = observations(2)
                        report_vals(5) = observations(4)

                        report_sdev(1) = volume_sdev(1)
                        report_sdev(2) = volume_sdev(4)*mag_factor
                        report_sdev(3) = obs_sdev(3)
                        report_sdev(4) = obs_sdev(2)
                        report_sdev(5) = obs_sdev(4)

                        Call get_moments(drifts(1:global_count,2),mean_value,sdev_value)
                        report_vals(6) = mean_value
                        report_sdev(6) = sdev_value



                    Endif

                    If (benchmark_mode .eq. 3) Then
                        report_vals(1) = volume_integrals(1)*volume_norm
                        report_vals(2) = volume_integrals(2)*volume_norm
                        report_vals(3) = volume_integrals(3)*volume_norm
                        report_vals(4) = observations(3)
                        report_vals(5) = observations(2)

                        report_sdev(1) = volume_sdev(1)*volume_norm
                        report_sdev(2) = volume_sdev(2)*volume_norm
                        report_sdev(3) = volume_sdev(3)*volume_norm
                        report_sdev(4) = obs_sdev(3)
                        report_sdev(5) = obs_sdev(2)

                        Call get_moments(drifts(1:global_count,2),mean_value,sdev_value)
                        report_vals(6) = mean_value
                        report_sdev(6) = sdev_value
                        fmtstr = '(ES14.6)'
                    Endif


                    If (benchmark_mode .eq. 4) Then
                        Do i = 1, 6
                            report_vals(i) = volume_integrals(i)*volume_norm
                            report_sdev(i) = volume_sdev(i)*volume_norm
                        Enddo
                        Do i = 4,6
                            report_vals(i) = report_vals(i)*mag_factor
                            report_sdev(i) = report_sdev(i)*mag_factor
                        Enddo


                        report_vals(7) = observations(3)
                        report_vals(8) = observations(2)
                        report_vals(9) = observations(4)

                        !Btheta can be either positive or negative.
                        !If it's negative, switch to positive so that
                        !relative error makes sense.
                        If (report_vals(9) .lt. 0.0d0) Then
                            report_vals(9) = -report_vals(9)
                        Endif

                        report_sdev(7) = obs_sdev(3)
                        report_sdev(8) = obs_sdev(2)
                        report_sdev(9) = obs_sdev(4)

                        Call get_moments(drifts(1:global_count,2),mean_value,sdev_value)
                        report_vals(10) = mean_value
                        report_sdev(10) = sdev_value
                        fmtstr = '(ES14.6)'
                    Endif


                    Do i = 1, num_rep
                        rel_diff = (report_vals(i)-suggested_vals(i))/suggested_vals(i)*100
                        Write(val_str,fmtstr)report_vals(i)
                        Write(sug_str,fmtstr)suggested_vals(i)
                        Write(rel_str,fmtstr)rel_diff
                        Write(sdev_str,fmtstr)report_sdev(i)
                        Write(funit,*)TRIM(report_names(i)), val_str,sug_str,rel_str,sdev_str
                    Enddo

                    Close(funit)


                Endif
            Endif
            DeAllocate(volume_integrals,volume_sdev, obs_sdev)
        Endif
        Endif
    End Subroutine Benchmark_Checkup

    Subroutine Assemble_Strips(inbuff)
        Real*8, Intent(In) :: inbuff(1:,my_r%min:,my_theta%min:,:)
        Real*8, Allocatable :: all_strips(:,:,:)
        Integer :: i,j, indst(3), ndata
        Integer :: rirqs(1:4), sirqs(1:4)

        Integer, Allocatable :: obs_inds(:)

        ! At the end of this routine, rank zero will have
        ! All four slices in memory
        Allocate(all_strips(1:n_phi,1:nobs,4))
        Allocate(obs_inds(1:nobs))
        obs_inds(1) = vr
        obs_inds(2) = vphi
        obs_inds(3) = tvar
        if (magnetism) obs_inds(4) = btheta

        indst(:) = 1
        ndata = n_phi*nobs !number of data points communicated during each send

        !Rank 0 Post receives as appropriate

        If (my_rank .eq. 0) Then
            Do j = 1, num_strips
                If (.not. have_strip(j)) Then
                    indst(3) = j
                    Call Ireceive(all_strips, rirqs(j), n_elements = ndata,source= strip_owners(j), &
                        &  tag=btags(j),grp = pfi%gcomm, indstart = indst)
                Endif
            Enddo
        Endif

        ! Ranks that own strips send them to rank 0 as appropriate
        Do j = 1, num_strips
            If (have_strip(j)) Then
                Do i = 1, nobs
                    all_strips(1:n_phi,i,1) = inbuff(1:n_phi,rvals(j),tvals(j),obs_inds(i))
                Enddo

                If (my_rank .ne. 0) Then
                    ! If this was rank 0, we only need to copy the strip into the strips array
                    indst(3) = j
                    Call ISend(all_strips, sirqs(j),n_elements = ndata, dest = 0, &
                        &  tag = btags(j), grp = pfi%gcomm)
                Endif
            Endif
        Enddo

        !Next, we wait on sends and receives to complete and compute the average
        If (my_rank .ne. 0) Then
            Do j = 1, num_strips
                If (have_strip(j)) Then
                    Call IWait(sirqs(j))
                Endif
            Enddo
        Else
            Do j = 1, num_strips
                If (.not. have_strip(j)) Then
                    Call IWait(rirqs(j))
                Endif
            Enddo



            strips(:,:) = 0.0d0
            Do j = 1, num_strips
                Do i = 1, nobs
                    strips(1:n_phi,i) = strips(1:n_phi,i)+all_strips(1:n_phi,i,j)
                Enddo
            Enddo
            strips = strips*(1.0d0/num_strips) ! convert to an average

        Endif
        DeAllocate(all_strips, obs_inds)

    End Subroutine Assemble_Strips

    Subroutine Point_Observations(time_in)
        Implicit None
        !Conduct point-wise observations of T, u_phi, and b_theta
        ! Observations taken where vr = 0 and dvr/dphi > 0
        ! Also calculates the drift velocity.
        ! There is a lot of (somewhat) confusing logic here related to
        ! dealing with the drift of vr=0 points needed to
        ! calculate the drift velocity properly
        Real*8, Intent(In) :: time_in
        Integer :: i, j, xind, im1
        Real*8 :: vrlast, vrnext, dvr, test,x, y2, y1, slope
        Real*8 :: domegadt2, domegadt
        Real*8 :: delta_time, delta_time_pair, xadd, thisx
        Real*8, Allocatable :: vsave(:,:)
        Logical :: adjustx
        Allocate(vsave(1:msymm,1:nobs))
        xind = 1
        im1 = n_phi
        Do i = 1, n_phi
            vrnext = strips(i,1)
            vrlast = strips(im1,1)
            dvr = vrnext-vrlast
            test = vrnext*vrlast
            If ( (test .lt. 0) .and. (xind .le. msymm) ) Then
                ! We have crossed a zero point
                If (dvr .gt. 0) Then
                    !dvr/dphi is positive at this zero point
                    x = -vrlast/dvr  ! zero crossing (relative to im1 = 0)
                    xnow(xind) = im1+x
                    If (xnow(xind) .gt. n_phi) xnow(xind) = xnow(xind)-n_phi ! x is always between 0 and n_phi
                    If (xnow(xind) .lt. 0) xnow(xind) = xnow(xind)+n_phi
                    Do j = 1, nobs
                        ! Linearly interpolate to find the values at the zero crossing of vr
                        ! Interpolated vr is zero by definition -- good sanity check to save
                        y1 = strips(im1,j)
                        y2 = strips(i,j)
                        slope = y2-y1
                        vsave(xind,j) = slope*x+y1
                    Enddo
                    xind = xind+1
                Endif
            Endif
            im1 = i
        Enddo

        If (xind .ne. (msymm+1)) Write(6,*)'ISSUE!'

        Do i = 1, nobs
            observations(i) = SUM(vsave(:,i))/msymm
        Enddo

        ! Check for wrapping around of the x coordinates
        adjustx = .false.
        Do i = 1, xind-1
            If (drift_sign .ge. 0) Then
                If (xnow(i) .le. xlast(i)) Then
                    adjustx = .true.
                Endif
            Else
                If (xnow(i) .ge. xlast(i)) Then
                    adjustx = .true.
                Endif
            Endif
        Enddo

        delta_time = time_in - drift_reference_time
        delta_time_pair = time_in - previous_time
        domegadt = 0.0d0
        domegadt2 = 0.0d0
        xadd = n_phi/dble(msymm)


        Do i = 1, xind-1
            thisx = xnow(i)
            If (adjustx) Then
                thisx = thisx + xadd*drift_sign
            Endif
            domegadt = domegadt+(thisx-xref(i))
            domegadt2 = domegadt2 + (thisx-xlast(i))
            xlast(i) = xnow(i)
            previous_time = time_in
        Enddo
        If (adjustx) Then
            xref(:) = xnow(:)
            drift_reference_time = time_in
            !Write(6,*)'Adjusting X'
        Endif

        domegadt = domegadt/(msymm*delta_time*n_phi)*two_pi
        domegadt2 = domegadt2/(msymm*delta_time_pair*n_phi)*two_pi
        drifts(numt_ind,1) = domegadt
        drifts(numt_ind,2) = domegadt2

        DeAllocate(vsave)
    End Subroutine Point_Observations


    Subroutine Write_Array(arr,filename)
        Implicit None
        Character*120, Optional, Intent(In) :: filename
        Integer :: i,j,sig = 314, nx,ny
        Real*8, Intent(In) :: arr(1:,1:)
        nx = size(arr,1)
        ny = size(arr,2)

        Open(unit=15,file=filename,form='unformatted', status='replace',access='stream')
        Write(15)sig
        Write(15)nx
        Write(15)ny
        Write(15)((arr(i,j),i=1,nx),j = 1, ny)

        Close(15)

    End Subroutine Write_Array
    Subroutine Get_Moments(arr,arr_mean,arr_sdev)
        Implicit None
        Real*8, Intent(In) :: arr(1:)
        Real*8, Intent(InOut) :: arr_mean, arr_sdev
        Integer :: i, n
        ! Computes the mean and standard deviation of arr
        arr_mean = 0.0d0
        n = size(arr)

        Do i = 1, n
            arr_mean = arr_mean+arr(i)
        Enddo
        arr_mean = arr_mean/n

        arr_sdev = 0.0d0
        Do i = 1, n
            arr_sdev = arr_sdev+(arr(i)-arr_mean)**2
        Enddo
        arr_sdev = arr_sdev/n
        arr_sdev = sqrt(arr_sdev)
    End Subroutine Get_Moments
End Module Benchmarking
