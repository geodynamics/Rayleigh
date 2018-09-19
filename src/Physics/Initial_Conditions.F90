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

Module Initial_Conditions
    Use ProblemSize
    Use Fields
    Use Parallel_Framework
    Use Fourier_Transform
    Use Legendre_Transforms, Only : Legendre_Transform
    Use SendReceive
    Use Math_Constants
    Use Checkpointing, Only : read_checkpoint, read_checkpoint_alt
    Use Controls
    Use Timers
    Use General_MPI, Only : BCAST2D
    Use ReferenceState, Only : s_conductive, heating_type,ref
    Use BoundaryConditions, Only : T_top, T_bottom, fix_tvar_Top, fix_tvar_bottom,&
         & fix_dtdr_top, fix_dtdr_bottom, dtdr_top, dtdr_bottom, &
         & C10_bottom, C11_bottom, C1m1_bottom
    Use ClockInfo, Only : Euler_Step
    Use TransportCoefficients, Only : kappa, dlnkappa
    Use Linear_Solve
    Use Math_Utility
    Use BufferedOutput
    Use Load_Balance, Only : l_lm_values, my_num_lm

    Implicit None
    Logical :: alt_check = .false.
    Integer :: init_type = 1
    Integer :: magnetic_init_type = 1
    Integer :: init_tag = 8989
    Integer :: restart_iter = 0
    Real*8 :: temp_amp = 1.0d0, temp_w = 0.3d0, mag_amp = 1.0d0
    Logical :: conductive_profile = .false.
    Logical :: rescale_velocity = .false.
    Logical :: rescale_bfield = .false.
    Logical :: rescale_pressure = .false.
    Logical :: rescale_tvar = .false.
    Logical :: rescale_entropy = .false.
    Real*8  :: velocity_scale = 1.0d0
    Real*8  :: bfield_scale = 1.0d0
    Real*8  :: tvar_scale = 1.0d0
    Real*8  :: pressure_scale = 1.0d0
    Real*8  :: mdelta = 0.0d0  ! mantle convection benchmark delta

    Namelist /Initial_Conditions_Namelist/ init_type, temp_amp, temp_w, restart_iter, &
            & magnetic_init_type,alt_check, mag_amp, conductive_profile, rescale_velocity, &
            & rescale_bfield, velocity_scale, bfield_scale, rescale_tvar, &
            & rescale_pressure, tvar_scale, pressure_scale, mdelta
Contains

    Subroutine Initialize_Fields()
        Implicit None
        Logical :: dbtrans, dbconfig
        Logical :: test_reduce = .true.
        ! When coming out of this routine, the RHS of the equation set should contain the field values.
        ! This setup is consistent with the program having just completed a time step

        ! wsp%p1b should contain the Adams-Bashforth terms from the previous (not current)
        ! values of WPST.  This means that wsp%p1b should be zero if this init is from scratch.
        ! If this init is from restart, wsp%p1b should contain the adams bashforth terms output
        ! as part of the checkpoint.

        ! Check control variables to see if we need want static or buffers

        If (my_rank .eq. 0) Then
            Call stdout%print(" -- Initializing Fields...")
            Call stdout%print(" ---- Specified parameters: ")
            If (conductive_profile) Then
                Call stdout%print(" ---- Conductive entropy profile is selected. ")
            Endif
        Endif
        dbtrans = .not. static_transpose
        dbconfig = .not. static_config


        Call wsp%init(field_count = wsfcount, config = 'p1b', &
            dynamic_transpose =dbtrans, dynamic_config = dbconfig, &
            hold_cargo = test_reduce, padding = pad_alltoall, num_cargo = nglobal_msgs)
        Call wsp%construct('p1b')    ! We will always start in p1b - should do wsp%set_config('p1b')
        wsp%p1b(:,:,:,:) = 0.0d0    ! All fields are zero initially

        ! Allocate the Equation Set RHS
        ! Set it to zero initially
        ! The equation set RHS's stays allocated throughout - it is effectively how we save the AB terms.
        Call Allocate_RHS(zero_rhs=.true.)


        !////////////////////////////////////////
        ! Read in checkpoint files as appropriate
        If (init_type .eq. -1) Then
            If (my_rank .eq. 0) Then
                Call stdout%print(" ---- Hydro Init Type    : RESTART ")
            Endif
        Endif
        If (magnetism .and. (magnetic_init_type .eq. -1) ) Then
            If (my_rank .eq. 0) Then
                Call stdout%print(" ---- Magnetic Init Type : RESTART ")
            Endif
        Endif
        If ( (init_type .eq. -1) .or. ( magnetism .and. (magnetic_init_type .eq. -1) ) ) Then
            Call restart_from_checkpoint(restart_iter)
        Endif


        !////////////////////////////////////
        ! Initialize the hydro variables

        If (init_type .eq. 1) Then
            call benchmark_init_hydro()
            If (my_rank .eq. 0) Then
                Call stdout%print(" ---- Hydro Init Type    : Benchmark (Christensen et al. 2001) ")
            Endif
        Endif

        If (init_type .eq. 2) Then
            Call Mantle_Benchmark_Init()
            If (my_rank .eq. 0) Then
                Call stdout%print(" ---- Hydro Init Type    : Mantle Benchmark (Arrial et al. 2014) ")
            Endif
        Endif
        If (init_type .eq. 3) Then
            call diffusion_init_hydro()
            If (my_rank .eq. 0) Then
                Call stdout%print(" ---- Hydro Init Type    : Diffusion ")
            Endif
        Endif

        If (init_type .eq. 6) Then
            call abenchmark_init_hydro()
            If (my_rank .eq. 0) Then
                Call stdout%print(" ---- Hydro Init Type    : Benchmark (Jones et al. 2011) ")
            Endif
        Endif
        If (init_Type .eq. 7) Then
            If (my_rank .eq. 0) Then
                Call stdout%print(" ---- Hydro Init Type    : Random Thermal Field ")
            Endif
            call random_thermal_init()
        Endif


        If (magnetism) Then
            ! Initialize the magnetic variables
            If (magnetic_init_type .eq. 1) Then
                call benchmark_insulating_init()
                If (my_rank .eq. 0) Then
                    Call stdout%print(" ---- Magnetic Init Type : Benchmark (Christensen et al. 2001) ")
                Endif
            Endif
            If (magnetic_init_type .eq. 7) Then
                call random_init_Mag()
                If (my_rank .eq. 0) Then
                    Call stdout%print(" ---- Magnetic Init Type : Random Field")
                Endif
            Endif
            If (magnetic_init_type .eq. 10) Then
                call Dipole_Field_Init()
                If (my_rank .eq. 0) Then
                    Call stdout%print(" ---- Magnetic Init Type : Dipole Field")
                Endif
            Endif
        Endif
        ! Fields are now initialized and loaded into the RHS.
        ! We are ready to enter the main loop
        If (my_rank .eq. 0) Then
            Call stdout%print(" -- Fields initialized.")
            Call stdout%print(" ")
        Endif
    End Subroutine Initialize_Fields

    Subroutine Restart_From_Checkpoint(iteration)
        Implicit None
        Integer, Intent(In) :: iteration
        type(SphericalBuffer) :: tempfield
        Integer :: fcount(3,2), rpars(1:2),prod
        Integer :: this_ell, lm
        Character*14 :: scstr
        Character*8 ::  scfmt ='(ES10.4)'
        !rpars(1) = 1 if hydro variables are to be read (0 otherwise)
        !rpars(2) = 1 if magnetic variables are to be read (0 otherwise)
        rpars(1:2) = 0

        If (magnetism) Then

            If (init_type .eq. -1) rpars(1) = 1
            If (magnetic_init_type .eq. -1) rpars(2) = 1

            !If both variable types are not read in, an euler_step is taken on restart
            prod = rpars(1)*rpars(2)
            if (prod .eq. 0) euler_step = .true.

        Else
            If ( init_type .eq. -1) Then
                rpars(1) = 1
            Endif
        Endif
        !//////////////////////


        ! This routine also reads in the relevant magnetic quantities
        ! They are overwritten later by whatever the magnetic initialization does
        fcount(:,:) = 4
        If (magnetism) Then
            fcount(:,:) = 6
        Endif

        Call tempfield%init(field_count = fcount, config = 'p1a')
        Call tempfield%construct('p1a')

        wsp%p1b(:,:,:,:) = 0.0d0
        tempfield%p1a(:,:,:,:) = 0.0d0

        Call StopWatch(cread_time)%StartClock()
        If (read_chk_type .eq. 2) Then
            Call Read_Checkpoint_Alt(tempfield%p1a,wsp%p1b,iteration,rpars)
        Else
            Call Read_Checkpoint(tempfield%p1a,wsp%p1b,iteration,rpars)
        Endif
        Call StopWatch(cread_time)%Increment()

        If (rescale_velocity) Then
            euler_step = .true.
            tempfield%p1a(:,:,:,wvar) = tempfield%p1a(:,:,:,wvar)*velocity_scale
            tempfield%p1a(:,:,:,zvar) = tempfield%p1a(:,:,:,zvar)*velocity_scale
            wsp%p1b(:,:,:,:) = 0.0d0

            If (my_rank .eq. 0) Then
                Write(scstr,scfmt)velocity_scale
                Call stdout%print(" Rescaling velocity field by: "//scstr)
            Endif
        Endif
        If (rescale_bfield) Then
            euler_step = .true.
            tempfield%p1a(:,:,:,cvar) = tempfield%p1a(:,:,:,cvar)*bfield_scale
            tempfield%p1a(:,:,:,avar) = tempfield%p1a(:,:,:,avar)*bfield_scale
            wsp%p1b(:,:,:,:) = 0.0d0
            If (my_rank .eq. 0) Then
                Write(scstr,scfmt)bfield_scale
                Call stdout%print(" Rescaling magnetic field by: "//scstr)
            Endif
        Endif
        If (rescale_pressure) Then
            ! We do not rescale the ell = 0 mode
            euler_step = .true.
            Do lm = 1, my_num_lm
                this_ell = l_lm_values(lm)
                If (this_ell .gt. 0) Then
                    tempfield%p1a(:,:,lm,pvar) = tempfield%p1a(:,:,lm,pvar)*pressure_scale
                Endif
            Enddo
            wsp%p1b(:,:,:,:) = 0.0d0
            If (my_rank .eq. 0) Then
                Write(scstr,scfmt)pressure_scale
                Call stdout%print(" Rescaling magnetic field by: "//scstr)
            Endif
        Endif
        If (rescale_tvar) Then
            ! We do not rescale the ell = 0 mode
            euler_step = .true.
            Do lm = 1, my_num_lm
                this_ell = l_lm_values(lm)
                If (this_ell .gt. 0) Then
                    tempfield%p1a(:,:,lm,tvar) = tempfield%p1a(:,:,lm,tvar)*tvar_scale
                Endif
            Enddo
            wsp%p1b(:,:,:,:) = 0.0d0
            If (my_rank .eq. 0) Then
                Write(scstr,scfmt)tvar_scale
                Call stdout%print(" Rescaling thermal field (ell > 0) by: "//scstr)
            Endif
        Endif


        Call Set_All_RHS(tempfield%p1a)
        Call tempfield%deconstruct('p1a')

    End Subroutine Restart_From_Checkpoint



    !///////////////////////////////////////////////////////////
    !       Random Perturbation Initializaton Routines
    Subroutine Generate_Random_Field(rand_amp, field_ind, infield,rprofile, &
                & ell0_profile)
        Implicit None
        Integer :: ncombinations, i, m, r, seed(1), mp,n, l, ind1, ind2
        Integer :: mode_count, my_mode_start, my_mode_end, fcount(3,2)
        Integer, Intent(In) :: field_ind
        Real*8, Intent(In) :: rand_amp
        Real*8, Intent(In), Optional :: rprofile(my_r%min:), ell0_profile(1:)
        Real*8, Allocatable :: rand(:,:), rfunc(:), lpow(:)
        Real*8 :: amp, phase, lmid, alpha,x

        type(SphericalBuffer), Intent(InOut) :: infield
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 1


        Allocate(rfunc(my_r%min: my_r%max))
        If (present(rprofile)) Then
            rfunc(:) = rprofile(:)
        Else

            Do r = my_r%min, my_r%max
                x = 2.0d0*pi*(radius(r)-r_inner)/(r_outer-r_inner)
                rfunc(r) = 0.5d0*(1.0d0-Cos(x))
            Enddo
        Endif


        ! We put our temporary field in spectral space
        Call tempfield%init(field_count = fcount, config = 's2b')
        Call tempfield%construct('s2b')


        !///////////////////////
        ncombinations = 0
        Do m = 0, l_max
            ncombinations = ncombinations+ (l_max-m+1)
        Enddo

        !Set up the random phases and amplitudes
        Allocate(rand(1:ncombinations*2,1))

        If (my_rank .eq. 0) Then
            Call system_clock(seed(1))
            Call random_seed()
            Call random_number(rand)

            Do i = 1, ncombinations
                rand(i,1) = 2*temp_amp*(rand(i,1)-0.5d0)        ! first half of rand contains the amplitude
            Enddo
            ! We leave the second half alone (contains phases)

            ! Send rand
            !Do n = 1, ncpu -1
            !    Call send(rand, dest = n,tag=init_tag, grp=pfi%gcomm)
            !Enddo
        ENDIF
        !Else
            ! receive rand
        !    Call receive(rand, source= 0,tag=init_tag,grp = pfi%gcomm)
        !Endif

            If (my_row_rank .eq. 0) Then
                ! Broadcast along the column
                Call BCAST2D(rand,grp = pfi%ccomm)
            Endif
            Call BCAST2D(rand,grp = pfi%rcomm)


        ! Everyone establishes their range of random phases
        mode_count = 0
        Do mp = 1, my_mp%max
            if (mp .eq. my_mp%min) then
                my_mode_start = mode_count+1
            endif
            m = m_values(mp)
            mode_count = mode_count + (l_max-m+1)
            if (mp .eq. my_mp%max) then
                my_mode_end = mode_count
            endif
        Enddo

        Allocate(lpow(0:l_max))
        lmid = l_max/2.0d0
        alpha = lmid/3.0d0
        Do l = 0, l_max
                lpow(l) = rand_amp*exp(- ((l-lmid)/alpha )**2)
        Enddo


        ind1 = my_mode_start
        ind2 = ind1+ncombinations
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            Do l = m, l_max
                tempfield%s2b(mp)%data(l,:,:,:) = 0.0d0
                amp = rand(ind1,1)*lpow(l)
                phase = rand(ind2,1)
                ind1 = ind1+1
                ind2 = ind2+1
                Do r = my_r%min, my_r%max
                    tempfield%s2b(mp)%data(l,r,1,1) = tempfield%s2b(mp)%data(l,r,1,1) + &
                        amp*rfunc(r)*phase  ! real part
                    tempfield%s2b(mp)%data(l,r,2,1) =  &
                         & tempfield%s2b(mp)%data(l,r,2,1) + &
                         & amp*rfunc(r)*(1.0d0-phase) ! imaginary part
                Enddo
            Enddo
            if (m .eq. 0) Then
                ! Ell = 0 modes have no imaginary component
                Do r = my_r%min, my_r%max
                    tempfield%s2b(mp)%data(0,r,:,1) = 0.0d0
                Enddo
                If (present(ell0_profile)) Then
                    ! replace the ell = 0 profile
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(0,r,1,1) = ell0_profile(r)*sqrt(4.0d0*pi)
                    Enddo
                Endif
            Endif
        Enddo
        DeAllocate(rfunc, lpow)

        Call tempfield%reform() ! goes to p1b

        If (chebyshev) Then
            ! we need to load the chebyshev coefficients, and not the physical representation into the RHS
            Call tempfield%construct('p1a')

            Call gridcp%To_Spectral(tempfield%p1b,tempfield%p1a)

            tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
            Call tempfield%deconstruct('p1a')
        Endif
        infield%p1b(:,:,:,field_ind) = tempfield%p1b(:,:,:,1)
        Call tempfield%deconstruct('p1b')
        !Call tempfield%obliterate()
    End Subroutine Generate_Random_Field

    Subroutine Random_Init_Mag()
        Implicit None
        Real*8 :: ampa, ampc, dr_fiducial
        Real*8, Allocatable :: zero_profile(:)
        Integer :: fcount(3,2)
        type(SphericalBuffer) :: a_and_c
        fcount(:,:) = 2

        dr_fiducial = (radius(1)-radius(N_r))/dble(n_r)
        ampa = dr_fiducial*mag_amp
        ampc = dr_fiducial*ampa

        ! Construct the streamfunction field buffer
        Call a_and_c%init(field_count = fcount, config = 'p1b')
        Call a_and_c%construct('p1b')

        Allocate(zero_profile(1:N_R))
        zero_profile = 0.0d0
        ! Randomize each field
        ! neither of the magnetic potentials has an ell=0 component (zero it out)
        Call Generate_Random_Field(ampa, 1, a_and_c,ell0_profile = zero_profile)
        Call Generate_Random_Field(ampc, 2, a_and_c,ell0_profile = zero_profile)
        DeAllocate(zero_profile)
        Call Set_RHS(aeq,a_and_c%p1b(:,:,:,1))
        Call Set_RHS(ceq,a_and_c%p1b(:,:,:,2))
        Call a_and_c%deconstruct('p1b')
        !Call a_and_c%obliterate()
    End Subroutine Random_Init_Mag

    Subroutine Random_Thermal_Init()
        ! Generates random initial thermal perturbations
        Implicit None
        Real*8 :: amp
        Real*8, Allocatable :: profile0(:)
        Integer :: fcount(3,2)
        type(SphericalBuffer) :: sbuffer
        fcount(:,:) = 1

        amp = temp_amp

        ! Construct the streamfunction field buffer
        Call sbuffer%init(field_count = fcount, config = 'p1b')
        Call sbuffer%construct('p1b')

        If (conductive_profile) Then
            Allocate(profile0(1:N_R))
            profile0(:) = 0.0d0
            If (allocated(s_conductive)) Then
                If (heating_type .eq. 0) Then
                    profile0(:) = t_bottom*s_conductive(:)
                Else
                    profile0(:) = s_conductive(:)
                Endif
            Else
                Allocate(s_conductive(1:N_R))
                !The conductive {S,T} profile depends on kappa and ref%heating, so do this here.
                Call Calculate_Conductive_Profile()
                profile0(:) = s_conductive(:)

            Endif
            ! Randomize the entropy
            Call Generate_Random_Field(amp, 1, sbuffer,ell0_profile = profile0)
            DeAllocate(profile0)

        Else



            ! Randomize the entropy
            Call Generate_Random_Field(amp, 1, sbuffer)
        Endif

        Call Set_RHS(teq,sbuffer%p1b(:,:,:,1))
        Call sbuffer%deconstruct('p1b')

    End Subroutine Random_Thermal_Init

    !//////////////////////////////////////////////////////////////////////////////////
    !  Diffusion Init (for linear solve development)
    !  Initializes the Toroidal Stream Function (Z)
    Subroutine Diffusion_Init_Hydro()
        Implicit None
        Real*8, Allocatable :: rfunc(:)
        Real*8 :: x
        Integer :: r, l, m, mp
        Integer :: fcount(3,2)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 1

        Allocate(rfunc(my_r%min: my_r%max))

        Do r = my_r%min, my_r%max
            x = 2.0d0*pi*(radius(r)-r_inner)/(r_outer-r_inner)
            rfunc(r) = (1-cos(x))*0.5d0

        Enddo
        !write(6,*)'rf max ', maxval(rfunc1)

        ! We put our temporary field in spectral space
        Call tempfield%init(field_count = fcount, config = 's2b')
        Call tempfield%construct('s2b')

        ! Set the ell = 0 temperature and the real part of Y44
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            tempfield%s2b(mp)%data(:,:,:,:) = 0.0d0
            Do l = m, l_max
                if ( (l .eq. 1) .and. (m .eq. 1) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc(r)
                    Enddo
                endif
            Enddo
        Enddo
        DeAllocate(rfunc)

        Call tempfield%reform() ! goes to p1b
        If (chebyshev) Then
            ! we need to load the chebyshev coefficients, and not the physical representation into the RHS
            Call tempfield%construct('p1a')

            Call gridcp%To_Spectral(tempfield%p1b,tempfield%p1a)

            tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
            Call tempfield%deconstruct('p1a')
        Endif
        ! Set Z (toroidal stream function).  Leave the other fields alone
        Call Set_RHS(zeq,tempfield%p1b(:,:,:,1))

        Call tempfield%deconstruct('p1b')
    End Subroutine Diffusion_Init_Hydro

    !//////////////////////////////////////////////////////////////////////////////////
    !  Benchmark Initialization Routines
    Subroutine Benchmark_Init_Hydro()
        Implicit None
        Real*8, Allocatable :: rfunc1(:), rfunc2(:)
        Real*8 :: x
        Integer :: r, l, m, mp
        Integer :: fcount(3,2)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 1

        Allocate(rfunc1(my_r%min: my_r%max))
        Allocate(rfunc2(my_r%min: my_r%max))

        Do r = my_r%min, my_r%max
            x = 2.0d0*radius(r)-r_inner-r_outer
            rfunc1(r) = 0.2d0*(1.0d0-3.0d0*x*x+3.0d0*x**4-x**6)
            rfunc2(r) = r_outer*r_inner/radius(r)-r_inner
        Enddo
        !write(6,*)'rf max ', maxval(rfunc1)

        ! We put our temporary field in spectral space
        Call tempfield%init(field_count = fcount, config = 's2b')
        Call tempfield%construct('s2b')

        ! Set the ell = 0 temperature and the real part of Y44
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            tempfield%s2b(mp)%data(:,:,:,:) = 0.0d0
            Do l = m, l_max
                if ( (l .eq. 4) .and. (m .eq. 4) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc1(r)
                    Enddo
                endif

                if ( (l .eq. 0) .and. (m .eq. 0) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc2(r)*sqrt(4.0d0*pi)
                    Enddo
                endif
            Enddo
        Enddo
        DeAllocate(rfunc1,rfunc2)

        Call tempfield%reform() ! goes to p1b
        If (chebyshev) Then
            ! we need to load the chebyshev coefficients, and not the physical representation into the RHS
            Call tempfield%construct('p1a')

            Call gridcp%To_Spectral(tempfield%p1b,tempfield%p1a)

            tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
            Call tempfield%deconstruct('p1a')
        Endif
        ! Set temperature.  Leave the other fields alone
        Call Set_RHS(teq,tempfield%p1b(:,:,:,1))

        Call tempfield%deconstruct('p1b')
    End Subroutine Benchmark_Init_Hydro

    Subroutine ABenchmark_Init_Hydro()
        Implicit None
        Real*8, Allocatable :: rfunc1(:), rfunc2(:)
        Real*8 :: norm
        Integer :: r, l, m, mp
        Integer :: fcount(3,2)
        Real*8 :: d, beta, denom, zeta_0, c0, c1, n_rho, bm_n,ee, delta_s
        Real*8, Allocatable :: zeta(:)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 1

        Allocate(rfunc1(my_r%min: my_r%max))
        Allocate(rfunc2(my_r%min: my_r%max))
        !!!!!!!!!
        ! Stuff to set up background entropy gradient (taken from Mark's benchmark init in ASH)
        Delta_S = t_bottom ! bottom_entropy
        If (my_rank .eq. 0) Write(6,*)'Delta_S is ', delta_s
        d = radius(1) - radius(N_R)
        beta = radius(N_R) / radius(1)
        n_rho = 5.0d0
        bm_n = 2.0d0
        denom = beta * exp(N_rho / bm_n) + 1.d0
        zeta_0 = (beta+1.d0)/denom

        c0 = (2.d0 * zeta_0 - beta - 1.d0) / (1.d0 - beta)

        denom = (1.d0 - beta)**2
        c1 = (1.d0+beta)*(1.d0-zeta_0)/denom

        Allocate(zeta(1:N_R))

        zeta = c0 + c1 * d / Radius

        ee = -1.d0*bm_n
        denom = zeta(1)**ee - zeta(N_R)**ee



        !!!!!!!

        norm = 2.0d0*Pi/(radius(1)-radius(N_R))
        Do r = my_r%min, my_r%max

            rfunc2(r) = Delta_S * (zeta(1)**ee - zeta(r)**ee) / denom

            rfunc1(r) = (1.0d0-Cos(norm*(radius(r)-radius(N_R))))*temp_amp
        Enddo

        DeAllocate(zeta)
        ! We put our temporary field in spectral space
        Call tempfield%init(field_count = fcount, config = 's2b')
        Call tempfield%construct('s2b')

        ! Set the ell = 0 temperature and the real part of Y_19^19    and Y_1_1
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            tempfield%s2b(mp)%data(:,:,:,:) = 0.0d0
            Do l = m, l_max
                if ( (l .eq. 19) .and. (m .eq. 19) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc1(r)
                    Enddo
                endif
                if ( (l .eq. 1) .and. (m .eq. 1) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc1(r)*0.1d0
                    Enddo
                endif
                if ( (l .eq. 0) .and. (m .eq. 0) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc2(r)*sqrt(4.0d0*pi)
                    Enddo
                endif
            Enddo
        Enddo
        DeAllocate(rfunc1,rfunc2)

        Call tempfield%reform() ! goes to p1b
        If (chebyshev) Then
            ! we need to load the chebyshev coefficients, and not the physical representation into the RHS
            Call tempfield%construct('p1a')
            Call gridcp%to_Spectral(tempfield%p1b,tempfield%p1a)
            tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
            Call tempfield%deconstruct('p1a')
        Endif

        ! Set temperature.  Leave the other fields alone
        Call Set_RHS(teq,tempfield%p1b(:,:,:,1))

        Call tempfield%deconstruct('p1b')
    End Subroutine ABenchmark_Init_Hydro




    ! Initial Condition for cubic steady state
    ! P.-A Arrial, N. Flyer, G.B. Wright, L.H. Kellogg, 2014
    ! Geosci. Model Dev., 7 2065-2076
    Subroutine Mantle_Benchmark_Init()
        Implicit None
        Real*8, Allocatable :: rfunc1(:), rfunc2(:)
        Real*8 :: x
        Integer :: r, l, m, mp
        Integer :: fcount(3,2)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 1

        Allocate(rfunc1(my_r%min: my_r%max))
        Allocate(rfunc2(my_r%min: my_r%max))

        Do r = my_r%min, my_r%max

            x = (radius(r)-r_outer)/(r_inner-r_outer)
            rfunc1(r) = (r_inner/radius(r))*x

            x = (radius(r)-r_inner)/(r_outer-r_inner)
            rfunc2(r) = sin(pi*x)*0.01d0
        Enddo
        !write(6,*)'rf max ', maxval(rfunc1)

        ! We put our temporary field in spectral space
        Call tempfield%init(field_count = fcount, config = 's2b')
        Call tempfield%construct('s2b')

        ! Set the ell = 0 temperature and the real part of Y44
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            tempfield%s2b(mp)%data(:,:,:,:) = 0.0d0
            Do l = m, l_max
                if ( (l .eq. 4) .and. (m .eq. 4) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc2(r)*(1.0d0-mdelta)*(5.0d0/7.0d0)*sqrt(2.0d0)
                    Enddo
                endif

                if ( (l .eq. 4) .and. (m .eq. 0) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc2(r)
                    Enddo
                endif

                if ( (l .eq. 0) .and. (m .eq. 0) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc1(r)*sqrt(4.0d0*pi)
                    Enddo
                endif
            Enddo
        Enddo
        DeAllocate(rfunc1,rfunc2)

        Call tempfield%reform() ! goes to p1b
        If (chebyshev) Then
            ! we need to load the chebyshev coefficients, and not the physical representation into the RHS
            Call tempfield%construct('p1a')

            Call gridcp%To_Spectral(tempfield%p1b,tempfield%p1a)

            tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
            Call tempfield%deconstruct('p1a')
        Endif
        ! Set temperature.  Leave the other fields alone
        Call Set_RHS(teq,tempfield%p1b(:,:,:,1))

        Call tempfield%deconstruct('p1b')
    End Subroutine Mantle_Benchmark_Init





    Subroutine Dipole_Field_Init()
        Implicit None


        Integer :: r, l, m, mp
        Integer :: fcount(3,2)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 1


        ! We put our temporary field in spectral space
        Call tempfield%init(field_count = fcount, config = 's2b')
        Call tempfield%construct('s2b')

        ! Set the ell = 1, m = 0 component of C streamfunction to fall off as 1/r
        ! Set other components to zero
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            tempfield%s2b(mp)%data(:,:,:,:) = 0.0d0
            If (m .eq. 0) Then
                Do r = my_r%min, my_r%max
                    tempfield%s2b(mp)%data(1,r,1,1) = C10_bottom/radius(r)
                Enddo
            Endif

            If (m .eq. 1) Then
                Do r = my_r%min, my_r%max
                    tempfield%s2b(mp)%data(1,r,1,1) = C11_bottom/radius(r)
                    tempfield%s2b(mp)%data(1,r,2,1) = C1m1_bottom/radius(r)
                Enddo
            Endif

        Enddo


        Call tempfield%reform() ! goes to p1b
        If (chebyshev) Then
            ! we need to load the chebyshev coefficients, and not the physical representation into the RHS
            Call tempfield%construct('p1a')

            Call gridcp%To_Spectral(tempfield%p1b,tempfield%p1a)

            tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
            Call tempfield%deconstruct('p1a')
        Endif

        Call Set_RHS(ceq,tempfield%p1b(:,:,:,1))

        Call tempfield%deconstruct('p1b')
    End Subroutine dipole_field_init

    !////////////////////////////////////
    !  Magnetic Initialization
    Subroutine Benchmark_Insulating_Init()
        Implicit None
        Real*8, Allocatable :: rfunc1(:), rfunc2(:)
        Real*8 :: nrm1, nrm2
        Integer :: r, l, m, mp, roff
        Integer :: fcount(3,2)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 2

        Allocate(rfunc1(my_r%min: my_r%max))
        Allocate(rfunc2(my_r%min: my_r%max))

        nrm2 = sqrt(pi/5.0d0)*(4.0d0/3.0d0)*5.0d0
        nrm1 = (5.0d0/8.0d0)*sqrt(pi/3.0d0)
        Do r = my_r%min, my_r%max

            rfunc1(r) = 8.0d0*r_outer-6.0d0*radius(r)        ! functional form for c_1_0
            rfunc1(r) = rfunc1(r)-2.0d0*(r_inner**4)/(radius(r)**3)
            rfunc1(r) = rfunc1(r)*nrm1*radius(r)**2

            rfunc2(r) = sin(pi*(radius(r)-r_inner))*radius(r)*nrm2  ! function form for a_2_0
        Enddo
        !write(6,*)'rf max ', maxval(rfunc1)

        ! We put our temporary field in spectral space
        Call tempfield%init(field_count = fcount, config = 's2b')
        Call tempfield%construct('s2b')
        roff= 2*my_r%delta
        ! Set the ell = 0 temperature and the real part of Y44
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            tempfield%s2b(mp)%data(:,:,:,:) = 0.0d0
            Do l = m, l_max
                if ( (l .eq. 1) .and. (m .eq. 0) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc1(r)
                    Enddo
                endif

                if ( (l .eq. 2) .and. (m .eq. 0) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,2) = rfunc2(r)
                    Enddo
                endif
            Enddo
        Enddo
        DeAllocate(rfunc1,rfunc2)

        Call tempfield%reform() ! goes to p1b
        If (chebyshev) Then
            ! we need to load the chebyshev coefficients, and not the physical representation into the RHS
            Call tempfield%construct('p1a')


            Call gridcp%to_Spectral(tempfield%p1b,tempfield%p1a)

            tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
            Call tempfield%deconstruct('p1a')
        Endif

        Call Set_RHS(ceq,tempfield%p1b(:,:,:,1))
        Call Set_RHS(aeq,tempfield%p1b(:,:,:,2))

        Call tempfield%deconstruct('p1b')
    End Subroutine benchmark_insulating_init

    Subroutine Load_Radial_Profile(profile_file,profile_out)
        Implicit None
        Character*120, Intent(In) :: profile_file
        Real*8, Intent(InOut) :: profile_out(1:)
        Real*8, Allocatable :: radius_in(:), profile_in(:), spy2(:)
        Real*8 :: min_r_in, max_r_in, min_p_in, max_p_in
        Real*8 :: splx, sply
        Integer :: nr_in, r
        ! Reads in a 1-D radial profile of some quantity,
        ! interpolates it (using cubic splines) to the current grid,
        ! and stores it in profile_out.

      Open(unit=892, file = profile_file, form = 'unformatted', status = 'old')
      Read(892)nr_in
      Allocate(profile_in(1:nr_in))
      Allocate(radius_in(1:nr_in))
      Read(892)(radius_in(r),r=1,nr_in)
      Read(892)(profile_in(r),r = 1, nr_in)
      Close(892)

     !--------------------------------------------------------------
      ! Interpolate onto our grid (assume custom_radius and custom_entropy are backwards
      min_r_in = radius_in(nr_in)
      max_r_in = radius_in(1)
      max_p_in = profile_in(nr_in)
      min_p_in = profile_in(1)

      Allocate(spy2(1:nr_in))
      spy2(1:nr_in) = 0.0d0
      profile_out(1:N_R) = 0.0d0
      Call Spline(radius_in, profile_in, nr_in, 2.0D30, 2.0D30, spy2)
      Do r = 1, N_R
         If ( (radius(r) .le. max_r_in) .and. (radius(r) .ge. min_r_in) ) Then
            splx = radius(r)
            Call Splint(radius_in, profile_in,spy2,nr_in, splx, sply)
            profile_out(r) = sply
         Endif
      Enddo

      ! Take care of any out of bounds radial values
      Do r = 1, N_R
         If (radius(r) .ge. max_r_in) Then
            profile_out(r) = min_p_in
         Endif
         If (radius(r) .le. min_r_in) Then
            profile_out(r) = max_p_in
         Endif
      Enddo

        DeAllocate(radius_in, profile_in, spy2)
    End Subroutine Load_Radial_Profile



    Subroutine Calculate_Conductive_Profile
        Implicit None
        Integer :: i, r, old_code
        Real*8 :: amp, grav_r_ref, dr, vhint, qadd2, ftest
        Real*8 :: lum_top, lum_bottom, sfactor, diff, r2dr, qadd, dsdr_mean
        Real*8, Allocatable :: tmp1d(:), tmp1d2(:), tmp1d3(:)
        Allocate(tmp1d(1:N_R),tmp1d2(1:N_R))
        tmp1d(:) =0.0d0
        tmp1d2(:) = 0.0d0
        !Calculates the conductive entropy profile

        !The conductive entropy profile is one that is
        !in agreement with the boundary conditions
        !and whose associated flux balances any
        !internal heat sources/sinks.

        If (heating_type .gt. 0) Then
            tmp1d = ref%heating*ref%density*ref%temperature*r_squared
            !tmp1d is zero otherwise - i.e., no heating
        Endif
        Call Indefinite_Integral(tmp1d,tmp1d2,radius)
        !tmp1d2(r) is now int_rmin_r Q rho T r^2 dr
        !tmp1d2 is also now r^2*F_conductive + A {A is undetermined}

        If (fix_dtdr_top .or. fix_dtdr_bottom) Then
            If (fix_dtdr_bottom) Then
                ftest = -dtdr_bottom*kappa(N_R)
                ftest = ftest*ref%density(N_R)*ref%temperature(N_R)*r_squared(N_R)
                tmp1d2 = tmp1d2-tmp1d2(N_R)+ftest
            Endif

            If (fix_dtdr_top .and. (.not. fix_dtdr_bottom)) Then
                ftest = -dtdr_top*kappa(1)
                ftest = ftest*ref%density(1)*ref%temperature(1)*r_squared(1)
                tmp1d2 = tmp1d2-tmp1d2(1)+ftest
            Endif
            tmp1d = -tmp1d2*OneOverRSquared/(kappa*ref%density*ref%temperature)
            !tmp1d is now dsdr_conductive
            Call Indefinite_Integral(tmp1d,tmp1d2,radius)
            !tmp1d2 is now s_conductive + A {A is undetermined}
            tmp1d2 = tmp1d2-tmp1d2(1) ! set to zero at the top ; adjust as needed
            If (fix_tvar_top) Then
                tmp1d2 = tmp1d2+T_Top
            Endif
            If (fix_tvar_bottom) Then
                tmp1d2 = tmp1d2 - tmp1d2(N_R)+T_Bottom
            Endif
            s_conductive = tmp1d2
        Else
            ! T is fixed at top and bottom - this is marginally more complicated
            Allocate(tmp1d3(1:N_R))
            tmp1d3 = OneOverRSquared/(kappa*ref%density*ref%temperature)
            tmp1d = -tmp1d2*tmp1d3
            Call Indefinite_Integral(tmp1d,tmp1d2,radius)
            Call Indefinite_Integral(tmp1d3,tmp1d,radius)
            amp = (tmp1d2(1)-tmp1d2(N_R))  - (T_top-T_bottom)
            amp = amp/ (tmp1d(1)-tmp1d(N_R))
            s_conductive = tmp1d2-amp*tmp1d
            s_conductive = s_conductive-s_conductive(1)+T_Top
            DeAllocate(tmp1d3)

        Endif




        DeAllocate(tmp1d,tmp1d2)
    End Subroutine Calculate_Conductive_Profile



    Subroutine Restore_InitialCondition_Defaults()
        Implicit None

        alt_check = .false.
        init_type = 1
        magnetic_init_type = 1
        init_tag = 8989
        restart_iter = 0
        temp_amp = 1.0d0
        temp_w   = 0.3d0
        mag_amp  = 1.0d0
        conductive_profile = .false.

    End Subroutine Restore_InitialCondition_Defaults
End Module Initial_Conditions
