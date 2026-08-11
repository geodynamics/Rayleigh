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
    Use Checkpointing, Only : read_checkpoint
    Use Generic_Input, Only : read_input
    Use Controls
    Use Controls, Only : spin_horizontal
    Use Spin_Conversions, Only : spin_state_is_slot, Convert_Slot_Pair_To_Spin
    Use Timers
    Use General_MPI, Only : BCAST2D
    Use PDE_Coefficients
    Use BoundaryConditions, Only : T_top, T_bottom, fix_tvar_Top, fix_tvar_bottom,&
         & fix_dtdr_top, fix_dtdr_bottom, dtdr_top, dtdr_bottom, &
         & C10_bottom, C11_bottom, C1m1_bottom
    Use ClockInfo, Only : Euler_Step
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
    Real*8  :: temp_add = 0.0d0 
    Character*120 :: t_init_file = '__nothing__'
    Character*120 :: w_init_file = '__nothing__'
    Character*120 :: p_init_file = '__nothing__'
    Character*120 :: z_init_file = '__nothing__'
    Character*120 :: c_init_file = '__nothing__'
    Character*120 :: a_init_file = '__nothing__'
    Character*120 :: chi_a_init_file(1:n_scalar_max) = '__nothing__'
    Character*120 :: chi_p_init_file(1:n_scalar_max) = '__nothing__'
    Character*120 :: custom_thermal_file = '__nothing__'

    Namelist /Initial_Conditions_Namelist/ init_type, temp_amp, temp_w, restart_iter, &
            & magnetic_init_type,alt_check, mag_amp, conductive_profile, rescale_velocity, &
            & rescale_bfield, velocity_scale, bfield_scale, rescale_tvar, &
            & rescale_pressure, tvar_scale, pressure_scale, mdelta, &
            & t_init_file, w_init_file, p_init_file, z_init_file, &
            & c_init_file, a_init_file, custom_thermal_file, chi_a_init_file, chi_p_init_file, temp_add
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
        
        ! Take care of IC's that require reading from a checkpoint (if init_type and magnetic_init_type = -1)
        ! If (init_type or/and magnetic_init_type = -2), read the checkpoint and add a user defined field to it for IC
        If ( (init_type .lt. 0) .or. ( magnetism .and. (magnetic_init_type .lt. 0) ) ) Then
            Call restart_from_checkpoint(restart_iter)
            If (init_type .eq. -2) Then
                If (my_rank .eq. 0) Then
                    Call stdout%print(" ---- Hydro Init Type    : ADD USER FILE TO CHECKPOINT ")
                Endif
               Call add_to_field(tvar, t_init_file)
            Endif
            If (magnetic_init_type .eq. -2) Then
                If (my_rank .eq. 0) Then
                    Call stdout%print(" ---- Magnetic Init Type : ADD USER FILE TO CHECKPOINT ")
                Endif
               Call add_to_field(cvar, c_init_file)
               Call add_to_field(avar, a_init_file)
            Endif
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

        if (init_type .eq. 8) then
            If (my_rank .eq. 0) Then
                call stdout%print(" ---- Hydro Init Type    : Input File ")
            Endif
            call file_init()
        Endif

        if (init_type .eq. 9) then
            If (my_rank .eq. 0) Then
                call stdout%print(" ---- Hydro Init Type Compressible    : Benchmark (Jones et al. 2011) ")
            Endif
            If (thermal_variable .eq. 2) Then
                call Compressible_Init_Hydro_Entropy()
            Else
                call Compressible_Init_Hydro_Temperature()
            Endif
        Endif

        if (init_type .eq. 42) Then
            call ellzero_init()
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
            If (magnetic_init_type .eq. 8) Then
                call magnetic_file_init()
                If (my_rank .eq. 0) Then
                    Call stdout%print(" ---- Magnetic Init Type : Input File ")
                Endif
            Endif
            If (magnetic_init_type .eq. 10) Then
                call Dipole_Field_Init()
                If (my_rank .eq. 0) Then
                    Call stdout%print(" ---- Magnetic Init Type : Dipole Field")
                Endif
            Endif
        Endif

        If ((newtonian_cooling) .and. (newtonian_cooling_profile_file .ne. '__nothing__')) Then
            Allocate(newtonian_cooling_profile(1:N_R))
            newtonian_cooling_profile(:) = 0.0d0
            Call Load_Radial_Profile(newtonian_cooling_profile_file,newtonian_cooling_profile)
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
        Integer :: this_ell, lm, numfields
        Character*14 :: scstr
        Character*8 ::  scfmt ='(ES10.4)'
        !rpars(1) = 1 if hydro variables are to be read (0 otherwise)
        !rpars(2) = 1 if magnetic variables are to be read (0 otherwise)
        rpars(1:2) = 0

        If (magnetism) Then

            If (init_type .eq. -2) rpars(1) = 1
            If (init_type .eq. -1) rpars(1) = 1
            If (magnetic_init_type .eq. -1) rpars(2) = 1
            If (magnetic_init_type .eq. -2) rpars(2) = 1

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
        numfields = n_equations + n_active_scalars + n_passive_scalars
        if (magnetism) then
          numfields = numfields + 2
        end if
        fcount(:,:) = numfields
        If (magnetism) Then
            fcount(:,:) = numfields
        Endif

        Call tempfield%init(field_count = fcount, config = 'p1a')
        Call tempfield%construct('p1a')

        wsp%p1b(:,:,:,:) = 0.0d0
        tempfield%p1a(:,:,:,:) = 0.0d0

        Call StopWatch(cread_time)%StartClock()

        Call Read_Checkpoint(tempfield%p1a,wsp%p1b,iteration,rpars)

        Call StopWatch(cread_time)%Increment()

        If (compressible .and. spin_horizontal) Then
            ! v13.1: faithful spin restart.  Convert the pair STATE and the pair
            ! AB history slot -> q+/- at read time, via a scratch buffer with the
            ! writer's 2*numfields field count (the parked 2-field Convert_AB_Pair
            ! violated the transpose plan sizing; the checkpoint WRITER already
            ! proves the 2*numfields p1a->s2a transpose works).  The solve-side
            ! RHS then receives native q+/- before the first solve -- the old
            ! one-shot converted only the transform-side buffer, so the first
            ! restart solve mangled slot coefficients as q+/- (pair-only,
            ! second-step-onset corruption).  No Euler re-prime needed.
            Block
                Type(SphericalBuffer) :: cvt
                Integer :: cfc(3,2), mp
                cfc(:,:) = numfields*2
                Call cvt%init(field_count=cfc, config='p1a')
                Call cvt%construct('p1a')
                cvt%p1a(:,:,:,:) = 0.0d0
                cvt%p1a(:,:,:,vtheta) = tempfield%p1a(:,:,:,vtheta)
                cvt%p1a(:,:,:,vphi)   = tempfield%p1a(:,:,:,vphi)
                cvt%p1a(:,:,:,numfields+vteq) = wsp%p1b(:,:,:,vteq)
                cvt%p1a(:,:,:,numfields+vpeq) = wsp%p1b(:,:,:,vpeq)
                Call cvt%reform()   ! p1a -> s2a (pure transpose)
                Call Convert_Slot_Pair_To_Spin(cvt%s2a, vtheta, vphi)
                Call Convert_Slot_Pair_To_Spin(cvt%s2a, numfields+vteq, numfields+vpeq)
                Call cvt%construct('s2b')
                Do mp = my_mp%min, my_mp%max
                    cvt%s2b(mp)%data(:,:,:,:) = cvt%s2a(mp)%data(:,:,:,:)
                Enddo
                Call cvt%deconstruct('s2a')
                cvt%config = 's2b'
                Call cvt%reform()   ! s2b -> p1b (pure transpose)
                tempfield%p1a(:,:,:,vtheta) = cvt%p1b(:,:,:,vtheta)
                tempfield%p1a(:,:,:,vphi)   = cvt%p1b(:,:,:,vphi)
                wsp%p1b(:,:,:,vteq) = cvt%p1b(:,:,:,numfields+vteq)
                wsp%p1b(:,:,:,vpeq) = cvt%p1b(:,:,:,numfields+vpeq)
                Call cvt%deconstruct('p1b')
            End Block
            spin_state_is_slot = .false.
            If (my_rank .eq. 0) Write(6,*) 'spin_horizontal: restart pair state+AB converted slot -> q+/- (v13.1)'
        Endif

        If (rescale_velocity) Then
            euler_step = .true.
            If (compressible) Then
                tempfield%p1a(:,:,:,vr) = tempfield%p1a(:,:,:,vr)*velocity_scale
                tempfield%p1a(:,:,:,vtheta) = tempfield%p1a(:,:,:,vtheta)*velocity_scale
                tempfield%p1a(:,:,:,vphi) = tempfield%p1a(:,:,:,vphi)*velocity_scale
            Else 
                tempfield%p1a(:,:,:,wvar) = tempfield%p1a(:,:,:,wvar)*velocity_scale
                tempfield%p1a(:,:,:,zvar) = tempfield%p1a(:,:,:,zvar)*velocity_scale
            Endif 
            
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
        Integer :: ncombinations, i, m, r, seed(1), mp, l, ind1, ind2
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
            Call Generate_Random_Field(amp, 1, sbuffer,ell0_profile = profile0)
            DeAllocate(profile0)            
                    
        Else If (trim(custom_thermal_file) .ne. '__nothing__') then
            Allocate(profile0(1:N_R))
            profile0(:) = 0.0d0
            
            Call Load_Radial_Profile(custom_thermal_file,profile0)

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

    subroutine file_init()
        ! initialize hydrodynamic variables from generic input files
        implicit none
        integer :: i
        integer :: fcount(3,2)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 1
        call tempfield%init(field_count = fcount, config = 'p1b')
        call tempfield%construct('p1b')

        if (trim(w_init_file) .ne. '__nothing__') then
            call read_input(w_init_file, 1, tempfield)
            call set_rhs(weq, tempfield%p1b(:,:,:,1))
        end if
 
        if (trim(p_init_file) .ne. '__nothing__') then
            call read_input(p_init_file, 1, tempfield)
            call set_rhs(peq, tempfield%p1b(:,:,:,1))
        end if
 
        if (trim(t_init_file) .ne. '__nothing__') then
            call read_input(t_init_file, 1, tempfield)
            call set_rhs(teq, tempfield%p1b(:,:,:,1))
        end if

        if (trim(z_init_file) .ne. '__nothing__') then
            call read_input(z_init_file, 1, tempfield)
            call set_rhs(zeq, tempfield%p1b(:,:,:,1))
        end if

        do i = 1, n_active_scalars
          if (trim(chi_a_init_file(i)) .ne. '__nothing__') then
              call read_input(chi_a_init_file(i), 1, tempfield)
              call set_rhs(chiaeq(i), tempfield%p1b(:,:,:,1))
          end if
        end do

        do i = 1, n_passive_scalars
          if (trim(chi_p_init_file(i)) .ne. '__nothing__') then
              call read_input(chi_p_init_file(i), 1, tempfield)
              call set_rhs(chipeq(i), tempfield%p1b(:,:,:,1))
          end if
        end do

        call tempfield%deconstruct('p1b')

    end subroutine file_init

    subroutine magnetic_file_init()
        ! initialize magnetic variables from generic input files
        implicit none
        integer :: fcount(3,2)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 1
        call tempfield%init(field_count = fcount, config = 'p1b')
        call tempfield%construct('p1b')

        if (trim(c_init_file) .ne. '__nothing__') then
            call read_input(c_init_file, 1, tempfield)
            call set_rhs(ceq, tempfield%p1b(:,:,:,1))
        end if
 
        if (trim(a_init_file) .ne. '__nothing__') then
            call read_input(a_init_file, 1, tempfield)
            call set_rhs(aeq, tempfield%p1b(:,:,:,1))
        end if
 
        call tempfield%deconstruct('p1b')

    end subroutine magnetic_file_init

    Subroutine Add_To_Field(field_index, field_file)
        ! initialize magnetic variables from generic input files
        Implicit None
        Integer, Intent(In) :: field_index
        Character*120, Intent(In) :: field_file
        Integer :: fcount(3,2)
        Type(SphericalBuffer) :: tempfield, tempfield2
        fcount(:,:) = 1
        Call tempfield%init(field_count = fcount, config = 'p1b')
        Call tempfield%construct('p1b')
        Call tempfield2%init(field_count = fcount, config = 'p1b')
        Call tempfield2%construct('p1b')

        Call get_rhs(field_index,tempfield2%p1b(:,:,:,1))

        If (trim(field_file) .ne. '__nothing__') Then
            Call read_input(field_file, 1, tempfield)

            tempfield%p1b = tempfield%p1b+tempfield2%p1b
            Call set_rhs(field_index, tempfield%p1b(:,:,:,1))
        End If


        Call tempfield%deconstruct('p1b')
        Call tempfield2%deconstruct('p1b')

    End Subroutine add_to_field

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


    !////////////////////////////////
    ! Compressible debugging
    Subroutine ellzero_init()
        Implicit None
        Real*8, Allocatable :: rfunc(:)
        Real*8 :: x, alpha, q, biga, bigb
        Integer :: r, l, m, mp
        Integer :: fcount(3,2)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 1

        Allocate(rfunc(my_r%min: my_r%max))

        alpha = 50.0d0
        q = -1.0d0 ! gas_gamma/prandtl_number
        biga = (-q/3.0)*(r_inner**3)
        bigb = -q*((1.0d0/6.0d0)*r_outer**2 +(1.0d0/3.0d0)*(r_inner**3)/r_outer)
        Do r = my_r%min, my_r%max
            x = 2.0d0*pi*(radius(r)-r_inner)/(r_outer-r_inner)
            rfunc(r) = 0.5d0*(1.0d0-Cos(x))
            x = (radius(r)-r_inner)/(r_outer-r_inner)
            rfunc(r) = exp(-alpha*(x-0.5d0)**2)
            rfunc(r) = (q/6.0)*radius(r)**2 -BigA/Radius(r) + BigB 

        Enddo
        rfunc = rfunc-exp(-alpha*(-0.5d0)**2)
        !write(6,*)'rf max ', maxval(rfunc1)

        ! We put our temporary field in spectral space
        Call tempfield%init(field_count = fcount, config = 's2b')
        Call tempfield%construct('s2b')

        ! Set the ell = 0 temperature
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            tempfield%s2b(mp)%data(:,:,:,:) = 0.0d0
            Do l = m, l_max

                if ( (l .eq. 0) .and. (m .eq. 0) ) Then
                    Do r = my_r%min, my_r%max
                        tempfield%s2b(mp)%data(l,r,1,1) = rfunc(r)*sqrt(4.0d0*pi)*temp_amp+sqrt(4.0d0*pi)*temp_add
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
        ! Set temperature.  Leave the other fields alone
        Call Set_RHS(teq,tempfield%p1b(:,:,:,1))

        Call tempfield%deconstruct('p1b')
    End Subroutine ellzero_init

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

    ! ==========================================================================
    ! Reference listing: Compressible_Init_Hydro, entropy branch
    ! (thermal_variable == 2). Companion to rayleigh_entropy_handdoc_partA.md,
    ! Edit 5. Written in the tree's own idiom (descending radius array:
    ! radius(1) = r_outer, radius(N_R) = r_inner; global radius/ref arrays on
    ! every rank; spectral loading skeleton retained verbatim from the current
    ! routine). Hand-transcription reference -- adapt names to taste.
    !
    ! Uses only: radius, ref%temperature, ref%density, ref%dT, ref%dlnrho,
    ! gas_gamma, bigz. GENERAL BACKGROUND (adiabatic NOT assumed): the entropy
    ! zero-point constant Cs is set at the inner radius, the background entropy
    ! profile Sbar(r) = bigz*ln( Tbar rhobar^{1-gamma} / Cs ) is computed from the
    ! reference arrays and loaded into the l=0 field (identically zero for an
    ! adiabatic background, so j2011 is unchanged), and gravity comes from the
    ! general hydrostatic identity g = -R[ dTbar/dr + Tbar dlnrho/dr ],
    ! R = (gamma-1)*bigz. The H-gate verifies this matches the momentum
    ! routine's gravity discretely.
    ! ==========================================================================

    Subroutine Compressible_Init_Hydro_Entropy()
      Implicit None
      Real*8, Allocatable :: rfunc1(:), Scond(:), dScond(:), Lprof(:)
      Real*8, Allocatable :: Stot(:), dStot(:)
      Real*8, Allocatable :: Jarr(:), finteg(:), gprof(:), Sbar(:)
      Real*8, Allocatable :: Fw(:), Iw(:)
      Real*8 :: norm, DeltaS, dr   ! (Cs comes from PDE_Coefficients -- do NOT redeclare locally, it would shadow the module variable)
      Real*8 :: anchor, mass_ref, danchor, m1, m2, Rg, pex
      Integer :: r, l, m, mp, it
      Integer :: fcount(3,2)
      type(SphericalBuffer) :: tempfield
      fcount(:,:) = 2

      DeltaS = 851225.7d0
      If (nulltest_deltas_zero) DeltaS = 0.0d0   ! adiabatic null diagnostic
      ! v13.2: honor the namelist knob.  temp_amp (namelist) multiplies the
      ! historical hardcoded seed 1e-4*DeltaS, so temp_amp=1.0 reproduces
      ! every existing run bit-identically; temp_amp=1e-4 gives the linear
      ! (gate-G) seed 1e-8*DeltaS.  (The old line overwrote the namelist
      ! value -- the knob was dead for the entropy init.)
      temp_amp = temp_amp*1.0d-4*DeltaS  ! seed amplitude in entropy units
      
      Allocate(rfunc1(my_r%min:my_r%max))
      Allocate(Scond(1:N_R), dScond(1:N_R), Lprof(1:N_R), Sbar(1:N_R))
      Allocate(Jarr(1:N_R), finteg(1:N_R), gprof(1:N_R))
      Allocate(Fw(1:N_R), Iw(1:N_R))
      
      ! SELF-TEST of Cheby_Antiderivative conventions (index order, wall
      ! row, domain scaling) in one shot: F=1 must give X = r - r_inner.
      Fw(:) = 1.0d0
      Call Cheby_Antiderivative(Fw, Iw)
      If (my_rank .eq. 0) Then
         Write(6,*) ' entropy init: antideriv self-test max err = ', &
              MaxVal(Abs(Iw(:) - (radius(:) - radius(N_R))))
         ! expect ~1e-12 * (r_outer - r_inner)-scale; anything larger
         ! means an index/scaling convention mismatch -- fix before use.
      Endif
      
      !------------------------------------------------------------------
      ! (0) entropy zero-point + background entropy profile (GENERAL).
      !     S = cv ln( T rho^{1-gamma} / Cs );  Cs fixes S(inner adiabat)=0.
      Cs = ref%temperature(N_R)*ref%density(N_R)**(1.0d0-gas_gamma)
      Do r = 1, N_R
         Sbar(r) = bigz*Log( ref%temperature(r) &
              *ref%density(r)**(1.0d0-gas_gamma)/Cs )
      Enddo
      If (my_rank .eq. 0) Then
         Write(6,*) ' entropy init: Cs = ', Cs
         Write(6,*) ' entropy init: max|Sbar|/DeltaS = ', &
              MaxVal(Abs(Sbar))/DeltaS
         ! For an adiabatic background (j2011) this is ~1e-12 -- Sbar = 0.
         ! Nonzero values are fine: that IS the non-adiabatic background,
         ! carried in the total-field S.
      Endif
      
      !------------------------------------------------------------------
      ! (1) conductive entropy profile:  rho T r^2 dS/dr = const
      !     S(r) = DeltaS * (1 - J(ri,r)/J(ri,ro)),  J = int dr/(rho T r^2)
      Do r = 1, N_R
         finteg(r) = 1.0d0/(ref%density(r)*ref%temperature(r)*radius(r)**2)
      Enddo
      ! SPECTRAL integral: Jarr = antiderivative of finteg, zeroed at the
      ! inner wall (implementation below: gridcp%dcheby + dgesv).
      Call Cheby_Antiderivative(finteg, Jarr)
      Jarr(:) = Jarr(:) - Jarr(N_R)
      Do r = 1, N_R
         Scond(r)  = DeltaS*(1.0d0 - Jarr(r)/Jarr(1))
         dScond(r) = -DeltaS*finteg(r)/Jarr(1)      ! analytic dS/dr
      Enddo
        
      !------------------------------------------------------------------
      ! (2) gravity from GENERAL background hydrostatics:
      !     g = -(1/rhobar) dpbar/dr = -R[ dTbar/dr + Tbar dlnrhobar/dr ]
      Do r = 1, N_R
         gprof(r) = -(gas_gamma-1.0d0)*bigz*( ref%dT(r) &
              + ref%temperature(r)*ref%dlnrho(r) )
      Enddo

      !------------------------------------------------------------------
      ! (3) hydrostatic total ln(rho) with S = Stot = Sbar + Scond, via the
      !     EXACT LINEARIZATION w = p^((gamma-1)/gamma):
      !        dw/dr = -((gamma-1)/gamma) (R*Cs)^(-1/gamma) g e^{-S/(gamma cv)}
      !     -> single anchor-independent quadrature; then
      !        p = w^(gamma/(gamma-1)),
      !        L = ( ln(p/(R*Cs)) )/gamma - S/(gamma*bigz).
      !     Mass Newton = scalar secant on the inner anchor (no
      !     re-integration: the integral is anchor-independent).
      Allocate(Stot(1:N_R), dStot(1:N_R))
      Do r = 1, N_R
         Stot(r)  = Sbar(r) + Scond(r)
         dStot(r) = dScond(r) + bigz*( ref%dT(r)/ref%temperature(r) &
              + (1.0d0-gas_gamma)*ref%dlnrho(r) )
      Enddo

      Rg  = (gas_gamma-1.0d0)*bigz
      pex = (gas_gamma-1.0d0)/gas_gamma
      Do r = 1, N_R
         Fw(r) = gprof(r)*Exp(-Stot(r)/(gas_gamma*bigz))
      Enddo
      Call Cheby_Antiderivative(Fw, Iw)
      Iw(:) = Iw(:) - Iw(N_R)

      ! reference mass (use the grid's spectral integration weights if
      ! available; trapezoid shown for self-containment):
      mass_ref = 0.0d0
      Do r = 1, N_R-1
         dr = radius(r) - radius(r+1)
         mass_ref = mass_ref + 0.5d0*(ref%density(r)*radius(r)**2 &
              + ref%density(r+1)*radius(r+1)**2)*dr
      Enddo
      
      anchor  = Log(ref%density(N_R))                ! first guess: rhobar
      danchor = 1.0d-3
      Do it = 1, 12
         Call W_To_Mass(anchor,          m1)
         If (Abs(m1/mass_ref - 1.0d0) .lt. 1.0d-14) Exit
         Call W_To_Mass(anchor+danchor,  m2)
         anchor = anchor - (m1 - mass_ref)*danchor/(m2 - m1)
      Enddo
      Call W_To_Mass(anchor, m1)
      If (my_rank .eq. 0) Then
         Write(6,*) ' entropy init: mass Newton dM/M = ', m1/mass_ref-1.0d0
         Write(6,*) ' entropy init: L(inner), L(outer) = ', &
              Lprof(N_R), Lprof(1)
      Endif

      !------------------------------------------------------------------
      ! (4) seed radial shape (unchanged from current routine)
      norm = 2.0d0*Pi/(radius(1)-radius(N_R))
      Do r = my_r%min, my_r%max
         rfunc1(r) = (1.0d0-Cos(norm*(radius(r)-radius(N_R))))*temp_amp
      Enddo
      
      !------------------------------------------------------------------
      ! (5) spectral loading -- the tree's existing skeleton with the new
      !     l=0 payloads (Stot, Lprof) and the entropy-unit seed.
      Call tempfield%init(field_count = fcount, config = 's2b')
      Call tempfield%construct('s2b')
      
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
                  tempfield%s2b(mp)%data(l,r,1,1) = Stot(r)*sqrt(4.0d0*pi)
                  tempfield%s2b(mp)%data(l,r,1,2) = Lprof(r)*sqrt(4.0d0*pi)
               Enddo
            endif
         Enddo
      Enddo
      
      Call tempfield%reform() ! goes to p1b
      If (chebyshev) Then
         ! load chebyshev coefficients, not the physical representation
         Call tempfield%construct('p1a')
         Call gridcp%to_Spectral(tempfield%p1b,tempfield%p1a)
         tempfield%p1b(:,:,:,:) = tempfield%p1a(:,:,:,:)
         Call tempfield%deconstruct('p1a')
      Endif
      
      Call Set_RHS(teq,tempfield%p1b(:,:,:,1))
      Call Set_RHS(rhoeq,tempfield%p1b(:,:,:,2))
      Call tempfield%deconstruct('p1b')

      DeAllocate(rfunc1, Scond, dScond, Lprof, Jarr, finteg, gprof)
      DeAllocate(Sbar, Stot, dStot, Fw, Iw)
    Contains

      Subroutine W_To_Mass(a_in, mass_out)
        ! Closed-form profile family from the precomputed integral Iw:
        ! fills Lprof for inner anchor a_in and returns the shell mass.
        Implicit None
        Real*8, Intent(In)  :: a_in
        Real*8, Intent(Out) :: mass_out
        Real*8 :: wi0, wr0, drm
        Integer :: rr
        wi0 = ( Rg*Cs*Exp(Stot(N_R)/bigz + gas_gamma*a_in) )**pex
        Do rr = 1, N_R
           wr0 = wi0 - pex*(Rg*Cs)**(-1.0d0/gas_gamma)*Iw(rr)
           Lprof(rr) = ( Log(wr0)/pex - Log(Rg*Cs) )/gas_gamma &
                - Stot(rr)/(gas_gamma*bigz)
        Enddo
        mass_out = 0.0d0
        Do rr = 1, N_R-1
           drm = radius(rr) - radius(rr+1)
           mass_out = mass_out + 0.5d0*( Exp(Lprof(rr))*radius(rr)**2 &
                + Exp(Lprof(rr+1))*radius(rr+1)**2 )*drm
        Enddo
      End Subroutine W_To_Mass
          
    End Subroutine Compressible_Init_Hydro_Entropy

    ! ---------------------------------------------------------------------
    ! Cheby_Antiderivative: solve X' = F with X(r_inner) = 0 at spectral
    ! precision, using the grid's OWN stored Chebyshev matrices
    ! (gridcp%dcheby(1)%data(point, coeff, order): order 0 = value synthesis
    ! T_n(r_k), order 1 = d/dr synthesis, already domain-scaled). Unknowns
    ! are the Chebyshev coefficients c of X:
    !     [ dcheby(:,:,1) ] c = F        (derivative rows, every node)
    !     [ dcheby(iwall,:,0) ] c = 0    (value row replacing the inner-wall
    !                                     derivative row -- the BC that makes
    !                                     the otherwise-singular d/dr matrix
    !                                     invertible)
    ! then X = dcheby(:,:,0) c. Because these are the SAME matrices the
    ! solver's radial derivative uses, the IC becomes discretely balanced
    ! with respect to the code's own d/dr -- the exact null criterion.
    !
    ! If the coefficient dimension (N_max, dealiased) is smaller than N_R,
    ! the stacked system is (N_R x N_max) overdetermined: use dgels instead
    ! of dgesv (smooth integrands live in the dealiased range; residual is
    ! machine-level). Square case shown.
    Subroutine Cheby_Antiderivative(F_in, X_out)
      Implicit None
      Real*8, Intent(In)  :: F_in(1:N_R)
      Real*8, Intent(Out) :: X_out(1:N_R)
      Real*8  :: Amat(1:N_R,1:N_R), rhs(1:N_R)
      Integer :: ipiv(1:N_R), info
      Amat(1:N_R,1:N_R) = gridcp%dcheby(1)%data(1:N_R,1:N_R,1)
      rhs(1:N_R)        = F_in(1:N_R)
      Amat(N_R,1:N_R)   = gridcp%dcheby(1)%data(N_R,1:N_R,0)  ! inner-wall value row
      rhs(N_R)          = 0.0d0
      Call dgesv(N_R, 1, Amat, N_R, ipiv, rhs, N_R, info)
      If ((info .ne. 0) .and. (my_rank .eq. 0)) Then
         Write(6,*) ' Cheby_Antiderivative: dgesv info = ', info
      Endif
      X_out(1:N_R) = MatMul( gridcp%dcheby(1)%data(1:N_R,1:N_R,0), rhs )
    End Subroutine Cheby_Antiderivative
    
    Subroutine Compressible_Init_Hydro_Temperature()
        Implicit None
        Real*8, Allocatable :: rfunc1(:), rfunc2(:)
        Real*8 :: norm
        Integer :: r, l, m, mp
        Integer :: fcount(3,2)
        type(SphericalBuffer) :: tempfield
        fcount(:,:) = 2
        Allocate(rfunc1(my_r%min: my_r%max))
        Allocate(rfunc2(my_r%min: my_r%max))
        temp_amp = 1.0d0
        norm = 2.0d0*Pi/(radius(1)-radius(N_R))
        Do r = my_r%min, my_r%max
            
            rfunc2(r) = (4234)*((radius(N_R)/radius(r) - radius(N_R)/radius(1))/(1-radius(N_R)/radius(1)))
    
            rfunc1(r) = (1.0d0-Cos(norm*(radius(r)-radius(N_R))))*temp_amp
        Enddo
    ! We put our temporary field in spectral space
        Call tempfield%init(field_count = fcount, config = 's2b')
        Call tempfield%construct('s2b')
    ! Set the ell = 0 temperature and the real part of Y_19^19 and Y_1_1
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
                        tempfield%s2b(mp)%data(l,r,1,1) = ((ref%temperature(r)) + rfunc2(r))*sqrt(4.0d0*pi)
                        tempfield%s2b(mp)%data(l,r,1,2) = log(ref%density(r))*sqrt(4.0d0*pi)
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
        Call Set_RHS(rhoeq,tempfield%p1b(:,:,:,2))
        Call tempfield%deconstruct('p1b')
    End Subroutine Compressible_Init_Hydro_Temperature

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


        Integer :: r, m, mp
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
                    tempfield%s2b(mp)%data(1,r,1,1) = C10_bottom*rmin/radius(r)
                Enddo
            Endif

            If (m .eq. 1) Then
                Do r = my_r%min, my_r%max
                    tempfield%s2b(mp)%data(1,r,1,1) = C11_bottom*rmin/radius(r)
                    tempfield%s2b(mp)%data(1,r,2,1) = C1m1_bottom*rmin/radius(r)
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
        Real*8 :: ov_root_fourpi
        type(SphericalBuffer) :: tempfield
        integer :: fcount(3,2)
        integer :: r, m, mp
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Reads in a 1-D radial profile of some quantity,
        ! interpolates it (using cubic splines) to the current grid,
        ! and stores it in profile_out.

        fcount(:,:) = 1
        call tempfield%init(field_count = fcount, config = 'p1b')
        call tempfield%construct('p1b')        ! read_input expects field to be in p1b configuration
        call read_input(profile_file, 1, tempfield) ! Read the data
        
        Call tempfield%construct('p1a')  ! transform from Chebyshev space to radial grid
        Call gridcp%From_Spectral(tempfield%p1b,tempfield%p1a)
        tempfield%config='p1a'
        Call tempfield%reform()  ! move from p1a to s2a
        
        ! Next move to rlm space
        
        ! Extract ell=0 profile

        ov_root_fourpi = 1.0d0/sqrt(four_pi)
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            If (m .eq. 0) Then

                Do r = my_r%min, my_r%max
                    profile_out(r) = tempfield%s2a(mp)%data(0,r,1,1)*ov_root_fourpi
                Enddo
            Endif
        Enddo
        
        
        ! clean up
        Call tempfield%deconstruct('s2a')
        Call tempfield%deconstruct('p1b')
        
        
 
    End Subroutine Load_Radial_Profile



    Subroutine Calculate_Conductive_Profile
        Implicit None
        Real*8 :: amp, ftest
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

        t_init_file = '__nothing__'
        w_init_file = '__nothing__'
        p_init_file = '__nothing__'
        z_init_file = '__nothing__'
        c_init_file = '__nothing__'
        a_init_file = '__nothing__'

    End Subroutine Restore_InitialCondition_Defaults
    Subroutine Convert_AB_Pair(abbuf)
        ! Faithful slot -> q+/- conversion of the horizontal pair's AB history.
        ! abbuf is in p1 layout (lm-distributed); route the pair through the
        ! transpose cycle to rlm space (all-l-local), convert with the
        ! certified machinery, and return via the standard s2b -> p1b leg.
        Implicit None
        Real*8, Intent(InOut) :: abbuf(:,:,:,:)
        Type(SphericalBuffer) :: abtmp
        Integer :: abcount(3,2), mp

        abcount(:,:) = 2
        Call abtmp%init(field_count = abcount, config = 'p1a')
        Call abtmp%construct('p1a')
        abtmp%p1a(:,:,:,1) = abbuf(:,:,:,vteq)
        abtmp%p1a(:,:,:,2) = abbuf(:,:,:,vpeq)
        Call abtmp%reform()          ! p1a -> s2a (pure transpose)

        Call Convert_Slot_Pair_To_Spin(abtmp%s2a, 1, 2)

        ! Return leg: s2a and s2b share the rlm layout; construct s2b, copy,
        ! and take the standard transpose back to p1b.
        Call abtmp%construct('s2b')
        Do mp = my_mp%min, my_mp%max
            abtmp%s2b(mp)%data(:,:,:,1:2) = abtmp%s2a(mp)%data(:,:,:,1:2)
        Enddo
        Call abtmp%deconstruct('s2a')
        abtmp%config = 's2b'
        Call abtmp%reform()          ! s2b -> p1b (pure transpose)

        abbuf(:,:,:,vteq) = abtmp%p1b(:,:,:,1)
        abbuf(:,:,:,vpeq) = abtmp%p1b(:,:,:,2)
        Call abtmp%deconstruct('p1b')
    End Subroutine Convert_AB_Pair

End Module Initial_Conditions
