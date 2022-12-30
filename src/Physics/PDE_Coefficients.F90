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

! This module defines the constant and nonconstant coefficients
! that appear in the system of PDEs solved by Rayleigh.
! These coefficients are broadly classified into two groups:
! 1.  The reference (or thermodynamic background) state
! 2.  Transport coefficients (nu, kappa, eta).

! This module is divided into four sections.  Search on any of the
! phrases below to jump to that section.
! I.  Variables describing the background reference state
! II.  Variables Related to the Transport Coefficients
! III.  Subroutines used to define the background reference state
! IV.  Subroutines used to define the transport coefficients

Module PDE_Coefficients
    Use ProblemSize
    Use Controls
    Use Math_Constants
    Use Math_Utility
    Use General_MPI, Only : BCAST2D
    Implicit None

    !///////////////////////////////////////////////////////////
    ! I.  Variables describing the background reference state

    ! The following object (handle for which is "ref") is used by the code to access the PDE coefficients
    Type ReferenceInfo
        Real*8, Allocatable :: Density(:)
        Real*8, Allocatable :: dlnrho(:)
        Real*8, Allocatable :: d2lnrho(:)

        Real*8, Allocatable :: Temperature(:)
        Real*8, Allocatable :: dlnT(:)

        Real*8, Allocatable :: dsdr(:)

        Real*8, Allocatable :: heating(:)

        Real*8 :: Coriolis_Coeff ! Multiplies z_hat x u in momentum eq.
        Real*8 :: Lorentz_Coeff ! Multiplies (Del X B) X B in momentum eq.
        Real*8, Allocatable :: Buoyancy_Coeff(:) ! Multiplies {S,T} in momentum eq. ..typically = gravity/cp
        Real*8, Allocatable :: chi_buoyancy_coeff(:,:) ! Multiplies Chis in momentum eq.
        Real*8, Allocatable :: dpdr_w_term(:) ! Multiplies d_by_dr{P/rho} in momentum eq.
        Real*8, Allocatable :: pressure_dwdr_term(:) ! Multiplies l(l+1)/r^2 (P/rho) in Div dot momentum eq.

        ! The following two terms are used to compute the ohmic and viscous heating
        Real*8, Allocatable :: ohmic_amp(:) ! Multiplied by {eta(r),H(r)}J^2 in dSdt eq.
        Real*8, Allocatable :: viscous_amp(:) ! Multiplied by {nu(r),N(r)}{e_ij terms) in dSdt eq.

    End Type ReferenceInfo

    Type(ReferenceInfo) :: ref
    ! Allow up to 50 active/passive scalar fields
    Integer, Parameter :: n_scalar_max = 50

    ! Version number for the "equation_coefficients" container that is output to the simulation directory
    ! (i.e., the human-obtainable version of "ref")
    Integer, Parameter  :: eqn_coeff_version = 1

    ! Which background state to use; default 1 (non-dimensional Boussinesq)
    Integer :: reference_type = 1 

    ! Nondimensional variables (reference_type = 1,3)
    Real*8 :: Rayleigh_Number         = 1.0d0
    Real*8 :: Ekman_Number            = 1.0d0
    Real*8 :: Prandtl_Number          = 1.0d0
    Real*8 :: Magnetic_Prandtl_Number = 1.0d0
    Real*8 :: gravity_power           = 0.0d0
    Real*8 :: Dissipation_Number      = 0.0d0
    Real*8 :: Modified_Rayleigh_Number = 0.0d0

    ! Nondimensional variables for the active/passive scalar fields
    Real*8 :: chi_a_rayleigh_number(1:n_scalar_max)          = 0.0d0
    Real*8 :: chi_a_prandtl_number(1:n_scalar_max)           = 1.0d0
    Real*8 :: chi_a_modified_rayleigh_number(1:n_scalar_max) = 0.0d0
    Real*8 :: chi_p_prandtl_number(1:n_scalar_max)           = 1.0d0

    ! Dimensional anelastic variables (reference_type = 2)
    Real*8 :: pressure_specific_heat  = 1.0d0 ! CP (not CV)
    Real*8 :: poly_n = 0.0d0 ! Polytropic index
    Real*8 :: poly_Nrho = 0.0d0 ! Number of density scale heights across domain
    Real*8 :: poly_mass = 0.0d0 ! Stellar mass; g(r) = G*poly_mass/r^2
    Real*8 :: poly_rho_i = 0.0d0 ! Density (g/cm^3) at the inner radius rmin
    Real*8 :: Angular_Velocity = -1.0d0 ! Frame rotation rate (sets Coriolis force)

    ! Custom reference-state variables (reference_type = 4)
    Integer, Parameter   :: max_ra_constants = 10 + 2*n_scalar_max
    Integer, Parameter   :: max_ra_functions = 14 + 2*n_scalar_max
    Integer              :: n_ra_constants
    Integer              :: n_ra_functions
    Logical              :: with_custom_reference = .false.
    Logical              :: override_constants = .false.
    Logical              :: override_constant(1:max_ra_constants) = .false. ! in namelist
    Integer              :: with_custom_constants(1:max_ra_constants) = 0   ! in namelist
    Integer              :: with_custom_functions(1:max_ra_functions) = 0   ! in namelist
    Real*8               :: ra_constants(1:max_ra_constants) = 0.0d0        ! in namelist
    Integer, Allocatable :: ra_constant_set(:)
    Integer, Allocatable :: ra_function_set(:)
    Logical, Allocatable :: use_custom_constant(:)
    Logical, Allocatable :: use_custom_function(:)
    Real*8, Allocatable  :: ra_functions(:,:)
    Logical              :: custom_reference_read = .false.
    Character*120        :: custom_reference_file ='nothing'    

    ! Internal heating variables
    Integer :: heating_type = 0 ! 0 means no reference heating.  > 0 selects optional reference heating
    Real*8  :: Luminosity = 0.0d0 ! Specifies the integral of the heating function
    Real*8  :: Heating_Integral = 0.0d0  ! Same as luminosity (for non-star watchers)
    Real*8  :: Heating_EPS = 1.0d-12  ! Small number to test whether luminosity was specified
    Logical :: adjust_reference_heating = .false. ! Flag used to decide if luminosity determined via boundary conditions
    Real*8, Allocatable :: s_conductive(:)
    
    ! Minimum time step based on rotation rate
    ! (determined by the reference state)
    Real*8 :: max_dt_rotation = 0.0d0

    ! Alter some of the above (I) by reading the main_input file
    Namelist /Reference_Namelist/ reference_type,poly_n, poly_Nrho, poly_mass, poly_rho_i, &
            & pressure_specific_heat, heating_type, luminosity, Angular_Velocity, &
            & Rayleigh_Number, Ekman_Number, Prandtl_Number, Magnetic_Prandtl_Number, &
            & gravity_power, custom_reference_file, &
            & Dissipation_Number, Modified_Rayleigh_Number, &
            & Heating_Integral, override_constants, override_constant, ra_constants, with_custom_constants, &
            & with_custom_functions, with_custom_reference, &
            & chi_a_rayleigh_number, chi_a_prandtl_number, &
            & chi_a_modified_rayleigh_number, chi_p_prandtl_number

    !///////////////////////////////////////////////////////////////////////////////////////
    ! II.  Variables Related to the Transport Coefficients

    Real*8, Allocatable :: nu(:), kappa(:), eta(:)
    Real*8, Allocatable :: dlnu(:), dlnkappa(:), dlneta(:)
    real*8, allocatable :: kappa_chi_a(:,:), kappa_chi_p(:,:)
    real*8, allocatable :: dlnkappa_chi_a(:,:), dlnkappa_chi_p(:,:)

    Real*8, Allocatable :: ohmic_heating_coeff(:)
    Real*8, Allocatable :: viscous_heating_coeff(:)

    Real*8, Allocatable :: W_Diffusion_Coefs_0(:), W_Diffusion_Coefs_1(:)
    Real*8, Allocatable :: dW_Diffusion_Coefs_0(:), dW_Diffusion_Coefs_1(:), dW_Diffusion_Coefs_2(:)
    Real*8, Allocatable :: S_Diffusion_Coefs_1(:), Z_Diffusion_Coefs_0(:), Z_Diffusion_Coefs_1(:)
    Real*8, Allocatable :: A_Diffusion_Coefs_1(:)
    real*8, allocatable :: chi_a_diffusion_coefs_1(:,:), chi_p_diffusion_coefs_1(:,:)

    Integer :: kappa_type = 1, nu_type = 1, eta_type = 1
    Real*8  :: nu_top = -1.0d0, kappa_top = -1.0d0,  eta_top = -1.0d0
    Real*8  :: nu_power = 0, eta_power = 0, kappa_power = 0
    Integer :: kappa_chi_a_type(1:n_scalar_max) = 1
    Real*8  :: kappa_chi_a_top(1:n_scalar_max) = -1.0d0
    Real*8  :: kappa_chi_a_power(1:n_scalar_max) = 0
    Integer :: kappa_chi_p_type(1:n_scalar_max) = 1
    Real*8  :: kappa_chi_p_top(1:n_scalar_max) = -1.0d0
    Real*8  :: kappa_chi_p_power(1:n_scalar_max) = 0

    Logical :: hyperdiffusion = .false.
    Real*8  :: hyperdiffusion_beta = 0.0d0
    Real*8  :: hyperdiffusion_alpha = 1.0d0

    ! Alter some of the above (II) by reading the main_input file
    Namelist /Transport_Namelist/ nu_type, kappa_type, eta_type, &
            & nu_power, kappa_power, eta_power, &
            & nu_top, kappa_top, eta_top, &
            & hyperdiffusion, hyperdiffusion_beta, hyperdiffusion_alpha, &
            & kappa_chi_a_type, kappa_chi_a_top, kappa_chi_a_power, &
            & kappa_chi_p_type, kappa_chi_p_top, kappa_chi_p_power

Contains

    !/////////////////////////////////////////////////////////////////
    ! III.  Subroutines used to define the background reference state

    Subroutine Initialize_Reference()
        Implicit None

        If (my_rank .eq. 0) Then
            Call stdout%print(" -- Initalizing Reference State...")
            Call stdout%print(" ---- Specified parameters:")
        Endif

        Call Allocate_Reference_State()

        If (reference_type .eq. 1) Then
            Call Constant_Reference()
        Endif

        If (reference_type .eq. 2) Then
            Call Polytropic_Reference()
        Endif

        If (reference_type .eq. 3) Then
            Call Polytropic_ReferenceND()
        Endif

        If (reference_type .eq. 4) Then
            Call Get_Custom_Reference()
        Endif

        If (my_rank .eq. 0) Then
            Call stdout%print(" -- Reference State initialized.")
            Call stdout%print(" ")
        Endif

        If (with_custom_reference) Call Augment_Reference()  

        If (rotation) Call Set_Rotation_dt()

    End Subroutine Initialize_Reference

    Subroutine Allocate_Reference_State
        Implicit None

        n_ra_constants = 10 + 2*(n_active_scalars + n_passive_scalars)
        n_ra_functions = 14 + 2*(n_active_scalars + n_passive_scalars)

        Allocate(ref%density(1:N_R))
        Allocate(ref%temperature(1:N_R))
        Allocate(ref%dlnrho(1:N_R))
        Allocate(ref%d2lnrho(1:N_R))
        Allocate(ref%dlnt(1:N_R))
        Allocate(ref%dsdr(1:N_R))
        Allocate(ref%Buoyancy_Coeff(1:N_R))
        Allocate(ref%dpdr_w_term(1:N_R))
        Allocate(ref%pressure_dwdr_term(1:N_R))
        Allocate(ref%ohmic_amp(1:N_R))
        Allocate(ref%viscous_amp(1:N_R))
        Allocate(ref%heating(1:N_R))
        Allocate(ref%chi_buoyancy_coeff(n_active_scalars,1:N_R))

        Allocate(ra_constant_set(1:n_ra_constants))
        ra_constant_set = 0
        Allocate(ra_function_set(1:n_ra_functions))
        ra_function_set = 0
        Allocate(use_custom_constant(1:n_ra_constants))
        use_custom_constant = .false.
        Allocate(use_custom_function(1:n_ra_functions))
        use_custom_function = .false.
        Allocate(ra_functions(1:N_R, 1:n_ra_functions))
        ra_functions(:,:) = Zero

        ref%density(:)            = Zero
        ref%temperature(:)        = Zero
        ref%dlnrho(:)             = Zero
        ref%d2lnrho(:)            = Zero
        ref%dlnt(:)               = Zero
        ref%dsdr(:)               = Zero
        ref%buoyancy_coeff(:)     = Zero
        ref%dpdr_w_term(:)        = Zero
        ref%pressure_dwdr_term(:) = Zero
        ref%ohmic_amp(:)          = Zero
        ref%viscous_amp(:)        = Zero
        ref%heating(:)            = Zero
        ref%chi_buoyancy_coeff(:,:) = Zero

        ref%Coriolis_Coeff = Zero
        ref%Lorentz_Coeff  = Zero

    End Subroutine Allocate_Reference_State

    Subroutine Set_Rotation_dt()
        Implicit None
        Real*8 :: rotational_timescale
        !Adjust the maximum timestep to account for rotation rate, if necessary.
        
        rotational_timescale = 1.0d0/ra_constants(1)
        
        ! Minimum sampling would require two time samples per rotational timescale.
        ! We specify 4 samples and further adjust by the CFL safety factor.
        max_dt_rotation = rotational_timescale*0.25d0*cflmax  
        
        ! We only modify the timestep if not running a benchmark.
        ! Those models require set timesteps to satisfy the automated
        ! benchmark test during continuous integration.
        If (benchmark_mode .eq. 0) max_time_step = Min(max_time_step,max_dt_rotation)
    
    End Subroutine Set_Rotation_dt

    Subroutine Constant_Reference()
        Implicit None
        Integer :: i,j
        Real*8 :: r_outer, r_inner, prefactor, amp, pscaling
        Character*12 :: dstring
        Character*8 :: dofmt = '(ES12.5)'

        viscous_heating = .false.  ! Turn this off for Boussinesq runs
        ohmic_heating = .false.
        If (my_rank .eq. 0) Then
            Call stdout%print(" ---- Reference type           : "//trim(" Boussinesq (Non-dimensional)"))
            Write(dstring,dofmt)Rayleigh_Number
            Call stdout%print(" ---- Rayleigh Number          : "//trim(dstring))
            Write(dstring,dofmt)Ekman_Number
            Call stdout%print(" ---- Ekman Number             : "//trim(dstring))
            Write(dstring,dofmt)Prandtl_Number
            Call stdout%print(" ---- Prandtl Number           : "//trim(dstring))
            If (magnetism) Then
                Write(dstring,dofmt)Magnetic_Prandtl_Number
                Call stdout%print(" ---- Magnetic Prandtl Number  : "//trim(dstring))
            Endif
        Endif

        ! Set the equation coefficients (except mostly the ones having to do with diffusivities)
        ! Do this by constants first, then functions (more or less by ascending index)

        ! Constants
        ra_constants(1) = 2.0d0/Ekman_Number
        ra_constants(2) = Rayleigh_Number/Prandtl_Number 
        If (rotation) Then
            ra_constants(3) = 1.0d0/Ekman_Number            
        Else
            ra_constants(3) = 1.0d0            
        Endif

        If (magnetism) Then
            ra_constants(4) = 1.0d0/(Magnetic_Prandtl_Number*Ekman_Number)
            eta_top     = 1.0d0/Magnetic_Prandtl_Number
        Else
            ra_constants(4) = 0.0d0
            eta_top     = 0.0d0
        Endif

        nu_top = 1.0d0
        kappa_top = 1.0d0/Prandtl_Number
        Do i = 1, n_active_scalars
            kappa_chi_a_top(i)   = 1.0d0/Chi_A_Prandtl_Number(i)
        Enddo
        Do i = 1, n_passive_scalars
            kappa_chi_p_top(i)   = 1.0d0/Chi_P_Prandtl_Number(i)
        Enddo
 
        ! Though it has no effect (viscous_heating = ohmic_heating = .false.), 
        ! turn off the viscous/ohmic dissipation terms (again) for logical consistency
        ra_constants(8) = 0.0d0
        ra_constants(9) = 0.0d0

        ! Functions
        ra_functions(:,1) = 1.0d0
        ra_functions(:,8) = 0.0d0
        ra_functions(:,9) = 0.0d0 

        ra_functions(:,2) = (radius(:)/radius(1))**gravity_power

        ra_functions(:,4) = 1.0d0
        ra_functions(:,10) = 0.0d0      

        ra_functions(:,14) = 0.0d0      

        ! Heating
        Call Initialize_Reference_Heating()
        If (heating_type .eq. 0) Then
            !Otherwise, s_conductive will be calculated in initial_conditions
            Allocate(s_conductive(1:N_R))
            r_outer = radius(1)
            r_inner = radius(N_R)
            prefactor = r_outer*r_inner/(r_inner-r_outer)
            Do i = 1, N_R
                s_conductive(i) = prefactor*(1.0d0/r_outer-1.0d0/radius(i))
            Enddo
        Endif

    End Subroutine Constant_Reference

    Subroutine Polytropic_ReferenceND()
        Implicit None
        Integer :: i
        Real*8 :: dtmp, otmp
        Real*8, Allocatable :: dtmparr(:), gravity(:)
        Character*12 :: dstring
        Character*8 :: dofmt = '(ES12.5)'

        If (my_rank .eq. 0) Then
            Call stdout%print(" ---- Reference type           : "//trim(" Polytrope (Non-dimensional)"))
            Write(dstring,dofmt)Modified_Rayleigh_Number
            Call stdout%print(" ---- Modified Rayleigh Number : "//trim(dstring))
            Write(dstring,dofmt)Ekman_Number
            Call stdout%print(" ---- Ekman Number             : "//trim(dstring))
            Write(dstring,dofmt)Prandtl_Number
            Call stdout%print(" ---- Prandtl Number           : "//trim(dstring))
            If (magnetism) Then
                Write(dstring,dofmt)Magnetic_Prandtl_Number
                Call stdout%print(" ---- Magnetic Prandtl Number  : "//trim(dstring))
            Endif
            Write(dstring,dofmt)poly_n
            Call stdout%print(" ---- Polytropic Index         : "//trim(dstring))
            Write(dstring,dofmt)poly_nrho
            Call stdout%print(" ---- Density Scaleheights     : "//trim(dstring))
        Endif

        If (aspect_ratio .lt. 0) Then
            aspect_ratio = rmax/rmin
        Endif
        Allocate(dtmparr(1:N_R), gravity(1:N_R))
        dtmparr(:) = 0.0d0

        Dissipation_Number = aspect_ratio*(exp(poly_Nrho/poly_n)-1.0D0)
        dtmp = 1.0D0/(1.0D0-aspect_ratio)

        ! Set the equation coefficients (except mostly the ones having to do with diffusivities)
        ! Do this by constants first, then functions (more or less by ascending index)

        ! Constants
        ra_constants(1) = 2.0d0
        ra_constants(2) = Modified_Rayleigh_Number 
        ra_constants(3) = 1.0d0

        nu_top = Ekman_Number
        kappa_top = Ekman_Number/Prandtl_Number
        If (magnetism) Then
            ra_constants(4) = Ekman_Number/Magnetic_Prandtl_Number
            eta_top = Ekman_Number/Magnetic_Prandtl_Number
        Else
            ra_constants(4) = 0.0d0
            eta_top     = 0.0d0
        Endif

        Do i = 1, n_active_scalars
            kappa_chi_a_top(i)   = Ekman_Number/Chi_A_Prandtl_Number(i)
        Enddo
        Do i = 1, n_passive_scalars
            kappa_chi_p_top(i)   = Ekman_Number/Chi_P_Prandtl_Number(i)
        Enddo
 
        ra_constants(8) = Dissipation_Number*Ekman_Number/Modified_Rayleigh_Number
        If (magnetism) Then
            ra_constants(9) = Dissipation_Number*Ekman_Number**2/ & 
                (Modified_Rayleigh_Number*Magnetic_Prandtl_Number**2)
        Endif ! if not magnetism, ra_constants(9) was initialized to zero

        ! Functions
        gravity = (rmax**2)*OneOverRSquared(:)
        ra_functions(:,4) = dtmp*Dissipation_Number*(dtmp*One_Over_R(:)-1.0D0)+1.0D0
        ra_functions(:,1) = ra_functions(:,4)**poly_n
        ra_functions(:,2) = ra_functions(:,1)*gravity

        !Compute the background temperature gradient : dTdr = -Dg,  d2Tdr2 = 2*D*g/r (for g ~1/r^2)
        dtmparr = -Dissipation_Number*gravity
        !Now, the logarithmic derivative of temperature
        ra_functions(:,10) = dtmparr/ra_functions(:,4)

        !And then logarithmic derivative of rho : dlnrho = n dlnT
        ra_functions(:,8) = poly_n*ra_functions(:,10)

        !Now, the second logarithmic derivative of rho :  d2lnrho = (n/T)*d2Tdr2 - n*(dlnT^2)
        ra_functions(:,9) = -poly_n*(ra_functions(:,10)**2)
        dtmparr = (poly_n/ra_functions(:,4))*(2.0d0*Dissipation_Number*gravity/radius) ! (n/T)*d2Tdr2
        ra_functions(:,9) = ra_functions(:,9)+dtmparr

        ! Loren, 12/22/2022: This is probably a bug...are stable layers disallowed under non-D anelastic?
        ra_functions(:,14) = 0.0d0

        ! Heating
        Call Initialize_Reference_Heating()

        DeAllocate(dtmparr, gravity)
    End Subroutine Polytropic_ReferenceND

    Subroutine Polytropic_Reference()
        Implicit None
        Integer :: i
        Real*8 :: zeta_0,  c0, c1, d
        Real*8 :: rho_c, P_c, T_c,denom
        Real*8 :: thermo_gamma, volume_specific_heat
        Real*8 :: beta
        Real*8 :: Gravitational_Constant = 6.67d-8 ! cgs units
        Real*8, Allocatable :: zeta(:), gravity(:)
        Real*8 :: One
        Real*8 :: InnerRadius, OuterRadius
        Character*12 :: dstring
        Character*8 :: dofmt = '(ES12.5)'
        If (my_rank .eq. 0) Then
            Call stdout%print(" ---- Reference type                : "//trim(" Polytrope (Dimensional)"))
            Write(dstring,dofmt)Angular_Velocity
            Call stdout%print(" ---- Angular Velocity (rad/s)      : "//trim(dstring))
            Write(dstring,dofmt)poly_rho_i
            Call stdout%print(" ---- Inner-Radius Density (g/cm^3) : "//trim(dstring))
            Write(dstring,dofmt)poly_mass
            Call stdout%print(" ---- Interior Mass  (g)            : "//trim(dstring))
            Write(dstring,dofmt)poly_n
            Call stdout%print(" ---- Polytropic Index              : "//trim(dstring))
            Write(dstring,dofmt)poly_nrho
            Call stdout%print(" ---- Density Scaleheights          : "//trim(dstring))
            Write(dstring,dofmt)pressure_specific_heat
            Call stdout%print(" ---- CP (erg g^-1 cm^-3 K^-1)      : "//trim(dstring))
        Endif

        ! Adiabatic, Polytropic Reference State (see, e.g., Jones et al. 2011)
        ! The following parameters are read from the input file.
        ! poly_n
        ! poly_Nrho
        ! poly_mass
        ! poly_rho_i
        ! Note that cp must also be specified.

        ! Set the equation coefficients (except mostly the ones having to do with diffusivities)
        ! Do this by constants first, then functions (more or less by ascending index)

        ! Constants
        ra_constants(1) = 2.0d0*Angular_Velocity
        ra_constants(2) = 1.0d0 
        ra_constants(3) = 1.0d0

        If (magnetism) Then
            ra_constants(4) = 1.0d0/four_pi
        Else
            ra_constants(4) = 0.0d0
        Endif
 
        ra_constants(8) = 1.0d0
        If (magnetism) Then
            ra_constants(9) = 1.0d0/four_pi
        Endif ! if not magnetism, ra_constants(9) was initialized to zero
        
        ! Functions

        ! Compute dimensional polytrope
        InnerRadius = Radius(N_r)
        OuterRadius = Radius(1)

        One = 1.0d0
        !-----------------------------------------------------------
        beta = InnerRadius/OuterRadius

        denom = beta * exp(poly_Nrho / poly_n) + 1.d0
        zeta_0 = (beta+1.d0)/denom

        c0 = (2.d0 * zeta_0 - beta - 1.d0) / (1.d0 - beta)

        denom = (1.d0 - beta)**2
        c1 = (1.d0+beta)*(1.d0-zeta_0)/denom

        !-----------------------------------------------------------
        ! allocate and define zeta
        ! also rho_c, T_c, P_c

        Allocate(zeta(N_R), gravity(1:N_R))

        d = OuterRadius - InnerRadius

        zeta = c0 + c1 * d / Radius

        rho_c = poly_rho_i / zeta(N_R)**poly_n

        denom = (poly_n+1.d0) * d * c1
        P_c = Gravitational_Constant * poly_mass * rho_c / denom

        T_c = (poly_n+1.d0) * P_c / (Pressure_Specific_Heat * rho_c)

        !-----------------------------------------------------------
        ! Initialize reference structure
        gravity = Gravitational_Constant * poly_mass / Radius**2

        ! The following is needed to calculate the entropy gradient
        thermo_gamma = 5.0d0/3.0d0
        volume_specific_heat = pressure_specific_heat / thermo_gamma

        ra_functions(:,1) = rho_c * zeta**poly_n
        ra_functions(:,8) = - poly_n * c1 * d / (zeta * Radius**2)
        ra_functions(:,9) = - ra_functions(:,8)*(2.0d0/Radius-c1*d/zeta/Radius**2)

        ra_functions(:,2) = ra_functions(:,1)*gravity/Pressure_Specific_Heat

        ra_functions(:,4) = T_c * zeta
        ra_functions(:,10) = -(c1*d/Radius**2)/zeta

        ra_functions(:,14) = volume_specific_heat * (ra_functions(:,10) - (thermo_gamma - 1.0d0) * ra_functions(:,8))

        Deallocate(zeta, gravity)

        ! Heating
        Call Initialize_Reference_Heating()     

    End Subroutine Polytropic_Reference

     Subroutine Initialize_Reference_Heating()
        Implicit None
        ! This is where a volumetric heating function Phi(r) = Q(r)/(rho(r)*T(r)) is computed
        ! This function appears in the entropy equation as
        ! dSdt = Phi(r)
        ! Phi(r) may represent internal heating of any type.  For stars, this heating would be
        ! associated with radiative diffusion of the background state and/or nuclear burning.

        If (heating_type .eq. 1) Then
            Call Constant_Entropy_Heating()
        Endif

        If (heating_type .eq. 4) Then
            Call  Constant_Energy_Heating()
        Endif

        ! Volumetric heating (Q ~ f_6) is normalized so that:
        ! Int_volume rho T f_6 dV = 1 

        !If luminosity or heating_integral has been set,
        !we set the integral of f_6 (i.e., c_10) using that

        !Otherwise, we will use the boundary thermal flux
        !To establish the integral
        If (heating_type .gt. 0) Then
            If (abs(luminosity) .gt. heating_eps) Then
                ra_constants(10) = luminosity
                adjust_reference_heating = .false.                

            Elseif (abs(Heating_Integral) .gt. heating_eps) Then
                ra_constants(10) = heating_integral
                adjust_reference_heating = .false.
            Else
                ! c_10 will be set later in BoundaryConditions.F90 (and f_6 rescaled)
                adjust_reference_heating = .true.
            Endif
        Endif
    End Subroutine Initialize_Reference_Heating

    Subroutine Augment_Reference()
        Implicit None
        Real*8, Allocatable :: temp_functions(:,:), temp_constants(:)

        If (my_rank .eq. 0) Then
            Call stdout%print('Reference state will be augmented.')
            Call stdout%print('Only heating and buoyancy may be modified.')
            Call stdout%print('Heating requires both c_10 and f_6 to be set.')
            Call stdout%print('Buoyancy requires both c_2 and f_2 to be set.')
            Call stdout%print('Reading from: '//Trim(custom_reference_file))
        Endif

        ! Before reading the custom reference file, guard against
        ! potential overwrites of all ra_functions and ra_constants which,
        ! in this case were determined by ref_types 1,2,3.  Only certain
        ! combinations are allowed to be overwritten.

        Allocate(temp_functions(1:n_r, 1:n_ra_functions))
        Allocate(temp_constants(1:n_ra_constants))
        temp_functions(:,:) = ra_functions(:,:)

        ! Note that ra_constants is allocated up to max_ra_constants
        temp_constants(:) = ra_constants(1:n_ra_constants)


        Call Read_Custom_Reference_File(custom_reference_file)

        If (use_custom_constant(10) .and. use_custom_function(6)) Then
            If (my_rank .eq. 0) Then
                Call stdout%print('Heating has been set to:')
                Call stdout%print('f_6*c_10')
                Call stdout%print(' ')
            Endif
            ref%heating(:) = ra_functions(:,6)/(ref%density*ref%temperature)*ra_constants(10)
            temp_functions(:,6) = ra_functions(:,6)
            temp_constants(10)  = ra_constants(10)
        Endif

        If (use_custom_constant(2) .and. use_custom_function(2)) Then
            If (my_rank .eq. 0) Then
                Call stdout%print('Buoyancy_coeff has been set to:')
                Call stdout%print('f_2*c_2')
                Call stdout%print(' ')
            Endif
            ref%buoyancy_coeff(:) = ra_constants(2)*ra_functions(:,2)
            temp_functions(:,2) = ra_functions(:,2)
            temp_constants(2) = ra_constants(2)
        Endif

        ra_constants(1:n_ra_constants) = temp_constants(:)
        ra_functions(:,:) = temp_functions(:,:)
        DeAllocate(temp_functions, temp_constants)

    End Subroutine Augment_Reference

    Subroutine Constant_Energy_Heating()
        Implicit None
        ! rho T dSdt = alpha : energy deposition per unit volume (Q ~ f_6) is constant
        ! dSdt = alpha/rhoT : entropy deposition per unit volume is ~ 1/Pressure
        ra_functions(:,6) = 1.0d0/shell_volume
    End Subroutine Constant_Energy_Heating

    Subroutine Constant_Entropy_Heating()
        Implicit None
        Real*8 :: integral, alpha

        ! dSdt = alpha : Entropy deposition per unit volume is constant
        ! rho T dSdt = rho T alpha : Energy deposition per unit volume (Q ~ f_6) is ~ Pressure

        ! Luminosity (c_10) is specified as an input
        ra_functions(:,6) = ra_functions(:,1)*ra_functions(:,4)
        Call Integrate_in_radius(ra_functions(:,6),integral) !Int_rmin_rmax rho T r^2 dr
        integral = integral*4.0d0*pi  ! Int_V temp dV

        ra_functions(:,6)= ra_functions(:,6)/integral
    End Subroutine Constant_Entropy_Heating

    Subroutine Integrate_in_radius(func,int_func)
        Implicit None
        Real*8, Intent(In) :: func(1:)
        Real*8, Intent(Out) :: int_func
        Integer :: i
        Real*8 :: rcube
        !compute integrate_r=rmin_r=rmax func*r^2 dr
        rcube = (rmax**3-rmin**3)*One_Third
        int_func = 0.0d0
        Do i = 1, n_r
            int_func = int_func+func(i)*radial_integral_weights(i)
        Enddo
        int_func = int_func*rcube

    End Subroutine Integrate_in_radius

    Subroutine Get_Custom_Reference()
        Implicit None
        Integer :: i, fi
        Character(len=2) :: ind
        Integer :: fi_to_check(4) = (/1, 2, 4, 6/)

        If (my_rank .eq. 0) Then
            Write(6,*)'Custom reference state specified.'
            Write(6,*)'Reading from: ', custom_reference_file
        Endif

        ! All we have to do for this type is read the custom_reference_file
        Call Read_Custom_Reference_File(custom_reference_file)

        ! Consistency check
        ! Issue: if it's an "error", the code should exit below; 
        ! if not, then ERROR --> WARNING
        Do i=1,4
            fi = fi_to_check(i)
            If (ra_function_set(fi) .eq. 0) Then
                If (my_rank .eq. 0) Then
                    Write(ind, '(I2)') fi
                    Call stdout%print('ERROR: function f_'//Adjustl(ind)//' must be set in the custom reference file')
                Endif
            Endif
        Enddo

    End Subroutine Get_Custom_Reference

    Subroutine Write_Equation_Coefficients_File(filename)
        Character*120, Intent(In), Optional :: filename
        Character*120 :: ref_file
        Integer :: i, k, pi_integer
        pi_integer = 314

        If (present(filename)) Then
            ref_file = Trim(my_path)//filename
        Else
            ref_file = Trim(my_path)//'equation_coefficients'
        Endif

        If (my_rank .eq. 0) Then
            Open(unit=15, file=ref_file, form='unformatted', status='replace', access='stream')

            Write(15) pi_integer
            Write(15) eqn_coeff_version
            ra_constant_set(:) = 1
            ra_function_set(:) = 1
            Write(15) (ra_constant_set(i), i=1,n_ra_constants)
            Write(15) (ra_function_set(i), i=1,n_ra_functions)
            Write(15) (ra_constants(i), i=1,n_ra_constants)
            Write(15) n_r
            Write(15) (radius(i), i=1,n_r)
            Do k=1, n_ra_functions
                Write(15) (ra_functions(i,k), i=1,n_r)
            Enddo

            Close(15)
        Endif

    End Subroutine Write_Equation_Coefficients_File

    Subroutine Read_Custom_Reference_File(filename)
        Character*120, Intent(In) :: filename
        Character*120 :: ref_file
        Integer :: pi_integer,nr_ref, eqversion
        Integer :: i, k, j, n_scalars
        Integer :: cset(1:n_ra_constants), fset(1:n_ra_functions)
        Real*8  :: input_constants(1:n_ra_constants)
        Real*8, Allocatable :: ref_arr_old(:,:), rtmp(:), rtmp2(:)
        Real*8, Allocatable :: old_radius(:)
        Character*12 :: dstring
        Character*8 :: dofmt = '(ES12.5)'
        Character(len=2) :: ind

        cset(:) = 0
        input_constants(:) = 0.0d0

        
        ref_file = Trim(my_path)//filename


        Open(unit=15,file=ref_file,form='unformatted', status='old',access='stream')

        !Verify Endianness
        Read(15)pi_integer
        If (pi_integer .ne. 314) Then
            close(15)
            Write(6,*)'Trying to convert.  pi = : ', pi_integer
            Open(unit=15,file=ref_file,form='unformatted', status='old', &
                 CONVERT = 'BIG_ENDIAN' , access='stream')
            Read(15)pi_integer
            If (pi_integer .ne. 314) Then
                Close(15)
                Write(6,*)'Trying to convert again.  pi is now: ', pi_integer
                Open(unit=15,file=ref_file,form='unformatted', status='old', &
                 CONVERT = 'LITTLE_ENDIAN' , access='stream')
                Read(15)pi_integer
                Write(6,*)'My final value of pi is: ', pi_integer
            Endif
        Endif

        If (pi_integer .eq. 314) Then

            ! Read in constants and their 'set' flags
            Read(15) eqversion
            Read(15) cset(1:n_ra_constants)
            Read(15) fset(1:n_ra_functions)
            Read(15) input_constants(1:n_ra_constants)
            
            ! Cset(i) is 1 if a constant(i) was set; it is 0 otherwise.
            ! The logic below deals with a constant set in both the reference
            ! file and in main_input.  Main_input values take precedence if
            ! override_constant(s) is set or if the reference file constant 
            ! was not set.
            Do i = 1, n_ra_constants
                If ( (.not. override_constants) .and. (.not. override_constant(i)) ) Then
                    ra_constants(i) = ra_constants(i) + cset(i)*(input_constants(i)-ra_constants(i))
                Endif
            Enddo

            ! determine which functions/constants were set by the user
            ra_function_set(:) = fset(:)
            Do i = 1, n_ra_constants
                If ((cset(i) .eq. 1) .or. override_constant(i) .or. override_constants) Then
                    ra_constant_set(i) = 1
                Endif
            Enddo

            ! Print the values of the constants
            Do k = 1, n_ra_constants
                If (my_rank .eq. 0) Then
                    Write(ind, '(I2)') k
                    Write(dstring,dofmt) ra_constants(k)
                    Call stdout%print('c_'//Adjustl(ind)//' = '//Trim(dstring))
                    !Write(6,*)'c: ', k, ra_constants(k)
                Endif
            Enddo

            ! Read the reference file's radial grid
            Read(15) nr_ref
            Allocate(ref_arr_old(1:nr_ref,1:n_ra_functions)) 
            Allocate(old_radius(1:nr_ref))

            Read(15)(old_radius(i),i=1,nr_ref)
            Do k = 1, n_ra_functions
                Read(15)(ref_arr_old(i,k) , i=1 , nr_ref)
            Enddo

            !Check to see if radius is tabulated in ascending or descending order.
            !If it is found to be in ascending order, reverse the radius and the 
            !input array of functions
            If (old_radius(1) .lt. old_radius(nr_ref)) Then

                If (my_rank .eq. 0) Write(6,*)'Reversing Radial Indices in Custom Ref File!'

                Allocate(rtmp(1:nr_ref))

                rtmp = old_radius
                Do i = 1, nr_ref
                    old_radius(i) = rtmp(nr_ref-i+1)
                Enddo

                Do k = 1, n_ra_functions
                    rtmp(:) = ref_arr_old(:,k)
                    Do i = 1, nr_ref
                        ref_arr_old(i,k) = rtmp(nr_ref-i+1)
                    Enddo
                Enddo

                DeAllocate(rtmp)

            Endif

            Close(15)
            custom_reference_read = .true.

            If (nr_ref .ne. n_r) Then
                !Interpolate onto the current radial grid if necessary
                !Note that the underlying assumption here is that same # of grid points
                ! means same grid - come back to this later for generality
                Allocate(   rtmp2(1:n_r))
                Allocate( rtmp(1:nr_ref))

                Do k = 1, n_ra_functions

                    rtmp(:) = ref_arr_old(:,k)
                    rtmp2(:) = 0.0d0

                    Call Spline_Interpolate(rtmp, old_radius, rtmp2, radius)

                    ra_functions(1:n_r,k) = rtmp2
                Enddo

                DeAllocate(rtmp,rtmp2)
            Else

                ! Bit redundant here, but may want to do filtering on ref_arr array
                ra_functions(1:n_r,1:n_ra_functions) = &
                       ref_arr_old(1:n_r,1:n_ra_functions)

                If (my_rank .eq. 0) Then
                    call stdout%print(" WARNING:  nr = nr_old.  Assuming grids are the same.")
                Endif
            Endif
            DeAllocate(ref_arr_old,old_radius)
            
            ! Finally, if the logarithmic derivatives of rho, T, nu, kappa, kappa_chi and eta were
            ! not specified, then we compute them here.
            ! only calculate the log derivative if the function was set, otherwise there
            ! are divide by zero issues
            If ((fset(8) .eq. 0) .and. (fset(1) .eq. 1)) Then
                Call log_deriv(ra_functions(:,1), ra_functions(:,8)) ! dlnrho
            Endif
            If ((fset(9) .eq. 0) .and. (fset(8) .eq. 1)) Then
                Call log_deriv(ra_functions(:,8), ra_functions(:,9), no_log=.true.) !d2lnrho
            Endif
            If ((fset(10) .eq. 0) .and. (fset(4) .eq. 1)) Then
                Call log_deriv(ra_functions(:,4), ra_functions(:,10)) !dlnT
            Endif
            If ((fset(11) .eq. 0) .and. (fset(3) .eq. 1)) Then
                Call log_deriv(ra_functions(:,3), ra_functions(:,11)) !dlnnu
            Endif
            If ((fset(12) .eq. 0) .and. (fset(5) .eq. 1)) Then
                Call log_deriv(ra_functions(:,5), ra_functions(:,12)) !dlnkappa
            Endif
            If ((fset(13) .eq. 0) .and. (fset(7) .eq. 1)) Then
                Call log_deriv(ra_functions(:,7), ra_functions(:,13)) !dlneta
            Endif
            n_scalars = n_active_scalars + n_passive_scalars
            do i = 0, (n_scalars - 1)
              If ((fset(16+i*2) .eq. 0) .and. (fset(15+i*2) .eq. 1)) Then
                  Call log_deriv(ra_functions(:,15+i*2), ra_functions(:,16+i*2)) !dlnkappa_chi
              Endif
            end do
        Else
            Write(6,*)'Error.  This file appears to be corrupt (check Endian convention).'
            Write(6,*)'Pi integer: ', pi_integer
        Endif

        ! only used if user wants to change reference_type=1,2,3
        If (with_custom_reference) Then
            Do i=1,n_ra_constants
                j = with_custom_constants(i)
                If ((j .gt. 0) .and. (j .le. n_ra_constants)) Then
                    If (ra_constant_set(j) .eq. 1) Then
                        use_custom_constant(j) = .true.
                    Else
                        If (my_rank .eq. 0) Then
                            Write(6,*)' '
                            Write(6,*)'You set with_custom_constant: ', j
                            Write(6,*)'But this constant was not set in either main_input or '
                            Write(6,*)'the custom reference file.  Selection will be ignored.'
                            Write(6,*)' '
                        Endif
                    Endif
                Endif
            Enddo

            Do i=1,n_ra_functions
                j = with_custom_functions(i)
                If ((j .gt. 0) .and. (j .le. n_ra_functions)) Then
                    use_custom_function(j) = .true.
                    If (ra_function_set(j) .eq. 1) Then
                        use_custom_function(j) = .true.
                    Else
                        If (my_rank .eq. 0) Then
                            Write(6,*)' '
                            Write(6,*)'You set with_custom_function: ', j
                            Write(6,*)'But this function was not set in either main_input or '
                            Write(6,*)'the custom reference file.  Selection will be ignored.'
                            Write(6,*)' '
                        Endif
                    Endif
                Endif
            Enddo
        Endif

    End Subroutine Read_Custom_Reference_File

    Subroutine Log_Deriv(arr1,arr2, no_log)
        Implicit None
        Real*8, Intent(In)    :: arr1(:)
        Real*8, Intent(InOut) :: arr2(:)
        Real*8, Allocatable   :: dtemp(:,:,:,:),dtemp2(:,:,:,:)
        Logical, Optional     :: no_log

        ! Computes logarithmic derivative of arr1 with respect to radius.
        ! Result is stored in arr2.
        ! Arr1 is assumed to be in physical space.
        ! Set no_log = .true. to take normal derivative.

        Allocate(dtemp(1:n_r,1,1,2))
        Allocate(dtemp2(1:n_r,1,1,2))

        dtemp(:,:,:,:) = 0.0d0
        dtemp2(:,:,:,:) = 0.0d0
        dtemp(1:n_r,1,1,1) = arr1(1:n_r)

        ! Transform to spectral space & de-alias
        Call gridcp%to_Spectral(dtemp,dtemp2)
        dtemp2((n_r*2)/3:n_r,1,1,1) = 0.0d0

        ! Take the derivative & de-alias
        Call gridcp%d_by_dr_cp(1,2,dtemp2,1)
        dtemp2((n_r*2)/3:n_r,1,1,2) = 0.0d0 

        !Transform back to physical space.
        Call gridcp%From_Spectral(dtemp2,dtemp)
        arr2(:) = dtemp(:,1,1,2)
        
        ! If desired, convert to logarithmic derivative (default)
        If (.not. present(no_log)) Then
            arr2(:) = arr2(:)/arr1(:)
        Endif

        DeAllocate(dtemp,dtemp2)

    End Subroutine log_deriv

    Subroutine Restore_Reference_Defaults
        Implicit None
        ! Restore all variables in this module to their default state.
        ! Deallocates all allocatable module variables.

        ! Which background state to use; default 1 (non-dimensional Boussinesq)
        reference_type = 1 

        ! Nondimensional variables (reference_type = 1,3)
        Rayleigh_Number         = 1.0d0
        Ekman_Number            = 1.0d0
        Prandtl_Number          = 1.0d0
        Magnetic_Prandtl_Number = 1.0d0
        gravity_power           = 0.0d0
        Dissipation_Number      = 0.0d0
        Modified_Rayleigh_Number = 0.0d0

        ! Nondimensional variables for the active/passive scalar fields
        chi_a_rayleigh_number(1:n_scalar_max)          = 0.0d0
        chi_a_prandtl_number(1:n_scalar_max)           = 1.0d0
        chi_a_modified_rayleigh_number(1:n_scalar_max) = 0.0d0
        chi_p_prandtl_number(1:n_scalar_max)           = 1.0d0

        ! Dimensional anelastic variables (reference_type = 2)
        pressure_specific_heat  = 1.0d0 ! CP (not CV)
        poly_n = 0.0d0 ! Polytropic index
        poly_Nrho = 0.0d0 ! Number of density scale heights across domain
        poly_mass = 0.0d0 ! Stellar mass; g(r) = G*poly_mass/r^2
        poly_rho_i = 0.0d0 ! Density (g/cm^3) at the inner radius rmin
        Angular_Velocity = -1.0d0 ! Frame rotation rate (sets Coriolis force)

        ! Custom reference-state variables (reference_type = 4)
        !       NOTE: n_ra_constants / n_ra_functions do not have default values (but maybe do, if you consider the Allocate_Reference_State() routine)
        with_custom_reference = .false.
        override_constants = .false.
        override_constant(1:max_ra_constants) = .false. ! in namelist
        with_custom_constants(1:max_ra_constants) = 0   ! in namelist
        with_custom_functions(1:max_ra_functions) = 0   ! in namelist
        ra_constants(1:max_ra_constants) = 0.0d0        ! in namelist
        custom_reference_read = .false.
        custom_reference_file ='nothing'  

        ! Internal heating variables
        heating_type = 0 ! 0 means no reference heating.  > 0 selects optional reference heating
        Luminosity = 0.0d0 ! Specifies the integral of the heating function
        Heating_Integral = 0.0d0  ! Same as luminosity (for non-star watchers)
        Heating_EPS = 1.0d-12  ! Small number to test whether luminosity was specified
        adjust_reference_heating = .false. ! Flag used to decide if luminosity determined via boundary conditions

        ! Minimum time step based on rotation rate
        ! (determined by the reference state)
        max_dt_rotation = 0.0d0

        ! Diffusion-coefficient (transport-namelist) variables
        kappa_type = 1
        nu_type = 1
        eta_type = 1
        nu_top = -1.0d0
        kappa_top = -1.0d0
        eta_top = -1.0d0
        nu_power = 0
        eta_power = 0
        kappa_power = 0

        kappa_chi_a_type(1:n_scalar_max) = 1
        kappa_chi_a_top(1:n_scalar_max) = -1.0d0
        kappa_chi_a_power(1:n_scalar_max) = 0
        kappa_chi_p_type(1:n_scalar_max) = 1
        kappa_chi_p_top(1:n_scalar_max) = -1.0d0
        kappa_chi_p_power(1:n_scalar_max) = 0

        hyperdiffusion = .false.
        hyperdiffusion_beta = 0.0d0
        hyperdiffusion_alpha = 1.0d0

        ! Deallocate "ref" object
        If (allocated(ref%Density)) DeAllocate(ref%density)
        If (allocated(ref%dlnrho)) DeAllocate(ref%dlnrho)
        If (allocated(ref%d2lnrho)) DeAllocate(ref%d2lnrho)

        If (allocated(ref%Temperature)) DeAllocate(ref%Temperature)
        If (allocated(ref%dlnT)) DeAllocate(ref%dlnT)

        If (allocated(ref%dsdr)) DeAllocate(ref%dsdr)

        If (allocated(ref%heating)) DeAllocate(ref%heating)

        !       NOTE: ref%Coriolis_Coeff and ref%Lorentz_Coeff have no default values (but maybe do, if you consider the Allocate_Reference_State() routine)
        If (allocated(ref%Buoyancy_Coeff)) DeAllocate(ref%Buoyancy_Coeff)
        If (allocated(ref%chi_buoyancy_coeff)) DeAllocate(ref%chi_buoyancy_coeff)
        If (allocated(ref%dpdr_w_term)) DeAllocate(ref%dpdr_w_term)
        If (allocated(ref%pressure_dwdr_term)) DeAllocate(ref%pressure_dwdr_term)

        If (allocated(ref%ohmic_amp)) DeAllocate(ref%ohmic_amp)
        If (allocated(ref%viscous_amp)) DeAllocate(ref%viscous_amp)

        ! Deallocate s_conductive
        If (allocated(s_conductive)) DeAllocate(s_conductive)
    
        ! Deallocate custom-reference stuff
        If (allocated(ra_constant_set)) DeAllocate(ra_constant_set)
        If (allocated(ra_function_set)) DeAllocate(ra_function_set)
        If (allocated(use_custom_constant)) DeAllocate(use_custom_constant)
        If (allocated(use_custom_function)) DeAllocate(use_custom_function)
        If (allocated(ra_functions)) DeAllocate(ra_functions)

        ! Deallocate transport-coefficient stuff
        If (allocated(nu)) DeAllocate(nu)
        If (allocated(kappa)) DeAllocate(kappa)
        If (allocated(eta)) DeAllocate(eta)
        If (allocated(dlnu)) DeAllocate(dlnu)
        If (allocated(dlnkappa)) DeAllocate(dlnkappa)
        If (allocated(dlneta)) DeAllocate(dlneta)
        If (allocated(kappa_chi_a)) DeAllocate(kappa_chi_a)
        If (allocated(kappa_chi_p)) DeAllocate(kappa_chi_p)
        If (allocated(dlnkappa_chi_a)) DeAllocate(dlnkappa_chi_a)
        If (allocated(dlnkappa_chi_p)) DeAllocate(dlnkappa_chi_p)

        If (allocated(ohmic_heating_coeff)) DeAllocate(ohmic_heating_coeff)
        If (allocated(viscous_heating_coeff)) DeAllocate(viscous_heating_coeff)

        If (allocated(W_Diffusion_Coefs_0)) DeAllocate(W_Diffusion_Coefs_0)
        If (allocated(W_Diffusion_Coefs_1)) DeAllocate(W_Diffusion_Coefs_1)
        If (allocated(dW_Diffusion_Coefs_0)) DeAllocate(dW_Diffusion_Coefs_0)
        If (allocated(dW_Diffusion_Coefs_1)) DeAllocate(dW_Diffusion_Coefs_1)
        If (allocated(dW_Diffusion_Coefs_2)) DeAllocate(dW_Diffusion_Coefs_2)
        If (allocated(S_Diffusion_Coefs_1)) DeAllocate(S_Diffusion_Coefs_1)
        If (allocated(Z_Diffusion_Coefs_0)) DeAllocate(Z_Diffusion_Coefs_0)
        If (allocated(Z_Diffusion_Coefs_1)) DeAllocate(Z_Diffusion_Coefs_1)
        If (allocated(A_Diffusion_Coefs_1)) DeAllocate(A_Diffusion_Coefs_1)
        If (allocated(chi_a_diffusion_coefs_1)) DeAllocate(chi_a_diffusion_coefs_1)
        If (allocated(chi_p_diffusion_coefs_1)) DeAllocate(chi_p_diffusion_coefs_1)

    End Subroutine Restore_Reference_Defaults


    !//////////////////////////////////////////////////////////////////////////
    ! IV.  Subroutines used to define the transport coefficients

    Subroutine Initialize_Transport_Coefficients()
        Implicit None
        Integer :: i
        Real*8, Allocatable :: temp_functions(:,:), temp_constants(:)
        Logical :: restore, need_custom

        ! In this routine, we may need to read in the custom file, which will overwrite
        ! ra_constants and ra_functions
        ! To be safe, we save the original ra_constants and ra_functions in temporary arrays
        ! We only modify the diffusion-related parts of the temporary arrays, then copy 
        ! everything over to ra_constants and ra_functions at the end
        Allocate(temp_functions(1:n_r, 1:n_ra_functions))
        Allocate(temp_constants(1:n_ra_constants))
        temp_functions(:,:) = ra_functions(:,:)
        ! Note that ra_constants is allocated up to max_ra_constants
        temp_constants(:) = ra_constants(1:n_ra_constants)      

        ! Figure out if we need to read anything from the custom file 
        ! (more complicated than it used to be because of the new scalar field coefficients)
        need_custom = .false. 
        If ( (nu_type .eq. 3) .or. (kappa_type .eq. 3) .or. (eta_type .eq. 3) ) Then
            need_custom = .true.
        EndIf

        Do i = 1, n_active_scalars
            If (kappa_chi_a_type(i) .eq. 3) Then
                need_custom = .true.
            Endif
        Enddo

        Do i = 1, n_passive_scalars
            If (kappa_chi_p_type(i) .eq. 3) Then
                need_custom = .true.
            Endif
        Enddo

        If ((.not. custom_reference_read) .and. need_custom) Then
            Call Read_Custom_Reference_File(custom_reference_file)
        EndIf

        Call Initialize_Diffusivity(temp_constants,temp_functions,nu_top,nu_type,nu_power,5,3,11)
        Call Initialize_Diffusivity(temp_constants,temp_functions,kappa_top,kappa_type,kappa_power,6,5,12)
        If (magnetism) Then
            Call Initialize_Diffusivity(temp_constants,temp_functions,eta_top,eta_type,eta_power,7,7,13)
        Else
            temp_functions(:,7) = 0.0d0 ! eta was already allocated, but never initialized since this
            temp_functions(:,13) = 0.0d0 ! run has magnetism = False. Explicitly set eta to zero
        Endif

        Do i = 1, n_active_scalars
            Call Initialize_Diffusivity(temp_constants,temp_functions,&
                kappa_chi_a_top(i),kappa_chi_a_type(i),kappa_chi_a_power(i),&
                11+(i-1)*2,15+(i-1)*2,16+(i-1)*2)
        Enddo

        Do i = 1, n_passive_scalars
            Call Initialize_Diffusivity(temp_constants,temp_functions,&
                kappa_chi_p_top(i),kappa_chi_p_type(i),kappa_chi_p_power(i),&
                11+(n_active_scalars+i-1)*2,15+(n_active_scalars+i-1)*2,16+(n_active_scalars+i-1)*2)
        Enddo

        ! copy temporary arrays to ra_constants and ra_functions
        ra_constants(1:n_ra_constants) = temp_constants(:)
        ra_functions(:,:) = temp_functions(:,:)
        DeAllocate(temp_functions, temp_constants)

    End Subroutine Initialize_Transport_Coefficients

    Subroutine Allocate_Transport_Coefficients()
        Implicit None

        Allocate(nu(1:N_r))
        Allocate(dlnu(1:N_r))
        Allocate(kappa(1:N_r))
        Allocate(dlnkappa(1:N_r))
        Allocate(kappa_chi_a(n_active_scalars,1:N_r))
        Allocate(dlnkappa_chi_a(n_active_scalars,1:N_r))
        Allocate(kappa_chi_p(n_passive_scalars,1:N_r))
        Allocate(dlnkappa_chi_p(n_passive_scalars,1:N_r))
        Allocate(eta(1:N_R))
        Allocate(dlneta(1:N_R))

    End Subroutine Allocate_Transport_Coefficients

    Subroutine Initialize_Diffusivity(temp_constants,temp_functions,xtop,xtype,xpower,ci,fi,dlnfi)
        Implicit None
        Real*8, Intent(InOut) :: xtop, temp_constants(:), temp_functions(:,:)
        Integer, Intent(In) :: xtype, ci, fi, dlnfi
        Real*8, Intent(In) :: xpower
        Character(len=2) :: ind

        ! This probably shouldn't be here anymore, since the diffusivity_types are independent
        ! of "reference_type .eq. 4". But leave on first pass (Loren, 12/22/2022)
        If (reference_type .eq. 4) Then
            If (xtop .le. 0) Then
                If (ra_constant_set(ci) .eq. 0) Then
                    If (my_rank .eq. 0) Then
                        Write(ind, '(I2)') ci
                        Call stdout%print('ERROR: constant c_'//Trim(Adjustl(ind))//' must be set in the custom reference file')
                    Endif
                Else
                    xtop = ra_constants(ci)
                Endif
            Endif
        Endif

        Select Case(xtype)
            Case(1)
                temp_constants(ci) = xtop
                temp_functions(:,fi) = 1.0d0
                temp_functions(:,dlnfi) = 0.0d0
            Case(2)
                temp_constants(ci) = xtop
                temp_functions(:,fi) = (temp_functions(:,1)/temp_functions(1,1))**xpower
                temp_functions(:,dlnfi) = xpower*temp_functions(:,8)
            Case(3)
                If ((ra_function_set(fi) .eq. 1) .and. (ra_constant_set(ci) .eq. 1)) Then
                    ! user specified both the constant and the function in the custom file
                    temp_constants(ci) = ra_constants(ci)                    
                    temp_functions(:,fi) = ra_functions(:,fi)
                    temp_functions(:,dlnfi) = ra_functions(:,dlnfi)
                ElseIf ((ra_function_set(fi) .eq. 1) .and. (ra_constant_set(ci) .eq. 0)) Then
                    ! User specified the function in the custom file, but not the constant
                    ! Get the constant from main_input (and normalize the function to its top value)
                    temp_constants(ci) = xtop
                    temp_functions(:,fi) = ra_functions(:,fi)/ra_functions(1,fi)
                    temp_functions(:,dlnfi) = ra_functions(:,dlnfi)
                Else
                    ! No function specified...user has messed up
                    If (my_rank .eq. 0) Then
                        Write(ind, '(I2)') fi
                        Call stdout%print('ERROR: function f_'//Adjustl(ind)//' must be set in the custom reference file')
                    EndIf
                EndIf

        End Select

    End Subroutine Initialize_Diffusivity

    Subroutine Restore_Transport_Defaults
        Implicit None

        If (Allocated(nu))           DeAllocate(nu)
        If (Allocated(kappa))        DeAllocate(kappa)
        If (Allocated(kappa_chi_a))  DeAllocate(kappa_chi_a)
        If (Allocated(kappa_chi_p))  DeAllocate(kappa_chi_p)
        If (Allocated(eta))          DeAllocate(eta)
        If (Allocated(dlnu))         DeAllocate(dlnu)
        If (Allocated(dlnkappa))     DeAllocate(dlnkappa)
        If (Allocated(dlnkappa_chi_a)) DeAllocate(dlnkappa_chi_a)
        If (Allocated(dlnkappa_chi_p)) DeAllocate(dlnkappa_chi_p)
        If (Allocated(dlneta))       DeAllocate(dlneta)

        If (allocated(W_Diffusion_Coefs_0) ) DeAllocate( W_Diffusion_Coefs_0)
        If (allocated(W_Diffusion_Coefs_1) ) DeAllocate( W_Diffusion_Coefs_1)

        If (allocated(dW_Diffusion_Coefs_0)) DeAllocate(dW_Diffusion_Coefs_0)
        If (allocated(dw_Diffusion_Coefs_1)) DeAllocate(dW_Diffusion_Coefs_1)
        If (allocated(dW_diffusion_coefs_2)) DeAllocate(dW_Diffusion_Coefs_2)

        If (allocated(S_Diffusion_Coefs_1) ) DeAllocate( S_Diffusion_Coefs_1)

        If (allocated(Z_Diffusion_Coefs_1) ) DeAllocate( Z_Diffusion_Coefs_1)
        If (allocated(Z_Diffusion_Coefs_0) ) DeAllocate( Z_Diffusion_Coefs_0)

        If (allocated(A_Diffusion_Coefs_1) ) DeAllocate( A_Diffusion_Coefs_1)

        kappa_type =1
        nu_type = 1
        eta_type = 1

        nu_top = 1.0d0
        kappa_top = 1.0d0
        eta_top = 1.0d0

        nu_power = 0
        eta_power = 0
        kappa_power = 0

        kappa_chi_a_type = 1
        kappa_chi_a_top = 1.0d0
        kappa_chi_a_power = 1.0d0
        kappa_chi_p_type = 1
        kappa_chi_p_top = 1.0d0
        kappa_chi_p_power = 1.0d0

    End Subroutine Restore_Transport_Defaults

    Subroutine Compute_Diffusion_Coefs()
        Implicit None
        ! These coefficients are nonzero only when nu and/or rho vary in radius
        ! They multiply derivatives of the field variables when
        ! constructing the diffusion terms in Sphere_Linear_Terms.F90.
        !////////////////////////////////////////+
        ! W Coefficients for W Equation
        Allocate(W_Diffusion_Coefs_0(1:N_R))
        Allocate(W_Diffusion_Coefs_1(1:N_R))
        W_Diffusion_Coefs_0 =     -nu*(4.0d0/3.0d0)*( dlnu*ref%dlnrho + ref%d2lnrho + ref%dlnrho/radius + &
            & 3.0d0*dlnu/radius )
        W_Diffusion_Coefs_1 = nu*(2.0d0*dlnu-ref%dlnrho/3.0d0)

        !/////////////////////////////////////
        ! W Coefficients for dWdr equation
        Allocate(DW_Diffusion_Coefs_0(1:N_R))
        Allocate(DW_Diffusion_Coefs_1(1:N_R))
        Allocate(DW_Diffusion_Coefs_2(1:N_R))
        DW_Diffusion_Coefs_2 = dlnu-ref%dlnrho
        DW_Diffusion_Coefs_1 = ref%d2lnrho+(2.0d0)/radius*ref%dlnrho+2.0d0/radius*dlnu+dlnu*ref%dlnrho
        DW_Diffusion_Coefs_0 = 2.0d0*ref%dlnrho/3.0d0+dlnu      !pulled out 2/r since that doesn't depend on rho or nu
            !include the factor of nu in these coefficients (and add minus sign for coefs 1 and 0)
        DW_Diffusion_Coefs_2 =  DW_Diffusion_Coefs_2*nu
        DW_Diffusion_Coefs_1 = -DW_Diffusion_Coefs_1*nu
        DW_Diffusion_Coefs_0 = -DW_Diffusion_Coefs_0*nu
        !//////////////////////////////////////// +
        ! S Coefficients for S Equation
        Allocate(S_Diffusion_Coefs_1(1:N_R))
        S_diffusion_Coefs_1 = kappa*(dlnkappa+ref%dlnrho+ref%dlnT)
        !//////////////////////////////////////// +
        ! chi Coefficients for chi Equation
        Allocate(chi_a_Diffusion_Coefs_1(n_active_scalars,1:N_R))
        chi_a_diffusion_Coefs_1 = kappa_chi_a*dlnkappa_chi_a
        Allocate(chi_p_Diffusion_Coefs_1(n_passive_scalars,1:N_R))
        chi_p_diffusion_Coefs_1 = kappa_chi_p*dlnkappa_chi_p
        !//////////////////////////////////////// +
        ! Z Coefficients for the Z Equation
        Allocate(Z_Diffusion_Coefs_0(1:N_R))
        Allocate(Z_Diffusion_Coefs_1(1:N_R))
        Z_Diffusion_Coefs_0 = -nu*( 2.0d0*dlnu/radius + ref%dlnrho*dlnu + &
            & ref%d2lnrho+2.0d0*ref%dlnrho/radius)
        Z_Diffusion_Coefs_1 = nu*(dlnu-ref%dlnrho)

        !////////////////////////////////////////
        ! A (vector potential) Coefficients
        If (magnetism) Then
            Allocate(A_Diffusion_Coefs_1(1:N_R))
            A_Diffusion_Coefs_1 = eta*dlneta
        Endif

    End Subroutine Compute_Diffusion_Coefs

    Subroutine Initialize_PDE_Coefficients()
        ! Sets Rayleigh's internal ReferenceInfo ("ref") structure from 
        ! ra_constants and ra_functions (i.e., the arrays more familiar to the user/coder)
        ! This routine is called at the end of everything the user/coder does to specify 
        ! ra_constants and ra_functions

        Implicit None
        Integer :: i
      
        ! Thermodynamic (historical "reference-state") variables
        ref%density(:) = ra_functions(:,1)
        ref%buoyancy_coeff(:) = ra_constants(2)*ra_functions(:,2)
        ref%temperature(:) = ra_functions(:,4)

        ref%dlnrho(:)  = ra_functions(:,8)
        ref%d2lnrho(:) = ra_functions(:,9)
        ref%dlnt(:) = ra_functions(:,10)
        ref%dsdr(:)     = ra_functions(:,14)

        ref%dpdr_w_term(:) = ra_constants(3)*ra_functions(:,1)
        ref%pressure_dwdr_term(:)= -ref%dpdr_w_term(:) 
        ref%viscous_amp(:) = 2.0*ra_constants(8)/ref%temperature(:)
        ref%lorentz_coeff = ra_constants(4)
        ref%ohmic_amp(:) = ra_constants(9)/(ref%density(:)*ref%temperature(:))

        Do i = 1, n_active_scalars
            ref%chi_buoyancy_coeff(i,:) = ra_constants(12+(i-1)*2)*ra_functions(:,2)
        Enddo

        ! Heating 
        ref%heating(:) = ra_functions(:,6)/(ra_functions(:,1)*ra_functions(:,4))
        ! Here we may rescale ref%heating to yield c_10 under integration. But if
        ! adjust_reference_heating is .True., then c_10 hasn't been set yet and will be
        ! (and ref%heating rescaled) in the BoundaryConditions.F90 file. 
        If (.not. adjust_reference_heating) Then
            ref%heating(:) = ra_constants(10)*ref%heating(:)
        Endif

        ! Coriolis coefficient 
        ! Before setting, allow angular_velocity to overwrite c_1
        If (angular_velocity .gt. 0) Then
            ra_constants(1) = 2.0d0*angular_velocity 
        Endif
        ref%coriolis_coeff = ra_constants(1)

        ! Diffusion coefficients
        ! Before setting, allow nu_top (and so on) to overwrite c_5 (and so on)
        ! NOTE: User should NOT set (e.g.) nu_top > 0 if f_3 (e.g.) has somehow been set to zero,
        ! since division by zero would then ensue here...
        Call Allocate_Transport_Coefficients

        If (nu_top .gt. 0.0d0) Then
            ra_constants(5) = nu_top/ra_functions(1,3)
        Endif
        nu(:) = ra_constants(5)*ra_functions(:,3)
        dlnu(:) = ra_functions(:,11)
        If (viscous_heating) Then
            Allocate(viscous_heating_coeff(1:N_R))
            viscous_heating_coeff(:) = ref%viscous_amp(:)*nu(:)
        Endif

        If (kappa_top .gt. 0.0d0) Then
            ra_constants(6) = kappa_top/ra_functions(1,5)
        Endif
        kappa(:) = ra_constants(6)*ra_functions(:,5)
        dlnkappa(:) = ra_functions(:,12)

        If (magnetism) Then
            If (eta_top .gt. 0.0d0) Then
                ra_constants(7) = eta_top/ra_functions(1,7)
            Endif
            eta(:) = ra_constants(7)*ra_functions(:,7)
            dlneta(:) = ra_functions(:,13)
            If (ohmic_heating) Then
                Allocate(ohmic_heating_coeff(1:N_R))
                ohmic_heating_coeff(:) = ref%ohmic_amp(:)*eta(:)
            Endif
        Endif

        Do i = 1, n_active_scalars
            If (kappa_chi_a_top(i) .gt. 0.0d0) Then
                ra_constants(11+(i-1)*2) = kappa_chi_a_top(i)/ra_functions(1,15+(i-1)*2)
            Endif
            kappa_chi_a(i,:) = ra_constants(11+(i-1)*2)*ra_functions(:,15+(i-1)*2)
            dlnkappa_chi_a(i,:) = ra_functions(:,16+(i-1)*2)
        Enddo

        Do i = 1, n_passive_scalars
            If (kappa_chi_p_top(i) .gt. 0.0d0) Then
                ra_constants(11+(n_active_scalars+i-1)*2) = kappa_chi_p_top(i)/ra_functions(1,15+(n_active_scalars+i-1)*2)
            Endif
            kappa_chi_p(i,:) = ra_constants(11+(n_active_scalars+i-1)*2)*ra_functions(:,15+(n_active_scalars+i-1)*2)
            dlnkappa_chi_p(i,:) = ra_functions(:,16+(n_active_scalars+i-1)*2)
        Enddo

        ! Finally, get the other "internal" diffusion coefficients
        Call Compute_Diffusion_Coefs()

    End Subroutine Initialize_PDE_Coefficients

End Module PDE_Coefficients
