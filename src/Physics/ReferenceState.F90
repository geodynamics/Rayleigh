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

! REFERENCE STATE MODULE
! Contains routines for initializing the reference state structure
! Reference state structure contains all information related to background
! stratification.  It DOES NOT contain transport variable (e.g. nu, kappa) information



Module ReferenceState
    Use ProblemSize
    Use Controls
    Use Math_Constants
    Use Math_Utility
    Use General_MPI, Only : BCAST2D
    Implicit None
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
        Real*8, Allocatable :: Buoyancy_Coeff(:)    ! Multiplies {S,T} in momentum eq. ..typically = gravity/cp
        Real*8, Allocatable :: dpdr_w_term(:)  ! multiplies d_by_dr{P/rho} in momentum eq.
        Real*8, Allocatable :: pressure_dwdr_term(:) !multiplies l(l+1)/r^2 (P/rho) in Div dot momentum eq.

        ! The following two terms are used to compute the ohmic and viscous heating
        Real*8, Allocatable :: ohmic_amp(:) !multiplied by {eta(r),H(r)}J^2 in dSdt eq.
        Real*8, Allocatable :: viscous_amp(:) !multiplied by {nu(r),N(r)}{e_ij terms) in dSdt eq.

        Real*8 :: script_N_Top ! If a nondimensional reference state is employed,
        Real*8 :: script_K_Top ! these are used in lieu of nu_top from the input file.
        Real*8 :: script_H_Top ! {N:nu, K:kappa, H:eta}

    End Type ReferenceInfo

    Integer, Parameter  :: eqn_coeff_version = 1

    ! Custom reference state variables
    Integer, Parameter  :: n_ra_constants = 10
    Integer, Parameter  :: n_ra_functions = 14
    Logical             :: with_custom_reference = .false.
    Logical             :: override_constants = .false.
    Logical             :: override_constant(1:n_ra_constants) = .false.
    Integer             :: ra_constant_set(1:n_ra_constants) = 0
    Integer             :: ra_function_set(1:n_ra_functions) = 0
    Logical             :: use_custom_constant(1:n_ra_constants) = .false.
    Logical             :: use_custom_function(1:n_ra_functions) = .false.
    Integer             :: with_custom_constants(1:n_ra_constants) = 0
    Integer             :: with_custom_functions(1:n_ra_functions) = 0
    Real*8              :: ra_constants(1:n_ra_constants) = 0.0d0
    Real*8, Allocatable :: ra_functions(:,:)
    Logical             :: custom_reference_read = .false.
    Character*120 :: custom_reference_file ='nothing'    

    Real*8, Allocatable :: s_conductive(:)

    Integer :: reference_type =1
    Integer :: heating_type = 0 ! 0 means no reference heating.  > 0 selects optional reference heating
    Integer :: cooling_type = 0 ! 0 means no reference cooling.  > 0 will ADD to any exisiting reference heating
    Real*8  :: Luminosity = 0.0d0 ! specifies the integral of the heating function
    Real*8  :: Heating_Integral = 0.0d0  !same as luminosity (for non-star watchers)
    Real*8  :: Heating_EPS = 1.0D-12  !Small number to test whether luminosity specified
    Logical :: adjust_reference_heating = .false.  ! Flag used to decide if luminosity determined via boundary conditions
    Real*8  :: heating_factor = 0.0d0, heating_r0 = 0.0d0  ! scaling and shifting factors for
    Real*8  :: cooling_factor = 0.0d0, cooling_r0 = 0.0d0  ! heating and cooling
    Type(ReferenceInfo) :: ref

    Real*8 :: pressure_specific_heat  = 1.0d0 ! CP (not CV)
    Real*8 :: poly_n = 0.0d0    !polytropic index
    Real*8 :: poly_Nrho = 0.0d0
    Real*8 :: poly_mass = 0.0d0
    Real*8 :: poly_rho_i =0.0d0
    Real*8 :: Angular_Velocity = 1.0d0

    !/////////////////////////////////////////////////////////////////////////////////////
    ! Nondimensional Parameters
    Real*8 :: Rayleigh_Number         = 1.0d0
    Real*8 :: Ekman_Number            = 1.0d0
    Real*8 :: Prandtl_Number          = 1.0d0
    Real*8 :: Magnetic_Prandtl_Number = 1.0d0
    Real*8 :: gravity_power           = 0.0d0
    Real*8 :: Dissipation_Number      = 0.0d0
    Real*8 :: Modified_Rayleigh_Number = 0.0d0
    Logical :: Dimensional_Reference = .false.  ! Changed depending on reference state specified

    !//////////////////////////////////////////////////////////////////////
    ! Development Section
    ! Everything below this line is related to in-development features.
    ! This is unsupported code.
    ! The code may well break if these variables are used in production.
    ! Email inquiries regarding these implementation or meaning of anything below this line
    ! will be politely ignored.
    Real*8, Allocatable :: paf_v2(:)
    Real*8, Allocatable :: paf_gv2(:)
    Real*8, Allocatable :: paf_p2(:)
    Real*8, Allocatable :: paf_gp2(:)

    Namelist /Reference_Namelist/ reference_type,poly_n, poly_Nrho, poly_mass,poly_rho_i, &
            & pressure_specific_heat, heating_type, luminosity, Angular_Velocity,     &
            & Rayleigh_Number, Ekman_Number, Prandtl_Number, Magnetic_Prandtl_Number, &
            & gravity_power, heating_factor, heating_r0, custom_reference_file,       &
            & cooling_type, cooling_r0, cooling_factor,        &
            & Dissipation_Number, Modified_Rayleigh_Number, Heating_Integral,         &
            & override_constants, override_constant, ra_constants, with_custom_constants, &
            & with_custom_functions, with_custom_reference
Contains

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
        Call Write_Reference()

    End Subroutine Initialize_Reference

    Subroutine Allocate_Reference_State
        Implicit None
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

        ref%Coriolis_Coeff = Zero
        ref%Lorentz_Coeff  = Zero
        ref%script_N_Top   = Zero
        ref%script_K_Top   = Zero
        ref%script_H_Top   = Zero

    End Subroutine Allocate_Reference_State

    Subroutine Constant_Reference()
        Implicit None
        Integer :: i
        Real*8 :: r_outer, r_inner, prefactor, amp, pscaling
        Character*12 :: dstring
        Character*8 :: dofmt = '(ES12.5)'
        ! devel variables (see note at top regarding devel)
        Real*8 :: pafk
        Real*8, Allocatable :: sink(:), cosk(:)
        Real*8, Allocatable :: array2d(:,:)
        Character*120 :: dvf
        Dimensional_Reference = .false.
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

        ref%density      = 1.0d0
        ref%dlnrho       = 0.0d0
        ref%d2lnrho      = 0.0d0
        ref%temperature  = 1.0d0
        ref%dlnT         = 0.0d0
        ref%dsdr         = 0.0d0

        amp = Rayleigh_Number/Prandtl_Number

        Do i = 1, N_R
            ref%Buoyancy_Coeff(i) = amp*(radius(i)/radius(1))**gravity_power
        Enddo

        pressure_specific_heat = 1.0d0
        Call initialize_reference_heating()
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

        If (.not. rotation) Then
            pscaling = 1.0d0
        Else
            pscaling = 1.0d0/Ekman_Number
        Endif

        !Define the various equation coefficients
        ref%dpdr_w_term(:)        =  ref%density*pscaling
        ref%pressure_dwdr_term(:) = -1.0d0*ref%density*pscaling
        ref%Coriolis_Coeff        =  2.0d0/Ekman_Number

        If (devel_physics) Then
            Allocate(paf_v2(1:N_R), paf_gv2(1:N_R))
            Allocate(paf_p2(1:N_R), paf_gp2(1:N_R))
            Allocate(sink(1:N_R), cosk(1:N_R) )
            pafk = 4.49d0/radius(1)  ! 4.49 is possibly a control parameter...possibly.
            sink = sin(pafk*radius)
            cosk = cos(pafk*radius)

            paf_p2(:) = sin(pafk*radius)/pafk/radius
            paf_p2 = paf_p2*paf_p2
            ! take derivative
            paf_gp2 = -2.0d0*paf_p2/radius
            paf_gp2 = paf_gp2+&
                      2.0d0*sin(pafk*radius)*cos(pafk*radius)/pafk/(radius**3)

            paf_v2 = sink*sink+(pafk**2)*r_squared*cosk*cosk
            paf_v2 = paf_v2 - two*pafk*radius*sink*cosk
            paf_v2 = paf_v2/(pafk**4)
            paf_v2 = paf_v2*OneOverRSquared*OneOverRSquared

            paf_gv2 = -two*(pafk**3)*r_squared*cosk*sink +two*(pafk**2)*radius*sink*sink
            paf_gv2 = paf_gv2/(pafk**4)
            paf_gv2 = paf_gv2*OneOverRSquared*OneOverRSquared

            paf_gv2 = paf_gv2-4.0d0*paf_v2*One_Over_R

            !renormalize (do gv2 first!)
            paf_gv2 = paf_gv2/maxval(paf_v2) *Rayleigh_Number
            paf_v2 = paf_v2/maxval(paf_v2) * Rayleigh_Number

            Allocate(array2d(1:N_R,5))
            array2d(:,1) = radius(:)
            array2d(:,2) = paf_v2(:)
            array2d(:,3) = paf_gv2(:)
            array2d(:,4) = paf_p2(:)
            array2d(:,5) = paf_gp2(:)
            dvf = 'paf.dat'
            Call Write_Profile(array2d,dvf)
            DeAllocate(array2d)
            DeAllocate(sink, cosk)
        Endif

        ref%script_N_top       = 1.0d0
        ref%script_K_top       = 1.0d0/Prandtl_Number
        ref%viscous_amp(1:N_R) = 2.0d0

        If (magnetism) Then
            ref%Lorentz_Coeff    = 1.0d0/(Magnetic_Prandtl_Number*Ekman_Number)
            ref%script_H_Top     = 1.0d0/Magnetic_Prandtl_Number
            ref%ohmic_amp(1:N_R) = ref%lorentz_coeff
        Else
            ref%Lorentz_Coeff    = 0.0d0
            ref%script_H_Top     = 0.0d0
            ref%ohmic_amp(1:N_R) = 0.0d0
        Endif

    End Subroutine Constant_Reference

    Subroutine Polytropic_ReferenceND()
        Implicit None
        Real*8 :: dtmp, otmp
        Real*8, Allocatable :: dtmparr(:), gravity(:)
        Character*12 :: dstring
        Character*8 :: dofmt = '(ES12.5)'
        Dimensional_Reference = .false.
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
        ref%temperature(:) = dtmp*Dissipation_Number*(dtmp*One_Over_R(:)-1.0D0)+1.0D0
        ref%density(:) = ref%temperature(:)**poly_n
        gravity = (rmax**2)*OneOverRSquared(:)
        ref%Buoyancy_Coeff = gravity*Modified_Rayleigh_Number*ref%density

        !Compute the background temperature gradient : dTdr = -Dg,  d2Tdr2 = 2*D*g/r (for g ~1/r^2)
        dtmparr = -Dissipation_Number*gravity
        !Now, the logarithmic derivative of temperature
        ref%dlnt = dtmparr/ref%temperature

        !And then logarithmic derivative of rho : dlnrho = n dlnT
        ref%dlnrho = poly_n*ref%dlnt

        !Now, the second logarithmic derivative of rho :  d2lnrho = (n/T)*d2Tdr2 - n*(dlnT^2)
        ref%d2lnrho = -poly_n*(ref%dlnT**2)

        dtmparr = (poly_n/ref%temperature)*(2.0d0*Dissipation_Number*gravity/radius) ! (n/T)*d2Tdr2
        ref%d2lnrho = ref%d2lnrho+dtmparr

        DeAllocate(dtmparr, gravity)

        ref%dsdr(:) = 0.0d0
        Call Initialize_Reference_Heating()

        ref%Coriolis_Coeff = 2.0d0
        ref%dpdr_w_term(:) = ref%density
        ref%pressure_dwdr_term(:) = -1.0d0*ref%density

        ref%script_N_top   = Ekman_Number
        ref%script_K_top   = Ekman_Number/Prandtl_Number
        ref%viscous_amp(1:N_R) = 2.0d0/ref%temperature(1:N_R)* &
                                 & Dissipation_Number/Modified_Rayleigh_Number

        If (magnetism) Then
            ref%Lorentz_Coeff    = Ekman_Number/(Magnetic_Prandtl_Number)
            ref%script_H_top     = Ekman_Number/Magnetic_Prandtl_Number

            otmp = (Dissipation_Number*Ekman_Number**2)/(Modified_Rayleigh_Number*Magnetic_Prandtl_Number**2)
            ref%ohmic_amp(1:N_R) = otmp/ref%density(1:N_R)/ref%temperature(1:N_R)
        Else
            ref%Lorentz_Coeff    = 0.0d0
            ref%script_H_Top     = 0.0d0
            ref%ohmic_amp(1:N_R) = 0.0d0
        Endif

    End Subroutine Polytropic_ReferenceND

    Subroutine Polytropic_Reference()
        Real*8 :: zeta_0,  c0, c1, d
        Real*8 :: rho_c, P_c, T_c,denom
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

        Dimensional_Reference = .true. ! This is actually the default

        ! Adiabatic, Polytropic Reference State (see, e.g., Jones et al. 2011)
        ! The following parameters are read from the input file.
        ! poly_n
        ! poly_Nrho
        ! poly_mass
        ! poly_rho_i

        ! Note that cp must also be specified.
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
        Gravity = Gravitational_Constant * poly_mass / Radius**2

        Ref%Density = rho_c * zeta**poly_n

        Ref%dlnrho = - poly_n * c1 * d / (zeta * Radius**2)
        Ref%d2lnrho = - Ref%dlnrho*(2.0d0/Radius-c1*d/zeta/Radius**2)

        Ref%Temperature = T_c * zeta
        Ref%dlnT = -(c1*d/Radius**2)/zeta

        Ref%dsdr = 0.d0

        Ref%Buoyancy_Coeff = gravity/Pressure_Specific_Heat*ref%density

        Deallocate(zeta, gravity)

        Call Initialize_Reference_Heating()

        ref%Coriolis_Coeff        = 2.0d0*Angular_velocity
        ref%dpdr_w_term(:)        = ref%density
        ref%pressure_dwdr_term(:) = -1.0d0*ref%density
        ref%viscous_amp(1:N_R)    = 2.0d0/ref%temperature(1:N_R)
        If (magnetism) Then
            ref%Lorentz_Coeff = 1.0d0/four_pi
            ref%ohmic_amp(1:N_R) = ref%lorentz_coeff/ref%density(1:N_R)/ref%temperature(1:N_R)
        Else
            ref%Lorentz_Coeff = 0.0d0
            ref%ohmic_amp(1:N_R) = 0.0d0
        Endif

    End Subroutine Polytropic_Reference

    Subroutine Initialize_Reference_Heating()
        Implicit None
        ! This is where a volumetric heating function Phi(r) is computed
        ! This function appears in the entropy equation as
        ! dSdt = Phi(r)
        ! Phi(r) may represent internal heating of any type.  For stars, this heating would be
        ! associated with temperature diffusion of the reference state and/or nuclear burning.

        If ( (heating_type .gt. 0) .or. (cooling_type .gt. 0) )Then
            If (.not. Allocated(ref%heating)) Allocate(ref%heating(1:N_R))
            ref%heating(:) = 0.0d0
        Endif

        If (heating_type .eq. 1) Then
            Call Constant_Entropy_Heating()
        Endif

        If (heating_type .eq. 2) Then
            Call Tanh_Reference_Heating()
        Endif

        If (heating_type .eq. 4) Then
            Call  Constant_Energy_Heating()
        Endif

        !///////////////////////////////////////////////////////////
        ! Next, compute a cooling function if desired and ADD it to
        ! whatever's in reference heating.
        If (cooling_type .eq. 2) Then
            Call Tanh_Reference_Cooling()
        Endif

        ! Heating Q_H and cooling Q_C are normalized so that:
        ! Int_volume rho T Q_H dV = 1 and Int_volume rho T Q_C dV = 1

        !If luminosity or heating_integral has been set,
        !we set the integral of ref%heating using that.

        !Otherwise, we will use the boundary thermal flux
        !To establish the integral
        If ( (heating_type .gt. 0) .or. (cooling_type .gt. 0) )Then
            adjust_reference_heating = .true.
            If (abs(Luminosity) .gt. heating_eps) Then
                adjust_reference_heating = .false.
                ref%heating = ref%heating*Luminosity
            Endif

            If (abs(Heating_Integral) .gt. heating_eps) Then
                adjust_reference_heating = .false.
                ref%heating = ref%heating*Heating_Integral
            Endif
        Endif
    End Subroutine Initialize_Reference_Heating

    Subroutine Augment_Reference()
        Implicit None

        If (my_rank .eq. 0) Then
            Call stdout%print('Reference state will be augmented.')
            Call stdout%print('Only heating and buoyancy may be modified.')
            Call stdout%print('Heating requires both c_10 and f_6 to be set.')
            Call stdout%print('Buoyancy requires both c_2 and f_2 to be set.')
            Call stdout%print('Reading from: '//Trim(custom_reference_file))
        Endif

        Call Read_Custom_Reference_File(custom_reference_file)

        If (use_custom_constant(10) .and. use_custom_function(6)) Then
            If (my_rank .eq. 0) Then
                Call stdout%print('Heating has been set to:')
                Call stdout%print('f_6*c_10')
                Call stdout%print(' ')
            Endif
            ref%heating(:) = ra_functions(:,6)/(ref%density*ref%temperature)*ra_constants(10)
        Endif

        If (use_custom_constant(2) .and. use_custom_function(2)) Then
            If (my_rank .eq. 0) Then
                Call stdout%print('Buoyancy_coeff has been set to:')
                Call stdout%print('f_2*c_2')
                Call stdout%print(' ')
            Endif
            ref%buoyancy_coeff(:) = ra_constants(2)*ra_functions(:,2)
        Endif

    End Subroutine Augment_Reference

    Subroutine Constant_Energy_Heating()
        Implicit None
        ! rho T dSdt = alpha : energy deposition per unit volume is constant
        ! dSdt = alpha/rhoT : entropy deposition per unit volume is ~ 1/Pressure
        ref%heating(:) = 1.0d0/shell_volume
        ref%heating = ref%heating/(ref%density*ref%temperature)
    End Subroutine Constant_Energy_Heating

    Subroutine Constant_Entropy_Heating()
        Implicit None
        Real*8 :: integral, alpha
        Real*8, Allocatable :: temp(:)

        ! dSdt = alpha : Entropy deposition per unit volume is constant
        ! rho T dSdt = rho T alpha : Energy deposition per unit volume is ~ Pressure

        ! Luminosity is specified as an input
        ! Phi(r) is set to alpha such that
        ! Integral_r=rinner_r=router (4*pi*alpha*rho(r)*T(r)*r^2 dr) = Luminosity
        Allocate(temp(1:N_R))

        temp = ref%density*ref%temperature
        Call Integrate_in_radius(temp,integral) !Int_rmin_rmax rho T r^2 dr
        integral = integral*4.0d0*pi  ! Int_V temp dV

        alpha = 1.0d0/integral

        ref%heating(:) = alpha
        DeAllocate(temp)
        !Note that in the boussinesq limit, alpha = Luminosity/Volume
    End Subroutine Constant_Entropy_Heating

    Subroutine Tanh_Reference_Heating()
        Implicit None
        Real*8 :: integral, alpha
        Real*8, Allocatable :: temp(:), x(:), temp2(:)
          Integer :: i
        Character*120 :: heating_file

        ! Luminosity is specified as an input
        ! Heating is set so that temp * 4 pi r^2 integrates to one Lsun
        ! Integral_r=rinner_r=router (4*pi*alpha*rho(r)*T(r)*r^2 dr) = Luminosity
        Allocate(temp(1:N_R))
        Allocate(temp2(1:N_R))
        Allocate(x(1:N_R))
        x = heating_factor*(radius-heating_r0)/ (maxval(radius)-minval(radius))
        ! x runs from zero to 1 if heating_r0 is min(radius)
        !Call tanh_profile(x,temp)
        Do i = 1, n_r
            temp2(i) = 0.5d0*(1.0d0-tanh(x(i))*tanh(x(i)))*heating_factor/(maxval(radius)-minval(radius))
        Enddo

        !temp2 = heating_factor*(1-(temp*temp))/ (maxval(radius)-minval(radius))

        temp2 = -temp2/(4*pi*radius*radius)

        Call Integrate_in_radius(temp2,integral)

        integral = integral*4.0d0*pi
        alpha = 1.0d0/integral

        ref%heating(:) = alpha*temp2/(ref%density*ref%temperature)
        If (my_rank .eq. 0) Then
            heating_file = Trim(my_path)//'reference_heating'
            Open(unit=15,file=heating_file,form='unformatted', status='replace',access='stream')
            Write(15)n_r
            Write(15)(radius(i),i=1,n_r)
            Write(15)(ref%heating(i), i = 1, n_r)
            Write(15)(temp(i), i = 1, n_r)
            Write(15)(temp2(i), i = 1, n_r)
            Write(15)(ref%density(i), i = 1, n_r)
            Write(15)(ref%temperature(i), i = 1, n_r)
        Endif

        DeAllocate(x,temp, temp2)
    End Subroutine Tanh_Reference_Heating

    Subroutine Tanh_Reference_Cooling()
        Implicit None
        Real*8 :: integral, alpha
        Real*8, Allocatable :: temp(:), x(:), temp2(:), cool_arr(:,:)
          Integer :: i
        Character*120 :: cooling_file

        ! Luminosity is specified as an input
        ! Heating is set so that temp * 4 pi r^2 integrates to one Lsun
        ! Integral_r=rinner_r=router (4*pi*alpha*rho(r)*T(r)*r^2 dr) = Luminosity
        Allocate(temp(1:N_R))
        Allocate(temp2(1:N_R))
        Allocate(x(1:N_R))
        x = cooling_factor*(radius-cooling_r0)/ (maxval(radius)-minval(radius)) ! x runs from zero to 1 if heating_r0 is min(radius)
        !Call tanh_profile(x,temp)
        Do i = 1, n_r
            temp2(i) = 0.5d0*(1.0d0-tanh(x(i))*tanh(x(i)))*cooling_factor/(maxval(radius)-minval(radius))
        Enddo

        !temp2 = heating_factor*(1-(temp*temp))/ (maxval(radius)-minval(radius))

        temp2 = -temp2/(4*pi*radius*radius)

        Call Integrate_in_radius(temp2,integral)

        integral = integral*4.0d0*pi
        alpha = 1.0d0/integral

        !/////////////////////
        ! The "Cooling" part comes in through the minus sign here - otherwise, this is identical to the heating function
        Do i = 1, n_r
            ref%heating(i) = ref%heating(i)-alpha*temp2(i)/(ref%density(i)*ref%temperature(i))
        Enddo
        !temp2 = ref%heating*ref%density*ref%temperature
        !Call Integrate_in_radius(temp2,integral)
        !write(6,*)'integral:  ', integral*4.0*pi

        cooling_file = Trim(my_path)//'reference_heating'
        Allocate(cool_arr(1:n_r,1:6))
        cool_arr(:,1) = radius
        cool_arr(:,2) = ref%heating
        cool_arr(:,3) = temp
        cool_arr(:,4) = temp2
        cool_arr(:,5) = ref%density
        cool_arr(:,6) = ref%temperature

        Call Write_Profile(cool_arr,cooling_file)
        DeAllocate(cool_arr)
        DeAllocate(x,temp, temp2)
    End Subroutine Tanh_Reference_Cooling

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

    Subroutine Write_Reference(filename)
        Implicit None
        Character*120, Optional, Intent(In) :: filename
        Character*120 :: ref_file
        Integer :: i,sig = 314
        if (present(filename)) then
            ref_file = Trim(my_path)//filename
        else
            ref_file = Trim(my_path)//'reference'
        endif

        If (my_rank .eq. 0) Then
            Open(unit=15,file=ref_file,form='unformatted', status='replace',access='stream')
            Write(15)sig
            Write(15)n_r
            Write(15)(radius(i),i=1,n_r)
            Write(15)(ref%density(i),i=1,n_r)
            Write(15)(ref%dlnrho(i),i=1,n_r)
            Write(15)(ref%d2lnrho(i),i=1,n_r)
            Write(15)(ref%temperature(i),i=1,n_r)
            Write(15)(ref%dlnT(i),i=1,n_r)
            Write(15)(ref%dsdr(i),i=1,n_r)
            Write(15)(ref%heating(i),i=1,n_r)
            Close(15)
        Endif
    End Subroutine Write_Reference

    Subroutine Write_Profile(arr,filename)
        Implicit None
        Character*120, Optional, Intent(In) :: filename
        Integer :: i,j,nq,sig = 314, nx
        Real*8, Intent(In) :: arr(1:,1:)
        nx = size(arr,1)
        nq = size(arr,2)
        If (my_rank .eq. 0) Then
            Open(unit=15,file=filename,form='unformatted', status='replace',access='stream')
            Write(15)sig
            Write(15)nx
            Write(15)nq
            Write(15)((arr(i,j),i=1,nx),j = 1, nq)

            Close(15)
        Endif
    End Subroutine Write_Profile

    Subroutine Get_Custom_Reference()
        Implicit None
        Integer :: i

        If (my_rank .eq. 0) Then
            Write(6,*)'Custom reference state specified.'
            Write(6,*)'Reading from: ', custom_reference_file
        Endif

        Call Read_Custom_Reference_File(custom_reference_file)

        Do i=1,7
            If (ra_function_set(i) .eq. 0) Then
                If (my_rank .eq. 0) Then
                    Write(6,*) "ERROR: function f_i must be set in the custom reference file",i
                Endif
            Endif
        Enddo

        ref%density(:) = ra_functions(:,1)
        ref%dlnrho(:)  = ra_functions(:,8)
        ref%d2lnrho(:) = ra_functions(:,9)
        ref%buoyancy_coeff(:) = ra_constants(2)*ra_functions(:,2)

        ref%temperature(:) = ra_functions(:,4)
        ref%dlnT(:) = ra_functions(:,10)

        ref%heating(:) = ra_functions(:,6)/(ref%density*ref%temperature)*ra_constants(10)
        
        ref%Coriolis_Coeff = ra_constants(1)
        ref%dpdr_w_term(:) = ra_constants(3)*ra_functions(:,1)
        ref%pressure_dwdr_term(:)= - ref%dpdr_w_term(:) 
        ref%viscous_amp(:) = 2.0/ref%temperature(:)*ra_constants(8)
        ref%Lorentz_Coeff = ra_constants(4)
        ref%ohmic_amp(:) = ra_constants(9)/(ref%density(:)*ref%temperature(:))

        ref%dsdr(:)     = ra_functions(:,14)

    End Subroutine Get_Custom_Reference

    Subroutine Write_Equation_Coefficients_File(filename)
        Character*120, Intent(In), Optional :: filename
        Character*120 :: ref_file
        Integer :: i, k, pi_integer
        pi_integer = 314

        If (present(filename)) Then
            ref_file = Trim(my_path)//filename
        Else
            ref_file = 'reference'
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
        Character*120, Intent(In), Optional :: filename
        Character*120 :: ref_file
        Integer :: pi_integer,nr_ref, eqversion
        Integer :: i, k, j
        Integer :: cset(1:n_ra_constants), fset(1:n_ra_functions)
        Real*8  :: input_constants(1:n_ra_constants)
        Real*8, Allocatable :: ref_arr_old(:,:), rtmp(:), rtmp2(:)
        Real*8, Allocatable :: old_radius(:)

        cset(:) = 0
        input_constants(:) = 0.0d0

        If (present(filename)) Then
            ref_file = Trim(my_path)//filename
        Else
            ref_file = 'reference'
        Endif

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
                    Write(6,*)'c: ', k, ra_constants(k)
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
                    call stdout%print("WARNING:  nr = nr_old.  Assuming grids are the same.")
                Endif
            Endif
            DeAllocate(ref_arr_old,old_radius)
            
            ! Finally, if the logarithmic derivatives of rho, T, nu, kappa, and eta were
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
        !Restore all values in this module to their default state.
        !Deallocates all allocatable module variables.
        reference_type =1
        heating_type = 0
        Luminosity =0.0d0
        heating_factor = 0.0d0
        heating_r0 = 0.0d0

        pressure_specific_heat = 0.0d0 ! CP (not CV)
        poly_n = 0.0d0
        poly_Nrho = 0.0d0
        poly_mass = 0.0d0
        poly_rho_i = 0.0d0

        Angular_Velocity = 1.0d0

        Rayleigh_Number         = 1.0d0
        Ekman_Number            = 1.0d0
        Prandtl_Number          = 1.0d0
        Magnetic_Prandtl_Number = 1.0d0
        gravity_power           = 0.0d0
        Dimensional_Reference = .true.
        custom_reference_file ='nothing'

        If (allocated(s_conductive)) DeAllocate(s_conductive)
        If (allocated(ref%Density)) DeAllocate(ref%density)
        If (allocated(ref%dlnrho)) DeAllocate(ref%dlnrho)
        If (allocated(ref%d2lnrho)) DeAllocate(ref%d2lnrho)
        If (allocated(ref%Temperature)) DeAllocate(ref%Temperature)
        If (allocated(ref%dlnT)) DeAllocate(ref%dlnT)
        If (allocated(ref%dsdr)) DeAllocate(ref%dsdr)
        If (allocated(ref%Buoyancy_Coeff)) DeAllocate(ref%Buoyancy_Coeff)
        If (allocated(ref%Heating)) DeAllocate(ref%Heating)

    End Subroutine Restore_Reference_Defaults

End Module ReferenceState
