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

Module TransportCoefficients
    Use ProblemSize
    Use ReferenceState
    Use Controls
    Use BoundaryConditions
    Use General_IO

    Implicit None
    Real*8, Allocatable :: nu(:), kappa(:), eta(:)
    Real*8, Allocatable :: dlnu(:), dlnkappa(:), dlneta(:)

    Real*8, Allocatable :: ohmic_heating_coeff(:)
    Real*8, Allocatable :: viscous_heating_coeff(:)
    !//////////

    Real*8, Allocatable :: W_Diffusion_Coefs_0(:), W_Diffusion_Coefs_1(:)
    Real*8, Allocatable :: dW_Diffusion_Coefs_0(:), dW_Diffusion_Coefs_1(:), dW_Diffusion_Coefs_2(:)
    Real*8, Allocatable :: S_Diffusion_Coefs_1(:), Z_Diffusion_Coefs_0(:), Z_Diffusion_Coefs_1(:)
    Real*8, Allocatable ::  A_Diffusion_Coefs_1(:)

    Integer :: kappa_type =1, nu_type = 1, eta_type = 1
    Real*8 :: nu_top = 1.0d0, kappa_top = 1.0d0, eta_top = 1.0d0
    Real*8 :: nu_power = 0, eta_power = 0, kappa_power = 0
    Real*8 :: eta_amp = 1.0d0

    Character*120 :: custom_eta_file = 'nothing'
    Character*120 :: custom_nu_file = 'nothing'
    Character*120 :: custom_kappa_file = 'nothing'

    Logical :: hyperdiffusion = .false.
    Real*8  :: hyperdiffusion_beta = 0.0d0
    Real*8  :: hyperdiffusion_alpha = 1.0d0


    Namelist /Transport_Namelist/ nu_type, kappa_type, eta_type, nu_power, kappa_power, eta_power, &
            & nu_top, kappa_top, eta_top, custom_nu_file, custom_eta_file, custom_kappa_file, &
              eta_amp, hyperdiffusion, hyperdiffusion_beta, hyperdiffusion_alpha


Contains

    Subroutine Compute_Diffusion_Coefs()
        ! These coefficients are nonzero only when nu and/or rho vary in radius
        ! The formulas here have been verified against my notes
        ! and against those implemented in ASH (which are slightly different from Brun et al. 2004
        ! due to sign errors in that paper)
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
        !DW_Diffusion_Coefs_0 = 2.0d0/radius+2.0d0*ref%dlnrho/3.0d0+dlnu
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

    Subroutine Initialize_Transport_Coefficients()
        Call Allocate_Transport_Coefficients
        If (.not. Dimensional_Reference) Then
            ! nu,kappa, and eta are based on the non-dimensionalization employed
            ! and are not read from the main_input file.
            nu_top    = ref%script_N_top
            kappa_top = ref%script_K_top
        Endif

        If (reference_type .eq. 4) Then
            nu_top = rayleigh_constants(5)*rayleigh_functions(1,3)
            kappa_top = rayleigh_constants(6)*rayleigh_functions(1,5)
        Endif

        Call Initialize_Nu()    ! Viscosity
        Call Initialize_Kappa() ! Thermal Diffusivity



        If (viscous_heating) Then
            Allocate(viscous_heating_coeff(1:N_R))
            viscous_heating_coeff(1:N_R) = ref%viscous_amp(1:N_R)*nu(1:N_R)
        Endif

        If (magnetism) Then
            If (.not. Dimensional_Reference) eta_top   = ref%script_H_top
            If (reference_type .eq. 4) eta_top = rayleigh_constants(7)*rayleigh_functions(1,7)
            Call Initialize_Eta()    ! Magnetic Diffusivity
            If (ohmic_heating) Then
                Allocate(ohmic_heating_coeff(1:N_R))
                ohmic_heating_coeff(1:N_R) = ref%ohmic_amp(1:N_R)*eta(1:N_R)
            Endif
        Endif


        Call Compute_Diffusion_Coefs()

        Call Transport_Dependencies()

        Call Write_Transport()

    End Subroutine Initialize_Transport_Coefficients

    Subroutine Write_Transport(filename)
        Implicit None
        Character*120, Optional, Intent(In) :: filename
        Character*120 :: transport_file
        Integer :: i, sig=314
        if (present(filename)) then
            transport_file = Trim(my_path)//filename
        else
            transport_file = Trim(my_path)//'transport'
        endif

        If (my_rank .eq. 0) Then
            Open(unit=15, file=transport_file, form='unformatted', status='replace', access='stream')
            Write(15) sig
            Write(15) n_r
            If (magnetism) Then
            ! Indicate whether or not there is magnetism
            ! (Determines size of output data file)
                Write(15) 1
            Else
                Write(15) 0
            Endif

            Write(15) (radius(i), i=1,n_r)

            Write(15) (nu(i),i=1, n_r)
            Write(15) (dlnu(i),i=1, n_r)

            Write(15) (kappa(i),i=1, n_r)
            Write(15) (dlnkappa(i),i=1, n_r)

            If (magnetism) Then
                Write(15) (eta(i),i=1,n_r)
                Write(15) (dlneta(i),i=1,n_r)
            Endif

            Close(15)
        Endif
    End Subroutine Write_Transport

    Subroutine Transport_Dependencies()
        Implicit None
        Real*8 :: fsun, lum_top, lum_bottom
        !Any odd boundary conditions that need the reference state
        ! or nu/kappa etc. can be set here
        ! As can any reference state quanitities, such as heating, that
        ! might depend on nu and kappa
        If (fix_tdt_bottom) Then
            ! Set the entropy gradient at the top based on Luminosity and kappa
            fsun = luminosity/four_pi/radius(1)/radius(1)
            dtdr_top = -fsun/kappa(1)/ref%density(1)/ref%temperature(1)
        Endif
        If (adjust_reference_heating ) Then
            !Renormalize the reference heating based on boundary fluxes
            lum_top    = 0.0d0
            lum_bottom = 0.0d0

            If ( fix_dtdr_top ) Then
                lum_top = -dtdr_top*kappa(1)*four_pi*(rmax**2)
                lum_top = lum_top*ref%density(1)*ref%temperature(1)
            Endif

            If ( fix_dtdr_bottom) Then
                lum_bottom = -dtdr_bottom*kappa(N_R)*four_pi*(rmin**2)
                lum_bottom = lum_bottom*ref%density(N_R)*ref%temperature(N_R)
            Endif

            ref%heating = ref%heating*(lum_top-lum_bottom)
            !If something has been set inconsistenly, this will result
            ! in zero reference heating
         Endif
    End Subroutine Transport_Dependencies


    Subroutine Allocate_Transport_Coefficients()
        Allocate(nu(1:N_r))
        Allocate(dlnu(1:N_r))
        Allocate(kappa(1:N_r))
        Allocate(dlnkappa(1:N_r))

        If (magnetism) Then
            Allocate(eta(1:N_R))
            Allocate(dlneta(1:N_R))
        Endif
    End Subroutine Allocate_Transport_Coefficients

    Subroutine Initialize_Nu()
        Select Case(nu_type)
            Case(1)    ! Constant nu
                nu(:) = nu_top
                dlnu(:) = 0.0d0
            Case(2)
                Call vary_with_density(nu,dlnu,nu_top, nu_power)
            Case(3)
                !Call get_custom_profile(nu,dlnu,custom_nu_file)
                nu(:) = rayleigh_constants(5)*rayleigh_functions(:,3)
                dlnu(:) = rayleigh_functions(:,11)
                nu_top = nu(1)

        End Select
    End Subroutine Initialize_Nu

    Subroutine Initialize_Kappa()
        Select Case(kappa_type)
            Case(1)    ! Constant Kappa
                kappa(:) = kappa_top
                dlnkappa(:) = 0.0d0
            Case(2)
                Call vary_with_density(kappa,dlnkappa,kappa_top, kappa_power)
            Case(3)
                !Call get_custom_profile(kappa,dlnkappa,custom_kappa_file)
                kappa(:) = rayleigh_constants(5)*rayleigh_functions(:,5)
                dlnkappa(:) = rayleigh_functions(:,12)
                kappa_top = kappa(1)
        End Select
    End Subroutine Initialize_Kappa

    Subroutine Initialize_Eta()
        Real*8, Allocatable :: tmp_arr(:,:)
        Character*120 :: eta_file = 'Eta_variation'
        Select Case(eta_type)
            Case(1)    ! Constant Eta
                eta(:) = eta_top
                dlneta(:) = 0.0d0
            Case(2)
                Call vary_with_density(eta,dlneta,eta_top, eta_power)
            Case(3)

                !Call get_custom_profile(eta,dlneta,custom_eta_file)
                !eta(:) = eta(:)*eta_top !this assume profile is 1 at the top

                eta(:) = rayleigh_constants(7)*rayleigh_functions(:,7)
                dlneta(:) = rayleigh_functions(:,13)
                eta_top = eta(1)

                If (my_rank .eq. 0) then
                    Allocate(tmp_arr(1:N_R,1:3))
                    tmp_arr(:,1) = radius(:)
                    tmp_arr(:,2) = eta(:)
                    tmp_arr(:,3) = dlneta(:)
                    Call Write_Profile(tmp_arr,eta_file)
                    DeAllocate(tmp_arr)
                Endif

        End Select
    End Subroutine Initialize_Eta

    Subroutine Get_Custom_Profile(coeff, dln, coeff_file)
        Real*8, Intent(InOut) :: coeff(:), dln(:)
        Real*8, Allocatable :: tmp_arr(:,:),dtemp(:,:,:,:),dtemp2(:,:,:,:)
        Integer :: dcheck(2)
        Character*120, Intent(In) :: coeff_file
        ! Reads density from a Rayleigh Profile File
        Allocate(tmp_arr(1:N_R,1:3))
        tmp_arr(:,:) = 0.0d0
        Call Read_Rayleigh_Array(coeff_file,tmp_arr,dims = dcheck)
        !Radius is assumed to be in tmp_arr(:,1)
        coeff(:) = tmp_arr(:,2)
        If (dcheck(2) >2) Then
            dln(:) = tmp_arr(:,3)
        Else
            Allocate(dtemp(1:n_r,1,1,2))
            Allocate(dtemp2(1:n_r,1,1,2))
            dtemp(:,:,:,:) = 0.0d0
            dtemp2(:,:,:,:) = 0.0d0
            dtemp(1:n_r,1,1,1) = coeff(1:n_r)
            Call gridcp%to_Spectral(dtemp,dtemp2)
            dtemp2((n_r*2)/3:n_r,1,1,1) = 0.0d0
            Call gridcp%d_by_dr_cp(1,2,dtemp2,1)
            dtemp2((n_r*2)/3:n_r,1,1,2) = 0.0d0  ! de-alias
            !transform back to physical
            Call gridcp%From_Spectral(dtemp2,dtemp)
            dln(:) = dtemp(:,1,1,2)/coeff
            DeAllocate(dtemp,dtemp2)
        Endif
        DeAllocate(tmp_arr)
    End Subroutine Get_Custom_Profile


    Subroutine Vary_With_Density(coeff, dln, coeff_top, coeff_power)
        Real*8, Intent(InOut) :: coeff(:), dln(:)
        Real*8, Intent(In) :: coeff_top, coeff_power
        ! Computes a transport coefficient and its logarithmic derivative
        ! using a density-dependent form for the coefficient:
        !        coeff = coeff_top*(rho/rho_top)**coeff_power
        coeff = coeff_top*(ref%density/ref%density(1))**coeff_power
        dln = coeff_power*ref%dlnrho
    End Subroutine Vary_With_Density

    Subroutine Restore_Transport_Defaults
        Implicit None

        If (Allocated(nu))       DeAllocate(nu)
        If (Allocated(kappa))    DeAllocate(kappa)
        If (Allocated(eta))      DeAllocate(eta)
        If (Allocated(dlnu))     DeAllocate(dlnu)
        If (Allocated(dlnkappa)) DeAllocate(dlnkappa)
        If (Allocated(dlneta))   DeAllocate(dlneta)

        If (allocated  (W_Diffusion_Coefs_0)) DeAllocate( W_Diffusion_Coefs_0)
        If (allocated  (W_Diffusion_Coefs_1)) DeAllocate( W_Diffusion_Coefs_1)

        If (allocated (dW_Diffusion_Coefs_0)) DeAllocate(dW_Diffusion_Coefs_0)
        If (allocated (dw_Diffusion_Coefs_1)) DeAllocate(dW_Diffusion_Coefs_1)
        If (allocated (dW_diffusion_coefs_2)) DeAllocate(dW_Diffusion_Coefs_2)

        If (allocated(S_Diffusion_Coefs_1)) DeAllocate(S_Diffusion_Coefs_1)

        If (allocated(Z_Diffusion_Coefs_1)) DeAllocate(Z_Diffusion_Coefs_1)
        If (allocated(Z_Diffusion_Coefs_0)) DeAllocate(Z_Diffusion_Coefs_0)

        If (allocated(A_Diffusion_Coefs_1)) DeAllocate(A_Diffusion_Coefs_1)


        kappa_type =1
        nu_type = 1
        eta_type = 1

        nu_top = 1.0d0
        kappa_top = 1.0d0
        eta_top = 1.0d0

        nu_power = 0
        eta_power = 0
        kappa_power = 0

        eta_amp = 1.0d0

        custom_eta_file = 'nothing'
        custom_nu_file = 'nothing'
        custom_kappa_file = 'nothing'

    End Subroutine Restore_Transport_Defaults
End Module TransportCoefficients
