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

#include "indices.F"
Module Diagnostics_Interface
    Use ProblemSize
    Use Controls
    Use Spherical_IO
    Use Fields
    Use Legendre_Polynomials, Only : gl_weights
    Use ReferenceState
    Use TransportCoefficients
    Use Math_Constants
    Use Diagnostics_Base

    Use Diagnostics_Second_Derivatives

    Use Diagnostics_Mean_Correction

    Use Diagnostics_Velocity_Field
    Use Diagnostics_Magnetic_Field
    Use Diagnostics_Energies

    Use Diagnostics_Thermodynamic_Gradients
    Use Diagnostics_Thermal_Energies
    Use Diagnostics_Thermal_Equation

    Use Diagnostics_Vorticity_Field
    Use Diagnostics_Current_Density

    Use Diagnostics_Linear_Forces
    Use Diagnostics_Inertial_Forces
    Use Diagnostics_Angular_Momentum
    Use Diagnostics_Lorentz_Forces

    Use Diagnostics_KE_Flux

    Use Diagnostics_TurbKE_Budget
    Use Diagnostics_Axial_Field
    Use Diagnostics_Induction
    Use Diagnostics_Poynting_Flux

    Use Diagnostics_Miscellaneous

    Use Diagnostics_Custom

    Implicit None


    !///////////////////////////////////
    !Real*8, Allocatable :: qty(:,:,:)   ! This variable holds each quantity that we output
    !Real*8, Allocatable :: tmp1(:,:,:)
    !Real*8, Allocatable :: rweights(:), tweights(:)

    !//////////////////////////////////

    Integer, private :: reboot_count = 0  ! Number of times diagnostics has been rebooted during this run

Contains

    !//////////////////////////////////////////////////////////////////////////////////////////
    ! When entering PS_OUTPUT, the indices below may be used to reference the 4th dimension of
    ! the buffer array.  That array is dimensioned as:
    !   buffer(1:n_phi+2, my_r%min:my_r%max, my_theta%min:my_theta%max,1:nvariables)
    !
    ! The extra 2 in the first index is needed for the in-place FFTs.  Care should be taken
    !   to only loop over 1 to n_phi.
    !
    ! Each index along the 4th dimension of buffer corresponds to a different variable.
    !  These indices (left) and the variables they correspond to (right) are given below.

    ! Field variables:
    !   vr      -- radial velocity
    !   vtheta  -- theta velocity
    !   vphi    -- phi velocity
    !   tvar    -- temperature or entropy
    !   pvar    -- pressure
    !   zvar    -- l(l+1)*Z/r^2  where Z is the toroidal streamfunction

    ! Radial Derivatives:
    !   dvrdr   -- d(v_r)/dr
    !   dvtdr   -- d(v_theta)/dr
    !   dvpdr   -- d(v_phi)/dr
    !   dtdr    -- d(temperature or entropy)/dr


    ! Theta Derivatives:
    !   dvrdt   -- d(v_r)/dtheta
    !   dvtdt   -- d(v_theta)/dtheta
    !   dvpdt   -- d(v_phi)/dtheta
    !   dtdt    -- d(temperature or entropy)/dtheta
    !

    ! Phi Derivatives:
    !   dvrdp   --  d(v_r)/dphi
    !   dvtdp   --  d(v_theta)/dphi
    !   dvpdp   --  d(v_phi)/dphi
    !   dtdp    --  d(temperature or entropy)/dphi


    ! If Magnetism is On, six additional variables are present:
    !   br      -- radial magnetic field
    !   btheta  -- theta magnetic field
    !   bphi    -- phi magnetic field
    !   curlbr      -- [Del x B]_r
    !   curlbtheta  -- [Del x B]_theta
    !   curlbphi    -- [Del x B]_phi


    ! If Induction Output is needed for this iteration,
    !   the buffer also holds the derivatives of each
    !   component of B.

    ! Radial Derivatives:
    !   dbrdr   -- d(b_r)/dr
    !   dbtdr   -- d(b_theta)/dr
    !   dbpdr   -- d(b_phi)/dr


    ! Theta Derivatives:
    !   dbrdt   -- d(b_r)/dtheta
    !   dbtdt   -- d(b_theta)/dtheta
    !   dbpdt   -- d(b_phi)/dtheta


    ! Phi Derivatives:
    !   dbrdp   --  d(b_r)/dphi
    !   dbtdp   --  d(b_theta)/dphi
    !   dbpdp   --  d(b_phi)/dphi



    Subroutine PS_Output(buffer,iteration, current_time)
        Implicit None
        Integer, Intent(In) :: iteration
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Intent(In) :: current_time
        Real*8 :: over_n_phi, tmp, tmp2, tmp3

        Integer :: p,t,r, nfields, bdims(1:4), pass_num, k


        If (time_to_output(iteration)) Then
            Call Begin_Outputting(iteration)

            !//////////////////////////////////////////////////////////////////////////
            !We go ahead and compute two averages over the buffer:
            !   1.  The ell=0 component of each field in the buffer (requires comm)
            !   2.  The m=0   component of each field in the buffer (requires no comm)
            bdims = shape(buffer)
            nfields = bdims(4)
            Allocate(ell0_values(my_r%min:my_r%max,1:nfields))
            Allocate(m0_values(my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields))
            Call ComputeEll0(buffer,ell0_values)
            Call ComputeM0(buffer,m0_values)
            Call Compute_Fluctuations(buffer)

            Call Initialize_Mean_Correction()


            IF (need_second_derivatives) THEN
                Call Compute_Second_Derivatives(buffer)
            ENDIF

            Allocate(qty(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
            Allocate(tmp1(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
            Allocate(tmp4(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
            Allocate(tmp1d(1:N_R))
            over_n_phi = 1.0d0/dble(n_phi)



            Call Mean_Correction(buffer)    ! Remove ell=0 component from radial and theta forces


            ! All requested Shell_Average quantities are computed twice
            ! During the first pass, ell = 0 and m = 0 averages are computed
            !   (for the shell_average quantities)
            ! During the second pass, all quantities are computed,
            !   file output is conducted, and the previously
            !   computed averages are used for moments in the shell_average output
            ! Compute_quantity returns false on the first pass for everything but shell_averages
            !////////////////////////
            Do pass_num = 1, 2
                ! Set the averaging flag, so that all quantities or only shell averages are computed
                Call Set_Avg_Flag(pass_num)

                Call Compute_Velocity_Components(buffer)
                Call Compute_Vorticity_Field(buffer)
                Call Compute_Velocity_Second_Derivatives()
                Call Compute_Axial_Field(buffer)                

                Call Compute_Thermodynamic_Gradients(buffer)
                Call Compute_Thermal_Energy(buffer)
                Call Compute_Thermal_Second_Derivatives()
                Call Compute_Thermal_Equation_Terms(buffer)



                Call Compute_Kinetic_Energy(buffer)
                Call Compute_Angular_Momentum_Balance(buffer)
                Call Compute_Inertial_Terms(buffer)
                Call Compute_Linear_Forces(buffer)

                Call Compute_KE_Flux(buffer)

                Call Compute_TurbulentKE_Budget(buffer)



                Call Compute_Misc_Diagnostics(buffer)
                Call Custom_Hydro_Diagnostics(buffer)

                !////// Magnetic Quantities
                If (magnetism) Then
                    Call Compute_BField_Components(buffer)
                    Call Compute_Magnetic_Second_Derivatives()
                    Call Compute_Lorentz_Forces(buffer)
                    Call Compute_J_Components(buffer)
                    Call Compute_Induction_Terms(buffer)
                    Call Compute_Magnetic_Diffusion(buffer)
                    Call Compute_Magnetic_Energy(buffer)
                    Call Compute_Poynting_Flux(buffer)
                    Call Custom_MHD_Diagnostics(buffer)
                Endif
                If (pass_num .eq. 1) Call Finalize_Averages()
            Enddo


            DeAllocate(qty,tmp1,tmp1d,tmp4)
            Call Complete_Output(iteration, current_time)

            DeAllocate(ell0_values,m0_values)
            Call DeAllocate_Fluctuations()
            IF (need_second_derivatives) THEN
                Call d2buffer%deconstruct('p3a')
                DeAllocate(d2_ell0,d2_m0,d2_fbuffer)
            ENDIF
            Call Finalize_Mean_Correction
        Endif  ! time_to_output(iteration)
    End Subroutine PS_Output

    Subroutine Read_Output_Namelist()
        Implicit None
        Character*120 :: input_file
        input_file = Trim(my_path)//'main_input'

        ! First read the main input file
        Open(unit=20, file=input_file, status="old", position="rewind")
        Read(unit=20, nml=output_namelist)
        Close(20)

    End Subroutine Read_Output_Namelist

    Subroutine Initialize_Diagnostics()
        Implicit None
        Integer :: i, isize
        Real*8 :: delr

        Allocate(tweights(1:n_theta))
        tweights(:) = gl_weights(:)/2.0d0

        Allocate(rweights(1:n_r))

        If (chebyshev) Then
            rweights = radial_integral_weights

        Else
            Do i = 2, n_r-1
                delr = (radius(i-1)-radius(i+1))/2.0d0
                rweights(i) = delr*radius(i)**2
            Enddo
            delr = ( radius(1)-radius(2) )/ 2.0d0
            rweights(1) = delr*radius(1)**2

            delr = (radius(n_r-1)-radius(n_r))/2.0d0
            rweights(n_r) = delr*radius(n_r)**2

            rweights = rweights/sum(rweights)

        Endif

        Call Initialize_Spherical_IO(radius,sintheta,rweights,tweights,costheta,my_path)

        Call Initialize_Diagnostic_Indices()
        !DeAllocate(tweights)  !<---- Used to deallocate these.  We now use these for the computing the ell0 components
        !DeAllocate(rweights)

        !Call Set_Spherical_IO_Integration_Weights(gl_weights, r_int_weights)


        Call Initialize_Second_Derivatives()

        Call Initialize_Diagnostics_Buffer()

    End Subroutine Initialize_Diagnostics


    Function count_digits(n) result(ndigits)
        !Counts the number of digits in the integer n
        Implicit None
        Integer, Intent(in) :: n
        Integer :: ndigits, tmp
        ndigits = 1
        tmp = abs(n)/10
            Do While( tmp .gt. 0)
                ndigits = ndigits+1
                tmp = tmp/10
            Enddo
    End Function count_digits
    Subroutine Reboot_Diagnostics(iteration,force_reboot)
        Implicit None
        ! Checks to see if a reboot file is found.  If so, reboot the diagnostics
        Integer, Intent(In) :: iteration
        Logical, Intent(In), Optional :: force_reboot

        Character*120 :: reboot_file
        Character*1 :: ndigstr
        Character*6 :: digfmt
        Character*4 :: suffix
        Logical :: reboot_now = .false.

        Integer :: ndigits

        If (present(force_reboot)) Then
            If (force_reboot) Then
                !This functionality is used in conjunction with a full reboot.
                reboot_count = 0
                reboot_now = .true.
            Endif
        Else
            If (MOD(iteration,diagnostic_reboot_interval) .eq. 0) Then
                if (my_rank .eq. 0) WRITE(6,*)'REBOOTING DIAGNOSTICS!'
                ! Find the name of the current reboot file.
                ndigits = count_digits(reboot_count)
                Write(ndigstr,'(i1.1)') ndigits
                digfmt = '(i'//ndigstr//'.'//ndigstr//')'
                If (my_rank .eq. 0) Write(suffix,digfmt)reboot_count
                reboot_file = 'reboot_diagnostics_'//Trim(suffix)


                INQUIRE(file = reboot_file, exist = reboot_now)
                If (reboot_now) reboot_count = reboot_count+1
            Endif
        Endif


        If (reboot_now) Then
            If (my_rank .eq. 0) Call stdout%print('Reboot file found.  Rebooting diagnostics.')
            Call   CleanUP_Spherical_IO()
            Call   Read_Output_Namelist()
            Call Initialize_Diagnostics()
            reboot_now = .false.
        Endif
    End Subroutine Reboot_Diagnostics

End Module Diagnostics_Interface
