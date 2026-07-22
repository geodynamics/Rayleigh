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

Module Sphere_Physical_Space
    Use Parallel_Framework
    Use Controls
    Use ProblemSize
    Use Fourier_Transform
    Use Spectral_Derivatives
    Use Fields
    Use Diagnostics_Interface, Only : PS_Output
    Use General_MPI, Only : global_max
    Use Timers
    Use ClockInfo
    Use PDE_Coefficients
    Use Math_Constants
    Use Benchmarking, Only : benchmark_checkup
    Implicit None

    Real*8, Allocatable :: tvar_eq(:,:,:)
    Real*8, Allocatable :: divu(:,:,:,:)
    Real*8, Allocatable :: gnu(:,:,:,:), gkappa(:,:,:,:)
    Real*8, Allocatable :: phi_visc(:,:,:), str(:,:,:,:)
    Real*8, Allocatable :: csquared(:,:,:)  ! sound speed squared
    Real*8, Allocatable :: work(:,:,:)
    Logical :: vars_allocated = .false.
    Integer, Parameter :: e_rr = 1, e_tt = 2, e_pp = 3
    Integer, Parameter :: e_rt = 4, e_rp = 5, e_tp = 6
    Logical :: debug = .false.

Contains
    
    Subroutine Physical_Space_Init()
        Implicit None
        Integer :: k, r, t
        Real*8, Allocatable :: cooling_profile(:)

        Allocate(cooling_profile(1:N_R))
        If (newtonian_cooling .and. (newtonian_cooling_profile_file .ne. '__nothing__')) Then
            cooling_profile(:) = newtonian_cooling_profile(:)
            If (my_rank .eq. 0) Then
                call stdout%print('Newtonian cooling is active.')
                call stdout%print('Cooling profile set from: '//TRIM(ADJUSTL(newtonian_cooling_profile_file)))
            Endif
        Else
            cooling_profile(:) = 1.0d0
        Endif
        
        ! Any persistant arrays needs for physical space routines can be
        ! initialized here.
        If (newtonian_cooling) Then
            Allocate(tvar_eq(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
            tvar_eq(:,:,:) = 0.0d0

            If (newtonian_cooling_type .eq. 1) Then
                ! No angular variation
                If (my_rank .eq. 0) call stdout%print('Newtonian cooling type = 1')
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do k =1, n_phi
                            tvar_eq(k,r,t) = newtonian_cooling_tvar_amp*newtonian_cooling_profile(r)
                        Enddo
                    Enddo
                Enddo
            Endif


            If (newtonian_cooling_type .eq. 2) Then
                ! Angular variation (ell=1,m=1, motivated by hot Jupiters)

                If (my_rank .eq. 0) call stdout%print('Newtonian cooling type = 2')
                Do t = my_theta%min, my_theta%max
                    Do r = my_r%min, my_r%max
                        Do k =1, n_phi
                            tvar_eq(k,r,t) = newtonian_cooling_tvar_amp*sintheta(t)*sinphi(k)
                            tvar_eq(k,r,t) = tvar_eq(k,r,t)*newtonian_cooling_profile(r)
                        Enddo
                    Enddo
                Enddo
            Endif
        Endif

        DeAllocate(cooling_profile)

    End Subroutine Physical_Space_Init

    Subroutine physical_space()
        Implicit None
        Integer :: i, t, r, k

        ! We aren't quite in physical space yet.
        ! 1st, get the phi derivatives
        Call StopWatch(dphi_time)%startclock()
        Call Phi_Derivatives()
        If (output_iteration) Then
            Call Diagnostics_Copy_and_Derivs()
        Endif
        Call StopWatch(dphi_time)%increment()

        
        !Brandon
        If (.not. vars_allocated .and. compressible) Then
            Allocate(divu(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:4))
            Allocate(phi_visc(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
            Allocate(gnu(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:4))
            Allocate(gkappa(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:4))
            Allocate(str(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:6))
            Allocate(csquared(1:n_phi+2,my_r%min:my_r%max,my_theta%min:my_theta%max))
            csquared = Zero
            Allocate(work(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
            vars_allocated = .true.
        Endif



        ! Next perform the FFT
        Call StopWatch(fft_time)%startclock()
        Call fft_to_physical(wsp%p3a,rsc = .true.)
        Call StopWatch(fft_time)%increment()


        Call StopWatch(pspace_time)%startclock()
        ! Convert all our terms of the form "sintheta var" to "var"
        Call StopWatch(sdiv_time)%startclock()
        Call sintheta_div(vtheta)    ! sintheta vtheta to vtheta etc.
        Call sintheta_div(vphi)
        Call sintheta_div(dvtdr)
        Call sintheta_div(dvpdr)
        Call sintheta_div(dtdt)
        Call sintheta_div(dvrdt)
        Call sintheta_div(dvpdp)
        Call sintheta_div(dvtdp)
        
        If (compressible) Then
            Call sintheta_div(drhodt)
            Call sintheta_div(dvtdt)
            Call sintheta_div(dvpdt)
            Call sintheta_div(d2vpdtdp)
            Call sintheta_div(d2vtdtdp)
            Call sintheta_div(d2vtdrdt)
            Call sintheta_div(d2vrdrdt)
        Endif

        do i = 1, n_active_scalars
          Call sintheta_div(dchiadt(i))  ! dchidt initially contains sin(theta) dsdtheta -- divide by sintheta
        end do
        do i = 1, n_passive_scalars
          Call sintheta_div(dchipdt(i))  ! dchidt initially contains sin(theta) dsdtheta -- divide by sintheta
        end do



        If (.not. compressible) Call Compute_dvtheta_by_dtheta()
        If (.not. compressible) Call Compute_dvphi_by_dtheta()

        If (magnetism) Then
            Call rsintheta_div(curlbtheta)
            Call rsintheta_div(curlbphi)
            Call rsintheta_div(Btheta)
            Call rsintheta_div(Bphi)
        Endif

        If (output_iteration) Then
            Call Diagnostics_Prep()
        Endif

        If (compressible) Then 
            Call sintheta_div(d2vtdt2)

            DO_IDX
                wsp%p3a(IDX,d2vtdt2) = wsp%p3a(IDX,d2vtdt2) - costheta(t)*wsp%p3a(IDX,dvtdt)
            END_DO 

            ! step 3:  divide by sin(theta), leaving d^2 vtheta/dtheta^2
            Call sintheta_div(d2vtdt2)
        Endif 



        Call StopWatch(sdiv_time)%increment()

        !////////////////////////////////////////////////////////////////////////
        !This is a good spot to do some simple diagnostic output while we debug the code
        !since velocity components, Pressure, and Temperature are all
        !in memory and in physical space at this point in time.

        Call ps_output(wsp%p3a, iteration,simulation_time)
        Call Benchmark_Checkup(wsp%p3a, iteration,simulation_time)
        !////////////////////////////////////////////////////////////////////////

        if (compressible) Call Compute_Sound_Speed
        Call Find_MyMinDT()    ! Piggyback CFL communication on transposes


        ! We are now ready to build the nonlinear terms
        Call wsp%construct('p3b')
        wsp%config = 'p3b'
        wsp%p3b(:,:,:,:) = 0.0d0
        !................................
        !Nonlinear Advection
        Call StopWatch(nl_time)%startclock()

        !Logic for computing compressible terms in physical space
        !Brandon 
        If (compressible) Then
            Call Compute_DivU()
            !Write(6, *) "DIVU,", Maxval(RHSP)
            Call Compute_Strain_Rate()
            !Write(6, *) "STR, ",Maxval(RHSP)
            Call Compute_Grad_Kappa()
            !Write(6, *) "KAPPA, ",Maxval(RHSP)
            Call Compute_Grad_Nu()
            !Write(6, *) "NU, " ,Maxval(RHSP)

            Call Compute_Phi_Visc()
            !Write(6, *) "PHI VISC, ", Maxval(RHSP)

            If (debug) Then
                Call Temperature_Diffusion()
                !Call Temperature_Heating()
            Else
                Call Temperature_Advection_Compressible()
                !Write(6, *) "TEP ADVECT COMP, ", Maxval(RHSP)
                Call Temperature_Compression()
                !Write(6, *) "TEMP COMPRESSIBLE, ", Maxval(RHSP)
                Call Temperature_Diffusion()
                !Write(6, *) "TEMP DIFFUSION,", Maxval(RHSP)
                !Call Temperature_Heating()
                Call Temperature_Viscous_Heating()
                !Write(6, *) "TEMP VISCSOUS HEAT,", Maxval(RHSP)

                Call Density_Advection()
                !Write(6, *) "DENSITY ADVECT", Maxval(RHSP)
                Call Density_Compression()
                !Write(6, *) "DENSITY COMPRESS", Maxval(RHSP)

                Call Velocity_Advection()
                !Write(6, *) "VELCOITY ADVECT", Maxval(RHSP)
                Call Velocity_Diffusion()
                !Write(6, *) "VELOCITY DIFFUSE", Maxval(RHSP)
                Call Pressure_Force()
                !Write(6, *) "PRESSURE FORCE", Maxval(RHSP)
                
                If (gravity) Call Compute_Gravity()
                If (rotation) Call Coriolis_Centrifugal()
            Endif
        Else
            Call Temperature_Advection()
            Call Volumetric_Heating()

            do i = 1, n_active_scalars
            Call chi_Advection(chiavar(i), dchiadr(i), dchiadt(i), dchiadp(i))
            Call chi_Source_function(chiavar(i), ref%chi_a_source(:,i))
            end do
            do i = 1, n_passive_scalars
            Call chi_Advection(chipvar(i), dchipdr(i), dchipdt(i), dchipdp(i))
            Call chi_Source_function(chipvar(i), ref%chi_p_source(:,i))
            end do
            
            If (viscous_heating) Call Compute_Viscous_Heating()

            Call Momentum_Advection_Radial()
            Call Momentum_Advection_Theta()
            Call Momentum_Advection_Phi()

        Endif

        If (magnetism) Then
            Call Compute_Ohmic_Heating()
            Call Compute_EMF()
        Endif

        Call StopWatch(nl_time)%increment()
        !...........................

        Call wsp%deconstruct('p3a')

        Call StopWatch(pspace_time)%increment()


        Call StopWatch(fft_time)%startclock()
        Call fft_to_spectral(wsp%p3b, rsc = .true.)
        Call StopWatch(fft_time)%increment()


        Call wsp%load_cargo(global_msgs)

        Call StopWatch(rtranspose_time)%startclock()
        Call wsp%reform()    ! Move to p2b
        Call StopWatch(rtranspose_time)%increment()
    End Subroutine Physical_Space

    Subroutine Compute_Sound_Speed()
        Implicit None
        Integer :: t, r, k
        Real*8 :: gfactor
        ! Checked:
        !       Fredy (8/22/19)
        !       Nick (8/27/19)
        gfactor = gas_gamma*(gas_gamma-1.0d0)*bigz
        DO_IDX
            csquared(IDX) = gfactor*FIELDSP(IDX,tvar)
        END_DO

    End Subroutine Compute_Sound_Speed

    Subroutine Compute_DivU()
        Implicit None
        Integer :: t, r, k
        ! Checked:
        !           Nick (8/27/19)

        !Note -- divu is allocated already before calling this routine
        DO_IDX
            divu(IDX,1) = wsp%p3a(IDX,vr)*two_over_r(r)+wsp%p3a(IDX,dvrdr) + &
                        One_Over_R(r)*( &
                        cottheta(t)*wsp%p3a(IDX,vtheta)+wsp%p3a(IDX,dvtdt) + &
                        csctheta(t)*wsp%p3a(IDX,dvpdp) )
        END_DO 

        ! Radial component of grad (div dot u)
        DO_IDX
            divu(IDX,2) = FIELDSP(IDX,d2vrdr2)                                 &
                          -OneOverRsquared(r)*(Two*FIELDSP(IDX,vr) +           &
                          cottheta(t)*FIELDSP(IDX,vtheta)+FIELDSP(IDX,dvtdt) + &
                          csctheta(t)*FIELDSP(IDX,dvpdp) ) +                   &
                          One_Over_R(r)*(Two*FIELDSP(IDX,dvrdr) +              &
                          cottheta(t)*FIELDSP(IDX,dvtdr)+FIELDSP(IDX,d2vtdrdt) + &
                          csctheta(t)*FIELDSP(IDX,d2vpdrdp) )            
        END_DO


        ! Theta component of grad (div dot  u)
        DO_IDX
            divu(IDX,3) =  One_Over_R(r)*FIELDSP(IDX,d2vrdrdt)          &
                           + OneOverRSquared(r)* (Two*FIELDSP(IDX,dvrdt) &
                           + FIELDSP(IDX,d2vtdt2)    &
                           + cottheta(t)*FIELDSP(IDX,dvtdt) &
                           - csctheta(t)*csctheta(t)*FIELDSP(IDX,vtheta) &
                           + csctheta(t)*FIELDSP(IDX,d2vpdtdp) &
                           -csctheta(t)*cottheta(t)*FIELDSP(IDX,dvpdp) )
        END_DO

        ! Phi component of grad (div dot u)
        DO_IDX
            divu(IDX,4) = OneOverRSquared(r)*csctheta(t)*(Two*FIELDSP(IDX,dvrdp) &
                         +FIELDSP(IDX,d2vtdtdp) &
                         +cottheta(t)*FIELDSP(IDX,dvtdp) &
                         +csctheta(t)*FIELDSP(IDX,d2vpdp2) ) &
                         +One_Over_R(r)*csctheta(t)*FIELDSP(IDX,d2vrdrdp) 
        END_DO


    End Subroutine Compute_DivU

    Subroutine Compute_dvtheta_by_dtheta()
        Implicit None
        Integer :: t, r, k

        DO_IDX
            wsp%p3a(IDX,dvtdt) = -wsp%p3a(IDX,vr)*(radius(r)*ref%dlnrho(r)+2.0d0) &
                                        - radius(r)*wsp%p3a(IDX,dvrdr) &
                                        - wsp%p3a(IDX,vtheta)*cottheta(t) &
                                        - wsp%p3a(IDX,dvpdp)*csctheta(t)
        END_DO

    End Subroutine Compute_dvtheta_by_dtheta

    Subroutine Compute_dvphi_by_dtheta()
        Implicit None
        Integer :: t, r,k
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            wsp%p3a(IDX,dvpdt) = radius(r)*wsp%p3a(IDX,zvar)+wsp%p3a(IDX,dvtdp)*csctheta(t) &
            -wsp%p3a(IDX,vphi)*cottheta(t)
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine Compute_dvphi_by_dtheta


    Subroutine chi_Advection(chivar, dchidr, dchidt, dchidp)
        Implicit None
        Integer :: chivar, dchidr, dchidt, dchidp
        Integer :: t,r,k
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        Do t = my_theta%min, my_theta%max
            Do r = my_r%min, my_r%max
                Do k =1, n_phi
                wsp%p3b(k,r,t,chivar) = -wsp%p3a(k,r,t,vr)*wsp%p3a(k,r,t,dchidr)     &
                                     - one_over_r(r)*(                           &
                                       wsp%p3a(k,r,t,dchidt)*wsp%p3a(k,r,t,vtheta) &
                                     + wsp%p3a(k,r,t,vphi)*wsp%p3a(k,r,t,dchidp)*csctheta(t) )

                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO
    End Subroutine chi_Advection

    Subroutine Temperature_Advection()
        Integer :: t,r,k
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        Do t = my_theta%min, my_theta%max
            Do r = my_r%min, my_r%max
                Do k =1, n_phi
                wsp%p3b(k,r,t,tvar) = -wsp%p3a(k,r,t,vr)*wsp%p3a(k,r,t,dtdr)     &
                                     - one_over_r(r)*(                           &
                                       wsp%p3a(k,r,t,dtdt)*wsp%p3a(k,r,t,vtheta) &
                                     + wsp%p3a(k,r,t,vphi)*wsp%p3a(k,r,t,dtdp)*csctheta(t) )

                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO

    End Subroutine Temperature_Advection

    Subroutine Volumetric_Heating()
        Implicit None
        Integer :: t,r,k
        If (heating_type .gt. 0) Then
            ! Added a volumetric heating to the energy equation
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    Do k =1, n_phi
                        wsp%p3b(k,r,t,tvar) = wsp%p3b(k,r,t,tvar)+ref%heating(r)
                    Enddo
                Enddo
            Enddo
            !$OMP END PARALLEL DO
        Endif


        If (newtonian_cooling) Then
            ! Added a volumetric heating to the energy equation
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    Do k =1, n_phi
                        wsp%p3b(k,r,t,tvar) = wsp%p3b(k,r,t,tvar) + &
                                      (tvar_eq(k,r,t) -wsp%p3b(k,r,t,tvar))/newtonian_cooling_time
                    Enddo
                Enddo
            Enddo
            !$OMP END PARALLEL DO
        Endif
    End Subroutine Volumetric_Heating


    Subroutine chi_Source_Function(chivar, source)
        Implicit None
        Integer, intent(in) :: chivar
        Real*8, intent(in) :: source(:)
        Integer :: t,r,k

        !$OMP PARALLEL DO PRIVATE(t,r,k)
        Do t = my_theta%min, my_theta%max
            Do r = my_r%min, my_r%max
                Do k =1, n_phi
                    wsp%p3b(k,r,t,chivar) = wsp%p3b(k,r,t,chivar)+source(r)
                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO

    End Subroutine chi_Source_Function

    Subroutine Compute_Viscous_Heating()
        Implicit None
        Integer :: t,r,k
        Real*8 :: tmp, tmp2
        Real*8, Allocatable :: htemp(:,:,:)

        Allocate(htemp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))

        ! Need to optimize these loops later, but for now, let's write this in
        ! easily debuggable way.

        !Contributions from E_rr, E_theta_theta     & E_phi_phi

        !$OMP PARALLEL DO PRIVATE(t,r,k,tmp,tmp2)
        Do t = my_theta%min, my_theta%max
            Do r = my_r%min, my_r%max
                Do k =1, n_phi
                    tmp = (wsp%p3a(k,r,t,dvpdp)*csctheta(t) +wsp%p3a(k,r,t,vr) &
                            +wsp%p3a(k,r,t,vtheta)*cottheta(t))*one_over_r(r)    !e_phi_phi
                    tmp2 = (wsp%p3a(k,r,t,dvtdt)+wsp%p3a(k,r,t,vr))*one_over_r(r) ! e_theta_theta
                    htemp(k,r,t) = wsp%p3a(k,r,t,dvrdr)*wsp%p3a(k,r,t,dvrdr)+tmp*tmp +tmp2*tmp2

                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO

        !E_r_phi
        !$OMP PARALLEL DO PRIVATE(t,r,k,tmp)
        Do t = my_theta%min, my_theta%max
            Do r = my_r%min, my_r%max
                Do k =1, n_phi
                    tmp = (wsp%p3a(IDX,dvrdp)*csctheta(t)- wsp%p3a(IDX,vphi))*one_over_r(r) &
                            +wsp%p3a(IDX,dvpdr) ! 2*e_r_phi

                    htemp(IDX) = htemp(IDX)+tmp*tmp*Half  ! +2 e_r_phi**2

                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO


        !E_r_theta
        !$OMP PARALLEL DO PRIVATE(t,r,k,tmp)
        Do t = my_theta%min, my_theta%max
            Do r = my_r%min, my_r%max
                Do k =1, n_phi
                    tmp = (wsp%p3a(IDX,dvrdt)-wsp%p3a(IDX,vtheta))*one_over_r(r) &
                            +wsp%p3a(IDX,dvtdr) ! 2*e_r_theta

                    htemp(IDX) = htemp(IDX)+tmp*tmp*Half   ! + 2+e_r_theta**2

                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO


        !E_phi_theta
        !$OMP PARALLEL DO PRIVATE(t,r,k,tmp)
        Do t = my_theta%min, my_theta%max
            Do r = my_r%min, my_r%max
                Do k =1, n_phi
                    tmp = (wsp%p3a(IDX,dvpdt) &
                            +wsp%p3a(IDX,dvtdp)*csctheta(t) &
                            -wsp%p3a(IDX,vphi)*cottheta(t) )*one_over_r(r)        ! 2*e_phi_theta

                    htemp(IDX) = htemp(IDX)+tmp*tmp*Half   ! + 2*e_phi_theta**2

                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO


        ! -1/3 (div dot v )**2
        !$OMP PARALLEL DO PRIVATE(t,r,k,tmp)
        Do t = my_theta%min, my_theta%max
            Do r = my_r%min, my_r%max
                Do k =1, n_phi
                    tmp = -wsp%p3a(IDX,vr)*ref%dlnrho(r)
                    htemp(IDX) = htemp(IDX)-tmp*tmp*one_third   ! + 2*e_phi_theta**2

                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO



        !$OMP PARALLEL DO PRIVATE(t,r,k)
        Do t = my_theta%min, my_theta%max
            Do r = my_r%min, my_r%max
                Do k =1, n_phi

                    wsp%p3b(k,r,t,tvar) = wsp%p3b(k,r,t,tvar)+viscous_heating_coeff(r)*htemp(k,r,t)

                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO

        DeAllocate(htemp)


    End Subroutine Compute_Viscous_Heating


    Subroutine Compute_Ohmic_Heating()
        Implicit None
        Integer :: t,r,k
        If (Ohmic_Heating) Then
            !We need a prefactor here for nondimensionalization

            !$OMP PARALLEL DO PRIVATE(t,r,k)
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    Do k =1, n_phi
                    wsp%p3b(k,r,t,tvar) = wsp%p3b(k,r,t,tvar) &
                                         + (wsp%p3a(k,r,t,curlbr)*wsp%p3a(k,r,t,curlbr) &
                                         + wsp%p3a(k,r,t,curlbtheta)*wsp%p3a(k,r,t,curlbtheta) &
                                         + wsp%p3a(k,r,t,curlbphi)*wsp%p3a(k,r,t,curlbphi))*ohmic_heating_coeff(r)
                    Enddo
                Enddo
            Enddo
            !$OMP END PARALLEL DO

        Endif
    End Subroutine Compute_Ohmic_Heating

    Subroutine Momentum_Advection_Radial()
        Implicit None
        Integer :: t,r,k

        ! Build -radius^2 [u dot grad u]_r

        If (momentum_advection) Then
            !$OMP PARALLEL DO PRIVATE(t,r,k)

            DO_IDX
                RHSP(IDX,wvar) = -FIELDSP(IDX,vr)*FIELDSP(IDX,dvrdr)*r_squared(r) &
                    - FIELDSP(IDX,vtheta) * ( FIELDSP(IDX,dvrdt)-FIELDSP(IDX,vtheta) )*radius(r)    &
                    - FIELDSP(IDX,vphi)*(FIELDSP(IDX,dvrdp)*csctheta(t)-FIELDSP(IDX,vphi) )*radius(r)
            END_DO

            !$OMP END PARALLEL DO
        Else
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,wvar) = 0.0d0
            END_DO
            !$OMP END PARALLEL DO
        Endif

        ! Add Coriolis Terms if so desired
        If (rotation) Then
        !    ! [- 2 z_hat cross u ]_r = 2 sintheta u_phi
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,wvar) = RHSP(IDX,wvar) + &
                    & ref%Coriolis_Coeff*sintheta(t)*FIELDSP(IDX,vphi)*R_squared(r)
            END_DO
            !$OMP END PARALLEL DO
        Endif


        ! Multiply advection/coriolis pieces by rho
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,wvar) = RHSP(IDX,wvar)*ref%density(r)
        END_DO
        !$OMP END PARALLEL DO


        If (magnetism .and. lorentz_forces) Then
            ! Add r_squared [JxB]_r
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,wvar)= RHSP(IDX,wvar) +r_squared(r)*ref%Lorentz_Coeff* &
                    (FIELDSP(IDX,curlbtheta)*FIELDSP(IDX,bphi)-FIELDSP(IDX,curlbphi)*FIELDSP(IDX,btheta))
            END_DO
            !$OMP END PARALLEL DO
        Endif

        ! Multiply by rho_*/rho = exp(s/c_P)  if pseudo-incompressible
        If (pseudo_incompressible) Then
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,wvar) = RHSP(IDX,wvar)*ref%exp_entropy(r)
            END_DO
            !$OMP END PARALLEL DO
        Endif


    End Subroutine Momentum_Advection_Radial

    Subroutine Compute_EMF()
        Implicit None
        Integer :: t,r,k

        ! Build the emf

        !$OMP PARALLEL DO PRIVATE(t,r,k)

        DO_IDX
            RHSP(IDX,emfr) = &
                  FIELDSP(IDX,vtheta) *  FIELDSP(IDX,bphi)  &
                - FIELDSP(IDX,vphi)     *  FIELDSP(IDX,btheta)
        END_DO

        !$OMP END PARALLEL DO
        
        !$OMP PARALLEL DO PRIVATE(t,r,k)

        DO_IDX
            RHSP(IDX,emftheta) = &
                - FIELDSP(IDX,vr) *  FIELDSP(IDX,bphi)  &
                + FIELDSP(IDX,vphi)   *  FIELDSP(IDX,br)
        END_DO

        !$OMP END PARALLEL DO

        !$OMP PARALLEL DO PRIVATE(t,r,k)

        DO_IDX
            RHSP(IDX,emfphi) = &
                  FIELDSP(IDX,vr)     *  FIELDSP(IDX,btheta)  &
                - FIELDSP(IDX,vtheta) *  FIELDSP(IDX,br)
        END_DO

        !$OMP END PARALLEL DO

        ! We need to divide by r/sintheta before taking the derivatives in the next space
        !$OMP PARALLEL DO PRIVATE(t,r,k)

        DO_IDX
            RHSP(IDX,emfphi) = RHSP(IDX,emfphi)*csctheta(t)*radius(r)
            RHSP(IDX,emftheta) = RHSP(IDX,emftheta)*csctheta(t)*radius(r)

        END_DO

        !$OMP END PARALLEL DO

    End Subroutine Compute_EMF

    Subroutine Momentum_Advection_Theta()
        Implicit None
        Integer :: t, r,k
        ! Build (radius/sintheta)[u dot grad u]_theta

        If (momentum_advection) Then
            ! First add all the terms that get multiplied by u_theta
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,pvar) = wsp%p3a(IDX,dvrdr)       &
                     + ( wsp%p3a(IDX,dvpdp)*csctheta(t)    & ! vphi/sintheta/r dvrdphi        !check this comment...
                     +   wsp%p3a(IDX,vtheta)*cottheta(t)   & !vtheta cot(theta)/r
                     +   wsp%p3a(IDX,vr) ) *one_over_r(r)                   &   !ur/r
                     +   wsp%p3a(IDX,vr)*ref%dlnrho(r) !ur dlnrho
            END_DO
            !$OMP END PARALLEL DO

            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,pvar) = -RHSP(IDX,pvar)*wsp%p3a(IDX,vtheta) & ! multiply by -u_theta
                    + wsp%p3a(IDX,vr  )*wsp%p3a(IDX,dvtdr)                         & ! vr dvthetadr
                    + wsp%p3a(IDX,vphi)*( wsp%p3a(IDX,dvtdp)*csctheta(t) & ! vphi/sintheta/r dvtheta dphi
                    - wsp%p3a(IDX,vphi )*cottheta(t) )*one_over_r(r)    ! vphi^2 cot(theta)/r

            END_DO
            !$OMP END PARALLEL DO
        Else
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,pvar) = 0.0d0
            END_DO
            !$OMP END PARALLEL DO
        Endif

        If (rotation) Then
            ! Add - the coriolis term (part of -RHS of theta)
            ! [2 z_hat cross u]_theta = -2 costheta u_phi

            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,pvar) = RHSP(IDX,pvar)- ref%Coriolis_Coeff*costheta(t)*FIELDSP(IDX,vphi)
            END_DO
            !$OMP END PARALLEL DO
        Endif

        ! Multiply advection/coriolis pieces by rho
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,pvar) = RHSP(IDX,pvar)*ref%density(r)
        END_DO
        !OMP END PARALLEL DO

        If (magnetism .and. lorentz_forces) Then
            ! Add -[JxB]_theta
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,pvar)= RHSP(IDX,pvar) &
                    - ref%Lorentz_Coeff*(FIELDSP(IDX,curlbphi)*FIELDSP(IDX,br)-FIELDSP(IDX,curlbr)*FIELDSP(IDX,bphi))
            END_DO
            !$OMP END PARALLEL DO
        Endif


        ! At this point, we have [u dot grad u]_theta
        ! Multiply by radius/sintheta so that we have r[u dot grad u]_theta/sintheta (getting ready for Z and dWdr RHS building)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,pvar) = RHSP(IDX,pvar)*radius(r)*csctheta(t)
        END_DO
        !$OMP END PARALLEL DO
        
                
        ! Multiply by rho_*/rho = exp(s/c_P)  if pseudo-incompressible
        If (pseudo_incompressible) Then
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,pvar) = RHSP(IDX,pvar)*ref%exp_entropy(r)
            END_DO
            !$OMP END PARALLEL DO
        Endif



    End Subroutine Momentum_Advection_Theta
    
    
    Subroutine Momentum_Advection_Phi()
        Implicit None
        Integer :: t, r, k
        ! Build (radius/sintheta)[u dot grad u]_phi


        If (momentum_advection) Then
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,zvar) = FIELDSP(IDX,vtheta)*(FIELDSP(IDX,zvar)  & ! terms multiplied by u_theta
                                        +FIELDSP(IDX,dvtdp)*csctheta(t)*one_over_r(r)) &
                    +FIELDSP(IDX,vr)*FIELDSP(IDX,dvpdr)    & ! radial advection
                    + FIELDSP(IDX,vphi) & ! terms multiplied by u_phi
                    * ( FIELDSP(IDX,dvpdp)*csctheta(t) + FIELDSP(IDX,vr))*one_over_r(r)
            END_DO
            !$OMP END PARALLEL DO

        Else
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,zvar) = 0.0d0
            END_DO
            !$OMP END PARALLEL DO
        Endif

        If (rotation) Then
            ! Add - Coriolis term (we are building -RHS of vphi)
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,zvar) = RHSP(IDX,zvar)                        &
                     + ref%Coriolis_Coeff*costheta(t)*FIELDSP(IDX,vtheta) &
                     + ref%Coriolis_Coeff*sintheta(t)*FIELDSP(IDX,vr)
            END_DO
            !$OMP END PARALLEL DO
        Endif

        ! Multiply advection/coriolis pieces by rho
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,zvar) = RHSP(IDX,zvar)*ref%density(r)
        END_DO
        !OMP END PARALLEL DO

        If (magnetism .and. lorentz_forces) Then
            ! Add -[JxB]_phi
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,zvar)= RHSP(IDX,zvar) - &
                    ref%Lorentz_Coeff*(FIELDSP(IDX,curlbr)*FIELDSP(IDX,btheta)-FIELDSP(IDX,curlbtheta)*FIELDSP(IDX,br))
            END_DO
            !$OMP END PARALLEL DO
        Endif




        ! At this point, we have [u dot grad u]_phi
        ! Multiply by radius/sintheta so that we have r[u dot grad u]_phi/sintheta (getting ready for Z and dWdr RHS building)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,zvar) = RHSP(IDX,zvar)*radius(r)*csctheta(t)
        END_DO
        !OMP END PARALLEL DO
        
        ! Multiply by rho_*/rho = exp(s/c_P)  if pseudo-incompressible
        If (pseudo_incompressible) Then
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,pvar) = RHSP(IDX,pvar)*ref%exp_entropy(r)
            END_DO
            !$OMP END PARALLEL DO
        Endif
        
    End Subroutine Momentum_Advection_Phi

    !/////////////////////////////////////////////////////
    ! COMPRESSIONAL ROUTINES START HERE  
    Subroutine Velocity_Advection()
        Implicit None
        Integer :: t, r,k

        ! Checked:
        !           Nick (8/20/19)

        ! Add advection terms to RHS of vr, vtheta, and vphi equations

        !Write(6, *) "Vr Max:min", Maxval(FIELDSP(:, :, :, vr)), Minval(FIELDSP(:, :, :, vr))
        !Write(6, *) "Vtheta Max:min", Maxval(FIELDSP(:, :, :, vtheta)), Minval(FIELDSP(:, :, :, vtheta))
        !Write(6, *) "Vphi Max:min", Maxval(FIELDSP(:, :, :, vphi)), Minval(FIELDSP(:, :, :, vphi))



        ! Radial component:  -[u dot grad u]_r
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vr) = -FIELDSP(IDX,vr)*FIELDSP(IDX,dvrdr) &
                - FIELDSP(IDX,vtheta) * ( FIELDSP(IDX,dvrdt)-FIELDSP(IDX,vtheta) )*One_Over_R(r)    &
                - FIELDSP(IDX,vphi)*(FIELDSP(IDX,dvrdp)*csctheta(t)-FIELDSP(IDX,vphi) )*One_Over_R(r)
        END_DO
        !$OMP END PARALLEL DO

        ! Theta component:  -[u dot grad u]_theta
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vtheta) = -FIELDSP(IDX,vr)*FIELDSP(IDX,dvtdr) &
                - One_Over_R(r)*FIELDSP(IDX,vtheta) * ( FIELDSP(IDX,dvtdt) + FIELDSP(IDX,vr) ) &
                - One_Over_R(r)*FIELDSP(IDX,vphi) & 
                * ( csctheta(t)*FIELDSP(IDX,dvtdp) - FIELDSP(IDX,vphi)*cottheta(t) ) ! - ...
        END_DO
        !$OMP END PARALLEL DO

        ! Phi component:  -[u dot grad u]_phi
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vphi) = -FIELDSP(IDX,vr)*FIELDSP(IDX,dvpdr) &
                - FIELDSP(IDX,vtheta)* One_Over_R(r)*( FIELDSP(IDX,dvpdt) & 
                + FIELDSP(IDX,vphi)*cottheta(t) ) & 
                - One_Over_R(r)*FIELDSP(IDX,vphi) & 
                * ( csctheta(t)*FIELDSP(IDX,dvpdp) + FIELDSP(IDX,vr) )
        END_DO
        !$OMP END PARALLEL DO



    End Subroutine Velocity_Advection

    Subroutine Pressure_Force()
        Implicit None
        Integer :: t, r,k
        Real*8 :: gfactor
        ! Add terms like du/dt = a Grad T + b T Grad lnrho
        ! Checked:
        !           Nick (8/20/19)


        
        gfactor = (1.0-gas_gamma)*bigz


        ! Radial component
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vr) = RHSP(IDX,vr) + gfactor*( FIELDSP(IDX,dtdr) &
                                        +FIELDSP(IDX,tvar)*FIELDSP(IDX,drhodr))
        END_DO
        !$OMP END PARALLEL DO

        If (Remove_Reference) Then     
            !$OMP PARALLEL DO PRIVATE(t,r,k)       
            DO_IDX
                RHSP(IDX,vr) = RHSP(IDX,vr) - gfactor*(ref%dT(r) +  ref%temperature(r)*ref%dlnrho(r))
            END_DO
            !$OMP END PARALLEL DO
        EndIf

        ! Theta component
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vtheta) = RHSP(IDX,vtheta) + gfactor*One_Over_R(r) & 
                                *( FIELDSP(IDX,dtdt) &
                                  +FIELDSP(IDX,tvar)*FIELDSP(IDX,drhodt) )
        END_DO
        !$OMP END PARALLEL DO

        ! Phi component
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vphi) = RHSP(IDX,vphi) +gfactor*One_Over_R(r)*csctheta(t) &
                                *( FIELDSP(IDX,dtdp) &
                                 + FIELDSP(IDX,tvar)*FIELDSP(IDX,drhodp) )
        END_DO
        !$OMP END PARALLEL DO


    End Subroutine Pressure_Force

    Subroutine Velocity_Diffusion()
        Implicit None
        Integer :: t, r,k
        Integer :: dvrdrdr,dvrdtdt,dvrdpdp      ! Maybe need to add this??
        Integer :: dvtdrdr,dvtdtdt,dvtdpdp
        Integer :: dvpdrdr,dvpdtdt,dvpdpdp
        ! These are pretty tough.
        ! See src_copy/Diagnostics/Diagnostics_Velocity_Diffusion lines 54, 87, and 115
        ! Compare what I did here against wikipedia Del^2{u}_{r,theta,phi} and see what's
        ! done there for inspiration.
        ! Rho is allowed to vary as r,theta, phi.   Nu only varies with radius.

        ! The goal for now is to compute nu * Del^2{u}_{r,theta,phi}.  There are other pieces too
        ! but those are a bit complicated and I need to do some careful math.
        
        ! Checked:
        !           Nick (8/27/2019)

        ! RADIAL DIFFUSION
        ! 1st:  nu*Del^2 {u_r}  (  NOT nu*[ Del^2{u} ]_r )
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vr) = RHSP(IDX,vr)+gnu(IDX,1)*(FIELDSP(IDX,d2vrdr2)+ &
                           Two_Over_R(r)*FIELDSP(IDX,dvrdr) -FIELDSP(IDX,hvr))
        END_DO
        !$OMP END PARALLEL DO

        !2nd:  Add geomtric terms to make this nu*[Del^2{u}]_r
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vr) = RHSP(IDX,vr) -2.0d0*OneOverRsquared(r)*gnu(IDX,1)*( &
                        FIELDSP(IDX,vr) + FIELDSP(IDX,dvtdt)+FIELDSP(IDX,vtheta)*cottheta(t) &
                        + csctheta(t)*FIELDSP(IDX,dvpdp) )
        END_DO
        !$OMP END PARALLEL DO

        !3rd:   Add radial component of gradient of div dot v (1/3 grad div.v)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vr) = RHSP(IDX,vr) + One_Third*gnu(IDX,1)*divu(IDX,2)
        END_DO
        !$OMP END PARALLEL DO


        !4th:
        
        ! Theta DIFFUSION
        ! 1st: nu*Del^2(u_t)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vtheta) = RHSP(IDX,vtheta)+gnu(IDX,1)*(FIELDSP(IDX,d2vtdr2)+ &
                           Two_Over_R(r)*FIELDSP(IDX,dvtdr) -FIELDSP(IDX,hvtheta))
        END_DO
        !$OMP END PARALLEL DO


        ! 2nd: Add geometric terms to make this nu*[Del^2{u}]_t
        DO_IDX
            RHSP(IDX,vtheta) = RHSP(IDX,vtheta) -OneOverRSquared(r)*gnu(IDX,1)*( &
                        csctheta(t)*csctheta(t)*FIELDSP(IDX,vtheta) &
                        -2.0d0*FIELDSP(IDX,dvrdt)+2.0d0*cottheta(t) &
                        *csctheta(t)*FIELDSP(IDX,dvpdp) )
        END_DO

        !3rd:   Add theta component of gradient of div dot v (1/3 grad div.v)

        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vtheta) = RHSP(IDX,vtheta) + One_Third*gnu(IDX,1)*divu(IDX,3)
        END_DO
        !$OMP END PARALLEL DO


        ! Phi DIFFUSION
        ! 1st: nu*Del^2(u_p)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vphi) = RHSP(IDX,vphi)+gnu(IDX,1)*(FIELDSP(IDX,d2vpdr2)+ &
                           Two_Over_R(r)*FIELDSP(IDX,dvpdr) -FIELDSP(IDX,hvphi))
        END_DO
        !$OMP END PARALLEL DO



        ! 2nd: Add geometric terms to make this nu*[Del^2{u}]_p
        DO_IDX
            RHSP(IDX,vphi) = RHSP(IDX,vphi)+csctheta(t)*OneOverRSquared(r) &
                        *(-csctheta(t)*FIELDSP(IDX,vphi)+2.0d0 &
                        *FIELDSP(IDX,dvrdp)+2.0d0*cottheta(t) &
                        *FIELDSP(IDX,dvtdp)  )*gnu(IDX,1)
        END_DO

        !$OMP END PARALLEL DO

        ! 3rd:  Add phi component of gradient of div dot v (1/3 grad div.v)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vphi) = RHSP(IDX,vphi) + One_Third*gnu(IDX,1)*divu(IDX,4)
        END_DO
        !$OMP END PARALLEL DO


        ! Finally, we compute 2/rho grad{rho nu} dot (e_ij - 1/3 divu del_ij)
        ! 2/rho grad {rho nu} = 2 grad(nu) + 2 nu grad(lnrho)
        ! Build each component and augment each RHS before building 
        ! next component

        ! I.  r-component of grad(mu)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            work(IDX) = Two*(gnu(IDX,2) + gnu(IDX,1)*FIELDSP(IDX,drhodr))
        END_DO
        !$OMP END PARALLEL DO

        !           ~~~~~~~~~~ r
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vr) = RHSP(IDX,vr) + work(IDX) * (STR(IDX,e_rr) - &
                           One_Third*divu(IDX,1))
        END_DO
        !$OMP END PARALLEL DO       

        !           ~~~~~~~~~~ theta
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vtheta) = RHSP(IDX,vtheta) + work(IDX) * STR(IDX,e_rt)
        END_DO
        !$OMP END PARALLEL DO     

        !           ~~~~~~~~~~ phi
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vphi) = RHSP(IDX,vphi) + work(IDX) * STR(IDX,e_rp)
        END_DO
        !$OMP END PARALLEL DO     

        ! II.  theta-component of grad(mu)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            work(IDX) = Two*One_Over_R(r)*(gnu(IDX,3) + gnu(IDX,1)*FIELDSP(IDX,drhodt))
        END_DO
        !$OMP END PARALLEL DO

        !           ~~~~~~~~~~ r
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vr) = RHSP(IDX,vr) + work(IDX) * STR(IDX,e_rt)
        END_DO
        !$OMP END PARALLEL DO       

        !           ~~~~~~~~~~ theta
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vtheta) = RHSP(IDX,vtheta) + work(IDX) * ( STR(IDX,e_tt) - &
                           One_Third*divu(IDX,1))
        END_DO
        !$OMP END PARALLEL DO     

        !           ~~~~~~~~~~ phi
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vphi) = RHSP(IDX,vphi) + work(IDX) * STR(IDX,e_tp)
        END_DO
        !$OMP END PARALLEL DO     

        ! III.  phi-component of grad(mu)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            work(IDX) = Two*One_Over_R(r)*csctheta(t)*(gnu(IDX,4) + gnu(IDX,1)*FIELDSP(IDX,drhodp))
        END_DO
        !$OMP END PARALLEL DO

        !           ~~~~~~~~~~ r
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vr) = RHSP(IDX,vr) + work(IDX) * STR(IDX,e_rp)
        END_DO
        !$OMP END PARALLEL DO       

        !           ~~~~~~~~~~ theta
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vtheta) = RHSP(IDX,vtheta) + work(IDX) * STR(IDX,e_tp)
        END_DO
        !$OMP END PARALLEL DO     

        !           ~~~~~~~~~~ phi
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vphi) = RHSP(IDX,vphi) + work(IDX) * (STR(IDX,e_pp) - &
                           One_Third*divu(IDX,1))
        END_DO
        !$OMP END PARALLEL DO     


    End Subroutine Velocity_Diffusion

    Subroutine Coriolis_Centrifugal()
        Implicit None
        Integer :: t, r,k
        ! Coriolis and centrifugal terms that appear on RHS of v-equations
        ! Checked:
        !           Nick (8/20/2019) 


        If (coriolis) Then
            ! Coriolis:  Radial
            ! [- z_hat cross u ]_r =  sintheta u_phi
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,vr) = RHSP(IDX,vr) + &
                    & ref%Coriolis_Coeff*sintheta(t)*FIELDSP(IDX,vphi) !*R_squared(r)
            END_DO
            !$OMP END PARALLEL DO

            ! Coriolis:  Theta
            ! [-z_hat cross u]_theta = costheta u_phi
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,vtheta) = RHSP(IDX,vtheta) + &
                    & ref%Coriolis_Coeff*costheta(t)*FIELDSP(IDX,vphi)
            END_DO
            !$OMP END PARALLEL DO

            ! Coriolis:  Phi  
            ! [- z_hat cross u]_phi = -u_theta costheta - u_r sintheta 
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,vphi) = RHSP(IDX,vphi) - &
                    & ref%Coriolis_Coeff*(costheta(t)*FIELDSP(IDX,vtheta) &
                    + sintheta(t)*FIELDSP(IDX,vr))
            END_DO
            !$OMP END PARALLEL DO
        Endif

        If (Centrifugal) Then
            ! Centrifugal Force:  
            ! Cylindrical radius s = r sin(theta)

            ! Centrifugal Radial:  s sin(theta)  
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,vr) = RHSP(IDX,vr) + &
                    ref%Centrifugal_Coeff*sintheta(t)*sintheta(t)*radius(r)
            END_DO
            !$OMP END PARALLEL DO

            ! Centrifugal Theta:  s cos(theta)  
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,vtheta) = RHSP(IDX,vtheta) + &
                    ref%Centrifugal_Coeff*costheta(t)*sintheta(t)*radius(r)
            END_DO
            !$OMP END PARALLEL DO
        Endif

        ! Centifugal force has no Phi-component
        !If (my_rank .eq. 0) Write(6,*) 'Centrifual: ', ref%centrifugal_Coeff

    End Subroutine Coriolis_Centrifugal

    Subroutine Compute_Gravity()
        Implicit None
        Integer :: t, r,k
        ! compressible gravity term that is solved explicitly 

        ! Gravity yerm:  Radial
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,vr) = RHSP(IDX,vr) - ref%gravity(r) 
        END_DO
        !$OMP END PARALLEL DO

    End Subroutine Compute_Gravity

    Subroutine Density_Advection()
        Implicit None
        Integer :: t, r,k
        ! density advection (dlnrho/dt = -v dot grad lnrho )
        ! rhovar is ln(rho), not rho

        ! Checked:
        !           Nick (8/20/19)
        !          
        !Write(6, *) "MIN:MAX dhroFIELD", Maxval(FIELDSP(:, :, :, drhodr)), Minval(FIELDSP(:, :, :, drhodr))
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,rhovar) = -FIELDSP(IDX,vr)*FIELDSP(IDX,drhodr)    &
                             - one_over_r(r)*(                     &
                             FIELDSP(IDX,drhodt)*FIELDSP(IDX,vtheta) &
                             + FIELDSP(IDX,vphi)*FIELDSP(IDX,drhodp)*csctheta(t) )
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine Density_Advection

    Subroutine Density_Compression()
        Implicit None
        Integer :: t, r,k
        ! dlnrho/dt = -div dot v
        ! Checked:
        !           Nick (8/20/19)
        !          

        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,rhovar) = RHSP(IDX,rhovar)-divu(IDX,1)
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine Density_Compression

    Subroutine Temperature_Advection_Compressible()
        Implicit None
        Integer :: t, r,k
        ! Temperature advection (dT/dt = -v dot grad T )
        ! Checked:
        !           Nick (8/20/19)
        !    

        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,tvar) = -FIELDSP(IDX,vr)*FIELDSP(IDX,dtdr)    &
                             - one_over_r(r)*(                     &
                             FIELDSP(IDX,dtdt)*FIELDSP(IDX,vtheta) &
                             + FIELDSP(IDX,vphi)*FIELDSP(IDX,dtdp)*csctheta(t) )
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine Temperature_Advection_Compressible


    Subroutine Temperature_Heating()
        Implicit None
        Integer :: t, r,k
        ! Real*8 :: delta
        ! delta = 0.01d0 ! square wave edges sharpeness
        ! Add the Q term to the temperature equation
        !  Just leave this one alone for now.
        !  I need to think about adding time-dependence into Q
        ! Note: Microwave is pulsed in a square wave. 
        !  when I get back from N.C.
        ! Added rounded square wave using atan and sine. 

        ! Still need to check

        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,tvar) = RHSP(IDX,tvar) + gas_gamma/Prandtl_Number
                !+(two*atan(sin(pulse_freq*simulation_time)/pulse_sharpness)/Pi + 1.0d0)/2.0d0 ! + Q
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine Temperature_Heating

    Subroutine Temperature_Viscous_Heating()
        Implicit None
        Integer :: t,r,k
        Real*8 :: zfactor
        zfactor = 1.0d0/bigz
        ! Add the PHI term to the temperature equation
        !$OMP PARALLEL DO PRIVATE(t,r,k)

        DO_IDX
            RHSP(IDX,tvar) = RHSP(IDX,tvar) +zfactor*Phi_Visc(IDX)
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine Temperature_Viscous_Heating

    Subroutine Temperature_Compression()
        Implicit None
        Integer :: t, r,k
        Real*8 :: gfactor
        ! Add the div dot v term to the temperature equation
        ! Checked:
        !           Nick (8/20/19)
        !    

        gfactor = 1.0D0-gas_gamma
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,tvar) = RHSP(IDX,tvar) + gfactor*FIELDSP(IDX,tvar)*divu(IDX,1)! + div V term
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine Temperature_Compression

    Subroutine Temperature_Diffusion()
        Implicit None
        Integer :: t, r,k
        Integer :: dtdtdt, dtdpdp
        Real*8  :: kcoeff
        Logical :: ddebug=.false.
        ! Add the diffusion term to the temperature equation
        ! Checked:
        !           Nick (8/20/19)
        !    

        kcoeff = gas_gamma/Prandtl_Number

        If (ddebug) Then
            DO_IDX
                RHSP(IDX,tvar) = kcoeff*gkappa(IDX,1)*( &
                    Two_Over_R(r)*FIELDSP(IDX,dtdr)+FIELDSP(IDX,d2tdr2) &
                    -FIELDSP(IDX,htvar) ) 
            END_DO
            !Write(6,*)'This branch'
        Else


        ! Kappa * del^2 T
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,tvar) = RHSP(IDX,tvar) +kcoeff*gkappa(IDX,1)*( &
                Two_Over_R(r)*FIELDSP(IDX,dtdr)+FIELDSP(IDX,d2tdr2) &
                -FIELDSP(IDX,htvar) ) ! This term is equivalent to those below
                !+OneOverRSquared(r)*(cottheta(t)*FIELDSP(IDX,dtdt) &
                !+FIELDSP(IDX,dtdtdt)+csctheta(t)*csctheta(t) &
                !*FIELDSP(IDX,dtdpdp))! + diffusion terms
        END_DO
        !$OMP END PARALLEL DO

        If (Remove_Reference) Then
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            DO_IDX
                RHSP(IDX,tvar) = RHSP(IDX,tvar) - kcoeff*gkappa(IDX,1)*(ref%d2T(r) + &
                2*ref%dT(r)/radius(r) + ref%dlnrho(r)*ref%dT(r)) &
                 - kcoeff*gkappa(IDX,2)*ref%dT(r)
            END_DO
            !$OMP END PARALLEL DO
        EndIf

        !kcoeff = 0.0d0
        ! Grad T dot Grad Kappa
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,tvar) = RHSP(IDX,tvar) +kcoeff*( &
                             FIELDSP(IDX,dtdr)*gkappa(IDX,2) + &  ! 2 is dkappa/dr
                             OneOverRSquared(r)*( &
                             FIELDSP(IDX,dtdt)*gkappa(IDX,3) + &  ! 3 is dkappa/dtheta
                             FIELDSP(IDX,dtdp)*gkappa(IDX,4)*csctheta(t)*csctheta(t))) ! 4 is dkappa/dphi
        END_DO
        !$OMP END PARALLEL DO



        ! Kapp Grad T dot Grad lnrho
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            RHSP(IDX,tvar) = RHSP(IDX,tvar) +kcoeff*gkappa(IDX,1)*( &
                             FIELDSP(IDX,dtdr)*FIELDSP(IDX,drhodr) + &
                             OneOverRSquared(r)*( &
                             FIELDSP(IDX,dtdt)*FIELDSP(IDX,drhodt) + & 
                             FIELDSP(IDX,dtdp)*FIELDSP(IDX,drhodp)*csctheta(t)*csctheta(t))) 
        END_DO
        !$OMP END PARALLEL DO
        Endif

    End Subroutine Temperature_Diffusion

    Subroutine Compute_Grad_Kappa()
        Implicit None
        Integer :: k,r,t
        ! Compute kappa and its 1st derivatives
        ! Kappa could be a function of position / depend on T, etc.
        ! For now, we set it (nondimensional kappa) to 1
        ! Checked:
        !           Nick (8/20/19)
 
        !nondimensional version
        !gkappa(:,:,:,1) = One
        !gkappa(:,:,:,2:4) = Zero

        Do r = my_r%min, my_r%max
            gkappa(:,r,:,1) = kappa(r)
            gkappa(:,r,:,2) = kappa(r)*dlnkappa(r)
            gkappa(:,r,:,3:4) = Zero
        Enddo
        
    End Subroutine Compute_Grad_Kappa

    Subroutine Compute_Grad_Nu()
        Implicit None
        Integer :: k,r,t
        ! Compute nu and its 1st derivatives
        ! Nu could be a function of position / depend on T, etc.
        ! For now, we set it (nondimensional nu) to 1
        ! Checked:
        !           Nick (8/20/19)
        !    

        !nondimensional version 
        !gnu(:,:,:,1) = One
        !gnu(:,:,:,2:4) = Zero

        Do r = my_r%min, my_r%max
            gnu(:,r,:,1) = nu(r)
            gnu(:,r,:,2) = nu(r)*dlnu(r)
            gnu(:,r,:,3:4) = Zero
        Enddo
        
    End Subroutine Compute_Grad_Nu

    Subroutine Compute_Strain_Rate()
        Implicit None
        Integer :: k,r,t
        ! We need this is a couple of places.  Let's just compute it once
        ! Checked:
        !       Fredy (8/22/19)
        !       Nick  (8/27/19)

        ! e_rr = du_r/dr
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            STR(IDX,e_rr) = FIELDSP(IDX,dvrdr)
        END_DO
        !$OMP END PARALLEL DO  

        ! e_theta_theta = 1/r( du_theta/dtheta + u_r)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            STR(IDX,e_tt) = One_Over_R(r)*(FIELDSP(IDX,vr)+FIELDSP(IDX,dvtdt))
        END_DO
        !$OMP END PARALLEL DO       

        ! e_phi_phi = 1/r( {1/sin(theta)}{du_phi/phi} + u_r + u_theta*cot(theta) )
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            STR(IDX,e_pp) = One_Over_R(r)*(FIELDSP(IDX,vr)+ &
                                           FIELDSP(IDX,dvpdp)*csctheta(t) + &
                                           FIELDSP(IDX,vtheta)*cottheta(t))
        END_DO
        !$OMP END PARALLEL DO      

        !e_theta_phi = 1/2[dvphi/dtheta - vphi*cotan(theta) +csc(theta)*dvtheta/dphi]
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            STR(IDX,e_tp) = Half*One_Over_R(r)*(FIELDSP(IDX,dvpdt) - &
                            FIELDSP(IDX,vphi)*cottheta(t) + &
                            csctheta(t)*FIELDSP(IDX,dvtdp) )
        END_DO
        !$OMP END PARALLEL DO  

        !e_r_phi = 1/2*[csc(theta)/r*du_r/dphi + du_phi/dr -u_phi/r]
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            STR(IDX,e_rp) = Half*(One_Over_R(r)*csctheta(t)*FIELDSP(IDX,dvrdp) &
                        + FIELDSP(IDX,dvpdr) - One_Over_R(r)*FIELDSP(IDX,vphi))
        END_DO
        !$OMP END PARALLEL DO  

        !e_r_theta = (1/2)*[ 1/r du_r/dtheta + du_theta/dr -u_theta/r]
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            STR(IDX,e_rt) = Half*One_Over_R(r)*(FIELDSP(IDX,dvrdt) - &
                            FIELDSP(IDX,vtheta)) + Half*FIELDSP(IDX,dvtdr)
        END_DO
        !$OMP END PARALLEL DO  

    End Subroutine Compute_Strain_Rate

    Subroutine Compute_Phi_Visc()
        Implicit None
        Integer :: k,r,t
        ! Computes Phi/rho   where
        ! Phi = 2*mu*{e_ij}{e_ij}-2/3*mu*(divu)^2
        ! mu  = rho*nu
        ! rho, nu, e_ij nondimensional
        ! Checked:
        !       Fredy (8/22/19)
        !       Nick  (8/27/19)

        ! Step 1.  divu term
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            Phi_Visc(IDX) = -One_Third*(divu(IDX,1))**2
        END_DO
        !$OMP END PARALLEL DO        

        ! Step 2.  Diagonal terms: e_rr^2 + e_tt^2 + e_pp^2
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            Phi_Visc(IDX) = Phi_Visc(IDX)+STR(IDX,e_rr)**2 + &
                            STR(IDX,e_tt)**2 + STR(IDX,e_pp)**2
        END_DO

        ! Step 3.  Off-diagonal terms
        ! 2*( e_r_theta^2 + e_r_phi^2 + e_theta_phi^2)
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            Phi_Visc(IDX) = Phi_Visc(IDX)+Two*(STR(IDX,e_rt)**2 + &
                            STR(IDX,e_rp)**2 + STR(IDX,e_tp)**2)
        END_DO
        !$OMP END PARALLEL DO  

        ! Step 4.  Multiply by 2*nu
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            Phi_Visc(IDX) = Phi_Visc(IDX)*Two*gnu(IDX,1)
        END_DO
        !$OMP END PARALLEL DO 

    End Subroutine Compute_Phi_Visc
    !END Compressional Routines
    !Some code can be recycled from routines below.
    !///////////////////////////////////////////////

    Subroutine Phi_Derivatives()
        Implicit None
        Integer :: i

        Call d_by_dphi(wsp%p3a,vr,dvrdp)
        Call d_by_dphi(wsp%p3a,vtheta,dvtdp)
        Call d_by_dphi(wsp%p3a,vphi,dvpdp)
        Call d_by_dphi(wsp%p3a,tvar,dtdp)
        do i = 1, n_active_scalars
          Call d_by_dphi(wsp%p3a,chiavar(i),dchiadp(i))
        end do
        do i = 1, n_passive_scalars
          Call d_by_dphi(wsp%p3a,chipvar(i),dchipdp(i))
        end do
    End Subroutine Phi_Derivatives
    Subroutine sintheta_div(ind)
        ! Divide by sintheta
        Implicit None
        Integer, Intent(In) :: ind
        Integer :: t,r,k
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            FIELDSP(IDX,ind) = FIELDSP(IDX,ind)*csctheta(t)
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine sintheta_div

    Subroutine rsintheta_div(ind)
        Implicit None
        !divide by rsintheta
        Integer, Intent(In) :: ind
        Integer :: t,r,k
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            FIELDSP(IDX,ind) = FIELDSP(IDX,ind)*csctheta(t)*one_over_r(r)
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine rsintheta_div

    Subroutine Find_MyMinDT()
        Implicit None
        Real*8 :: ovt2, ovht2, ovrt2
        Integer :: r
        Call StopWatch(ts_time)%startclock()

        ovt2 = 0.0d0    ! "over t squared"
        If (compressible) Then
            Do r = my_r%min, my_r%max
                ovht2 = Maxval(csquared(:,r,:) + wsp%p3a(:,r,:,vtheta)**2+wsp%p3a(:,r,:,vphi)**2) &
                                    *OneOverRSquared(r)*l_l_plus1(l_max) ! horizontal
                ovt2  = Max(ovt2, ovht2)
                ovrt2 = Maxval(csquared(:,r,:) + wsp%p3a(:,r,:,vr)**2)/(delta_r(r)**2)    ! radial
                ovt2  = Max(ovt2,ovrt2)
            Enddo

        Else            
            Do r = my_r%min, my_r%max
                ovht2 = Maxval(wsp%p3a(:,r,:,vtheta)**2+wsp%p3a(:,r,:,vphi)**2) &
                                    *OneOverRSquared(r)*l_l_plus1(l_max) ! horizontal
                ovt2  = Max(ovt2, ovht2)
                ovrt2 = Maxval(wsp%p3a(:,r,:,vr)**2)/(delta_r(r)**2)    ! radial
                ovt2  = Max(ovt2,ovrt2)
            Enddo
        Endif 

        If (magnetism) Then
            ! Check on alfven speed as well
            Do r = my_r%min, my_r%max
                ovht2 = Maxval(wsp%p3a(:,r,:,btheta)**2+wsp%p3a(:,r,:,bphi)**2) &
                                *OneOverRSquared(r)*l_l_plus1(l_max)/(ref%density(r))*ref%Lorentz_Coeff ! horizontal
                ovt2  = Max(ovt2, ovht2)
                ovrt2 = Maxval(wsp%p3a(:,r,:,br)**2)/(delta_r(r)**2)/(ref%density(r))*ref%Lorentz_Coeff    ! radial
                ovt2  = Max(ovt2,ovrt2)
            Enddo
        Endif

        global_msgs(1) = ovt2


        Call StopWatch(ts_time)%increment()
    End Subroutine Find_MyMinDT


    !/////////////////////////////////////////////////////
    ! Support routines for getting additional diagnostic fields sorted out


    Subroutine Diagnostics_Copy_and_Derivs()
        Implicit None
        !Copy everything from out auxiliary output buffer into the main buffer

        If (.not. compressible) wsp%p3a(:,:,:,dpdr) = cobuffer%p3a(:,:,:,dpdr_cb)
        If (.not. compressible) wsp%p3a(:,:,:,dpdt) = cobuffer%p3a(:,:,:,dpdt_cb)
        If (magnetism) Then
            wsp%p3a(:,:,:,dbrdr) = cobuffer%p3a(:,:,:,dbrdr_cb)
            wsp%p3a(:,:,:,dbtdr) = cobuffer%p3a(:,:,:,dbtdr_cb)
            wsp%p3a(:,:,:,dbpdr) = cobuffer%p3a(:,:,:,dbpdr_cb)
            wsp%p3a(:,:,:,dbpdt) = cobuffer%p3a(:,:,:, avar_cb)
            wsp%p3a(:,:,:,dbrdt) = cobuffer%p3a(:,:,:,dbrdt_cb)
        Endif

        !Everything we need is in main buffer - reset the auxiliary buffer
        Call cobuffer%deconstruct('p3a')
        cobuffer%config = 'p1a'

        !Take phi derivatives
        Call d_by_dphi(wsp%p3a,pvar,dpdp)
        If (magnetism) Then
            Call d_by_dphi(wsp%p3a,br,dbrdp)
            Call d_by_dphi(wsp%p3a,btheta,dbtdp)
            Call d_by_dphi(wsp%p3a,bphi,dbpdp)
        Endif


    End Subroutine Diagnostics_Copy_and_Derivs

    Subroutine Diagnostics_Prep()
        Implicit None
        Integer :: t,r,k
        !convert d/dr(p/rho) to dpdr
        
        if (.not. compressible) Then
            Call sintheta_div(dpdt)
            DO_IDX
                wsp%p3a(IDX,dpdr) = wsp%p3a(IDX,dpdr)*ref%density(r)+ &
                                    & wsp%p3a(IDX,pvar)*ref%dlnrho(r)
            END_DO
        Endif 

        If (magnetism) Then

            Call rsintheta_div(dbtdp)
            Call rsintheta_div(dbpdp)

            Call sintheta_div(dbrdt) !these do not have the one over r factor
            Call sintheta_div(dbpdr)
            Call sintheta_div(dbtdr)

            Call Compute_dbtheta_by_dtheta()
            Call Compute_dbphi_by_dtheta()

        Endif



    End Subroutine Diagnostics_Prep

    Subroutine Compute_dbtheta_by_dtheta()
        Implicit None
        Integer :: t, r,k

        DO_IDX
            wsp%p3a(IDX,dbtdt) = - wsp%p3a(IDX,br)*2.0d0 &
                                 - radius(r)*wsp%p3a(IDX,dbrdr) &
                                 - wsp%p3a(IDX,btheta)*cottheta(t) &
                                 - wsp%p3a(IDX,dbpdp)*csctheta(t)
        END_DO

    End Subroutine Compute_dbtheta_by_dtheta

    Subroutine Compute_dbphi_by_dtheta()
        Implicit None
        Integer :: t, r,k
        !Note: the A streamfunction was stored in dbpdt earlier.  We overwrite it with actual d B_phi d_theta now
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
            wsp%p3a(IDX,dbpdt) = radius(r)*wsp%p3a(IDX,dbpdt)+wsp%p3a(IDX,dbtdp)*csctheta(t) &
            -wsp%p3a(IDX,bphi)*cottheta(t)
        END_DO
        !$OMP END PARALLEL DO
    End Subroutine Compute_dbphi_by_dtheta


End Module Sphere_Physical_Space
