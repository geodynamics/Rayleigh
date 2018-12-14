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
    Use ReferenceState
    Use TransportCoefficients
    Use Math_Constants
    Use Benchmarking, Only : benchmark_checkup
    Use Stable_Plugin
    Implicit None

Contains
    Subroutine physical_space()
        Implicit None
        Integer :: t, r,k
        ! We aren't quite in physical space yet.
        ! 1st, get the phi derivatives
        Call StopWatch(dphi_time)%startclock()
        Call Phi_Derivatives()
        If (output_iteration) Then
            Call Diagnostics_Copy_and_Derivs()
        Endif
        Call StopWatch(dphi_time)%increment()



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



        Call Compute_dvtheta_by_dtheta()
        Call Compute_dvphi_by_dtheta()

        If (magnetism) Then
            Call rsintheta_div(curlbtheta)
            Call rsintheta_div(curlbphi)
            Call rsintheta_div(Btheta)
            Call rsintheta_div(Bphi)
        Endif

        If (output_iteration) Then
            Call Diagnostics_Prep()
        Endif





        Call StopWatch(sdiv_time)%increment()

        !////////////////////////////////////////////////////////////////////////
        !This is a good spot to do some simple diagnostic output while we debug the code
        !since velocity components, Pressure, and Temperature are all
        !in memory and in physical space at this point in time.

        Call ps_output(wsp%p3a, iteration,simulation_time)
        Call Benchmark_Checkup(wsp%p3a, iteration,simulation_time)
        !////////////////////////////////////////////////////////////////////////


        Call Find_MyMinDT()    ! Piggyback CFL communication on transposes


        ! We are now ready to build the nonlinear terms
        Call wsp%construct('p3b')
        wsp%config = 'p3b'
        wsp%p3b(:,:,:,:) = 0.0d0
        !................................
        !Nonlinear Advection
        Call StopWatch(nl_time)%startclock()

        Call Temperature_Advection()
        Call Volumetric_Heating()
        If (viscous_heating) Call Compute_Viscous_Heating()


        Call Momentum_Advection_Radial()
        Call Momentum_Advection_Theta()
        Call Momentum_Advection_Phi()


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

    Subroutine Compute_dvtheta_by_dtheta()
        Implicit None
        Integer :: t, r,k

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

        !================================ STABLE =================
        If (stable_flag) Then
            ! Add -v dot grad SU
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    Do k =1, n_phi
                    wsp%p3b(k,r,t,tvar) = wsp%p3b(k,r,t,tvar) &
                                        -wsp%p3a(k,r,t,vr)*grad_su%p3a(k,r,t,1)     &
                                         - one_over_r(r)*wsp%p3a(k,r,t,vtheta)*grad_su%p3a(k,r,t,2)

                    Enddo
                Enddo
            Enddo
            !$OMP END PARALLEL DO

            !Add - U_mean dot grad s
            !$OMP PARALLEL DO PRIVATE(t,r,k)
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    Do k =1, n_phi
                    wsp%p3b(k,r,t,tvar) = wsp%p3b(k,r,t,tvar) &
                                         -Vmean_r(r,t)*wsp%p3a(k,r,t,dtdr)    &
                                         -one_over_r(r)* &
                                         ( Vmean_theta(r,t) *wsp%p3a(k,r,t,dtdt) &
                                          +Vmean_phi(r,t)*csctheta(t)*wsp%p3a(k,r,t,dtdp) )

                    Enddo
                Enddo
            Enddo
            !$OMP END PARALLEL DO

        Endif
        !=============================== STABLE ==================
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
    End Subroutine Volumetric_Heating

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

        !==================== STABLE ==========================
        If (stable_flag)  Then
        !$OMP PARALLEL DO PRIVATE(t,r,k)

        DO_IDX
        ! Add terms due to mean flow to the momentum equation
        RHSP(IDX,wvar) = RHSP(IDX,wvar) - &
            (Vmean_r(IDXM)*FIELDSP(IDX,dvrdr) + FIELDSP(IDX,vr)*dVmean_rdr(IDXM))*r_squared(r) &
            - ( Vmean_theta(IDXM)*FIELDSP(IDX,dvrdt) + FIELDSP(IDX,vtheta)*dVmean_rdt(IDXM)    &
                +  Vmean_phi(IDXM)*FIELDSP(IDX,dvrdp)*csctheta(t) &
                -  2.d0*(Vmean_theta(IDXM)*FIELDSP(IDX,vtheta) + Vmean_phi(IDXM)*FIELDSP(IDX,vphi)) &
                ) * radius(r)
        !print*, &
        !    (Vmean_r(IDXM)*FIELDSP(IDX,dvrdr) + FIELDSP(IDX,vr)*dVmean_rdr(IDXM))*r_squared(r)
        END_DO

        !$OMP END PARALLEL DO
        Endif
        !=======================================================

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

                !==================== STABLE ==========================
                ! Add terms to the MHD Induction equation

        If (STABLE_flag) Then

           !$OMP PARALLEL DO PRIVATE(t,r,k)

            DO_IDX
            RHSP(IDX,emfr) = RHSP(IDX,emfr) + &
                                  Vmean_theta(IDXM) *  FIELDSP(IDX,bphi)  &
                - Vmean_phi(IDXM)   *  FIELDSP(IDX,btheta)
           END_DO

           !$OMP END PARALLEL DO

           !$OMP PARALLEL DO PRIVATE(t,r,k)

           DO_IDX
                        RHSP(IDX,emftheta) = RHSP(IDX,emftheta) + &
                                  Vmean_phi(IDXM) *  FIELDSP(IDX,br)  &
                - Vmean_r(IDXM)   *  FIELDSP(IDX,bphi)
           END_DO

           !$OMP END PARALLEL DO

           !$OMP PARALLEL DO PRIVATE(t,r,k)

           DO_IDX
                        RHSP(IDX,emfphi) = RHSP(IDX,emfphi) + &
                                  Vmean_r(IDXM) *  FIELDSP(IDX,btheta)  &
                - Vmean_theta(IDXM)  *  FIELDSP(IDX,br)
           END_DO

           !$OMP END PARALLEL DO

            ! This is only for the axisymmetric benchmark
            If (Poloidal_Source == 1) Then
                Call Alpha_Bmean()
                ! after computing Bmean, add it to the emf

                !$OMP PARALLEL DO PRIVATE(t,r,k)
                DO_IDX
                    RHSP(IDX,emfphi) = RHSP(IDX,emfphi) + &
                    kindy_alpha(r)*kindy_gtheta(t)*Bmean(t) &
                    / (1.0d0 + (Bmean(t)/5.0d4)**2)
                END_DO
                !$OMP END PARALLEL DO

            EndIf

        EndIf

                !==================== STABLE ==========================


        ! We need to divide by r/sintheta before taking the derivatives in the next space
        !$OMP PARALLEL DO PRIVATE(t,r,k)

        DO_IDX
            RHSP(IDX,emfphi) = RHSP(IDX,emfphi)*csctheta(t)*radius(r)
            RHSP(IDX,emftheta) = RHSP(IDX,emftheta)*csctheta(t)*radius(r)

        END_DO

        !$OMP END PARALLEL DO

    End Subroutine Compute_EMF

    !============================= STABLE ====================================
    Subroutine Alpha_Bmean()
       ! This is only called for Poloidal_Source = 1
       ! It's an axisymmetric, non-local alpha effect designed to test the
       ! code with the axisymmetric flux-transport dynamo benchmark of
       ! Jouve et al (2008)
       ! This routine computes Bmean(my_theta%min:my_theta%max)

       Real*8, allocatable :: buff(:)
       Real*8 :: rc, rr, onp, beta
       Integer :: r1, r2, r, k, t

       ! hardwired in for the benchmark
       rc = 0.7d0 * Radius(1)

       onp = 1.d0 / N_Phi

       ! first find the two radial values that bracket rc
       r1 = N_R
       rr = Radius(r1)
       do while (rr <= rc)
          r1 = r1 - 1
          rr = Radius(r1)
       enddo
       r2 = r1 + 1

       ! Interpolate in radius but really treat this as a
       ! sum with appropriate weights
       Allocate(buff(my_theta%min:my_theta%max))
       buff(:) = 0.d0

       beta = (rc - Radius(r2))/(Radius(r1)-Radius(r2))

       do t=my_theta%min,my_theta%max
          do r=my_r%min,my_r%max
             If (r == r1) Then
                do k = 1,N_phi
                   buff(t) = buff(t) + onp*beta*FIELDSP(IDX,bphi)
                enddo
             ElseIf (r == r2) Then
                do k = 1,N_phi
                   buff(t) = buff(t) + onp*(1.d0-beta)*FIELDSP(IDX,bphi)
                enddo
             EndIf
          enddo
       enddo

       ! sum over all columns (radii)
       Call dallsum1d(buff,Bmean,pfi%ccomm)

       ! I believe all the radial nodes (columns) know the result
       ! check here; Yes, it do
       !do t=my_theta%min,my_theta%max
       !   If (t == 50) Then
       !      print*,'Bmean check: ',my_column_rank,Bmean(t)
       !   EndIf
       !enddo

       Deallocate(buff)

    End Subroutine Alpha_Bmean

    !============================= STABLE ====================================


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
        !==================== STABLE ==========================
        If (STABLE_flag) Then
        ! Add terms due to mean flow to the momentum equation
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
              RHSP(IDX,pvar) = RHSP(IDX,pvar) + &
               Vmean_r(IDXM)*wsp%p3a(IDX,dvtdr)+wsp%p3a(IDX,vr)*dVmean_tdr(IDXM) + &
               ( Vmean_theta(IDXM)*wsp%p3a(IDX,dvtdt) + wsp%p3a(IDX,vtheta)*dVmean_tdt(IDXM) + &
                 Vmean_phi(IDXM)*wsp%p3a(IDX,dvtdp)*csctheta(t) + &
                 Vmean_theta(IDXM)*wsp%p3a(IDX,vr) + wsp%p3a(IDX,vtheta)*Vmean_r(IDXM) -  &
                 2.d0*cottheta(t)*Vmean_phi(IDXM)*wsp%p3a(IDX,vphi) )*one_over_r(r)
        END_DO
        !$OMP END PARALLEL DO
        EndIF
        !======================================================

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

        !==================== STABLE ==========================
        If (STABLE_flag) Then
        ! Add terms due to mean flow to the momentum equation
        !$OMP PARALLEL DO PRIVATE(t,r,k)
        DO_IDX
              RHSP(IDX,zvar) = RHSP(IDX,zvar) + &
               Vmean_r(IDXM)*wsp%p3a(IDX,dvpdr)+wsp%p3a(IDX,vr)*dVmean_pdr(IDXM) + &
               ( Vmean_theta(IDXM)*wsp%p3a(IDX,dvpdt) + wsp%p3a(IDX,vtheta)*dVmean_pdt(IDXM) + &
                 Vmean_phi(IDXM)*wsp%p3a(IDX,dvpdp)*csctheta(t) + &
                 Vmean_phi(IDXM)*wsp%p3a(IDX,vr) + wsp%p3a(IDX,vphi)*Vmean_r(IDXM) +  &
                 cottheta(t)*(Vmean_phi(IDXM)*wsp%p3a(IDX,vtheta)+Vmean_theta(IDXM)*wsp%p3a(IDX,vphi)) &
               ) * one_over_r(r)
        END_DO
        !$OMP END PARALLEL DO
        EndIf
        !======================================================

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
    End Subroutine Momentum_Advection_Phi
    Subroutine Phi_Derivatives()
        Implicit None
        Integer :: r,t,k


        Call d_by_dphi(wsp%p3a,vr,dvrdp)
        Call d_by_dphi(wsp%p3a,vtheta,dvtdp)
        Call d_by_dphi(wsp%p3a,vphi,dvpdp)
        Call d_by_dphi(wsp%p3a,tvar,dtdp)
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
        Do r = my_r%min, my_r%max
            ovht2 = Maxval(wsp%p3a(:,r,:,vtheta)**2+wsp%p3a(:,r,:,vphi)**2) &
                                *OneOverRSquared(r)*l_l_plus1(l_max) ! horizontal
            ovt2  = Max(ovt2, ovht2)
            ovrt2 = Maxval(wsp%p3a(:,r,:,vr)**2)/(delta_r(r)**2)    ! radial
            ovt2  = Max(ovt2,ovrt2)
        Enddo
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
        Integer :: t,r,k
        !Copy everything from out auxiliary output buffer into the main buffer

        wsp%p3a(:,:,:,dpdr) = cobuffer%p3a(:,:,:,dpdr_cb)
        wsp%p3a(:,:,:,dpdt) = cobuffer%p3a(:,:,:,dpdt_cb)

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
        Call sintheta_div(dpdt)
        !convert d/dr(p/rho) to dpdr
        DO_IDX
            wsp%p3a(IDX,dpdr) = wsp%p3a(IDX,dpdr)*ref%density(r)+ &
                                & wsp%p3a(IDX,pvar)*ref%dlnrho(r)
        END_DO

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
