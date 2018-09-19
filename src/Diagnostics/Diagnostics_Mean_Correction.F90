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


Module Diagnostics_Mean_Correction
    Use Diagnostics_Base
    Use Spherical_IO
    Use Diagnostics_ADotGradB

Contains

    Subroutine Initialize_Mean_Correction()
        Implicit None
        Logical :: compute_full_full = .false.
        Logical :: compute_fluct_fluct = .false.
        Logical :: compute_mean_mean = .false.
        ncorrect = 0

        !//////////////////////////////////////////////////
        ! We only correct the radial forces
        ! (we remove the ell=0 mean).  This should not be necessary
        ! for the theta/phi forces


        !////////////////////////////////////////////////////////////////////////////////
        ! Coriolis Corrections (phi is currently untouched in Diagnostics_Linear_Forces
        If (compute_quantity(coriolis_force_r) .or. &
            compute_quantity(coriolis_mforce_r) ) Then
            ncorrect = ncorrect+1
            cforce_r = ncorrect
        Endif

        !If (compute_quantity(coriolis_force_theta) .or. &
        !    compute_quantity(coriolis_mforce_theta) ) Then
        !    ncorrect     = ncorrect+1
        !    cforce_theta = ncorrect
        !Endif

        !If (compute_quantity(coriolis_force_phi) .or. &
        !    compute_quantity(coriolis_mforce_phi) ) Then
        !    ncorrect   = ncorrect+1
        !    cforce_phi = ncorrect
        !Endif


        !//////////////////////////////////////////////////////////////////////////////////
        ! Advective terms
        !//////////////////////////////////////////////
        ! Advective terms
        compute_full_full = .false.
        If (compute_quantity(v_grad_v_r))     compute_full_full = .true.
        If (compute_quantity(v_grad_v_theta)) compute_full_full = .true.
        If (compute_quantity(v_grad_v_phi))   compute_full_full = .true.
        If (compute_quantity(advec_work))     compute_full_full = .true.


        ! -- v v
        If (compute_full_full ) Then

            If (compute_quantity(v_grad_v_r) .or. compute_quantity(advec_work)) Then
                ncorrect=ncorrect+1
                aforce_r = ncorrect
            Endif

            !If (compute_quantity(v_grad_v_theta) .or. compute_quantity(advec_work)) Then
            !    ncorrect=ncorrect+1
            !    aforce_theta = ncorrect
            !Endif

            !If (compute_quantity(v_grad_v_phi) .or. compute_quantity(advec_work)) Then
            !    ncorrect=ncorrect+1
            !    aforce_phi = ncorrect
            !Endif

        Endif

        ! ---- v' v'
        compute_fluct_fluct = .false.
        If (compute_quantity(vp_grad_vp_r))     compute_fluct_fluct = .true.
        If (compute_quantity(vp_grad_vp_theta)) compute_fluct_fluct = .true.
        If (compute_quantity(vp_grad_vp_phi))   compute_fluct_fluct = .true.
        If (compute_quantity(advec_work_ppp))   compute_fluct_fluct = .true.
        If (compute_quantity(advec_work_mpp))   compute_fluct_fluct = .true.

        If (compute_fluct_fluct) Then

            If (compute_quantity(vp_grad_vp_r) .or. compute_quantity(advec_work_ppp) &
                .or. compute_quantity(advec_work_mpp)) Then

                ncorrect=ncorrect+1
                aforcepp_r = ncorrect
            Endif

            !If (compute_quantity(vp_grad_vp_theta) .or. compute_quantity(advec_work_ppp) &
            !    .or. compute_quantity(advec_work_mpp)) Then
            !    ncorrect=ncorrect+1
            !    aforcepp_theta = ncorrect
            !Endif

            !If (compute_quantity(vp_grad_vp_phi) .or. compute_quantity(samom_advec_pp) &
            !    .or. compute_quantity(advec_work_ppp) .or. compute_quantity(advec_work_mpp)) Then
            !    ncorrect=ncorrect+1
            !    aforcepp_phi = ncorrect
            !Endif
        Endif


        compute_mean_mean = .false.
        If (compute_quantity(vm_grad_vm_r))     compute_mean_mean = .true.
        If (compute_quantity(vm_grad_vm_theta)) compute_mean_mean = .true.
        If (compute_quantity(vm_grad_vm_phi))   compute_mean_mean = .true.
        If (compute_quantity(advec_work_mmm))   compute_mean_mean = .true.


        If (compute_mean_mean .or. compute_quantity(advec_work_mmm)) Then


            If (compute_quantity(vm_grad_vm_r) .or. compute_quantity(advec_work_mmm)) Then
                ncorrect = ncorrect+1
                aforcemm_r = ncorrect
            Endif

            !If (compute_quantity(vm_grad_vm_theta) .or. compute_quantity(advec_work_mmm)) Then
            !    ncorrect = ncorrect+1
            !    aforcemm_theta = ncorrect
            !Endif

            !If (compute_quantity(vm_grad_vm_phi) .or. compute_quantity(samom_advec_mm) &
            !    .or. compute_quantity(advec_work_mmm)) Then

            !    ncorrect = ncorrect+1
            !    aforcemm_phi = ncorrect

            !Endif
        Endif


        !/////////////////////////////////////////////////////
        ! Viscous Force
        If (compute_quantity(viscous_force_r) .or. compute_quantity(visc_work) .or. &
            compute_quantity(viscous_mforce_r) ) Then
            ncorrect = ncorrect+1
            vforce_r = ncorrect
        Endif


        !///////////////////////////////////////////////////////
        ! Lorentz Force
        If (compute_quantity(j_cross_b_r) .or. compute_quantity(mag_work)) Then
            ncorrect = ncorrect+1
            lforce_r = ncorrect
        Endif

        If (compute_quantity(jm_cross_bm_r) .or. compute_quantity(mag_work_mmm)) Then
            ncorrect = ncorrect+1
            lforcemm_r = ncorrect
        Endif

        If (compute_quantity(jp_cross_bp_r) .or. compute_quantity(mag_work_ppp) &
            .or. compute_quantity(mag_work_mpp)) Then
            ncorrect = ncorrect+1
            lforcepp_r = ncorrect

        Endif

        !///////////////////////////////////////////////////////////////////////////////////
        If (ncorrect .gt. 0) Then
            Allocate(mean_3dbuffer(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max,1:ncorrect))
            Allocate(mean_ell0buffer(my_r%min:my_r%max,1:ncorrect))
        Endif
    End Subroutine Initialize_Mean_Correction

    Subroutine Finalize_Mean_Correction()
        Implicit None
        If (ncorrect .gt. 0) Then
            DeAllocate(mean_3dbuffer)
            DeAllocate(mean_ell0buffer)
        Endif
    End Subroutine Finalize_Mean_Correction

    Subroutine Mean_Correction(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Allocatable :: cbuffer(:,:,:,:)
        Integer :: r,k, t
        Logical :: compute_full_full = .false.
        Logical :: compute_fluct_fluct = .false.
        Logical :: compute_mean_mean = .false.
        logical :: compute_mean_correct =.false.
        Real*8 :: amp,del2u, estress
        Real*8, Allocatable :: mu_visc(:), dmudr(:), ovstheta(:), ovs2theta(:)

        compute_mean_correct = .false.
        !///////////////////////////////////////////////////
        ! Coriolis forces
        If (compute_quantity(coriolis_force_r) .or. &
            compute_quantity(coriolis_mforce_r) ) Then
            Call Compute_Coriolis_Force_r(buffer)
            DO_PSI
                mean_3dbuffer(PSI,cforce_r) = qty(PSI)
            END_DO
            compute_mean_correct=.true.
        Endif

        !If (compute_quantity(coriolis_force_theta) .or. &
        !    compute_quantity(coriolis_mforce_theta) ) Then
        !    Call Compute_Coriolis_Force_theta(buffer)
        !    DO_PSI
        !        mean_3dbuffer(PSI,cforce_theta) = qty(PSI)
        !    END_DO
        !Endif

        !If (compute_quantity(coriolis_force_phi) .or. &
        !    compute_quantity(coriolis_mforce_phi) ) Then
        !    Call Compute_Coriolis_Force_phi(buffer)
        !    DO_PSI
        !        mean_3dbuffer(PSI,cforce_phi) = qty(PSI)
        !    END_DO
        !Endif

        !//////////////////////////////////////////////
        ! Advective terms

        ! -- v v
        Allocate(cbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))

        compute_full_full = .false.
        If (compute_quantity(v_grad_v_r))     compute_full_full = .true.
        !If (compute_quantity(v_grad_v_theta)) compute_full_full = .true.
        !If (compute_quantity(v_grad_v_phi))   compute_full_full = .true.
        If (compute_quantity(advec_work))     compute_full_full = .true.



        If (compute_full_full ) Then

            Call ADotGradB(buffer,buffer,cbuffer,aindices=vindex,bindices=vindex)

            If (compute_quantity(v_grad_v_r) .or. compute_quantity(advec_work)) Then
                DO_PSI
                    mean_3dbuffer(PSI,aforce_r) = cbuffer(PSI,1)*ref%density(r)
                END_DO
            Endif

            !If (compute_quantity(v_grad_v_theta) .or. compute_quantity(advec_work)) Then
            !    DO_PSI
            !        mean_3dbuffer(PSI,aforce_theta) = cbuffer(PSI,2)*ref%density(r)
            !    END_DO
            !Endif

            !If (compute_quantity(v_grad_v_phi) .or. compute_quantity(advec_work)) Then
            !    DO_PSI
            !        mean_3dbuffer(PSI,aforce_phi) = cbuffer(PSI,3)*ref%density(r)
            !    END_DO
            !Endif
            compute_mean_correct=.true.
        Endif


        ! -- v' v'
        compute_fluct_fluct = .false.
        If (compute_quantity(vp_grad_vp_r))     compute_fluct_fluct = .true.
        If (compute_quantity(advec_work_ppp))     compute_fluct_fluct = .true.
        If (compute_quantity(advec_work_mpp))     compute_fluct_fluct = .true.

        If (compute_fluct_fluct) Then
            Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices=vindex,bindices=vindex)

            If (compute_quantity(vp_grad_vp_r) .or. compute_quantity(advec_work_ppp) &
                .or. compute_quantity(advec_work_mpp)) Then

                DO_PSI
                    mean_3dbuffer(PSI,aforcepp_r) = cbuffer(PSI,1)*ref%density(r)
                END_DO
            Endif

            !If (compute_quantity(vp_grad_vp_theta) .or. compute_quantity(advec_work_ppp) &
            !    .or. compute_quantity(advec_work_mpp)) Then
            !    DO_PSI
            !        mean_3dbuffer(PSI,aforcepp_theta) = cbuffer(PSI,2)*ref%density(r)
            !    END_DO
            !Endif

            !If (compute_quantity(vp_grad_vp_phi) .or. compute_quantity(samom_advec_pp) &
            !    .or. compute_quantity(advec_work_ppp) .or. compute_quantity(advec_work_mpp)) Then
            !    DO_PSI
            !        mean_3dbuffer(PSI,aforcepp_phi) = cbuffer(PSI,3)*ref%density(r)
            !    END_DO
            !Endif
            compute_mean_correct=.true.
        Endif

        ! -- <v> <v>
        compute_mean_mean = .false.
        If (compute_quantity(vm_grad_vm_r))     compute_mean_mean = .true.
        !If (compute_quantity(vm_grad_vm_theta)) compute_mean_mean = .true.
        !If (compute_quantity(vm_grad_vm_phi))   compute_mean_mean = .true.
        If (compute_quantity(advec_work_mmm))   compute_mean_mean = .true.


        If (compute_mean_mean .or. compute_quantity(advec_work_mmm)) Then


            If (compute_quantity(vm_grad_vm_r) .or. compute_quantity(advec_work_mmm)) Then
                DO_PSI
                    mean_3dbuffer(PSI,aforcemm_r) = cbuffer(PSI,1)*ref%density(r)
                END_DO
            Endif

            !If (compute_quantity(vm_grad_vm_theta) .or. compute_quantity(advec_work_mmm)) Then
            !    DO_PSI
            !        mean_3dbuffer(PSI,aforcemm_theta) = cbuffer(PSI,2)*ref%density(r)
            !    END_DO
            !Endif

            !If (compute_quantity(vm_grad_vm_phi) .or. compute_quantity(samom_advec_mm) &
            !    .or. compute_quantity(advec_work_mmm)) Then
            !    DO_PSI
            !        mean_3dbuffer(PSI,aforcemm_phi) = cbuffer(PSI,3)*ref%density(r)
            !    END_DO
            !Endif
            compute_mean_correct=.true.
        Endif


        DeAllocate(cbuffer)



        !///////////////////////////////////////////////////////////////////////
        ! Viscous Force

        Allocate(ovstheta(1:N_theta), ovs2theta(1:N_theta)) ! 1/sin; 1/sin^2
        ovstheta = 1.0d0/sintheta
        ovs2theta = 1.0d0/sin2theta

        !Compute the dynamic viscosity mu=rho*nu (nu is kinematic viscosity)
        Allocate(mu_visc(1:N_R), dmudr(1:N_R))

        mu_visc = ref%density*nu
        dmudr = mu_visc*(ref%dlnrho+dlnu)





        If (compute_quantity(viscous_force_r) .or. compute_quantity(visc_work) .or. &
            compute_quantity(viscous_mforce_r) ) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_r}
                del2u = DDBUFF(PSI,dvrdrdr)+Two_Over_R(r)*buffer(PSI,dvrdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvrdtdt)+cottheta(t)*buffer(PSI,dvrdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u-2.0d0*OneOverRsquared(r)*( &
                        buffer(PSI,vr) + &
                        buffer(PSI,dvtdt)+buffer(PSI,vtheta)*cottheta(t) + &
                        ovstheta(t)*buffer(PSI,dvpdp) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*(buffer(PSI,dvrdr)*ref%dlnrho(r)+ &
                        buffer(PSI,vr)*ref%d2lnrho(r) )



                ! Finally, add the piece due to the gradient of mu
                estress = buffer(PSI,dvrdr)-One_Third*buffer(PSI,vr)*ref%dlnrho(r)

                mean_3dbuffer(PSI,vforce_r) = 2.0d0*dmudr(r)*estress + mu_visc(r)*del2u


            END_DO

            compute_mean_correct=.true.
        Endif


        DeAllocate(mu_visc, dmudr)
        DeAllocate(ovstheta,ovs2theta)



        !//////////////////////////////////////////////////////////
        ! Lorentz Forces
        If (compute_quantity(j_cross_b_r) .or. compute_quantity(mag_work)) Then
            DO_PSI
                mean_3dbuffer(PSI, lforce_r) = (buffer(PSI,curlbtheta)*buffer(PSI,bphi)- &
                         & buffer(PSI,btheta)*buffer(PSI,curlbphi) ) *ref%Lorentz_Coeff
            END_DO
            compute_mean_correct=.true.
        Endif

        If (compute_quantity(jm_cross_bm_r) .or. compute_quantity(mag_work_mmm)) Then
            DO_PSI2
                mean_3dbuffer(1:n_phi,PSI2, lforcemm_r) = ( m0_values(PSI2,curlbtheta)*m0_values(PSI2,bphi)- &
                                  & m0_values(PSI2,btheta)*m0_values(PSI2,curlbphi) )*ref%Lorentz_Coeff
            END_DO2
            compute_mean_correct=.true.
        Endif

        If (compute_quantity(jp_cross_bp_r) .or. compute_quantity(mag_work_ppp) &
            .or. compute_quantity(mag_work_mpp)) Then

            DO_PSI
                mean_3dbuffer(PSI, lforcepp_r) = ( fbuffer(PSI,curlbtheta)*fbuffer(PSI,bphi)- &
                         & fbuffer(PSI,btheta)*fbuffer(PSI,curlbphi) )*ref%Lorentz_Coeff
            END_DO
            compute_mean_correct=.true.
        Endif




        !//////////////////////////////////////////////////////
        ! Perform the averaging
        if (compute_mean_correct) then
            !Write(6,*)'Allocated: ', allocated(mean_3dbuffer), allocated(mean_ell0buffer)

            Call ComputeEll0(mean_3dbuffer,mean_ell0buffer)
        endif
    End Subroutine Mean_Correction

    !/////////////////////////////////////////
    ! Coriolis Force Routines

    Subroutine Compute_Coriolis_Force_r(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8 :: coriolis_term
        Integer :: r,k, t
        coriolis_term = ref%Coriolis_Coeff
        DO_PSI
            qty(PSI) = ref%density(r)*coriolis_term*sintheta(t)*buffer(PSI,vphi)
        END_DO
    End Subroutine Compute_Coriolis_Force_r

    Subroutine Compute_Coriolis_Force_theta(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8 :: coriolis_term
        Integer :: r,k, t
        coriolis_term = ref%Coriolis_Coeff
        DO_PSI
            qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*buffer(PSI,vphi)
        END_DO
    End Subroutine Compute_Coriolis_Force_theta

    Subroutine Compute_Coriolis_Force_Phi(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8 :: coriolis_term
        Integer :: r,k, t
        coriolis_term = ref%Coriolis_Coeff
        DO_PSI
            qty(PSI) = - (coriolis_term*costheta(t)*buffer(PSI,vtheta) &
                       + coriolis_term*sintheta(t)*buffer(PSI,vr))*ref%density(r)
        END_DO
    End Subroutine Compute_Coriolis_Force_Phi

    !/////////////////////////////////////////
    ! Advection Routines
    !Subroutine Advection_Radial_Mean()

    !End Subroutine Advection_Radial_mean


    !//////////////////////////////////////////
    ! Lorentz Force routines


End Module Diagnostics_Mean_Correction
