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
Module Diagnostics_Linear_Forces
    Use Diagnostics_Base
    Implicit None

Contains

    Subroutine Compute_Linear_Forces(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Call Compute_Buoyancy_Force(buffer)
        Call Compute_Coriolis_Force(buffer)
        Call Compute_Viscous_Force(buffer)
        Call Compute_Pressure_Force(buffer)
    End Subroutine Compute_Linear_Forces

    Subroutine Compute_Buoyancy_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        ! Buoyancy forces as they contribute to the ell .ne. 0 components of the momentum equation

        ! NOTE:  The ell=0 component of the r-momentum equation is entirely described by
        !        hydrostatic balance between the pressure and entropy perturbations
        !        (the reference state is assumed to also be in hydrostatic balance).
        !        As such, the ell=0 buoancy force is uninteresting from the point of
        !        view of the flow.  We explicitly separate the ell=0 component for this
        !        term (as with the pressure term).
        ! -- full buoyancy
        If (compute_quantity(buoyancy_force) .or. compute_quantity(buoy_work)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*(buffer(PSI,tvar)-&
                           & ell0_values(r,tvar))
            END_DO
            If (compute_quantity(buoyancy_force)) Call Add_Quantity(qty)
            If (compute_quantity(buoy_work)) Then
                DO_PSI
                    qty(PSI)=buffer(PSI,vr)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(buoyancy_force_ell0)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*ell0_values(r,tvar)
            END_DO
            Call Add_Quantity(qty)

        Endif

        ! -- fluctuating buoyancy (ell = 0, m =0 already subtracted)
        If (compute_quantity(buoyancy_pforce) .or. compute_quantity(buoy_work_pp)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*fbuffer(PSI,tvar)
            END_DO
            If (compute_quantity(buoyancy_pforce)) Call Add_Quantity(qty)
            If (compute_quantity(buoy_work_pp)) Then
                DO_PSI
                    qty(PSI)=fbuffer(PSI,vr)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        ! -- mean buoyancy
        If (compute_quantity(buoyancy_mforce) .or. compute_quantity(buoy_work_mm)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*(m0_values(PSI2,tvar)-&
                           & ell0_values(r,tvar))
            END_DO
            If (compute_quantity(buoyancy_mforce)) Call Add_Quantity(qty)
            If (compute_quantity(buoy_work_mm)) Then
                DO_PSI
                    qty(PSI)=m0_values(PSI2,vr)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

    End Subroutine Compute_Buoyancy_Force

    Subroutine Compute_Coriolis_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: coriolis_term
        coriolis_term = ref%Coriolis_Coeff

        If(compute_quantity(Coriolis_Force_r)) Then
            DO_PSI
                qty(PSI) = mean_3dbuffer(PSI,cforce_r)-mean_ell0buffer(r,cforce_r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If(compute_quantity(Coriolis_pForce_r)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*sintheta(t)*fbuffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If(compute_quantity(Coriolis_mForce_r)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*sintheta(t)*m0_values(PSI2,vphi) &
                           - mean_ell0buffer(r,cforce_r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_Force_Theta)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*buffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_pForce_Theta)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*fbuffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_mForce_Theta)) Then
            DO_PSI
                qty(PSI) = ref%density(r)*coriolis_term*costheta(t)*m0_values(PSI2,vphi) ! &
                           !-mean_ell0buffer(r,cforce_theta)
            END_DO
            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(Coriolis_Force_Phi)) Then
            DO_PSI
                qty(PSI) = - (coriolis_term*costheta(t)*buffer(PSI,vtheta) &
                           + coriolis_term*sintheta(t)*buffer(PSI,vr))*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_pForce_Phi)) Then
            DO_PSI
                qty(PSI) = - (coriolis_term*costheta(t)*fbuffer(PSI,vtheta) &
                           + coriolis_term*sintheta(t)*fbuffer(PSI,vr))*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(Coriolis_mForce_Phi) .or. compute_quantity(samom_coriolis)) Then
            DO_PSI
                qty(PSI) = - (coriolis_term*costheta(t)*m0_values(PSI2,vtheta) &
                           + coriolis_term*sintheta(t)*m0_values(PSI2,vr))*ref%density(r)
            END_DO
            If (compute_quantity(Coriolis_mForce_Phi)) Call Add_Quantity(qty)
            If (compute_quantity(samom_coriolis)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*radius(r)*sintheta(t)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

    End Subroutine Compute_Coriolis_Force


    Subroutine Compute_Viscous_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: mypi, amp,del2u, estress
        Real*8, Allocatable :: mu_visc(:), dmudr(:), ovstheta(:), ovs2theta(:)

        Allocate(ovstheta(1:N_theta), ovs2theta(1:N_theta)) ! 1/sin; 1/sin^2
        ovstheta = 1.0d0/sintheta
        ovs2theta = 1.0d0/sin2theta

        !Compute the dynamic viscosity mu=rho*nu (nu is kinematic viscosity)
        Allocate(mu_visc(1:N_R), dmudr(1:N_R))

        mu_visc = ref%density*nu
        dmudr = mu_visc*(ref%dlnrho+dlnu)


        !////////////////////////////////////////////////////////

        ! r-direction; Full
        !If (compute_quantity(viscous_force_r) .or. compute_quantity(visc_work)) Then

        !    DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_r}
        !        del2u = DDBUFF(PSI,dvrdrdr)+Two_Over_R(r)*buffer(PSI,dvrdr)
        !        del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvrdtdt)+cottheta(t)*buffer(PSI,dvrdt))
        !        del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
        !        del2u = del2u-2.0d0*OneOverRsquared(r)*( &
        !                buffer(PSI,vr) + &
        !                buffer(PSI,dvtdt)+buffer(PSI,vtheta)*cottheta(t) + &
        !                ovstheta(t)*buffer(PSI,dvpdp) )

                ! Move onto the non del squared terms (compressibility term)
        !        del2u = del2u-One_Third*(buffer(PSI,dvrdr)*ref%dlnrho(r)+ &
        !                buffer(PSI,vr)*ref%d2lnrho(r) )



                ! Finally, add the piece due to the gradient of mu
        !        estress = buffer(PSI,dvrdr)-One_Third*buffer(PSI,vr)*ref%dlnrho(r)

        !        qty(PSI) = 2.0d0*dmudr(r)*estress + mu_visc(r)*del2u


        !    END_DO

        !    If (compute_quantity(viscous_force_r)) Call Add_Quantity(qty)
        !    If (compute_quantity(visc_work)) Then
        !        DO_PSI
        !            tmp1(PSI)=buffer(PSI,vr)*qty(PSI)
        !        END_DO
        !    Endif
        !Endif

        ! r-direction; Full
        If (compute_quantity(viscous_force_r) .or. compute_quantity(visc_work)) Then
            DO_PSI
                qty(PSI) = mean_3dbuffer(PSI,vforce_r)-mean_ell0buffer(r,vforce_r)
            END_DO
            If (compute_quantity(viscous_force_r)) Call Add_Quantity(qty)
            If (compute_quantity(visc_work)) Then
                DO_PSI
                    tmp1(PSI)=buffer(PSI,vr)*qty(PSI)
                END_DO
            Endif
        Endif

        !Theta-direction; Full
        If (compute_quantity(viscous_force_theta) .or. compute_quantity(visc_work)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_theta}
                del2u = DDBUFF(PSI,dvtdrdr)+Two_Over_R(r)*buffer(PSI,dvtdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvtdtdt)+cottheta(t)*buffer(PSI,dvtdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_theta
                del2u = del2u +OneOverRSquared(r)*( 2.0d0*buffer(PSI,dvrdt) - &
                        ovs2theta(t)*(   buffer(PSI,vtheta) + &
                        2.0d0*costheta(t)*buffer(PSI,dvpdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*buffer(PSI,dvrdt)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(buffer(PSI,dvrdt)-buffer(PSI,vtheta) ) +buffer(PSI,dvtdr)

                qty(PSI) = dmudr(r)*estress +mu_visc(r)*del2u
            END_DO

            If (compute_quantity(viscous_force_theta)) Call Add_Quantity(qty)
            If (compute_quantity(visc_work)) Then
                DO_PSI
                    tmp1(PSI)=tmp1(PSI)+buffer(PSI,vtheta)*qty(PSI)
                END_DO
            Endif
        Endif


        !Phi-direction
        If (compute_quantity(viscous_force_phi) .or. compute_quantity(visc_work)) Then

            DO_PSI
                del2u = DDBUFF(PSI,dvpdrdr)+Two_Over_R(r)*buffer(PSI,dvpdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvpdtdt)+cottheta(t)*buffer(PSI,dvpdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvpdpdp)*ovs2theta(t)

                ! Add geometric terms here
                ! del2u -
                !Add geometric terms to make this { Del^2{u} }_phi
                del2u = del2u +OneOverRSquared(r)*( 2.0d0*buffer(PSI,dvrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   buffer(PSI,vphi) - &
                        2.0d0*costheta(t)*buffer(PSI,dvtdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*ovstheta(t)*buffer(PSI,dvrdp)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(ovstheta(t)*buffer(PSI,dvrdp)-buffer(PSI,vphi) )+ &
                          buffer(PSI,dvpdr)

                qty(PSI) =dmudr(r)*estress + mu_visc(r)*del2u
            END_DO

            If (compute_quantity(viscous_force_phi)) Call Add_Quantity(qty)
            If (compute_quantity(visc_work)) Then
                DO_PSI
                    tmp1(PSI)=tmp1(PSI)+buffer(PSI,vphi)*qty(PSI)
                END_DO
                Call Add_Quantity(tmp1)
            Endif
        Endif


        !.............................................................

        ! r-direction; fluctuating
        If (compute_quantity(viscous_pforce_r) .or. compute_quantity(visc_work_pp)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_r}
                del2u = d2_fbuffer(PSI,dvrdrdr)+Two_Over_R(r)*fbuffer(PSI,dvrdr)
                del2u = del2u+OneOverRSquared(r)*(d2_fbuffer(PSI,dvrdtdt)+cottheta(t)*fbuffer(PSI,dvrdt))
                del2u = del2u+OneOverRSquared(r)*d2_fbuffer(PSI,dvrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u-2.0d0*OneOverRsquared(r)*( &
                        fbuffer(PSI,vr) + &
                        fbuffer(PSI,dvtdt)+fbuffer(PSI,vtheta)*cottheta(t) + &
                        ovstheta(t)*fbuffer(PSI,dvpdp) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*(fbuffer(PSI,dvrdr)*ref%dlnrho(r)+ &
                        fbuffer(PSI,vr)*ref%d2lnrho(r) )



                ! Finally, add the piece due to the gradient of mu
                estress = fbuffer(PSI,dvrdr)+One_Third*fbuffer(PSI,vr)*ref%dlnrho(r)

                qty(PSI) = 2.0d0*dmudr(r)*estress + mu_visc(r)*del2u


            END_DO

            If (Compute_quantity(viscous_pforce_r)) Call Add_Quantity(qty)
            If (compute_quantity(visc_work_pp)) Then
                DO_PSI
                    tmp1(PSI)=fbuffer(PSI,vr)*qty(PSI)
                END_DO
            Endif
        Endif

        !Theta-direction; Fluctuating
        If (compute_quantity(viscous_pforce_theta) .or. compute_quantity(visc_work_pp)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_theta}
                del2u = d2_fbuffer(PSI,dvtdrdr)+Two_Over_R(r)*fbuffer(PSI,dvtdr)
                del2u = del2u+OneOverRSquared(r)*(d2_fbuffer(PSI,dvtdtdt)+cottheta(t)*fbuffer(PSI,dvtdt))
                del2u = del2u+OneOverRSquared(r)*d2_fbuffer(PSI,dvtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_theta
                del2u = del2u +OneOverRSquared(r)*( 2.0d0*fbuffer(PSI,dvrdt) - &
                        ovs2theta(t)*(   fbuffer(PSI,vtheta) + &
                        2.0d0*costheta(t)*fbuffer(PSI,dvpdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*fbuffer(PSI,dvrdt)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(fbuffer(PSI,dvrdt)-fbuffer(PSI,vtheta) ) +fbuffer(PSI,dvtdr)

                qty(PSI) = dmudr(r)*estress +mu_visc(r)*del2u
            END_DO


            If (Compute_quantity(viscous_pforce_theta)) Call Add_Quantity(qty)
            If (compute_quantity(visc_work_pp)) Then
                DO_PSI
                    tmp1(PSI)=tmp1(PSI)+fbuffer(PSI,vtheta)*qty(PSI)
                END_DO
            Endif
        Endif

        !Phi-direction (fluctuating)
        If (compute_quantity(viscous_pforce_phi) .or. compute_quantity(visc_work_pp)) Then

            DO_PSI
                del2u = d2_fbuffer(PSI,dvpdrdr)+Two_Over_R(r)*fbuffer(PSI,dvpdr)
                del2u = del2u+OneOverRSquared(r)*(d2_fbuffer(PSI,dvpdtdt)+cottheta(t)*fbuffer(PSI,dvpdt))
                del2u = del2u+OneOverRSquared(r)*d2_fbuffer(PSI,dvpdpdp)*ovs2theta(t)

                ! Add geometric terms here
                ! del2u -
                !Add geometric terms to make this { Del^2{u} }_phi
                del2u = del2u +OneOverRSquared(r)*( 2.0d0*fbuffer(PSI,dvrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   fbuffer(PSI,vphi) - &
                        2.0d0*costheta(t)*fbuffer(PSI,dvtdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*ovstheta(t)*fbuffer(PSI,dvrdp)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(ovstheta(t)*fbuffer(PSI,dvrdp)-fbuffer(PSI,vphi) )+ &
                          fbuffer(PSI,dvpdr)

                qty(PSI) =dmudr(r)*estress + mu_visc(r)*del2u
            END_DO

            If (Compute_quantity(viscous_pforce_phi)) Call Add_Quantity(qty)
            If (compute_quantity(visc_work_pp)) Then
                DO_PSI
                    tmp1(PSI)=tmp1(PSI)+fbuffer(PSI,vphi)*qty(PSI)
                END_DO
                Call Add_Quantity(tmp1)
            Endif
        Endif

        ! r-direction; mean
        If (compute_quantity(viscous_mforce_r) .or. compute_quantity(visc_work_mm)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_r}
                del2u = d2_m0(PSI2,dvrdrdr)+Two_Over_R(r)*m0_values(PSI2,dvrdr)
                del2u = del2u+OneOverRSquared(r)*(d2_m0(PSI2,dvrdtdt)+cottheta(t)*m0_values(PSI2,dvrdt))
                del2u = del2u+OneOverRSquared(r)*d2_m0(PSI2,dvrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u-2.0d0*OneOverRsquared(r)*( &
                        m0_values(PSI2,vr) + &
                        m0_values(PSI2,dvtdt)+m0_values(PSI2,vtheta)*cottheta(t) + &
                        ovstheta(t)*m0_values(PSI2,dvpdp) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*(m0_values(PSI2,dvrdr)*ref%dlnrho(r)+ &
                        m0_values(PSI2,vr)*ref%d2lnrho(r) )



                ! Finally, add the piece due to the gradient of mu
                estress = m0_values(PSI2,dvrdr)+One_Third*m0_values(PSI2,vr)*ref%dlnrho(r)

                qty(PSI) = 2.0d0*dmudr(r)*estress + mu_visc(r)*del2u
                qty(PSI) = qty(PSI)-mean_ell0buffer(r,vforce_r)

            END_DO


            If (Compute_quantity(viscous_mforce_r)) Call Add_Quantity(qty)
            If (compute_quantity(visc_work_mm)) Then
                DO_PSI
                    tmp1(PSI)=m0_values(PSI2,vr)*qty(PSI)
                END_DO
            Endif
        Endif





        !Theta-direction; Mean
        If (compute_quantity(viscous_mforce_theta) .or. compute_quantity(visc_work_mm)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_theta}
                del2u = d2_m0(PSI2,dvtdrdr)+Two_Over_R(r)*m0_values(PSI2,dvtdr)
                del2u = del2u+OneOverRSquared(r)*(d2_m0(PSI2,dvtdtdt)+cottheta(t)*m0_values(PSI2,dvtdt))
                del2u = del2u+OneOverRSquared(r)*d2_m0(PSI2,dvtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_theta
                del2u = del2u +OneOverRSquared(r)*( 2.0d0*m0_values(PSI2,dvrdt) - &
                        ovs2theta(t)*(   m0_values(PSI2,vtheta) + &
                        2.0d0*costheta(t)*m0_values(PSI2,dvpdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*m0_values(PSI2,dvrdt)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(m0_values(PSI2,dvrdt)-m0_values(PSI2,vtheta) ) +m0_values(PSI2,dvtdr)

                qty(PSI) = dmudr(r)*estress +mu_visc(r)*del2u
            END_DO

            If (Compute_quantity(viscous_mforce_theta)) Call Add_Quantity(qty)
            If (compute_quantity(visc_work_mm)) Then
                DO_PSI
                    tmp1(PSI)=tmp1(PSI)+m0_values(PSI2,vtheta)*qty(PSI)
                END_DO
            Endif
        Endif





        !Phi-direction (mean)
        If (compute_quantity(viscous_mforce_phi) .or. compute_quantity(samom_diffusion) .or. compute_quantity(visc_work_mm)) Then

            DO_PSI
                del2u = d2_m0(PSI2,dvpdrdr)+Two_Over_R(r)*m0_values(PSI2,dvpdr)
                del2u = del2u+OneOverRSquared(r)*(d2_m0(PSI2,dvpdtdt)+cottheta(t)*m0_values(PSI2,dvpdt))
                del2u = del2u+OneOverRSquared(r)*d2_m0(PSI2,dvpdpdp)*ovs2theta(t)

                ! Add geometric terms here
                ! del2u -
                !Add geometric terms to make this { Del^2{u} }_phi
                del2u = del2u +OneOverRSquared(r)*( 2.0d0*m0_values(PSI2,dvrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   m0_values(PSI2,vphi) - &
                        2.0d0*costheta(t)*m0_values(PSI2,dvtdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*ovstheta(t)*m0_values(PSI2,dvrdp)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(ovstheta(t)*m0_values(PSI2,dvrdp)-m0_values(PSI2,vphi) )+ &
                          m0_values(PSI2,dvpdr)

                qty(PSI) =dmudr(r)*estress + mu_visc(r)*del2u
            END_DO

            If (compute_quantity(viscous_mforce_phi)) Call Add_Quantity(qty)
            If (compute_quantity(visc_work_mm)) Then
                DO_PSI
                    tmp1(PSI)=tmp1(PSI)+m0_values(PSI2,vphi)*qty(PSI)
                END_DO
                Call Add_Quantity(tmp1)
            Endif
            If (compute_quantity(samom_diffusion)) Then
                    DO_PSI
                        qty(PSI) = qty(PSI)*radius(r)*sintheta(t)
                    END_DO
                    Call Add_Quantity(qty)
            Endif
        Endif
        DeAllocate(mu_visc, dmudr)
        DeAllocate(ovstheta,ovs2theta)
    End Subroutine Compute_Viscous_Force

Subroutine Compute_Pressure_Force(buffer)
    Implicit None
    Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
    Real*8  :: pfactor(my_r%min:my_r%max)
    Integer :: r,k, t
    ! This routine essentially just compute the pressure gradient, but
    ! scaled by the reference pressure term.  This is typically either 1
    ! or 1/Ekman number.  We still treat it as a function of radius to
    ! keep things general.
    pfactor(my_r%min:my_r%max) = ref%dpdr_w_term(my_r%min:my_r%max) &
                                /ref%density(my_r%min:my_r%max)

    !////////////////////////////////////////
    !       Pressure Force


    ! NOTE:  The ell=0 component of the r-momentum equation is entirely described by
    !        hydrostatic balance between the pressure and entropy perturbations
    !        (the reference state is assumed to be in hydrostatic balance).
    !        As such, the ell=0 pressure force is uninteresting from the point of
    !        view of the flow.  We explicitly separate the ell=0 component for this
    !        term (as with the buoyancy term).

    ! pressure full
    ! r
    If (compute_quantity(pressure_force_r) .or. compute_quantity(press_work)) Then
        DO_PSI
            qty(PSI) = -(buffer(PSI,dpdr)-ell0_values(r,dpdr))*pfactor(r) + &
                        (buffer(PSI,pvar)-ell0_values(r,pvar))*pfactor(r) * &
                        ref%dlnrho(r)
        END_DO
        If (compute_quantity(pressure_force_r)) Call Add_Quantity(qty)
        If (compute_quantity(press_work)) Then
            DO_PSI
                tmp1(PSI)=qty(PSI)*buffer(PSI,vr)
            END_DO
        Endif
    Endif

    ! theta
    If (compute_quantity(pressure_force_theta) .or. compute_quantity(press_work)) Then
        DO_PSI
            qty(PSI) = -buffer(PSI,dpdt)*pfactor(r)*One_Over_R(r)
        END_DO
        If (compute_quantity(pressure_force_theta)) Call Add_Quantity(qty)
        If (compute_quantity(press_work)) Then
            DO_PSI
                tmp1(PSI)=tmp1(PSI)+qty(PSI)*buffer(PSI,vtheta)
            END_DO
        Endif
    Endif

    ! phi
    If (compute_quantity(pressure_force_phi) .or. compute_quantity(press_work)) Then
        DO_PSI
            qty(PSI) = -buffer(PSI,dpdp)*pfactor(r)*One_Over_R(r)*csctheta(t)
        END_DO
        If (compute_quantity(pressure_force_phi)) Call Add_Quantity(qty)
        If (compute_quantity(press_work)) Then
            DO_PSI
                tmp1(PSI)=tmp1(PSI)+qty(PSI)*buffer(PSI,vphi)
            END_DO
            Call Add_Quantity(tmp1)
        Endif
    Endif


    !fluctuating pressure
    ! r
    If (compute_quantity(pressure_pforce_r) .or. compute_quantity(press_work_pp)) Then
        DO_PSI
            qty(PSI) = -fbuffer(PSI,dpdr)*pfactor(r) + &
                        fbuffer(PSI,pvar)*pfactor(r) * &
                        ref%dlnrho(r)
        END_DO
        If (compute_quantity(pressure_pforce_r)) Call Add_Quantity(qty)
        If (compute_quantity(press_work_pp)) Then
            DO_PSI
                tmp1(PSI)=qty(PSI)*fbuffer(PSI,vr)
            END_DO
        Endif
    Endif

    ! theta
    If (compute_quantity(pressure_pforce_theta) .or. compute_quantity(press_work_pp)) Then
        DO_PSI
            qty(PSI) = -fbuffer(PSI,dpdt)*pfactor(r)*One_Over_R(r)
        END_DO
        If (compute_quantity(pressure_pforce_theta)) Call Add_Quantity(qty)
        If (compute_quantity(press_work_pp)) Then
            DO_PSI
                tmp1(PSI)=tmp1(PSI)+qty(PSI)*fbuffer(PSI,vtheta)
            END_DO
        Endif
    Endif

    ! phi
    If (compute_quantity(pressure_pforce_phi) .or. compute_quantity(press_work_pp)) Then
        DO_PSI
            qty(PSI) = -fbuffer(PSI,dpdp)*pfactor(r)*One_Over_R(r)*csctheta(t)
        END_DO
        If (compute_quantity(pressure_pforce_phi)) Call Add_Quantity(qty)
        If (compute_quantity(press_work_pp)) Then
            DO_PSI
                tmp1(PSI)=tmp1(PSI)+qty(PSI)*fbuffer(PSI,vphi)
            END_DO
            Call Add_Quantity(tmp1)
        Endif
    Endif


    ! Mean pressure

    ! r
    If (compute_quantity(pressure_mforce_r).or. compute_quantity(press_work_mm)) Then
        DO_PSI
            qty(PSI) = -( m0_values(PSI2,dpdr) - ell0_values(r,dpdr) )*pfactor(r) + &
                        ( m0_values(PSI2,pvar) - ell0_values(r,pvar) )*pfactor(r) * &
                        ref%dlnrho(r)
        END_DO

        If (compute_quantity(pressure_mforce_r)) Call Add_Quantity(qty)

        If (compute_quantity(press_work_mm)) Then
            DO_PSI
                tmp1(PSI)=qty(PSI)*m0_values(PSI2,vr)
            END_DO
        Endif
    Endif



    ! Theta
    If (compute_quantity(pressure_mforce_theta) .or. compute_quantity(press_work_mm)) Then
        DO_PSI
            qty(PSI) = -m0_values(PSI2,dpdt)*pfactor(r)*One_Over_R(r)
        END_DO
        If (compute_quantity(pressure_mforce_theta)) Call Add_Quantity(qty)

        If (compute_quantity(press_work_mm)) Then
            DO_PSI
                tmp1(PSI)=tmp1(PSI)+qty(PSI)*m0_values(PSI2,vtheta)
            END_DO
        Endif
    Endif

    ! Phi
    If (compute_quantity(pressure_mforce_phi) .or. compute_quantity(press_work_mm)) Then
        DO_PSI
            qty(PSI) = -m0_values(PSI2,dpdp)*pfactor(r)*One_Over_R(r)*csctheta(t)
        END_DO

        If (compute_quantity(pressure_mforce_phi)) Call Add_Quantity(qty)

        If (compute_quantity(press_work_mm)) Then
            DO_PSI
                tmp1(PSI)=tmp1(PSI)+qty(PSI)*m0_values(PSI2,vphi)
            END_DO
            Call Add_Quantity(tmp1)
        Endif

    Endif



    ! Spherically symmetric component
    If (compute_quantity(pressure_force_ell0_r)) Then
        DO_PSI
            qty(PSI) = -(ell0_values(r,dpdr))*pfactor(r) + &
                        (ell0_values(r,pvar))*pfactor(r) * &
                        ref%dlnrho(r)
        END_DO
        Call Add_Quantity(qty)
    Endif



End Subroutine Compute_Pressure_Force

End Module Diagnostics_Linear_Forces
