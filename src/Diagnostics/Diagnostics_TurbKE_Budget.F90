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
Module Diagnostics_TurbKE_Budget
    Use Diagnostics_Base
    Implicit None


Contains


    !=======================================================================================!
    !                                                                                       !
    !            Turbulent Kinetic Energy Budget                                            !
    !                                                                                       !
    !    The following code snippets compute terms (that once averaged in longitude)        !
    !   appear on the right hand side of the turbulent kinetic energy equation.             !
    !                                                                                       !
    !    (d/dt) [(rho_bar/2) <|u'|**2>] = B_T + P_T - Phi_T - div.F_T                       !
    !                                                                                       !
    !    F_T = F_TP + F_TV + F_MA + F_TA                Transport Flux                      !
    !                                                                                       !
    !    B_T = (g rho_bar / c_p) <S' w'>                Buoyant Production                  !
    !    P_T = -rho_bar <u'u'> : <e>                    Shear Production                    !
    !    Phi_T = 2 rho_bar nu * [<e':e'> - (1/3) <(div.u')**2>]    Viscous Dissipation      !
    !    F_TP = <P'u'> +                                Pressure Transport                  !
    !    F_TV = -<sigma'.u'>                            Viscous Transport                   !
    !    F_MA = (1/2) rho_bar <|u'|**2> <u>             Advection by the Mean               !
    !    F_TA = (1/2) rho_bar <|u'|**2 u'>              Turbulent Advection                 !
    !                                                                                       !
    !   The means in the Reynolds decomposition are longitudinal averages.                  !
    !                                                                                       !
    !   Note: The computation of the turbulent viscous transport requires MOST of the       !
    !    spatial second derivatives of the velocity components.  The ones needed are        !
    !    as follows (in Einstein notation for derivatives):                                 !
    !                                                                                       !
    !    v_r,rr   ;  v_r,tt   ;   v_r,pp   ;   v_r,rt   ;   v_r,rp                          !
    !    v_t,rr   ;  v_t,tt   ;   v_t,pp   ;   v_t,tr   ;   v_t,tp                          !
    !    v_p,rr   ;  v_p,tt   ;   v_p,pp   ;   v_p,pr   ;   v_r,pt                          !
    !                                                                                       !
    !   where r = radius, t = theta (colatitude), and p = phi (longitude). The missing      !
    !   combinations are few (v_r,tp ; v_t,rp ; v_p,rt).                                    !
    !      I presume that the indices will will be identified using the following scheme    !
    !   to access this information in the buffers:                                          !
    !                                                                                       !
    !    u_r,rr -- dvrdrdr ,  u_r,rt -- dvrdrdt and so on                                   !
    !                                                                                       !
    !   For cross derivatives, r precedes t which precedes p. So, its dvrdrdt NOT dvrdtdr.  !
    !                                                                                       !
    !=======================================================================================!

    Subroutine Compute_TurbulentKE_Budget(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        Real*8 :: htmp1, htmp2, htmp3             ! temporary variables for use if needed
        Real*8 :: one_over_rsin, ctn_over_r        ! spherical trig
        Real*8 :: Err,Ett,Epp, Ert,Erp,Etp        ! variables to store the components of the rate of strain
        Real*8 :: Lap_r, Lap_t, Lap_p            ! variables to store Laplacians
        Real*8 :: divu                          ! divergence of the velocity
        Real*8 :: mu, dmudr                ! the dynamic viscosity and its radial derivativ


        !-------------------------------------
        ! PRODUCTION AND DISSIPATION TERMS
        !-------------------------------------

        ! Buoyant Production of turbulent kinetic energy.
        !       B_T = (g rho_bar / c_p) S' u'_r
        If (compute_quantity(production_buoyant_pKE)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r)*fbuffer(PSI,tvar)*fbuffer(PSI,vr)    ! Assuming ref%Buoyancy_Coeff = g rho / cp
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Shear Production of turbulent kinetic energy.
        !       P_T = -rho_bar u'u' : <E>
        If (compute_quantity(production_shear2_pKE)) Then
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)
                ctn_over_r = one_over_r(r) * cottheta(t)

                ! Compute elements of the mean rate of strain tensor E_ij
                Err = m0_values(PSI2,dvrdr)
                Ett = one_over_r(r) * (m0_values(PSI2,dvtdt) + m0_values(PSI2,vr))
                Epp = one_over_rsin * m0_values(PSI2,dvpdp) + ctn_over_r*m0_values(PSI2,vtheta)    &
                    + one_over_r(r) * m0_values(PSI2,vr)

                ! The diagonal elements, e.g.,  Ert = E_rt
                Ert = one_over_r(r) * (m0_values(PSI2,dvrdt) - m0_values(PSI2,vtheta))        &
                    + m0_values(PSI2,dvtdr)
                Ert = 0.5D0*Ert
                Erp = m0_values(PSI2,dvpdr) + one_over_rsin * m0_values(PSI2,dvrdp)        &
                    - one_over_r(r) * m0_values(PSI2,vphi)
                Erp = 0.5D0*Erp
                Etp = one_over_rsin * m0_values(PSI2,dvtdp) - ctn_over_r * m0_values(PSI2,vphi)    &
                    + one_over_r(r) * m0_values(PSI2,dvpdt)
                Etp = 0.5D0*Etp

                ! Compute u'_i u'_j E_ij
                DO k = 1, n_phi
                    ! Compute diagonal elements of the double contraction
                    htmp1 = fbuffer(PSI,vr)**2 * Err
                    htmp2 = fbuffer(PSI,vtheta)**2 * Ett
                    htmp3 = fbuffer(PSI,vphi)**2 * Epp
                    qty(PSI) = htmp1 + htmp2 + htmp3

                    ! Compute off-diagonal elements
                    htmp1 = fbuffer(PSI,vr)*fbuffer(PSI,vtheta)*Ert
                    htmp2 = fbuffer(PSI,vr)*fbuffer(PSI,vphi)*Erp
                    htmp3 = fbuffer(PSI,vtheta)*fbuffer(PSI,vphi)*Etp
                    qty(PSI) = qty(PSI) + 2D0*(htmp1 + htmp2 + htmp3)

                    qty(PSI) = -ref%density(r) * qty(PSI)
                ENDDO        ! End of phi loop
            END_DO2        ! End of theta & r loop
            Call Add_Quantity(qty)
        Endif


        !Viscous Dissipation of turbulent kinetic energy.
        !    Phi_T = 2 rho_bar nu [ e': e' - (1/3) (div.u')**2 ]
        If (compute_quantity(dissipation_viscous_pKE)) Then
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)
                ctn_over_r = one_over_r(r) * cottheta(t)
                mu = nu(r) * ref%density(r)
                DO k = 1, n_phi
                    ! Compute elements of the turbulent rate of strain tensor e'_ij
                    err = fbuffer(PSI,dvrdr)
                    ett = one_over_r(r) * (fbuffer(PSI,dvtdt) + fbuffer(PSI,vr))
                    epp = one_over_rsin * fbuffer(PSI,dvpdp) + ctn_over_r * fbuffer(PSI,vtheta)    &
                        + one_over_r(r) * fbuffer(PSI,vr)

                    ! The diagonal elements of the turbulent rate of strain tensor, e.g.,  ert = e'_rt
                    ert = one_over_r(r) * (fbuffer(PSI,dvrdt) - fbuffer(PSI,vtheta))        &
                        + fbuffer(PSI,dvtdr)
                    ert = 0.5D0 * ert 
                    erp = fbuffer(PSI,dvpdr) + one_over_rsin * fbuffer(PSI,dvrdp)        &
                        - one_over_r(r) * fbuffer(PSI,vphi)
                    erp = 0.5D0 * erp 
                    etp = one_over_rsin * fbuffer(PSI,dvtdp) - ctn_over_r * fbuffer(PSI,vphi)    &
                        + one_over_r(r) * fbuffer(PSI,dvpdt)
                    etp = 0.5D0 * etp 

                    ! First compute e'_ij e'_ij
                    qty(PSI) = err*err + ett*ett + epp*epp            ! Diagonal: e'_ii e'_ii
                    qty(PSI) = qty(PSI) + 2D0*(ert*ert + erp*erp + etp*etp)    ! + Off-Diagonal

                    ! Next add -(1/3) (div.u)^2 & multiply by 2 rho_bar nu
                    divu = -ref%dlnrho(r) * fbuffer(PSI,vr)               ! Assume anelasticity
                    qty(PSI) = 2D0*mu * (qty(PSI) - one_third*divu**2)    ! Turbulent Viscous Dissipation        
                ENDDO        ! End of phi loop
            END_DO2        ! End of theta & r loop
            Call Add_Quantity(qty)
        Endif


        !-------------------------------------
        ! TRANSPORT TERMS
        !-------------------------------------

        ! Pressure Transport of turbulent kinetic energy.
        !    -div.F_TP = -div . (P' u')
        If (compute_quantity(transport_pressure_pKE)) Then
            DO_PSI
                htmp1 = fbuffer(PSI,vr)*fbuffer(PSI,dpdr)
                htmp2 = fbuffer(PSI,vtheta)*fbuffer(PSI,dpdt)*one_over_r(r)
                htmp3 = fbuffer(PSI,vphi)*fbuffer(PSI,dpdp)*one_over_r(r)/sintheta(t)

                qty(PSI) = -(htmp1 + htmp2 + htmp3)
                qty(PSI) = qty(PSI) + ref%dlnrho(r)*fbuffer(PSI,pvar)*fbuffer(PSI,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Viscous Transport of turbulent kinetic energy.
        !    -div.F_TV = div.(sigma'.u') = (div.sigma').u' + Phi_T
        If (compute_quantity(transport_viscous_pKE)) Then
            !Write(6,*)dvrdrdr, shape(d2_fbuffer)
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)
                ctn_over_r = one_over_r(r) * cottheta(t)
                mu = nu(r) * ref%density(r)
                dmudr =  mu * (ref%dlnrho(r) + dlnu(r))

                DO k = 1, n_phi
                    ! Compute elements of the turbulent rate of strain tensor e'_ij
                    err = fbuffer(PSI,dvrdr)
                    ett = one_over_r(r) * (fbuffer(PSI,dvtdt) + fbuffer(PSI,vr))
                    epp = one_over_rsin * fbuffer(PSI,dvpdp) + ctn_over_r * fbuffer(PSI,vtheta)    &
                        + one_over_r(r) * fbuffer(PSI,vr)

                    ! The diagonal elements, e.g.,  ert = e'_rt
                    ert = one_over_r(r) * (fbuffer(PSI,dvrdt) - fbuffer(PSI,vtheta))        &
                        + fbuffer(PSI,dvtdr)
                    ert = 0.5D0 * ert
                    erp = fbuffer(PSI,dvpdr) + one_over_rsin * fbuffer(PSI,dvrdp)        &
                        - one_over_r(r) * fbuffer(PSI,vphi)
                    erp = 0.5D0 * erp
                    etp = one_over_rsin * fbuffer(PSI,dvtdp) - ctn_over_r * fbuffer(PSI,vphi)    &
                        + one_over_r(r) * fbuffer(PSI,dvpdt)
                    etp = 0.5D0 * etp

                    ! First compute e'_ij e'_ij
                    qty(PSI) = err*err + ett*ett + epp*epp            ! Diagonal: e'_ii e'_ii
                    qty(PSI) = qty(PSI) + 2D0*(ert*ert + erp*erp + etp*etp)    ! + Off-Diagonal

                    ! Next add -(1/3) (div.u)^2 & multiply by 2 rho_bar nu
                    divu = -ref%dlnrho(r) * fbuffer(PSI,vr)               ! Assume anelasticity
                    qty(PSI) = 2D0*mu * (qty(PSI) - one_third*divu**2)    ! Turbulent Viscous Dissipation


                    ! Compute the Laplacians of the velocity components
                    htmp1 = d2_fbuffer(PSI,dvrdrdr) + 2D0*one_over_r(r) * fbuffer(PSI,dvrdr)
                    htmp2 = one_over_r(r)**2 * (d2_fbuffer(PSI,dvrdtdt) + cottheta(t)*fbuffer(PSI,dvrdt))
                    htmp3 = one_over_rsin**2 * d2_fbuffer(PSI,dvrdpdp)
                    Lap_r = htmp1 + htmp2 + htmp3                    ! Lap u'_r

                    htmp1 = d2_fbuffer(PSI,dvtdrdr) + 2D0*one_over_r(r) * fbuffer(PSI,dvtdr)
                    htmp2 = one_over_r(r)**2 * (d2_fbuffer(PSI,dvtdtdt) + cottheta(t)*fbuffer(PSI,dvtdt))
                    htmp3 = one_over_rsin**2 * d2_fbuffer(PSI,dvtdpdp)
                    Lap_t = htmp1 + htmp2 + htmp3                    ! Lap u'_theta

                    htmp1 = d2_fbuffer(PSI,dvpdrdr) + 2D0*one_over_r(r) * fbuffer(PSI,dvpdr)
                    htmp2 = one_over_r(r)**2 * (d2_fbuffer(PSI,dvpdtdt) + cottheta(t)*fbuffer(PSI,dvpdt))
                    htmp3 = one_over_rsin**2 * d2_fbuffer(PSI,dvpdpdp)
                    Lap_p = htmp1 + htmp2 + htmp3                    ! Lap u'_phi


                    ! Compute vector Laplacian u'
                    htmp1 = fbuffer(PSI,vr) + fbuffer(PSI,dvtdt) + cottheta(t) * fbuffer(PSI,vtheta)
                    htmp2 = one_over_r(r)*one_over_rsin * fbuffer(PSI,dvpdp)
                    Lap_r = Lap_r - 2D0*(one_over_r(r)**2 * htmp1 + htmp2)        ! r component (Lap u')_r

                    htmp1 = one_over_rsin * fbuffer(PSI,vtheta) + 2D0*ctn_over_r * fbuffer(PSI,dvpdp)
                    htmp2 = 2D0*one_over_r(r)**2 * fbuffer(PSI,dvrdt)
                    Lap_t = Lap_t + htmp2 - one_over_rsin * htmp1            ! theta component (Lap u')_theta

                    htmp1 = 2D0*ctn_over_r * fbuffer(PSI,dvtdp) - one_over_rsin * fbuffer(PSI,vphi)
                    htmp2 = 2D0*one_over_r(r) * fbuffer(PSI,dvrdp)
                    Lap_p = Lap_p + one_over_rsin * (htmp1 + htmp2)            ! phi component (Lap u')_phi

         
                    ! Compute  grad (div.u') using anelasticity (thus avoiding the need for second derivatives)
                    htmp1 = -ref%dlnrho(r) * fbuffer(PSI,dvrdr) - ref%d2lnrho(r) * fbuffer(PSI,vr)  ! [Grad (div.u')]_r
                    htmp2 = -ref%dlnrho(r)*one_over_r(r) * fbuffer(PSI,dvrdt) ! [Grad (div.u')]_theta
                    htmp3 = -ref%dlnrho(r)*one_over_rsin * fbuffer(PSI,dvrdp) ! [Grad (div.u')]_phi 


                    ! Compute viscous force: div . sigma'
                    htmp1 = 2D0*dmudr*(err - one_third*divu) + mu * (Lap_r + one_third*htmp1) ! (div.sigma')_r
                    htmp2 = 2D0*dmudr*ert + mu*(Lap_t + one_third*htmp2)                ! (div.sigma')_theta
                    htmp3 = 2D0*dmudr*erp + mu*(Lap_p + one_third*htmp3)                ! (div.sigma')_phi

                    qty(PSI) = qty(PSI) + htmp1*fbuffer(PSI,vr) + htmp2*fbuffer(PSI,vtheta)
                    qty(PSI) = qty(PSI) + htmp3*fbuffer(PSI,vphi) ! + (div.sigma).u'
                ENDDO        ! End of phi loop
            END_DO2        ! End of r and theta loop
            Call Add_Quantity(qty)
        Endif


        ! Turbulent Advective Transport of turbulent kinetic energy.
        !    -div.F_TA = -div . (1/2 rho_bar |u'|**2 u' ) = -(1/2) rho_bar (u'. grad) |u'|**2 = -rho_bar u'_i (u'. grad) u'_i
        If (compute_quantity(transport_turbadvect_pKE)) Then
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)

                DO k = 1, n_phi
                    ! Compute u'_r (u'.grad) u'_r
                        htmp1 = fbuffer(PSI,vr) * fbuffer(PSI,dvrdr)
                    htmp2 = one_over_r(r) * fbuffer(PSI,vtheta) * fbuffer(PSI,dvrdt)
                    htmp3 = one_over_rsin * fbuffer(PSI,vphi) * fbuffer(PSI,dvrdp)
                        qty(PSI) = fbuffer(PSI,vr)*(htmp1 + htmp2 + htmp3)        ! u'_r (u'. grad) u'_r

                    ! Compute u'_t (u'.grad) u'_t
                        htmp1 = fbuffer(PSI,vr) * fbuffer(PSI,dvtdr)
                    htmp2 = one_over_r(r) * fbuffer(PSI,vtheta) * fbuffer(PSI,dvtdt)
                    htmp3 = one_over_rsin * fbuffer(PSI,vphi) * fbuffer(PSI,dvtdp)
                        qty(PSI) = qty(PSI) + fbuffer(PSI,vtheta)*(htmp1+htmp2+htmp3)    ! +  u'_t (u'. grad) u'_t

                    ! Compute u'_p (u'.grad) u'_p
                        htmp1 = fbuffer(PSI,vr) * fbuffer(PSI,dvpdr)
                    htmp2 = one_over_r(r) * fbuffer(PSI,vtheta) * fbuffer(PSI,dvpdt)
                    htmp3 = one_over_rsin * fbuffer(PSI,vphi) * fbuffer(PSI,dvpdp)
                        qty(PSI) = qty(PSI) + fbuffer(PSI,vphi)*(htmp1+htmp2+htmp3)    ! +  u'_p (u'. grad) u'_p

                    qty(PSI) = -ref%density(r)*qty(PSI)
                ENDDO        ! End phi loop
            END_DO2        ! End r and theta loops
            Call Add_Quantity(qty)
        Endif


        ! Mean Advective Transport of turbulent kinetic energy.
        !    -div.F_MA = -div . (1/2 rho_bar |u'|**2 U ) = -(1/2) rho_bar (U.grad) |u'|**2 = -rho_bar u'_i (U . grad) u'_i
        If (compute_quantity(transport_meanadvect_pKE)) Then
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)

                DO k = 1, n_phi
                    ! Compute u'_r (U.grad) u'_r
                        htmp1 = m0_values(PSI2,vr) * fbuffer(PSI,dvrdr)
                    htmp2 = one_over_r(r) * m0_values(PSI2,vtheta) * fbuffer(PSI,dvrdt)
                    htmp3 = one_over_rsin * m0_values(PSI2,vphi) * fbuffer(PSI,dvrdp)
                        qty(PSI) = fbuffer(PSI,vr)*(htmp1+htmp2+htmp3)        ! u'_r (U . grad) u'_r

                    ! Compute u'_t (U.grad) u'_t
                        htmp1 = m0_values(PSI2,vr) * fbuffer(PSI,dvtdr)
                    htmp2 = one_over_r(r) * m0_values(PSI2,vtheta) * fbuffer(PSI,dvtdt)
                    htmp3 = one_over_rsin * m0_values(PSI2,vphi) * fbuffer(PSI,dvtdp)
                        qty(PSI) = qty(PSI) + fbuffer(PSI,vtheta)*(htmp1+htmp2+htmp3)    ! +  u'_t (U . grad) u'_t

                    ! Compute u'_p (U.grad) u'_p
                        htmp1 = m0_values(PSI2,vr) * fbuffer(PSI,dvpdr)
                    htmp2 = one_over_r(r) * m0_values(PSI2,vtheta) * fbuffer(PSI,dvpdt)
                    htmp3 = one_over_rsin * m0_values(PSI2,vphi) * fbuffer(PSI,dvpdp)
                        qty(PSI) = qty(PSI) + fbuffer(PSI,vphi)*(htmp1+htmp2+htmp3)    ! +  u'_p (U . grad) u'_p

                    qty(PSI) = -ref%density(r)*qty(PSI)
                ENDDO        ! End phi loop
            END_DO2        ! End r and theta loops
            Call Add_Quantity(qty)
        Endif


        !-------------------------------------
        ! FLUXES
        !-------------------------------------

        ! Radial Pressure Flux of turbulent kinetic energy.
        !    r_hat . F_TP = P' u'_r
        If (compute_quantity(rflux_pressure_pKE)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr) * fbuffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Colatitudinal Pressure Flux of turbulent kinetic energy.
        !    theta_hat . F_TP = P' u'_theta
        If (compute_quantity(thetaflux_pressure_pKE)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vtheta) * fbuffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Azimuthal Pressure Flux of turbulent kinetic energy.
        !    phi_hat . F_TP = P' u'_phi
        If (compute_quantity(phiflux_pressure_pKE)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vphi) * fbuffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Radial Viscous Flux of turbulent kinetic energy.
        !    r_hat . F_TV = -r_hat . sigma'.u'
        If (compute_quantity(rflux_viscous_pKE)) Then
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)
                ctn_over_r = one_over_r(r) * cottheta(t)
                    mu = nu(r) * ref%density(r)

                DO k = 1, n_phi
                    ! Compute elements of the turbulent rate of strain tensor e'_ij
                    err = fbuffer(PSI,dvrdr)

                    ! The diagonal elements, e.g.,  ert = e'_rt
                    ert = one_over_r(r) * (fbuffer(PSI,dvrdt) - fbuffer(PSI,vtheta))        &
                        + fbuffer(PSI,dvtdr)
                    ert = 0.5D0 * ert
                    erp = fbuffer(PSI,dvpdr) + one_over_rsin * fbuffer(PSI,dvrdp)        &
                        - one_over_r(r) * fbuffer(PSI,vphi)
                    erp = 0.5D0 * erp

                    ! Radial component of the viscous stess contracted w/ the velocity: r_hat . sigma' . u'
                    divu = -ref%dlnrho(r) * fbuffer(PSI,vr)                       ! Assume anelasticity
                    htmp1 = 2D0 * (err*fbuffer(PSI,vr) + ert*fbuffer(PSI,vtheta) + erp*fbuffer(PSI,vphi))
                    htmp2 = 2D0 * one_third * divu * fbuffer(PSI,vr)
                    qty(PSI) = mu * (htmp2 - htmp1)

                ENDDO
            END_DO2
            Call Add_Quantity(qty)
        Endif


        ! Colatitudinal Viscous Flux of turbulent kinetic energy.
        !    theta_hat . F_TV = -theta_hat . sigma'.u'
        If (compute_quantity(thetaflux_viscous_pKE)) Then
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)
                ctn_over_r = one_over_r(r) * cottheta(t)
                    mu = nu(r) * ref%density(r)

                DO k = 1, n_phi
                    ! Compute elements of the turbulent rate of strain tensor e'_ij
                    err = fbuffer(PSI,dvrdr)
                    ett = one_over_r(r) * (fbuffer(PSI,dvtdt) + fbuffer(PSI,vr))

                    ! The diagonal elements, e.g.,  ert = e'_rt
                    ert = one_over_r(r) * (fbuffer(PSI,dvrdt) - fbuffer(PSI,vtheta))        &
                        + fbuffer(PSI,dvtdr)
                    ert = 0.5D0 * ert
                    etp = one_over_rsin * fbuffer(PSI,dvtdp) - ctn_over_r * fbuffer(PSI,vphi)    &
                        + one_over_r(r) * fbuffer(PSI,dvpdt)
                    etp = 0.5D0 * etp

                    ! Colatitudinal component of the viscous stess contracted w/ the velocity: theta_hat . sigma' . u'
                    divu = -ref%dlnrho(r) * fbuffer(PSI,vr)                     ! Assume anelasticity
                    htmp1 = 2D0 * (ett*fbuffer(PSI,vtheta) + ert*fbuffer(PSI,vr) + etp*fbuffer(PSI,vphi))
                    htmp2 = 2D0 * one_third * divu * fbuffer(PSI,vtheta)
                    qty(PSI) = mu * (htmp2 - htmp1)

                ENDDO
            END_DO2
            Call Add_Quantity(qty)
        Endif


        ! Azimuthal Viscous Flux of turbulent kinetic energy.
        !    phi_hat . F_TV = -phi_hat . sigma'.u'
        If (compute_quantity(phiflux_viscous_pKE)) Then
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)
                ctn_over_r = one_over_r(r) * cottheta(t)
                mu = nu(r) * ref%density(r)

                DO k = 1, n_phi
                    ! Compute elements of the turbulent rate of strain tensor
                    ! e'_ij
                    erp = 0.5D0*(one_over_rsin*fbuffer(PSI,dvrdp) - one_over_r(r)*fbuffer(PSI,vphi)     &
                        + fbuffer(PSI,dvpdr))
                    etp = 0.5D0*(one_over_rsin*fbuffer(PSI,dvtdp) - ctn_over_r*fbuffer(PSI,vphi)        &
                        + one_over_r(r)*fbuffer(PSI,dvpdt))
                    epp = one_over_rsin*fbuffer(PSI,dvpdp) + one_over_r(r)*fbuffer(PSI,vr)              &
                        + ctn_over_r*fbuffer(PSI,vtheta)

                    ! Zonal component of the viscous stess contracted w/ the velocity: phi_hat . sigma' . u'
                    divu = -ref%dlnrho(r) * fbuffer(PSI,vr)                       ! Assume anelasticity
                    htmp1 = erp*fbuffer(PSI,vr) + etp*fbuffer(PSI,vtheta) + epp*fbuffer(PSI,vphi)
                    htmp2 = one_third * divu * fbuffer(PSI,vphi)
                    qty(PSI) = 2D0*mu * (htmp2 - htmp1)

                ENDDO
            END_DO2
            Call Add_Quantity(qty)
        Endif


        ! Radial Turbulent Advective Flux of turbulent kinetic energy.
        !    r_hat . F_TA = 1/2 rho_bar |u'|**2 u'_r
        If (compute_quantity(rflux_turbadvect_pKE)) Then
            DO_PSI
                htmp1 = fbuffer(PSI,vr)**2 + fbuffer(PSI,vtheta)**2 + fbuffer(PSI,vphi)**2
                qty(PSI) = 0.5D0*ref%density(r) * htmp1 * fbuffer(PSI,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Colatitudinal Turbulent Advective Flux of turbulent kinetic energy.
        !    theta_hat . F_TA = 1/2 rho_bar |u'|**2 u'_theta
        If (compute_quantity(thetaflux_turbadvect_pKE)) Then
            DO_PSI
                htmp1 = fbuffer(PSI,vr)**2 + fbuffer(PSI,vtheta)**2 + fbuffer(PSI,vphi)**2
                qty(PSI) = 0.5D0*ref%density(r) * htmp1 * fbuffer(PSI,vtheta)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Azimuthal Turbulent Advective Flux of turbulent kinetic energy.
        !    phi_hat . F_TA = 1/2 rho_bar |u'|**2 u'_phi
        If (compute_quantity(phiflux_turbadvect_pKE)) Then
            DO_PSI
                htmp1 = fbuffer(PSI,vr)**2 + fbuffer(PSI,vtheta)**2 + fbuffer(PSI,vphi)**2
                qty(PSI) = 0.5D0*ref%density(r) * htmp1 * fbuffer(PSI,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Radial Mean Advective Flux of turbulent kinetic energy.
        !    r_hat . F_MA = 1/2 rho_bar |u'|**2 U_r
        If (compute_quantity(rflux_meanadvect_pKE)) Then
            DO_PSI
                htmp1 = fbuffer(PSI,vr)**2 + fbuffer(PSI,vtheta)**2 + fbuffer(PSI,vphi)**2
                qty(PSI) = 0.5D0*ref%density(r) * htmp1 * m0_values(PSI2,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Colatitudinal Mean Advective Flux of turbulent kinetic energy.
        !    theta_hat . F_MA = 1/2 rho_bar |u'|**2 U_theta
        If (compute_quantity(thetaflux_meanadvect_pKE)) Then
            DO_PSI
                htmp1 = fbuffer(PSI,vr)**2 + fbuffer(PSI,vtheta)**2 + fbuffer(PSI,vphi)**2
                qty(PSI) = 0.5D0*ref%density(r) * htmp1 * m0_values(PSI2,vtheta)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Azimuthal Mean Advective Flux of turbulent kinetic energy.
        !    phi_hat . F_MA = 1/2 rho_bar |u'|**2 U_phi
        If (compute_quantity(phiflux_meanadvect_pKE)) Then
            DO_PSI
                htmp1 = fbuffer(PSI,vr)**2 + fbuffer(PSI,vtheta)**2 + fbuffer(PSI,vphi)**2
                qty(PSI) = 0.5D0*ref%density(r) * htmp1 * m0_values(PSI2,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_TurbulentKE_Budget

End Module Diagnostics_TurbKE_Budget
