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
!////////////////////// Diagnostics KE Flux ///////////////////////
!
!       This module handles calculation of terms associated with the kinetic energy flux
!///////////////////////////////////////////////////////////////////
Module Diagnostics_KE_Flux
    Use Diagnostics_Base

Contains
    Subroutine Compute_KE_Flux(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: htmp1, htmp2, htmp3             ! temporary variables for use if needed
        Real*8 :: one_over_rsin, ctn_over_r        ! spherical trig
        Real*8 :: Err,Ett,Epp, Ert,Erp,Etp        ! variables to store the components of the rate of strain
        Real*8 :: Lap_r, Lap_t, Lap_p            ! variables to store Laplacians
        Real*8 :: mu, dmudr                ! the dynamic viscosity and its radial derivativ
        ! First, we compute the flux of total KE in each direction
        If (compute_quantity(ke_flux_radial) .or. compute_quantity(ke_flux_theta) &
            .or. compute_quantity(ke_flux_phi)) Then

            DO_PSI
                tmp1(PSI) = buffer(PSI,vr)**2 + buffer(PSI,vtheta)**2 + buffer(PSI,vphi)**2
            END_DO
            DO_PSI
                tmp1(PSI) = tmp1(PSI)*half*ref%density(r)
            END_DO

            If (compute_quantity(ke_flux_radial)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,vr)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ke_flux_theta)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,vtheta)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ke_flux_phi)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,vphi)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        ! Flux of mean KE due to mean flow
        If (compute_quantity(mke_mflux_radial) .or. compute_quantity(mke_mflux_theta) &
            .or. compute_quantity(mke_mflux_phi)) Then

            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vr)**2 + m0_values(PSI2,vtheta)**2 + m0_values(PSI2,vphi)**2
            END_DO
            DO_PSI
                tmp1(PSI) = tmp1(PSI)*half*ref%density(r)
            END_DO

            If (compute_quantity(mke_mflux_radial)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vr)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(mke_mflux_theta)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vtheta)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(mke_mflux_phi)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vphi)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        ! Flux of fluctuating KE due to mean flow
        If (compute_quantity(pke_mflux_radial) .or. compute_quantity(pke_mflux_theta) &
            .or. compute_quantity(pke_mflux_phi)) Then

            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr)**2 + fbuffer(PSI,vtheta)**2 + fbuffer(PSI,vphi)**2
            END_DO
            DO_PSI
                tmp1(PSI) = tmp1(PSI)*half*ref%density(r)
            END_DO

            If (compute_quantity(pke_mflux_radial)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vr)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(pke_mflux_theta)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vtheta)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(pke_mflux_phi)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vphi)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        ! Flux of fluctuating KE due to fluctuating flow
        If (compute_quantity(pke_pflux_radial) .or. compute_quantity(pke_pflux_theta) &
            .or. compute_quantity(pke_pflux_phi)) Then

            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr)**2 + fbuffer(PSI,vtheta)**2 + fbuffer(PSI,vphi)**2
            END_DO
            DO_PSI
                tmp1(PSI) = tmp1(PSI)*half*ref%density(r)
            END_DO

            If (compute_quantity(pke_pflux_radial)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,vr)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(pke_pflux_theta)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,vtheta)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(pke_pflux_phi)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,vphi)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif


        !/////////////////////////////////////////////////////////////////
        ! Now move on to the viscous fluxes
        !First, the full fluxes

        If (compute_quantity(visc_flux_r)) Then
            DO_PSI
                tmp1(PSI) = 2*buffer(PSI,vr) * buffer(PSI,dvrdr) ! 2 vr dot e_r_r
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,vtheta) * &              ! 2 vtheta dot e_r_theta
                ( buffer(PSI,dvtdr) -one_over_r(r)*(buffer(PSI,vtheta) - buffer(PSI,dvrdt)) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,vphi) * &              ! 2 vphi dot e_r_phi
                ( buffer(PSI,dvpdr) -one_over_r(r)*(buffer(PSI,vphi) - csctheta(t)*buffer(PSI,dvrdp)) )
            END_DO

            DO_PSI                                              ! -2 vr dot div(v)/3
                tmp1(PSI) = tmp1(PSI)+2*buffer(PSI,vr)*(buffer(PSI,vr)*ref%dlnrho(r)*one_third )
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*nu(r)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(visc_flux_theta)) Then
            DO_PSI
                tmp1(PSI) = buffer(PSI,vr) * & ! 2 vr dot e_r_theta
                ( buffer(PSI,dvtdr) -one_over_r(r)*(buffer(PSI,vtheta) - buffer(PSI,dvrdt)) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+2*buffer(PSI,vtheta) * & ! 2 vtheta dot e_theta_theta
                one_over_r(r)*(buffer(PSI,vr) + buffer(PSI,dvtdt))
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,vphi) * & ! 2 vphi dot e_theta_phi
                one_over_r(r)*(csctheta(t)*buffer(PSI,dvtdp) + buffer(PSI,dvpdt) - buffer(PSI,vphi)*cottheta(t) )
            END_DO

            DO_PSI                                              ! -2 vtheta dot div(v)/3
                tmp1(PSI) = tmp1(PSI)+2*buffer(PSI,vtheta)*(buffer(PSI,vr)*ref%dlnrho(r)*one_third )
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*nu(r)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)


        Endif

        If (compute_quantity(visc_flux_phi)) Then
            DO_PSI
                tmp1(PSI) = buffer(PSI,vr) * &  ! 2 vr dot e_r_phi
                ( buffer(PSI,dvpdr) -one_over_r(r)*(buffer(PSI,vphi) - csctheta(t)*buffer(PSI,dvrdp)) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,vtheta) * & ! 2 vtheta dot e_theta_phi
                 one_over_r(r)*(csctheta(t)*buffer(PSI,dvtdp) + buffer(PSI,dvpdt) - buffer(PSI,vphi)*cottheta(t) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,vphi) * & ! 2 vphi dot e_phi_phi
                Two_Over_R(r)*(buffer(PSI,vr) + cottheta(t)*buffer(PSI,vtheta) + csctheta(t)*buffer(PSI,dvpdp))
            END_DO

            DO_PSI                                              ! -2 vphi dot div(v)/3
                tmp1(PSI) = tmp1(PSI)+2*buffer(PSI,vphi)*(buffer(PSI,vr)*ref%dlnrho(r)*one_third )
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*nu(r)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)

        Endif

        !----
        If (compute_quantity(visc_fluxpp_r)) Then
            DO_PSI
                tmp1(PSI) = 2*fbuffer(PSI,vr) * fbuffer(PSI,dvrdr) ! 2 vr dot e_r_r
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+fbuffer(PSI,vtheta) * &              ! 2 vtheta dot e_r_theta
                ( fbuffer(PSI,dvtdr) -one_over_r(r)*(fbuffer(PSI,vtheta) - fbuffer(PSI,dvrdt)) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+fbuffer(PSI,vphi) * &              ! 2 vphi dot e_r_phi
                ( fbuffer(PSI,dvpdr) -one_over_r(r)*(fbuffer(PSI,vphi) - csctheta(t)*fbuffer(PSI,dvrdp)) )
            END_DO

            DO_PSI                                              ! -2 vr dot div(v)/3
                tmp1(PSI) = tmp1(PSI)+2*fbuffer(PSI,vr)*(fbuffer(PSI,vr)*ref%dlnrho(r)*one_third )
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*nu(r)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(visc_fluxpp_theta)) Then
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr) * & ! 2 vr dot e_r_theta
                ( fbuffer(PSI,dvtdr) -one_over_r(r)*(fbuffer(PSI,vtheta) - fbuffer(PSI,dvrdt)) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+2*fbuffer(PSI,vtheta) * & ! 2 vtheta dot e_theta_theta
                one_over_r(r)*(fbuffer(PSI,vr) + fbuffer(PSI,dvtdt))
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+fbuffer(PSI,vphi) * & ! 2 vphi dot e_theta_phi
                one_over_r(r)*(csctheta(t)*fbuffer(PSI,dvtdp) + fbuffer(PSI,dvpdt) - fbuffer(PSI,vphi)*cottheta(t) )
            END_DO

            DO_PSI                                              ! -2 vtheta dot div(v)/3
                tmp1(PSI) = tmp1(PSI)+2*fbuffer(PSI,vtheta)*(fbuffer(PSI,vr)*ref%dlnrho(r)*one_third )
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*nu(r)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)


        Endif

        If (compute_quantity(visc_fluxpp_phi)) Then
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr) * &  ! 2 vr dot e_r_phi
                ( fbuffer(PSI,dvpdr) -one_over_r(r)*(fbuffer(PSI,vphi) - csctheta(t)*fbuffer(PSI,dvrdp)) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+fbuffer(PSI,vtheta) * & ! 2 vtheta dot e_theta_phi
                 one_over_r(r)*(csctheta(t)*fbuffer(PSI,dvtdp) + fbuffer(PSI,dvpdt) - fbuffer(PSI,vphi)*cottheta(t) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+fbuffer(PSI,vphi) * & ! 2 vphi dot e_phi_phi
                Two_Over_R(r)*(fbuffer(PSI,vr) + cottheta(t)*fbuffer(PSI,vtheta) + csctheta(t)*fbuffer(PSI,dvpdp))
            END_DO

            DO_PSI                                              ! -2 vphi dot div(v)/3
                tmp1(PSI) = tmp1(PSI)+2*fbuffer(PSI,vphi)*(fbuffer(PSI,vr)*ref%dlnrho(r)*one_third )
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*nu(r)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)

        Endif



        !----
        If (compute_quantity(visc_fluxmm_r)) Then
            DO_PSI
                tmp1(PSI) = 2*m0_values(PSI2,vr) * m0_values(PSI2,dvrdr) ! 2 vr dot e_r_r
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+m0_values(PSI2,vtheta) * &              ! 2 vtheta dot e_r_theta
                ( m0_values(PSI2,dvtdr) -one_over_r(r)*(m0_values(PSI2,vtheta) - m0_values(PSI2,dvrdt)) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+m0_values(PSI2,vphi) * &              ! 2 vphi dot e_r_phi
                ( m0_values(PSI2,dvpdr) -one_over_r(r)*(m0_values(PSI2,vphi) - csctheta(t)*m0_values(PSI2,dvrdp)) )
            END_DO

            DO_PSI                                              ! -2 vr dot div(v)/3
                tmp1(PSI) = tmp1(PSI)+2*m0_values(PSI2,vr)*(m0_values(PSI2,vr)*ref%dlnrho(r)*one_third )
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*nu(r)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(visc_fluxmm_theta)) Then
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vr) * & ! 2 vr dot e_r_theta
                ( m0_values(PSI2,dvtdr) -one_over_r(r)*(m0_values(PSI2,vtheta) - m0_values(PSI2,dvrdt)) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+2*m0_values(PSI2,vtheta) * & ! 2 vtheta dot e_theta_theta
                one_over_r(r)*(m0_values(PSI2,vr) + m0_values(PSI2,dvtdt))
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+m0_values(PSI2,vphi) * & ! 2 vphi dot e_theta_phi
                one_over_r(r)*(csctheta(t)*m0_values(PSI2,dvtdp) + m0_values(PSI2,dvpdt) - m0_values(PSI2,vphi)*cottheta(t) )
            END_DO

            DO_PSI                                              ! -2 vtheta dot div(v)/3
                tmp1(PSI) = tmp1(PSI)+2*m0_values(PSI2,vtheta)*(m0_values(PSI2,vr)*ref%dlnrho(r)*one_third )
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*nu(r)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)


        Endif

        If (compute_quantity(visc_fluxmm_phi)) Then
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vr) * &  ! 2 vr dot e_r_phi
                ( m0_values(PSI2,dvpdr) -one_over_r(r)*(m0_values(PSI2,vphi) - csctheta(t)*m0_values(PSI2,dvrdp)) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+m0_values(PSI2,vtheta) * & ! 2 vtheta dot e_theta_phi
                 one_over_r(r)*(csctheta(t)*m0_values(PSI2,dvtdp) + m0_values(PSI2,dvpdt) - m0_values(PSI2,vphi)*cottheta(t) )
            END_DO

            DO_PSI
                tmp1(PSI) = tmp1(PSI)+m0_values(PSI2,vphi) * & ! 2 vphi dot e_phi_phi
                Two_Over_R(r)*(m0_values(PSI2,vr) + cottheta(t)*m0_values(PSI2,vtheta) + csctheta(t)*m0_values(PSI2,dvpdp))
            END_DO

            DO_PSI                                              ! -2 vphi dot div(v)/3
                tmp1(PSI) = tmp1(PSI)+2*m0_values(PSI2,vphi)*(m0_values(PSI2,vr)*ref%dlnrho(r)*one_third )
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*nu(r)*ref%density(r)
            END_DO

            Call Add_Quantity(qty)

        Endif


        !----



        ! Pressure transport terms
        If (compute_quantity(press_flux_r)) Then
            DO_PSI
                qty(PSI)=-buffer(PSI,vr)*buffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(press_flux_theta)) Then
            DO_PSI
                qty(PSI)=-buffer(PSI,vr)*buffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(press_flux_phi)) Then
            DO_PSI
                qty(PSI)=-buffer(PSI,vphi)*buffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(press_fluxpp_r)) Then
            DO_PSI
                qty(PSI)=-fbuffer(PSI,vr)*fbuffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(press_fluxpp_theta)) Then
            DO_PSI
                qty(PSI)=-fbuffer(PSI,vr)*fbuffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(press_fluxpp_phi)) Then
            DO_PSI
                qty(PSI)=-fbuffer(PSI,vphi)*fbuffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(press_fluxmm_r)) Then
            DO_PSI
                qty(PSI)=-m0_values(PSI2,vr)*m0_values(PSI2,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(press_fluxmm_theta)) Then
            DO_PSI
                qty(PSI)=-m0_values(PSI2,vr)*m0_values(PSI2,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(press_fluxmm_phi)) Then
            DO_PSI
                qty(PSI)=-m0_values(PSI2,vphi)*m0_values(PSI2,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Shear Production of turbulent kinetic energy.
        !       P_T = -rho_bar uu : E
        If (compute_quantity(production_shear_KE)) Then
            DO_PSI
                one_over_rsin = one_over_r(r) * csctheta(t)
                ctn_over_r = one_over_r(r) * cottheta(t)

                ! Compute elements of the mean rate of strain tensor E_ij
                Err = buffer(PSI,dvrdr)
                Ett = one_over_r(r) * (buffer(PSI,dvtdt) + buffer(PSI,vr))
                Epp = one_over_rsin * buffer(PSI,dvpdp) + ctn_over_r*buffer(PSI,vtheta)    &
                    + one_over_r(r) * buffer(PSI,vr)

                ! Twice the diagonal elements, e.g.,  Ert = 2 * E_rt
                Ert = one_over_r(r) * (buffer(PSI,dvrdt) - buffer(PSI,vtheta))        &
                    + buffer(PSI,dvtdr)
                Erp = buffer(PSI,dvpdr) + one_over_rsin * buffer(PSI,dvrdp)        &
                    - one_over_r(r) * buffer(PSI,vphi)
                Etp = one_over_rsin * buffer(PSI,dvtdp) - ctn_over_r * buffer(PSI,vphi)    &
                    + one_over_r(r) * buffer(PSI,dvpdt)

                ! Compute u'_i u'_j E_ij
                !DO k = 1, n_phi
                    ! Compute diagonal elements of the double contraction
                    htmp1 = buffer(PSI,vr)**2 * Err
                        htmp2 = buffer(PSI,vtheta)**2 * Ett
                    htmp3 = buffer(PSI,vphi)**2 * Epp
                            qty(PSI) = htmp1 + htmp2 + htmp3

                    ! Compute off-diagonal elements
                    htmp1 = buffer(PSI,vr)*buffer(PSI,vtheta)*Ert
                    htmp2 = buffer(PSI,vr)*buffer(PSI,vphi)*Erp
                    htmp3 = buffer(PSI,vtheta)*buffer(PSI,vphi)*Etp
                    qty(PSI) = qty(PSI) + htmp1 + htmp2 + htmp3

                    qty(PSI) = -ref%density(r) * qty(PSI)
                !ENDDO        ! End of phi loop
             END_DO       ! End of theta & r loop
                Call Add_Quantity(qty)
        Endif



        ! Shear Production of mean kinetic energy.
        !       P_T = -rho_bar <u><u> : E
        If (compute_quantity(production_shear_mKE)) Then
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)
                ctn_over_r = one_over_r(r) * cottheta(t)

                ! Compute elements of the mean rate of strain tensor E_ij
                Err = m0_values(PSI2,dvrdr)
                Ett = one_over_r(r) * (m0_values(PSI2,dvtdt) + m0_values(PSI2,vr))
                Epp = one_over_rsin * m0_values(PSI2,dvpdp) + ctn_over_r*m0_values(PSI2,vtheta)    &
                    + one_over_r(r) * m0_values(PSI2,vr)

                ! Twice the diagonal elements, e.g.,  Ert = 2 * E_rt
                Ert = one_over_r(r) * (m0_values(PSI2,dvrdt) - m0_values(PSI2,vtheta))        &
                    + m0_values(PSI2,dvtdr)
                Erp = m0_values(PSI2,dvpdr) + one_over_rsin * m0_values(PSI2,dvrdp)        &
                    - one_over_r(r) * m0_values(PSI2,vphi)
                Etp = one_over_rsin * m0_values(PSI2,dvtdp) - ctn_over_r * m0_values(PSI2,vphi)    &
                    + one_over_r(r) * m0_values(PSI2,dvpdt)

                ! Compute u'_i u'_j E_ij
                DO k = 1, n_phi
                    ! Compute diagonal elements of the double contraction
                    htmp1 = m0_values(PSI2,vr)**2 * Err
                        htmp2 = m0_values(PSI2,vtheta)**2 * Ett
                    htmp3 = m0_values(PSI2,vphi)**2 * Epp
                            qty(PSI) = htmp1 + htmp2 + htmp3

                    ! Compute off-diagonal elements
                    htmp1 = m0_values(PSI2,vr)*m0_values(PSI2,vtheta)*Ert
                    htmp2 = m0_values(PSI2,vr)*m0_values(PSI2,vphi)*Erp
                    htmp3 = m0_values(PSI2,vtheta)*m0_values(PSI2,vphi)*Etp
                    qty(PSI) = qty(PSI) + htmp1 + htmp2 + htmp3

                    qty(PSI) = -ref%density(r) * qty(PSI)
                ENDDO        ! End of phi loop
                END_DO2        ! End of theta & r loop
                Call Add_Quantity(qty)
        Endif
        ! Shear Production of turbulent kinetic energy.
        !       P_T = -rho_bar u'u' : <E>
        If (compute_quantity(production_shear_pKE)) Then
            DO_PSI2
                one_over_rsin = one_over_r(r) * csctheta(t)
                ctn_over_r = one_over_r(r) * cottheta(t)

                ! Compute elements of the mean rate of strain tensor E_ij
                Err = m0_values(PSI2,dvrdr)
                Ett = one_over_r(r) * (m0_values(PSI2,dvtdt) + m0_values(PSI2,vr))
                Epp = one_over_rsin * m0_values(PSI2,dvpdp) + ctn_over_r*m0_values(PSI2,vtheta)    &
                    + one_over_r(r) * m0_values(PSI2,vr)

                ! Twice the diagonal elements, e.g.,  Ert = 2 * E_rt
                Ert = one_over_r(r) * (m0_values(PSI2,dvrdt) - m0_values(PSI2,vtheta))        &
                    + m0_values(PSI2,dvtdr)
                Erp = m0_values(PSI2,dvpdr) + one_over_rsin * m0_values(PSI2,dvrdp)        &
                    - one_over_r(r) * m0_values(PSI2,vphi)
                Etp = one_over_rsin * m0_values(PSI2,dvtdp) - ctn_over_r * m0_values(PSI2,vphi)    &
                    + one_over_r(r) * m0_values(PSI2,dvpdt)

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
                    qty(PSI) = qty(PSI) + htmp1 + htmp2 + htmp3

                    qty(PSI) = -ref%density(r) * qty(PSI)
                ENDDO        ! End of phi loop
                END_DO2        ! End of theta & r loop
                Call Add_Quantity(qty)
        Endif


    End Subroutine Compute_KE_Flux

End Module Diagnostics_KE_Flux

