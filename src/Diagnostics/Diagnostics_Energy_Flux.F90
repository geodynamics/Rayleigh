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
!////////////////////// Diagnostics Energy Flux ///////////////////////
!
!       This module handles calculation of terms comprising
!       the RADIAL energy flux.  These are:
!       1.  Enthalpy Flux
!       2.  Kinetic Energy Flux
!       3.  Conductive Heat Flux
!       4.  Poynting Flux
!       5.  Viscous flux of KE
!
!///////////////////////////////////////////////////////////////////
Module Diagnostics_Energy_Flux
    Use Diagnostics_Base
    Use Diagnostics_ADotGradB
    Implicit None

Contains

    Subroutine Compute_Energy_Flux(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8 :: dt_by_dp, dt_by_ds, tpert
        Integer :: r,k, t
        Real*8 :: dr, qadd, fpr2dr

        !First, the radial viscous flux of energy
        If (compute_quantity(visc_flux_r)) Then

            !Radial contribution (mod rho*nu)
            DO_PSI
                qty(PSI) = buffer(PSI,vr)*ref%dlnrho(r)/3.0d0+buffer(PSI,dvrdr)
                qty(PSI) = qty(PSI)*buffer(PSI,vr)*2.0d0
            END_DO

            !Theta contribution (mod rho*nu)
            DO_PSI
                tmp1(PSI) = (2.0d0/3.0d0)*buffer(PSI,vr)*ref%dlnrho(r)
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,dvtdr)-buffer(PSI,vtheta)/radius(r)
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,dvrdt)/radius(r)
                tmp1(PSI) = tmp1(PSI)*buffer(PSI,vtheta)
            END_DO


            DO_PSI
                qty(PSI) = qty(PSI)+tmp1(PSI)
            END_DO

            !phi contribution (mod rho*nu)
            DO_PSI
                tmp1(PSI) = (2.0d0/3.0d0)*buffer(PSI,vr)*ref%dlnrho(r)
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,dvpdr)-buffer(PSI,vphi)/radius(r)
                tmp1(PSI) = tmp1(PSI)+buffer(PSI,dvrdp)/radius(r)/sintheta(t)
                tmp1(PSI) = tmp1(PSI)*buffer(PSI,vphi)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)+tmp1(PSI)
            END_DO

            !Multiply by rho and nu
            DO_PSI
                qty(PSI) = qty(PSI)*nu(r)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !Conductive Flux
        If (compute_quantity(cond_flux_r)) Then
            DO_PSI
                qty(PSI) = -ref%density(r) &
                  & *ref%temperature(r)*kappa(r)*buffer(PSI,dtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !Radial KE Flux
        If (compute_quantity(ke_flux_radial)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vr)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vtheta)**2

            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)*buffer(1:n_phi,:,:,vr)
            DO_PSI
                qty(PSI) = qty(PSI)*ref%density(r)*0.5d0
            END_DO
            Call Add_Quantity(qty)
        Endif

        !Enthalpy Flux
        If (compute_quantity(Enth_flux_radial)) Then

            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    dt_by_ds = ref%temperature(r)/pressure_specific_heat
                    dt_by_dp = 1.0d0/pressure_specific_heat/ref%density(r)
                    Do k = 1, n_phi
                        tpert = dt_by_ds*buffer(PSI,tvar) &
                         & + dt_by_dp*buffer(PSI,pvar)
                        tpert = tpert*ref%density(r) ! This is now T'*rho_bar
                        qty(PSI) = tpert*buffer(PSI,vr)*pressure_specific_heat
                    Enddo
                Enddo
            Enddo
            Call Add_Quantity(qty)
        Endif


        !Thermal Energy Flux (vr {rho_bar T_bar S}  OR vr T)
        If (compute_quantity(thermalE_flux_radial)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vr)*buffer(PSI,tvar)* &
                           & ref%density(r)*ref%temperature(r)
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Volume Heating
        If (compute_quantity(vol_heating)) Then
            If (allocated(ref%heating)) Then
                DO_PSI
                    qty(PSI) = ref%heating(r)* &
                        & ref%density(r)*ref%temperature(r)
                END_DO
            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif

        !The "Flux" associated with the volume heating
        If (compute_quantity(vol_heat_flux)) Then
            tmp1d(:) = 0.0d0
            If (heating_type .gt. 0) Then

                ! Note that radial_integral_weights give int{f r^2}/int(r^2}
                Do r = N_R-1, 1,-1
                    qadd = ref%heating(r)*ref%density(r)*ref%temperature(r) ! the heat
                    qadd = qadd+ ref%heating(r+1)*ref%density(r+1)*ref%temperature(r+1)
                    qadd = qadd*half
                    fpr2dr = (r_squared(r) + r_squared(r+1))*Half*four_pi
                    fpr2dr = fpr2dr*(radius(r)-radius(r+1))
                    tmp1d(r) = tmp1d(r+1)+qadd*fpr2dr
                Enddo
                tmp1d = tmp1d(1)-tmp1d
                tmp1d = tmp1d/four_pi/r_squared
            Endif
            DO_PSI
                qty(PSI) = tmp1d(r)
            END_DO
            Call Add_Quantity(qty)
        Endif



        !/////////// Radial component of ExB
        If (magnetism) Then
            If (compute_quantity(ecrossb_r)) Then
                ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
                DO_PSI
                    tmp1(PSI) = buffer(PSI,vphi)*buffer(PSI,br)- &
                                & buffer(PSI,vr)*buffer(PSI,bphi)
                    tmp1(PSI) = tmp1(PSI)-eta(r)*buffer(PSI,curlbtheta)
                END_DO
                DO_PSI
                    qty(PSI) = tmp1(PSI)*buffer(PSI,bphi) ! E_theta B_phi
                END_DO

                !Next, compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
                DO_PSI
                    tmp1(PSI) = buffer(PSI,vr)*buffer(PSI,btheta)- &
                                & buffer(PSI,vtheta)*buffer(PSI,br)
                    tmp1(PSI) = tmp1(PSI)-eta(r)*buffer(PSI,curlbphi)
                END_DO

                DO_PSI
                    qty(PSI) = qty(PSI)-tmp1(PSI)*buffer(PSI,btheta) ! E_phi B_theta
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

    End Subroutine Compute_Energy_Flux



End Module Diagnostics_Energy_Flux
