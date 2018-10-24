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
Module Diagnostics_Poynting_Flux
    Use Diagnostics_Base
    Implicit None


Contains

    Subroutine Compute_Poynting_Flux(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        !/////////// Radial component of ExB (Full)

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

        !/////////// Theta component of ExB (Full)
        If (compute_quantity(ecrossb_theta)) Then
            ! compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = buffer(PSI,vr)*buffer(PSI,btheta)- &
                            & buffer(PSI,vtheta)*buffer(PSI,br)
                tmp1(PSI) = tmp1(PSI)-eta(r)*buffer(PSI,curlbphi)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*buffer(PSI,br) ! E_phi B_r
            END_DO

            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = buffer(PSI,vtheta)*buffer(PSI,bphi)- &
                            & buffer(PSI,vphi)*buffer(PSI,btheta)
                tmp1(PSI) = tmp1(PSI)-eta(r)*buffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*buffer(PSI,bphi) ! E_r B_phi
            END_DO
            Call Add_Quantity(qty)
        Endif

        !/////////// Phi component of ExB (Full)
        If (compute_quantity(ecrossb_phi)) Then
            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = buffer(PSI,vtheta)*buffer(PSI,bphi)- &
                            & buffer(PSI,vphi)*buffer(PSI,btheta)
                tmp1(PSI) = tmp1(PSI)-eta(r)*buffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*buffer(PSI,btheta) ! E_r B_theta
            END_DO

            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = buffer(PSI,vphi)*buffer(PSI,br)- &
                            & buffer(PSI,vr)*buffer(PSI,bphi)
                tmp1(PSI) = tmp1(PSI)-eta(r)*buffer(PSI,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*buffer(PSI,br) ! E_theta B_r
            END_DO


            Call Add_Quantity(qty)
        Endif

        !///////////////////////////////////////////////////////////////////////
        ! Next, we have terms of the formed from triple products of fluctuations
        ! Radial

        If (compute_quantity(ecrossb_ppp_r)) Then
            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vphi)*fbuffer(PSI,br)- &
                            & fbuffer(PSI,vr)*fbuffer(PSI,bphi)
                tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*fbuffer(PSI,bphi) ! E_theta B_phi
            END_DO

            !Next, compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr)*fbuffer(PSI,btheta)- &
                            & fbuffer(PSI,vtheta)*fbuffer(PSI,br)
                tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbphi)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*fbuffer(PSI,btheta) ! E_phi B_theta
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta
        If (compute_quantity(ecrossb_ppp_theta)) Then
            ! compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr)*fbuffer(PSI,btheta)- &
                            & fbuffer(PSI,vtheta)*fbuffer(PSI,br)
                tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbphi)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*fbuffer(PSI,br) ! E_phi B_r
            END_DO

            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vtheta)*fbuffer(PSI,bphi)- &
                            & fbuffer(PSI,vphi)*fbuffer(PSI,btheta)
                tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*fbuffer(PSI,bphi) ! E_r B_phi
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi
        If (compute_quantity(ecrossb_ppp_phi)) Then
            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vtheta)*fbuffer(PSI,bphi)- &
                            & fbuffer(PSI,vphi)*fbuffer(PSI,btheta)
                tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*fbuffer(PSI,btheta) ! E_r B_theta
            END_DO

            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vphi)*fbuffer(PSI,br)- &
                            & fbuffer(PSI,vr)*fbuffer(PSI,bphi)
                tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*fbuffer(PSI,br) ! E_theta B_r
            END_DO


            Call Add_Quantity(qty)
        Endif

        !///////////////////////////////////////////////////////////////////////
        ! Next, we have terms of the formed from triple products of means
        ! Radial

        If (compute_quantity(ecrossb_mmm_r)) Then
            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vphi)*m0_values(PSI2,br)- &
                            & m0_values(PSI2,vr)*m0_values(PSI2,bphi)
                tmp1(PSI) = tmp1(PSI)-eta(r)*m0_values(PSI2,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*m0_values(PSI2,bphi) ! E_theta B_phi
            END_DO

            !Next, compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vr)*m0_values(PSI2,btheta)- &
                            & m0_values(PSI2,vtheta)*m0_values(PSI2,br)
                tmp1(PSI) = tmp1(PSI)-eta(r)*m0_values(PSI2,curlbphi)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*m0_values(PSI2,btheta) ! E_phi B_theta
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta
        If (compute_quantity(ecrossb_mmm_theta)) Then
            ! compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vr)*m0_values(PSI2,btheta)- &
                            & m0_values(PSI2,vtheta)*m0_values(PSI2,br)
                tmp1(PSI) = tmp1(PSI)-eta(r)*m0_values(PSI2,curlbphi)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*m0_values(PSI2,br) ! E_phi B_r
            END_DO

            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vtheta)*m0_values(PSI2,bphi)- &
                            & m0_values(PSI2,vphi)*m0_values(PSI2,btheta)
                tmp1(PSI) = tmp1(PSI)-eta(r)*m0_values(PSI2,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*m0_values(PSI2,bphi) ! E_r B_phi
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi
        If (compute_quantity(ecrossb_mmm_phi)) Then
            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vtheta)*m0_values(PSI2,bphi)- &
                            & m0_values(PSI2,vphi)*m0_values(PSI2,btheta)
                tmp1(PSI) = tmp1(PSI)-eta(r)*m0_values(PSI2,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*m0_values(PSI2,btheta) ! E_r B_theta
            END_DO

            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vphi)*m0_values(PSI2,br)- &
                            & m0_values(PSI2,vr)*m0_values(PSI2,bphi)
                tmp1(PSI) = tmp1(PSI)-eta(r)*m0_values(PSI2,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*m0_values(PSI2,br) ! E_theta B_r
            END_DO


            Call Add_Quantity(qty)
        Endif

        !///////////////////////////////////////////////////////////////////////
        ! Next, we have terms of the form (v' x B') x < B >
        ! Radial

        If (compute_quantity(ecrossb_ppm_r)) Then
            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vphi)*fbuffer(PSI,br)- &
                            & fbuffer(PSI,vr)*fbuffer(PSI,bphi)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*m0_values(PSI2,bphi) ! E_theta B_phi
            END_DO

            !Next, compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr)*fbuffer(PSI,btheta)- &
                            & fbuffer(PSI,vtheta)*fbuffer(PSI,br)
               ! tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbphi)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*m0_values(PSI2,btheta) ! E_phi B_theta
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta
        If (compute_quantity(ecrossb_ppm_theta)) Then
            ! compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr)*fbuffer(PSI,btheta)- &
                            & fbuffer(PSI,vtheta)*fbuffer(PSI,br)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbphi)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*m0_values(PSI2,br) ! E_phi B_r
            END_DO

            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vtheta)*fbuffer(PSI,bphi)- &
                            & fbuffer(PSI,vphi)*fbuffer(PSI,btheta)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*m0_values(PSI2,bphi) ! E_r B_phi
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi
        If (compute_quantity(ecrossb_ppm_phi)) Then
            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vtheta)*fbuffer(PSI,bphi)- &
                            & fbuffer(PSI,vphi)*fbuffer(PSI,btheta)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*m0_values(PSI2,btheta) ! E_r B_theta
            END_DO

            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vphi)*fbuffer(PSI,br)- &
                            & fbuffer(PSI,vr)*fbuffer(PSI,bphi)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*m0_values(PSI2,br) ! E_theta B_r
            END_DO


            Call Add_Quantity(qty)
        Endif

        !///////////////////////////////////////////////////////////////////////
        ! Next, we have terms of the form (v' x <B> ) x B'
        ! Radial

        If (compute_quantity(ecrossb_pmp_r)) Then
            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vphi)*m0_values(PSI2,br)- &
                            & fbuffer(PSI,vr)*m0_values(PSI2,bphi)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*m0_values(PSI2,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*fbuffer(PSI,bphi) ! E_theta B_phi
            END_DO

            !Next, compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr)*m0_values(PSI2,btheta)- &
                            & fbuffer(PSI,vtheta)*m0_values(PSI2,br)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*m0_values(PSI2,curlbphi)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*fbuffer(PSI,btheta) ! E_phi B_theta
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta
        If (compute_quantity(ecrossb_pmp_theta)) Then
            ! compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vr)*m0_values(PSI2,btheta)- &
                            & fbuffer(PSI,vtheta)*m0_values(PSI2,br)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbphi)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*fbuffer(PSI,br) ! E_phi B_r
            END_DO

            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vtheta)*m0_values(PSI2,bphi)- &
                            & fbuffer(PSI,vphi)*m0_values(PSI2,btheta)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*fbuffer(PSI,bphi) ! E_r B_phi
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi
        If (compute_quantity(ecrossb_pmp_phi)) Then
            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vtheta)*m0_values(PSI2,bphi)- &
                            & fbuffer(PSI,vphi)*m0_values(PSI2,btheta)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*fbuffer(PSI,btheta) ! E_r B_theta
            END_DO

            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = fbuffer(PSI,vphi)*m0_values(PSI2,br)- &
                            & fbuffer(PSI,vr)*m0_values(PSI2,bphi)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*fbuffer(PSI,br) ! E_theta B_r
            END_DO


            Call Add_Quantity(qty)
        Endif

        !///////////////////////////////////////////////////////////////////////
        ! Finally, we have terms of the form (<v> X B') X B'
        ! Radial

        If (compute_quantity(ecrossb_mpp_r)) Then
            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vphi)*fbuffer(PSI,br)- &
                            & m0_values(PSI2,vr)*fbuffer(PSI,bphi)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*fbuffer(PSI,bphi) ! E_theta B_phi
            END_DO

            !Next, compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vr)*fbuffer(PSI,btheta)- &
                            & m0_values(PSI2,vtheta)*fbuffer(PSI,br)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbphi)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*fbuffer(PSI,btheta) ! E_phi B_theta
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta
        If (compute_quantity(ecrossb_mpp_theta)) Then
            ! compute [vxB]_phi -eta (j_phi) {i.e., E_phi}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vr)*fbuffer(PSI,btheta)- &
                            & m0_values(PSI2,vtheta)*fbuffer(PSI,br)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbphi)
            END_DO
            DO_PSI
                qty(PSI) = tmp1(PSI)*fbuffer(PSI,br) ! E_phi B_r
            END_DO

            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vtheta)*fbuffer(PSI,bphi)- &
                            & m0_values(PSI2,vphi)*fbuffer(PSI,btheta)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*fbuffer(PSI,bphi) ! E_r B_phi
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi
        If (compute_quantity(ecrossb_mpp_phi)) Then
            !Next, compute [vxB]_r -eta (j_r) {i.e., E_r}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vtheta)*fbuffer(PSI,bphi)- &
                            & m0_values(PSI2,vphi)*fbuffer(PSI,btheta)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbr)
            END_DO

            DO_PSI
                qty(PSI) = tmp1(PSI)*fbuffer(PSI,btheta) ! E_r B_theta
            END_DO

            ! compute [vxB]_theta -eta (j_theta) {i.e., E_theta}
            DO_PSI
                tmp1(PSI) = m0_values(PSI2,vphi)*fbuffer(PSI,br)- &
                            & m0_values(PSI2,vr)*fbuffer(PSI,bphi)
                !tmp1(PSI) = tmp1(PSI)-eta(r)*fbuffer(PSI,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = qty(PSI)-tmp1(PSI)*fbuffer(PSI,br) ! E_theta B_r
            END_DO


            Call Add_Quantity(qty)
        Endif
    End Subroutine Compute_Poynting_Flux


End Module Diagnostics_Poynting_Flux
