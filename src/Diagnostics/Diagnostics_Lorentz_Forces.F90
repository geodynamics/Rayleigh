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

Module Diagnostics_Lorentz_Forces

    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_Lorentz_Forces(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        ! Full JxB terms
        If (compute_quantity(j_cross_b_r) .or. compute_quantity(mag_work)) Then
            DO_PSI
                !qty(PSI) = (buffer(PSI,curlbtheta)*buffer(PSI,bphi)- &
                !         & buffer(PSI,btheta)*buffer(PSI,curlbphi) ) *ref%Lorentz_Coeff
                qty(PSI) = mean_3dbuffer(PSI,lforce_r) - mean_ell0buffer(r,lforce_r)
            END_DO
            If (compute_quantity(j_cross_b_r)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work)) Then
                DO_PSI
                    tmp1(PSI) = qty(PSI)*buffer(PSI,vr)
                END_DO
            Endif
        Endif

        If (compute_quantity(j_cross_b_theta) .or. compute_quantity(mag_work)) Then
            DO_PSI
                qty(PSI) = ( buffer(PSI,br)*buffer(PSI,curlbphi)- &
                         & buffer(PSI,curlbr)*buffer(PSI,bphi) )*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(j_cross_b_theta)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*buffer(PSI,vtheta)
                END_DO
            Endif
        Endif

        If (compute_quantity(j_cross_b_phi) .or. compute_quantity(mag_work)) Then
            DO_PSI
                qty(PSI) = ( buffer(PSI,curlbr)*buffer(PSI,btheta)- &
                           & buffer(PSI,br)*buffer(PSI,curlbtheta) )*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(j_cross_b_phi)) Call Add_Quantity(qty)

            If (compute_quantity(mag_work)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*buffer(PSI,vphi)
                END_DO
                Call Add_Quantity(tmp1)
            Endif
        Endif

        !///////////////////////////////////////////////
        !                   <J> x <B> terms

        If (compute_quantity(jm_cross_bm_r) .or. compute_quantity(mag_work_mmm)) Then
            DO_PSI2
                !qty(1:n_phi,PSI2) = ( m0_values(PSI2,curlbtheta)*m0_values(PSI2,bphi)- &
                !                  & m0_values(PSI2,btheta)*m0_values(PSI2,curlbphi) )*ref%Lorentz_Coeff
                qty(1:n_phi, PSI2) = mean_3dbuffer(1:n_phi,PSI2,lforcemm_r) - mean_ell0buffer(r, lforcemm_r)
            END_DO2
            If (compute_quantity(jm_cross_bm_r)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work_mmm)) Then
                DO_PSI
                    tmp1(PSI) = qty(PSI)*m0_values(PSI2,vr)
                END_DO
            Endif
        Endif

        If (compute_quantity(jm_cross_bm_theta) .or. compute_quantity(mag_work_mmm)) Then
            DO_PSI2
                qty(1:n_phi,PSI2) = ( m0_values(PSI2,br)*m0_values(PSI2,curlbphi)- &
                                  & m0_values(PSI2,curlbr)*m0_values(PSI2,bphi) )*ref%Lorentz_Coeff
            END_DO2
            If (compute_quantity(jm_cross_bm_theta)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work_mmm)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*m0_values(PSI2,vtheta)
                END_DO
            Endif
        Endif

        If (compute_quantity(jm_cross_bm_phi) .or. compute_quantity(samom_lorentz_mm) &
            .or. compute_quantity(mag_work_mmm)) Then
            DO_PSI2
                qty(1:n_phi,PSI2) = ( m0_values(PSI2,curlbr)*m0_values(PSI2,btheta)- &
                                  & m0_values(PSI2,br)*m0_values(PSI2,curlbtheta) )*ref%Lorentz_Coeff
            END_DO2
            If (compute_quantity(jm_cross_bm_phi)) Call Add_Quantity(qty)

            If (compute_quantity(mag_work_mmm)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*m0_values(PSI2,vphi)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(samom_lorentz_mm)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*radius(r)*sintheta(t)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif


        !////////////////////////////////////////////////////////
        !  J' X B' terms
        If (compute_quantity(jp_cross_bp_r) .or. compute_quantity(mag_work_ppp) &
            .or. compute_quantity(mag_work_mpp)) Then
            DO_PSI
                !qty(PSI) = ( fbuffer(PSI,curlbtheta)*fbuffer(PSI,bphi)- &
                !         & fbuffer(PSI,btheta)*fbuffer(PSI,curlbphi) )*ref%Lorentz_Coeff
                qty(PSI) = mean_3dbuffer(PSI,lforcepp_r) - mean_ell0buffer(r,lforcepp_r)
            END_DO
            If (compute_quantity(jp_cross_bp_r)) Call Add_Quantity(qty)

            If (compute_quantity(mag_work_mpp)) Then
                DO_PSI
                    tmp1(PSI) = qty(PSI)*m0_values(PSI2,vr)
                END_DO
            Endif

            If (compute_quantity(mag_work_ppp)) Then
                DO_PSI
                    tmp4(PSI) = qty(PSI)*fbuffer(PSI,vr)
                END_DO
            Endif

        Endif

        If (compute_quantity(jp_cross_bp_theta) .or. compute_quantity(mag_work_ppp) &
            .or. compute_quantity(mag_work_mpp)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,br)*fbuffer(PSI,curlbphi)- &
                         & fbuffer(PSI,curlbr)*fbuffer(PSI,bphi) )*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(jp_cross_bp_theta)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work_mpp)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*m0_values(PSI2,vtheta)
                END_DO
            Endif

            If (compute_quantity(mag_work_ppp)) Then
                DO_PSI
                    tmp4(PSI) = tmp4(PSI)+qty(PSI)*fbuffer(PSI,vtheta)
                END_DO
            Endif
        Endif

        If (compute_quantity(jp_cross_bp_phi) .or. compute_quantity(samom_lorentz_pp) &
            .or. compute_quantity(mag_work_ppp) .or. compute_quantity(mag_work_mmm)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,curlbr)*fbuffer(PSI,btheta)- &
                           & fbuffer(PSI,br)*fbuffer(PSI,curlbtheta) )*ref%Lorentz_Coeff
            END_DO

            If (compute_quantity(jp_cross_bp_phi))            Call Add_Quantity(qty)
            If (compute_quantity(mag_work_mpp)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*m0_values(PSI2,vphi)
                END_DO
                Call Add_Quantity(tmp1)
            Endif

            If (compute_quantity(mag_work_ppp)) Then
                DO_PSI
                    tmp4(PSI) = tmp4(PSI)+qty(PSI)*fbuffer(PSI,vphi)
                END_DO
                Call Add_Quantity(tmp4)
            Endif

            If (compute_quantity(samom_lorentz_pp)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*radius(r)*sintheta(t)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        !//////////////////////////////////////////////////////////
        !               J' x <B> terms
        If (compute_quantity(jp_cross_bm_r) .or. compute_quantity(mag_work_ppm)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,curlbtheta)*m0_values(PSI2,bphi)- &
                         & m0_values(PSI2,btheta)*fbuffer(PSI,curlbphi) )*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(jp_cross_bm_r)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work_ppm)) Then
                DO_PSI
                    tmp1(PSI) = qty(PSI)*fbuffer(PSI,vr)
                END_DO
            Endif
        Endif

        If (compute_quantity(jp_cross_bm_theta) .or. compute_quantity(mag_work_ppm)) Then
            DO_PSI
                qty(PSI) = ( m0_values(PSI2,br)*fbuffer(PSI,curlbphi)- &
                         & fbuffer(PSI,curlbr)*m0_values(PSI2,bphi) )*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(jp_cross_bm_theta)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work_ppm)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*fbuffer(PSI,vtheta)
                END_DO
            Endif
        Endif

        If (compute_quantity(jp_cross_bm_phi) .or. compute_quantity(mag_work_ppm)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,curlbr)*m0_values(PSI2,btheta)- &
                           & m0_values(PSI2,br)*fbuffer(PSI,curlbtheta) )*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(jp_cross_bm_phi)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work_ppm)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*fbuffer(PSI,vphi)
                END_DO
                Call Add_Quantity(tmp1)
            Endif
        Endif

        !//////////////////////////////////////////////////////////////
        !               <J> x B' terms
        !////////////////////////////////////////////////////////
        !  J' X B' terms
        If (compute_quantity(jm_cross_bp_r) .or. compute_quantity(mag_work_pmp)) Then
            DO_PSI
                qty(PSI) = ( m0_values(PSI2,curlbtheta)*fbuffer(PSI,bphi)- &
                         & fbuffer(PSI,btheta)*m0_values(PSI2,curlbphi) )*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(jm_cross_bp_r)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work_pmp)) Then
                DO_PSI
                    tmp1(PSI) = qty(PSI)*fbuffer(PSI,vr)
                END_DO
            Endif
        Endif

        If (compute_quantity(jm_cross_bp_theta) .or. compute_quantity(mag_work_pmp)) Then
            DO_PSI
                qty(PSI) = ( fbuffer(PSI,br)*m0_values(PSI2,curlbphi)- &
                         & m0_values(PSI2,curlbr)*fbuffer(PSI,bphi) )*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(jm_cross_bp_theta)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work_pmp)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*fbuffer(PSI,vtheta)
                END_DO
            Endif
        Endif

        If (compute_quantity(jm_cross_bp_phi) .or. compute_quantity(mag_work_pmp)) Then
            DO_PSI
                qty(PSI) = ( m0_values(PSI2,curlbr)*fbuffer(PSI,btheta)- &
                           & fbuffer(PSI,br)*m0_values(PSI2,curlbtheta) )*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(jm_cross_bp_phi)) Call Add_Quantity(qty)
            If (compute_quantity(mag_work_pmp)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*fbuffer(PSI,vphi)
                END_DO
                Call Add_Quantity(tmp1)
            Endif
        Endif

    End Subroutine Compute_Lorentz_Forces


End Module Diagnostics_Lorentz_Forces
