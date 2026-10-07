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

Module Diagnostics_Curl_Momentum
    Use Diagnostics_Base
    Use Spectral_Derivatives
    Use Finite_Difference, Only : d_by_dx3d3
    Use Structures
    Use Load_Balance, Only : l_lm_values, my_lm_min
    Implicit None

    ! Internal slot indices of Viscous_Force_Derivatives, per force set (1=full,
    ! 2=fluctuating, 3=mean); -1 when not needed.  The output indices (vfd_*)
    ! are in Diagnostics_Base.
    Integer :: nvftrans = 0, nvfwork = 0, nvfback = 0
    Integer :: vf_fr(3), vf_sft(3), vf_sfp(3), vf_q(3)
    Integer :: vf_sft_dr(3), vf_sfp_dr(3), vf_q_dr(3), vf_q_d2r(3)

Contains

    Subroutine Compute_Curl_Momentum_Forces(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Call Compute_Curl_Advection_Force(buffer)
        Call Compute_Curl_Buoyancy_Force(buffer)
        Call Compute_Curl_Magnetic_Force(buffer)
        Call Compute_Curl_Coriolis_Force(buffer)
        Call Compute_Curl_Pressure_Force(buffer)
        Call Compute_Curl_Viscous_Force(buffer)
    End Subroutine

    Subroutine Compute_Curl_Advection_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t
        Real*8 :: vgv_abs_r, vgv_abs_t, vgv_abs_p

        Real*8  :: pfactor(my_r%min:my_r%max)
        pfactor(my_r%min:my_r%max) = ref%dpdr_w_term(my_r%min:my_r%max) &
                                        /ref%density(my_r%min:my_r%max)
        
        !!!!!!!!!!!!!!!!!!!!!!!!! Advection Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        If (compute_quantity(curl_v_grad_v_r) .or. compute_quantity(curl_v_grad_v_r_squared)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvpdrdt) * buffer(PSI,vr) * one_over_r(r) * ref%density(r) &
                + DDBUFF(PSI,dvpdtdt) * buffer(PSI,vtheta) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvpdr) * buffer(PSI,dvrdt) * one_over_r(r) * ref%density(r) &
                + buffer(PSI,dvpdt) * buffer(PSI,dvtdt) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvpdt) * buffer(PSI,vr) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvrdt) * buffer(PSI,vphi) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,vphi) * buffer(PSI,vtheta) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + DDBUFF(PSI,dvpdtdp) * buffer(PSI,vphi) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvpdp) * buffer(PSI,dvpdt) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - DDBUFF(PSI,dvtdpdp) * buffer(PSI,vphi) * ref%density(r) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                - DDBUFF(PSI,dvtdrdp) * buffer(PSI,vr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                - DDBUFF(PSI,dvtdtdp) * buffer(PSI,vtheta) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r) &
                - buffer(PSI,dvpdp) * buffer(PSI,dvtdp) * ref%density(r) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                - buffer(PSI,dvrdp) * buffer(PSI,dvtdr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                - buffer(PSI,dvrdp) * buffer(PSI,vtheta) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvtdp) * buffer(PSI,dvtdt) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvtdp) * buffer(PSI,vr) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvpdr) * buffer(PSI,vr) * costheta(t) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + buffer(PSI,dvtdt) * buffer(PSI,vphi) * costheta(t) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r) &
                + buffer(PSI,vphi) * buffer(PSI,vr) * costheta(t) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r) &
                + 2 * buffer(PSI,dvpdp) * buffer(PSI,vphi) * costheta(t) * ref%density(r) * csctheta(t) * csctheta(t) &
                    * one_over_r(r) * one_over_r(r) &
                + 2 * buffer(PSI,dvpdt) * buffer(PSI,vtheta) * costheta(t) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r)
            END_DO
            If (compute_quantity(curl_v_grad_v_r)) Call Add_Quantity(qty)
            If (compute_quantity(curl_v_grad_v_r_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif            
            
        Endif
    
        If (compute_quantity(curl_v_grad_v_theta) .or. compute_quantity(curl_v_grad_v_theta_squared)) Then
            DO_PSI
                qty(PSI) = -DDBUFF(PSI,dvpdrdr) * buffer(PSI,vr) * ref%density(r) &
                - buffer(PSI,dvpdr) * buffer(PSI,dvrdr) * ref%density(r) &
                - DDBUFF(PSI,dvpdrdt) * buffer(PSI,vtheta) * one_over_r(r) * ref%density(r) &
                - buffer(PSI,dvpdr) * buffer(PSI,vr) * ref%density(r) * ref%dlnrho(r) &
                - buffer(PSI,dvpdt) * buffer(PSI,dvtdr) * one_over_r(r) * ref%density(r) &
                - buffer(PSI,dvrdr) * buffer(PSI,vphi) * one_over_r(r) * ref%density(r) &
                - 2 * buffer(PSI,dvpdr) * buffer(PSI,vr) * one_over_r(r) * ref%density(r) &
                + DDBUFF(PSI,dvrdpdp) * buffer(PSI,vphi) * ref%density(r) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                + DDBUFF(PSI,dvrdrdp) * buffer(PSI,vr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + DDBUFF(PSI,dvrdtdp) * buffer(PSI,vtheta) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r) &
                + buffer(PSI,dvpdp) * buffer(PSI,dvrdp) * ref%density(r) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                + buffer(PSI,dvrdp) * buffer(PSI,dvrdr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + buffer(PSI,dvrdt) * buffer(PSI,dvtdp) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - DDBUFF(PSI,dvpdrdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%density(r) &
                - buffer(PSI,dvpdp) * buffer(PSI,dvpdr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                - buffer(PSI,dvpdr) * buffer(PSI,vtheta) * cottheta(t) * one_over_r(r) * ref%density(r) &
                - buffer(PSI,dvpdt) * buffer(PSI,vtheta) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                - buffer(PSI,dvtdr) * buffer(PSI,vphi) * cottheta(t) * one_over_r(r) * ref%density(r) &
                - buffer(PSI,vphi) * buffer(PSI,vr) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                - 2 * buffer(PSI,dvpdp) * buffer(PSI,vphi) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - 2 * buffer(PSI,dvtdp) * buffer(PSI,vtheta) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvpdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                - buffer(PSI,vphi) * buffer(PSI,vtheta) * cottheta(t) * one_over_r(r) * ref%density(r) * ref%dlnrho(r)
            END_DO
            If (compute_quantity(curl_v_grad_v_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_v_grad_v_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif
        
        If (compute_quantity(curl_v_grad_v_phi) .or. compute_quantity(curl_v_grad_v_phi_squared)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvtdrdr) * buffer(PSI,vr) * ref%density(r) &
                + buffer(PSI,dvrdr) * buffer(PSI,dvtdr) * ref%density(r) &
                + DDBUFF(PSI,dvtdrdt) * buffer(PSI,vtheta) * one_over_r(r) * ref%density(r) &
                + buffer(PSI,dvrdr) * buffer(PSI,vtheta) * one_over_r(r) * ref%density(r) &
                + buffer(PSI,dvtdr) * buffer(PSI,dvtdt) * one_over_r(r) * ref%density(r) &
                + buffer(PSI,dvtdr) * buffer(PSI,vr) * ref%density(r) * ref%dlnrho(r) &
                - DDBUFF(PSI,dvrdrdt) * buffer(PSI,vr) * one_over_r(r) * ref%density(r) &
                - DDBUFF(PSI,dvrdtdt) * buffer(PSI,vtheta) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvrdr) * buffer(PSI,dvrdt) * one_over_r(r) * ref%density(r) &
                - buffer(PSI,dvrdt) * buffer(PSI,dvtdt) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + 2 * buffer(PSI,dvpdt) * buffer(PSI,vphi) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + 2 * buffer(PSI,dvtdr) * buffer(PSI,vr) * one_over_r(r) * ref%density(r) &
                + 2 * buffer(PSI,dvtdt) * buffer(PSI,vtheta) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + DDBUFF(PSI,dvtdrdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + buffer(PSI,dvpdr) * buffer(PSI,dvtdp) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + buffer(PSI,dvtdt) * buffer(PSI,vtheta) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                + buffer(PSI,vr) * buffer(PSI,vtheta) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                - DDBUFF(PSI,dvrdtdp) * buffer(PSI,vphi) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvpdt) * buffer(PSI,dvrdp) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - cottheta(t) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) * buffer(PSI,vphi) * buffer(PSI,vphi) &
                - 2 * buffer(PSI,dvpdr) * buffer(PSI,vphi) * cottheta(t) * one_over_r(r) * ref%density(r) &
                + buffer(PSI,dvrdp) * buffer(PSI,vphi) * costheta(t) * ref%density(r) * csctheta(t) * csctheta(t) &
                    * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvtdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%density(r) * ref%dlnrho(r)
            END_DO
            If (compute_quantity(curl_v_grad_v_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_v_grad_v_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        If (compute_quantity(curl_v_grad_v_abs)) Then
            DO_PSI
                vgv_abs_r = DDBUFF(PSI,dvpdrdt) * buffer(PSI,vr) * one_over_r(r) &
                + DDBUFF(PSI,dvpdtdt) * buffer(PSI,vtheta) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvpdr) * buffer(PSI,dvrdt) * one_over_r(r) &
                + buffer(PSI,dvpdt) * buffer(PSI,dvtdt) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvpdt) * buffer(PSI,vr) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvrdt) * buffer(PSI,vphi) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,vphi) * buffer(PSI,vtheta) * one_over_r(r) * one_over_r(r) &
                + DDBUFF(PSI,dvpdtdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvpdp) * buffer(PSI,dvpdt) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - DDBUFF(PSI,dvtdpdp) * buffer(PSI,vphi) * csctheta(t) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - DDBUFF(PSI,dvtdrdp) * buffer(PSI,vr) * csctheta(t) * one_over_r(r) &
                - DDBUFF(PSI,dvtdtdp) * buffer(PSI,vtheta) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvpdp) * buffer(PSI,dvtdp) * csctheta(t) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvrdp) * buffer(PSI,dvtdr) * csctheta(t) * one_over_r(r) &
                - buffer(PSI,dvrdp) * buffer(PSI,vtheta) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvtdp) * buffer(PSI,dvtdt) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvtdp) * buffer(PSI,vr) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvpdr) * buffer(PSI,vr) * costheta(t) * csctheta(t) * one_over_r(r) &
                + buffer(PSI,dvtdt) * buffer(PSI,vphi) * costheta(t) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,vphi) * buffer(PSI,vr) * costheta(t) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                + 2 * buffer(PSI,dvpdp) * buffer(PSI,vphi) * costheta(t) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                + 2 * buffer(PSI,dvpdt) * buffer(PSI,vtheta) * costheta(t) * csctheta(t) * one_over_r(r) * one_over_r(r)

                vgv_abs_t = -DDBUFF(PSI,dvpdrdr) * buffer(PSI,vr) &
                - buffer(PSI,dvpdr) * buffer(PSI,dvrdr) &
                - DDBUFF(PSI,dvpdrdt) * buffer(PSI,vtheta) * one_over_r(r) &
                - buffer(PSI,dvpdr) * buffer(PSI,vr) * ref%dlnrho(r) &
                - buffer(PSI,dvpdt) * buffer(PSI,dvtdr) * one_over_r(r) &
                - buffer(PSI,dvrdr) * buffer(PSI,vphi) * one_over_r(r) &
                - 2 * buffer(PSI,dvpdr) * buffer(PSI,vr) * one_over_r(r) &
                + DDBUFF(PSI,dvrdpdp) * buffer(PSI,vphi) * csctheta(t) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                + DDBUFF(PSI,dvrdrdp) * buffer(PSI,vr) * csctheta(t) * one_over_r(r) &
                + DDBUFF(PSI,dvrdtdp) * buffer(PSI,vtheta) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvpdp) * buffer(PSI,dvrdp) * csctheta(t) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                + buffer(PSI,dvrdp) * buffer(PSI,dvrdr) * csctheta(t) * one_over_r(r) &
                + buffer(PSI,dvrdt) * buffer(PSI,dvtdp) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - DDBUFF(PSI,dvpdrdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) &
                - buffer(PSI,dvpdp) * buffer(PSI,dvpdr) * csctheta(t) * one_over_r(r) &
                - buffer(PSI,dvpdr) * buffer(PSI,vtheta) * cottheta(t) * one_over_r(r) &
                - buffer(PSI,dvpdt) * buffer(PSI,vtheta) * one_over_r(r) * ref%dlnrho(r) &
                - buffer(PSI,dvtdr) * buffer(PSI,vphi) * cottheta(t) * one_over_r(r) &
                - buffer(PSI,vphi) * buffer(PSI,vr) * one_over_r(r) * ref%dlnrho(r) &
                - 2 * buffer(PSI,dvpdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - 2 * buffer(PSI,dvtdp) * buffer(PSI,vtheta) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvpdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%dlnrho(r) &
                - buffer(PSI,vphi) * buffer(PSI,vtheta) * cottheta(t) * one_over_r(r) * ref%dlnrho(r)

                vgv_abs_p = DDBUFF(PSI,dvtdrdr) * buffer(PSI,vr) &
                + buffer(PSI,dvrdr) * buffer(PSI,dvtdr) &
                + DDBUFF(PSI,dvtdrdt) * buffer(PSI,vtheta) * one_over_r(r) &
                + buffer(PSI,dvrdr) * buffer(PSI,vtheta) * one_over_r(r) &
                + buffer(PSI,dvtdr) * buffer(PSI,dvtdt) * one_over_r(r) &
                + buffer(PSI,dvtdr) * buffer(PSI,vr) * ref%dlnrho(r) &
                - DDBUFF(PSI,dvrdrdt) * buffer(PSI,vr) * one_over_r(r) &
                - DDBUFF(PSI,dvrdtdt) * buffer(PSI,vtheta) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvrdr) * buffer(PSI,dvrdt) * one_over_r(r) &
                - buffer(PSI,dvrdt) * buffer(PSI,dvtdt) * one_over_r(r) * one_over_r(r) &
                + 2 * buffer(PSI,dvpdt) * buffer(PSI,vphi) * one_over_r(r) * one_over_r(r) &
                + 2 * buffer(PSI,dvtdr) * buffer(PSI,vr) * one_over_r(r) &
                + 2 * buffer(PSI,dvtdt) * buffer(PSI,vtheta) * one_over_r(r) * one_over_r(r) &
                + DDBUFF(PSI,dvtdrdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) &
                + buffer(PSI,dvpdr) * buffer(PSI,dvtdp) * csctheta(t) * one_over_r(r) &
                + buffer(PSI,dvtdt) * buffer(PSI,vtheta) * one_over_r(r) * ref%dlnrho(r) &
                + buffer(PSI,vr) * buffer(PSI,vtheta) * one_over_r(r) * ref%dlnrho(r) &
                - DDBUFF(PSI,dvrdtdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - buffer(PSI,dvpdt) * buffer(PSI,dvrdp) * csctheta(t) * one_over_r(r) * one_over_r(r) &
                - cottheta(t) * one_over_r(r) * ref%dlnrho(r) * buffer(PSI,vphi) * buffer(PSI,vphi) &
                - 2 * buffer(PSI,dvpdr) * buffer(PSI,vphi) * cottheta(t) * one_over_r(r) &
                + buffer(PSI,dvrdp) * buffer(PSI,vphi) * costheta(t) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                + buffer(PSI,dvtdp) * buffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%dlnrho(r)

                qty(PSI) = ref%density(r) * sqrt(vgv_abs_r * vgv_abs_r + vgv_abs_t * vgv_abs_t + vgv_abs_p * vgv_abs_p)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vp_grad_vp_r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvpdrdt) * fbuffer(PSI,vr) * one_over_r(r) * ref%density(r) &
                + d2_fbuffer(PSI,dvpdtdt) * fbuffer(PSI,vtheta) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + fbuffer(PSI,dvpdr) * fbuffer(PSI,dvrdt) * one_over_r(r) * ref%density(r) &
                + fbuffer(PSI,dvpdt) * fbuffer(PSI,dvtdt) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + fbuffer(PSI,dvpdt) * fbuffer(PSI,vr) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + fbuffer(PSI,dvrdt) * fbuffer(PSI,vphi) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - fbuffer(PSI,vphi) * fbuffer(PSI,vtheta) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + d2_fbuffer(PSI,dvpdtdp) * fbuffer(PSI,vphi) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + fbuffer(PSI,dvpdp) * fbuffer(PSI,dvpdt) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - d2_fbuffer(PSI,dvtdpdp) * fbuffer(PSI,vphi) * ref%density(r) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                - d2_fbuffer(PSI,dvtdrdp) * fbuffer(PSI,vr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                - d2_fbuffer(PSI,dvtdtdp) * fbuffer(PSI,vtheta) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r) &
                - fbuffer(PSI,dvpdp) * fbuffer(PSI,dvtdp) * ref%density(r) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                - fbuffer(PSI,dvrdp) * fbuffer(PSI,dvtdr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                - fbuffer(PSI,dvrdp) * fbuffer(PSI,vtheta) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - fbuffer(PSI,dvtdp) * fbuffer(PSI,dvtdt) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - fbuffer(PSI,dvtdp) * fbuffer(PSI,vr) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + fbuffer(PSI,dvpdr) * fbuffer(PSI,vr) * costheta(t) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + fbuffer(PSI,dvtdt) * fbuffer(PSI,vphi) * costheta(t) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r) &
                + fbuffer(PSI,vphi) * fbuffer(PSI,vr) * costheta(t) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r) &
                + 2 * fbuffer(PSI,dvpdp) * fbuffer(PSI,vphi) * costheta(t) * ref%density(r) * csctheta(t) * csctheta(t) &
                    * one_over_r(r) * one_over_r(r) &
                + 2 * fbuffer(PSI,dvpdt) * fbuffer(PSI,vtheta) * costheta(t) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vp_grad_vp_theta)) Then
            DO_PSI
                qty(PSI) = -d2_fbuffer(PSI,dvpdrdr) * fbuffer(PSI,vr) * ref%density(r) &
                - fbuffer(PSI,dvpdr) * fbuffer(PSI,dvrdr) * ref%density(r) &
                - d2_fbuffer(PSI,dvpdrdt) * fbuffer(PSI,vtheta) * one_over_r(r) * ref%density(r) &
                - fbuffer(PSI,dvpdr) * fbuffer(PSI,vr) * ref%density(r) * ref%dlnrho(r) &
                - fbuffer(PSI,dvpdt) * fbuffer(PSI,dvtdr) * one_over_r(r) * ref%density(r) &
                - fbuffer(PSI,dvrdr) * fbuffer(PSI,vphi) * one_over_r(r) * ref%density(r) &
                - 2 * fbuffer(PSI,dvpdr) * fbuffer(PSI,vr) * one_over_r(r) * ref%density(r) &
                + d2_fbuffer(PSI,dvrdpdp) * fbuffer(PSI,vphi) * ref%density(r) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                + d2_fbuffer(PSI,dvrdrdp) * fbuffer(PSI,vr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + d2_fbuffer(PSI,dvrdtdp) * fbuffer(PSI,vtheta) * csctheta(t) * ref%density(r) * one_over_r(r) &
                    * one_over_r(r) &
                + fbuffer(PSI,dvpdp) * fbuffer(PSI,dvrdp) * ref%density(r) * csctheta(t) * csctheta(t) * one_over_r(r) &
                    * one_over_r(r) &
                + fbuffer(PSI,dvrdp) * fbuffer(PSI,dvrdr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + fbuffer(PSI,dvrdt) * fbuffer(PSI,dvtdp) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - d2_fbuffer(PSI,dvpdrdp) * fbuffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%density(r) &
                - fbuffer(PSI,dvpdp) * fbuffer(PSI,dvpdr) * csctheta(t) * one_over_r(r) * ref%density(r) &
                - fbuffer(PSI,dvpdr) * fbuffer(PSI,vtheta) * cottheta(t) * one_over_r(r) * ref%density(r) &
                - fbuffer(PSI,dvpdt) * fbuffer(PSI,vtheta) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                - fbuffer(PSI,dvtdr) * fbuffer(PSI,vphi) * cottheta(t) * one_over_r(r) * ref%density(r) &
                - fbuffer(PSI,vphi) * fbuffer(PSI,vr) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                - 2 * fbuffer(PSI,dvpdp) * fbuffer(PSI,vphi) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - 2 * fbuffer(PSI,dvtdp) * fbuffer(PSI,vtheta) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - fbuffer(PSI,dvpdp) * fbuffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                - fbuffer(PSI,vphi) * fbuffer(PSI,vtheta) * cottheta(t) * one_over_r(r) * ref%density(r) * ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vp_grad_vp_phi)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvtdrdr) * fbuffer(PSI,vr) * ref%density(r) &
                + fbuffer(PSI,dvrdr) * fbuffer(PSI,dvtdr) * ref%density(r) &
                + d2_fbuffer(PSI,dvtdrdt) * fbuffer(PSI,vtheta) * one_over_r(r) * ref%density(r) &
                + fbuffer(PSI,dvrdr) * fbuffer(PSI,vtheta) * one_over_r(r) * ref%density(r) &
                + fbuffer(PSI,dvtdr) * fbuffer(PSI,dvtdt) * one_over_r(r) * ref%density(r) &
                + fbuffer(PSI,dvtdr) * fbuffer(PSI,vr) * ref%density(r) * ref%dlnrho(r) &
                - d2_fbuffer(PSI,dvrdrdt) * fbuffer(PSI,vr) * one_over_r(r) * ref%density(r) &
                - d2_fbuffer(PSI,dvrdtdt) * fbuffer(PSI,vtheta) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - fbuffer(PSI,dvrdr) * fbuffer(PSI,dvrdt) * one_over_r(r) * ref%density(r) &
                - fbuffer(PSI,dvrdt) * fbuffer(PSI,dvtdt) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + 2 * fbuffer(PSI,dvpdt) * fbuffer(PSI,vphi) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + 2 * fbuffer(PSI,dvtdr) * fbuffer(PSI,vr) * one_over_r(r) * ref%density(r) &
                + 2 * fbuffer(PSI,dvtdt) * fbuffer(PSI,vtheta) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                + d2_fbuffer(PSI,dvtdrdp) * fbuffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + fbuffer(PSI,dvpdr) * fbuffer(PSI,dvtdp) * csctheta(t) * one_over_r(r) * ref%density(r) &
                + fbuffer(PSI,dvtdt) * fbuffer(PSI,vtheta) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                + fbuffer(PSI,vr) * fbuffer(PSI,vtheta) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) &
                - d2_fbuffer(PSI,dvrdtdp) * fbuffer(PSI,vphi) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - fbuffer(PSI,dvpdt) * fbuffer(PSI,dvrdp) * csctheta(t) * ref%density(r) * one_over_r(r) * one_over_r(r) &
                - cottheta(t) * one_over_r(r) * ref%density(r) * ref%dlnrho(r) * fbuffer(PSI,vphi) * fbuffer(PSI,vphi) &
                - 2 * fbuffer(PSI,dvpdr) * fbuffer(PSI,vphi) * cottheta(t) * one_over_r(r) * ref%density(r) &
                + fbuffer(PSI,dvrdp) * fbuffer(PSI,vphi) * costheta(t) * ref%density(r) * csctheta(t) * csctheta(t) &
                    * one_over_r(r) * one_over_r(r) &
                + fbuffer(PSI,dvtdp) * fbuffer(PSI,vphi) * csctheta(t) * one_over_r(r) * ref%density(r) * ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vm_grad_vm_r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvpdrdt)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                + d2_m0(PSI2,dvpdtdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + m0_values(PSI2,dvpdr)*m0_values(PSI2,dvrdt)*one_over_r(r)*ref%density(r) &
                + m0_values(PSI2,dvpdt)*m0_values(PSI2,dvtdt)*one_over_r(r)**2*ref%density(r) &
                + m0_values(PSI2,dvpdt)*m0_values(PSI2,vr)*one_over_r(r)**2*ref%density(r) &
                + m0_values(PSI2,dvrdt)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                - m0_values(PSI2,vphi)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*m0_values(PSI2,dvpdr)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                + cottheta(t)*m0_values(PSI2,dvtdt)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*m0_values(PSI2,vphi)*m0_values(PSI2,vr)*one_over_r(r)**2*ref%density(r) &
                + 2*cottheta(t)*m0_values(PSI2,dvpdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vm_grad_vm_theta)) Then
            DO_PSI
                qty(PSI) = -d2_m0(PSI2,dvpdrdr)*m0_values(PSI2,vr)*ref%density(r) &
                - m0_values(PSI2,dvpdr)*m0_values(PSI2,dvrdr)*ref%density(r) &
                - d2_m0(PSI2,dvpdrdt)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r) &
                - m0_values(PSI2,dvpdr)*m0_values(PSI2,vr)*ref%density(r)*ref%dlnrho(r) &
                - m0_values(PSI2,dvpdt)*m0_values(PSI2,dvtdr)*one_over_r(r)*ref%density(r) &
                - m0_values(PSI2,dvrdr)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r) &
                - 2*m0_values(PSI2,dvpdr)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                - cottheta(t)*m0_values(PSI2,dvpdr)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r) &
                - cottheta(t)*m0_values(PSI2,dvtdr)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r) &
                - m0_values(PSI2,dvpdt)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - m0_values(PSI2,vphi)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - cottheta(t)*m0_values(PSI2,vphi)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r)*ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vm_grad_vm_phi)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvtdrdr)*m0_values(PSI2,vr)*ref%density(r) &
                + m0_values(PSI2,dvrdr)*m0_values(PSI2,dvtdr)*ref%density(r) &
                + d2_m0(PSI2,dvtdrdt)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r) &
                + m0_values(PSI2,dvrdr)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r) &
                + m0_values(PSI2,dvtdr)*m0_values(PSI2,dvtdt)*one_over_r(r)*ref%density(r) &
                + m0_values(PSI2,dvtdr)*m0_values(PSI2,vr)*ref%density(r)*ref%dlnrho(r) &
                - d2_m0(PSI2,dvrdrdt)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                - d2_m0(PSI2,dvrdtdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                - m0_values(PSI2,dvrdr)*m0_values(PSI2,dvrdt)*one_over_r(r)*ref%density(r) &
                - m0_values(PSI2,dvrdt)*m0_values(PSI2,dvtdt)*one_over_r(r)**2*ref%density(r) &
                + 2*m0_values(PSI2,dvpdt)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + 2*m0_values(PSI2,dvtdr)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                + 2*m0_values(PSI2,dvtdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + m0_values(PSI2,dvtdt)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                + m0_values(PSI2,vr)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - cottheta(t)*m0_values(PSI2,vphi)**2*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - 2*cottheta(t)*m0_values(PSI2,dvpdr)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vp_grad_vm_r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvpdrdt)*fbuffer(PSI,vr)*one_over_r(r)*ref%density(r) &
                + d2_m0(PSI2,dvpdtdt)*fbuffer(PSI,vtheta)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,dvpdt)*m0_values(PSI2,vr)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,dvrdt)*m0_values(PSI2,dvpdr)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvtdt)*m0_values(PSI2,dvpdt)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,vphi)*m0_values(PSI2,dvrdt)*one_over_r(r)**2*ref%density(r) &
                - fbuffer(PSI,vphi)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,dvpdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,vphi)*m0_values(PSI2,dvtdt)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,vphi)*m0_values(PSI2,vr)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,vr)*m0_values(PSI2,dvpdr)*one_over_r(r)*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,vtheta)*m0_values(PSI2,dvpdt)*one_over_r(r)**2*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvrdp)*m0_values(PSI2,dvtdr)*one_over_r(r)*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvtdp)*m0_values(PSI2,dvtdt)*one_over_r(r)**2*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvtdp)*m0_values(PSI2,vr)*one_over_r(r)**2*ref%density(r) &
                + costheta(t)*csctheta(t)**2*fbuffer(PSI,dvpdp)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vp_grad_vm_theta)) Then
            DO_PSI
                qty(PSI) = -d2_m0(PSI2,dvpdrdr)*fbuffer(PSI,vr)*ref%density(r) &
                - fbuffer(PSI,dvrdr)*m0_values(PSI2,dvpdr)*ref%density(r) &
                - d2_m0(PSI2,dvpdrdt)*fbuffer(PSI,vtheta)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,dvpdr)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,dvtdr)*m0_values(PSI2,dvpdt)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,vphi)*m0_values(PSI2,dvrdr)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,vr)*m0_values(PSI2,dvpdr)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,vr)*m0_values(PSI2,dvpdr)*ref%density(r)*ref%dlnrho(r) &
                + csctheta(t)*fbuffer(PSI,dvrdp)*m0_values(PSI2,dvrdr)*one_over_r(r)*ref%density(r) &
                + csctheta(t)*fbuffer(PSI,dvtdp)*m0_values(PSI2,dvrdt)*one_over_r(r)**2*ref%density(r) &
                - cottheta(t)*fbuffer(PSI,dvpdr)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r) &
                - cottheta(t)*fbuffer(PSI,vphi)*m0_values(PSI2,dvtdr)*one_over_r(r)*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvpdp)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvtdp)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                - fbuffer(PSI,vphi)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - fbuffer(PSI,vtheta)*m0_values(PSI2,dvpdt)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - cottheta(t)*fbuffer(PSI,vphi)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r)*ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vp_grad_vm_phi)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvtdrdr)*fbuffer(PSI,vr)*ref%density(r) &
                + fbuffer(PSI,dvrdr)*m0_values(PSI2,dvtdr)*ref%density(r) &
                + d2_m0(PSI2,dvtdrdt)*fbuffer(PSI,vtheta)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvpdt)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,dvtdr)*m0_values(PSI2,dvtdt)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvtdr)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvtdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,vphi)*m0_values(PSI2,dvpdt)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,vr)*m0_values(PSI2,dvtdr)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,vr)*m0_values(PSI2,dvtdr)*ref%density(r)*ref%dlnrho(r) &
                + fbuffer(PSI,vtheta)*m0_values(PSI2,dvrdr)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,vtheta)*m0_values(PSI2,dvtdt)*one_over_r(r)**2*ref%density(r) &
                - d2_m0(PSI2,dvrdrdt)*fbuffer(PSI,vr)*one_over_r(r)*ref%density(r) &
                - d2_m0(PSI2,dvrdtdt)*fbuffer(PSI,vtheta)*one_over_r(r)**2*ref%density(r) &
                - fbuffer(PSI,dvrdt)*m0_values(PSI2,dvrdr)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,dvtdt)*m0_values(PSI2,dvrdt)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,vtheta)*m0_values(PSI2,dvtdt)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                + fbuffer(PSI,vtheta)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - cottheta(t)*fbuffer(PSI,dvpdr)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r) &
                - cottheta(t)*fbuffer(PSI,vphi)*m0_values(PSI2,dvpdr)*one_over_r(r)*ref%density(r) &
                - cottheta(t)*fbuffer(PSI,vphi)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r)*ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vm_grad_vp_r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvpdrdt)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                + d2_fbuffer(PSI,dvpdtdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,dvpdr)*m0_values(PSI2,dvrdt)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvpdt)*m0_values(PSI2,dvtdt)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,dvrdt)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,vr)*m0_values(PSI2,dvpdt)*one_over_r(r)**2*ref%density(r) &
                - fbuffer(PSI,vtheta)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,dvpdr)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,dvpdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,dvtdt)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,vr)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + cottheta(t)*fbuffer(PSI,vtheta)*m0_values(PSI2,dvpdt)*one_over_r(r)**2*ref%density(r) &
                + csctheta(t)*d2_fbuffer(PSI,dvpdtdp)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + csctheta(t)*fbuffer(PSI,dvpdp)*m0_values(PSI2,dvpdt)*one_over_r(r)**2*ref%density(r) &
                - csctheta(t)*d2_fbuffer(PSI,dvtdrdp)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                - csctheta(t)*d2_fbuffer(PSI,dvtdtdp)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvrdp)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                - csctheta(t)**2*d2_fbuffer(PSI,dvtdpdp)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + costheta(t)*csctheta(t)**2*fbuffer(PSI,dvpdp)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vm_grad_vp_theta)) Then
            DO_PSI
                qty(PSI) = -d2_fbuffer(PSI,dvpdrdr)*m0_values(PSI2,vr)*ref%density(r) &
                - fbuffer(PSI,dvpdr)*m0_values(PSI2,dvrdr)*ref%density(r) &
                - d2_fbuffer(PSI,dvpdrdt)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,dvpdr)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,dvpdr)*m0_values(PSI2,vr)*ref%density(r)*ref%dlnrho(r) &
                - fbuffer(PSI,dvpdt)*m0_values(PSI2,dvtdr)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,dvrdr)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,vr)*m0_values(PSI2,dvpdr)*one_over_r(r)*ref%density(r) &
                + csctheta(t)*d2_fbuffer(PSI,dvrdrdp)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                + csctheta(t)*d2_fbuffer(PSI,dvrdtdp)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + csctheta(t)**2*d2_fbuffer(PSI,dvrdpdp)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                - cottheta(t)*fbuffer(PSI,dvtdr)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r) &
                - cottheta(t)*fbuffer(PSI,vtheta)*m0_values(PSI2,dvpdr)*one_over_r(r)*ref%density(r) &
                - csctheta(t)*d2_fbuffer(PSI,dvpdrdp)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvpdp)*m0_values(PSI2,dvpdr)*one_over_r(r)*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvpdp)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvtdp)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                - fbuffer(PSI,dvpdt)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - fbuffer(PSI,vr)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - cottheta(t)*fbuffer(PSI,vtheta)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - csctheta(t)*fbuffer(PSI,dvpdp)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r)*ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_vm_grad_vp_phi)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvtdrdr)*m0_values(PSI2,vr)*ref%density(r) &
                + fbuffer(PSI,dvtdr)*m0_values(PSI2,dvrdr)*ref%density(r) &
                + d2_fbuffer(PSI,dvtdrdt)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvpdt)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,dvrdr)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvtdr)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvtdr)*m0_values(PSI2,vr)*ref%density(r)*ref%dlnrho(r) &
                + fbuffer(PSI,dvtdt)*m0_values(PSI2,dvtdr)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvtdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,vphi)*m0_values(PSI2,dvpdt)*one_over_r(r)**2*ref%density(r) &
                + fbuffer(PSI,vr)*m0_values(PSI2,dvtdr)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,vtheta)*m0_values(PSI2,dvtdt)*one_over_r(r)**2*ref%density(r) &
                - d2_fbuffer(PSI,dvrdrdt)*m0_values(PSI2,vr)*one_over_r(r)*ref%density(r) &
                - d2_fbuffer(PSI,dvrdtdt)*m0_values(PSI2,vtheta)*one_over_r(r)**2*ref%density(r) &
                - fbuffer(PSI,dvrdr)*m0_values(PSI2,dvrdt)*one_over_r(r)*ref%density(r) &
                - fbuffer(PSI,dvrdt)*m0_values(PSI2,dvtdt)*one_over_r(r)**2*ref%density(r) &
                + csctheta(t)*d2_fbuffer(PSI,dvtdrdp)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r) &
                + csctheta(t)*fbuffer(PSI,dvtdp)*m0_values(PSI2,dvpdr)*one_over_r(r)*ref%density(r) &
                + fbuffer(PSI,dvtdt)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                + fbuffer(PSI,vr)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - cottheta(t)*fbuffer(PSI,dvpdr)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r) &
                - cottheta(t)*fbuffer(PSI,vphi)*m0_values(PSI2,dvpdr)*one_over_r(r)*ref%density(r) &
                - csctheta(t)*d2_fbuffer(PSI,dvrdtdp)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                - csctheta(t)*fbuffer(PSI,dvrdp)*m0_values(PSI2,dvpdt)*one_over_r(r)**2*ref%density(r) &
                + costheta(t)*csctheta(t)**2*fbuffer(PSI,dvrdp)*m0_values(PSI2,vphi)*one_over_r(r)**2*ref%density(r) &
                + csctheta(t)*fbuffer(PSI,dvtdp)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r)*ref%dlnrho(r) &
                - cottheta(t)*fbuffer(PSI,vphi)*m0_values(PSI2,vphi)*one_over_r(r)*ref%density(r)*ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Curl_Advection_Force

    Subroutine Compute_Curl_Buoyancy_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t
        
        !!!!!!!!!!!!!!!!!!!!!!!!!! Buoyancy Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

        If (compute_quantity(curl_buoyancy_force_theta) .or. compute_quantity(curl_buoyancy_force_theta_squared)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r) * (csctheta(t) * &
                            one_over_r(r) * buffer(PSI,dtdp))
            END_DO
            If (compute_quantity(curl_buoyancy_force_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_buoyancy_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        If (compute_quantity(curl_buoyancy_force_phi) .or. compute_quantity(curl_buoyancy_force_phi_squared)) Then
            DO_PSI
                qty(PSI) = -ref%Buoyancy_Coeff(r) * ( one_over_r(r) * buffer(PSI,dtdt))
            END_DO
            If (compute_quantity(curl_buoyancy_force_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_buoyancy_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif
        
        If (compute_quantity(curl_buoyancy_force_abs)) Then
            DO_PSI
                qty(PSI) = one_over_r(r) * ((ref%Buoyancy_Coeff(r) * csctheta(t) * buffer(PSI,dtdp))**2 + &
                           (-ref%Buoyancy_Coeff(r) * buffer(PSI,dtdt))**2)**(0.5)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_buoyancy_pforce_theta)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r) * (csctheta(t) * &
                            one_over_r(r) * fbuffer(PSI,dtdp))
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_buoyancy_pforce_phi)) Then
            DO_PSI
                qty(PSI) = -ref%Buoyancy_Coeff(r) * ( one_over_r(r) * fbuffer(PSI,dtdt))
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_buoyancy_mforce_phi)) Then
            DO_PSI
                qty(PSI) = -ref%Buoyancy_Coeff(r) * ( one_over_r(r) * m0_values(PSI2,dtdt))
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Curl_Buoyancy_Force

    Subroutine Compute_Curl_Magnetic_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t
        Real*8 :: jxb_abs_r, jxb_abs_t, jxb_abs_p

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Magnetic Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        If (compute_quantity(curl_j_cross_b_r) .or. compute_quantity(curl_j_cross_b_r_squared)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbpdrdt)*buffer(PSI,br)*one_over_r(r)*ref%Lorentz_Coeff &
                + DDBUFF(PSI,dbpdtdp)*buffer(PSI,bphi)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + DDBUFF(PSI,dbpdtdt)*buffer(PSI,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - DDBUFF(PSI,dbtdpdp)*buffer(PSI,bphi)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                - DDBUFF(PSI,dbtdrdp)*buffer(PSI,br)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - DDBUFF(PSI,dbtdtdp)*buffer(PSI,btheta)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + buffer(PSI,bphi)*buffer(PSI,br)*cottheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - buffer(PSI,bphi)*buffer(PSI,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*buffer(PSI,bphi)*buffer(PSI,dbpdp)*costheta(t)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                + buffer(PSI,bphi)*buffer(PSI,dbrdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + buffer(PSI,bphi)*buffer(PSI,dbtdt)*cottheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + buffer(PSI,br)*buffer(PSI,dbpdr)*cottheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + buffer(PSI,br)*buffer(PSI,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - buffer(PSI,br)*buffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*buffer(PSI,btheta)*buffer(PSI,dbpdt)*cottheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - buffer(PSI,btheta)*buffer(PSI,dbrdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + buffer(PSI,dbpdp)*buffer(PSI,dbpdt)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - buffer(PSI,dbpdp)*buffer(PSI,dbtdp)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                + buffer(PSI,dbpdr)*buffer(PSI,dbrdt)*one_over_r(r)*ref%Lorentz_Coeff &
                + buffer(PSI,dbpdt)*buffer(PSI,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - buffer(PSI,dbrdp)*buffer(PSI,dbtdr)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - buffer(PSI,dbtdp)*buffer(PSI,dbtdt)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(curl_j_cross_b_r)) Call Add_Quantity(qty)
            If (compute_quantity(curl_j_cross_b_r_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif            
        
        Endif
        
        
        If (compute_quantity(curl_j_cross_b_theta) .or. compute_quantity(curl_j_cross_b_theta_squared)) Then
            DO_PSI
                qty(PSI) = -DDBUFF(PSI,dbpdrdr)*buffer(PSI,br)*ref%Lorentz_Coeff &
                - buffer(PSI,dbpdr)*buffer(PSI,dbrdr)*ref%Lorentz_Coeff &
                - DDBUFF(PSI,dbpdrdt)*buffer(PSI,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                - buffer(PSI,bphi)*buffer(PSI,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - buffer(PSI,dbpdt)*buffer(PSI,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - 2*buffer(PSI,br)*buffer(PSI,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + DDBUFF(PSI,dbrdpdp)*buffer(PSI,bphi)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                + DDBUFF(PSI,dbrdrdp)*buffer(PSI,br)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + DDBUFF(PSI,dbrdtdp)*buffer(PSI,btheta)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + buffer(PSI,dbpdp)*buffer(PSI,dbrdp)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                + buffer(PSI,dbrdp)*buffer(PSI,dbrdr)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + buffer(PSI,dbrdt)*buffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - DDBUFF(PSI,dbpdrdp)*buffer(PSI,bphi)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - buffer(PSI,bphi)*buffer(PSI,dbtdr)*cottheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - buffer(PSI,btheta)*buffer(PSI,dbpdr)*cottheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - buffer(PSI,dbpdp)*buffer(PSI,dbpdr)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - 2*buffer(PSI,bphi)*buffer(PSI,dbpdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - 2*buffer(PSI,btheta)*buffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(curl_j_cross_b_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_j_cross_b_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif            

        Endif
        
        If (compute_quantity(curl_j_cross_b_phi) .or. compute_quantity(curl_j_cross_b_phi_squared)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbtdrdr)*buffer(PSI,br)*ref%Lorentz_Coeff &
                + buffer(PSI,dbrdr)*buffer(PSI,dbtdr)*ref%Lorentz_Coeff &
                + DDBUFF(PSI,dbtdrdt)*buffer(PSI,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                + buffer(PSI,btheta)*buffer(PSI,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + buffer(PSI,dbtdr)*buffer(PSI,dbtdt)*one_over_r(r)*ref%Lorentz_Coeff &
                - DDBUFF(PSI,dbrdrdt)*buffer(PSI,br)*one_over_r(r)*ref%Lorentz_Coeff &
                - DDBUFF(PSI,dbrdtdt)*buffer(PSI,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - buffer(PSI,dbrdr)*buffer(PSI,dbrdt)*one_over_r(r)*ref%Lorentz_Coeff &
                - buffer(PSI,dbrdt)*buffer(PSI,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*buffer(PSI,bphi)*buffer(PSI,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*buffer(PSI,br)*buffer(PSI,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + 2*buffer(PSI,btheta)*buffer(PSI,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + DDBUFF(PSI,dbtdrdp)*buffer(PSI,bphi)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + buffer(PSI,dbpdr)*buffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - DDBUFF(PSI,dbrdtdp)*buffer(PSI,bphi)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - buffer(PSI,dbpdt)*buffer(PSI,dbrdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - 2*buffer(PSI,bphi)*buffer(PSI,dbpdr)*cottheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + buffer(PSI,bphi)*buffer(PSI,dbrdp)*costheta(t)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            If (compute_quantity(curl_j_cross_b_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_j_cross_b_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif            

        Endif
        
        If (compute_quantity(curl_j_cross_b_abs)) Then
            DO_PSI
                jxb_abs_r = DDBUFF(PSI,dbpdrdt)*buffer(PSI,br)*one_over_r(r) &
                + DDBUFF(PSI,dbpdtdp)*buffer(PSI,bphi)*csctheta(t)*one_over_r(r)**2 &
                + DDBUFF(PSI,dbpdtdt)*buffer(PSI,btheta)*one_over_r(r)**2 &
                - DDBUFF(PSI,dbtdpdp)*buffer(PSI,bphi)*csctheta(t)**2*one_over_r(r)**2 &
                - DDBUFF(PSI,dbtdrdp)*buffer(PSI,br)*csctheta(t)*one_over_r(r) &
                - DDBUFF(PSI,dbtdtdp)*buffer(PSI,btheta)*csctheta(t)*one_over_r(r)**2 &
                + buffer(PSI,bphi)*buffer(PSI,br)*cottheta(t)*one_over_r(r)**2 &
                - buffer(PSI,bphi)*buffer(PSI,btheta)*one_over_r(r)**2 &
                + 2*buffer(PSI,bphi)*buffer(PSI,dbpdp)*costheta(t)*csctheta(t)**2*one_over_r(r)**2 &
                + buffer(PSI,bphi)*buffer(PSI,dbrdt)*one_over_r(r)**2 &
                + buffer(PSI,bphi)*buffer(PSI,dbtdt)*cottheta(t)*one_over_r(r)**2 &
                + buffer(PSI,br)*buffer(PSI,dbpdr)*cottheta(t)*one_over_r(r) &
                + buffer(PSI,br)*buffer(PSI,dbpdt)*one_over_r(r)**2 &
                - buffer(PSI,br)*buffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)**2 &
                + 2*buffer(PSI,btheta)*buffer(PSI,dbpdt)*cottheta(t)*one_over_r(r)**2 &
                - buffer(PSI,btheta)*buffer(PSI,dbrdp)*csctheta(t)*one_over_r(r)**2 &
                + buffer(PSI,dbpdp)*buffer(PSI,dbpdt)*csctheta(t)*one_over_r(r)**2 &
                - buffer(PSI,dbpdp)*buffer(PSI,dbtdp)*csctheta(t)**2*one_over_r(r)**2 &
                + buffer(PSI,dbpdr)*buffer(PSI,dbrdt)*one_over_r(r) &
                + buffer(PSI,dbpdt)*buffer(PSI,dbtdt)*one_over_r(r)**2 &
                - buffer(PSI,dbrdp)*buffer(PSI,dbtdr)*csctheta(t)*one_over_r(r) &
                - buffer(PSI,dbtdp)*buffer(PSI,dbtdt)*csctheta(t)*one_over_r(r)**2

                jxb_abs_t = -DDBUFF(PSI,dbpdrdr)*buffer(PSI,br) - buffer(PSI,dbpdr)*buffer(PSI,dbrdr) &
                - DDBUFF(PSI,dbpdrdt)*buffer(PSI,btheta)*one_over_r(r) &
                - buffer(PSI,bphi)*buffer(PSI,dbrdr)*one_over_r(r) &
                - buffer(PSI,dbpdt)*buffer(PSI,dbtdr)*one_over_r(r) &
                - 2*buffer(PSI,br)*buffer(PSI,dbpdr)*one_over_r(r) &
                + DDBUFF(PSI,dbrdpdp)*buffer(PSI,bphi)*csctheta(t)**2*one_over_r(r)**2 &
                + DDBUFF(PSI,dbrdrdp)*buffer(PSI,br)*csctheta(t)*one_over_r(r) &
                + DDBUFF(PSI,dbrdtdp)*buffer(PSI,btheta)*csctheta(t)*one_over_r(r)**2 &
                + buffer(PSI,dbpdp)*buffer(PSI,dbrdp)*csctheta(t)**2*one_over_r(r)**2 &
                + buffer(PSI,dbrdp)*buffer(PSI,dbrdr)*csctheta(t)*one_over_r(r) &
                + buffer(PSI,dbrdt)*buffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)**2 &
                - DDBUFF(PSI,dbpdrdp)*buffer(PSI,bphi)*csctheta(t)*one_over_r(r) &
                - buffer(PSI,bphi)*buffer(PSI,dbtdr)*cottheta(t)*one_over_r(r) &
                - buffer(PSI,btheta)*buffer(PSI,dbpdr)*cottheta(t)*one_over_r(r) &
                - buffer(PSI,dbpdp)*buffer(PSI,dbpdr)*csctheta(t)*one_over_r(r) &
                - 2*buffer(PSI,bphi)*buffer(PSI,dbpdp)*csctheta(t)*one_over_r(r)**2 &
                - 2*buffer(PSI,btheta)*buffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)**2

                jxb_abs_p = DDBUFF(PSI,dbtdrdr)*buffer(PSI,br) + buffer(PSI,dbrdr)*buffer(PSI,dbtdr) &
                + DDBUFF(PSI,dbtdrdt)*buffer(PSI,btheta)*one_over_r(r) &
                + buffer(PSI,btheta)*buffer(PSI,dbrdr)*one_over_r(r) &
                + buffer(PSI,dbtdr)*buffer(PSI,dbtdt)*one_over_r(r) &
                - DDBUFF(PSI,dbrdrdt)*buffer(PSI,br)*one_over_r(r) &
                - DDBUFF(PSI,dbrdtdt)*buffer(PSI,btheta)*one_over_r(r)**2 &
                - buffer(PSI,dbrdr)*buffer(PSI,dbrdt)*one_over_r(r) &
                - buffer(PSI,dbrdt)*buffer(PSI,dbtdt)*one_over_r(r)**2 &
                + 2*buffer(PSI,bphi)*buffer(PSI,dbpdt)*one_over_r(r)**2 &
                + 2*buffer(PSI,br)*buffer(PSI,dbtdr)*one_over_r(r) &
                + 2*buffer(PSI,btheta)*buffer(PSI,dbtdt)*one_over_r(r)**2 &
                + DDBUFF(PSI,dbtdrdp)*buffer(PSI,bphi)*csctheta(t)*one_over_r(r) &
                + buffer(PSI,dbpdr)*buffer(PSI,dbtdp)*csctheta(t)*one_over_r(r) &
                - DDBUFF(PSI,dbrdtdp)*buffer(PSI,bphi)*csctheta(t)*one_over_r(r)**2 &
                - buffer(PSI,dbpdt)*buffer(PSI,dbrdp)*csctheta(t)*one_over_r(r)**2 &
                - 2*buffer(PSI,bphi)*buffer(PSI,dbpdr)*cottheta(t)*one_over_r(r) &
                + buffer(PSI,bphi)*buffer(PSI,dbrdp)*costheta(t)*csctheta(t)**2*one_over_r(r)**2

                qty(PSI) = ref%Lorentz_Coeff * sqrt(jxb_abs_r*jxb_abs_r + jxb_abs_t*jxb_abs_t + jxb_abs_p*jxb_abs_p)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jp_cross_bp_r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbpdrdt)*fbuffer(PSI,br)*one_over_r(r)*ref%Lorentz_Coeff &
                + d2_fbuffer(PSI,dbpdtdp)*fbuffer(PSI,bphi)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + d2_fbuffer(PSI,dbpdtdt)*fbuffer(PSI,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbtdpdp)*fbuffer(PSI,bphi)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbtdrdp)*fbuffer(PSI,br)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbtdtdp)*fbuffer(PSI,btheta)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,bphi)*fbuffer(PSI,br)*cottheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,bphi)*fbuffer(PSI,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*fbuffer(PSI,bphi)*fbuffer(PSI,dbpdp)*costheta(t)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,bphi)*fbuffer(PSI,dbrdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,bphi)*fbuffer(PSI,dbtdt)*cottheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,br)*fbuffer(PSI,dbpdr)*cottheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,br)*fbuffer(PSI,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,br)*fbuffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*fbuffer(PSI,btheta)*fbuffer(PSI,dbpdt)*cottheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,btheta)*fbuffer(PSI,dbrdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdp)*fbuffer(PSI,dbpdt)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbpdp)*fbuffer(PSI,dbtdp)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdr)*fbuffer(PSI,dbrdt)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdt)*fbuffer(PSI,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbrdp)*fbuffer(PSI,dbtdr)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbtdp)*fbuffer(PSI,dbtdt)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jp_cross_bp_theta)) Then
            DO_PSI
                qty(PSI) = -d2_fbuffer(PSI,dbpdrdr)*fbuffer(PSI,br)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbpdr)*fbuffer(PSI,dbrdr)*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbpdrdt)*fbuffer(PSI,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,bphi)*fbuffer(PSI,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbpdt)*fbuffer(PSI,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - 2*fbuffer(PSI,br)*fbuffer(PSI,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + d2_fbuffer(PSI,dbrdpdp)*fbuffer(PSI,bphi)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                + d2_fbuffer(PSI,dbrdrdp)*fbuffer(PSI,br)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + d2_fbuffer(PSI,dbrdtdp)*fbuffer(PSI,btheta)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdp)*fbuffer(PSI,dbrdp)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbrdp)*fbuffer(PSI,dbrdr)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbrdt)*fbuffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbpdrdp)*fbuffer(PSI,bphi)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,bphi)*fbuffer(PSI,dbtdr)*cottheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,btheta)*fbuffer(PSI,dbpdr)*cottheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbpdp)*fbuffer(PSI,dbpdr)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - 2*fbuffer(PSI,bphi)*fbuffer(PSI,dbpdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - 2*fbuffer(PSI,btheta)*fbuffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jp_cross_bp_phi)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbtdrdr)*fbuffer(PSI,br)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbrdr)*fbuffer(PSI,dbtdr)*ref%Lorentz_Coeff &
                + d2_fbuffer(PSI,dbtdrdt)*fbuffer(PSI,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,btheta)*fbuffer(PSI,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbtdr)*fbuffer(PSI,dbtdt)*one_over_r(r)*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbrdrdt)*fbuffer(PSI,br)*one_over_r(r)*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbrdtdt)*fbuffer(PSI,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbrdr)*fbuffer(PSI,dbrdt)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbrdt)*fbuffer(PSI,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*fbuffer(PSI,bphi)*fbuffer(PSI,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*fbuffer(PSI,br)*fbuffer(PSI,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + 2*fbuffer(PSI,btheta)*fbuffer(PSI,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + d2_fbuffer(PSI,dbtdrdp)*fbuffer(PSI,bphi)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdr)*fbuffer(PSI,dbtdp)*csctheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbrdtdp)*fbuffer(PSI,bphi)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbpdt)*fbuffer(PSI,dbrdp)*csctheta(t)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - 2*fbuffer(PSI,bphi)*fbuffer(PSI,dbpdr)*cottheta(t)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,bphi)*fbuffer(PSI,dbrdp)*costheta(t)*csctheta(t)**2*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jm_cross_bm_r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbpdrdt)*m0_values(PSI2,br)*one_over_r(r)*ref%Lorentz_Coeff &
                + d2_m0(PSI2,dbpdtdt)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + m0_values(PSI2,bphi)*m0_values(PSI2,dbrdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + m0_values(PSI2,br)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + m0_values(PSI2,dbpdr)*m0_values(PSI2,dbrdt)*one_over_r(r)*ref%Lorentz_Coeff &
                + m0_values(PSI2,dbpdt)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - m0_values(PSI2,bphi)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + cottheta(t)*m0_values(PSI2,bphi)*m0_values(PSI2,br)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + cottheta(t)*m0_values(PSI2,bphi)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + cottheta(t)*m0_values(PSI2,br)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + 2*cottheta(t)*m0_values(PSI2,btheta)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jm_cross_bm_theta)) Then
            DO_PSI
                qty(PSI) = -d2_m0(PSI2,dbpdrdr)*m0_values(PSI2,br)*ref%Lorentz_Coeff &
                - m0_values(PSI2,dbpdr)*m0_values(PSI2,dbrdr)*ref%Lorentz_Coeff &
                - d2_m0(PSI2,dbpdrdt)*m0_values(PSI2,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                - m0_values(PSI2,bphi)*m0_values(PSI2,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - m0_values(PSI2,dbpdt)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - 2*m0_values(PSI2,br)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - cottheta(t)*m0_values(PSI2,bphi)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - cottheta(t)*m0_values(PSI2,btheta)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jm_cross_bm_phi)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbtdrdr)*m0_values(PSI2,br)*ref%Lorentz_Coeff &
                + m0_values(PSI2,dbrdr)*m0_values(PSI2,dbtdr)*ref%Lorentz_Coeff &
                + d2_m0(PSI2,dbtdrdt)*m0_values(PSI2,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                + m0_values(PSI2,btheta)*m0_values(PSI2,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + m0_values(PSI2,dbtdr)*m0_values(PSI2,dbtdt)*one_over_r(r)*ref%Lorentz_Coeff &
                - d2_m0(PSI2,dbrdrdt)*m0_values(PSI2,br)*one_over_r(r)*ref%Lorentz_Coeff &
                - d2_m0(PSI2,dbrdtdt)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - m0_values(PSI2,dbrdr)*m0_values(PSI2,dbrdt)*one_over_r(r)*ref%Lorentz_Coeff &
                - m0_values(PSI2,dbrdt)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*m0_values(PSI2,bphi)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*m0_values(PSI2,br)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + 2*m0_values(PSI2,btheta)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - 2*cottheta(t)*m0_values(PSI2,bphi)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jp_cross_bm_r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbpdrdt)*m0_values(PSI2,br)*one_over_r(r)*ref%Lorentz_Coeff &
                + d2_fbuffer(PSI,dbpdtdt)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,bphi)*m0_values(PSI2,dbrdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdr)*m0_values(PSI2,dbrdt)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdt)*m0_values(PSI2,br)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdt)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,bphi)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + cottheta(t)*fbuffer(PSI,bphi)*m0_values(PSI2,br)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + cottheta(t)*fbuffer(PSI,bphi)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + cottheta(t)*fbuffer(PSI,dbpdr)*m0_values(PSI2,br)*one_over_r(r)*ref%Lorentz_Coeff &
                + csctheta(t)*d2_fbuffer(PSI,dbpdtdp)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*d2_fbuffer(PSI,dbtdrdp)*m0_values(PSI2,br)*one_over_r(r)*ref%Lorentz_Coeff &
                - csctheta(t)*d2_fbuffer(PSI,dbtdtdp)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbrdp)*m0_values(PSI2,dbrdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbtdp)*m0_values(PSI2,br)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbtdp)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)**2*d2_fbuffer(PSI,dbtdpdp)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*cottheta(t)*fbuffer(PSI,dbpdt)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + costheta(t)*csctheta(t)**2*fbuffer(PSI,dbpdp)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jp_cross_bm_theta)) Then
            DO_PSI
                qty(PSI) = -d2_fbuffer(PSI,dbpdrdr)*m0_values(PSI2,br)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbpdr)*m0_values(PSI2,dbrdr)*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbpdrdt)*m0_values(PSI2,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,bphi)*m0_values(PSI2,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbpdt)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - 2*fbuffer(PSI,dbpdr)*m0_values(PSI2,br)*one_over_r(r)*ref%Lorentz_Coeff &
                + csctheta(t)*d2_fbuffer(PSI,dbrdrdp)*m0_values(PSI2,br)*one_over_r(r)*ref%Lorentz_Coeff &
                + csctheta(t)*d2_fbuffer(PSI,dbrdtdp)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + csctheta(t)*fbuffer(PSI,dbrdp)*m0_values(PSI2,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + csctheta(t)*fbuffer(PSI,dbtdp)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + csctheta(t)**2*d2_fbuffer(PSI,dbrdpdp)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - cottheta(t)*fbuffer(PSI,bphi)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - cottheta(t)*fbuffer(PSI,dbpdr)*m0_values(PSI2,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                - csctheta(t)*d2_fbuffer(PSI,dbpdrdp)*m0_values(PSI2,bphi)*one_over_r(r)*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbpdp)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbtdp)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jp_cross_bm_phi)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbtdrdr)*m0_values(PSI2,br)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbtdr)*m0_values(PSI2,dbrdr)*ref%Lorentz_Coeff &
                + d2_fbuffer(PSI,dbtdrdt)*m0_values(PSI2,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,bphi)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,btheta)*m0_values(PSI2,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,btheta)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdr)*m0_values(PSI2,dbpdt)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdt)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbtdr)*m0_values(PSI2,dbtdt)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbtdt)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbrdrdt)*m0_values(PSI2,br)*one_over_r(r)*ref%Lorentz_Coeff &
                - d2_fbuffer(PSI,dbrdtdt)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbpdt)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbrdt)*m0_values(PSI2,dbrdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbrdt)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*fbuffer(PSI,dbtdr)*m0_values(PSI2,br)*one_over_r(r)*ref%Lorentz_Coeff &
                + csctheta(t)*d2_fbuffer(PSI,dbtdrdp)*m0_values(PSI2,bphi)*one_over_r(r)*ref%Lorentz_Coeff &
                + csctheta(t)*fbuffer(PSI,dbtdp)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - cottheta(t)*fbuffer(PSI,bphi)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - cottheta(t)*fbuffer(PSI,dbpdr)*m0_values(PSI2,bphi)*one_over_r(r)*ref%Lorentz_Coeff &
                - csctheta(t)*d2_fbuffer(PSI,dbrdtdp)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbrdp)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + costheta(t)*csctheta(t)**2*fbuffer(PSI,dbrdp)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jm_cross_bp_r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbpdrdt)*fbuffer(PSI,br)*one_over_r(r)*ref%Lorentz_Coeff &
                + d2_m0(PSI2,dbpdtdt)*fbuffer(PSI,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,br)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbrdt)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbrdt)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbtdt)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,btheta)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + cottheta(t)*fbuffer(PSI,br)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + cottheta(t)*fbuffer(PSI,br)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + cottheta(t)*fbuffer(PSI,dbtdt)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + csctheta(t)*fbuffer(PSI,dbpdp)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + csctheta(t)*fbuffer(PSI,dbrdp)*m0_values(PSI2,dbrdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbrdp)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbrdp)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + 2*cottheta(t)*fbuffer(PSI,btheta)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + costheta(t)*csctheta(t)**2*fbuffer(PSI,dbpdp)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jm_cross_bp_theta)) Then
            DO_PSI
                qty(PSI) = -d2_m0(PSI2,dbpdrdr)*fbuffer(PSI,br)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbrdr)*m0_values(PSI2,dbpdr)*ref%Lorentz_Coeff &
                - d2_m0(PSI2,dbpdrdt)*fbuffer(PSI,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbrdr)*m0_values(PSI2,bphi)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbtdr)*m0_values(PSI2,dbpdt)*one_over_r(r)*ref%Lorentz_Coeff &
                - 2*fbuffer(PSI,br)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + csctheta(t)*fbuffer(PSI,dbtdp)*m0_values(PSI2,dbrdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - cottheta(t)*fbuffer(PSI,btheta)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - cottheta(t)*fbuffer(PSI,dbtdr)*m0_values(PSI2,bphi)*one_over_r(r)*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbpdp)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbpdp)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbtdp)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - csctheta(t)*fbuffer(PSI,dbtdp)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_jm_cross_bp_phi)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbtdrdr)*fbuffer(PSI,br)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbrdr)*m0_values(PSI2,dbtdr)*ref%Lorentz_Coeff &
                + d2_m0(PSI2,dbtdrdt)*fbuffer(PSI,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,bphi)*m0_values(PSI2,dbpdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,btheta)*m0_values(PSI2,dbtdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdt)*m0_values(PSI2,bphi)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbpdt)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbrdr)*m0_values(PSI2,btheta)*one_over_r(r)*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbtdt)*m0_values(PSI2,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + fbuffer(PSI,dbtdt)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - d2_m0(PSI2,dbrdrdt)*fbuffer(PSI,br)*one_over_r(r)*ref%Lorentz_Coeff &
                - d2_m0(PSI2,dbrdtdt)*fbuffer(PSI,btheta)*one_over_r(r)**2*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbpdr)*m0_values(PSI2,dbpdt)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbrdr)*m0_values(PSI2,dbrdt)*one_over_r(r)*ref%Lorentz_Coeff &
                - fbuffer(PSI,dbtdt)*m0_values(PSI2,dbrdt)*one_over_r(r)**2*ref%Lorentz_Coeff &
                + 2*fbuffer(PSI,br)*m0_values(PSI2,dbtdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - cottheta(t)*fbuffer(PSI,bphi)*m0_values(PSI2,dbpdr)*one_over_r(r)*ref%Lorentz_Coeff &
                - cottheta(t)*fbuffer(PSI,dbpdr)*m0_values(PSI2,bphi)*one_over_r(r)*ref%Lorentz_Coeff
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Curl_Magnetic_Force

    Subroutine Compute_Curl_Coriolis_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t

        !!!!!!!!!!!!!!!!!!!!!!!!!!!! Coriolis Force !!!!!!!!!!!!!!!!!!!!!!!!!!

        If (compute_quantity(curl_coriolis_force_r) .or. compute_quantity(curl_coriolis_force_r_squared)) Then
            DO_PSI
                qty(PSI) = - ref%Coriolis_Coeff * ref%density(r) * one_over_r(r) * &
                            (-sintheta(t) * buffer(PSI,vtheta) + cottheta(t) * costheta(t) * buffer(PSI,vtheta) + &
                            costheta(t) *  buffer(PSI,dvtdt) + 2 * costheta(t) *  buffer(PSI,vr) + sintheta(t) * &
                            buffer(PSI,dvrdt) + cottheta(t) *  buffer(PSI,dvpdp))
            END_DO
            If (compute_quantity(curl_coriolis_force_r)) Call Add_Quantity(qty)
            If (compute_quantity(curl_coriolis_force_r_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        
        Endif

        If (compute_quantity(curl_coriolis_force_theta) .or. compute_quantity(curl_coriolis_force_theta_squared)) Then
            DO_PSI
                qty(PSI) = ref%Coriolis_Coeff * ref%density(r) * (one_over_r(r) * (buffer(PSI,dvpdp) + &
                            costheta(t) *  buffer(PSI,vtheta) + sintheta(t) * buffer(PSI,vr)) + costheta(t) * &
                            buffer(PSI,dvtdr) + sintheta(t) * buffer(PSI,dvrdr) + ref%dlnrho(r) * costheta(t) * &
                            buffer(PSI,vtheta) + ref%dlnrho(r) * sintheta(t) * buffer(PSI,vr))
            END_DO
            If (compute_quantity(curl_coriolis_force_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_coriolis_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
                
        Endif

        If (compute_quantity(curl_coriolis_force_phi) .or. compute_quantity(curl_coriolis_force_phi_squared)) Then
            DO_PSI
                qty(PSI) = ref%Coriolis_Coeff * ref%density(r) * (ref%dlnrho(r) * costheta(t) * buffer(PSI,vphi) + &
                            costheta(t) * buffer(PSI,dvpdr) - one_over_r(r) * sintheta(t) * buffer(PSI,dvpdt)) 

            END_DO
            If (compute_quantity(curl_coriolis_force_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_coriolis_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        
        Endif

        If (compute_quantity(curl_coriolis_force_abs)) Then
            DO_PSI
                qty(PSI) =((- ref%Coriolis_Coeff * ref%density(r) * one_over_r(r) * &
                            (-sintheta(t) * buffer(PSI,vtheta) + cottheta(t) * costheta(t) * buffer(PSI,vtheta) + &
                            costheta(t) *  buffer(PSI,dvtdt) + 2 * costheta(t) *  buffer(PSI,vr) + sintheta(t) * &
                            buffer(PSI,dvrdt) + cottheta(t) *  buffer(PSI,dvpdp)))**2 + &
                           (ref%Coriolis_Coeff * ref%density(r) * (one_over_r(r) * (buffer(PSI,dvpdp) + &
                            costheta(t) *  buffer(PSI,vtheta) + sintheta(t) * buffer(PSI,vr)) + costheta(t) * &
                            buffer(PSI,dvtdr) + sintheta(t) * buffer(PSI,dvrdr) + ref%dlnrho(r) * costheta(t) * &
                            buffer(PSI,vtheta) + ref%dlnrho(r) * sintheta(t) * buffer(PSI,vr)))**2 + &
                           (ref%Coriolis_Coeff * ref%density(r) * (ref%dlnrho(r) * costheta(t) * buffer(PSI,vphi) + &
                            costheta(t) * buffer(PSI,dvpdr) - one_over_r(r) * sintheta(t) * buffer(PSI,dvpdt)))**2)**(0.5)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_coriolis_pforce_r)) Then
            DO_PSI
                qty(PSI) = -costheta(t)*fbuffer(PSI,dvtdt)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                - cottheta(t)*fbuffer(PSI,dvpdp)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                - fbuffer(PSI,dvrdt)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t) &
                - csctheta(t)*fbuffer(PSI,vtheta)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                - 2*costheta(t)*fbuffer(PSI,vr)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                + 2*fbuffer(PSI,vtheta)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_coriolis_pforce_theta)) Then
            DO_PSI
                qty(PSI) = costheta(t)*fbuffer(PSI,dvtdr)*ref%Coriolis_Coeff*ref%density(r) &
                + fbuffer(PSI,dvpdp)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                + fbuffer(PSI,dvrdr)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t) &
                + costheta(t)*fbuffer(PSI,vtheta)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                + costheta(t)*fbuffer(PSI,vtheta)*ref%Coriolis_Coeff*ref%density(r)*ref%dlnrho(r) &
                + fbuffer(PSI,vr)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t) &
                + fbuffer(PSI,vr)*ref%Coriolis_Coeff*ref%density(r)*ref%dlnrho(r)*sintheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_coriolis_pforce_phi)) Then
            DO_PSI
                qty(PSI) = costheta(t)*fbuffer(PSI,dvpdr)*ref%Coriolis_Coeff*ref%density(r) &
                + costheta(t)*fbuffer(PSI,vphi)*ref%Coriolis_Coeff*ref%density(r)*ref%dlnrho(r) &
                - fbuffer(PSI,dvpdt)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_coriolis_mforce_r)) Then
            DO_PSI
                qty(PSI) = -costheta(t)*m0_values(PSI2,dvtdt)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                - m0_values(PSI2,dvrdt)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t) &
                - csctheta(t)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                - 2*costheta(t)*m0_values(PSI2,vr)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                + 2*m0_values(PSI2,vtheta)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_coriolis_mforce_theta)) Then
            DO_PSI
                qty(PSI) = costheta(t)*m0_values(PSI2,dvtdr)*ref%Coriolis_Coeff*ref%density(r) &
                + m0_values(PSI2,dvrdr)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t) &
                + costheta(t)*m0_values(PSI2,vtheta)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r) &
                + costheta(t)*m0_values(PSI2,vtheta)*ref%Coriolis_Coeff*ref%density(r)*ref%dlnrho(r) &
                + m0_values(PSI2,vr)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t) &
                + m0_values(PSI2,vr)*ref%Coriolis_Coeff*ref%density(r)*ref%dlnrho(r)*sintheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_coriolis_mforce_phi)) Then
            DO_PSI
                qty(PSI) = costheta(t)*m0_values(PSI2,dvpdr)*ref%Coriolis_Coeff*ref%density(r) &
                + costheta(t)*m0_values(PSI2,vphi)*ref%Coriolis_Coeff*ref%density(r)*ref%dlnrho(r) &
                - m0_values(PSI2,dvpdt)*one_over_r(r)*ref%Coriolis_Coeff*ref%density(r)*sintheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Curl_Coriolis_Force

    Subroutine Compute_Curl_Pressure_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t
        ! NOTE: pfactor is assumed constant in the derivation below
        Real*8  :: pfactor(my_r%min:my_r%max)
        pfactor(my_r%min:my_r%max) = ref%dpdr_w_term(my_r%min:my_r%max) &
                                        /ref%density(my_r%min:my_r%max)

        !!!!!!!!!!!!!!!!!!!!!!!!!!! Pressure Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! curl_pressure_force_r is always 0
         
        If (compute_quantity(curl_pressure_force_theta) .or. compute_quantity(curl_pressure_force_theta_squared)) Then
            DO_PSI
                qty(PSI) = pfactor(r) * &
                            One_Over_R(r) * csctheta(t) * ref%dlnrho(r) * buffer(PSI,dpdp)
            END_DO
            If (compute_quantity(curl_pressure_force_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_pressure_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(curl_pressure_force_phi) .or. compute_quantity(curl_pressure_force_phi_squared)) Then
            DO_PSI
                qty(PSI) = - pfactor(r) * &
                              One_Over_R(r) * ref%dlnrho(r) * buffer(PSI,dpdt)
            END_DO
            If (compute_quantity(curl_pressure_force_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_pressure_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        ! curl_pressure_pforce_r is always 0.
        
        If (compute_quantity(curl_pressure_pforce_theta)) Then
            DO_PSI
                qty(PSI) = pfactor(r) * &
                            One_Over_R(r) * csctheta(t) * ref%dlnrho(r) * fbuffer(PSI,dpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_pressure_pforce_phi)) Then
            DO_PSI
                qty(PSI) = - pfactor(r) * &
                              One_Over_R(r) * ref%dlnrho(r) * fbuffer(PSI,dpdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Mean pressure curl: P -> <P>. <P> is axisymmetric (m0_values has no
        ! phi index), so curl_pressure_mforce_theta ~ d<P>/dphi is identically 0.
        ! Only curl_pressure_mforce_phi ~ d<P>/dtheta survives.
        If (compute_quantity(curl_pressure_mforce_phi)) Then
            DO_PSI
                qty(PSI) = - pfactor(r) * &
                              One_Over_R(r) * ref%dlnrho(r) * m0_values(PSI2,dpdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Curl_Pressure_Force

    Subroutine Compute_Curl_Viscous_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t

        !!!!!!!!!!!!!!!!!!!!!!!!!!! Viscous Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! The derivatives of the viscous force F are computed by
        ! Viscous_Force_Derivatives (see the comments there), so that
        !    (curl F)_r     = VFDBUFF(vfd_curl_r)
        !    (curl F)_theta = (1/r) (dF_r/dphi / sin(theta) - F_phi) - dF_phi/dr
        !    (curl F)_phi   = dF_theta/dr + (1/r) (F_theta - dF_r/dtheta)
        ! The second index of each vfd_* array selects the full (1), fluctuating (2)
        ! or mean (3) force.

        If (compute_quantity(curl_viscous_force_r) .or. compute_quantity(curl_viscous_force_r_squared)) Then
            DO_PSI
                qty(PSI) = VFDBUFF(PSI,vfd_curl_r(1))
            END_DO
            If (compute_quantity(curl_viscous_force_r)) Call Add_Quantity(qty)
            If (compute_quantity(curl_viscous_force_r_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(curl_viscous_force_theta) .or. compute_quantity(curl_viscous_force_theta_squared)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*(csctheta(t)*VFDBUFF(PSI,vfd_r_dp(1)) - &
                                vforce_buffer(PSI,vf_p)) - &
                            VFDBUFF(PSI,vfd_p_dr(1))
            END_DO
            If (compute_quantity(curl_viscous_force_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_viscous_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(curl_viscous_force_phi) .or. compute_quantity(curl_viscous_force_phi_squared)) Then
            DO_PSI
                qty(PSI) = VFDBUFF(PSI,vfd_t_dr(1)) + &
                            One_Over_R(r)*(vforce_buffer(PSI,vf_t) - &
                                VFDBUFF(PSI,vfd_r_dt(1)))
            END_DO
            If (compute_quantity(curl_viscous_force_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_viscous_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(curl_viscous_force_abs)) Then
            DO_PSI
                qty(PSI) = sqrt((VFDBUFF(PSI,vfd_curl_r(1)))**2 + &
                                (One_Over_R(r)*(csctheta(t)*VFDBUFF(PSI,vfd_r_dp(1)) - &
                                    vforce_buffer(PSI,vf_p)) - &
                                VFDBUFF(PSI,vfd_p_dr(1)))**2 + &
                                (VFDBUFF(PSI,vfd_t_dr(1)) + &
                                One_Over_R(r)*(vforce_buffer(PSI,vf_t) - &
                                    VFDBUFF(PSI,vfd_r_dt(1))))**2)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Fluctuating and mean viscous force curls: identical in form to
        ! curl_viscous_force_r/theta/phi above, just using the pforce/mforce
        ! vforce_buffer offsets and the second/third set of VFDBUFF indices.
        If (compute_quantity(curl_viscous_pforce_r)) Then
            DO_PSI
                qty(PSI) = VFDBUFF(PSI,vfd_curl_r(2))
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_viscous_pforce_theta)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*(csctheta(t)*VFDBUFF(PSI,vfd_r_dp(2)) - &
                                vforce_buffer(PSI,vfp_p)) - &
                            VFDBUFF(PSI,vfd_p_dr(2))
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_viscous_pforce_phi)) Then
            DO_PSI
                qty(PSI) = VFDBUFF(PSI,vfd_t_dr(2)) + &
                            One_Over_R(r)*(vforce_buffer(PSI,vfp_t) - &
                                VFDBUFF(PSI,vfd_r_dt(2)))
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_viscous_mforce_r)) Then
            DO_PSI
                qty(PSI) = VFDBUFF(PSI,vfd_curl_r(3))
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_viscous_mforce_theta)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*(csctheta(t)*VFDBUFF(PSI,vfd_r_dp(3)) - &
                                vforce_buffer(PSI,vfm_p)) - &
                            VFDBUFF(PSI,vfd_p_dr(3))
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(curl_viscous_mforce_phi)) Then
            DO_PSI
                qty(PSI) = VFDBUFF(PSI,vfd_t_dr(3)) + &
                            One_Over_R(r)*(vforce_buffer(PSI,vfm_t) - &
                                VFDBUFF(PSI,vfd_r_dt(3)))
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Curl_Viscous_Force

    Subroutine Vforce_Derivative_Logic(check, need)
        ! Trigger-code logic shared between the once-at-startup buffer-sizing
        ! pass (check => Sometimes_Compute) and the per-iteration recheck of
        ! whether Viscous_Force_Derivatives needs to run this iteration
        ! (check => Compute_Quantity).
        ! need(c,s) is true if component c (1=r, 2=theta, 3=phi) of the curl
        ! of force set s (1=full, 2=fluctuating, 3=mean) is required.
        Implicit None
        Procedure(Quantity_Check_If) :: check
        Logical, Intent(Out) :: need(3,3)

        need(1,1) = check(curl_viscous_force_r) .or. check(curl_viscous_force_r_squared)
        need(2,1) = check(curl_viscous_force_theta) .or. check(curl_viscous_force_theta_squared)
        need(3,1) = check(curl_viscous_force_phi) .or. check(curl_viscous_force_phi_squared)
        If (check(curl_viscous_force_abs)) need(:,1) = .true.

        need(1,2) = check(curl_viscous_pforce_r)
        need(2,2) = check(curl_viscous_pforce_theta)
        need(3,2) = check(curl_viscous_pforce_phi)

        need(1,3) = check(curl_viscous_mforce_r)
        need(2,3) = check(curl_viscous_mforce_theta)
        need(3,3) = check(curl_viscous_mforce_phi)

    End Subroutine Vforce_Derivative_Logic

    Function Vforce_Derivatives_Needed() result(needed)
        ! Per-iteration determination of whether Viscous_Force_Derivatives needs to
        ! run, via compute_quantity (this iteration's menu) rather than
        ! sometimes_compute, using the same trigger-code logic as
        ! Initialize_Viscous_Force_Derivatives (Vforce_Derivative_Logic).
        Implicit None
        Logical :: needed
        Logical :: need(3,3)

        Call Vforce_Derivative_Logic(Compute_Quantity, need)
        needed = any(need)

    End Function Vforce_Derivatives_Needed

    Subroutine Initialize_Viscous_Force_Derivatives()
        Implicit None
        Integer :: s, nder, nback, n3
        Integer :: vfdfcount(3,2) ! buffer sizes
        Logical :: need(3,3)

        ! Fields Legendre transformed on the way in (p3b slots), per force set s
        ! (see Viscous_Force_Derivatives):
        vf_fr(:) = -1     ! F_r
        vf_sft(:) = -1    ! sin(theta) F_theta
        vf_sfp(:) = -1    ! sin(theta) F_phi
        vf_q(:) = -1      ! Q = r omega_r
        ! Radial derivatives (work array slots):
        vf_sft_dr(:) = -1
        vf_sfp_dr(:) = -1
        vf_q_dr(:) = -1
        vf_q_d2r(:) = -1
        ! Outputs (VFDBUFF slots); the first four are the only fields transposed back:
        vfd_r_dt(:) = -1   ! F_r -> sin(theta) dF_r/dtheta -> dF_r/dtheta
        vfd_t_dr(:) = -1   ! dF_theta/dr
        vfd_p_dr(:) = -1   ! dF_phi/dr
        vfd_curl_r(:) = -1 ! (curl F)_r
        vfd_r_dp(:) = -1   ! dF_r/dphi (FFT only, added in p3a)

        Call Vforce_Derivative_Logic(Sometimes_Compute, need)

        nvftrans = 0
        Do s = 1, 3
            If (need(3,s)) Then
                Call next_slot(nvftrans, vf_fr(s))
                Call next_slot(nvftrans, vf_sft(s))
            Endif
            If (need(2,s)) Call next_slot(nvftrans, vf_sfp(s))
            If (need(1,s)) Call next_slot(nvftrans, vf_q(s))
        Enddo
        If (nvftrans .eq. 0) Return

        nder = nvftrans
        Do s = 1, 3
            If (need(3,s)) Call next_slot(nder, vf_sft_dr(s))
            If (need(2,s)) Call next_slot(nder, vf_sfp_dr(s))
            If (need(1,s)) Then
                Call next_slot(nder, vf_q_dr(s))
                Call next_slot(nder, vf_q_d2r(s))
            Endif
        Enddo
        nvfwork = nder

        nback = 0
        Do s = 1, 3
            If (need(3,s)) Then
                Call next_slot(nback, vfd_r_dt(s))
                Call next_slot(nback, vfd_t_dr(s))
            Endif
            If (need(2,s)) Call next_slot(nback, vfd_p_dr(s))
            If (need(1,s)) Call next_slot(nback, vfd_curl_r(s))
        Enddo
        nvfback = nback

        n3 = nback
        Do s = 1, 3
            If (need(2,s)) Call next_slot(n3, vfd_r_dp(s))
        Enddo

        ! size the buffers at each config stage
        vfdfcount(1,1) = nback    ! config 1a
        vfdfcount(2,1) = nback    ! 2a
        vfdfcount(3,1) = n3       ! 3a
        vfdfcount(3,2) = nvftrans ! 3b
        vfdfcount(2,2) = nvftrans ! 2b
        vfdfcount(1,2) = nvftrans ! 1b

        Call d_vforce_buffer%init(field_count = vfdfcount, config = 'p3b')

        Call d_vforce_buffer%construct('p3a')

        Call d_vforce_buffer%deconstruct('p3a')

    Contains

        Subroutine next_slot(counter, slot)
            Integer, Intent(InOut) :: counter
            Integer, Intent(Out) :: slot
            counter = counter+1
            slot = counter
        End Subroutine next_slot

    End Subroutine Initialize_Viscous_Force_Derivatives

    Subroutine Viscous_Force_Derivatives(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: s, r, k, t, mp, m, imi, lm, l
        Real*8 :: ll1
        Real*8, Allocatable :: mu_visc(:), dmudr(:)
        Real*8, Allocatable :: work(:,:,:,:), work2(:,:,:,:)
        Type(rmcontainer3D), Allocatable :: ddtemp(:)

        ! Computes the derivatives of the viscous force F needed for its curl,
        ! for each requested force set s (1=full, 2=fluctuating, 3=mean).
        !
        ! F_r is a smooth scalar but F_theta and F_phi are not, 
        ! so, as in Compute_Second_Derivatives, they are never
        ! Legendre transformed directly:
        !
        !    (curl F)_theta, (curl F)_phi need dF_r/dtheta, dF_r/dphi, dF_theta/dr and
        !    dF_phi/dr.  F_r and sin(theta) F_h are smooth scalars of degree <= l_max,
        !    so they are transformed exactly; dF_r/dtheta = [sin(theta) dF_r/dtheta]/sin(theta)
        !    (computed spectrally), dF_h/dr = d(sin(theta) F_h)/dr / sin(theta) and
        !    dF_r/dphi is computed with FFTs only.
        !
        !    (curl F)_r needs dF_phi/dtheta, a theta derivative of a horizontal
        !    component.  Instead we use the identity (assuming mu = mu(r))
        !        (curl F)_r = (mu/r) Del^2(Q) + dmu/dr d(Q/r)/dr
        !    where Q = r vort_r = dv_phi/dtheta + cot(theta) v_phi - dv_theta/dphi / sin(theta)
        !    is a smooth scalar built from first derivatives, and
        !        Del^2(Q) = d2Q/dr2 + (2/r) dQ/dr - l(l+1) Q/r^2
        !    is assembled per mode in p1a, with no division by sin(theta).
        !
        ! Derivation of the (curl F)_r identity.  Notation: for any vector A,
        !        rhat.curl(A) = (1/r) C[A],   C[A] = (1/sin(theta)) [d(sin(theta) A_phi)/dtheta - dA_theta/dphi]
        ! where the angular operator C involves no r, so Q = C[v].
        ! 1) With S = e - (1/3) div(v) I (e the strain rate tensor), 2 div(e) = Del^2 v + grad(div v)
        !    and grad(mu) = mu' rhat:
        !        F = div(2 mu S) = mu [Del^2 v + (1/3) grad(div v)] + 2 mu' S.rhat
        !    (Viscous_Force uses div(v) = -v_r dlnrho/dr; the identity holds for any v.)
        ! 2) rhat.curl(mu A) = mu rhat.curl(A), since grad(mu) x A has no r component.
        !    curl(grad(div v)) = 0, and curl(Del^2 v) = Del^2 omega, so the mu term is
        !    mu (Del^2 omega)_r.  For divergence-free A,
        !        (Del^2 A)_r = Del^2 A_r - (2/r^2) A_r - (2/r^2) r div_h(A_h)
        !    with r div_h(A_h) = -(1/r) d(r^2 A_r)/dr, so
        !        (Del^2 A)_r = Del^2 A_r + (2/r) dA_r/dr + (2/r^2) A_r = (1/r) Del^2(r A_r)
        !    (the last step is the product rule Del^2(fg) = f Del^2 g + g Del^2 f + 2 grad f.grad g
        !    with f = r, grad r = rhat, Del^2 r = 2/r).  With A = omega: (mu/r) Del^2(Q).
        ! 3) The (1/3) div(v) rhat part of S.rhat is radial, so has no radial curl.  The rest is
        !        2 e.rhat = d(v)/dr + grad(v_r) - v_h/r
        !    (componentwise: 2 e_rtheta = dv_theta/dr - v_theta/r + (1/r) dv_r/dtheta, etc.), and
        !        C[d(v)/dr]    = dQ/dr      (C involves no r)
        !        C[grad(v_r)]  = 0          (curl of a gradient)
        !        C[v_h/r]      = Q/r
        !    so rhat.curl(2 mu' S.rhat) = mu' (1/r) (dQ/dr - Q/r) = mu' d(Q/r)/dr.
        !
        ! Outline:
        ! 1) Load F_r, sin(theta) F_theta, sin(theta) F_phi and Q (as needed) into p3b
        ! 2) Transform to p1b, copy to the work array and take radial derivatives there
        ! 3) Fill p1a with F_r, d(sin(theta) F_h)/dr and (curl F)_r
        ! 4) In s2a: F_r -> sin(theta) dF_r/dtheta
        ! 5) Transform to p3a, take dF_r/dphi, FFT to physical space, divide out sin(theta)

        Allocate(mu_visc(1:N_R), dmudr(1:N_R))
        mu_visc = ref%density*nu
        dmudr = mu_visc*(ref%dlnrho+dlnu)

        !///////////////////////////////////////////////////////////
        ! Step 1:  Load the fields to be transformed
        call d_vforce_buffer%construct('p3b')
        d_vforce_buffer%config = 'p3b'
        d_vforce_buffer%p3b = 0.0d0

        Do s = 1, 3
            If (vf_fr(s) .gt. 0) Then
                DO_PSI
                    d_vforce_buffer%p3b(PSI,vf_fr(s)) = vforce_buffer(PSI,vf_set(1,s))
                END_DO
            Endif
            If (vf_sft(s) .gt. 0) Then
                DO_PSI
                    d_vforce_buffer%p3b(PSI,vf_sft(s)) = vforce_buffer(PSI,vf_set(2,s))*sintheta(t)
                END_DO
            Endif
            If (vf_sfp(s) .gt. 0) Then
                DO_PSI
                    d_vforce_buffer%p3b(PSI,vf_sfp(s)) = vforce_buffer(PSI,vf_set(3,s))*sintheta(t)
                END_DO
            Endif
        Enddo

        ! Q = r vort_r, from the first derivatives of the matching velocity set
        If (vf_q(1) .gt. 0) Then
            DO_PSI
                d_vforce_buffer%p3b(PSI,vf_q(1)) = buffer(PSI,dvpdt) + cottheta(t)*buffer(PSI,vphi) &
                                                 - csctheta(t)*buffer(PSI,dvtdp)
            END_DO
        Endif
        If (vf_q(2) .gt. 0) Then
            DO_PSI
                d_vforce_buffer%p3b(PSI,vf_q(2)) = fbuffer(PSI,dvpdt) + cottheta(t)*fbuffer(PSI,vphi) &
                                                 - csctheta(t)*fbuffer(PSI,dvtdp)
            END_DO
        Endif
        If (vf_q(3) .gt. 0) Then
            DO_PSI
                d_vforce_buffer%p3b(PSI,vf_q(3)) = m0_values(PSI2,dvpdt) + cottheta(t)*m0_values(PSI2,vphi) &
                                                 - csctheta(t)*m0_values(PSI2,dvtdp)
            END_DO
        Endif

        !///////////////////////////////////////////////////////////
        ! Step 2:  Transform to p1b, copy to the work array and take radial derivatives
        Call fft_to_spectral(d_vforce_buffer%p3b, rsc = .true.)
        call d_vforce_buffer%reform() ! move to p2b
        call d_vforce_buffer%construct('s2b')
        call Legendre_Transform(d_vforce_buffer%p2b, d_vforce_buffer%s2b)
        call d_vforce_buffer%deconstruct('p2b')
        d_vforce_buffer%config = 's2b'

        call d_vforce_buffer%reform() ! move to p1b

        Allocate(work(1:size(d_vforce_buffer%p1b,1),1:2,1:size(d_vforce_buffer%p1b,3),1:nvfwork))
        work(:,:,:,:) = 0.0d0
        if (chebyshev) then
            call gridcp%to_Spectral(d_vforce_buffer%p1b(:,:,:,1:nvftrans), work(:,:,:,1:nvftrans))
            call gridcp%dealias_buffer(work(:,:,:,1:nvftrans))
        else
            work(:,:,:,1:nvftrans) = d_vforce_buffer%p1b(:,:,:,1:nvftrans)
        end if

        Do s = 1, 3
            Call radial_derivative(vf_sft(s), vf_sft_dr(s), 1)
            Call radial_derivative(vf_sfp(s), vf_sfp_dr(s), 1)
            Call radial_derivative(vf_q(s), vf_q_dr(s), 1)
            Call radial_derivative(vf_q(s), vf_q_d2r(s), 2)
        Enddo

        if (chebyshev) then
            Allocate(work2(1:size(work,1),1:2,1:size(work,3),1:nvfwork))
            call gridcp%from_spectral(work, work2)
            work = work2
            DeAllocate(work2)
        end if

        !///////////////////////////////////////////////////////////
        ! Step 3:  Fill p1a with the fields to be transposed back
        call d_vforce_buffer%construct('p1a')
        d_vforce_buffer%config = 'p1a'

        Do s = 1, 3
            If (vfd_r_dt(s) .gt. 0) d_vforce_buffer%p1a(:,:,:,vfd_r_dt(s)) = work(:,:,:,vf_fr(s))
            If (vfd_t_dr(s) .gt. 0) d_vforce_buffer%p1a(:,:,:,vfd_t_dr(s)) = work(:,:,:,vf_sft_dr(s))
            If (vfd_p_dr(s) .gt. 0) d_vforce_buffer%p1a(:,:,:,vfd_p_dr(s)) = work(:,:,:,vf_sfp_dr(s))
            If (vfd_curl_r(s) .gt. 0) Then
                ! (curl F)_r = (mu/r) [Q'' + (2/r) Q' - l(l+1) Q/r^2] + mu' [Q'/r - Q/r^2]
                Do lm = 1, size(work,3)
                    l = l_lm_values(my_lm_min+lm-1)
                    ll1 = l*(l+1.0d0)
                    Do imi = 1, 2
                        Do r = 1, size(work,1)
                            d_vforce_buffer%p1a(r,imi,lm,vfd_curl_r(s)) = &
                                mu_visc(r)*One_Over_R(r)*( work(r,imi,lm,vf_q_d2r(s)) &
                                    + Two_Over_R(r)*work(r,imi,lm,vf_q_dr(s)) &
                                    - ll1*OneOverRSquared(r)*work(r,imi,lm,vf_q(s)) ) &
                              + dmudr(r)*( One_Over_R(r)*work(r,imi,lm,vf_q_dr(s)) &
                                    - OneOverRSquared(r)*work(r,imi,lm,vf_q(s)) )
                        Enddo
                    Enddo
                Enddo
            Endif
        Enddo
        DeAllocate(work)
        call d_vforce_buffer%deconstruct('p1b')

        !///////////////////////////////////////////////////////////
        ! Step 4:  F_r -> sin(theta) dF_r/dtheta
        call d_vforce_buffer%reform() ! now in s2a

        Allocate(ddtemp(my_mp%min:my_mp%max))
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            Allocate(ddtemp(mp)%data(m:l_max,my_r%min:my_r%max,1:2))
        Enddo
        Do s = 1, 3
            If (vfd_r_dt(s) .gt. 0) Then
                ! sin(theta) dF_r/dtheta needs modes up to l_max+1, which s2a cannot
                ! hold.  F_r has no l_max content, but truncate it anyway so that
                ! the result fits exactly and remains divisible by sin(theta).
                Do mp = my_mp%min, my_mp%max
                    d_vforce_buffer%s2a(mp)%data(l_max,:,:,vfd_r_dt(s)) = 0.0d0
                Enddo
                Call d_by_dtheta(d_vforce_buffer%s2a, vfd_r_dt(s), ddtemp)
                DO_IDX2
                    d_vforce_buffer%s2a(mp)%data(IDX2,vfd_r_dt(s)) = ddtemp(mp)%data(IDX2)
                END_DO
            Endif
        Enddo
        Do mp = my_mp%min, my_mp%max
            DeAllocate(ddtemp(mp)%data)
        Enddo
        DeAllocate(ddtemp)

        call d_vforce_buffer%construct('p2a')
        call Legendre_Transform(d_vforce_buffer%s2a, d_vforce_buffer%p2a)
        call d_vforce_buffer%deconstruct('s2a')
        d_vforce_buffer%config = 'p2a'

        !///////////////////////////////////////////////////////////
        ! Step 5:  dF_r/dphi (FFT only), FFT to physical space, divide out sin(theta)
        call d_vforce_buffer%reform() ! move to p3a

        If (size(d_vforce_buffer%p3a,4) .gt. nvfback) Then
            d_vforce_buffer%p3a(:,:,:,nvfback+1:) = 0.0d0
            Do s = 1, 3
                If (vfd_r_dp(s) .gt. 0) Then
                    DO_PSI
                        VFDBUFF(PSI,vfd_r_dp(s)) = vforce_buffer(PSI,vf_set(1,s))
                    END_DO
                Endif
            Enddo
            Call fft_to_spectral_rsc(d_vforce_buffer%p3a(:,:,:,nvfback+1:))
            Do s = 1, 3
                If (vfd_r_dp(s) .gt. 0) Call d_by_dphi(d_vforce_buffer%p3a, vfd_r_dp(s), vfd_r_dp(s))
            Enddo
        Endif

        call FFT_To_Physical(d_vforce_buffer%p3a, rsc=.true.)

        Do s = 1, 3
            Call divide_sintheta(vfd_r_dt(s))
            Call divide_sintheta(vfd_t_dr(s))
            Call divide_sintheta(vfd_p_dr(s))
        Enddo

        ! d_vforce_buffer%p3a (VFDBUFF) now holds, in physical space, for each
        ! force set s (1=full, 2=fluctuating, 3=mean) that needs it:
        !    vfd_r_dt(s)   -- dF_r/dtheta      (needed for (curl F)_phi)
        !    vfd_t_dr(s)   -- dF_theta/dr      (needed for (curl F)_phi)
        !    vfd_p_dr(s)   -- dF_phi/dr        (needed for (curl F)_theta)
        !    vfd_curl_r(s) -- (curl F)_r
        !    vfd_r_dp(s)   -- dF_r/dphi        (needed for (curl F)_theta)
        ! The first four are packed into slots [1 : nvfback] (the only fields
        ! transposed back), and vfd_r_dp(s) follows in [nvfback+1 : ...].
        ! Indices are -1 for fields not needed (see Initialize_Viscous_Force_Derivatives).

        DeAllocate(mu_visc, dmudr)

    Contains

        Subroutine radial_derivative(fin, fout, dorder)
            Integer, Intent(In) :: fin, fout, dorder
            If (fout .le. 0) Return
            if (chebyshev) then
                call gridcp%d_by_dr_cp(fin, fout, work, dorder)
            else
                call d_by_dx3d3(fin, fout, work, dorder)
            end if
        End Subroutine radial_derivative

        Subroutine divide_sintheta(f)
            Integer, Intent(In) :: f
            If (f .le. 0) Return
            DO_PSI
                VFDBUFF(PSI,f) = VFDBUFF(PSI,f)*csctheta(t)
            END_DO
        End Subroutine divide_sintheta

        Integer Function vf_set(c, iset)
            ! vforce_buffer index of component c (1=r, 2=theta, 3=phi) of force set iset
            Integer, Intent(In) :: c, iset
            Integer :: idx(3,3)
            idx(:,1) = [vf_r,  vf_t,  vf_p ]
            idx(:,2) = [vfp_r, vfp_t, vfp_p]
            idx(:,3) = [vfm_r, vfm_t, vfm_p]
            vf_set = idx(c, iset)
        End Function vf_set

    End Subroutine Viscous_Force_Derivatives

End Module Diagnostics_Curl_Momentum
