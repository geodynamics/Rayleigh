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
!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_VELOCITY_FIELD
!               This module computes the components of the velocity field
!               and their derivatives.  Zonal means and fluctuations about
!               those means are also computed.
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Velocity_Field
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_Velocity_Components(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        !/////////////////////////////////////////
        ! 1. terms involving radial velocity
        If (compute_quantity(v_r)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(vp_r)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,vr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(vm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- radial derivatives of v_r
        If (compute_quantity(dv_r_dr)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvrdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_r_dr)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvrdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_r_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- theta derivatives of v_r
        If (compute_quantity(dv_r_dt)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvrdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_r_dt)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvrdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_r_dt)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- phi derivatives of v_r
        If (compute_quantity(dv_r_dp)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvrdp)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_r_dp)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvrdp)
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_r_dp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/r d(v_r)/dtheta}
        If (compute_quantity(dv_r_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvrdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_r_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvrdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_r_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/(r sin[theta]) d(v_r)/dphi
        If (compute_quantity(dv_r_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvrdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_r_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvrdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_r_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvrdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !/////////////////////////////////////////
        ! 2. terms involving theta velocity
        If (compute_quantity(v_theta)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vtheta)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(vp_theta)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,vtheta)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(vm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vtheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- radial derivatives of v_theta
        If (compute_quantity(dv_theta_dr)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvtdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_theta_dr)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvtdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_theta_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- theta derivatives of v_theta
        If (compute_quantity(dv_theta_dt)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvtdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_theta_dt)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvtdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_theta_dt)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- phi derivatives of v_theta
        If (compute_quantity(dv_theta_dp)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvtdp)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_theta_dp)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvtdp)
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_theta_dp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/r d(v_theta)/dtheta}
        If (compute_quantity(dv_theta_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvtdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_theta_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvtdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_theta_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/(r sin[theta]) d(v_theta)/dphi
        If (compute_quantity(dv_theta_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvtdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_theta_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvtdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_theta_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif



        !/////////////////////////////////////////
        ! 3. terms involving phi velocity
        If (compute_quantity(v_phi)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(vp_phi)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,vphi)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(vm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- radial derivatives of v_phi
        If (compute_quantity(dv_phi_dr)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvpdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_phi_dr)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvpdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_phi_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- theta derivatives of v_phi
        If (compute_quantity(dv_phi_dt)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvpdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_phi_dt)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvpdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_phi_dt)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- phi derivatives of v_phi
        If (compute_quantity(dv_phi_dp)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dvpdp)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_phi_dp)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dvpdp)
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_phi_dp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/r d(v_phi)/dtheta}
        If (compute_quantity(dv_phi_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvpdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_phi_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvpdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_phi_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/(r sin[theta]) d(v_phi)/dphi
        If (compute_quantity(dv_phi_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvpdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvp_phi_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvpdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dvm_phi_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvpdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !/////////////////////////////////////////////////////////////
        ! Finally, if desired, we compute the mass flux.

        If (compute_quantity(rhov_r)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vr)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(rhovp_r)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(rhovm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(rhov_theta)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vtheta)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(rhovp_theta)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vtheta)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(rhovm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vtheta)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(rhov_phi)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vphi)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(rhovp_phi)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vphi)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(rhovm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vphi)*ref%density(r)
            END_DO
            Call Add_Quantity(qty)
        Endif


    End Subroutine Compute_Velocity_Components

    Subroutine Compute_Velocity_Second_Derivatives()
        Implicit None
        Integer :: r,k, t

        !//////////////////////////////////////////
        ! Radial second derivatives of

        ! Radial velocity
        If (compute_quantity(dv_r_d2r)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvrdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_r_d2r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvrdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_r_d2r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvrdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(dv_theta_d2r)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvtdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_theta_d2r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvtdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_theta_d2r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvtdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(dv_phi_d2r)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvpdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_phi_d2r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvpdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_phi_d2r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvpdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !//////////////////////////////////////////
        ! Theta second derivatives of

        ! Radial velocity
        If (compute_quantity(dv_r_d2t)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvrdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_r_d2t)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvrdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_r_d2t)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvrdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(dv_theta_d2t)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvtdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_theta_d2t)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvtdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_theta_d2t)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvtdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(dv_phi_d2t)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvpdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_phi_d2t)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvpdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_phi_d2t)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvpdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! phi second derivatives of

        ! Radial velocity
        If (compute_quantity(dv_r_d2p)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvrdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_r_d2p)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvrdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_r_d2p)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvrdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(dv_theta_d2p)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvtdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_theta_d2p)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvtdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_theta_d2p)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvtdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(dv_phi_d2p)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvpdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_phi_d2p)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvpdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_phi_d2p)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvpdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !//////////////////////////////////////////
        ! r-Theta second derivatives of

        ! Radial velocity
        If (compute_quantity(dv_r_d2rt)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvrdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_r_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvrdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_r_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvrdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(dv_theta_d2rt)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvtdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_theta_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvtdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_theta_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvtdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(dv_phi_d2rt)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvpdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_phi_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvpdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_phi_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvpdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! r-Phi second derivatives of

        ! Radial velocity
        If (compute_quantity(dv_r_d2rp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvrdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_r_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvrdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_r_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvrdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(dv_theta_d2rp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvtdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_theta_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvtdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_theta_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvtdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(dv_phi_d2rp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvpdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_phi_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvpdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_phi_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvpdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! Theta-Phi second derivatives of

        ! Radial velocity
        If (compute_quantity(dv_r_d2tp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvrdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_r_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvrdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_r_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvrdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(dv_theta_d2tp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvtdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_theta_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvtdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_theta_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvtdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(dv_phi_d2tp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dvpdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvp_phi_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dvpdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dvm_phi_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dvpdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Velocity_Second_Derivatives

End Module Diagnostics_Velocity_Field
