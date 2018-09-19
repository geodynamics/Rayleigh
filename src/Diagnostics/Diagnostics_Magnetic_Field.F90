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
!               DIAGNOSTICS_Magnetic_FIELD
!               This module computes the components of the magnetic field
!               and their derivatives.  Zonal means and fluctuations about
!               those means are also computed.
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Magnetic_Field
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_BField_Components(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        !/////////////////////////////////////////
        ! 1. terms involving radial Bfield
        If (compute_quantity(b_r)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,br)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(bp_r)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,br)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(bm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,br)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- radial derivatives of b_r
        If (compute_quantity(db_r_dr)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dbrdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_r_dr)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dbrdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_r_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- theta derivatives of b_r
        If (compute_quantity(db_r_dt)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dbrdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_r_dt)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dbrdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_r_dt)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- phi derivatives of b_r
        If (compute_quantity(db_r_dp)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dbrdp)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_r_dp)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dbrdp)
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_r_dp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/r d(b_r)/dtheta}
        If (compute_quantity(db_r_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dbrdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_r_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dbrdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_r_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbrdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/(r sin[theta]) d(b_r)/dphi
        If (compute_quantity(db_r_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dbrdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_r_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dbrdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_r_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbrdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !/////////////////////////////////////////
        ! 2. terms involving theta Bfield
        If (compute_quantity(b_theta)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,btheta)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(bp_theta)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,btheta)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(bm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,btheta)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- radial derivatives of b_theta
        If (compute_quantity(db_theta_dr)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dbtdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_theta_dr)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dbtdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_theta_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- theta derivatives of b_theta
        If (compute_quantity(db_theta_dt)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dbtdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_theta_dt)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dbtdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_theta_dt)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- phi derivatives of b_theta
        If (compute_quantity(db_theta_dp)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dbtdp)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_theta_dp)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dbtdp)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_theta_dp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !-- -- {1/r d(b_theta)/dtheta}
        If (compute_quantity(db_theta_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dbtdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_theta_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dbtdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_theta_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbtdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/(r sin[theta]) d(b_theta)/dphi
        If (compute_quantity(db_theta_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dbtdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_theta_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dbtdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_theta_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbtdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif



        !/////////////////////////////////////////
        ! 3. terms involving phi Bfield
        If (compute_quantity(b_phi)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,bphi)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(bp_phi)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,bphi)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(bm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,bphi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- radial derivatives of b_phi
        If (compute_quantity(db_phi_dr)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dbpdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_phi_dr)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dbpdr)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_phi_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- theta derivatives of b_phi
        If (compute_quantity(db_phi_dt)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dbpdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_phi_dt)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dbpdt)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_phi_dt)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbpdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- phi derivatives of b_phi
        If (compute_quantity(db_phi_dp)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,dbpdp)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_phi_dp)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,dbpdp)
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_phi_dp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/r d(b_phi)/dtheta}
        If (compute_quantity(db_phi_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dbpdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_phi_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dbpdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_phi_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbpdt)*one_over_r(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !-- -- {1/(r sin[theta]) d(b_phi)/dphi
        If (compute_quantity(db_phi_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dbpdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbp_phi_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dbpdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(dbm_phi_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dbpdp)*one_over_r(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif


    End Subroutine Compute_BField_Components

    Subroutine Compute_Magnetic_Second_Derivatives()
        Implicit None
        Integer :: r,k, t

        !//////////////////////////////////////////
        ! Radial second derivatives of

        ! Radial velocity
        If (compute_quantity(db_r_d2r)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbrdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_r_d2r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbrdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_r_d2r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbrdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(db_theta_d2r)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbtdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_theta_d2r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbtdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_theta_d2r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbtdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(db_phi_d2r)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbpdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_phi_d2r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbpdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_phi_d2r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbpdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !//////////////////////////////////////////
        ! Theta second derivatives of

        ! Radial velocity
        If (compute_quantity(db_r_d2t)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbrdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_r_d2t)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbrdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_r_d2t)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbrdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(db_theta_d2t)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbtdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_theta_d2t)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbtdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_theta_d2t)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbtdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(db_phi_d2t)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbpdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_phi_d2t)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbpdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_phi_d2t)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbpdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! phi second derivatives of

        ! Radial velocity
        If (compute_quantity(db_r_d2p)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbrdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_r_d2p)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbrdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_r_d2p)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbrdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(db_theta_d2p)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbtdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_theta_d2p)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbtdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_theta_d2p)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbtdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(db_phi_d2p)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbpdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_phi_d2p)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbpdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_phi_d2p)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbpdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !//////////////////////////////////////////
        ! r-Theta second derivatives of

        ! Radial velocity
        If (compute_quantity(db_r_d2rt)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbrdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_r_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbrdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_r_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbrdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(db_theta_d2rt)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbtdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_theta_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbtdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_theta_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbtdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(db_phi_d2rt)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbpdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_phi_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbpdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_phi_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbpdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! r-Phi second derivatives of

        ! Radial velocity
        If (compute_quantity(db_r_d2rp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbrdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_r_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbrdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_r_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbrdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(db_theta_d2rp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbtdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_theta_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbtdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_theta_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbtdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(db_phi_d2rp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbpdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_phi_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbpdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_phi_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbpdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! Theta-Phi second derivatives of

        ! Radial velocity
        If (compute_quantity(db_r_d2tp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbrdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_r_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbrdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_r_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbrdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta velocity
        If (compute_quantity(db_theta_d2tp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbtdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_theta_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbtdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_theta_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbtdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi velocity
        If (compute_quantity(db_phi_d2tp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dbpdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbp_phi_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dbpdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(dbm_phi_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dbpdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Magnetic_Second_Derivatives

End Module Diagnostics_Magnetic_Field
