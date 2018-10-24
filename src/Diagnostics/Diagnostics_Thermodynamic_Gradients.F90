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
Module Diagnostics_Thermodynamic_Gradients
    Use Diagnostics_Base
    Implicit None

Contains

    Subroutine Compute_Thermodynamic_Gradients(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        ! Here we compute the gradient of pressure and entropy/temperature


        !////////////////////////////////////////
        !       Entropy

        !  Entropy: field
        If (compute_quantity(entropy)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Entropy:  radial derivatives
        If (compute_quantity(entropy_dr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_dr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Entropy:  theta derivatives
        If (compute_quantity(entropy_dtheta)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_dtheta)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_dtheta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Entropy:  phi derivatives
        If (compute_quantity(entropy_dphi)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_dphi)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_dphi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Entropy:  (1/r)*theta derivatives
        If (compute_quantity(entropy_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dtdt)*One_Over_R(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dtdt)*One_Over_R(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dtdt)*One_Over_R(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Entropy:  (1/{r sintheta}) * phi derivatives
        If (compute_quantity(entropy_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dtdp)*csctheta(t)*One_Over_R(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dtdp)*csctheta(t)*One_Over_R(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dtdp)*csctheta(t)*One_Over_R(r)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !////////////////////////////////////////
        !       Pressure

        !  pressure: field
        If (compute_quantity(pressure)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,pvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! pressure:  radial derivatives
        If (compute_quantity(pressure_dr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_dr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! pressure:  theta derivatives
        If (compute_quantity(pressure_dtheta)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dpdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_dtheta)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dpdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_dtheta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dpdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! pressure:  phi derivatives
        If (compute_quantity(pressure_dphi)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_dphi)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_dphi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! pressure:  (1/r)*theta derivatives
        If (compute_quantity(pressure_dtr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dpdt)*One_Over_R(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_dtr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dpdt)*One_Over_R(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_dtr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dpdt)*One_Over_R(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! pressure:  (1/{r sintheta}) * phi derivatives
        If (compute_quantity(pressure_dprs)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dpdp)*One_Over_R(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_dprs)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dpdp)*One_Over_R(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_dprs)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dpdp)*One_Over_R(r)*csctheta(t)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !Pressure:  d_by_dr(P/rho_bar)
        If (compute_quantity(rhopressure_dr)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dpdr)-buffer(PSI,pvar)*ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(rhopressurep_dr)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dpdr)-fbuffer(PSI,pvar)*ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(rhopressurem_dr)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dpdr)-m0_values(PSI2,pvar)*ref%dlnrho(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Thermodynamic_Gradients

    Subroutine Compute_Thermal_Second_Derivatives()
        Implicit None
        Integer :: r,k, t

        !//////////////////////////////////////////
        ! Radial second derivatives of
        ! Entropy/Temperature
        If (compute_quantity(entropy_d2r)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dtdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_d2r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dtdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_d2r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dtdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Pressure
        If (compute_quantity(pressure_d2r)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dpdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_d2r)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dpdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_d2r)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dpdrdr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! Theta second derivatives of

        ! Entropy/Temperature
        If (compute_quantity(entropy_d2t)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dtdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_d2t)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dtdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_d2t)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dtdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Pressure
        If (compute_quantity(pressure_d2t)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dpdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_d2t)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dpdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_d2t)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dpdtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! Phi second derivatives of

        ! Entropy/Temperature
        If (compute_quantity(entropy_d2p)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dtdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_d2p)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dtdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_d2p)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dtdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Pressure
        If (compute_quantity(pressure_d2p)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dpdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_d2p)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dpdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_d2p)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dpdpdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! r-Theta second derivatives of

        ! Entropy/Temperature
        If (compute_quantity(entropy_d2rt)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dtdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dtdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dtdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Pressure
        If (compute_quantity(pressure_d2rt)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dpdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dpdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_d2rt)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dpdrdt)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! r-Phi second derivatives of

        ! Entropy/Temperature
        If (compute_quantity(entropy_d2rp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dtdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dtdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dtdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Pressure
        If (compute_quantity(pressure_d2rp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dpdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dpdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_d2rp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dpdrdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !//////////////////////////////////////////
        ! theta-Phi second derivatives of

        ! Entropy/Temperature
        If (compute_quantity(entropy_d2tp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dtdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_p_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dtdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(entropy_m_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dtdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Pressure
        If (compute_quantity(pressure_d2tp)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dpdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_p_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dpdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(pressure_m_d2tp)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dpdtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Thermal_Second_Derivatives
End Module Diagnostics_Thermodynamic_Gradients
