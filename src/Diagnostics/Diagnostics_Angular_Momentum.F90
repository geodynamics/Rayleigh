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
Module Diagnostics_Angular_Momentum
    Use Diagnostics_Base
    Implicit None

!///////////////////////////////////////////////////////////////////////////////////////////////
! In this module, we compute FLUXES for the mean-angular momentum equation
! Everything that governs the evolution of r sin(theta) <v>_phi is here
! SOURCE terms, which are arrived at by multiplying force densities by a moment arm,
! are computed alongside the forces in their respective routines
Contains

    Subroutine Compute_Angular_Momentum_Balance(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)

        Call Compute_Angular_Momentum_Sources(buffer)
        Call Compute_Angular_Momentum_Fluxes(buffer)

    End Subroutine Compute_Angular_Momentum_Balance


    Subroutine Compute_Angular_Momentum_Sources(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

    End Subroutine Compute_Angular_Momentum_Sources



    Subroutine Compute_Angular_Momentum_Fluxes(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t



        If (compute_quantity(famom_fluct_r)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*radius(r)*sintheta(t) &
                    & *(fbuffer(PSI,vr)*fbuffer(PSI,vphi))
            END_DO

            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(famom_fluct_theta)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*radius(r)*sintheta(t) &
                    & *(fbuffer(PSI,vtheta)*fbuffer(PSI,vphi))
            END_DO

            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(famom_dr_r)) Then

            DO_PSI
                    qty(PSI) = ref%density(r)*radius(r)*sintheta(t) &
                        & *(m0_values(PSI2,vr)*m0_values(PSI2,vphi))
            END_DO

            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(famom_dr_theta)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*radius(r)*sintheta(t) &
                    & *(m0_values(PSI2,vtheta)*m0_values(PSI2,vphi))
            END_DO

            Call Add_Quantity(qty)
        Endif

        !These need to be adjusted to handle non-dimensionalization
        If (compute_quantity(famom_mean_r)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*((radius(r)*sintheta(t))**2) &
                    & *(m0_values(PSI2,vr)*Angular_Velocity)
            END_DO

            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(famom_mean_theta)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*((radius(r)*sintheta(t))**2) &
                    & *(m0_values(PSI2,vtheta)*Angular_Velocity)
            END_DO

            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(famom_diff_r)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*sintheta(t)*nu(r)*(m0_values(PSI2,vphi)-radius(r)*m0_values(PSI2,dvpdr))
            END_DO

            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(famom_diff_theta)) Then

            DO_PSI
                qty(PSI) = ref%density(r)*nu(r)*(costheta(t)*m0_values(PSI2,vphi)-sintheta(t)*m0_values(PSI2,dvpdt))
            END_DO

            Call Add_Quantity(qty)
        Endif

        If (magnetism) Then
            If (compute_quantity(famom_maxstr_r)) Then
                DO_PSI
                    qty(PSI) = -radius(r)*sintheta(t)*fbuffer(PSI,br)*fbuffer(PSI,bphi)*ref%Lorentz_Coeff
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(famom_maxstr_theta)) Then
                DO_PSI
                    qty(PSI) = -radius(r)*sintheta(t)*fbuffer(PSI,btheta)*fbuffer(PSI,bphi)*ref%Lorentz_Coeff
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(famom_magtor_r)) Then
                DO_PSI
                    qty(PSI) = -radius(r)*sintheta(t)*m0_values(PSI2,br)*m0_values(PSI2,bphi)*ref%Lorentz_Coeff
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(famom_magtor_theta)) Then
                DO_PSI
                    qty(PSI) = -radius(r)*sintheta(t)*m0_values(PSI2,btheta)*m0_values(PSI2,bphi)*ref%Lorentz_Coeff
                END_DO
                Call Add_Quantity(qty)
            Endif


        Endif

    End Subroutine Compute_Angular_Momentum_Fluxes

End Module Diagnostics_Angular_Momentum
