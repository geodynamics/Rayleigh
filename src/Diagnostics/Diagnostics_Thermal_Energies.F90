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
!////////////////////// Diagnostics Thermal Energies ///////////////////////
!

!
!///////////////////////////////////////////////////////////////////
Module Diagnostics_Thermal_Energies
    Use Diagnostics_Base
    Use Diagnostics_ADotGradB
    Implicit None

Contains

    Subroutine Compute_Thermal_Energy(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8 :: dt_by_ds, dt_by_dp
        Integer :: r,k, t

        If (compute_quantity(thermal_energy_full) .or. compute_quantity(thermal_energy_sq)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,tvar)*ref%density(r)*ref%temperature(r)
            END_DO
            If (compute_quantity(thermal_energy_full)) Call Add_Quantity(qty)
            If (compute_quantity(thermal_energy_sq)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)**2
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(thermal_energy_p) .or. compute_quantity(thermal_energyp_sq)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,tvar)*ref%density(r)*ref%temperature(r)
            END_DO
            If (compute_quantity(thermal_energy_p)) Call Add_Quantity(qty)
            If (compute_quantity(thermal_energyp_sq)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)**2
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif
        If (compute_quantity(thermal_energy_m) .or. compute_quantity(thermal_energym_sq)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,tvar)*ref%density(r)*ref%temperature(r)
            END_DO
            If (compute_quantity(thermal_energy_m)) Call Add_Quantity(qty)
            If (compute_quantity(thermal_energym_sq)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)**2
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        !////////////////////////////////////////////////////
        ! Enthalpy density  (c_p rho_bar T' = c_p*rho_bar*(dt/ds*s + dt/dp*p)

        If (compute_quantity(enthalpy_full) .or. compute_quantity(enthalpy_sq)) Then
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    dt_by_ds = ref%temperature(r)/pressure_specific_heat
                    dt_by_dp = 1.0d0/pressure_specific_heat/ref%density(r)
                    Do k = 1, n_phi
                        qty(PSI) = pressure_specific_heat*ref%density(r) * &
                            (dt_by_ds*buffer(PSI,tvar) + dt_by_dp*buffer(PSI,pvar))
                    Enddo
                Enddo
            Enddo

            If (compute_quantity(enthalpy_full)) Call Add_Quantity(qty)

            If (compute_quantity(enthalpy_sq)) Then
                DO_PSI
                    qty(PSI)=qty(PSI)**2
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif
        If (compute_quantity(enthalpy_p) .or. compute_quantity(enthalpyp_sq)) Then
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    dt_by_ds = ref%temperature(r)/pressure_specific_heat
                    dt_by_dp = 1.0d0/pressure_specific_heat/ref%density(r)
                    Do k = 1, n_phi
                        qty(PSI) = pressure_specific_heat*ref%density(r) * &
                            (dt_by_ds*fbuffer(PSI,tvar) + dt_by_dp*fbuffer(PSI,pvar))
                    Enddo
                Enddo
            Enddo

            If (compute_quantity(enthalpy_p)) Call Add_Quantity(qty)

            If (compute_quantity(enthalpyp_sq)) Then
                DO_PSI
                    qty(PSI)=qty(PSI)**2
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(enthalpy_m) .or. compute_quantity(enthalpym_sq)) Then
            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    dt_by_ds = ref%temperature(r)/pressure_specific_heat
                    dt_by_dp = 1.0d0/pressure_specific_heat/ref%density(r)
                    Do k = 1, n_phi
                        qty(PSI) = pressure_specific_heat*ref%density(r) * &
                            (dt_by_ds*m0_values(PSI2,tvar) + dt_by_dp*m0_values(PSI2,pvar))
                    Enddo
                Enddo
            Enddo

            If (compute_quantity(enthalpy_m)) Call Add_Quantity(qty)

            If (compute_quantity(enthalpym_sq)) Then
                DO_PSI
                    qty(PSI)=qty(PSI)**2
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif



    End Subroutine Compute_Thermal_Energy

End MOdule Diagnostics_Thermal_Energies
