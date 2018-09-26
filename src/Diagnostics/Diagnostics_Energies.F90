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
Module Diagnostics_Energies
    Use Diagnostics_Base
    Implicit None

Contains

    Subroutine Compute_Kinetic_Energy(buffer)
        Implicit None
        Real*8, Allocatable :: dfact(:)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Allocate(dfact(1:N_R))
        !We do a few less multiplies this way, but more importantly,
        ! this leaves less room for typos!

        !Kinetic_energy_factor is a factor that can correct for the 1/2 on
        ! the off chance someone employs an odd nondimensionalization.

        !It's hard to imagine such a situation, but there is parity with
        ! the magnetic_energy_factor this way.

        !Note that kinetic_energy_factor is 1/2 by default
        dfact = ref%density*Half


        ! Energies associated with the full (fluctuating + mean) velocity field
        If (compute_quantity(kinetic_energy)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vr)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vtheta)**2
            DO_PSI
                qty(PSI) = qty(PSI)*dfact(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(radial_ke)) Then
            DO_PSI
                qty(PSI) = dfact(r)*buffer(PSI,vr)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(theta_ke)) Then
            DO_PSI
                qty(PSI) = dfact(r)*buffer(PSI,vtheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(phi_ke)) Then
            DO_PSI
                qty(PSI) = dfact(r)*buffer(PSI,vphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Energies associated with the mean velocity field
        If (compute_quantity(mkinetic_energy)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr    )**2
                qty(PSI) = qty(PSI)+m0_values(PSI2,vtheta)**2
                qty(PSI) = qty(PSI)+m0_values(PSI2,vphi  )**2
                qty(PSI) = qty(PSI)*dfact(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(radial_mke)) Then
            DO_PSI
                qty(PSI) = dfact(r)*m0_values(PSI2,vr)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(theta_mke)) Then
            DO_PSI
                qty(PSI) = dfact(r)*m0_values(PSI2,vtheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(phi_mke)) Then
            DO_PSI
                qty(PSI) = dfact(r)*m0_values(PSI2,vphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Energies associated with the fluctuating velocity field
        If (compute_quantity(pkinetic_energy)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,vphi)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+fbuffer(1:n_phi,:,:,vr)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+fbuffer(1:n_phi,:,:,vtheta)**2
            DO_PSI
                qty(PSI) = qty(PSI)*dfact(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(radial_pke)) Then
            DO_PSI
                qty(PSI) = dfact(r)*fbuffer(PSI,vr)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(theta_pke)) Then
            DO_PSI
                qty(PSI) = dfact(r)*fbuffer(PSI,vtheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(phi_pke)) Then
            DO_PSI
                qty(PSI) = dfact(r)*fbuffer(PSI,vphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! And now the squared fields (without the (1/2) density factor)
        ! Squared full velocity field
        If (compute_quantity(vsq)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,vphi)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vr)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,vtheta)**2

            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(radial_vsq)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vr)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(theta_vsq)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vtheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(phi_vsq)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Squared mean velocity field
        If (compute_quantity(mvsq)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr    )**2
                qty(PSI) = qty(PSI)+m0_values(PSI2,vtheta)**2
                qty(PSI) = qty(PSI)+m0_values(PSI2,vphi  )**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(radial_mvsq)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(theta_mvsq)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vtheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(phi_mvsq)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Squared fluctuating velocity field
        If (compute_quantity(pvsq)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,vphi)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+fbuffer(1:n_phi,:,:,vr)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+fbuffer(1:n_phi,:,:,vtheta)**2

            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(radial_pvsq)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(theta_pvsq)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vtheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(phi_pvsq)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif


        DeAllocate(dfact)

    End Subroutine Compute_Kinetic_Energy


    Subroutine Compute_Magnetic_Energy(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8 :: mfact
        Integer :: r,k, t

        !This changes depending on the non-dimensionalization (dimensionalization) employed
        mfact = half*ref%Lorentz_Coeff


        ! Energies associated with the full (fluctuating + mean) magnetic field
        If (compute_quantity(magnetic_energy)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,bphi)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,br)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+buffer(1:n_phi,:,:,btheta)**2
            DO_PSI
                qty(PSI) = qty(PSI)*mfact
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(radial_me)) Then
            DO_PSI
                qty(PSI) = mfact*buffer(PSI,br)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(theta_me)) Then
            DO_PSI
                qty(PSI) = mfact*buffer(PSI,btheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(phi_me)) Then
            DO_PSI
                qty(PSI) = mfact*buffer(PSI,bphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Energies associated with the mean magnetic field
        If (compute_quantity(mmagnetic_energy)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,br    )**2
                qty(PSI) = qty(PSI)+m0_values(PSI2,btheta)**2
                qty(PSI) = qty(PSI)+m0_values(PSI2,bphi  )**2
                qty(PSI) = qty(PSI)*mfact
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(radial_mme)) Then
            DO_PSI
                qty(PSI) = mfact*m0_values(PSI2,br)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(theta_mme)) Then
            DO_PSI
                qty(PSI) = mfact*m0_values(PSI2,btheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(phi_mme)) Then
            DO_PSI
                qty(PSI) = mfact*m0_values(PSI2,bphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif


        ! Energies associated with the fluctuating magnetic field
        If (compute_quantity(pmagnetic_energy)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,bphi)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+fbuffer(1:n_phi,:,:,br)**2
            qty(1:n_phi,:,:) = qty(1:n_phi,:,:)+fbuffer(1:n_phi,:,:,btheta)**2
            DO_PSI
                qty(PSI) = qty(PSI)*mfact
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(radial_pme)) Then
            DO_PSI
                qty(PSI) = mfact*fbuffer(PSI,br)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(theta_pme)) Then
            DO_PSI
                qty(PSI) = mfact*fbuffer(PSI,btheta)**2
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(phi_pme)) Then
            DO_PSI
                qty(PSI) = mfact*fbuffer(PSI,bphi)**2
            END_DO
            Call Add_Quantity(qty)
        Endif


    End Subroutine Compute_Magnetic_Energy


End Module Diagnostics_Energies
