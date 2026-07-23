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
Module Diagnostics_Scalars
    Use Diagnostics_Base
    Use Controls, Only: n_active_scalars, n_passive_scalars
    Implicit None

Contains

    Subroutine Compute_Scalars(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t, ii, ind, scoff, n_scalars
        Integer :: chivar, dchidr, dchidt, dchidp, d2chidr2
        
        n_scalars = n_active_scalars + n_passive_scalars
        Do ii = 1, n_scalars
            if (ii .le. n_active_scalars) then
                 ind   = ii
                 scoff = a_scalar_offset + (ind-1)*scalar_skip
                 chivar   = chiavar(ind)
                 dchidr   = dchiadr(ind)
                 dchidt   = dchiadt(ind)
                 dchidp   = dchiadp(ind)
                 d2chidr2 = d2chiadr2(ind)
            else
                 ind = ii - n_active_scalars
                 scoff = p_scalar_offset + (ind-1)*scalar_skip
                 chivar   = chipvar(ind)
                 dchidr   = dchipdr(ind)
                 dchidt   = dchipdt(ind)
                 dchidp   = dchipdp(ind)
                 d2chidr2 = d2chipdr2(ind)
            end if

            !  Chi: field
            If (compute_quantity(chi+scoff)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,chivar)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(chi_p+scoff)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,chivar)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(chi_m+scoff)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,chivar)
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! Chi:  radial derivatives
            If (compute_quantity(chi_dr+scoff)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,dchidr)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(chi_p_dr+scoff)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,dchidr)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(chi_m_dr+scoff)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,dchidr)
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! Chi:  theta derivatives
            If (compute_quantity(chi_dtheta+scoff)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,dchidt)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(chi_p_dtheta+scoff)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,dchidt)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(chi_m_dtheta+scoff)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,dchidt)
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! Chi:  phi derivatives
            If (compute_quantity(chi_dphi+scoff)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,dchidp)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(chi_p_dphi+scoff)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,dchidp)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(chi_m_dphi+scoff)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,dchidp)
                END_DO
                Call Add_Quantity(qty)
            Endif

            ! only valid for active scalars
            if (ii .le. n_active_scalars) then
                ! NOTE:  The ell=0 component of the r-momentum equation is entirely described by
                !        hydrostatic balance between the pressure and entropy and composition perturbations
                !        (the reference state is assumed to also be in hydrostatic balance).
                !        As such, the ell=0 buoancy force is uninteresting from the point of
                !        view of the flow.  We explicitly separate the ell=0 component for this
                !        term (as with the entropy and pressure terms).
                ! -- full buoyancy
                If (compute_quantity(chi_buoyancy_force+scoff) .or. compute_quantity(chi_buoy_work+scoff)) Then
                    DO_PSI
                        qty(PSI) = ref%chi_buoyancy_coeff(ii,r)*(buffer(PSI,chivar)-&
                                & ell0_values(r,chivar))
                    END_DO
                    If (compute_quantity(chi_buoyancy_force+scoff)) Call Add_Quantity(qty)
                    If (compute_quantity(chi_buoy_work+scoff)) Then
                        DO_PSI
                            qty(PSI)=buffer(PSI,vr)*qty(PSI)
                        END_DO
                        Call Add_Quantity(qty)
                    Endif
                Endif

                ! -- fluctuating buoyancy (ell = 0, m =0 already subtracted)
                If (compute_quantity(chi_buoyancy_pforce+scoff) .or. compute_quantity(chi_buoy_work_pp+scoff)) Then
                    DO_PSI
                        qty(PSI) = ref%chi_buoyancy_coeff(ii,r)*fbuffer(PSI,chivar)
                    END_DO
                    If (compute_quantity(chi_buoyancy_pforce+scoff)) Call Add_Quantity(qty)
                    If (compute_quantity(chi_buoy_work_pp+scoff)) Then
                        DO_PSI
                            qty(PSI)=buffer(PSI,vr)*qty(PSI)
                        END_DO
                        Call Add_Quantity(qty)
                    Endif
                Endif

                ! -- mean buoyancy
                If (compute_quantity(chi_buoyancy_mforce+scoff) .or. compute_quantity(chi_buoy_work_mm+scoff)) Then
                    DO_PSI
                        qty(PSI) = ref%chi_buoyancy_coeff(ii,r)*(m0_values(PSI2,chivar)-&
                                & ell0_values(r,chivar))
                    END_DO
                    If (compute_quantity(chi_buoyancy_mforce+scoff)) Call Add_Quantity(qty)
                    If (compute_quantity(chi_buoy_work_mm+scoff)) Then
                        DO_PSI
                            qty(PSI)=buffer(PSI,vr)*qty(PSI)
                        END_DO
                        Call Add_Quantity(qty)
                    Endif
                Endif

                If (compute_quantity(chi_buoyancy_force_ell0+scoff)) Then
                    DO_PSI
                        qty(PSI) = ref%chi_buoyancy_coeff(ii,r)*ell0_values(r,chivar)
                    END_DO
                    Call Add_Quantity(qty)
                Endif

                ! Chi Buoyant Production of turbulent kinetic energy.
                If (compute_quantity(chi_production_buoyant_pKE+scoff)) Then
                    DO_PSI
                        qty(PSI) = ref%chi_buoyancy_coeff(ii,r)*fbuffer(PSI,chivar)*fbuffer(PSI,vr)
                    END_DO
                    Call Add_Quantity(qty)
                Endif
            endif
            
        Enddo
    End Subroutine Compute_Scalars

End Module Diagnostics_Scalars
