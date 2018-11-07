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
!               DIAGNOSTICS_CURRENT_DENSITY
!               This module computes the components of del x v and enstrophy.
!               Zonal means and fluctuations about those means are
!               also computed (if desired).
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Vorticity_Field
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_Vorticity_Field(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8, Allocatable :: ens(:,:,:), ens_pm(:,:,:), ens_mm(:,:,:), ens_pp(:,:,:)
        Real*8, Allocatable :: kin_hel(:,:,:), kin_hel_pm(:,:,:), kin_hel_mm(:,:,:)
        Real*8, Allocatable :: kin_hel_pp(:,:,:), kin_hel_mp(:,:,:)
        Logical :: compute_efluct = .false., compute_emean = .false.


        If (compute_quantity(enstrophy)) Then
            Allocate(ens(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif

        If (compute_quantity(enstrophy_mm)) Then
            Allocate(ens_mm(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
            compute_emean = .true.
        Endif

        If (compute_quantity(enstrophy_pm)) Then
            Allocate(ens_pm(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
            compute_emean = .true.
            compute_efluct = .true.
        Endif

        If (compute_quantity(enstrophy_pp)) Then
            Allocate(ens_pp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
            compute_efluct = .true.
        Endif

        If (compute_quantity(kin_helicity)) Then
            Allocate(kin_hel(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif

        If (compute_quantity(kin_helicity_mm)) Then
            Allocate(kin_hel_mm(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif

        If (compute_quantity(kin_helicity_mp)) Then
            Allocate(kin_hel_mp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif

        If (compute_quantity(kin_helicity_pm)) Then
            Allocate(kin_hel_pm(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif

        If (compute_quantity(kin_helicity_pp)) Then
            Allocate(kin_hel_pp(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif
        !/////////////////////////////////////////
        ! 1. terms involving radial vorticity
        If (compute_quantity(kin_helicity_r) .or. compute_quantity(kin_helicity) .or. &
            compute_quantity(vort_r) .or. &
            compute_quantity(enstrophy) .or. compute_quantity(vort_r_sq)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( buffer(PSI,dvpdt) + &
                           cottheta(t)*buffer(PSI,vphi) - &
                           csctheta(t)*buffer(PSI,dvtdp) )
            END_DO
            If (compute_quantity(kin_helicity_r) .or. compute_quantity(kin_helicity)) Then
                DO_PSI
                   tmp1(PSI) = buffer(PSI,vr)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_r)) Then
                   Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity)) Then
                   kin_hel = tmp1
                Endif
            Endif
            If (compute_quantity(vort_r)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(enstrophy)) Then
                ens = qty**2
            Endif
            If (compute_quantity(vort_r_sq)) Then
                qty = qty**2
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(kin_helicity_pp_r) .or. compute_quantity(kin_helicity_mp_r) .or. &
            compute_quantity(kin_helicity_pp) .or. compute_quantity(kin_helicity_mp) .or. &
            compute_quantity(vortp_r) .or. compute_efluct .or. compute_quantity(vortp_r_sq)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( fbuffer(PSI,dvpdt)+ &
                           cottheta(t)*fbuffer(PSI,vphi) - &
                           csctheta(t)*fbuffer(PSI,dvtdp))

            END_DO
            If (compute_quantity(kin_helicity_pp_r) .or. compute_quantity(kin_helicity_pp)) Then
                DO_PSI
                   tmp1(PSI) = fbuffer(PSI,vr)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_pp_r)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_pp)) Then
                    kin_hel_pp = tmp1
                Endif
            Endif
            If (compute_quantity(kin_helicity_mp_r) .or. compute_quantity(kin_helicity_mp)) Then
                DO_PSI
                   tmp1(PSI) = m0_values(PSI2,vr)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_mp_r)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_mp)) Then
                    kin_hel_mp = tmp1
                Endif
            Endif
            If (compute_quantity(vortp_r)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(enstrophy_pp)) Then
                ens_pp = qty**2
            Endif
            If (compute_quantity(enstrophy_pm)) Then
                tmp1 = qty
            Endif
            If (compute_quantity(vortp_r_sq)) Then
                qty = qty**2
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(kin_helicity_mm_r) .or. compute_quantity(kin_helicity_pm_r) .or. &
            compute_quantity(kin_helicity_mm) .or. compute_quantity(kin_helicity_pm) .or. &
            compute_quantity(vortm_r) .or. compute_emean .or. compute_quantity(vortm_r_sq)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( m0_values(PSI2,dvpdt) +&
                           cottheta(t)*m0_values(PSI2,vphi) - &
                           csctheta(t)*m0_values(PSI2,dvtdp) )
            END_DO
            If (compute_quantity(kin_helicity_mm_r) .or. compute_quantity(kin_helicity_mm)) Then
                DO_PSI
                   tmp1(PSI) = m0_values(PSI2,vr)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_mm_r)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_mm)) Then
                    kin_hel_mm = tmp1
                Endif
            Endif
            If (compute_quantity(kin_helicity_pm_r) .or. compute_quantity(kin_helicity_pm)) Then
                DO_PSI
                   tmp1(PSI) = fbuffer(PSI,vr)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_pm_r)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_pm)) Then
                    kin_hel_pm = tmp1
                Endif
            Endif
            If (compute_quantity(vortm_r)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(enstrophy_mm)) Then
                ens_mm = qty**2
            Endif
            If (compute_quantity(enstrophy_pm)) Then
                ens_pm = qty*tmp1
            Endif
            If (compute_quantity(vortm_r_sq)) Then
                qty = qty**2
                Call Add_Quantity(qty)
            Endif
        Endif

        !/////////////////////////////////////////////////
        ! 2. terms involving theta vorticity
        If (compute_quantity(kin_helicity_theta) .or. compute_quantity(kin_helicity) .or. &
            compute_quantity(vort_theta) .or. &
            compute_quantity(enstrophy) .or. compute_quantity(vort_theta_sq)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( csctheta(t)*buffer(PSI,dvrdp) - &
                           buffer(PSI,vphi) )-buffer(PSI,dvpdr)
            END_DO
            If (compute_quantity(kin_helicity_theta) .or. compute_quantity(kin_helicity)) Then
                DO_PSI
                   tmp1(PSI) = buffer(PSI,vtheta)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_theta)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity)) Then
                    kin_hel = kin_hel + tmp1
                Endif
            Endif
            If (compute_quantity(vort_theta)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(enstrophy)) Then
                ens = ens+qty**2
            Endif
            If (compute_quantity(vort_theta_sq)) Then
                qty = qty**2
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(kin_helicity_pp_theta) .or. compute_quantity(kin_helicity_mp_theta) .or. &
            compute_quantity(kin_helicity_pp) .or. compute_quantity(kin_helicity_mp) .or. &
            compute_quantity(vortp_theta) .or. compute_efluct .or. compute_quantity(vortp_theta_sq)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( csctheta(t)*fbuffer(PSI,dvrdp) - &
                           fbuffer(PSI,vphi) )-fbuffer(PSI,dvpdr)
            END_DO
            If (compute_quantity(kin_helicity_pp_theta) .or. compute_quantity(kin_helicity_pp)) Then
                DO_PSI
                   tmp1(PSI) = fbuffer(PSI,vtheta)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_pp_theta)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_pp)) Then
                    kin_hel_pp = kin_hel_pp + tmp1
                Endif
            Endif
            If (compute_quantity(kin_helicity_mp_theta) .or. compute_quantity(kin_helicity_mp)) Then
                DO_PSI
                   tmp1(PSI) = m0_values(PSI2,vtheta)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_mp_theta)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_mp)) Then
                    kin_hel_mp = kin_hel_mp + tmp1
                Endif
            Endif
            If (compute_quantity(vortp_theta)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(enstrophy_pp)) Then
                ens_pp = ens_pp+qty**2
            Endif
            If (compute_quantity(enstrophy_pm)) Then
                tmp1 = qty
            Endif
            If (compute_quantity(vortp_theta_sq)) Then
                qty = qty**2
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(kin_helicity_mm_theta) .or. compute_quantity(kin_helicity_pm_theta) .or. &
            compute_quantity(kin_helicity_mm) .or. compute_quantity(kin_helicity_pm) .or. &
            compute_quantity(vortm_theta) .or. compute_emean .or. compute_quantity(vortm_theta_sq)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*( csctheta(t)*m0_values(PSI2,dvrdp) - &
                           m0_values(PSI2,vphi) )-m0_values(PSI2,dvpdr)
            END_DO

            If (compute_quantity(kin_helicity_mm_theta) .or. compute_quantity(kin_helicity_mm)) Then
                DO_PSI
                   tmp1(PSI) = m0_values(PSI2,vtheta)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_mm_theta)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_mm)) Then
                    kin_hel_mm = kin_hel_mm + tmp1
                Endif
            Endif
            If (compute_quantity(kin_helicity_pm_theta) .or. compute_quantity(kin_helicity_pm)) Then
                DO_PSI
                   tmp1(PSI) = fbuffer(PSI,vtheta)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_pm_theta)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_pm)) Then
                    kin_hel_pm = kin_hel_pm + tmp1
                Endif
            Endif
            If (compute_quantity(vortm_theta)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(enstrophy_mm)) Then
                ens_mm = ens_mm+qty**2
            Endif
            If (compute_quantity(enstrophy_pm)) Then
                ens_pm = ens_pm+qty*tmp1
            Endif
            If (compute_quantity(vortm_theta_sq)) Then
                qty = qty**2
                Call Add_Quantity(qty)
            Endif
        Endif

        !///////////////////////////////////////////
        ! 3. terms involving phi vorticity
        If (compute_quantity(kin_helicity_phi) .or. compute_quantity(kin_helicity) .or. &
            compute_quantity(vort_phi) .or. &
            compute_quantity(enstrophy) .or. compute_quantity(vort_phi_sq)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,dvtdr)+One_Over_R(r)*(buffer(PSI,vtheta) -buffer(PSI,dvrdt))
            END_DO
            If (compute_quantity(kin_helicity_phi) .or. compute_quantity(kin_helicity)) Then
                DO_PSI
                   tmp1(PSI) = buffer(PSI,vphi)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_phi)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity)) Then
                    kin_hel = kin_hel + tmp1
                    Call Add_Quantity(kin_hel)
                Endif
            Endif
            If (compute_quantity(vort_phi)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(enstrophy)) Then
                ens = ens+qty**2
                Call Add_Quantity(ens)
            Endif
            If (compute_quantity(vort_phi_sq)) Then
                qty = qty**2
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(kin_helicity_pp_phi) .or. compute_quantity(kin_helicity_mp_phi) .or. &
            compute_quantity(kin_helicity_pp) .or. compute_quantity(kin_helicity_mp) .or. &
            compute_quantity(vortp_phi) .or. compute_efluct .or. compute_quantity(vortp_phi_sq)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,dvtdr)+One_Over_R(r)*(fbuffer(PSI,vtheta) -fbuffer(PSI,dvrdt))
            END_DO
            If (compute_quantity(kin_helicity_pp_phi) .or. compute_quantity(kin_helicity_pp)) Then
                DO_PSI
                   tmp1(PSI) = fbuffer(PSI,vphi)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_pp_phi)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_pp)) Then
                    kin_hel_pp = kin_hel_pp + tmp1
                    Call Add_Quantity(kin_hel_pp)
                Endif
            Endif
            If (compute_quantity(kin_helicity_mp_phi) .or. compute_quantity(kin_helicity_mp)) Then
                DO_PSI
                   tmp1(PSI) = m0_values(PSI2,vphi)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_mp_phi)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_mp)) Then
                    kin_hel_mp = kin_hel_mp + tmp1
                    Call Add_Quantity(kin_hel_mp)
                Endif
            Endif
            If (compute_quantity(vortp_phi)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(enstrophy_pp)) Then
                ens_pp = ens_pp+qty**2
                Call Add_Quantity(ens_pp)
            Endif
            If (compute_quantity(enstrophy_pm)) Then
                tmp1 = qty
            Endif
            If (compute_quantity(vortp_phi_sq)) Then
                qty = qty**2
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(kin_helicity_mm_phi) .or. compute_quantity(kin_helicity_pm_phi) .or. &
            compute_quantity(kin_helicity_mm) .or. compute_quantity(kin_helicity_pm) .or. &
            compute_quantity(vortm_phi) .or. compute_emean .or. compute_quantity(vortm_phi_sq)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,dvtdr)+One_Over_R(r)*(m0_values(PSI2,vtheta) -m0_values(PSI2,dvrdt))
            END_DO
            If (compute_quantity(kin_helicity_mm_phi) .or. compute_quantity(kin_helicity_mm)) Then
                DO_PSI
                   tmp1(PSI) = m0_values(PSI2,vphi)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_mm_phi)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_mm)) Then
                    kin_hel_mm = kin_hel_mm + tmp1
                    Call Add_Quantity(kin_hel_mm)
                Endif
            Endif
            If (compute_quantity(kin_helicity_pm_phi) .or. compute_quantity(kin_helicity_pm)) Then
                DO_PSI
                   tmp1(PSI) = fbuffer(PSI,vphi)*qty(PSI)
                END_DO
                If (compute_quantity(kin_helicity_pm_phi)) Then
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(kin_helicity_pm)) Then
                    kin_hel_pm = kin_hel_pm + tmp1
                    Call Add_Quantity(kin_hel_pm)
                Endif
            Endif
            If (compute_quantity(vortm_phi)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(enstrophy_mm)) Then
                ens_mm = ens_mm+qty**2
                Call Add_Quantity(ens_mm)
            Endif
            If (compute_quantity(enstrophy_pm)) Then
                ens_pm = ens_pm+qty*tmp1
                Call Add_Quantity(ens_pm)
            Endif
            If (compute_quantity(vortm_phi_sq)) Then
                qty = qty**2
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(zstream)) Then

            DO_PSI
                qty(PSI) = buffer(PSI,zvar)
            END_DO
            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(enstrophy))    DeAllocate(ens)
        If (compute_quantity(enstrophy_mm)) DeAllocate(ens_mm)
        If (compute_quantity(enstrophy_pm)) DeAllocate(ens_pm)
        If (compute_quantity(enstrophy_pp)) DeAllocate(ens_pp)

        If (compute_quantity(kin_helicity))    DeAllocate(kin_hel)
        If (compute_quantity(kin_helicity_mm)) DeAllocate(kin_hel_mm)
        If (compute_quantity(kin_helicity_pm)) DeAllocate(kin_hel_pm)
        If (compute_quantity(kin_helicity_mp)) DeAllocate(kin_hel_mp)
        If (compute_quantity(kin_helicity_pp)) DeAllocate(kin_hel_pp)

    End Subroutine Compute_Vorticity_Field

End Module Diagnostics_Vorticity_Field
