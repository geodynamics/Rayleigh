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
        Integer :: r,k, t, ii, ind, scoff
        Integer :: chivar, dchidr, dchidt, dchidp, d2chidr2
        
        Do ii = 1, n_active_scalars + n_passive_scalars
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
            ! compute the field
            If (compute_quantity(chi+scoff)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,chivar)
                END_DO
                Call Add_Quantity(qty)
            Endif
            ! compute the dr gradient
            If (compute_quantity(chi_dr+scoff)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,dchidr)
                END_DO
                Call Add_Quantity(qty)
            Endif
            ! compute the dtheta gradient
            If (compute_quantity(chi_dtheta+scoff)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,dchidt)
                END_DO
                Call Add_Quantity(qty)
            Endif
            ! compute the dphi gradient
            If (compute_quantity(chi_dphi+scoff)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,dchidp)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Enddo
    End Subroutine Compute_Scalars

End Module Diagnostics_Scalars
