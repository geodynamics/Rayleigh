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

Module Diagnostics_Miscellaneous
    Use Diagnostics_Base
    Implicit None

    ! Place diagnostics that don't clearly belong in another module here.
Contains

    Subroutine Compute_Misc_Diagnostics(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: mypi, amp
        !////////////////////////////////////////////////////////
        ! Diagnostics for verifying output is working
        If (compute_quantity(diagnostic1)) Then
            mypi = acos(-1.0d0)
            DO_PSI
                qty(PSI) = sin(k*2.0d0*mypi/n_phi) &
                    & *(sintheta(t)**2)*radius(r)
            END_DO
            Write(6,*)'Diagnostic1!', my_rank
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(diagnostic2)) Then
            mypi = acos(-1.0d0)
            DO_PSI
                qty(PSI) = sin(k*4.0d0*mypi/n_phi) &
                    & *(sintheta(t)*costheta(t))*radius(r)**2
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(test_y11)) Then
            amp = -sqrt(3.0d0*over_eight_pi)
            DO_PSI
                qty(PSI) = 2*amp*sintheta(t)*cos((k-1)*two_pi/n_phi)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(test_y22) .or. compute_quantity(test_y22_sq)) Then
            amp = 0.25d0*sqrt(60.0d0*over_eight_pi)
            DO_PSI
                qty(PSI) = 3*amp*(sintheta(t)**2)*cos(2*(k-1)*two_pi/n_phi)
            END_DO
            If (compute_quantity(test_y22)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(test_y22_sq)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif


    End Subroutine Compute_Misc_Diagnostics

End Module Diagnostics_Miscellaneous
