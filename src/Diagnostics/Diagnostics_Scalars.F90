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
    Implicit None

Contains

    Subroutine Compute_Scalars(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t, ii, scoff


        !////////////////////////////////////////
        !       Scalar
        !  Scalar: field
        Do ii = 1, n_scalar_fields
            scoff = (ii-1)*scalar_offset
            ! compute_quantity also sets the current quantity code in the background as a side effect.
            ! In this case, the code will be set to scalar+scoff
            If (compute_quantity(scalar+scoff)) Then  
                DO_PSI
                    qty(PSI) = buffer(PSI,chivar) ! need to figure out what to do with chivar.  Maybe chivar(ii)?
                END_DO
                Call Add_Quantity(qty)
            Endif

        Enddo
    End Subroutine Compute_Scalars

End Module Diagnostics_Scalars
