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

Module Math_Constants
    Implicit None
    Real*8 :: One_Third = 1.0d0/3.0d0
    Real*8 :: Pi  = 3.1415926535897932384626433832795028841972d0
    Real*8 :: four_pi, over_eight_pi, two_pi
    Real*8 :: Half = 0.5d0
    Real*8 :: Zero = 0.0d0
    Real*8 :: one=1.0d0
    Real*8 :: two=2.0d0
Contains
    Subroutine Set_Math_Constants()
        Implicit None
        four_pi = 4.0d0*pi
        two_pi = 2.0d0*pi
        over_eight_pi = 1.0d0/(8.0d0*pi)
    End Subroutine Set_Math_Constants
End Module Math_Constants
