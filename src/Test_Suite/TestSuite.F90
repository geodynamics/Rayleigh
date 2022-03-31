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

Module TestSuite
    Use Test_SHT
    Use Test_Cheby
    Use ProblemSize
    Implicit None
    Logical :: test_mode = .false.
    Logical :: test_legendre = .false.
    Logical :: test_full_sht = .false.
    Logical :: test_chebyshev = .false.
    Namelist /Test_Namelist/ test_mode, ntest_transform_loops, test_legendre, test_full_sht, test_chebyshev
Contains
    Subroutine Test_Lib()
        If (my_rank .eq. 0) Then
            write(6,*)'Initiating library function tests.'
            write(6,*)
        Endif
        If (test_legendre) Call Test_Legendre_Transforms()
        If (test_full_sht) Call Test_Spherical_Transforms()
        If (test_chebyshev) Call Test_Chebyshev_Transforms()
    End Subroutine Test_Lib


End Module TestSuite
