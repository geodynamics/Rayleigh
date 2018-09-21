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

Module MakeDir

    Implicit None

Contains

    !Wrapper to C routine that calls mkdir.
    !System calls are often inaccessible through Fortran, but not through C.
    Subroutine Make_Directory(dirname, ecode)
        character(*), Intent(In) :: dirname
        Integer, Intent(Out) :: ecode
        Integer :: cmkdir
        ecode = cmkdir(Trim(dirname)//char(0))
    End Subroutine Make_Directory

End Module MakeDir
