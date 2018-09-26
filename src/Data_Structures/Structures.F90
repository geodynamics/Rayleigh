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

Module Structures
    !//////////////////////////////////////////////////////////
    ! This module contains various data structures that do not
    !   have associated methods.
    Type, Public :: rmcontainer1d
        Real*8, Allocatable :: data(:)
    End Type rmcontainer1d

    Type, Public :: rmcontainer
        Real*8, Allocatable :: data(:,:)
    End Type rmcontainer

    Type, Public :: rmcontainer3d
        Real*8, Allocatable :: data(:,:,:)
    End Type rmcontainer3d

    Type, Public :: rmcontainer4d
        Real*8, Allocatable :: data(:,:,:,:)
    End Type rmcontainer4d

End Module Structures
