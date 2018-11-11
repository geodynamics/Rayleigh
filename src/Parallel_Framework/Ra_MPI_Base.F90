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

Module Ra_MPI_Base
    Use MPI

    Type communicator
        Integer :: comm    ! The mpi handle for this group
        Integer :: np    ! The number of processors in this group
        Integer :: rank ! A processes's local rank within this group
    End Type communicator

    Private ::  mpi_null_delete_fn, mpi_dup_fn, mpi_null_copy_fn
    Public :: mpi_wtime,mpi_wtick
End Module Ra_MPI_Base
