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

Module MPI_LAYER
    Use RA_MPI_BASE
    Use All_to_All
    Implicit None
    !Public :: Standard_Transpose




Contains


    Function Init_Main_Group( err) result(grp)
        Type(communicator) :: grp
        Integer, Intent(out) :: err
        ! This sets a communicators (grp's) info to that of mpi_comm_world
         Call mpi_init(err)
        grp%comm = mpi_comm_world
        Call mpi_comm_size(grp%comm, grp%np, err)
        Call mpi_comm_rank(grp%comm, grp%rank, err)

    End Function Init_Main_Group

    Function Init_SubGroup(ingrp,ncpus,ierr) result(grp)
        ! This routine is for multi-run jobs.  It splits ingrp into
        ! an size(ncpus) subgroups, with the ith subgroup containing
        ! ncpus(i) ranks.
        !
        ! ingrp is a communicator with X ranks
        ! ncpus is an integer array arbitrary size <= X
        ! Sum(ncpus) must equal X (number of ranks in ingrp)

        Type(communicator) :: grp
        Integer, Intent(out) :: ierr
        Integer, Intent(In) :: ncpus(1:)
        Type(communicator), Intent(In) :: ingrp
        Integer :: i,grank,gcolor, nclusters, mn_rank, mx_rank
        grank = ingrp%rank
        mn_rank = 0
        nclusters = size(ncpus)

        Do i = 1, nclusters
             mx_rank = ncpus(i)-1+mn_rank
             if ( (grank .ge. mn_rank) .and. (grank .le. mx_rank) ) Then
                gcolor = i-1
             Endif
             mn_rank = mx_rank+1
        Enddo
         Call mpi_comm_split(ingrp%comm, gcolor, ingrp%rank, grp%comm, ierr)

        Call mpi_comm_size(grp%comm, grp%np, ierr)
        Call mpi_comm_rank(grp%comm, grp%rank, ierr)

    End Function Init_SubGroup



    Subroutine RowColSplit(grp,rgrp,cgrp,nprow, err)
        ! Take one group and split it using a row/column decomposition
        Type(communicator), Intent(InOut) :: grp, rgrp, cgrp
        Integer, Intent(out) :: err
        Integer, Intent(In) ::  nprow
        Integer :: row_rank, col_rank

        row_rank = mod(grp%rank,nprow)
        col_rank = grp%rank/nprow

        Call mpi_comm_split(grp%comm, col_rank, grp%rank, rgrp%comm, err)
        Call mpi_comm_split(grp%comm, row_rank, grp%rank, cgrp%comm, err)


        Call mpi_comm_size(rgrp%comm, rgrp%np, err)
        Call mpi_comm_rank(rgrp%comm, rgrp%rank, err)

        Call mpi_comm_size(cgrp%comm, cgrp%np, err)
        Call mpi_comm_rank(cgrp%comm, cgrp%rank, err)

        If (cgrp%rank .ne. col_rank) Write(6,*)'Error - ', cgrp%rank, col_rank
        If (rgrp%rank .ne. row_rank) Write(6,*)'Error - ', rgrp%rank, row_rank
    End Subroutine RowColSplit

    Subroutine Exit_Comm_Lib(err)
        Integer, Intent(out) :: err

        Call mpi_finalize(err)

    End Subroutine Exit_Comm_Lib

    Subroutine Barrier(grp)
        Implicit None
        Type(communicator) :: grp
        Integer :: ierr
        Call MPI_Barrier(grp%comm,ierr)
    End Subroutine Barrier

End Module MPI_LAYER
