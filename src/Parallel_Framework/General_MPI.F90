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

!////////////////////////////////////////////////////////////////////////////
!   MODULE:  GENERAL_MPI
!
!   DESCRIPTION:  Provides wrappers for several common, but cumbersome to code,
!               MPI global operations.
!
!   MEMBER SUBROUTINES:
!       Global_MAX
!       Global_IMAX
!       DSUM1D
!       DSUM2D
!       DALLSUM1D
!       DALLSUM2D
!       BCAST2D
!////////////////////////////////////////////////////////////////////////////
Module General_MPI
    Use RA_MPI_BASE
    Implicit None
Contains

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  Global_MAX
    !
    ! DESCRIPTION:  Computes the maximum value of double-precision variable sendbuf
    !                 across all processes in grp.  Results are stored in recvbuf.
    !                 uses MPI_ALLREDUCE
    !
    ! INPUTS:
    !            sendbuf - Double precision variable to compute maximum value of across grp
    !                grp - MPI communicator across which the REDUCE is conducted
    !                       (optional; default = MPI_COMM_WORLD)
    ! OUTPUTS:
    !
    !           recvbuf - Double precision variable that stores the result
    !                       of the REDUCE operation
    !/////////////////////////////////////////////////////////////////////

    Subroutine Global_Max(sendbuf, recvbuf, grp)
        Real*8, Intent(In)  :: sendbuf
        Real*8, Intent(Out) :: recvbuf
        Type(communicator), Optional :: grp
        Integer :: icount,  comm
        Integer :: MPI_err

        icount = 1

        If (Present(grp)) Then
            comm = grp%comm
        Else
            comm = MPI_COMM_WORLD
        End If


        Call MPI_ALLREDUCE(sendbuf, recvbuf, icount, MPI_DOUBLE_PRECISION, MPI_MAX, comm, MPI_err)

    End Subroutine Global_Max

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  Global_IMAX
    !
    ! DESCRIPTION:  Computes the maximum value of 4-byte integer variable sendbuf
    !                 across all processes in grp.  Results are stored in recvbuf.
    !                 uses MPI_ALLREDUCE
    !
    ! INPUTS:
    !            sendbuf - Integer variable (4-byte) to compute maximum value of across grp
    !                grp - MPI communicator across which the REDUCE is conducted
    !                       (optional; default = MPI_COMM_WORLD)
    ! OUTPUTS:
    !
    !           recvbuf - Integer variable (4-byte) that stores the result
    !                       of the REDUCE operation
    !/////////////////////////////////////////////////////////////////////
    Subroutine Global_IMax(sendbuf, recvbuf, grp)
        Integer*4, Intent(In)  :: sendbuf
        Integer*4, Intent(Out) :: recvbuf
        Type(communicator), Optional :: grp
        Integer :: icount,  comm, MPI_err

        icount = 1

        If (Present(grp)) Then
            comm = grp%comm
        Else
            comm = MPI_COMM_WORLD
        End If


        Call MPI_ALLREDUCE(sendbuf, recvbuf, icount, MPI_INTEGER, MPI_MAX, comm, MPI_err)

    End Subroutine Global_IMax

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  DSUM1D
    !
    ! DESCRIPTION:  Performs MPI REDUCE (SUM) of 1-D array
    !                 sendbuf across grp.  Results are stored in rank ddest.
    !
    ! INPUTS:
    !            sendbuf - 1-D double-precision array to SUM across grp
    !                grp - MPI communicator across which the REDUCE is conducted
    !                       (optional; default = MPI_COMM_WORLD)
    !              ddest - MPI rank within grp that stores the results of the REDUCE operation
    ! OUTPUTS:
    !
    !           recvbuf - 1-D double-precision array that stores the result
    !                       of the REDUCE on rank ddest
    !/////////////////////////////////////////////////////////////////////
    Subroutine DSUM1D(sendbuf, recvbuf, grp, ddest)
        Real*8, Intent(In)  :: sendbuf(:)
        Real*8, Intent(Out) :: recvbuf(:)
        Type(communicator), Optional :: grp
        Integer, Intent(In), Optional :: ddest
        Integer :: icount,  comm, MPI_err, dest
        If (present(ddest)) then
            dest = ddest
        Else
            dest = 0
        Endif

        icount = size(sendbuf)

        If (Present(grp)) Then
            comm = grp%comm
        Else
            comm = MPI_COMM_WORLD
        End If

        Call MPI_REDUCE(sendbuf, recvbuf, icount, MPI_DOUBLE_PRECISION, MPI_SUM, dest, comm, MPI_err)

    End Subroutine DSUM1D

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  DALLSUM1D
    !
    ! DESCRIPTION:  Performs MPI AllReduce (SUM) of 1-D array
    !                 sendbuf across grp
    !
    ! INPUTS:
    !            sendbuf - 1-D double-precision array to SUM across grp
    !                grp - MPI communicator across which the ALLREDUCE is conducted
    !                    (optional; default = MPI_COMM_WORLD)
    ! OUTPUTS:
    !
    !           recvbuf - 2-D double-precision array that stores the result
    !                       of the ALLREDUCE
    !/////////////////////////////////////////////////////////////////////
    Subroutine DALLSUM1D(sendbuf, recvbuf, grp)
        Real*8 :: sendbuf(1:)
        Real*8, Intent(Out) :: recvbuf(1:)
        Type(communicator), Optional :: grp
        Integer :: icount,  comm, MPI_err

        icount = size(sendbuf)

        If (Present(grp)) Then
            comm = grp%comm
        Else
            comm = MPI_COMM_WORLD
        End If

        Call MPI_ALLREDUCE(sendbuf, recvbuf, icount, MPI_DOUBLE_PRECISION, MPI_SUM, comm, MPI_err)

    End Subroutine DALLSUM1D

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  DSUM3D
    !
    ! DESCRIPTION:  Performs MPI REDUCE (SUM) of 3-D array
    !                 sendbuf across grp.  Results are stored in rank ddest.
    !
    ! INPUTS:
    !            sendbuf - 3-D double-precision array to SUM across grp
    !                grp - MPI communicator across which the REDUCE is conducted
    !                       (optional; default = MPI_COMM_WORLD)
    !              ddest - MPI rank within grp that stores the results of the REDUCE operation
    ! OUTPUTS:
    !
    !           recvbuf - 3-D double-precision array that stores the result
    !                       of the REDUCE on rank ddest
    !/////////////////////////////////////////////////////////////////////
    Subroutine DSUM3D(sendbuf, recvbuf, grp, ddest)
        Real*8, Intent(In)  :: sendbuf(:,:,:)
        Real*8, Intent(Out) :: recvbuf(:,:,:)
        Type(communicator), Optional  :: grp
        Integer, Intent(In), Optional :: ddest
        Integer :: icount,  comm, MPI_err, dest

        If (present(ddest)) then
            dest = ddest
        Else
            dest = 0
        Endif

        icount = size(sendbuf)

        If (Present(grp)) Then
            comm = grp%comm
        Else
            comm = MPI_COMM_WORLD
        End If


        Call MPI_REDUCE(sendbuf, recvbuf, icount, MPI_DOUBLE_PRECISION, MPI_SUM, dest, comm, MPI_err)

    End Subroutine DSUM3D


    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  DSUM2D
    !
    ! DESCRIPTION:  Performs MPI REDUCE (SUM) of 2-D array
    !                 sendbuf across grp.  Results are stored in rank ddest.
    !
    ! INPUTS:
    !            sendbuf - 2-D double-precision array to SUM across grp
    !                grp - MPI communicator across which the REDUCE is conducted
    !                       (optional; default = MPI_COMM_WORLD)
    !              ddest - MPI rank within grp that stores the results of the REDUCE operation
    ! OUTPUTS:
    !
    !           recvbuf - 2-D double-precision array that stores the result
    !                       of the REDUCE on rank ddest
    !/////////////////////////////////////////////////////////////////////
    Subroutine DSUM2D(sendbuf, recvbuf, grp, ddest)
        Real*8, Intent(In)  :: sendbuf(:,:)
        Real*8, Intent(Out) :: recvbuf(:,:)
        Type(communicator), Optional  :: grp
        Integer, Intent(In), Optional :: ddest
        Integer :: icount,  comm, MPI_err, dest

        If (present(ddest)) then
            dest = ddest
        Else
            dest = 0
        Endif

        icount = size(sendbuf)

        If (Present(grp)) Then
            comm = grp%comm
        Else
            comm = MPI_COMM_WORLD
        End If


        Call MPI_REDUCE(sendbuf, recvbuf, icount, MPI_DOUBLE_PRECISION, MPI_SUM, dest, comm, MPI_err)

    End Subroutine DSUM2D

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  DALLSUM2D
    !
    ! DESCRIPTION:  Performs MPI AllReduce (SUM) of 2-D array
    !                 sendbuf across grp
    !
    ! INPUTS:
    !            sendbuf - 2-D double-precision array to SUM across grp
    !                grp - MPI communicator across which the ALLREDUCE is conducted
    !                    (optional; default = MPI_COMM_WORLD)
    ! OUTPUTS:
    !
    !           recvbuf - 2-D double-precision array that stores the result
    !                       of the ALLREDUCE
    !/////////////////////////////////////////////////////////////////////
    Subroutine DALLSUM2D(sendbuf, recvbuf, grp)
        Real*8, Intent(IN)  :: sendbuf(1:,1:)
        Real*8, Intent(Out) :: recvbuf(1:,1:)
        Type(communicator), Optional :: grp
        Integer :: icount,  comm, MPI_err


        icount = size(sendbuf)

        If (Present(grp)) Then
            comm = grp%comm
        Else
            comm = MPI_COMM_WORLD
        End If


        Call MPI_ALLREDUCE(sendbuf, recvbuf, icount, MPI_DOUBLE_PRECISION, MPI_SUM, comm, MPI_err)

    End Subroutine DALLSUM2D

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  BCAST2D
    !
    ! DESCRIPTION:  Performs MPI Broadcast of 2-D dimensional array buff across grp
    !              from source rank broot.
    !
    ! INPUTS:
    !            buff - 2-D double-precision array to be broadcast
    !             grp - MPI communicator across which the broadcast is conducted
    !                    (optional; default = MPI_COMM_WORLD)
    !           broot - MPI rank within grp that serves as the source of
    !                    the broadcast (optional; default = 0)
    !/////////////////////////////////////////////////////////////////////
    Subroutine BCAST2D(buff, grp, broot)
        Real*8, INTENT(INOUT) :: buff(:,:)
        Type(communicator), INTENT(IN), Optional :: grp
        Integer, Intent(In), Optional :: broot
        Integer :: icount,  comm, MPI_err, root

        If (present(broot)) then
            root = broot
        Else
            root = 0
        Endif

        icount = size(buff)

        If (Present(grp)) Then
            comm = grp%comm
        Else
            comm = MPI_COMM_WORLD
        End If

        Call MPI_BCAST(buff, icount, MPI_DOUBLE_PRECISION,  root, comm, MPI_err)

    End Subroutine BCAST2D

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  BCAST1D
    !
    ! DESCRIPTION:  Performs MPI Broadcast of 1-D dimensional array buff across grp
    !              from source rank broot.
    !
    ! INPUTS:
    !            buff -  4-byte-integer array to be broadcast
    !             grp - MPI communicator across which the broadcast is conducted
    !                    (optional; default = MPI_COMM_WORLD)
    !           broot - MPI rank within grp that serves as the source of
    !                    the broadcast (optional; default = 0)
    !/////////////////////////////////////////////////////////////////////
    Subroutine BCAST1D(buff, grp, broot)
        Integer*4, INTENT(INOUT) :: buff(:)
        Type(communicator), INTENT(IN), Optional :: grp
        Integer, Intent(In), Optional :: broot
        Integer :: icount,  comm, MPI_err, root

        If (present(broot)) then
            root = broot
        Else
            root = 0
        Endif

        icount = size(buff)

        If (Present(grp)) Then
            comm = grp%comm
        Else
            comm = MPI_COMM_WORLD
        End If

        Call MPI_BCAST(buff, icount, MPI_INTEGER,  root, comm, MPI_err)

    End Subroutine BCAST1D

End Module General_MPI
