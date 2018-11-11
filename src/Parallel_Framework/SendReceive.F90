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

Module SendReceive
    Use RA_MPI_BASE
    Implicit None
    Private
    Integer :: mpi_err
    Public :: send, receive
    Interface Send
        Module Procedure D_Send_5D, D_Send_4D, D_Send_3D, D_Send_2D, D_Send_1D
    End Interface

    Interface Receive
        Module Procedure D_Receive_4D, D_Receive_3D, D_Receive_2D, D_Receive_1D
        Module Procedure D_Receive_5D
    End Interface

Contains

    Subroutine D_Send_5D(x, n_elements, dest, tag, grp, indstart)
    Real*8, Intent(in)  :: x(1:,1:,1:,1:,1:)

    Integer, Optional :: dest, n_elements, tag,indstart(1:5)
     Integer :: istart, kstart, jstart,lstart, mstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
     If (Present(indstart)) Then
        istart = indstart(1)
        jstart = indstart(2)
        kstart = indstart(3)
        lstart = indstart(4)
        mstart = indstart(5)
    Else
        istart = 1
        jstart = 1
        kstart = 1
        lstart = 1
        mstart = 1
    Endif
    Call mpi_send(x(istart,jstart,kstart,lstart,mstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2,  mpi_err)
    !write(6,*)'zs ', p
    End Subroutine D_Send_5D



    Subroutine D_Send_4D(x, n_elements, dest, tag, grp, indstart)
    Real*8, Intent(in)  :: x(:,:,:,:)

    Integer, Optional :: dest, n_elements, tag,indstart(1:4)
     Integer :: istart, kstart, jstart,lstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
     If (Present(indstart)) Then
        istart = indstart(1)
        jstart = indstart(2)
        kstart = indstart(3)
        lstart = indstart(4)
    Else
        istart = 1
        jstart = 1
        kstart = 1
        lstart = 1
    Endif
    Call mpi_send(x(istart,jstart,kstart,lstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2,  mpi_err)
    !write(6,*)'zs ', p
    End Subroutine D_Send_4D

    Subroutine D_Send_3D(x, n_elements, dest, tag, grp, indstart)
    Real*8, Intent(in)  :: x(1:,1:,1:)

    Integer, Optional :: dest, n_elements, tag,indstart(1:3)
     Integer ::  istart, kstart, jstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
     If (Present(indstart)) Then
        istart = indstart(1)
        jstart = indstart(2)
        kstart = indstart(3)
    Else
        istart = 1
        jstart = 1
        kstart = 1
    Endif
    Call mpi_send(x(istart,jstart,kstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2,  mpi_err)
    !write(6,*)'zs ', p
    End Subroutine D_Send_3D

    Subroutine D_Send_2D(x, n_elements, dest, tag, grp, indstart)
    Real*8, Intent(in)  :: x(:,:)

    Integer, Optional :: dest, n_elements, tag,indstart(1:2)
     Integer :: istart, jstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
     If (Present(indstart)) Then
        istart = indstart(1)
        jstart = indstart(2)
    Else
        istart = 1
        jstart = 1
    Endif
    Call mpi_send(x(istart,jstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2,  mpi_err)

    End Subroutine D_Send_2D

    Subroutine D_Send_1D(x, n_elements, dest, tag, grp, indstart)
    Real*8, Intent(in)  :: x(:)

    Integer, Optional :: dest, n_elements, tag,indstart(1)
     Integer :: istart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(dest)) Then
       p = dest
    Else
       p = 0
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = p
    End If
     If (Present(indstart)) Then
        istart = indstart(1)
    Else
        istart = 1
    Endif
    Call mpi_send(x(istart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2,  mpi_err)

    End Subroutine D_Send_1D



    Subroutine D_Receive_5D(x, n_elements, source, tag, grp,indstart)
        Real*8, Intent(out)  :: x(1:,1:,1:,1:,1:)

    Integer, Optional :: source, n_elements, tag,indstart(1:5)
     Integer :: istart,jstart,kstart,lstart, mstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, mstatus(MPI_STATUS_SIZE)

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

    If (Present(indstart)) Then
        istart = indstart(1)
        jstart = indstart(2)
        kstart = indstart(3)
        lstart = indstart(4)
        mstart = indstart(5)
    Else
        istart = 1
        jstart = 1
        kstart = 1
        lstart = 1
        mstart = 1
    Endif

    Call mpi_recv(x(istart,jstart,kstart,lstart,mstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, mstatus, mpi_err)


    End Subroutine D_Receive_5D



    Subroutine D_Receive_4D(x, n_elements, source, tag, grp,indstart)
        Real*8, Intent(out)  :: x(:,:,:,:)

    Integer, Optional :: source, n_elements, tag,indstart(1:4)
     Integer :: istart,jstart,kstart,lstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, mstatus(MPI_STATUS_SIZE)

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

    If (Present(indstart)) Then
        istart = indstart(1)
        jstart = indstart(2)
        kstart = indstart(3)
        lstart = indstart(4)
    Else
        istart = 1
        jstart = 1
        kstart = 1
        lstart = 1
    Endif

    Call mpi_recv(x(istart,jstart,kstart,lstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, mstatus, mpi_err)


    End Subroutine D_Receive_4D

    Subroutine D_Receive_3D(x, n_elements, source, tag, grp,indstart)
        Real*8, Intent(out)  :: x(1:,1:,1:)

    Integer, Optional :: source, n_elements, tag,indstart(1:3)
     Integer :: istart,jstart,kstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, mstatus(MPI_STATUS_SIZE)

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

    If (Present(indstart)) Then
        istart = indstart(1)
        jstart = indstart(2)
        kstart = indstart(3)
    Else
        istart = 1
        jstart = 1
        kstart = 1
    Endif

    Call mpi_recv(x(istart,jstart,kstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, mstatus, mpi_err)


    End Subroutine D_Receive_3D

    Subroutine D_Receive_2D(x, n_elements, source, tag, grp,indstart)
        Real*8, Intent(out)  :: x(:,:)

    Integer, Optional :: source, n_elements, tag,indstart(1:2)
     Integer :: istart,jstart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, mstatus(MPI_STATUS_SIZE)

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

    If (Present(indstart)) Then
        istart = indstart(1)
        jstart = indstart(2)
    Else
        istart = 1
        jstart = 1
    Endif

    Call mpi_recv(x(istart,jstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, mstatus, mpi_err)


    End Subroutine D_Receive_2D

    Subroutine D_Receive_1D(x, n_elements, source, tag, grp,indstart)
        Real*8, Intent(out)  :: x(:)

    Integer, Optional :: source, n_elements, tag,indstart(1)
     Integer :: istart
    Type(communicator), optional :: grp
    Integer :: p, n, comm2, tag2, mstatus(MPI_STATUS_SIZE)

    If (Present(n_elements)) Then
       n = n_elements
    Else
       n = Size(x)
    End If

    If (Present(source)) Then
       p = source
    Else
       p = MPI_ANY_SOURCE
    End If

    If (Present(grp)) Then
       comm2 = grp%comm

    Else
       comm2 = MPI_COMM_WORLD
    End If

    If (Present(tag)) Then
       tag2 = tag
    Else
       tag2 = MPI_ANY_TAG
    End If

    If (Present(indstart)) Then
        istart = indstart(1)
    Else
        istart = 1
    Endif

    Call mpi_recv(x(istart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, mstatus, mpi_err)


    End Subroutine D_Receive_1D

End Module SendReceive
