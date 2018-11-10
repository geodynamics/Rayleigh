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
!   MODULE:  ISENDRECEIVE
!
!   DESCRIPTION:  Provides interface to and wrappers for calling
!                   MPI_ISEND/MPI_IRECEIVE
!
!   MEMBER SUBROUTINES:
!       IWait, IWaitALL
!       D_ISend_1D, D_ISend_2D, D_ISend_3D, D_ISend_4D, D_ISend_5D
!       Z_ISend_1D, Z_ISend_2D, Z_ISend_3D
!       D_IReceive_1D, D_IReceive_2D, D_IReceive_3D, D_IReceive_4D, D_IReceive_5D
!       Z_IReceive_1D, Z_IReceive_2D, Z_IReceive_3D
!
!////////////////////////////////////////////////////////////////////////////

Module ISendReceive
    Use RA_MPI_BASE

    !////////////////////////////////////////////////////////////////////////
    !
    ! Interface:   ISend
    !
    ! Description:  Simple interface for MPI_ISend
    !
    !////////////////////////////////////////////////////////////////////////
    Interface ISend
        Module Procedure D_ISend_1D, D_ISend_2D, D_ISend_3D, D_ISend_4D
        Module Procedure D_ISend_5D
        Module Procedure Z_ISend_1D, Z_ISend_2D, Z_ISend_3D
    End Interface

    !////////////////////////////////////////////////////////////////////////
    !
    ! Interface:   IReceive
    !
    ! Description:  Simple interface for MPI_IReceive
    !
    !////////////////////////////////////////////////////////////////////////
     Interface IReceive
        Module Procedure D_IReceive_1D, D_IReceive_2D, D_IReceive_3D, D_IReceive_4D
        Module Procedure D_IReceive_5D
        Module Procedure Z_IReceive_1D, Z_IReceive_2D, Z_IReceive_3D
    End Interface

  Integer :: mpi_err

Contains

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  IWait
    !
    ! DESCRIPTION:  Wrapper to MPI_WAIT - blocks program until the process
    !                 associated with irq has completed.
    !
    ! INPUTS:
    !            irq - Integer variable whose value is associated with
    !                    a non-blocking MPI process
    !
    !/////////////////////////////////////////////////////////////////////
    Subroutine IWait(irq)
        Implicit None
        Integer :: irq, status(MPI_STATUS_SIZE), mpi_err
        Call MPI_WAIT(irq,status,mpi_err)
    End Subroutine IWait

    !/////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  IWaitAll
    !
    ! DESCRIPTION:  Wrapper to MPI_WAIT - blocks program until the processes
    !                 associated with irq have completed.
    !
    ! INPUTS:
    !              n - number of requests to wait on/size of array irq
    !            irq - Integer array whose values are associated with
    !                    multiple non-blocking MPI processes
    !
    !/////////////////////////////////////////////////////////////////////
    Subroutine IWaitAll(n,irq)
        Integer :: irq(:)
        Integer, Intent(In) :: n
        Integer :: mpi_err
        Integer, Allocatable :: istat(:,:)
        Allocate(istat(MPI_STATUS_SIZE,1:n))
        Call MPI_WAITALL(n,irq,istat,mpi_err)
        DeAllocate(istat)
    End Subroutine IWaitAll


    !//////////////////////////////////////////////////////////////////////
    !  NOTE:        The send and receive routines follow the same calling pattern
    !           with only minor differences related to the input/output
    !           parameters.  We provide a general description of their
    !           characteristics below
    !//////////////////////////////////////////////////////////////////////
    !  SUBROUTINE:  (X)_ISend_(Y)D
    !               X:  Denotes the variable type for the send and
    !                     receive buffers.  D indicates double precision.
    !                     Z indicates double-complex precision.
    !               Y:  The dimensionality of the send and receive buffers (1,2,3, etc.)
    !
    !  DESCRIPTION:  Performs a non-blocking send of array 'x' to MPI rank
    !                   'dest' in MPI communictory group 'grp'
    !  INPUTS:
    !             x          - Array of dimension Y and of datatype X to be sent to dest
    !             indstart   - Integer array of dimension Y indicating where within x
    !                          the send should initiate from.
    !                          (optional; default value is 1,1,.. - (send initiates at x(1,1,...)
    !             n_elements - Number of elements from x to broadcast
    !                          (optional; default value is size(x))
    !             dest       - The MPI rank (within communicator grp) to which x is sent
    !                          (optional; defaults to 0 if unspecified)
    !             tag        - The MPI tag for the send.
    !                          (optional;  defaults to dest if unspecified)
    !             grp        - The MPI communicator associated with this send
    !                          (optional; defaults to MPI_COMM_WORLD if unspecified)
    !
    !  OUTPUTS:   irq - (Integer) MPI request handle for use with MPI_WAIT
    !////////////////////////////////////////////////////////////////////////////////////////

    !//////////////////////////////////////////////////////////////////////
    !  SUBROUTINES:  (X)_IReceive_(Y)D
    !               X:  Denotes the variable type for the send and
    !                     receive buffers.  D indicates double precision.
    !                     Z indicates double-complex precision.
    !               Y:  The dimensionality of the send and receive buffers (1,2,3, etc.)
    !
    !  DESCRIPTION:  Performs a non-blocking receive of array 'x' from MPI rank
    !                   'source' in MPI communictory group 'grp'
    !  INPUTS:
    !             indstart   - Integer array of dimension Y indicating where within x
    !                          the receive should initiate from.
    !                          (optional; default value is 1,1,.. - (send initiates at x(1,1,...)
    !             n_elements - Number of elements within x to receive
    !                          (optional; default value is size(x))
    !             source     - The MPI rank (within communicator grp) from which x is received
    !                          (optional; defaults to 0 if unspecified)
    !             tag        - The MPI tag for the receive.
    !                          (optional;  defaults to source if unspecified)
    !             grp        - The MPI communicator associated with this send
    !                          (optional; defaults to MPI_COMM_WORLD if unspecified)
    !
    !  OUTPUTS:   irq - (Integer) MPI request handle for use with MPI_WAIT
    !               x -  Array of dimension Y and of datatype X to be received from source

    !Subroutine D_IReceive_4D(x, irq,n_elements, source, tag, grp,indstart)


    Subroutine D_ISend_1D(x, irq,n_elements, istart,dest, tag, grp)
        Real(8), Intent(In)  :: x(:)
        Integer, Intent(In), Optional :: dest, n_elements, tag,istart
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, ione

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

        If (Present(istart)) Then
            ione = istart
        Else
            ione = 1
        End if

        Call mpi_isend(x(ione), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)

    End Subroutine D_ISend_1D

    Subroutine D_ISend_4D(x, irq,n_elements, dest, tag, grp, indstart)
        Real*8, Intent(in)  :: x(1:,1:,1:,1:)
        Integer, Intent(In), Optional :: dest, n_elements, tag,indstart(1:4)
        Type(communicator), Intent(In), optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, istart, kstart, jstart,lstart

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
        Call mpi_isend(x(istart,jstart,kstart,lstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)

    End Subroutine D_ISend_4D

    Subroutine D_IReceive_4D(x, irq,n_elements, source, tag, grp,indstart)
        Real*8, Intent(Out)  :: x(1:,1:,1:,1:)
        Integer, Intent(In), Optional :: source, n_elements, tag,indstart(1:4)
        Type(communicator), Intent(In), Optional :: grp
        Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE)
        Integer, Intent(Out) :: irq
        Integer :: istart,jstart,kstart,lstart

        If (Present(n_elements)) Then
           n = n_elements
        Else
           n = Size(x)
        End If

        If (Present(source)) Then
           p = source
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

        Call mpi_irecv(x(istart,jstart,kstart,lstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

    End Subroutine D_IReceive_4D

    Subroutine D_ISend_5D(x, irq,n_elements, dest, tag, grp, indstart)
        Real*8, Intent(In)  :: x(1:,1:,1:,1:,1:)
        Integer, Intent(In), Optional :: dest, n_elements, tag,indstart(1:5)
        Type(communicator), optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2
        Integer :: istart, kstart, jstart,lstart, mstart

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
        Call mpi_isend(x(istart,jstart,kstart,lstart,mstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)

    End Subroutine D_ISend_5D

    Subroutine D_IReceive_5D(x, irq,n_elements, source, tag, grp,indstart)
        Real*8, Intent(Out)  :: x(1:,1:,1:,1:,1:)
        Integer, Intent(In), Optional :: source, n_elements, tag,indstart(1:5)
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE)
        Integer :: istart,jstart,kstart,lstart,mstart

        If (Present(n_elements)) Then
           n = n_elements
        Else
           n = Size(x)
        End If

        If (Present(source)) Then
           p = source
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

        Call mpi_irecv(x(istart,jstart,kstart,lstart,mstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

    End Subroutine D_IReceive_5D

    Subroutine Z_ISend_1D(x, irq,n_elements, istart,dest, tag, grp)
        Complex*16, Intent(In)  :: x(:)
        Integer, Intent(In), Optional :: dest, n_elements, tag,istart
        Type(communicator),Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, ione

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

        If (Present(istart)) Then
            ione = istart
        Else
            ione = 1
        End if

        Call mpi_isend(x(ione), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq,mpi_err)

    End Subroutine Z_ISend_1D

    Subroutine D_ISend_2D(x, irq,n_elements, dest, tag, grp,indstart)
        Real(8), Intent(In)  :: x(1:,1:)
        Type(communicator), Intent(In), optional :: grp
        Integer, Intent(In), Optional :: dest, n_elements, tag, indstart(2)
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2
        Integer :: ione, jone

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
            ione = indstart(1)
            jone = indstart(2)
        Else
            ione = 1
            jone = 1
        Endif
        Call mpi_isend(x(ione,jone), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)

    End Subroutine D_ISend_2D

    Subroutine D_ISend_3D(x, irq,n_elements, dest, tag, grp, indstart)
        Real*8, Intent(in)  :: x(1:,1:,1:)
        Integer, Intent(In), Optional :: dest, n_elements, tag,indstart(1:3)
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2
        Integer :: istart, kstart, jstart

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
        Call mpi_isend(x(istart,jstart,kstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq,mpi_err)

    End Subroutine D_ISend_3D

    Subroutine Z_ISend_2D(x, irq,n_elements, dest, tag, grp)
        Complex*16, Intent(In)  :: x(:,:)
        Integer, Intent(In), Optional :: dest, n_elements, tag
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
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

        Call mpi_isend(x(1,1), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq,mpi_err)

    End Subroutine Z_ISend_2D

    Subroutine Z_ISend_3D(x, irq,n_elements, dest, tag, grp, indstart)
        Complex*16, Intent(In)  :: x(:,:,:)
        Integer, Intent(In), Optional :: dest, n_elements, tag,indstart(1:3)
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, istart, kstart, jstart

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
        Call mpi_isend(x(istart,jstart,kstart), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq,mpi_err)

      End Subroutine Z_ISend_3D

    Subroutine D_IReceive_1D(x, irq, n_elements, istart, source, tag, grp)
        Real(8), Intent(Out)  :: x(:)
        Integer, Intent(In), Optional :: source, n_elements, tag, istart
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE), ione

        If (Present(n_elements)) Then
            n = n_elements
        Else
            n = Size(x)
        End If

        If (Present(source)) Then
            p = source
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


        If (Present(istart)) Then
            ione = istart
        Else
            ione = 1
        Endif

        Call mpi_irecv(x(ione), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

    End Subroutine D_IReceive_1D

    Subroutine Z_IReceive_1D(x, irq, n_elements, istart, source, tag, grp)
        Complex*16, Intent(Out)  :: x(:)
        Integer, Intent(In), Optional :: source, n_elements, tag, istart
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE), ione

        If (Present(n_elements)) Then
           n = n_elements
        Else
           n = Size(x)
        End If

        If (Present(source)) Then
           p = source
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


        If (Present(istart)) Then
            ione = istart
        Else
            ione = 1
        Endif

        Call mpi_irecv(x(ione), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq, mpi_err)

    End Subroutine Z_IReceive_1D

    Subroutine D_IReceive_2D(x, irq, n_elements, source, tag, grp, indstart)
        Real(8), Intent(Out):: x(1:,1:)
        Integer, Intent(In), Optional :: source, n_elements, tag, indstart(1:2)
        Type(communicator),  Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE), ione, jone

        If (Present(n_elements)) Then
            n = n_elements
        Else
            n = Size(x)
        End If

        If (Present(source)) Then
            p = source
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

        If (present(indstart)) Then
            ione = indstart(1)
            jone = indstart(2)
        Else
            ione = 1
            jone = 1
        Endif


        Call mpi_irecv(x(ione,jone), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

    End Subroutine D_IReceive_2D

    Subroutine D_IReceive_3D(x, irq,n_elements, source, tag, grp,indstart)
        Real*8, Intent(Out)  :: x(1:,1:,1:)
        Integer, Intent(In), Optional :: source, n_elements, tag,indstart(1:3)
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE)
        Integer :: istart,jstart,kstart

        If (Present(n_elements)) Then
           n = n_elements
        Else
           n = Size(x)
        End If

        If (Present(source)) Then
           p = source
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

        Call mpi_irecv(x(istart,jstart,kstart), n, MPI_DOUBLE_PRECISION, p, tag2, comm2, irq, mpi_err)

    End Subroutine D_IReceive_3D

    Subroutine Z_IReceive_2D(x, irq, n_elements, source, tag, grp, indstart)
        Complex*16, Intent(Out)  :: x(:,:)
        Integer, Intent(In), Optional :: source, n_elements, tag, indstart(1:2)
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE), istart, jstart

        If (Present(n_elements)) Then
           n = n_elements
        Else
           n = Size(x)
        End If

        If (Present(source)) Then
           p = source
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

        Call mpi_irecv(x(istart,jstart), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq, mpi_err)

    End Subroutine Z_IReceive_2D

    Subroutine Z_IReceive_3D(x, irq,n_elements, source, tag, grp,indstart)
        Complex*16, Intent(Out)  :: x(:,:,:)
        Integer, Intent(In), Optional :: source, n_elements, tag,indstart(1:3)
        Type(communicator), Intent(In), Optional :: grp
        Integer, Intent(Out) :: irq
        Integer :: p, n, comm2, tag2, status(MPI_STATUS_SIZE)
        Integer :: istart, jstart, kstart

        If (Present(n_elements)) Then
           n = n_elements
        Else
           n = Size(x)
        End If

        If (Present(source)) Then
           p = source
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

        Call mpi_irecv(x(istart,jstart,kstart), n, MPI_DOUBLE_COMPLEX, p, tag2, comm2, irq, mpi_err)

    End Subroutine Z_IReceive_3D

End Module ISendReceive
