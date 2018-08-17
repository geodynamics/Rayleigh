!////////////////////////////////////////////////////////////////////////////
!   MODULE:  ALL_TO_ALL
!
!   DESCRIPTION:  Provides the Standard_Transpose interface for calling MPI_ALLTOALL
!                   and MPI_ALLTOALLV
!
!   INTERFACES:
!       Standard_Transpose
!
!   MEMBER SUBROUTINES:
!       Z_Transpose_v_1D
!       D_Transpose_v_1D
!       D_Transpose_choose_1D
!////////////////////////////////////////////////////////////////////////////

Module All_To_All
    Use MPI_BASE

    Implicit None

    Private

    Public :: Standard_Transpose


    !////////////////////////////////////////////////////////////////////////
    !
    ! Interface:   Standard_Transpose
    !
    ! Description:  Simple interface for MPI AlltoAllv and AlltoALL wrappers
    !
    !////////////////////////////////////////////////////////////////////////
    Interface Standard_Transpose
        Module Procedure Z_Transpose_v_1D, D_Transpose_v_1D, D_Transpose_choose_1D
    End Interface

Contains

    !////////////////////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  Z_Tranpose_v_1D
    !
    ! DESCRIPTION:  Performs an AlltoALLv for a 1-D, double-complex array across grp
    !
    ! INPUTS:
    !            send_buf - Double-Complex 1-D array that serves as the source for
    !                         the alltoallv operation
    !                 grp - MPI communicator across which the REDUCE is conducted 
    !                       (optional; default = MPI_COMM_WORLD)
    !            send_displ - array of send displacements for mpi_alltoallv
    !            recv_displ - array of receive displacements for mpi_alltoallv
    !            send_count - array of send counts for mpi_alltoallv
    !            recv_count - array of receive counts for mpi_alltoallv
    !
    ! OUTPUTS:
    !
    !            recvbuf - Double-Complex 1-D array that serves as the receptor for
    !                        the alltoallv operation
    !////////////////////////////////////////////////////////////////////////////////////
    Subroutine Z_Transpose_v_1D(send_buf, recv_buf, send_count, send_displ, recv_count, recv_displ, grp)
        Complex*16, Intent(In)  :: send_buf(:)
        Complex*16, Intent(Out) :: recv_buf(:)

        Integer, Intent(in) :: send_count(:), send_displ(:)
        Integer, Intent(in) :: recv_count(:), recv_displ(:)
        Type(communicator), Intent(in) :: grp
        Integer :: MPI_err

        call MPI_ALLTOALLv(send_buf, send_count, send_displ, MPI_DOUBLE_COMPLEX, recv_buf, &
            & recv_count, recv_displ, MPI_DOUBLE_COMPLEX, grp%comm, MPI_err)
    End Subroutine Z_Transpose_v_1D

    !////////////////////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  D_Tranpose_v_1D
    !
    ! DESCRIPTION:  Performs an AlltoALLv for a 1-D, double-precision array across grp
    !
    ! INPUTS:
    !            send_buf - Double-precision 1-D array that serves as the source for
    !                         the alltoallv operation
    !                 grp - MPI communicator across which the REDUCE is conducted 
    !                       (optional; default = MPI_COMM_WORLD)
    !            send_displ - array of send displacements for mpi_alltoallv
    !            recv_displ - array of receive displacements for mpi_alltoallv
    !            send_count - array of send counts for mpi_alltoallv
    !            recv_count - array of receive counts for mpi_alltoallv
    !
    ! OUTPUTS:
    !
    !            recvbuf - Double-precision 1-D array that serves as the receptor for
    !                        the alltoallv operation
    !////////////////////////////////////////////////////////////////////////////////////

    Subroutine D_Transpose_v_1D(send_buf, recv_buf, send_count, send_displ, recv_count, recv_displ, grp)
        Real*8, Intent(In)  :: send_buf(:)
        Real*8, Intent(Out) :: recv_buf(:)

        Integer, Intent(in) :: send_count(:), send_displ(:)
        Integer, Intent(in) :: recv_count(:), recv_displ(:)
        Type(communicator), Intent(in) :: grp
        Integer :: MPI_err

        call MPI_ALLTOALLv(send_buf, send_count, send_displ, MPI_DOUBLE_PRECISION, recv_buf, &
            & recv_count, recv_displ, MPI_DOUBLE_PRECISION, grp%comm, MPI_err)
    End Subroutine D_Transpose_v_1D

    !////////////////////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  D_Tranpose_choose_1D
    !
    ! DESCRIPTION:  Performs an AlltoALLv OR AlltoALL for a 1-D, 
    !                 double-precision array across grp
    !
    ! INPUTS:
    !            send_buf - Double-precision 1-D array that serves as the source for
    !                         the alltoallv operation
    !                 grp - MPI communicator across which the REDUCE is conducted 
    !                       (optional; default = MPI_COMM_WORLD)
    !            send_displ - array of send displacements for mpi_alltoallv
    !            recv_displ - array of receive displacements for mpi_alltoallv
    !            send_count - array of send counts for mpi_alltoallv
    !            recv_count - array of receive counts for mpi_alltoallv
    !            normal - Boolean flag that controls use of AlltoAll (TRUE) or AlltoAllv (False)
    !
    ! OUTPUTS:
    !
    !            recvbuf - Double-precision 1-D array that serves as the receptor for
    !                        the alltoallv operation
    !////////////////////////////////////////////////////////////////////////////////////
    Subroutine D_Transpose_choose_1D(send_buf, recv_buf, send_count, send_displ, recv_count, recv_displ, grp, normal)
        Real*8, Intent(In)  :: send_buf(:)
        Real*8, Intent(Out) :: recv_buf(:)
        Integer, Intent(in) :: send_count(:), send_displ(:)
        Integer, Intent(in) :: recv_count(:), recv_displ(:)
        Type(communicator), Intent(in) :: grp
        Integer :: MPI_err
        Logical, Intent(in) :: normal
        If (normal) Then
            call MPI_ALLTOALL(send_buf, send_count(1), MPI_DOUBLE_PRECISION, recv_buf, &
                & recv_count(1), MPI_DOUBLE_PRECISION, grp%comm, MPI_err)
        Else
            call MPI_ALLTOALLv(send_buf, send_count, send_displ, MPI_DOUBLE_PRECISION, recv_buf, &
                & recv_count, recv_displ, MPI_DOUBLE_PRECISION, grp%comm, MPI_err)
        Endif
    End Subroutine D_Transpose_choose_1D

End Module All_To_All
