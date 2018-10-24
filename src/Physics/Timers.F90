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

Module Timers
    Use Timing
    Use Parallel_Framework
    Use SendReceive
    Use Controls, Only : my_path
    Implicit None
    Integer, Parameter :: loop_time = 1, legendre_time = 2, fft_time = 3, solve_time = 4
    Integer, Parameter :: rtranspose_time = 5, ctranspose_time = 6
    Integer, Parameter :: rlmb_time = 7, rlma_time = 8, pspace_time = 9,psolve_time = 10
    Integer, Parameter :: dphi_time = 11, nl_time = 12, sdiv_time = 13, ts_time = 14
    Integer, Parameter :: ar_time = 15, seteq_time = 16, init_time = 17, cread_time = 18, cwrite_time = 19
    Integer, Parameter :: walltime = 20  ! This contains local elapsed time since mpi was initialized

    Integer, Parameter :: ntimers = 21
    Type(Timer), Allocatable :: StopWatch(:)
    Real*8 :: timer_ticklength !Length of 1 tick in seconds
Contains
    Subroutine Initialize_Timers()
        Integer :: i
        Allocate(StopWatch(1:ntimers))
        Do i = 1, ntimers
            Call StopWatch(i)%init()
        Enddo
        call get_ticklength(timer_ticklength)
    End Subroutine Initialize_Timers

    Subroutine Finalize_Timing(nr_in,lmax_in,niter)
        Implicit None
        Real*8, Allocatable :: mytimes(:), rowtimes(:,:), alltimes(:,:), buff(:), buff2(:,:)
        Integer :: row_rank, column_rank,colnp,rownp, i,j,np,p,offset
        Integer, Allocatable :: row_ranks(:), column_ranks(:)
        Integer :: timing_tag = 432
        Integer, Intent(In) :: nr_in, lmax_in,niter
        Character*120 :: timing_file
        Character*4   :: nr_string, lmax_string
        Character*4   :: row_string,col_string
        ! Gather times into a single array
        Allocate(mytimes(1:ntimers))


        Do i = 1, ntimers
            mytimes(i) = StopWatch(i)%elapsed
        Enddo

        ! Everyone sends down their row to rank zero
        ! Row rank zeros send to column rank zero
        row_rank = pfi%rcomm%rank
        column_rank = pfi%ccomm%rank
        rownp = pfi%rcomm%np
        colnp = pfi%ccomm%np

        If (row_rank .eq. 0) Then
            Allocate(rowtimes(1:ntimers,1:rownp))
            Allocate(buff(1:ntimers))
            rowtimes(:,1) = mytimes(:)
            Do p = 1, rownp-1
                Call receive(buff, source= p,tag=timing_tag,grp = pfi%rcomm)
                rowtimes(:,p+1) = buff(:)
            Enddo
            DeAllocate(buff)
            If (column_rank .eq. 0) Then
                Allocate(buff2(1:ntimers,1:rownp))
                Allocate(alltimes(1:ntimers, 1:rownp*colnp))
                alltimes(1:ntimers,1:rownp) = rowtimes(1:ntimers,1:rownp)
                offset = 1+rownp
                Do p = 1, colnp -1
                    Call receive(buff2, source= p,tag=timing_tag,grp = pfi%ccomm)
                    alltimes(1:ntimers,offset:offset+rownp-1) = buff2(1:ntimers,1:rownp)
                    offset = offset+rownp
                Enddo
                DeAllocate(buff2)
            Else
                Call send(rowtimes, dest = 0,tag=timing_tag, grp=pfi%ccomm)
            Endif
        Else
             Call send(mytimes, dest = 0,tag=timing_tag, grp=pfi%rcomm)
        Endif
        DeAllocate(mytimes)


        If ( (column_rank .eq. 0) .and. (row_rank .eq. 0) ) Then
            np = colnp*rownp
            Allocate(column_ranks(1:np))
            Allocate(row_ranks(1:np))
            offset = 1
            Do i = 0, colnp -1
                Do j = 0, rownp-1
                    column_ranks(offset) = i
                    row_ranks(offset) = j
                    offset = offset+1
                Enddo
            Enddo

            write(nr_string,'(i4.4)') nr_in
            write(lmax_string,'(i4.4)') lmax_in
            write(row_string,'(i4.4)') rownp
            write(col_string,'(i4.4)') colnp
            timing_file = 'Timings/lmax'//TRIM(lmax_string)//'_nr'//TRIM(nr_string)//'_ncol'//TRIM(col_string)
            timing_file = Trim(my_path)//TRIM(timing_file)//'_nrow'//TRIM(row_string)
         Open(unit=15,file=timing_file,status='replace', ACCESS="STREAM")
         Write(15)colnp
            Write(15)rownp
            Write(15)ntimers
            Write(15)nr_in
            Write(15)lmax_in
            Write(15)niter
         Write(15)(column_ranks(i),i=1,np)
         Write(15)(row_ranks(i),i=1,np)
            Write(15)((alltimes(i,j),i=1,ntimers),j=1,np)
         Close(15)


            DeAllocate(column_ranks,row_ranks,alltimes)
        Endif

    End Subroutine Finalize_Timing

End Module Timers
