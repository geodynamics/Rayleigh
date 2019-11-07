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

Module Checkpointing
    Use Timers, Only : stopwatch
    Use ProblemSize
    Use Parallel_Framework
    Use Spherical_Buffer
    Use Linear_Solve, Only : get_all_rhs
    Use SendReceive
    Use ISendReceive
    Use Controls
    Use RA_MPI_BASE
    Use Chebyshev_Polynomials_Alt
    Use Physical_IO_Buffer

    Use BufferedOutput
    ! Simple Checkpointing Module
    ! Uses MPI-IO to split writing of files amongst rank zero processes from each row
    Implicit None
    Type(SphericalBuffer) :: chktmp, chktmp2
    Integer, private :: numfields = 4 ! 6 for MHD
    Integer,private,Allocatable :: mode_count(:)
    Integer,private :: nlm_total, checkpoint_tag = 425
    Integer, Allocatable, Private :: lmstart(:)
    Integer :: buffsize, buffsize2 ! These cannot be long (8 byte) - MPI-IO goes haywire.
    Integer*8, Private :: my_check_disp,full_disp        ! These need to be long ! for writing checkpoints
    Integer*8, Private :: my_check_disp2
    Integer*8, Private :: my_in_disp, full_in_disp        ! for reading checkpoints
    Character*3 :: wchar = 'W', pchar = 'P', tchar = 'T', zchar = 'Z', achar = 'A', cchar = 'C'
    Character*120 :: checkpoint_prefix ='nothing'
    Character*6 :: auto_fmt = '(i2.2)'  ! Format code for quicksaves
    Character*3 :: checkpoint_suffix(12)
    Integer :: checkpoint_iter = 0
    Real*8  :: checkpoint_dt, checkpoint_newdt
    Real*8  :: checkpoint_time
    !//////////////////////////
    ! New variables for new checkpointing style
    !Integer :: Noutputs_per_row, Nradii_per_output
    !Integer, Allocatable :: Nradii_at_rank(:), Rstart_at_rank(:)
    Logical :: I_Will_Output = .False.

    !///////////////////////////////////////////////////////////
    ! These variables are used for determining if it's time for a checkpoint
    Logical :: ItIsTimeForACheckpoint = .false.
    Logical :: ItIsTimeForAQuickSave = .false.
    Integer :: quicksave_num = -1
    Real*8  :: checkpoint_t0 = 0.0d0
    Real*8  :: checkpoint_elapsed = 0.0d0  ! Time elapsed since checkpoint_t0
    Real*8  :: quicksave_seconds = -1  ! Time between quick saves

    Type(Cheby_Transform_Interface) :: cheby_info

    Type(io_buffer_physical) :: checkpoint_buffer, checkpoint_inbuffer

Contains

    !//////////////////////////////////////////////////////////
    ! Modifications related to the new Memory-Friendly Checkpointing Style
    Subroutine Initialize_Checkpointing()
        Implicit None
        Integer :: nfs(6)
        Integer :: p, np, nl, m, mp, rextra

        checkpoint_t0 = stopwatch(walltime)%elapsed
        If (check_frequency .gt. 0) Then
            !this is for backwards compatibility
            !If specified, use check_frequency
            checkpoint_interval = check_frequency
        Endif

        If (num_quicksaves .gt. 100) Then
            !Maximum number of quick saves is 100 (hopefully far more than needed).
            num_quicksaves = 100
        Endif

        If (quicksave_minutes .gt. 0) Then
            quicksave_seconds = quicksave_minutes*60
            quicksave_interval = -1
        Endif

        if (magnetism) Then
            numfields = 6
            checkpoint_suffix(1:12) = &
                (/ 'W  ', 'P  ', 'T  ', 'Z  ', 'C  ', 'A  ', 'WAB', 'PAB', 'TAB', 'ZAB', 'CAB', 'AAB' /)
        Else
            checkpoint_suffix(1:8) = &
                (/ 'W  ', 'P  ', 'T  ', 'Z  ', 'WAB', 'PAB', 'TAB', 'ZAB' /)
        Endif
        nfs(:) = numfields*2
        Call chktmp%init(field_count = nfs, config = 'p1a')            ! This structure hangs around through the entire run

        Call checkpoint_buffer%Init(mpi_tag=checkpoint_tag,spectral=.true., cache_spectral = .true., &
                                    spec_comp = .true.)

 
        If (my_row_rank .lt. 1) Then
            I_Will_Output = .true.
        Endif

        np = pfi%rcomm%np
        Allocate(mode_count(0:np-1))
        mode_count(:)  = 0    ! This is how many total l-m combinations rank p of a row owns



        nlm_total = 0            ! This is the total number of l-m combinations
        Do p = 0, np -1
            Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
                m = m_values(mp)
                nl = l_max-m+1
                mode_count(p) = mode_count(p)+nl
                nlm_total = nlm_total+nl
            Enddo
        Enddo
        if (I_Will_Output) Then
            Allocate(lmstart(0:l_max))
            lmstart(0) = 1
            do m = 0, l_max-1
                nl = l_max-m+1
                lmstart(m+1) = lmstart(m)+nl
            enddo
            np = pfi%ccomm%np
            my_check_disp = 0
            Do p = 1, my_column_rank
                my_check_disp = my_check_disp+pfi%all_1p(p-1)%delta
            enddo
            my_check_disp = my_check_disp*nlm_total  ! for old Checkpointing style
            buffsize = nlm_total*my_r%delta  ! for old
            full_disp = N_r
            full_disp = full_disp*nlm_total
        Endif

    End Subroutine Initialize_Checkpointing

    Subroutine gen_data()
        Implicit None
        Integer :: mp, m, r, l, imi
        Character*120 :: testfile='test_data'
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            chktmp%s2a(mp)%data(:,:,:,1) = 0.0d0
            

            Do imi = 1, 2
            Do r = my_r%min, my_r%max
            Do l = m, l_max
                chktmp%s2a(mp)%data(l,r,imi,1) = (m+1)*(1.0d0+imi)*(r + l)
            Enddo
            Enddo
            Enddo

        Enddo
        Call checkpoint_buffer%cache_data_spectral(chktmp%s2a,1)
        Call checkpoint_buffer%write_data(filename=testfile)
    End Subroutine gen_data

    Subroutine Check_Data
        Implicit None
        Type(IO_Buffer_Physical) :: inbuffer(3)
        Character*120 :: files(3)
        Integer :: i, lmaxes(3), mp, m, r, l, imi
        Real*8 :: maxdiff, diff
        files(1) = 'test_data_32'
        files(2) = 'test_data_64'
        files(3) = 'test_data_128'

        lmaxes(1) = 31
        lmaxes(2) = 63
        lmaxes(3) = 127
        Write(6,*)'lmax is: ', l_max
        Do i = 1,3
            Call inbuffer(i)%Init(mpi_tag=checkpoint_tag,spectral=.true., cache_spectral = .true., &
                spec_comp = .true., lmax_in = lmaxes(i), mode = 2)
            Call inbuffer(i)%read_data(filename=files(i))

            Do mp = my_mp%min, my_mp%max
                chktmp%s2a(mp)%data(:,:,:,1) = 0.0d0
            Enddo

            Call inbuffer(i)%grab_data_spectral(chktmp%s2a,1)
            maxdiff = 0.0d0
            Do mp = my_mp%min, my_mp%max
                m = m_values(mp)

                Do imi = 1, 2
                Do r = my_r%min, my_r%max
                Do l = m, l_max
                    If (l .le. lmaxes(i)) then
                        diff  =  chktmp%s2a(mp)%data(l,r,imi,1)- (m+1)*(1.0d0+imi)*(r + l)
                    Else
                        diff  =  chktmp%s2a(mp)%data(l,r,imi,1) ! should be zero
                    Endif
                    diff = abs(diff)
                    if (diff .gt. maxdiff) maxdiff = diff
                Enddo
                Enddo
                Enddo

            Enddo
            Write(6,*)'i/maxdiff', i, maxdiff
        Enddo
    End Subroutine Check_Data

    Subroutine Write_Checkpoint(abterms,iteration,dt,new_dt,elapsed_time)
        Implicit None
        Real*8, Intent(In) :: abterms(:,:,:,:), dt, new_dt, elapsed_time
        Integer, Intent(In) :: iteration
        Integer :: mp, m, i
        Character*120 :: autostring, iterstring, cfile, checkfile
 
        Call chktmp%construct('p1a')
        chktmp%config = 'p1a'
        !Copy the RHS into chtkmp
        Call Get_All_RHS(chktmp%p1a)
        chktmp%p1a(:,:,:,numfields+1:numfields*2) = abterms(:,:,:,1:numfields)
        !Now we want to move from p1a to s2a (rlm space)
        Call chktmp%reform()

        ! This is where we cache, index by index...
        ! Don't forget to initialize the buffer:  Call Full_3D_Buffer%Init(mpi_tag=54) 

        If (ItIsTimeForAQuickSave) Then
            write(autostring,auto_fmt) (quicksave_num+1) !quick save number starts at 1
            checkpoint_prefix = 'Checkpoints/quicksave_'//trim(autostring)
        Else
            write(iterstring,int_out_fmt) iteration
            checkpoint_prefix = 'Checkpoints/'//trim(iterstring)
        Endif


        Do i = 1, numfields*2
            checkfile = Trim(my_path)//trim(checkpoint_prefix)//'_'//trim(checkpoint_suffix(i))
            !Write(6,*)'Checkpoint_suffix: ', checkpoint_suffix(i), checkfile
            Call checkpoint_buffer%cache_data_spectral(chktmp%s2a,i)
            Call checkpoint_buffer%write_data(filename=checkfile)

        Enddo
        !        Call Gen_Data()
        Call Check_Data()

        Call chktmp%deconstruct('s2a')

        If (my_rank .eq. 0) Then
            ! rank 0 writes out a file with the grid, etc.
            ! This file should contain everything that needs to be known

            cfile = Trim(my_path)//trim(checkpoint_prefix)//'_'//'grid_etc'

            open(unit=15,file=cfile,form='unformatted', status='replace')
            Write(15)n_r
            Write(15)grid_type
            Write(15)l_max
            Write(15)dt
            Write(15)new_dt
            Write(15)(radius(i),i=1,N_R)
            Write(15)elapsed_time
            Write(15)iteration
            Close(15)

            open(unit=15,file=Trim(my_path)//'Checkpoints/last_checkpoint',form='formatted', status='replace')
            If (ItIsTimeForAQuickSave) Then
                Write(15,int_minus_out_fmt)-iteration
                Write(15,'(i2.2)')(quicksave_num+1)
            Else
                Write(15,int_out_fmt)iteration
            Endif
            Close(15)

            open(unit=15,file=Trim(my_path)//'Checkpoints/checkpoint_log',form='formatted', status='unknown', &
                position='Append')
            If (ItIsTimeForAQuickSave) Then
                Write(iterstring,int_out_fmt)iteration
                Write(autostring,auto_fmt)quicksave_num+1
                Write(15,*)iterstring, ' ', autostring
            Else
                Write(15,int_out_fmt)iteration
            Endif
            Close(15)

        Endif

    End Subroutine Write_Checkpoint

    !~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Subroutine Read_Checkpoint_Buffer(fields, abterms,iteration,read_pars)
        Implicit None
        Integer, Intent(In) :: iteration, read_pars(1:2)
        Real*8, Intent(InOut) :: fields(:,:,:,:), abterms(:,:,:,:)
        Integer :: n_r_old, l_max_old, grid_type_old, nr_read
        Integer :: i, ierr, m, p, np, mp, lb,ub, f,  r, ind
        Integer :: old_pars(5), fcount(3,2)
        Integer :: last_iter, last_auto
        Integer :: read_magnetism = 0, read_hydro = 0
        Real*8 :: dt_pars(3),dt,new_dt
        Real*8, Allocatable :: old_radius(:), radius_old(:)
        Real*8, Allocatable :: tempfield1(:,:,:,:), tempfield2(:,:,:,:)
        Character*120 :: autostring, cfile, checkfile, dstring, iterstring

        read_hydro = read_pars(1)
        read_magnetism = read_pars(2)
        checkpoint_iter = iteration

        If (my_rank .eq. 0) Then
            old_pars(4) = checkpoint_iter
            old_pars(5) = -1
            If (checkpoint_iter .eq. 0) Then
                open(unit=15,file=Trim(my_path)//'Checkpoints/last_checkpoint',form='formatted', status='old')
                read(15,int_minus_in_fmt)last_iter
                If (last_iter .lt. 0) Then  !Indicates a quicksave
                    Read(15,'(i2.2)')last_auto
                    old_pars(4) = -last_iter
                    old_pars(5) = last_auto
                    Write(autostring,auto_fmt)last_auto
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
                Else
                    !Not a quicksave
                    old_pars(4) = last_iter
                    Write(iterstring,int_in_fmt) last_iter
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
                Endif

                Close(15)
            ElseIf (checkpoint_iter .lt. 0) Then
                !User has specified a particular quicksave file
                    last_auto = -checkpoint_iter
                    old_pars(5) = last_auto
                    Write(autostring,auto_fmt) last_auto
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
            Else
                Write(iterstring,int_in_fmt) iteration
                checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
            Endif

            !process zero reads all the old info and broadcasts to all other ranks
            cfile = Trim(checkpoint_prefix)//'_'//'grid_etc'
            open(unit=15,file=cfile,form='unformatted', status='old')
            Read(15)n_r_old
            Read(15)grid_type_old
            Read(15)l_max_old
            Read(15)dt
            Read(15)new_dt
            Allocate(old_radius(1:N_r_old))
            Read(15)(old_radius(i),i=1,N_R_old)
            Read(15)Checkpoint_time
            If (checkpoint_iter .lt. 0) Then
                ! We're loading a quicksave file
                ! Need to retrieve iteration from the grid_etc file because
                ! iteration specified in main_input was a low, negative number
                Read(15)Checkpoint_iter
                old_pars(4) = Checkpoint_iter
            Endif
            Close(15)

            write(dstring,sci_note_fmt)checkpoint_time
            call stdout%print(' ------ Checkpoint time is: '//trim(dstring))
            old_pars(1) = n_r_old
            old_pars(2) = grid_type_old
            old_pars(3) = l_max_old
            dt_pars(1) = dt
            dt_pars(2) = new_dt
            dt_pars(3) = checkpoint_time

            If (l_max_old .lt. l_max) Then
                    Write(6,*)' '
                    Write(6,*)'#####################################################################'
                    Write(6,*)'# '
                    Write(6,*)'#  Checkpoint horizontal resolution is lower than current resolution.'
                    Write(6,*)'#  The old solution will be interpolated onto horizontal grid with '
                    Write(6,*)'#  higher resolution corresponding to the new l_max.'
                    Write(6,*)'#  Old l_max: ', l_max_old
                    Write(6,*)'#  New l_max: ', l_max
                    Write(6,*)'# '
                    Write(6,*)'#####################################################################'
                    Write(6,*)' '
            Endif
            If (l_max_old .gt. l_max) Then
                    Write(6,*)' '
                    Write(6,*)'#####################################################################'
                    Write(6,*)'# '
                    Write(6,*)'#  Checkpoint horizontal resolution is higher than current resolution.'
                    Write(6,*)'#  The old SPH expansion will be truncated at the new l_max.'
                    Write(6,*)'#  This might not be a good idea.'
                    Write(6,*)'#  Old l_max: ', l_max_old
                    Write(6,*)'#  New l_max: ', l_max
                    Write(6,*)'# '
                    Write(6,*)'#####################################################################'
                    Write(6,*)' '
            Endif
        Endif

        If (my_row_rank .eq. 0) Then    !2-D broadcast pattern
            Call MPI_Bcast(old_pars,5, MPI_INTEGER, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(old_pars,5, MPI_INTEGER, 0, pfi%rcomm%comm, ierr)

        n_r_old       = old_pars(1)
        grid_type_old = old_pars(2)
        l_max_old     = old_pars(3)
        checkpoint_iter = old_pars(4)
        last_auto = old_pars(5)

        If (last_auto .ne. -1) Then
            !The prefix should be formed using quicksave
            Write(autostring,auto_fmt)last_auto
            checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
        Else
            !The prefix should reflect that this is a normal checkpoint file
            Write(iterstring,int_in_fmt) checkpoint_iter
            checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
        Endif

        !///////// Later we only want to do this if the grid is actually different
        If (my_rank .ne. 0) Then
            Allocate(old_radius(1:n_r_old))
        Endif

        If (my_row_rank .eq. 0) Then
            Call MPI_Bcast(old_radius,n_r_old, MPI_DOUBLE_PRECISION, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(old_radius,n_r_old, MPI_DOUBLE_PRECISION, 0, pfi%rcomm%comm, ierr)

        If (my_row_rank .eq. 0) Then
            Call MPI_Bcast(dt_pars,3, MPI_DOUBLE_PRECISION, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(dt_pars,3, MPI_DOUBLE_PRECISION, 0, pfi%rcomm%comm, ierr)

        checkpoint_dt    = dt_pars(1)
        checkpoint_newdt = dt_pars(2)
        checkpoint_Time  = dt_pars(3)

        Call chktmp%construct('s2b')
        chktmp%config = 's2b'
        Do mp = my_mp%min, my_mp%max
            chktmp%s2b(mp)%data(:,:,:,:) = 0.0d0
        Enddo


        ! Initialize the input buffer  (still need to deal properly with nr_old .ne. nr)

        Call checkpoint_inbuffer%Init(mpi_tag=checkpoint_tag,spectral=.true., cache_spectral = .true., &
                                    spec_comp = .true., lmax_in = l_max_old, mode = 2)
        Do i = 1, numfields*2
            checkfile = trim(checkpoint_prefix)//'_'//trim(checkpoint_suffix(i))
            Write(6,*)'Grabbing data from file: ', checkfile
            Call checkpoint_inbuffer%read_data(filename=checkfile)
            Call checkpoint_inbuffer%grab_data_spectral(chktmp%s2b,i)
        Enddo


        Call chktmp%reform()    ! move to p1b


        ! NOW, if n_r_old and grid_type_old are the same, we can copy chtkmp%p1b into abterms and
        ! fields.  Otherwise, we need to interpolate onto the current grid
        If  ((n_r_old .ne. n_r) .or. (grid_type_old .ne. grid_type) ) Then
            ! Interpolate
            ! We will assume the user kept the same radial domain bounds.
            ! If they  have not, this will end badly.
            If (my_rank .eq. 0) Then
                Write(6,*)'Grid has changed.  Interpolating onto new grid.'
                Write(6,*)'Old grid_type:     ', grid_type_old
                Write(6,*)'Current grid_type: ', grid_type
                Write(6,*)'Old N_R:           ', n_r_old
                Write(6,*)'Current N_R:           ', n_r
            Endif
            If (chebyshev) Then
                If (n_r_old .lt. n_r) Then

                    ! The fields are OK - they are already in chebyshev space
                    fields(:,:,:,1:numfields) = chktmp%p1b(:,:,:,1:numfields)

                    ! The AB terms are stored in physical space (in radius).
                    ! They need to be transformed, coefficients copied, and transformed back..
                    ! First, we need to initialize the old chebyshev grid.
                    Allocate(radius_old(1:n_r_old))
                    Call cheby_info%Init(radius_old,rmin,rmax)  ! We assume that rmax and rmin do not change
                    fcount(:,:) = numfields
                    Call chktmp2%init(field_count = fcount, config = 'p1a')
                    Call chktmp2%construct('p1a')
                    chktmp2%p1a(:,:,:,:) = 0.0d0
                    ! Allocate tempfield1, tempfield2
                    lb = lbound(chktmp%p1b,3)
                    ub = ubound(chktmp%p1b,3)
                    allocate(tempfield1(1:n_r_old,1:2,lb:ub,1))
                    allocate(tempfield2(1:n_r_old,1:2,lb:ub,1))

                    Do i = 1, numfields
                        tempfield1(:,:,:,:) = 0.0d0
                        tempfield2(:,:,:,:) = 0.0d0
                        tempfield1(1:n_r_old,:,:,1) = chktmp%p1b(1:n_r_old,:,:,numfields+i)
                        call cheby_info%tospec4d(tempfield1,tempfield2)
                        chktmp2%p1a(1:n_r_old,:,:,i) = tempfield2(1:n_r_old,:,:,1)
                    Enddo
                    DeAllocate(tempfield1,tempfield2)


                    Call chktmp2%construct('p1b')
                    !Normal transform(p1a,p1b)
                    Call gridcp%From_Spectral(chktmp2%p1a,chktmp2%p1b)

                    abterms(:,:,:,1:numfields) = chktmp2%p1b(:,:,:,1:numfields)
                    Call cheby_info%destroy()
                    Call chktmp2%deconstruct('p1a')
                    Call chktmp2%deconstruct('p1b')
                    Deallocate(radius_old)
                Endif
            Else
                Write(6,*)'Interpolation for FD not supported yet'
            Endif
        Else

            ! Interpolation is complete, now we just copy into the other arrays
            fields(:,:,:,1:numfields) = chktmp%p1b(:,:,:,1:numfields)
            abterms(:,:,:,1:numfields) = chktmp%p1b(:,:,:,numfields+1:numfields*2)

        Endif
        Call chktmp%deconstruct('p1b')
        DeAllocate(old_radius)

    End Subroutine Read_Checkpoint_Buffer

    Subroutine IsItTimeForACheckpoint(iter)
        Implicit None
        Integer, Intent(In) :: iter
        ItIsTimeForACheckpoint = .false.
        ItIsTimeForAQuickSave = .false.
        If (Mod(iter,checkpoint_interval) .eq. 0) Then
            ItIsTimeForACheckpoint = .true.
            checkpoint_t0 = checkpoint_elapsed      ! quicksaves not written
            checkpoint_elapsed = 0.0d0
            !If the long interval check is satisfied, nothing,
            ! nothing related to the short interval is executed.
        Else

            !Check for quick-save status.  This will be based on either iteration #
            ! OR on the time since the last checkpoint

            If (quicksave_interval .gt. 0) Then
                If (Mod(iter,quicksave_interval) .eq. 0) Then
                    ItIsTimeForACheckpoint = .true.
                    ItIsTimeForAQuickSave = .true.
                    quicksave_num = quicksave_num+1
                    quicksave_num = Mod(quicksave_num,num_quicksaves)

                Endif
            Endif

            If (quicksave_seconds .gt. 0) Then
                checkpoint_elapsed = global_msgs(2) - checkpoint_t0
                If (checkpoint_elapsed .gt. quicksave_seconds) Then

                    checkpoint_t0 = global_msgs(2)
                    checkpoint_elapsed = 0.0d0
                    ItIsTimeForACheckpoint = .true.
                    ItIsTimeForAQuickSave = .true.
                    quicksave_num = quicksave_num+1
                    quicksave_num = Mod(quicksave_num,num_quicksaves)

                Endif

            Endif

        Endif
    End Subroutine



    !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


    Subroutine Read_Checkpoint(fields, abterms,iteration,read_pars)
        Implicit None
        Integer, Intent(In) :: iteration, read_pars(1:2)
        Real*8, Intent(InOut) :: fields(:,:,:,:), abterms(:,:,:,:)
        Integer :: n_r_old, l_max_old, grid_type_old, nr_read
        Integer :: i, ierr, nlm_total_old, m, nl,p, np, mxread
        Integer :: maxl, dim2,offset, nl_load,lstart,mp, offset_index
        Integer :: old_pars(5)
        Integer, Allocatable :: lmstart_old(:)
        Real*8, Allocatable :: old_radius(:), radius_old(:)
        Real*8, Allocatable :: rowstrip(:,:), myarr(:,:), sendarr(:,:)
        Real*8 :: dt_pars(3),dt,new_dt
        Real*8, Allocatable :: tempfield1(:,:,:,:), tempfield2(:,:,:,:)
        Character*120 :: iterstring
        Character*120 :: autostring
        Character*120 :: dstring
        Character*120 :: cfile
        Integer :: fcount(3,2)
        Integer :: lb,ub, f, imi, r, ind
        Integer :: last_iter, last_auto
        Real*8 :: mxvt, mxvp
        Integer :: read_magnetism = 0, read_hydro = 0

        read_hydro = read_pars(1)
        read_magnetism = read_pars(2)

        dim2 = tnr*numfields*2
        checkpoint_iter = iteration

        If (my_rank .eq. 0) Then
            old_pars(4) = checkpoint_iter
            old_pars(5) = -1
            If (checkpoint_iter .eq. 0) Then
                open(unit=15,file=Trim(my_path)//'Checkpoints/last_checkpoint',form='formatted', status='old')
                read(15,int_minus_in_fmt)last_iter
                If (last_iter .lt. 0) Then  !Indicates a quicksave
                    Read(15,'(i2.2)')last_auto
                    old_pars(4) = -last_iter
                    old_pars(5) = last_auto
                    Write(autostring,auto_fmt)last_auto
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
                Else
                    !Not a quicksave
                    old_pars(4) = last_iter
                    Write(iterstring,int_in_fmt) last_iter
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
                Endif

                Close(15)
            ElseIf (checkpoint_iter .lt. 0) Then
                !User has specified a particular quicksave file
                    last_auto = -checkpoint_iter
                    old_pars(5) = last_auto
                    Write(autostring,auto_fmt) last_auto
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
            Else
                Write(iterstring,int_in_fmt) iteration
                checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
            Endif


            !process zero reads all the old info and broadcasts to all other ranks
            cfile = Trim(checkpoint_prefix)//'_'//'grid_etc'
            open(unit=15,file=cfile,form='unformatted', status='old')
            Read(15)n_r_old
            Read(15)grid_type_old
            Read(15)l_max_old
            Read(15)dt
            Read(15)new_dt
            Allocate(old_radius(1:N_r_old))
            Read(15)(old_radius(i),i=1,N_R_old)
            Read(15)Checkpoint_time
            If (checkpoint_iter .lt. 0) Then
                ! We're loading a quicksave file
                ! Need to retrieve iteration from the grid_etc file because
                ! iteration specified in main_input was a low, negative number
                Read(15)Checkpoint_iter
                old_pars(4) = Checkpoint_iter
            Endif
            Close(15)

            write(dstring,sci_note_fmt)checkpoint_time
            call stdout%print(' ------ Checkpoint time is: '//trim(dstring))
            old_pars(1) = n_r_old
            old_pars(2) = grid_type_old
            old_pars(3) = l_max_old
            dt_pars(1) = dt
            dt_pars(2) = new_dt
            dt_pars(3) = checkpoint_time

            If (l_max_old .lt. l_max) Then
                    Write(6,*)' '
                    Write(6,*)'#####################################################################'
                    Write(6,*)'# '
                    Write(6,*)'#  Checkpoint horizontal resolution is lower than current resolution.'
                    Write(6,*)'#  The old solution will be interpolated onto horizontal grid with '
                    Write(6,*)'#  higher resolution corresponding to the new l_max.'
                    Write(6,*)'#  Old l_max: ', l_max_old
                    Write(6,*)'#  New l_max: ', l_max
                    Write(6,*)'# '
                    Write(6,*)'#####################################################################'
                    Write(6,*)' '
            Endif
            If (l_max_old .gt. l_max) Then
                    Write(6,*)' '
                    Write(6,*)'#####################################################################'
                    Write(6,*)'# '
                    Write(6,*)'#  Checkpoint horizontal resolution is higher than current resolution.'
                    Write(6,*)'#  The old SPH expansion will be truncated at the new l_max.'
                    Write(6,*)'#  This might not be a good idea.'
                    Write(6,*)'#  Old l_max: ', l_max_old
                    Write(6,*)'#  New l_max: ', l_max
                    Write(6,*)'# '
                    Write(6,*)'#####################################################################'
                    Write(6,*)' '
            Endif
        Endif

        If (my_row_rank .eq. 0) Then    !2-D broadcast pattern
            Call MPI_Bcast(old_pars,5, MPI_INTEGER, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(old_pars,5, MPI_INTEGER, 0, pfi%rcomm%comm, ierr)

        n_r_old       = old_pars(1)
        grid_type_old = old_pars(2)
        l_max_old     = old_pars(3)
        checkpoint_iter = old_pars(4)
        last_auto = old_pars(5)

        If (last_auto .ne. -1) Then
            !The prefix should be formed using quicksave
            Write(autostring,auto_fmt)last_auto
            checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
        Else
            !The prefix should reflect that this is a normal checkpoint file
            Write(iterstring,int_in_fmt) checkpoint_iter
            checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
        Endif


        !///////// Later we only want to do this if the grid is actually different
        If (my_rank .ne. 0) Then
            Allocate(old_radius(1:n_r_old))
        Endif

        If (my_row_rank .eq. 0) Then
            Call MPI_Bcast(old_radius,n_r_old, MPI_DOUBLE_PRECISION, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(old_radius,n_r_old, MPI_DOUBLE_PRECISION, 0, pfi%rcomm%comm, ierr)

        If (my_row_rank .eq. 0) Then
            Call MPI_Bcast(dt_pars,3, MPI_DOUBLE_PRECISION, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(dt_pars,3, MPI_DOUBLE_PRECISION, 0, pfi%rcomm%comm, ierr)

        checkpoint_dt    = dt_pars(1)
        checkpoint_newdt = dt_pars(2)
        checkpoint_Time  = dt_pars(3)

        Call chktmp%construct('s2b')
        chktmp%config = 's2b'
        Do mp = my_mp%min, my_mp%max
            chktmp%s2b(mp)%data(:,:,:,:) = 0.0d0
        Enddo

        !//////////////////////////////////
        ! Determine if we actually read anything
        ! The old nr might be lower than the current nr
        ! (currently does not support old nr being larger...)

        mxread = n_r_old-my_r%min+1
        nr_read = min(my_r%delta,mxread)  !  Was tnr, not nr....
        if(my_r%min .gt. n_r_old) nr_read = 0


        ! Rank zero from each row participates in the read
        ! Rank zero within each row is responsible for reading in all l-m combinations
        ! for the subset of radii that its row is responsible for
        If (my_row_rank .eq. 0) Then
            !If (nr_read .gt. 0) Then    ! SMALL BUG HERE RELATED TO MPI_IO Logic... -- Revist this



                ! Column zero reads in the old checkpoint no matter what
                ! the old dimensions are. We'll broadcast back to all members of the row later

                ! Interpolation/truncation of the old spherical harmonic basis is easy
                ! Interpolation up in radius is easy, but down is more difficult to program
                ! and unlikely to be used.  I will write a serial version of the input
                ! for those cases where interpolation down is required.

                nlm_total_old = 0            ! This is the total number of l-m combinations in the CHECKPOINT FILE
                Do m = 0, l_max_old
                    nl = l_max_old-m+1
                    nlm_total_old = nlm_total_old+nl
                Enddo

                np = pfi%ccomm%np
                my_in_disp = 0
                Do p = 1, my_column_rank
                    my_in_disp = my_in_disp+pfi%all_1p(p-1)%delta
                enddo

                my_in_disp = my_in_disp*nlm_total_old
                full_in_disp = nlm_total_old*n_r_old

                ! tnr is the number of radial points for this row x 2 (for complex values)
                ! In addition to holding each of the 4 (or 6 in MHD) fields, the rowstrip
                ! array will also hold the Adams-Bashforth arrays associated with each field
                ! (hence the ADDITIONAL factor of 2 below).
                Allocate( rowstrip(1:nlm_total_old, 1:tnr*numfields*2))
                rowstrip(:,:) = 0

                If (read_hydro .eq. 1) Then
                    Call Read_Field(rowstrip,1,wchar, iteration,nr_read,nlm_total_old)
                    Call Read_Field(rowstrip,2,pchar, iteration,nr_read,nlm_total_old)
                    Call Read_Field(rowstrip,3,tchar, iteration,nr_read,nlm_total_old)
                    Call Read_Field(rowstrip,4,zchar, iteration,nr_read,nlm_total_old)
                Endif
                offset_index = 4
                If (magnetism) Then
                    If (read_magnetism .eq. 1) Then
                        Call Read_Field(rowstrip,5,cchar, iteration,nr_read,nlm_total_old)
                        Call Read_Field(rowstrip,6,achar, iteration,nr_read,nlm_total_old)
                    Endif
                    offset_index = 6
                Endif

                If (read_hydro .eq. 1) Then
                    Call Read_Field(rowstrip,offset_index+1,'WAB', iteration,nr_read,nlm_total_old)
                    Call Read_Field(rowstrip,offset_index+2,'PAB', iteration,nr_read,nlm_total_old)
                    Call Read_Field(rowstrip,offset_index+3,'TAB', iteration,nr_read,nlm_total_old)
                    Call Read_Field(rowstrip,offset_index+4,'ZAB', iteration,nr_read,nlm_total_old)
                Endif
                If (magnetism) Then
                    If (read_magnetism .eq. 1) Then
                        Call Read_Field(rowstrip,offset_index+5,'CAB', iteration,nr_read,nlm_total_old)
                        Call Read_Field(rowstrip,offset_index+6,'AAB', iteration,nr_read,nlm_total_old)
                    Endif
                Endif

                ! Now the head of each row owns all modes of each field at the
                ! radii owned by that row.  The different modes now need to be
                ! distributed to their respect owners.
                ! This is where horizontal interpolation or truncation is done (implicitly)
                ! by only loading  appropriate modes into the send arrays

                !/////////////////////////////
                !  Do some book-keeping related to the old l_max
                !  This array tells us where in the first dimension of rowstrip, ell values for mode m start
                np = pfi%rcomm%np

                Allocate(lmstart_old(0:l_max_old))
                lmstart_old(0) = 1
                do m = 0, l_max_old-1
                    nl = l_max_old-m+1
                    lmstart_old(m+1) = lmstart_old(m)+nl
                enddo

                maxl = min(l_max,l_max_old)    ! Take care to only read in modes that are common
                                                        ! to the checkpoint & and the current simulation

                ! First, each row-head pulls out their own modes
                Do mp = my_mp%min, my_mp%max
                    m = m_values(mp)
                    chktmp%s2b(mp)%data(:,:,:,:) = 0.0d0
                    If (m .le. l_max_old) Then
                        nl = maxl-m+1
                        lstart = lmstart_old(m)
                        ind = 1
                        Do f = 1, numfields*2
                        Do imi = 1, 2
                        Do r = my_r%min, my_r%max
                            chktmp%s2b(mp)%data(m:maxl,r,imi,f) = rowstrip(lstart:lstart+nl-1,ind)
                            ind = ind+1
                        Enddo
                        Enddo
                        Enddo
                    Endif
                Enddo


                !/////////////////////////////////////////////////////////////
                ! Next send each of the other processors in the row their information
                Do p = 1, np -1
                        Allocate(sendarr(1:mode_count(p),1:dim2))
                        sendarr(:,:) = 0.0d0
                        offset = 1
                        Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
                            m = m_values(mp)
                            nl = l_max-m+1
                            If (m .le. l_max_old) Then
                                nl_load = maxl-m+1
                                lstart = lmstart_old(m)
                                sendarr(offset:offset+nl_load-1,:)    = rowstrip(lstart:lstart+nl_load-1,:)
                            Endif
                            offset = offset+nl
                        Enddo
                   Call send(sendarr, dest= p,tag=checkpoint_tag,grp = pfi%rcomm)
                        DeAllocate(sendarr)
                Enddo
                !/////////

                DeAllocate(lmstart_old)
                DeAllocate(rowstrip)
            !Endif
        Else
            ! Receive my modes
            !If (nr_read .gt. 0) Then
                Allocate(myarr(1:mode_count(my_row_rank),1:dim2))
                Call receive(myarr, source= 0,tag=checkpoint_tag,grp = pfi%rcomm)
                offset =1
                Do mp = my_mp%min, my_mp%max
                    m = m_values(mp)
                    nl = l_max-m+1
                        ind = 1
                    Do f = 1, numfields*2
                    Do imi = 1, 2
                    Do r = my_r%min, my_r%max
                        chktmp%s2b(mp)%data(m:l_max,r,imi,f) = myarr(offset:offset+nl-1,ind)
                        ind = ind+1
                    Enddo
                    Enddo
                    Enddo
                    offset = offset+nl
                Enddo
                DeAllocate(myarr)
            !Endif
        Endif



        Call chktmp%reform()    ! move to p1b



        ! NOW, if n_r_old and grid_type_old are the same, we can copy chtkmp%p1b into abterms and
        ! fields.  Otherwise, we need to interpolate onto the current grid
        If  ((n_r_old .ne. n_r) .or. (grid_type_old .ne. grid_type) ) Then
            ! Interpolate
            ! We will assume the user kept the same radial domain bounds.
            ! If they  have not, this will end badly.
            If (my_rank .eq. 0) Then
                Write(6,*)'Grid has changed.  Interpolating onto new grid.'
                Write(6,*)'Old grid_type:     ', grid_type_old
                Write(6,*)'Current grid_type: ', grid_type
                Write(6,*)'Old N_R:           ', n_r_old
                Write(6,*)'Current N_R:           ', n_r
            Endif
            If (chebyshev) Then
                If (n_r_old .lt. n_r) Then

                    ! The fields are OK - they are already in chebyshev space
                    fields(:,:,:,1:numfields) = chktmp%p1b(:,:,:,1:numfields)

                    ! The AB terms are stored in physical space (in radius).
                    ! They need to be transformed, coefficients copied, and transformed back..
                    ! First, we need to initialize the old chebyshev grid.
                    Allocate(radius_old(1:n_r_old))
                    Call cheby_info%Init(radius_old,rmin,rmax)  ! We assume that rmax and rmin do not change
                    fcount(:,:) = numfields
                    Call chktmp2%init(field_count = fcount, config = 'p1a')
                    Call chktmp2%construct('p1a')
                    chktmp2%p1a(:,:,:,:) = 0.0d0
                    ! Allocate tempfield1, tempfield2
                    lb = lbound(chktmp%p1b,3)
                    ub = ubound(chktmp%p1b,3)
                    allocate(tempfield1(1:n_r_old,1:2,lb:ub,1))
                    allocate(tempfield2(1:n_r_old,1:2,lb:ub,1))

                    Do i = 1, numfields
                        tempfield1(:,:,:,:) = 0.0d0
                        tempfield2(:,:,:,:) = 0.0d0
                        tempfield1(1:n_r_old,:,:,1) = chktmp%p1b(1:n_r_old,:,:,numfields+i)
                        call cheby_info%tospec4d(tempfield1,tempfield2)
                        chktmp2%p1a(1:n_r_old,:,:,i) = tempfield2(1:n_r_old,:,:,1)
                    Enddo
                    DeAllocate(tempfield1,tempfield2)


                    Call chktmp2%construct('p1b')
                    !Normal transform(p1a,p1b)
                    Call gridcp%From_Spectral(chktmp2%p1a,chktmp2%p1b)

                    abterms(:,:,:,1:numfields) = chktmp2%p1b(:,:,:,1:numfields)
                    Call cheby_info%destroy()
                    Call chktmp2%deconstruct('p1a')
                    Call chktmp2%deconstruct('p1b')
                    Deallocate(radius_old)
                Endif
            Else
                Write(6,*)'Interpolation for FD not supported yet'
            Endif
        Else

            ! Interpolation is complete, now we just copy into the other arrays
            fields(:,:,:,1:numfields) = chktmp%p1b(:,:,:,1:numfields)
            abterms(:,:,:,1:numfields) = chktmp%p1b(:,:,:,numfields+1:numfields*2)

        Endif
        Call chktmp%deconstruct('p1b')
        DeAllocate(old_radius)

    End Subroutine Read_Checkpoint

    Subroutine Read_Field(arr,ind,tag,iter,nread,nlm)
        Implicit None
        Integer, Intent(In) :: ind, iter,nread, nlm
        Real*8, Intent(InOut) :: arr(1:,1:)
        Character*3, Intent(In) :: tag
        Character*120 :: cfile
        Integer ::bsize_in
        integer ierr, funit , v_offset1, v_offset2
        integer(kind=MPI_OFFSET_KIND) disp1,disp2
        Integer :: mstatus(MPI_STATUS_SIZE)

        cfile = Trim(checkpoint_prefix)//'_'//trim(tag)

        v_offset1 = (ind-1)*tnr+1
        v_offset2 = v_offset1+my_r%delta

        call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, &
        MPI_MODE_RDONLY, &
        MPI_INFO_NULL, funit, ierr)
            if (ierr .ne. 0) Then
                Write(6,*)'Error opening: ', pfi%ccomm%rank
            Endif
        bsize_in = nread*nlm
        If (ierr .eq. 0) Then

            disp1 = my_in_disp*8
            disp2 = (my_in_disp+full_in_disp)*8
            call MPI_FILE_SET_VIEW(funit, disp1, MPI_DOUBLE_PRECISION, &
                MPI_DOUBLE_PRECISION, 'native', &
                MPI_INFO_NULL, ierr)
            if (ierr .ne. 0) Then
                Write(6,*)'Error setting view1: ', pfi%ccomm%rank
            Endif
            If (nread .gt. 0) Then
            call MPI_FILE_READ(funit, arr(1,v_offset1), bsize_in, MPI_DOUBLE_PRECISION, &
            mstatus, ierr)
            if (ierr .ne. 0) Then
                Write(6,*)'Error reading1: ', pfi%ccomm%rank
            Endif
            Endif
            call MPI_FILE_SET_VIEW(funit, disp2, MPI_DOUBLE_PRECISION, &
                MPI_DOUBLE_PRECISION, 'native', &
                MPI_INFO_NULL, ierr)
            if (ierr .ne. 0) Then
                Write(6,*)'Error setting view2: ', pfi%ccomm%rank
            Endif
            If (nread .gt. 0) Then
            call MPI_FILE_READ(funit, arr(1,v_offset2), bsize_in, MPI_DOUBLE_PRECISION, &
            mstatus, ierr)
            if (ierr .ne. 0) Then
                Write(6,*)'Error reading2: ', pfi%ccomm%rank
            Endif
            Endif

        Else
            If (my_rank .eq. 0) Then
                Write(6,*)'File read error.  Associated array will contain only zeroes: ', cfile
            Endif
        Endif
        call MPI_FILE_CLOSE(funit, ierr)
        if (ierr .ne. 0) Then
                Write(6,*)'Error closing file: ', pfi%ccomm%rank
        Endif
    End Subroutine Read_Field



End Module Checkpointing
