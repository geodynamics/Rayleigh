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

    Integer :: checkpoint_iter = 0
    Real*8 :: checkpoint_dt, checkpoint_newdt
    Real*8 :: checkpoint_time
    !//////////////////////////
    ! New variables for new checkpointing style
    Integer :: Noutputs_per_row, Nradii_per_output
    Integer, Allocatable :: Nradii_at_rank(:), Rstart_at_rank(:)
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
        Endif
        nfs(:) = numfields*2
        Call chktmp%init(field_count = nfs, config = 'p1a')            ! This structure hangs around through the entire run

        !////////////////////////

        !//////////////////////////////////////////////////
        ! Revisions for mem-friendly I/O

        !Write(6,*)'NPOUT: ', npout
        Noutputs_per_row = npout! # of processors outputting within row
        If (Noutputs_per_row .gt. my_r%delta) Then
            Noutputs_per_row = my_r%delta
        Endif
        If (Noutputs_per_row .gt. nprow) Then
                Noutputs_per_row = nprow
        Endif
        Nradii_per_output = my_r%delta/Noutputs_per_row    ! Number of radii each processor outputs
        rextra = Mod(my_r%delta,Noutputs_per_row)           ! Remainder
        Allocate(nradii_at_rank(0:Noutputs_per_row-1))    ! Number of radii output by each rank
        Allocate(rstart_at_rank(0:Noutputs_per_row-1))    ! First radial index (local) this rank receives
        Do p = 0, Noutputs_per_row -1
            Nradii_at_rank(p) = Nradii_per_output
            If (rextra .gt. 0) Then
                If (p .lt. rextra) Then
                    Nradii_at_rank(p) = Nradii_at_rank(p)+1
                Endif
            Endif
        Enddo
        rstart_at_rank(0) = 1
        Do p = 1, Noutputs_per_row-1
            rstart_at_rank(p) = rstart_at_rank(p-1)+nradii_at_rank(p-1)
        Enddo
        ! Determine if each rank will output
        If (my_row_rank .lt. Noutputs_per_row) Then
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
            my_check_disp2 = (my_check_disp+rstart_at_rank(my_row_rank)-1)*nlm_total    ! for new checkpointing style
            my_check_disp = my_check_disp*nlm_total  ! for old Checkpointing style
            buffsize2 = nlm_total*nradii_at_rank(my_row_rank)    ! for new
            buffsize = nlm_total*my_r%delta  ! for old
            full_disp = N_r
            full_disp = full_disp*nlm_total
        Endif

    End Subroutine Initialize_Checkpointing



    Subroutine Write_Checkpoint(abterms,iteration,dt,new_dt,elapsed_time)
        Implicit None
        Real*8, Intent(In) :: abterms(:,:,:,:)
        Real*8, Intent(In) :: dt, new_dt
        Integer, Intent(In) :: iteration
        Integer :: mp, m, offset,nl,p,np, f, imi, r, ind
        Integer :: dim2, lstart, i, offset_index
        Real*8, Allocatable :: myarr(:,:), rowstrip(:,:)
        Real*8, Intent(In) :: elapsed_time
        Character*2 :: autostring
        Character*8 :: iterstring
        Character*120 :: cfile
        np = pfi%rcomm%np

        Call chktmp%construct('p1a')
        chktmp%config = 'p1a'
        !Copy the RHS into chtkmp
        Call Get_All_RHS(chktmp%p1a)
        chktmp%p1a(:,:,:,numfields+1:numfields*2) = abterms(:,:,:,1:numfields)
        !Now we want to move from p1a to s2a (rlm space)
        Call chktmp%reform()

        ! Next, each process stripes their s2a array into a true 2-D array
        dim2 = tnr*numfields*2
        Allocate(myarr(1:mode_count(my_row_rank),1:dim2))
        offset =1
        Do mp = my_mp%min, my_mp%max
                m = m_values(mp)
                nl = l_max-m+1
                ind = 1
                Do f = 1, numfields*2
                Do imi = 1, 2
                Do r = my_r%min, my_r%max
                myarr(offset:offset+nl-1,ind) = chktmp%s2a(mp)%data(m:l_max,r,imi,f)
                ind = ind+1
                Enddo
                Enddo
                Enddo

                offset = offset+nl
        Enddo
        Call chktmp%deconstruct('s2a')

        ! Everyone sends to then 0 process of each row, who organizes the data into one large strip for output.


        If (my_row_rank .ne. 0) Then
            ! Send myarr

            Call Send(myarr, dest = 0, tag = checkpoint_tag, grp = pfi%rcomm)
            DeAllocate(myarr)
        Else
            Allocate( rowstrip(1:nlm_total, 1:tnr*numfields*2))
            ! first, copy myarr into the larger array
            offset = 1
            Do mp = my_mp%min, my_mp%max
                m = m_values(mp)
                nl = l_max-m+1
                lstart = lmstart(m)
                rowstrip(lstart:lstart+nl-1,:) = myarr(offset:offset+nl-1,:)
                offset = offset+nl
            Enddo
            DeAllocate(myarr)
            Do p = 1, np -1
                    Allocate(myarr(1:mode_count(p),1:dim2))
                    ! Receive
                    Call receive(myarr, source= p,tag=checkpoint_tag,grp = pfi%rcomm)
                    offset = 1
                    Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
                        m = m_values(mp)
                        nl = l_max-m+1
                        lstart = lmstart(m)
                        rowstrip(lstart:lstart+nl-1,:) = myarr(offset:offset+nl-1,:)
                        offset = offset+nl
                    Enddo
                    DeAllocate(myarr)
            Enddo
                If (ItIsTimeForAQuickSave) Then
                    write(autostring,'(i2.2)') (quicksave_num+1) !quick save number starts at 1
                    checkpoint_prefix = 'Checkpoints/quicksave_'//trim(autostring)
                Else
                    write(iterstring,'(i8.8)') iteration
                    checkpoint_prefix = 'Checkpoints/'//trim(iterstring)
                Endif

                Call Write_Field(rowstrip,1,wchar, iteration)
                Call Write_Field(rowstrip,2,pchar, iteration)
                Call Write_Field(rowstrip,3,tchar, iteration)
                Call Write_Field(rowstrip,4,zchar, iteration)
                offset_index = 4
                If (magnetism) Then
                    Call Write_Field(rowstrip,5,cchar, iteration)
                    Call Write_Field(rowstrip,6,achar, iteration)
                    offset_index = 6
                Endif

                Call Write_Field(rowstrip,offset_index+1,'WAB', iteration)
                Call Write_Field(rowstrip,offset_index+2,'PAB', iteration)
                Call Write_Field(rowstrip,offset_index+3,'TAB', iteration)
                Call Write_Field(rowstrip,offset_index+4,'ZAB', iteration)
                If (magnetism) Then
                    Call Write_Field(rowstrip,offset_index+5,'CAB', iteration)
                    Call Write_Field(rowstrip,offset_index+6,'AAB', iteration)
                Endif
            DeAllocate(rowstrip)
            If (my_column_rank .eq. 0) Then
                ! row/column 0 writes out a file with the grid, etc.
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
                    !Write(iterstring,'(i8.8)')iteration
                    !Write(autostring,'(i2.2)')quicksave_num+1
                    !Write(15,*)iterstring, ' ', autostring
                    Write(15,'(i9.8)')-iteration
                    Write(15,'(i2.2)')(quicksave_num+1)
                Else
                    Write(15,'(i8.8)')iteration
                Endif
                Close(15)

                open(unit=15,file=Trim(my_path)//'Checkpoints/checkpoint_log',form='formatted', status='unknown', &
                    position='Append')
                If (ItIsTimeForAQuickSave) Then
                    !Write(15,'(i9.8)')-iteration
                    !Write(15,'(i2.2)')(quicksave_num+1)
                    Write(iterstring,'(i8.8)')iteration
                    Write(autostring,'(i2.2)')quicksave_num+1
                    Write(15,*)iterstring, ' ', autostring

                Else
                    Write(15,'(i8.8)')iteration
                Endif
                Close(15)

            Endif
        Endif


    End Subroutine Write_Checkpoint

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
        Character*8 :: iterstring
        Character*2 :: autostring
        Character*12 :: dstring
        Character*8 :: dofmt = '(ES12.5)'
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
        !Write(iterstring,'(i8.8)') iteration
        If (my_rank .eq. 0) Then
            old_pars(4) = checkpoint_iter
            old_pars(5) = -1
            If (checkpoint_iter .eq. 0) Then
                open(unit=15,file=Trim(my_path)//'Checkpoints/last_checkpoint',form='formatted', status='old')
                read(15,'(i9.8)')last_iter
                If (last_iter .lt. 0) Then  !Indicates a quicksave
                    Read(15,'(i2.2)')last_auto
                    old_pars(4) = -last_iter
                    old_pars(5) = last_auto
                    Write(autostring,'(i2.2)')last_auto
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
                Else
                    !Not a quicksave
                    old_pars(4) = last_iter
                    Write(iterstring,'(i8.8)') last_iter
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
                Endif

                Close(15)
            ElseIf (checkpoint_iter .lt. 0) Then
                !User has specified a particular quicksave file
                    last_auto = -checkpoint_iter
                    old_pars(5) = -checkpoint_iter
                    Write(autostring,'(i2.2)')-checkpoint_iter
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
            Else
                Write(iterstring,'(i8.8)') iteration
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

            write(dstring,dofmt)checkpoint_time
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
            Write(autostring,'(i2.2)')last_auto
            checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
        Else
            !The prefix should reflect that this is a normal checkpoint file
            Write(iterstring,'(i8.8)') checkpoint_iter
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

    Subroutine Write_Field(arr,ind,tag,iter)
                Implicit None
                Integer, Intent(In) :: ind, iter
                Real*8, Intent(In) :: arr(1:,1:)
                Character*2 :: autostring
                Character*8 :: iterstring
                Character*3, Intent(In) :: tag
                Character*120 :: cfile

                integer ierr, funit , v_offset1, v_offset2
                integer(kind=MPI_OFFSET_KIND) disp1,disp2
                Integer :: mstatus(MPI_STATUS_SIZE)



                cfile = Trim(my_path)//Trim(checkpoint_prefix)//'_'//Trim(tag)
                 ! We have to be careful here.  Each processor does TWO writes.
                ! The first write places the real part of the field into the file.
                ! The view then changes and advances to the appropriate location of the
                ! imaginary part.  This step is crucial for checkpoints to work with
                ! Different processor configurations.
                 v_offset1 = (ind-1)*tnr+1
                v_offset2 = v_offset1+my_r%delta

                call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, funit, ierr)
                if (ierr .ne. 0) Then
                    Write(6,*)'Error Opening File: ', pfi%ccomm%rank
                Endif

                disp1 = my_check_disp*8
                disp2 = (my_check_disp+full_disp)*8


                Call MPI_File_Seek(funit,disp1,MPI_SEEK_SET,ierr)
                If (ierr .ne. 0) Write(6,*)'Error Seeking 1: ', pfi%ccomm%rank


                Call MPI_FILE_WRITE(funit, arr(1,v_offset1), buffsize, MPI_DOUBLE_PRECISION, &
                        mstatus, ierr)
                If (ierr .ne. 0) Write(6,*)'Error Writing 1: ', pfi%ccomm%rank


                Call MPI_File_Seek(funit,disp2,MPI_SEEK_SET,ierr)
                If (ierr .ne. 0) Write(6,*)'Error Seeking 2: ', pfi%ccomm%rank


                Call MPI_FILE_WRITE(funit, arr(1,v_offset2), buffsize, MPI_DOUBLE_PRECISION, &
                        mstatus, ierr)
                If (ierr .ne. 0) Write(6,*)'Error Writing 2: ', pfi%ccomm%rank


                Call MPI_FILE_CLOSE(funit, ierr)
                if (ierr .ne. 0) Write(6,*)'Error Closing File: ', pfi%ccomm%rank





    End Subroutine Write_Field





    Subroutine Read_Field(arr,ind,tag,iter,nread,nlm)
        Implicit None
        Integer, Intent(In) :: ind, iter,nread, nlm

        Real*8, Intent(InOut) :: arr(1:,1:)
        Character*8 :: iterstring
        Character*3, Intent(In) :: tag
        Character*120 :: cfile
        Integer ::bsize_in
        integer ierr, funit , v_offset1, v_offset2
        integer(kind=MPI_OFFSET_KIND) disp1,disp2
        Integer :: mstatus(MPI_STATUS_SIZE)
        !write(iterstring,'(i8.8)') iter
        !cfile = Trim(my_path)//'Checkpoints/'//trim(iterstring)//'_'//trim(tag)
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
       ! If (nread .gt. 0) Then     !!! SET_VIEW OR _READ ARE BLOCKING APPRENTLY?...
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
       ! Endif
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

    Subroutine Write_Spectral_Field3D(arrin,ind,tag)
        ! Parallel Writing Routine For Fields in Spectral rlm configuration
        ! ISends and IReceives are used
        Implicit None
        Integer :: var_offset, offset, new_off


        Integer :: nrirq, nsirq, irq_ind, rtag, stag
        Integer, Allocatable :: rirqs(:), sirqs(:)
        Integer :: rone, rtwo, my_nrad, num_el
        Integer :: i, mp, m, lstart, nl,p,r
        Integer :: indstart(2)
        Real*8, Intent(In) :: arrin(1:,1:)
        Real*8, Allocatable :: arr(:,:,:), tarr(:)


                Integer, Intent(In) :: ind
                Character*8 :: iterstring
                Character*3, Intent(In) :: tag
                Character*120 :: cfile

                integer ierr, funit
                integer(kind=MPI_OFFSET_KIND) disp1,disp2
                Integer :: mstatus(MPI_STATUS_SIZE)
            cfile = Trim(checkpoint_prefix)//'_'//trim(tag)


        var_offset = (ind-1)*tnr
        If (I_Will_Output) Then
            my_nrad = nradii_at_rank(my_row_rank)
            Allocate( arr(1:nlm_total,1:my_nrad,2))
            Allocate(tarr(1:nlm_total*my_nrad))

            nrirq = nprow-1
            Allocate(rirqs(1:nrirq))
            nsirq = Noutputs_per_row-1
            Allocate(sirqs(1:nsirq))

            Do i = 1, 2        ! 2 passess.  Real and imaginary.   Could do one, but worried about Memory.
                rone = var_offset+rstart_at_rank(my_row_rank)+(i-1)*my_r%delta
                rtwo = rone+my_nrad-1
                ! Post receives for everyone in my row
                ! tag depends on i & p
                offset = 1
                irq_ind = 1
                Do p = 0, nprow-1
                        If (p .ne. my_row_rank) Then
                            rtag = p*(i+1)
                            num_el = mode_count(p)*my_nrad
                            Call IReceive(tarr, rirqs(irq_ind),num_el, offset,p, rtag, pfi%rcomm)
                            irq_ind = irq_ind+1
                        Else
                            ! File my stuff into tarr here
                            new_off = offset
                            Do r = rone, rtwo
                                Do m = 1, mode_count(p)
                                    tarr(new_off) = arrin(m,r)
                                    new_off = new_off+1
                                Enddo
                            Enddo
                        Endif

                        offset = offset+mode_count(p)*my_nrad
                Enddo


                ! Post sends to all output processes in my row
                irq_ind = 1
                Do p = 0, Noutputs_per_Row-1
                    If (p .ne. my_row_rank) Then
                        num_el = mode_count(my_row_rank)*nradii_at_rank(p)
                        stag = my_row_rank*(i+1)
                        indstart(1) = 1
                        indstart(2) = var_offset+rstart_at_rank(p)+(i-1)*my_r%delta
                        Call ISend(arrin, sirqs(irq_ind),num_el, p, stag, pfi%rcomm, indstart)
                        irq_ind = irq_ind+1
                    Endif
                Enddo

                ! Block and wait on receive irqs

                Call IWaitAll(nrirq, rirqs)

                ! File Data into Main Array
                offset = 1
                Do p = 0, nprow-1
                    Do r = 1, my_nrad
                    Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
                        m = m_values(mp)
                        nl = l_max-m+1
                        lstart = lmstart(m)
                        arr(lstart:lstart+nl-1,r,i) = tarr(offset:offset+nl-1)
                        offset = offset+nl
                    Enddo
                    Enddo
                Enddo

                Call IWaitAll(nsirq, sirqs)


            Enddo
            DeAllocate(sirqs)
            DeAllocate(tarr)
            DeAllocate(rirqs)
            !//////////////////////////////
            !  MPI Write




             ! We have to be careful here.  Each processor does TWO writes.
            ! The first write places the real part of the field into the file.
            ! The view then changes and advances to the appropriate location of the
            ! imaginary part.  This step is crucial for checkpoints to work with
            ! Different processor configurations.   The Real stuff sits at the
            ! Beginning of the file.  The imaginary stuff sits at the end.


            !///// NEED TO SET THESE DISPLACEMENTS

            Call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, &
                MPI_MODE_WRONLY + MPI_MODE_CREATE, &
            MPI_INFO_NULL, funit, ierr)
            disp1 = my_check_disp2*8
            disp2 = (my_check_disp2+full_disp)*8

            Call MPI_FILE_SET_VIEW(funit, disp1, MPI_DOUBLE_PRECISION, &     ! Real part
              MPI_DOUBLE_PRECISION, 'native', &
              MPI_INFO_NULL, ierr)
            Call MPI_FILE_WRITE(funit, arr(1,1,1), buffsize2, MPI_DOUBLE_PRECISION, &
              mstatus, ierr)

            Call MPI_FILE_SET_VIEW(funit, disp2, MPI_DOUBLE_PRECISION, &     ! Imaginary part
              MPI_DOUBLE_PRECISION, 'native', &
              MPI_INFO_NULL, ierr)
            Call MPI_FILE_WRITE(funit, arr(1,1,2), buffsize2, MPI_DOUBLE_PRECISION, &
              mstatus, ierr)

            Call MPI_FILE_CLOSE(funit, ierr)


            !////////////////////////////////
            DeAllocate(arr)
        Else
            ! This rank does not output
            ! Post sends to all output processes in my row
            nsirq = Noutputs_per_row
            Allocate(sirqs(1:nsirq))

            Do i = 1, 2
                irq_ind = 1
                Do p = 0, Noutputs_per_Row-1
                    num_el = mode_count(my_row_rank)*nradii_at_rank(p)
                    stag = my_row_rank*(i+1)
                    indstart(1) = 1
                    indstart(2) = var_offset+rstart_at_rank(p)+(i-1)*my_r%delta
                    Call ISend(arrin, sirqs(irq_ind),num_el, p, stag, pfi%rcomm, indstart)
                    irq_ind = irq_ind+1
                Enddo


                Call IWaitAll(nsirq, sirqs)


            Enddo
            DeAllocate(sirqs)
        Endif

    End Subroutine Write_Spectral_Field3D


    Subroutine Write_Checkpoint_Alt(abterms,iteration,dt,new_dt,elapsed_time)
        ! This uses the memory friendly framework
        Implicit None
        Real*8, Intent(In) :: abterms(:,:,:,:)
        Real*8, Intent(In) :: dt, new_dt, elapsed_time
        Integer, Intent(In) :: iteration
        Integer :: mp, m, offset,nl,np
        Integer :: dim2, i, offset_index, r, imi,f,ind
        Real*8, Allocatable :: myarr(:,:)
        Character*2 :: autostring
        Character*8 :: iterstring
        Character*120 :: cfile
        np = pfi%rcomm%np

        Call chktmp%construct('p1a')
        chktmp%config = 'p1a'
        !Copy the RHS into chtkmp
        Call Get_All_RHS(chktmp%p1a)
        chktmp%p1a(:,:,:,numfields+1:numfields*2) = abterms(:,:,:,1:numfields)
        !Now we want to move from p1a to s2a (rlm space)
        Call chktmp%reform()

        ! Next, each process stripes their s2a array into a true 2-D array
        dim2 = tnr*numfields*2
        Allocate(myarr(1:mode_count(my_row_rank),1:dim2))
        offset =1
        Do mp = my_mp%min, my_mp%max
                m = m_values(mp)
                nl = l_max-m+1
                ind = 1
                Do f = 1, numfields*2
                Do imi = 1, 2
                Do r = my_r%min, my_r%max
                    myarr(offset:offset+nl-1,ind) = chktmp%s2a(mp)%data(m:l_max,r,imi,f)
                    ind = ind+1
                Enddo
                Enddo
                Enddo
                offset = offset+nl
        Enddo
        Call chktmp%deconstruct('s2a')
        If (ItIsTimeForAQuickSave) Then
            write(autostring,'(i2.2)') (quicksave_num+1) !quick save number starts at 1
            checkpoint_prefix = trim(my_path)//'Checkpoints/quicksave_'//trim(autostring)
        Else
            write(iterstring,'(i8.8)') iteration
            checkpoint_prefix = trim(my_path)//'Checkpoints/'//trim(iterstring)
        Endif

        Call Write_Spectral_Field3D(myarr,1,wchar)
        Call Write_Spectral_Field3D(myarr,2,pchar)
        Call Write_Spectral_Field3D(myarr,3,tchar)
        Call Write_Spectral_Field3D(myarr,4,zchar)
        offset_index = 4
        If (magnetism) Then
            Call Write_Spectral_Field3D(myarr,5,cchar)
            Call Write_Spectral_Field3D(myarr,6,achar)
            offset_index = 6
        Endif

        Call Write_Spectral_Field3D(myarr,offset_index+1,'WAB')
        Call Write_Spectral_Field3D(myarr,offset_index+2,'PAB')
        Call Write_Spectral_Field3D(myarr,offset_index+3,'TAB')
        Call Write_Spectral_Field3D(myarr,offset_index+4,'ZAB')
        If (magnetism) Then
            Call Write_Spectral_Field3D(myarr,offset_index+5,'CAB')
            Call Write_Spectral_Field3D(myarr,offset_index+6,'AAB')
        Endif
      DeAllocate(myarr)

    If (my_column_rank .eq. 0) Then
    If (my_row_rank .eq. 0) Then
        ! row/column 0 writes out a file with the grid, etc.
        ! This file should contain everything that needs to be known
        write(iterstring,'(i8.8)') iteration
        cfile = Trim(checkpoint_prefix)//'_'//'grid_etc'
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
        Close(15)



        open(unit=15,file=Trim(my_path)//'Checkpoints/last_checkpoint',form='formatted', status='replace')
        If (ItIsTimeForAQuickSave) Then
            Write(15,'(i9.8)')-iteration
            Write(15,'(i2.2)')(quicksave_num+1)
        Else
            Write(15,'(i8.8)')iteration
        Endif
        Close(15)

        open(unit=15,file=Trim(my_path)//'Checkpoints/checkpoint_log',form='formatted', status='unknown', &
            position='Append')
        If (ItIsTimeForAQuickSave) Then
            Write(iterstring,'(i8.8)')iteration
            Write(autostring,'(i2.2)')quicksave_num+1
            Write(15,*)iterstring, ' ', autostring

        Else
            Write(15,'(i8.8)')iteration
        Endif
        Close(15)


        Endif
        Endif



    End Subroutine Write_Checkpoint_Alt

    Subroutine Read_Spectral_Field3D(arrin,ind,tag)
        ! Parallel Reading Routine For Fields in Spectral rlm configuration
        ! ISends and IReceives are used
        ! DOES NOT SUPPORT CHANGES IN PROBLEM SIZE (LMAX OR NR) CURRENTLY
        ! IN A HURRY FOR INCITE PROFILING >>>> WILL FIX THIS IN JULY

        ! NOTE - MAY NEED TO CHECK MPI GROUP FOR FILE READ AND WRITE.
        ! SHOULD PROBABLY BE GCOMM NOW -

        Implicit None
        Integer :: var_offset, offset, new_off


        Integer :: nrirq, nsirq, irq_ind, rtag, stag
        Integer, Allocatable :: rirqs(:), sirqs(:)
        Integer :: rone, rtwo, my_nrad, num_el
        Integer :: i, mp, m, lstart, nl,p,r
        Integer :: indstart(2)
        Real*8, Intent(InOut) :: arrin(1:,1:)
        Real*8, Allocatable :: arr(:,:,:), tarr(:)


        Integer, Intent(In) :: ind
        Character*8 :: iterstring
        Character*3, Intent(In) :: tag
        Character*120 :: cfile

        integer ierr, funit
        integer(kind=MPI_OFFSET_KIND) disp1,disp2
        Integer :: mstatus(MPI_STATUS_SIZE)
        !write(iterstring,'(i8.8)') iter
        cfile = trim(checkpoint_prefix)//'_'//trim(tag)


        var_offset = (ind-1)*tnr
        If (I_Will_Output) Then
            my_nrad = nradii_at_rank(my_row_rank)
            Allocate( arr(1:nlm_total,1:my_nrad,2))
            Allocate(tarr(1:nlm_total*my_nrad))

            Call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, &
                MPI_MODE_RDONLY, &
            MPI_INFO_NULL, funit, ierr)
            disp1 = my_check_disp2*8
            disp2 = (my_check_disp2+full_disp)*8

            Call MPI_FILE_SET_VIEW(funit, disp1, MPI_DOUBLE_PRECISION, &     ! Real part
              MPI_DOUBLE_PRECISION, 'native', &
              MPI_INFO_NULL, ierr)
            Call MPI_FILE_READ(funit, arr(1,1,1), buffsize2, MPI_DOUBLE_PRECISION, &
              mstatus, ierr)

            Call MPI_FILE_SET_VIEW(funit, disp2, MPI_DOUBLE_PRECISION, &     ! Imaginary part
              MPI_DOUBLE_PRECISION, 'native', &
              MPI_INFO_NULL, ierr)
            Call MPI_FILE_READ(funit, arr(1,1,2), buffsize2, MPI_DOUBLE_PRECISION, &
              mstatus, ierr)

            Call MPI_FILE_CLOSE(funit, ierr)


            !////////////////////////////////
            nrirq = Noutputs_per_Row-1
            Allocate(rirqs(1:nrirq))
            nsirq = nprow-1
            Allocate(sirqs(1:nsirq))
            Do i = 1, 2
                rone = var_offset+rstart_at_rank(my_row_rank)+(i-1)*my_r%delta
                rtwo = rone+my_nrad-1
                ! File Data into Temporary Array
                offset = 1
                Do p = 0, nprow-1
                    Do r = 1, my_nrad
                    Do mp = pfi%all_3s(p)%min, pfi%all_3s(p)%max
                        m = m_values(mp)
                        nl = l_max-m+1
                        lstart = lmstart(m)
                        !arr(lstart:lstart+nl-1,r,i) = tarr(offset:offset+nl-1)
                        tarr(offset:offset+nl-1) = arr(lstart:lstart+nl-1,r,i)
                        offset = offset+nl
                    Enddo
                    Enddo
                Enddo
                ! Post receives
                irq_ind = 1
                Do p = 0, Noutputs_per_Row-1
                        If (p .ne. my_row_rank) Then
                            rtag = p*(i+1)
                            num_el = mode_count(my_row_rank)*nradii_at_rank(p)
                            indstart(1) = 1
                            indstart(2) = var_offset+rstart_at_rank(p)+(i-1)*my_r%delta
                            Call IReceive(arrin, rirqs(irq_ind),num_el, p, rtag, pfi%rcomm, indstart)
                            irq_ind = irq_ind+1
                        Endif
                Enddo



                ! Post sends
                offset = 1
                irq_ind = 1
                Do p = 0, nprow-1
                    If (p .ne. my_row_rank) Then
                        stag = my_row_rank*(i+1)
                        num_el = mode_count(p)*my_nrad
                        Call ISend(tarr, sirqs(irq_ind),num_el, offset,p, stag, pfi%rcomm)
                        irq_ind = irq_ind+1
                    Else
                        new_off = offset
                        Do r = rone, rtwo
                            Do m = 1, mode_count(p)
                                !tarr(new_off) = arrin(m,r)
                                arrin(m,r) = tarr(new_off)
                                new_off = new_off+1
                            Enddo
                        Enddo
                    Endif
                    offset = offset+mode_count(p)*my_nrad
                Enddo
                Call IWaitAll(nsirq, sirqs)
                Call IWaitAll(nrirq, rirqs)
            Enddo
            DeAllocate(rirqs,sirqs)
            DeAllocate(tarr)
            DeAllocate(arr)
        Else
            ! This rank does not output
            ! Post sends to all output processes in my row
            nrirq = Noutputs_per_row
            Allocate(rirqs(1:nrirq))

            Do i = 1, 2
                irq_ind = 1
                Do p = 0, Noutputs_per_Row-1
                    num_el = mode_count(my_row_rank)*nradii_at_rank(p)
                    rtag = p*(i+1)
                    indstart(1) = 1
                    indstart(2) = var_offset+rstart_at_rank(p)+(i-1)*my_r%delta
                    Call IReceive(arrin, rirqs(irq_ind),num_el, p, rtag, pfi%rcomm, indstart)
                    irq_ind = irq_ind+1
                Enddo


                Call IWaitAll(nrirq, rirqs)


            Enddo
            DeAllocate(rirqs)
        Endif

    End Subroutine Read_Spectral_Field3D

    Subroutine Read_Checkpoint_Alt(fields, abterms,iteration,read_pars)
        Implicit None
        !///////////////////////////////////////
        ! DOES NOT YET PROPERLY HANDLE CHANGE IN RESOLUTION
        ! NEW MEMORY FRIENDLY VERSION FOR MIRA
        Integer, Intent(In) :: iteration, read_pars(1:2)
        Real*8, Intent(InOut) :: fields(:,:,:,:), abterms(:,:,:,:)
        Integer :: n_r_old, l_max_old, grid_type_old
        Integer :: i, ierr,  m, nl, r, f, imi, ind
        Integer :: dim2,offset, mp, offset_index
        Integer :: old_pars(5)

        Real*8, Allocatable :: old_radius(:)
        Real*8, Allocatable ::  myarr(:,:)
        Real*8 :: dt_pars(3),dt,new_dt
        Character*8 :: iterstring
        Character*2 :: autostring
        Character*120 :: cfile
        Integer :: read_hydro = 0, read_magnetism = 0
        Integer :: last_iter, last_auto

        read_hydro = read_pars(1)
        read_magnetism = read_pars(2)

        dim2 = tnr*numfields*2
        checkpoint_iter = iteration
        Write(iterstring,'(i8.8)') iteration
        If (my_rank .eq. 0) Then
            old_pars(4) = checkpoint_iter
            old_pars(5) = -1
            If (checkpoint_iter .eq. 0) Then
                open(unit=15,file=Trim(my_path)//'Checkpoints/last_checkpoint',form='formatted', status='old')
                read(15,'(i9.8)')last_iter
                If (last_iter .lt. 0) Then  !Indicates a quicksave
                    Read(15,'(i2.2)')last_auto
                    old_pars(4) = -last_iter
                    old_pars(5) = last_auto
                    Write(autostring,'(i2.2)')last_auto
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
                Else
                    !Not a quicksave
                    old_pars(4) = last_iter
                    Write(iterstring,'(i8.8)') last_iter
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
                Endif

                Close(15)
            ElseIf (checkpoint_iter .lt. 0) Then
                !User has specified a particular quicksave file
                    last_auto = -checkpoint_iter
                    old_pars(5) = -checkpoint_iter
                    Write(autostring,'(i2.2)')-checkpoint_iter
                    checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
            Else
                Write(iterstring,'(i8.8)') iteration
                checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
            Endif


            !process zero reads all the old info and broadcasts to all other ranks
            cfile = trim(checkpoint_prefix)//'_'//'grid_etc'
            open(unit=15,file=cfile,form='unformatted', status='old')
            Read(15)n_r_old
            Read(15)grid_type_old
            Read(15)l_max_old
            Read(15)dt
            Read(15)new_dt
            Allocate(old_radius(1:N_r_old))
            Read(15)(old_radius(i),i=1,N_R)
            Read(15)Checkpoint_time
            If (checkpoint_iter .lt. 0) Then
                ! We're loading a quicksave file
                ! Need to retrieve iteration from the grid_etc file because
                ! iteration specified in main_input was a low, negative number
                Read(15)Checkpoint_iter
                old_pars(4) = Checkpoint_iter
            Endif

            Close(15)
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
        Call MPI_Bcast(old_pars,5, MPI_INTEGER, 0, pfi%gcomm%comm, ierr)

        n_r_old = old_pars(1)
        grid_type_old = old_pars(2)
        l_max_old = old_pars(3)
        checkpoint_iter = old_pars(4)
        last_auto = old_pars(5)
        If (last_auto .ne. -1) Then
            !The prefix should be formed using quicksave
            Write(autostring,'(i2.2)')last_auto
            checkpoint_prefix = Trim(my_path)//'Checkpoints/quicksave_'//Trim(autostring)
        Else
            !The prefix should reflect that this is a normal checkpoint file
            Write(iterstring,'(i8.8)') checkpoint_iter
            checkpoint_prefix = Trim(my_path)//'Checkpoints/'//Trim(iterstring)
        Endif


        !///////// Later we only want to do this if the grid is actually different
        If (my_rank .ne. 0) Then
            Allocate(old_radius(1:n_r_old))
        Endif
        !Call MPI_Bcast(old_radius,n_r_old, MPI_DOUBLE_PRECISION, 0, pfi%gcomm%comm, ierr)
        !Call MPI_Bcast(dt_pars,3, MPI_DOUBLE_PRECISION, 0, pfi%gcomm%comm, ierr)

        If (my_row_rank .eq. 0) Then
            Call MPI_Bcast(old_radius,n_r_old, MPI_DOUBLE_PRECISION, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(old_radius,n_r_old, MPI_DOUBLE_PRECISION, 0, pfi%rcomm%comm, ierr)

        If (my_row_rank .eq. 0) Then
            Call MPI_Bcast(dt_pars,3, MPI_DOUBLE_PRECISION, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(dt_pars,3, MPI_DOUBLE_PRECISION, 0, pfi%rcomm%comm, ierr)


        checkpoint_dt = dt_pars(1)
        checkpoint_newdt = dt_pars(2)
        checkpoint_time = dt_pars(3)

        !////////////////////////////
        ! Next, each process stripes their s2a array into a true 2-D array
        dim2 = tnr*numfields*2
        Allocate(myarr(1:mode_count(my_row_rank),1:dim2))
        myarr(:,:) = 0.0d0
        If (read_hydro .eq. 1) Then
            Call Read_Spectral_Field3D(myarr,1,wchar)
            Call Read_Spectral_Field3D(myarr,2,pchar)
            Call Read_Spectral_Field3D(myarr,3,tchar)
            Call Read_Spectral_Field3D(myarr,4,zchar)
        Endif
        offset_index = 4
        If (magnetism) Then
            If (read_magnetism .eq. 1) Then
                Call Read_Spectral_Field3D(myarr,5,cchar)
                Call Read_Spectral_Field3D(myarr,6,achar)
            Endif
            offset_index = 6
        Endif

        If (read_hydro .eq. 1) Then
            Call Read_Spectral_Field3D(myarr,offset_index+1,'WAB')
            Call Read_Spectral_Field3D(myarr,offset_index+2,'PAB')
            Call Read_Spectral_Field3D(myarr,offset_index+3,'TAB')
            Call Read_Spectral_Field3D(myarr,offset_index+4,'ZAB')
        Endif
        If (magnetism) Then
            If (read_magnetism .eq. 1) Then
                Call Read_Spectral_Field3D(myarr,offset_index+5,'CAB')
                Call Read_Spectral_Field3D(myarr,offset_index+6,'AAB')
            Endif
        Endif


        Call chktmp%construct('s2b')
        chktmp%config = 's2b'
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


        !//////////////////////////
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
        Endif

        ! Interpolation is complete, now we just copy into the other arrays
        fields(:,:,:,1:numfields) = chktmp%p1b(:,:,:,1:numfields)
        abterms(:,:,:,1:numfields) = chktmp%p1b(:,:,:,numfields+1:numfields*2)

        Call chktmp%deconstruct('p1b')
        DeAllocate(old_radius)

    End Subroutine Read_Checkpoint_Alt

!///////////////////////////////////////////////////////////
! Routines below the dashed-double line are "staged" for deletion
!==================================================================
    Subroutine Initialize_Checkpointing_old()
        Implicit None
        Integer :: nfs(6)
        Integer :: p, np, nl, m, mp
        if (magnetism) Then
            numfields = 6
        Endif
        nfs(:) = numfields*2
        Call chktmp%init(field_count = nfs, config = 'p1a')            ! This structure hangs around through the entire run

        !////////////////////////
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
        if (my_row_rank .eq. 0) Then
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
            my_check_disp = my_check_disp*nlm_total !*2
            buffsize = nlm_total*my_r%delta
            full_disp = N_r
            full_disp = full_disp*nlm_total
        Endif
    End Subroutine Initialize_Checkpointing_old

    Subroutine Write_Field_Orig(arr,ind,tag,iter)
                Implicit None
                Integer, Intent(In) :: ind, iter
                Real*8, Intent(In) :: arr(1:,1:)
                Character*8 :: iterstring
                Character*3, Intent(In) :: tag
                Character*120 :: cfile

                integer ierr, funit , v_offset1, v_offset2
                integer(kind=MPI_OFFSET_KIND) disp1,disp2
                Integer :: mstatus(MPI_STATUS_SIZE)
             write(iterstring,'(i8.8)') iter
            cfile = Trim(my_path)//'Checkpoints/'//trim(iterstring)//'_'//trim(tag)


                 ! We have to be careful here.  Each processor does TWO writes.
                ! The first write places the real part of the field into the file.
                ! The view then changes and advances to the appropriate location of the
                ! imaginary part.  This step is crucial for checkpoints to work with
                ! Different processor configurations.
                 v_offset1 = (ind-1)*tnr+1
                v_offset2 = v_offset1+my_r%delta

                call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, &
                       MPI_MODE_WRONLY + MPI_MODE_CREATE, &
                       MPI_INFO_NULL, funit, ierr)
                if (ierr .ne. 0) Then
                    Write(6,*)'Error Opening File: ', pfi%ccomm%rank
                Endif

                disp1 = my_check_disp*8
                disp2 = (my_check_disp+full_disp)*8

                call MPI_FILE_SET_VIEW(funit, disp1, MPI_DOUBLE_PRECISION, &     ! Real part
                           MPI_DOUBLE_PRECISION, 'native', &
                           MPI_INFO_NULL, ierr)
                if (ierr .ne. 0) Then
                    Write(6,*)'Error Setting View 1: ', pfi%ccomm%rank
                Endif
                call MPI_FILE_WRITE(funit, arr(1,v_offset1), buffsize, MPI_DOUBLE_PRECISION, &
                        mstatus, ierr)
                if (ierr .ne. 0) Then
                    Write(6,*)'Error Writing 1: ', pfi%ccomm%rank
                Endif
                call MPI_FILE_SET_VIEW(funit, disp2, MPI_DOUBLE_PRECISION, &     ! Imaginary part
                           MPI_DOUBLE_PRECISION, 'native', &
                           MPI_INFO_NULL, ierr)
                if (ierr .ne. 0) Then
                    Write(6,*)'Error Setting View 2: ', pfi%ccomm%rank
                Endif
                call MPI_FILE_WRITE(funit, arr(1,v_offset2), buffsize, MPI_DOUBLE_PRECISION, &
                        mstatus, ierr)
                if (ierr .ne. 0) Then
                    Write(6,*)'Error Writing 2: ', pfi%ccomm%rank
                Endif
                call MPI_FILE_CLOSE(funit, ierr)
                if (ierr .ne. 0) Then
                    Write(6,*)'Error Closing File: ', pfi%ccomm%rank
                Endif




    End Subroutine Write_Field_Orig
End Module Checkpointing
