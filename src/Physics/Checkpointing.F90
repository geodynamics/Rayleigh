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
    Use Controls
    Use RA_MPI_BASE
    Use Chebyshev_Polynomials_Alt
    Use Parallel_IO
    Use MakeDir
    Use PDE_Coefficients, Only : Write_Equation_Coefficients_File
    Use BufferedOutput
    Use Load_Balance, Only : my_num_lm

    ! Simple Checkpointing Module
    ! Uses MPI-IO to split writing of files amongst rank zero processes from each row
    Implicit None
    Type(SphericalBuffer) :: chktmp, chktmp2, bctmp
    Integer, private :: numfields
    Integer, private :: check_err_off = 100  ! Checkpoint errors report in range 100-200.
    Integer, private :: checkpoint_version = 2
    Integer,private :: checkpoint_tag = 425
    Character*120 :: checkpoint_prefix ='nothing'
    Character*6 :: auto_fmt = '(i2.2)'  ! Format code for quicksaves
    Character*6, allocatable :: checkpoint_suffix(:)
    Integer :: checkpoint_iter = 0
    Real*8  :: checkpoint_dt, checkpoint_newdt
    Real*8  :: checkpoint_time
    Real*8, Allocatable :: boundary_mask(:,:,:,:) ! Copy of the boundary values array
    

    !///////////////////////////////////////////////////////////
    ! These variables are used for determining if it's time for a checkpoint
    Logical :: ItIsTimeForACheckpoint = .false.
    Logical :: ItIsTimeForAQuickSave = .false.
    Integer :: quicksave_num = -1
    Real*8  :: checkpoint_t0 = 0.0d0
    Real*8  :: checkpoint_elapsed = 0.0d0  ! Time elapsed since checkpoint_t0
    Real*8  :: quicksave_seconds = -1  ! Time between quick saves

    Type(Cheby_Transform_Interface) :: cheby_info

    Type(io_buffer) :: checkpoint_buffer, checkpoint_inbuffer

Contains

    !//////////////////////////////////////////////////////////
    ! Modifications related to the new Memory-Friendly Checkpointing Style
    Subroutine Initialize_Checkpointing()
        Implicit None
        Integer :: nfs(6), i, j
        Integer, Allocatable :: gpars(:,:)
        Character*4 :: sstring

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

        numfields = 4 + n_active_scalars + n_passive_scalars
        if (magnetism) then
          numfields = numfields + 2
        end if
        allocate(checkpoint_suffix(numfields*2))

        checkpoint_suffix(1:4) = (/ 'W     ', 'P     ', 'T     ', 'Z     '/)
        i = 4
        if (magnetism) then
           checkpoint_suffix(5:6) = (/'C     ', 'A     '/)
           i = 6
        end if
        do j = 1, n_active_scalars
          write(sstring,auto_fmt) (j)
          checkpoint_suffix(i+j) = 'Xa'//sstring
        end do
        i = i + n_active_scalars
        do j = 1, n_passive_scalars
          write(sstring,auto_fmt) (j)
          checkpoint_suffix(i+j) = 'Xp'//sstring
        end do
        i = i + n_passive_scalars

        do j = 1, i
          checkpoint_suffix(i+j) = trim(checkpoint_suffix(j))//'AB'
        end do

        nfs(:) = numfields*2
        Call chktmp%init(field_count = nfs, config = 'p1a') ! Persistent throughout run

        Allocate(boundary_mask(2, 2, my_num_lm, numfields))
        boundary_mask(:,:,:,:) = 0.0d0

        nfs(:) = 1
        Call bctmp%init(field_count = nfs, config = 'p1a') ! For the boundary mask

        Allocate(gpars(1:5,1:5))
        gpars(:,:) = -1
        Call checkpoint_buffer%Init(gpars, mpi_tag=checkpoint_tag, &
                spectral=.true., cache_spectral = .true., spec_comp = .true.)
        DeAllocate(gpars)
    End Subroutine Initialize_Checkpointing

    Subroutine Write_Checkpoint(abterms,iteration,dt,new_dt,elapsed_time, input_file)
        Implicit None
        Real*8, Intent(In) :: abterms(:,:,:,:), dt, new_dt, elapsed_time
        Integer, Intent(In) :: iteration
        Integer :: i, ecode, endian_tag
        Character*120 :: autostring, iterstring, cfile, checkfile
        Character*256, Intent(Out) :: input_file 
        Character*120 :: coeff_file

        endian_tag=314
        Call chktmp%construct('p1a')
        chktmp%config = 'p1a'
        !Copy the RHS into chtkmp
        Call Get_All_RHS(chktmp%p1a)
        chktmp%p1a(:,:,:,numfields+1:numfields*2) = abterms(:,:,:,1:numfields)
        !Now we want to move from p1a to s2a (rlm space)
        Call chktmp%reform()

        If (ItIsTimeForAQuickSave) Then
            write(autostring,auto_fmt) (quicksave_num+1) !quick save number starts at 1
            checkpoint_prefix = 'Checkpoints/quicksave_'//trim(autostring)
        Else
            write(iterstring,int_out_fmt) iteration
            checkpoint_prefix = 'Checkpoints/'//trim(iterstring)
        Endif
        input_file = TRIM(my_path)//trim(checkpoint_prefix)//'/main_input'
        coeff_file = TRIM(checkpoint_prefix)//'/equation_coefficients'
        If (my_rank .eq. 0) Call Make_Directory(Trim(my_path)//checkpoint_prefix,ecode)


        ! Cache and write data, index by index.
        Do i = 1, numfields*2
            checkfile = Trim(my_path)//trim(checkpoint_prefix)//'/'//trim(checkpoint_suffix(i))
            Call checkpoint_buffer%cache_data_spectral(chktmp%s2a,i)
            Call checkpoint_buffer%write_data(filename=checkfile, clear_existing = .true.)

        Enddo

        Call chktmp%deconstruct('s2a')
        
        !/////////////////////////////////////////////////////////////////////
        ! Now write out the boundary values the boundary values
        Call bctmp%construct('p1a')
        bctmp%config = 'p1a'

        Do i = 1, numfields
            bctmp%p1a(  2*(i-1)+1,:,:,1) = boundary_mask(1,:,:,i)
            bctmp%p1a(  2*(i-1)+2,:,:,1) = boundary_mask(2,:,:,i)
        Enddo

        Call bctmp%reform() ! move to s2a

        checkfile = TRIM(my_path)//TRIM(checkpoint_prefix)//'/boundary_conditions'
        Call checkpoint_buffer%cache_data_spectral(bctmp%s2a,1)
        Call checkpoint_buffer%write_data(filename=checkfile, clear_existing = .true.)

        Call bctmp%deconstruct('s2a')
        !/////////////////////////////////////////////////////////////////////

        If (my_rank .eq. 0) Then
            ! rank 0 writes out a file with the grid, etc.
            ! This file should contain everything that needs to be known
            Call Write_Equation_Coefficients_File(coeff_file)

            cfile = Trim(my_path)//trim(checkpoint_prefix)//'/'//'grid_etc'

            Open(unit=15,file=cfile,form='unformatted', status='replace', access='stream')
            Write(15)endian_tag
            Write(15)checkpoint_version
            Write(15)n_r
            Write(15)grid_type
            Write(15)l_max
            Write(15)dt
            Write(15)new_dt
            Write(15)(radius(i),i=1,N_R)
            Write(15)elapsed_time
            Write(15)iteration
            Close(15)

            Open(unit=15,file=Trim(my_path)//'Checkpoints/last_checkpoint',form='formatted', status='replace')
            If (ItIsTimeForAQuickSave) Then
                Write(15,int_minus_out_fmt)-iteration
                Write(15,'(i2.2)')(quicksave_num+1)
            Else
                Write(15,int_out_fmt)iteration
            Endif
            Close(15)

            Open(unit=15,file=Trim(my_path)//'Checkpoints/checkpoint_log',form='formatted', status='unknown', &
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

    Subroutine Read_Checkpoint(fields, abterms,iteration,read_pars)
        Implicit None
        Integer, Intent(In) :: iteration, read_pars(1:2)
        Real*8, Intent(InOut) :: fields(:,:,:,:), abterms(:,:,:,:)
        Integer :: n_r_old, l_max_old, grid_type_old
        Integer :: i, ierr, mp, lb,ub, ab_offset
        Integer :: old_pars(7 + nsubmax), fcount(3,2), version
        Integer :: last_iter, last_auto, endian_tag, funit
        Integer*8 :: found_bytes, expected_bytes, n_r_old_big, l_max_old_big
        Integer :: read_magnetism = 0, read_hydro = 0
        Integer, Allocatable :: rinds(:), gpars(:,:), read_var(:)
        Real*8 :: dt_pars(3),dt,new_dt
        Real*8, Allocatable :: old_radius(:), radius_old_loc(:)
        Real*8, Allocatable :: tempfield1(:,:,:,:), tempfield2(:,:,:,:)
        Character*120 :: autostring, dstring, iterstring, access_type
        Character*256 :: grid_file, checkfile, input_file
        Character*1 :: under_slash 
        Character*13 :: szstr
        Logical :: legacy_format, fexist
        Integer :: ncheby_old(1:nsubmax)
        Integer :: idom, n_r_loc, n_r_old_loc, irmin, irmax, irmin_old, irmax_old
        Logical :: first_time_interpolating = .True.

        read_hydro = read_pars(1)
        read_magnetism = read_pars(2)
        checkpoint_iter = iteration
        
        under_slash='/'      ! '_' underscore for legacy-format checkpoints
        access_type='stream' ! 'sequential' for legacy-format 
        old_pars(6) = 1      ! 2 for legacy-format 
        legacy_format=.false.

        allocate(read_var(numfields*2))
        read_var(:) = 0
        If (magnetism) Then
            ! hydro, magnetic, or both sets of field can be read
            ab_offset = 7 + n_active_scalars + n_passive_scalars
            read_var(1:4)   = read_hydro
            read_var(ab_offset:ab_offset+3)  = read_hydro
            read_var(5:6)   = read_magnetism
            read_var(ab_offset+4:ab_offset+5) = read_magnetism
        Else
            read_var(:) = 1
        Endif

        If (my_rank .eq. 0) Then
            old_pars(4) = checkpoint_iter
            old_pars(5) = -1
            If (checkpoint_iter .eq. 0) Then
                Open(newunit=funit,file=Trim(my_path)//'Checkpoints/last_checkpoint',form='formatted', status='old')
                Read(funit,int_minus_in_fmt)last_iter
                If (last_iter .lt. 0) Then  !Indicates a quicksave
                    Read(funit,'(i2.2)')last_auto
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

                Close(funit)
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
           
            grid_file = Trim(checkpoint_prefix)//under_slash//'grid_etc'
            Inquire(File=grid_file, Exist=fexist)
            If (.not. fexist) Then
                ! The user may be running from legacy-format checkpoints.
                ! Prior to 1.0, checkpoints did not have a dedicated directory
                legacy_format =.true.
                under_slash='_'
                access_type='sequential'
                old_pars(6)=2
                grid_file = Trim(checkpoint_prefix)//under_slash//'grid_etc'
                Inquire(File=grid_file, Exist=fexist)       
                If (fexist) Then
                    Call stdout%print(' ')
                    Call stdout%print(' ------ Legacy checkpoint format detected. ')
                    Call stdout%print(' ')
                Endif     
            Endif
            If (fexist) Then
                Open(newunit=funit,file=grid_file,form='unformatted', status='old',&
                     access=access_type, iostat=ierr)
                If (ierr .ne. 0) Then
                    Call stdout%print(' ')
                    Call stdout%print(' ------ Error: Failed to open '//TRIM(grid_file)//'.')
                    Call stdout%print(' ')
                    ierr = 2 ! grid file could not be opened
                Endif
            Else
                Call stdout%print(' ')
                Call stdout%print(' ------ Error: Could not find '//Trim(checkpoint_prefix)//'/grid_etc')
                Call stdout%print(' ------        or '//TRIM(grid_file)//'.')
                Call stdout%print(' ')
                ierr = 1 ! grid_file did not exist
            Endif

            If (ierr .eq. 0) Then
                endian_tag=0
                If (.not. legacy_format) Read(funit)endian_tag
                If (legacy_format .or. (endian_tag .eq. 314)) Then
                    If (.not. legacy_format) Then
                        Read(funit)version
                        Read(funit)n_r_old
                    Else
                        Read(funit)n_r_old
                    Endif
                    Read(funit)grid_type_old
                    Read(funit)l_max_old
                    Read(funit)dt
                    Read(funit)new_dt
                    Allocate(old_radius(1:N_r_old))
                    Read(funit)(old_radius(i),i=1,N_R_old)
                    Read(funit)Checkpoint_time
                    If (checkpoint_iter .lt. 0) Then
                        ! We're loading a quicksave file
                        ! Need to retrieve iteration from the grid_etc file because
                        ! iteration specified in main_input was a low, negative number
                        Read(funit)Checkpoint_iter
                        old_pars(4) = Checkpoint_iter
                    Endif
                    Close(funit)
                Else
                    Call stdout%print(' ')
                    Call stdout%print(' ------ Error: Endian-check failed when reading '//TRIM(grid_file)//'.')
                    Call stdout%print(' ')
                    ierr = 3
                Endif
            Endif

            ! Get the old ncheby
            If (ndomains .gt. 1) Then ! get the old ncheby from Checkpoint
                ! main_input
                input_file = TRIM(my_path)//trim(checkpoint_prefix)//'/main_input'
                Call Get_Old_Ncheby(input_file, ncheby_old)
            Else ! for ndomains = 1, user need not specify ncheby
                ncheby_old(1) = n_r_old
            Endif


            If (ierr .eq. 0) Then
                ! Verify that all checkpoint files exist and  have the correct size.
                n_r_old_big = n_r_old
                l_max_old_big = l_max_old
                expected_bytes = n_r_old_big*((l_max_old_big+1)**2 + l_max_old_big+1)*8 !TODO: Make this work for general precision
                write(6,*)'check: ', endian_tag, version, n_r_old
                Do i = 1, numfields*2
                    If (read_var(i) .eq. 1) Then
                        checkfile = trim(checkpoint_prefix)//under_slash//trim(checkpoint_suffix(i))
                        Inquire(File=checkfile, Exist=fexist,Size=found_bytes)
                        If (.not. fexist ) Then
                            Call stdout%print(' ------ Error: '//TRIM(checkfile)//' does not exist.')
                            Call stdout%print(' ')
                            ierr = 4
                        Else If (found_bytes .ne. expected_bytes) Then
                            Call stdout%print(' ------ Error: '//TRIM(checkfile)//' is the wrong size and may be corrupted.')
                            Write(szstr,'(i13)')expected_bytes
                            Call stdout%print('               Expected size (bytes): '//TRIM(szstr))
                            Write(szstr,'(i13)')found_bytes
                            Call stdout%print('               Actual size   (bytes): '//TRIM(szstr))
                            Call stdout%print(' ')
                            ierr = 5
                        Endif
                    Endif
                Enddo
            Endif

            old_pars(1) = n_r_old
            old_pars(2) = grid_type_old
            old_pars(3) = l_max_old
            old_pars(7) = ierr
            old_pars(8:8+nsubmax-1) = ncheby_old(1:nsubmax)
            dt_pars(1) = dt
            dt_pars(2) = new_dt
            dt_pars(3) = checkpoint_time

        Endif

        If (my_row_rank .eq. 0) Then    !2-D broadcast pattern
            Call MPI_Bcast(old_pars,7+nsubmax, MPI_INTEGER, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(old_pars,7+nsubmax, MPI_INTEGER, 0, pfi%rcomm%comm, ierr)

        n_r_old       = old_pars(1)
        grid_type_old = old_pars(2)
        l_max_old     = old_pars(3)
        checkpoint_iter = old_pars(4)
        last_auto = old_pars(5)

        If (old_pars(6) .eq. 2) Then
            under_slash='_'
            legacy_format=.true.
        Endif

        ncheby_old(1:nsubmax) = old_pars(8:8+nsubmax-1)

        If (old_pars(7) .ne. 0) Then
            ! Something is wrong with this checkpoint.
            If (my_rank .eq. 0) Then
                Call stdout%print(' ------ Exiting...')
                Call stdout%print(' ')
            Endif
            Call pfi%exit(check_err_off+old_pars(7))
        Endif


        !////////////////////////////////////////////////////////////////
        ! The grid file was read successfully and all primary checkpoint 
        ! files have the correct size.
        If (my_rank .eq. 0) Then
            Write(dstring,sci_note_fmt)checkpoint_time
            Call stdout%print(' ------ Checkpoint time is: '//trim(dstring))
            If (l_max_old .lt. l_max) Then
                    Call stdout%print(' ')
                    Call stdout%print('------  Checkpoint horizontal resolution is lower than current resolution.')
                    Call stdout%print('------  The old solution will be interpolated onto horizontal grid with ')
                    Call stdout%print('------  higher resolution corresponding to the new l_max.')
                    Write(szstr,'(i13)')l_max_old
                    Call stdout%print('------  Old l_max: '//TRIM(szstr))
                    Write(szstr,'(i13)')l_max
                    Call stdout%print('------  New l_max: '//TRIM(szstr))
                    Call stdout%print(' ')
            Endif
            If (l_max_old .gt. l_max) Then
                    Call stdout%print(' ')
                    Call stdout%print('------  Checkpoint horizontal resolution is higher than current resolution.')
                    Call stdout%print('------  The old SPH expansion will be truncated at the new l_max.')
                    Call stdout%print('------  This might not be a good idea.')
                    Write(szstr,'(i13)')l_max_old
                    Call stdout%print('------  Old l_max: '//TRIM(szstr))
                    Write(szstr,'(i13)')l_max
                    Call stdout%print('------  New l_max: '//TRIM(szstr))

                    Call stdout%print(' ')
            Endif

        Endif

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
            Call MPI_Bcast(old_radius,n_r_old, MPI_REAL8, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(old_radius,n_r_old, MPI_REAL8, 0, pfi%rcomm%comm, ierr)

        If (my_row_rank .eq. 0) Then
            Call MPI_Bcast(dt_pars,3, MPI_REAL8, 0, pfi%ccomm%comm, ierr)
        Endif
        Call MPI_Bcast(dt_pars,3, MPI_REAL8, 0, pfi%rcomm%comm, ierr)

        checkpoint_dt    = dt_pars(1)
        checkpoint_newdt = dt_pars(2)
        checkpoint_Time  = dt_pars(3)

        Call chktmp%construct('s2b')
        chktmp%config = 's2b'
        Do mp = my_mp%min, my_mp%max
            chktmp%s2b(mp)%data(:,:,:,:) = 0.0d0
        Enddo

        ! Initialize the input buffer  (still need to deal properly with nr_old .ne. nr)
        If (n_r_old .eq. n_r) Then
            Allocate(gpars(1:5,1:5))
            gpars(:,:) = -1
            Call checkpoint_inbuffer%Init(gpars, mpi_tag=checkpoint_tag, &
                        spectral=.true., cache_spectral = .true., spec_comp = .true., &
                        lmax_in = l_max_old, mode = 2)
            DeAllocate(gpars)
        Else
            Allocate(rinds(1:n_r_old))
            Do i = 1, n_r_old
                rinds(i) = i
            Enddo
            Allocate(gpars(1:size(rinds),1:5))
            gpars(:,:) = -1
            gpars(:,1) = rinds(:)
            gpars(1,5) = size(rinds) 
            Call checkpoint_inbuffer%Init(gpars, mpi_tag=checkpoint_tag,spectral=.true., & 
                           cache_spectral = .true., spec_comp = .true., &
                           lmax_in = l_max_old, mode = 2)
            DeAllocate(rinds)
            DeAllocate(gpars)
        Endif
        Do i = 1, numfields*2
            If (read_var(i) .eq. 1) Then
                checkfile = trim(checkpoint_prefix)//under_slash//trim(checkpoint_suffix(i))
                Call checkpoint_inbuffer%read_data(filename=checkfile)
                Call checkpoint_inbuffer%grab_data_spectral(chktmp%s2b,i)
            Endif
        Enddo

        Call chktmp%reform()    ! move to p1b

        If (.not. legacy_format) Then
            ! Load the boundary values array

            Call bctmp%construct('s2b')
            bctmp%config = 's2b'

            Do mp = my_mp%min, my_mp%max
                bctmp%s2b(mp)%data(:,:,:,:) = 0.0d0
            Enddo

            checkfile = trim(checkpoint_prefix)//'/boundary_conditions'
            Call checkpoint_inbuffer%read_data(filename=checkfile)
            Call checkpoint_inbuffer%grab_data_spectral(bctmp%s2b,1)

            Call bctmp%reform() ! move to p1b

            Do i = 1, numfields
                boundary_mask(1,:,:,i) = bctmp%p1b(  2*(i-1)+1,:,:,1)
                boundary_mask(2,:,:,i) = bctmp%p1b(  2*(i-1)+2,:,:,1)
            Enddo

            Call bctmp%deconstruct('p1b')

        Endif

        ! Loop over the domains to set the Chebyshev coefficients, 
        ! possibly interpolating in radius for each subdomain
        ! The fields are easy, since they are stored in spectral (Chebyshev) space
        ! The AB terms are a bit more complicated, since they are currently stored in physical (radius) space
        ! When we change the checkpointing format, should also store AB terms in cheby-space
        irmin_old = n_r_old
        irmax_old = n_r_old - ncheby_old(1) + 1
        irmin = n_r
        irmax = n_r - ncheby(1) + 1
        Do idom = 1, ndomains
            ! NOW, if n_r_old and grid_type_old are the same, we can copy chtkmp%p1b into abterms and
            ! fields.  Otherwise, we need to interpolate onto the current grid
            ! And by "interpolate," we mean set the lowest-order Chebyshev coefficients of 
            ! the new checkpoint using the old checkpoint, and zero out the higher-order coefficients. 
            ! For abterms, we will transform to spectral space, set the lowest-order coefficients and
            ! and transform back

            ! Loop from inner domain to outer domain
            n_r_old_loc = ncheby_old(idom)
            n_r_loc = ncheby(idom)
            
            ! Decrement the rmin/rmax indices (which go down with increasing radius/idom
            ! for each subdomain after the first
            If (idom .gt. 1) Then
                irmin_old = irmin_old - ncheby_old(idom - 1)
                irmax_old = irmax_old - ncheby_old(idom)
                irmin = irmin - ncheby(idom - 1)
                irmax = irmax - ncheby(idom)
            Endif

            ! Interpolate if necessary (in each domain separately)
            If  ((n_r_old_loc .ne. n_r_loc) .or. (grid_type_old .ne. grid_type) ) Then
                ! Interpolate
                ! We will assume the user kept the same radial domain bounds.
                ! If they  have not, this will end badly.
                If (my_rank .eq. 0) Then
                    Call stdout%print(' ')
                    If (ndomains .gt. 1) Then ! if more than one domain, specify
                                              ! which one we are in
                        Write(szstr,'(i13)')idom
                        Call stdout%print('------ In domain:     '//TRIM(szstr))
                    Endif
                    Call stdout%print('------ Radial grid has changed.')
                    Call stdout%print('------ Interpolating onto new grid.')
                    Write(szstr,'(i13)')grid_type_old
                    Call stdout%print('------ Old grid_type:     '//TRIM(szstr))
                    Write(szstr,'(i13)')grid_type
                    Call stdout%print('------ Current grid_type: '//TRIM(szstr))
                    Write(szstr,'(i13)')n_r_old_loc
                    Call stdout%print('------ Old N_R:           '//TRIM(szstr))
                    Write(szstr,'(i13)')n_r_loc
                    Call stdout%print('------ Current N_R:       '//TRIM(szstr))
                    Call stdout%print(' ')
                Endif

                If (n_r_old_loc .lt. n_r_loc) Then
                    ! The fields are OK - they are already in chebyshev space
                    ! Just set the lowest-order coefficients from the checkpoint
                    fields(irmax:irmax+n_r_old_loc-1,:,:,1:numfields) = chktmp%p1b(irmax_old:irmin_old,:,:,1:numfields)

                    ! The AB terms are stored in physical space (in radius).
                    ! They need to be transformed, lowest-order coefficients copied, and transformed back.
                    ! First, we need to initialize the old chebyshev grid.
                    Allocate(radius_old_loc(1:n_r_old_loc))
                    Call cheby_info%Init(radius_old_loc,radius(irmin),radius(irmax))  ! We assume that the domain bounds do not change
                    fcount(:,:) = numfields
                    If (first_time_interpolating) Then ! Can only initialize chktmp2 once
                        Call chktmp2%init(field_count = fcount, config = 'p1a')
                        first_time_interpolating = .False.
                    Endif
                    Call chktmp2%construct('p1a')
                    chktmp2%p1a(:,:,:,:) = 0.0d0
                    ! Allocate tempfield1, tempfield2
                    lb = lbound(chktmp%p1b,3)
                    ub = ubound(chktmp%p1b,3)
                    Allocate(tempfield1(1:n_r_old_loc,1:2,lb:ub,1))
                    Allocate(tempfield2(1:n_r_old_loc,1:2,lb:ub,1))

                    ! For each field, put the physical AB terms from the checkpoint in tempfield1
                    ! use "tospec4d" to put the spectral tempfield1 coefficients in tempfield2
                    ! copy the tempfield2 spectral coefficients into chktmp2 
                    ! (thus interpolating, since chktmp2 will have zeros for the higher-order coefficients)
                    Do i = 1, numfields 
                        tempfield1(:,:,:,:) = 0.0d0
                        tempfield2(:,:,:,:) = 0.0d0
                        tempfield1(:,:,:,1) = chktmp%p1b(irmax_old:irmin_old,:,:,numfields+i)
                        Call cheby_info%tospec4d(tempfield1,tempfield2)
                        chktmp2%p1a(irmax:irmax+n_r_old_loc-1,:,:,i) = tempfield2(:,:,:,1)
                    Enddo
                    DeAllocate(tempfield1,tempfield2)

                    ! chktmp2 now has the interpolated abterms from the checkpoint in spectral space
                    ! use "From_Spectral" to transform these fields into physical space and copy into abterms
                    Call chktmp2%construct('p1b')
                    !Normal transform(p1a,p1b)
                    Call gridcp%From_Spectral(chktmp2%p1a,chktmp2%p1b)
                    abterms(irmax:irmin,:,:,1:numfields) = chktmp2%p1b(irmax:irmin,:,:,1:numfields)

                    ! clean up
                    Call cheby_info%destroy()
                    Call chktmp2%deconstruct('p1a')
                    Call chktmp2%deconstruct('p1b')
                    Deallocate(radius_old_loc)
                Endif

            Else ! n_r_old_loc = n_r_loc
                ! no interpolation is needed, we just copy the checkpoint into fields and abterms
                fields(irmax:irmin,:,:,1:numfields) = chktmp%p1b(irmax_old:irmax_old+n_r_loc-1,:,:,1:numfields)
                abterms(irmax:irmin,:,:,1:numfields) = chktmp%p1b(irmax_old:irmax_old+n_r_loc-1,:,:,numfields+1:numfields*2)

            Endif
        Enddo

        DeAllocate(old_radius)
        Call chktmp%deconstruct('p1b')

    End Subroutine Read_Checkpoint

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

    Subroutine Store_BC_Mask(bvals)
        Implicit None
        Real*8, Intent(In) :: bvals(:,:,:,:)
        boundary_mask(:,:,:,:) = bvals(:,:,:,:)
    End Subroutine Store_BC_Mask

    Subroutine Load_BC_Mask(bvals)
        Implicit None
        Real*8, Intent(Out) :: bvals(:,:,:,:)
        bvals(:,:,:,:) = boundary_mask(:,:,:,:)
    End Subroutine Load_BC_Mask

    Subroutine Get_Old_Ncheby(input_file, ncheby_old)
        Implicit None
        Character*120, Intent(In) :: input_file
        Integer, Intent(Out) :: ncheby_old(1:nsubmax)
        Integer :: n_r_tmp, n_theta_tmp, nprow_tmp, npcol_tmp
        Real*8  :: rmin_tmp, rmax_tmp
        Integer :: npout_tmp, precise_bounds_tmp, grid_type_tmp, l_max_tmp, n_l_tmp
        Real*8 :: aspect_ratio_tmp, shell_depth_tmp
        !Integer, Parameter :: nsubmax = 256
        Integer :: ncheby_tmp(1:nsubmax)
        Real*8  :: domain_bounds_tmp(1:nsubmax+1)
        Integer :: dealias_by_tmp(1:nsubmax)
        Integer :: n_uniform_domains_tmp
        Logical :: uniform_bounds_tmp

        ! Don't overwrite the old problemsize namelist!
        n_r_tmp = n_r
        n_theta_tmp = n_theta
        nprow_tmp = nprow
        npcol_tmp = npcol
        rmin_tmp = rmin
        rmax_tmp = rmax
        npout_tmp = npout
        precise_bounds_tmp = precise_bounds
        grid_type_tmp = grid_type
        l_max_tmp = l_max
        n_l_tmp = n_l
        aspect_ratio_tmp = aspect_ratio
        shell_depth_tmp = shell_depth
        ncheby_tmp(1:nsubmax) = ncheby(1:nsubmax) 
        domain_bounds_tmp(1:nsubmax+1)  = domain_bounds(1:nsubmax+1) 
        dealias_by_tmp(1:nsubmax)  = dealias_by(1:nsubmax) 
        n_uniform_domains_tmp = n_uniform_domains
        uniform_bounds_tmp = uniform_bounds

        ! Read problemsize_namelist of input_file
        Open(unit=20, file=input_file, status="old", position="rewind")
        Read(unit=20, nml=problemsize_namelist)
        Close(20)

        ! "ncheby_old" is now in ncheby 
        ncheby_old(1:nsubmax) = ncheby(1:nsubmax)

        ! Reset the problemsize_namelist
        n_r = n_r_tmp
        n_theta = n_theta_tmp
        nprow = nprow_tmp
        npcol = npcol_tmp
        rmin = rmin_tmp
        rmax = rmax_tmp
        npout = npout_tmp
        precise_bounds = precise_bounds_tmp
        grid_type = grid_type_tmp
        l_max = l_max_tmp
        n_l = n_l_tmp
        aspect_ratio = aspect_ratio_tmp
        shell_depth = shell_depth_tmp
        ncheby(1:nsubmax) = ncheby_tmp(1:nsubmax) 
        domain_bounds(1:nsubmax+1)  = domain_bounds_tmp(1:nsubmax+1) 
        dealias_by(1:nsubmax)  = dealias_by_tmp(1:nsubmax) 
        n_uniform_domains = n_uniform_domains_tmp
        uniform_bounds = uniform_bounds_tmp
    End Subroutine Get_Old_Ncheby

End Module Checkpointing
