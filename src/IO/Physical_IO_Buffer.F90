Module Physical_IO_Buffer
    Use RA_MPI_Base
    Use Parallel_Framework
    Use Structures
    Use ISendReceive
    Type, Public :: IO_Buffer_Physical


        ! MPI-related variables
        Integer :: col_rank, row_rank, rank  ! column, row, and global ranks
        Integer :: orank                     ! rank in output communicator
        Integer :: tag                       ! MPI tag used in message exchange        
    
        !General variables, true for each column and row rank
        Integer :: nr            ! number of r-values in potentially subsampled array
        Integer :: ntheta        ! number of theta values in (potentially subsampled) array
        Integer :: nphi          ! number of phi values in the (potentially subsampled) array
        Integer :: ncache        ! number of fields (or fields/time indices) in the array
        Integer :: cache_index=1 ! current index within cache

        Integer :: nr_local     ! number of r-values subsampled locally (possibly 0)
        Integer :: ntheta_local     ! number of theta-values subsampled locally (possibly 0)
        Integer :: npts         ! number of points I output

        Integer, Allocatable :: r_global(:), theta_global(:)
        Integer, Allocatable :: r_local(:), theta_local(:), phi_local(:)
        
        Integer, Allocatable :: phi_ind(:)  ! phi indices in (potentially subsampled) array 
        Integer, Allocatable :: r_local_ind(:) ! r indices in subsampled array (r-my_r%min+1)

        !**Row-specific** variables related to the cascade operation
        !Everything here is true within a given row, but may different between
        !rows
        Integer, Allocatable :: ntheta_at_column(:)  ! number of theta points for each row rank
        Integer, Allocatable :: theta_inds_at_column(:,:) ! theta indices at each column
        Integer, Allocatable :: npts_at_column(:)    ! number of points (nr x ntheta x nphi) for each row rank (per field)
        !Variables only needed by output ranks
        Integer, Allocatable :: nr_out_at_column(:) ! number of radii output by row rank X in this row
        Integer, Allocatable :: nr_out_at_row(:)    ! number of radii output by column rank Y's row
        Integer, Allocatable :: nrecv_from_column(:)  ! number of points, per field, received from column X
        Integer :: nr_out = 0   ! Number of radii this process outputs


        Type(rmcontainer4d), Allocatable :: recv_buffers(:)

        ! Variables that (may) help efficiency of subsampling
        Logical :: simple   = .false.   ! no subsampling (3-D output)
        Logical :: general  = .false.   ! subsampling in all 3 dimensions (point-probe output)
        Logical :: r_general = .false.   ! subsampling in radius only (Shell Slices)
        Logical :: rp_general = .false. 

        Integer :: cascade_type  ! 1 for full cache, 2 for single cache index

        ! Buffer-specific variables
        Real*8, Allocatable :: cache(:,:,:,:)   ! storage space
        Logical :: output_rank = .false.        ! True if this rank performs I/O
        Type(communicator) :: ocomm             ! communicator for parallel output
        Integer :: nout_cols = 0                ! Number of columns that participate in I/O for this subsample

        Integer(kind=MPI_OFFSET_KIND) :: base_disp
        Integer(kind=MPI_OFFSET_KIND) :: qdisp
        Integer :: buffsize
        Integer :: nbytes = 8
        Real*8, Allocatable :: collated_data(:,:,:,:)

    Contains
        Procedure :: init => Initialize_Physical_IO_Buffer
        Procedure :: cache_data
        Procedure :: allocate_cache
        Procedure :: collate
        Procedure :: simple_cascade
        Procedure :: Initialize_IO_MPI
        Procedure :: Load_Balance_IO
        Procedure :: Allocate_Receive_Buffers
        Procedure :: DeAllocate_Receive_Buffers
        Procedure :: write_data
    End Type IO_Buffer_Physical

Contains

    Subroutine Initialize_Physical_IO_Buffer(self,r_indices, theta_indices, &
                                             phi_indices,ncache, mpi_tag)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer, Intent(In), Optional :: r_indices(1:), theta_indices(1:), phi_indices(1:)
        Integer, Intent(In), Optional :: ncache, mpi_tag
        Integer, Allocatable :: tmp(:)
        Integer :: r, ii, ind, my_min, my_max

        !//////////////////////////////////////////////
        ! Establish my identity
        self%col_rank = pfi%ccomm%rank
        self%row_rank = pfi%rcomm%rank
        self%rank     = pfi%gcomm%rank



        self%cascade_type = 1

        !///////////////////////////////////////////////////////////
        ! Note how many fields/time-steps will be stored and output
        If (.not. present(ncache)) Then
            self%ncache = 1
        Else
            self%ncache = ncache
            !Write(6,*)'ncahce is: ', self%ncache
        Endif

        !/////////////////////////////////////////////////
        ! Set a unique tag for message exchange
        If (.not. present(mpi_tag)) Then
            self%tag = 1
        Else
            self%tag = mpi_tag
        Endif

        !///////////////////////////////////////////////////////////
        ! Establish how this buffer will be subsampled, if at all

        If (present(r_indices)) Then
            self%nr = size(r_indices)
            Allocate(self%r_global(1:self%nr))
            self%r_global(:) =r_indices(:)
        Endif

        If (present(theta_indices)) Then
            self%ntheta = size(theta_indices)
            Allocate(self%theta_global(1:self%ntheta))
        Else
            self%ntheta = pfi%n2p
        Endif

        If (present(phi_indices)) Then
            self%nphi = size(phi_indices)
            Allocate(self%phi_local(1:self%nphi))
            self%phi_local(:) = phi_indices(:)
        Endif


        If (.not. present(r_indices)) Then
            If (.not. present(theta_indices)) Then
                If (.not. present(phi_indices)) Then
                    self%simple = .true.
                Endif
            Endif
        Else
            If (.not. present(theta_indices)) Then
                If (.not. present(phi_indices)) Then
                    self%r_general=.true.
                Else
                    self%rp_general = .true.
                Endif
            Else
                self%general = .true.
            Endif
        Endif



        Call self%Load_Balance_IO()
        Call self%Initialize_IO_MPI()
        Call self%Allocate_Cache()

    End Subroutine Initialize_Physical_IO_Buffer


    Subroutine Load_Balance_IO(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer :: m, n, p,  nextra, shsize, nout_cols
        Integer :: my_min, my_max
        Integer :: r, rmin, rmax
        Integer :: t, tmax, tmin
        Integer, Allocatable :: tmp(:)

        Allocate(self%ntheta_at_column(0:pfi%nprow-1))
        Allocate(self%npts_at_column(0:pfi%nprow-1))
        Allocate(self%nr_out_at_column(0:pfi%nprow-1))
        Allocate(self%nr_out_at_row(0:pfi%npcol-1))
        Allocate(self%nrecv_from_column(0:pfi%nprow-1))

        self%nr_out_at_column(:)  = 0
        self%nr_out_at_row(:)     = 0
        self%npts_at_column(:)    = 0
        self%ntheta_at_column(:)  = 0
        self%nrecv_from_column(:) = 0


        If (self%simple) Then
            self%ntheta = pfi%n2p
            self%nphi   = pfi%n3p
            self%nr       = pfi%n1p
            self%nr_local = pfi%all_1p(self%col_rank)%delta
            self%ntheta_local = pfi%all_2p(self%row_rank)%delta

            Do p = 0, pfi%nprow-1
                n = pfi%all_2p(p)%delta
                self%ntheta_at_column(p) = n
                self%npts_at_column(p) = self%nphi*self%nr_local*n
            Enddo

            
            Do p = 0, pfi%npcol-1
                self%nr_out_at_row(p) = pfi%all_1p(p)%delta
            Enddo

        Endif

        If (.not. self%simple) Then
            ! In this case, at minimum, a subset of radial indices
            ! have been specified
            self%nr_local = 0
            my_min = pfi%all_1p(self%col_rank)%min
            my_max = pfi%all_1p(self%col_rank)%max

            Allocate(tmp(1:(my_max-my_min)))

            Do p = 0, pfi%npcol-1
                rmin = pfi%all_1p(p)%min
                rmax = pfi%all_1p(p)%max
                n = 1
                Do r = 1, self%nr
                    m = self%r_global(r)
                    If ((m .ge. rmin ) .and. (m .le. rmax)) Then
                        self%nr_out_at_row(p) = n
                        If (p .eq. self%col_rank) tmp(n) = m-my_min+1
                        n = n+1
                    Endif
                Enddo
            Enddo
            self%nr_local = self%nr_out_at_row(self%col_rank)
            Allocate(self%r_local(1:self%nr_local))
            self%r_local(1:self%nr_local) = tmp(1:self%nr_local)
            DeAllocate(tmp)
            !write(6,*)'r_local_ind: ', self%r_local(1:self%nr_local), my_min, my_max


            ! Phi
            If ((.not. self%rp_general) .and. (.not. self%general) ) Then
                self%nphi   = pfi%n3p
            Endif

            ! Theta
            If (self%general) Then

                self%ntheta_local = 0
                my_min = pfi%all_2p(self%row_rank)%min
                my_max = pfi%all_2p(self%row_rank)%max

                Allocate(tmp(1:(my_max-my_min)))

                Do p = 0, pfi%nprow-1
                    tmin = pfi%all_2p(p)%min
                    tmax = pfi%all_2p(p)%max
                    n = 1
                    Do t = 1, self%ntheta
                        m = self%theta_global(t)
                        If ((m .ge. tmin ) .and. (m .le. tmax)) Then
                            self%ntheta_at_column(p) = n
                            If (p .eq. self%row_rank) tmp(n) = m-my_min+1
                            n = n+1
                        Endif
                    Enddo
                Enddo
                self%ntheta_local = self%ntheta_at_column(self%row_rank)
                Allocate(self%theta_local(1:self%ntheta_local))
                self%theta_local(1:self%ntheta_local) = tmp(1:self%ntheta_local)
                DeAllocate(tmp)

            Else

                self%ntheta_local = pfi%all_2p(self%row_rank)%delta
                Do p = 0, pfi%nprow-1
                    self%ntheta_at_column(p) = pfi%all_2p(p)%delta
                Enddo

            Endif

            Do p = 0, pfi%nprow-1
                n = self%ntheta_at_column(p)
                self%npts_at_column(p) = self%nphi*self%nr_local*n
            Enddo


        Endif



        nout_cols = pfi%output_columns
        If (nout_cols .gt. self%nr_local) Then
            nout_cols = self%nr_local
        Endif
        self%nout_cols = nout_cols
        If (self%row_rank .lt. nout_cols) self%output_rank = .true.


        self%npts = self%npts_at_column(self%row_rank)

        !Determine how many radii each rank in this row outputs

        If (self%nout_cols .gt. 0) Then
            n = self%nr_local / self%nout_cols
            nextra = Mod(self%nr_local,self%nout_cols)

            Do p = 0, self%nout_cols-1
                self%nr_out_at_column(p) = n
                If (p .lt. nextra) Then
                    self%nr_out_at_column(p) = n+1
                Endif
            Enddo
        Endif
        self%nr_out = self%nr_out_at_column(self%row_rank)
        

        If (self%output_rank) Then



            Do p = 0, pfi%nprow-1
                self%nrecv_from_column(p) = self%nr_out* &
                                            self%nphi*self%ntheta_at_column(p) 
            Enddo


            ! Determine offsets (in bytes) for MPI-IO
            shsize = self%ntheta*self%nphi*self%nbytes  ! size in bytes of a single shell
            self%qdisp = self%nr*shsize

            self%base_disp = 0
            Do p = 0, self%col_rank-1
                n = self%nr_out_at_row(p)*shsize
                self%base_disp = self%base_disp+n
            Enddo

            Do p = 0, self%row_rank-1
                n = self%nr_out_at_column(p)*shsize
                self%base_disp = self%base_disp+n
            Enddo


            self%buffsize = self%nr_out*self%nphi*self%ntheta 
        Endif


    End Subroutine Load_Balance_IO

    Subroutine Allocate_Receive_Buffers(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer :: p, np,nr,nt
        If (self%output_rank) Then
            np = self%nphi
            nr = self%nr_out
            Do p = 0, pfi%nprow-1
                nt = self%ntheta_at_column(p)
                
                If (self%cascade_type .eq. 1) Then
                    Allocate(self%recv_buffers(p)%data(1:np,1:nt,1:self%ncache,1:nr))
                Else
                    Allocate(self%recv_buffers(p)%data(1:np,1:nt,1:nr,self%ncache))
                Endif
            Enddo
        Endif
    End Subroutine Allocate_Receive_Buffers

    Subroutine DeAllocate_Receive_Buffers(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer :: p
        If (self%output_rank) Then
            Do p = 0, pfi%nprow -1
                DeAllocate(self%recv_buffers(p)%data)
            Enddo
        Endif
    End Subroutine DeAllocate_Receive_Buffers

    Subroutine Initialize_IO_MPI(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer :: p
        Integer :: nout_cols
        Integer :: color, ierr



        !Initialize Multi-column I/O
        !If (self%simple) Then
        !
        !    nout_cols = pfi%output_columns
        !    If (nout_cols .gt. pfi%my_1p%delta) Then
        !        nout_cols = pfi%my_1p%delta
        !    Endif
        !    self%nout_cols = nout_cols
        !    If (self%row_rank .lt. nout_cols) self%output_rank = .true.
        !
        !Endif

        color = 0
        If (self%output_rank) color = 1

        Call mpi_comm_split(pfi%gcomm%comm, color, self%rank, self%ocomm%comm, ierr)
        Call mpi_comm_size(self%ocomm%comm, self%ocomm%np, ierr)
        Call mpi_comm_rank(self%ocomm%comm, self%ocomm%rank, ierr)

        self%orank = self%ocomm%rank
        If (self%output_rank) Then
            !Write(6,*)'my_rank: ', pfi%gcomm%rank, ' output_rank: ', self%ocomm%rank
            Allocate(self%recv_buffers(0:pfi%nprow-1))
        Endif


    End Subroutine Initialize_IO_MPI

    Subroutine Allocate_Cache(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        If (self%cascade_type .eq. 1) Then
            Allocate(self%cache(1:self%nphi, 1:self%ntheta_local, &
                            1:self%ncache, 1:self%nr_local))
        Else
            Allocate(self%cache(1:self%nphi, 1:self%ntheta_local, &
                            1:self%nr_local,1))
        Endif
        self%cache(:,:,:,:) = 0.0d0
    End Subroutine Allocate_Cache

    Subroutine Cache_Data(self,vals)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Real*8, Intent(In) :: vals(1:,1:,1:)
        Integer :: r, t, p
        
        If (self%cascade_type .eq. 1) Then
            ! Need r_simple etc. logic here.  OK for now
            If (self%simple) Then
                Do t = 1, self%ntheta_local
                    Do r = 1, self%nr_local
                        Do p = 1, self%nphi
                            self%cache(p,t,self%cache_index,r) = vals(p,r,t)
                        Enddo
                    Enddo
                Enddo
            Endif

            If (self%r_general) Then
                !If (pfi%gcomm%rank .eq. 0) Write(6,*)'caching index: ', self%cache_index, &
                !        MOD(self%ncache,self%cache_index+1), self%ncache
                Do t = 1, self%ntheta_local
                    Do r = 1, self%nr_local
                        Do p = 1, self%nphi
                            self%cache(p,t,self%cache_index,r) = vals(p,self%r_local(r),t)
                        Enddo
                    Enddo
                Enddo
            Endif
        Else
            Do t = 1, self%ntheta_local
                Do r = 1, self%nr_local
                    Do p = 1, self%nphi
                        self%cache(p,t,r,self%cache_index) = vals(p,r,t)
                    Enddo
                Enddo
            Enddo

        Endif
        self%cache_index = MOD(self%cache_index, self%ncache)+1
    End Subroutine Cache_Data

    Subroutine Collate(self,ind)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer, Intent(In), Optional :: ind

        !If (self%simple) call self%simple_cascade()
        Call self%simple_cascade()
        


    End Subroutine Collate

    Subroutine Simple_Cascade(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer :: p, n, nn, rstart, rend,  nrirq, nsirq
        Integer, Allocatable :: rirqs(:), sirqs(:)
        Integer :: inds(4) = (/1,1,1,1/)
        Integer :: cache_ind, i, tstart,tend, r, t
        
        Call self%Allocate_Receive_Buffers()

        If (self%output_rank) Then
            !Post receives
            Allocate(rirqs(1:pfi%nprow-1))
            nrirq = 0
            Do p = 0, pfi%nprow-1
                If (p .ne. self%row_rank) Then
                    n = self%nrecv_from_column(p)*self%ncache
                    If (n .gt. 0) Then
                        nrirq =nrirq+1
                        Call IReceive(self%recv_buffers(p)%data, rirqs(nrirq),n_elements = n, &
                                &  source= p,tag = self%tag, grp = pfi%rcomm)	            
                    Endif
                Endif
            Enddo
        Endif

        nsirq = self%nout_cols
        If (self%output_rank) nsirq = nsirq-1
        Allocate(sirqs(1:nsirq))
        nn = 1
        rstart = 1
        inds(1) = 1
        inds(2) = 1
        inds(3) = 1
        Do p = 0, self%nout_cols-1
            inds(4) = rstart
            If (p .ne. self%row_rank) Then
                n = self%nr_out_at_column(p)*self%ncache*self%nphi*self%ntheta_local
                !If ((self%row_rank .eq. self%nout_cols) ) Then
                !    Write(6,*)'sending to ',p, rstart, self%cache(1,1,1,rstart), self%nr_out_at_column(p)
                !Endif
                If (n .gt. 0) Then
                    Call ISend(self%cache, sirqs(nn),n_elements = n, dest = p, tag = self%tag, & 
                        grp = pfi%rcomm, indstart = inds)
                    nn = nn+1
                Endif

            Else
                rend = rstart+self%nr_out-1
                self%recv_buffers(p)%data(:,:,:,:) = self%cache(:,:,:,rstart:rend)
            Endif
            rstart = rstart+self%nr_out_at_column(p)
            !If(self%row_rank .ge. self%nout_cols) Then
            !    Write(6,*)'non I/O: ', maxval(self%cache)
            !Endif
        Enddo


        If (self%output_rank) Then
            Call IWaitAll(nrirq, rirqs)
            DeAllocate(rirqs)
        Endif
        Call IWaitAll(nsirq,sirqs)
        DeAllocate(sirqs)


        If (self%output_rank) Then
            Allocate(self%collated_data(self%nphi,self%ntheta,self%nr_local,self%ncache))
            tstart = 1
            Do p = 0, pfi%nprow-1
                If (self%nrecv_from_column(p) .gt. 0) Then
                    tend = tstart+self%ntheta_at_column(p)-1
                    Do i = 1, self%ncache
                        Do r =1, self%nr_out
                            self%collated_data(:,tstart:tend,r,i) = self%recv_buffers(p)%data(:,:,i,r)
                        Enddo
                    Enddo
                    tstart = tend+1
                Endif
            Enddo            
        Endif

        Call self%deallocate_receive_buffers()

    End Subroutine Simple_Cascade

    Subroutine Write_Data(self,disp,file_unit,filename)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer, Intent(In), Optional :: file_unit
        Character*120, Intent(In), Optional :: filename
        Integer :: funit
        Integer :: write_mode
        Logical :: error
        Integer :: i, ierr, cache_start, cache_end, cache_ind
		Integer(kind=MPI_OFFSET_KIND), Intent(In), Optional :: disp
        Integer(kind=MPI_OFFSET_KIND) :: my_disp, hdisp
		Integer :: mstatus(MPI_STATUS_SIZE)


        hdisp = 0
        If (present(disp)) hdisp = disp

        ! The file can be opened previously or opened by this routine

        ! The write can be conducted in various ways:
        !   1.)  Full write -- everything stored in the output ranks
        !   2.)  Single cache item write -- specified index is cascaded and written
        !   3.)  Iterative cache write -- All cache items written, but communicated one-at-a-time

        write_mode = 1

        error = .false.
        ! The full write
        If (present(file_unit)) Then
            funit = file_unit
        Else
            If (present(filename)) Then
                If (self%output_rank) Then
		            Call MPI_FILE_OPEN(self%ocomm%comm, filename, & 
                           MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                           MPI_INFO_NULL, funit, ierr) 
                Endif
            Else
                error = .true.
            Endif
        Endif 

        If (.not. error) Then
            If (write_mode .eq. 1) Then
                Call self%collate()
                cache_start = 1
                cache_end = self%ncache
            Endif
            If (write_mode .eq. 2) Then
                cache_start = cache_ind
                cache_end = cache_ind
            Endif

            my_disp = hdisp + self%base_disp
            Do i = cache_start, cache_end
                cache_ind = i
                If (write_mode .eq. 3 ) Call self%collate(i)
                If (write_mode .gt. 1 ) cache_ind = 1
                If (self%output_rank) Then

                    !Write(6,*)'info: ', my_disp, self%buffsize
        
                    Call MPI_File_Seek(funit,my_disp,MPI_SEEK_SET,ierr)                    
                    Call MPI_FILE_WRITE(funit, self%collated_data(1,1,1,cache_ind), &
                           self%buffsize, MPI_DOUBLE_PRECISION, mstatus, ierr)
                    my_disp = my_disp+self%qdisp

                Endif
            Enddo

            If (self%output_rank) Then
                DeAllocate(self%collated_data)
                If (present(filename)) Then
			        Call MPI_FILE_CLOSE(funit, ierr) 
                Endif
            Endif
        Endif
    End Subroutine Write_Data


End Module Physical_IO_Buffer
