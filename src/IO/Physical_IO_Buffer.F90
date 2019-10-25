Module Physical_IO_Buffer
    Use RA_MPI_Base
    Use Parallel_Framework
    Use Structures
    Use ISendReceive
    Use MPI_Layer
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
        Logical :: phi_general = .false.
        Logical :: theta_general = .false. ! equatorial slices

        Logical :: r_spec = .false.
        Logical :: t_spec = .false.
        Logical :: p_spec = .false.

        ! Averaging variables
        Logical :: sum_r =.false.        ! Perform a weighted sum in r
        Logical :: sum_theta = .false.   ! Perform a weighted sum in theta
        Logical :: sum_phi = .false.     ! Straight average in phi (weighted sum not supported)
        Logical :: weighted_sum = .false. 
        Real*8, Allocatable :: theta_weights(:)  ! For summing over theta points
        Real*8, Allocatable :: radial_weights(:)
        Real*8 :: phi_weight = 0.0d0

        ! Caching variables
        Integer :: ncache_per_rec = 1
        Integer :: nrec = 1
        Integer :: rec_skip = 0 ! number of bytes to skip between cached records
        Integer :: time_index = 1

        Integer :: cascade_type  ! 1 for full cache, 2 for single cache index

        ! Time stamps:
        Integer, Allocatable :: iter(:)
        Real*8, Allocatable :: time(:)
        Logical :: write_timestamp = .false.

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
        Procedure :: timestamp
        Procedure :: collate
        Procedure :: Initialize_IO_MPI
        Procedure :: Load_Balance_IO
        Procedure :: Allocate_Receive_Buffers
        Procedure :: DeAllocate_Receive_Buffers
        Procedure :: write_data
    End Type IO_Buffer_Physical

Contains


    Subroutine Initialize_Physical_IO_Buffer(self,r_indices, theta_indices, &
                                             phi_indices,ncache, mpi_tag, &
                                             cascade, sum_weights_theta, nrec, &
                                             skip, write_timestamp, &
                                             averaging_axes)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer, Intent(In), Optional :: r_indices(1:), theta_indices(1:), phi_indices(1:)
        Integer, Intent(In), Optional :: ncache, mpi_tag, cascade, nrec, skip
        Real*8, Intent(In), Optional :: sum_weights_theta(:)
        Logical, Intent(In), Optional :: write_timestamp
        Integer, Intent(In), Optional :: averaging_axes(3)
        Integer, Allocatable :: tmp(:)
        Integer :: r, ii, ind, my_min, my_max

        !//////////////////////////////////////////////
        ! Establish my identity
        self%col_rank = pfi%ccomm%rank
        self%row_rank = pfi%rcomm%rank
        self%rank     = pfi%gcomm%rank

        If (present(averaging_axes)) Then
            If (averaging_axes(2) .eq. 1) self%sum_r = .true.
            If (averaging_axes(3) .eq. 1) self%sum_theta = .true.
            If (averaging_axes(1) .eq. 1) self%sum_phi = .true.
        Endif

        self%weighted_sum = (self%sum_r .or. self%sum_theta .or. self%sum_phi)

        self%cascade_type = 1
        if (present(cascade)) self%cascade_type = cascade

        !///////////////////////////////////////////////////////////
        ! Note how many fields/time-steps will be stored and output
        If (.not. present(ncache)) Then
            self%ncache = 1
        Else
            self%ncache = ncache
            !Write(6,*)'ncahce is: ', self%ncache
        Endif

        If (.not. present(nrec)) Then
            self%nrec = 1
        Else
            self%nrec = nrec
        Endif
        self%ncache_per_rec = self%ncache/self%nrec

        If (present(skip)) self%rec_skip = skip

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
            self%r_spec = .true.
        Else
            self%nr = pfi%n1p
        Endif

        If (present(theta_indices)) Then
            self%ntheta = size(theta_indices)
            Allocate(self%theta_global(1:self%ntheta))
            self%theta_global(:) = theta_indices(:)
            self%t_spec = .true.
            If (present(sum_weights_theta)) Then
                Allocate(self%theta_weights(1:self%ntheta))
                self%theta_weights(:) = sum_weights_theta(:)
                self%sum_theta = .true.
            Endif
        Else
            self%ntheta = pfi%n2p
        Endif

        If (present(phi_indices)) Then
            !Write(6,*)'here...phi', present(phi_indices), phi_indices
            self%nphi = size(phi_indices)
            Allocate(self%phi_local(1:self%nphi))
            self%phi_local(:) = phi_indices(:)
            self%p_spec = .true.
            !Write(6,*)'phi: ', phi_indices
        Else
            self%nphi = pfi%n3p
            If (self%sum_phi) Then
                self%nphi=1
                self%phi_weight = 1.0d0/pfi%n3p
            Endif
        Endif

        If ((.not. self%r_spec ) .and. (.not. self%t_spec) &
             .and. (.not. self%p_spec) ) self%simple = .true.

        If ( (self%r_spec ) .and. (.not. self%t_spec) &
             .and. (.not. self%p_spec) ) self%r_general = .true.

        If ( (self%t_spec ) .and. (.not. self%r_spec) &
             .and. (.not. self%p_spec) ) self%theta_general = .true.


        If ((.not. self%r_spec ) .and. (.not. self%t_spec) &
             .and. (self%p_spec) ) self%phi_general = .true.

        If (self%r_spec .and. self%t_spec &
             .and. self%p_spec ) self%general = .true.

       ! Write(6,*)self%simple, self%r_general, self%phi_general, self%general



        
        Call self%Load_Balance_IO()
        Call self%Initialize_IO_MPI()


        If (self%general) Write(6,*)'check: ', self%nr_local, self%ntheta_local, self%nphi, self%buffsize

        If (present(write_timestamp)) Then
            If ((self%output_rank) .and. (self%ocomm%rank .eq. 0) ) Then
                self%write_timestamp = write_timestamp
            Endif
            If (write_timestamp) self%rec_skip = 12
        Endif
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

        If (self%t_spec) Then

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
                n = pfi%all_2p(p)%delta
                self%ntheta_at_column(p) = n
            Enddo
        Endif

        If (self%r_spec) Then

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
        Else
            self%nr_local = pfi%all_1p(self%col_rank)%delta
            Do p = 0, pfi%npcol-1
                self%nr_out_at_row(p) = pfi%all_1p(p)%delta
            Enddo
        Endif

        !If Phi indices are specified, no extra logic is required
        !since phi is always in-process

        Do p = 0, pfi%nprow-1
            n = self%ntheta_at_column(p)
            self%npts_at_column(p) = self%nphi*self%nr_local*n
        Enddo
        self%npts = self%npts_at_column(self%row_rank)

        nout_cols = pfi%output_columns
        If (nout_cols .gt. self%nr_local) Then
            nout_cols = self%nr_local
        Endif
        self%nout_cols = nout_cols
        If (self%row_rank .lt. nout_cols) self%output_rank = .true.


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
            If (self%sum_theta) shsize = shsize/self%ntheta
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
            if (self%sum_theta) self%buffsize = self%buffsize/self%ntheta

        Endif

        !If (self%row_rank .eq. 0) Write(6,*)'c: ', self%col_rank, self%nr_out_at_column(:), self%nr_out_at_row(self%col_rank)

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
                    !If (p .eq. 0) Then
                    !    Write(6,*)'allocating: ', np, nt, nr
                    !Endif
                    Allocate(self%recv_buffers(p)%data(1:np,1:nt,1:nr,1))
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

        color = 0
        If (self%output_rank) color = 1

        Call mpi_comm_split(pfi%gcomm%comm, color, self%rank, self%ocomm%comm, ierr)
        Call mpi_comm_size(self%ocomm%comm, self%ocomm%np, ierr)
        Call mpi_comm_rank(self%ocomm%comm, self%ocomm%rank, ierr)

        self%orank = self%ocomm%rank
        If (self%output_rank) Then
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
                            1:self%nr_local,self%ncache))
        Endif
        self%cache(:,:,:,:) = 0.0d0
        If (self%write_timestamp) Then
            Allocate(self%iter(1:self%nrec))
            Allocate(self%time(1:self%nrec))
        Endif
    End Subroutine Allocate_Cache

    Subroutine Timestamp(self,ival,tval)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer, Intent(In) :: ival
        Real*8, Intent(In) :: tval
        If (self%write_timestamp) Then
            !Write(6,*)'True!', self%output_rank, self%ocomm%rank,self%nrec
            self%iter(self%time_index) = ival
            self%time(self%time_index) = tval
            self%time_index = MOD(self%time_index, self%nrec)+1
        Endif
    End Subroutine Timestamp

    Subroutine Cache_Data(self,vals)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Real*8, Intent(In) :: vals(1:,1:,1:)
        Integer :: r, t, p
        
        If (self%cascade_type .eq. 1) Then
            ! If cascade_type is 1, we may be performing a weighted sum
            If (self%weighted_sum) Then

                ! sum_phi sum_theta sum_r

                ! sum_phi and sum_theta
                
                If (self%sum_phi) Then  ! Put this one last (logical ordering)
                    Do t = 1, self%ntheta_local
                        Do r = 1, self%nr_local
                            self%cache(1,t,self%cache_index,r) = SUM(vals(:,r,t))*self%phi_weight
                        Enddo
                    Enddo          
                Endif

            Else


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
                    Do t = 1, self%ntheta_local
                        Do r = 1, self%nr_local
                            Do p = 1, self%nphi
                                self%cache(p,t,self%cache_index,r) = vals(p,self%r_local(r),t)
                            Enddo
                        Enddo
                    Enddo
                Endif

                If (self%phi_general) Then
                    Do t = 1, self%ntheta_local
                        Do r = 1, self%nr_local
                            Do p = 1, self%nphi
                                self%cache(p,t,self%cache_index,r) = vals(self%phi_local(p),r,t)
                            Enddo
                        Enddo
                    Enddo
                Endif

                If (self%theta_general) Then        ! We need to get these together...
                    Do t = 1, self%ntheta_local
                        Do r = 1, self%nr_local
                            Do p = 1, self%nphi
                                self%cache(p,t,self%cache_index,r) = vals(p,r,self%theta_local(t))
                            Enddo
                        Enddo
                    Enddo
                Endif

                If (self%general) Then
                    Do t = 1, self%ntheta_local
                        Do r = 1, self%nr_local
                            Do p = 1, self%nphi
                                self%cache(p,t,self%cache_index,r) = &
                                    & vals(self%phi_local(p),self%r_local(r),self%theta_local(t))
                            Enddo
                        Enddo
                    Enddo
                Endif

            Endif

        Else
            If (self%simple) then
                Do t = 1, self%ntheta_local
                    Do r = 1, self%nr_local
                        Do p = 1, self%nphi
                            self%cache(p,t,r,self%cache_index) = vals(p,r,t)
                        Enddo
                    Enddo
                Enddo
            Endif

            If (self%r_general) Then
                Do t = 1, self%ntheta_local
                    Do r = 1, self%nr_local
                        Do p = 1, self%nphi
                            self%cache(p,t,r,self%cache_index) = vals(p,self%r_local(r),t)
                        Enddo
                    Enddo
                Enddo
            Endif

            If (self%phi_general) Then
                !Write(6,*)'cascade type 2'
                Do t = 1, self%ntheta_local
                    Do r = 1, self%nr_local
                        Do p = 1, self%nphi
                            self%cache(p,t,r,self%cache_index) = vals(self%phi_local(p),r,t)
                        Enddo
                    Enddo
                Enddo
            Endif

        Endif

        self%cache_index = MOD(self%cache_index, self%ncache)+1
    End Subroutine Cache_Data



    Subroutine Collate(self,cache_ind)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer, Intent(In), Optional :: cache_ind
        Integer :: p, n, nn, rstart, rend,  nrirq, nsirq
        Integer, Allocatable :: rirqs(:), sirqs(:)
        Integer :: inds(4)
        Integer :: i, tstart,tend, r, t, ncache
        Real*8, Allocatable :: data_copy(:,:,:,:)
        
        Call self%Allocate_Receive_Buffers()

        ncache =1
        If (self%cascade_type .eq.1) ncache = self%ncache


        If (self%output_rank) Then
            !Post receives
            Allocate(rirqs(1:pfi%nprow-1))
            nrirq = 0
            Do p = 0, pfi%nprow-1
                If (p .ne. self%row_rank) Then
                    n = self%nrecv_from_column(p)*ncache
                    If (n .gt. 0) Then
                        !Write(6,*)self%rank, ' receiving ', n, ' from rank ', p
                        nrirq =nrirq+1
                        Call IReceive(self%recv_buffers(p)%data, rirqs(nrirq),n_elements = n, &
                                &  source= p,tag = self%tag, grp = pfi%rcomm)	            
                    Endif
                Endif
            Enddo
        Endif

        nsirq = self%nout_cols      ! maximum possible sends (point probes can alter this)
        If (self%output_rank) nsirq = nsirq-1
        Allocate(sirqs(1:nsirq))
        nn = 1
        rstart = 1
        inds(1) = 1
        inds(2) = 1
        Do p = 0, self%nout_cols-1
            if (self%cascade_type .eq. 1) Then
                inds(3) = 1
                inds(4) = rstart
            Else
                inds(3) = rstart
                inds(4) = cache_ind
            Endif
            If (p .ne. self%row_rank) Then
                n = self%nr_out_at_column(p)*ncache*self%nphi*self%ntheta_local
                !Write(6,*)self%rank, n
                If (n .gt. 0) Then
                    !Write(6,*)self%rank, ' sending ', n, ' to rank ', p
                    Call ISend(self%cache, sirqs(nn),n_elements = n, dest = p, tag = self%tag, & 
                        grp = pfi%rcomm, indstart = inds)
                    nn = nn+1
                Endif

            Else
                rend = rstart+self%nr_out-1
                If (self%cascade_type .eq. 1) Then
                    self%recv_buffers(p)%data(:,:,:,:) = self%cache(:,:,:,rstart:rend)
                Else
                    !Write(6,*)'rstart/rend', self%row_rank, rstart,rend
                    self%recv_buffers(p)%data(:,:,:,1) =  self%cache(:,:,rstart:rend,cache_ind)
                Endif
            Endif
            rstart = rstart+self%nr_out_at_column(p)
            !If(self%row_rank .ge. self%nout_cols) Then
            !    Write(6,*)'non I/O: ', maxval(self%cache)
            !Endif
        Enddo
        nsirq = nn-1  ! actual number of sends undertaken

        If (self%output_rank) Then
            !WRite(6,*)'outputting', self%rank
            !WRite(6,*)'shapes: ', nrirq, shape(rirqs)
            if (nrirq .gt. 0) Call IWaitAll(nrirq, rirqs)
            DeAllocate(rirqs)
        Endif
        If (nsirq .gt. 0) Call IWaitAll(nsirq,sirqs)
        DeAllocate(sirqs)

        
        If (self%output_rank) Then
            tstart = 1
            !write(6,*)(shape(self%collated_data))
            Do p = 0, pfi%nprow-1
                If (self%nrecv_from_column(p) .gt. 0) Then
                    tend = tstart+self%ntheta_at_column(p)-1
                    If (self%cascade_type .eq. 1) Then
                        Do i = 1, self%ncache
                            Do r =1, self%nr_out
                                self%collated_data(:,tstart:tend,r,i) = self%recv_buffers(p)%data(:,:,i,r)
                            Enddo
                        Enddo
                    Else
                        Do r =1, self%nr_out
                            self%collated_data(:,tstart:tend,r,1) = self%recv_buffers(p)%data(:,:,r,1)
                        Enddo                        
                    Endif
                    tstart = tend+1
                Endif
            Enddo            
        Endif

        Call self%deallocate_receive_buffers()

        If (self%output_rank .and. self%sum_theta) Then

            Allocate(data_copy(self%nphi,self%ntheta,self%nr_local,self%ncache))
            data_copy(:,:,:,:) = self%collated_data(:,:,:,:)
            DeAllocate(self%collated_data)
            Allocate(self%collated_data(self%nphi,1,self%nr_local,self%ncache))
            self%collated_data = 0.0d0

            Do t = 1, self%ntheta
                self%collated_data(:,1,:,:) = self%collated_data(:,1,:,:) + &
                        data_copy(:,t,:,:)*self%theta_weights(t)
            Enddo   

            DeAllocate(data_copy)

        Endif

    End Subroutine Collate

    Subroutine Write_Data(self,disp,file_unit,filename, mode)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer, Intent(In), Optional :: file_unit, mode
        Character*120, Intent(In), Optional :: filename
        Integer :: funit
        Integer :: write_mode
        Logical :: error
        Integer :: i, ierr, j, cache_start, cache_end, cache_ind
		Integer(kind=MPI_OFFSET_KIND), Intent(In), Optional :: disp
        Integer(kind=MPI_OFFSET_KIND) :: my_disp, hdisp, tdisp
		Integer :: mstatus(MPI_STATUS_SIZE)


        hdisp = 0
        If (present(disp)) hdisp = disp
      
        ! The file can be opened previously or opened by this routine

        ! The write can be conducted in various ways:
        !   1.)  Full write -- everything stored in the output ranks
        !   2.)  Single cache item write -- specified index is cascaded and written
        !   3.)  Iterative cache write -- All cache items written, but communicated one-at-a-time

        If (present(mode)) Then
            If ((mode .gt. 0) .and. (mode .le. 3)) write_mode = mode
        Else
            write_mode = 1
        Endif
        !Write(6,*)'write_mode is: ', write_mode

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

            If (self%output_rank) Then
                If (write_mode .eq. 1) Then
                    Allocate(self%collated_data(self%nphi,self%ntheta,self%nr_local,self%ncache))
                Else
                    Allocate(self%collated_data(self%nphi,self%ntheta,self%nr_local,1))
                Endif
            Endif

            If (write_mode .eq. 1) Then
                Call self%collate()
                cache_start = 1
                cache_end = self%ncache
            Endif
            If (write_mode .eq. 2) Then
                cache_start = cache_ind
                cache_end = cache_ind
            Endif
            If (write_mode .eq. 3) Then
                cache_start = 1
                cache_end = self%ncache
            Endif

            my_disp = hdisp + self%base_disp
            Do j = 1, self%nrec
                cache_start = (j-1)*self%ncache_per_rec+1
                cache_end   = j*self%ncache_per_rec
                !Write(6,*)'cache: ', cache_start, cache_end
                Do i = cache_start, cache_end
                    cache_ind = i
                    
                    If (write_mode .eq. 3 ) Call self%collate(i)
                    If (write_mode .gt. 1 ) cache_ind = 1  ! mode 2 not used, but I think this is incorrect
                    If (self%output_rank) Then

                        !Write(6,*)'info: ', my_disp, self%buffsize
                        !If (self%general) Write(6,*)my_disp        
                        !Write(6,*)my_disp
                        Call MPI_File_Seek(funit,my_disp,MPI_SEEK_SET,ierr)    
                        !IF (self%reduced) Then
                        ! Write the reduced buffer
                        ! ELSE     
                        !Write(6,*)my_disp, hdisp, funit, shape(self%collated_data)
                        !         self%collated_data(:,:,:,cache_ind)           
                        Call MPI_FILE_WRITE(funit, self%collated_data(1,1,1,cache_ind), &
                               self%buffsize, MPI_DOUBLE_PRECISION, mstatus, ierr)
                        my_disp = my_disp+self%qdisp
                        !if (ierr .ne. 0) Write(6,*)'error!: ', self%rank, self%col_rank, self%row_rank
                    Endif
                Enddo

                If (self%write_timestamp) Then  ! output_rank 0 will be the timestamp writer
                    tdisp = my_disp !  This works because output rank 0 always has my_disp 0 !-12
                    Call MPI_File_Seek(funit,tdisp,MPI_SEEK_SET,ierr)
                    !Write(6,*)'time stammping'
                    Call MPI_FILE_WRITE(funit, self%time(j), 1, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                    Call MPI_FILE_WRITE(funit, self%iter(j), 1, & 
                           MPI_INTEGER, mstatus, ierr)
                Endif


                my_disp = my_disp+self%rec_skip
            Enddo
            If (self%output_rank) Then
                !Write(6,*)'buff: ', self%orank, self%qdisp, hdisp, self%base_disp, my_disp
                DeAllocate(self%collated_data)

                If (present(filename)) Then
                    Write(6,*)'closing...'
			        Call MPI_FILE_CLOSE(funit, ierr) 
                Endif
            Endif
        Endif
    End Subroutine Write_Data


End Module Physical_IO_Buffer
