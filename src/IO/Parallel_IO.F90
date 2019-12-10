Module Parallel_IO
    Use RA_MPI_Base
    Use Parallel_Framework
    Use Structures
    Use ISendReceive
    Use MPI_Layer
    Use Fourier_Transform
    Use Legendre_Transforms, Only : Legendre_Transform
    Use Spherical_Buffer

    Type, Public :: io_buffer

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

        Integer :: nlm_local, nlm   ! Number of spherical harmonic modes locally and globally
        Integer :: nlm_out          ! Number of modes to output (same as nlm for shell_spectra)
        Integer :: lmax             ! Maximum degree ell in spherical harmonic expansion

        Integer :: nlm_in =1           ! number of modes to be read in
        Integer :: lmax_in =1          ! trunction degree l for input data

        Integer, Allocatable :: r_global(:), theta_global(:)
        Integer, Allocatable :: r_local(:), theta_local(:), phi_local(:)
        
        Integer, Allocatable :: phi_ind(:)  ! phi indices in (potentially subsampled) array 
        Integer, Allocatable :: r_local_ind(:) ! r indices in subsampled array (r-my_r%min+1)

        Integer, Allocatable :: l_values(:)   ! l-values subsampled from full spectrum
        Integer              :: n_l_samp

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
        Integer, Allocatable :: nsend_to_column(:)    ! number of points, per field, send to column X
        Integer :: nr_out = 0   ! Number of radii this process outputs

        Integer, Allocatable :: nlm_at_column(:)  ! Number of spherical harmonic modes in each column

        Type(rmcontainer4d), Allocatable :: recv_buffers(:)
        Type(SphericalBuffer) :: spectral_buffer

        ! Variables that (may) help efficiency of subsampling
        Logical :: simple   = .false.   ! no subsampling (3-D output)
        Logical :: general  = .false.   ! subsampling in all 3 dimensions (point-probe output)
        Logical :: r_general = .false.   ! subsampling in radius only (Shell Slices)
        Logical :: rp_general = .false. 
        Logical :: phi_general = .false.
        Logical :: theta_general = .false. ! equatorial slices

        Logical :: spectral = .false.
        Logical :: r_general_spectral = .false.
        Logical :: cache_spectral = .false.  ! cache in spectral configuration
        Logical :: spec_comp = .false.   ! rather than 2-d array, write out the "triangle" only

        Logical :: r_spec = .false.
        Logical :: t_spec = .false.
        Logical :: p_spec = .false.
        Logical :: l_spec = .false.


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
        Integer :: nwrites=1
        Integer :: rec_skip = 0 ! number of bytes to skip between cached records
        Integer :: time_index = 1

        Integer :: write_mode = 1  ! 1 for full cache, single cache index otherwise

        ! Time stamps:
        Integer, Allocatable :: iter(:)
        Real*8, Allocatable :: time(:)
        Logical :: write_timestamp = .false.

        ! Buffer-specific variables
        Real*8, Allocatable :: cache(:,:,:,:)   ! storage space
        Logical :: output_rank = .false.        ! True if this rank performs I/O
        Type(communicator) :: ocomm             ! communicator for parallel output
        Integer :: nout_cols = 0                ! Number of columns that participate in I/O for this subsample

        Integer(kind=MPI_OFFSET_KIND) :: base_disp, base_disp_in
        Integer(kind=MPI_OFFSET_KIND) :: qdisp, qdisp_in
        Integer(kind=MPI_OFFSET_KIND), Allocatable :: file_disp(:), file_disp_in(:)
        Integer :: buffsize   ! size to write out
        Integer :: in_buffsize ! size to read in
        Integer :: nbytes = 8
        Real*8, Pointer :: collated_data(:,:,:,:)
        Real*8, Allocatable :: buffer(:)
        Integer, Allocatable :: buffer_disp(:), buffer_indisp(:)
        Integer :: io_buffer_size, in_buffer_size
        Logical, Allocatable :: communicate(:)
        Integer, Allocatable :: ind(:)

    Contains
        Procedure :: init => Initialize_Physical_IO_Buffer
        Procedure :: Initialize_Sampled_Grid
        Procedure :: cache_data
        Procedure :: cache_data_spectral
        Procedure :: allocate_cache
        Procedure :: timestamp
        Procedure :: Initialize_IO_MPI
        Procedure :: Load_Balance_IO
        Procedure :: Allocate_Receive_Buffers
        Procedure :: DeAllocate_Receive_Buffers
        Procedure :: write_data
        Procedure :: Spectral_Prep

        Procedure :: gather_data
        Procedure :: cascade
        Procedure :: collate_physical
        Procedure :: collate_spectral
        Procedure :: set_displacements

        Procedure :: read_data
        Procedure :: grab_data_spectral
        Procedure :: distribute_data
        Procedure :: decollate_spectral
        Procedure :: despectral_prep
        Procedure :: decascade

    End Type io_buffer

Contains

    Subroutine Initialize_Physical_IO_Buffer(self,grid_pars, &
                                             ncache, mpi_tag, &
                                             cascade, averaging_weights, nrec, &
                                             skip, write_timestamp, &
                                             averaging_axes, spectral, mode, &
                                             l_values, cache_spectral, spec_comp, &
                                             lmax_in)
        Implicit None
        Class(io_buffer) :: self
        !Integer, Intent(In) :: r_indices(1:), theta_indices(1:), phi_indices(1:)
        Integer, Intent(In) :: grid_pars(1:,1:)
        Integer, Intent(In), Optional :: l_values(1:)
        Integer, Intent(In), Optional :: ncache, mpi_tag, cascade, nrec, skip
        Real*8 , Intent(In), Optional :: averaging_weights(1:,1:)
        Logical, Intent(In), Optional :: write_timestamp, spectral, cache_spectral, spec_comp
        Integer, Intent(In), Optional :: averaging_axes(3)
        Integer, Intent(In), Optional :: mode
        Integer, Intent(In), Optional :: lmax_in
        Integer :: r, ii, ind, my_min, my_max, adim(2), in_lmax
        Real*8, Allocatable :: avg_weights(:,:)

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

        ! The write can be conducted in various ways:
        !   1.)  Full write -- everything stored in the output ranks
        !   2.)  Single cache item write -- specified index is cascaded and written
        !   3.)  Iterative cache write -- All cache items written, but communicated one-at-a-time

        If (present(mode)) Then
            If ((mode .gt. 0) .and. (mode .le. 3)) self%write_mode = mode
        Else
            self%write_mode = 1
        Endif

        If (present(spectral)) self%spectral = spectral
        If (present(cache_spectral)) self%cache_spectral = cache_spectral

        !///////////////////////////////////////////////////////////
        ! Note how many fields/time-steps will be stored and output
        If (.not. present(ncache)) Then
            self%ncache = 1
        Else
            self%ncache = ncache
        Endif
        If (self%spectral) self%ncache=self%ncache*2  ! real/imaginary


        If (.not. present(nrec)) Then
            self%nrec = 1
        Else
            self%nrec = nrec
        Endif
        self%ncache_per_rec = self%ncache/self%nrec
        self%nwrites = self%nrec*self%ncache

        If (present(skip)) self%rec_skip = skip
        If (present(spec_comp)) self%spec_comp = spec_comp

        If (present(averaging_weights)) Then
            adim = shape(averaging_weights)
            Allocate(avg_weights(adim(1),adim(2)))
            avg_weights(:,:) = averaging_weights(:,:)
        Else  ! No averaging to be done.  Minimal-size array, init to -1
            Allocate(avg_weights(3,3))
            avg_weights(:,:) = -1
        Endif
        in_lmax = -1
        If (present(lmax_in)) in_lmax = lmax_in


        !/////////////////////////////////////////////////
        ! Set a unique tag for message exchange
        If (.not. present(mpi_tag)) Then
            self%tag = 1
        Else
            self%tag = mpi_tag
        Endif

        !///////////////////////////////////////////////////////////
        ! Establish how this buffer will be subsampled, if at all

        Call self%Initialize_Sampled_Grid(grid_pars, avg_weights, in_lmax)
        DeAllocate(avg_weights)

        !///////////////////////


        Call self%Load_Balance_IO()

        If (self%output_rank) Call self%Set_Displacements()
        
        Call self%Initialize_IO_MPI()

        If (present(write_timestamp)) Then
            If ((self%output_rank) .and. (self%ocomm%rank .eq. 0) ) Then
                self%write_timestamp = write_timestamp
            Endif
            If (write_timestamp) self%rec_skip = 12
        Endif
        Call self%Allocate_Cache()

    End Subroutine Initialize_Physical_IO_Buffer

    Subroutine Initialize_Sampled_Grid(self,grid_pars, averaging_weights, lmax_in)
        Implicit None
        Class(io_buffer) :: self
        Integer, Intent(In) :: grid_pars(1:,1:), lmax_in
        Real*8, Intent(In) :: averaging_weights(1:,1:)

        If (grid_pars(1,1) .ne. -1) Then
            self%nr = grid_pars(1,5)
            Allocate(self%r_global(1:self%nr))
            self%r_global(:) =grid_pars(1:self%nr,1)
            self%r_spec = .true.
        Else
            self%nr = pfi%n1p
        Endif

        If (grid_pars(1,2) .ne. -1) Then
            self%ntheta = grid_pars(2,5)
            Allocate(self%theta_global(1:self%ntheta))
            self%theta_global(:) = grid_pars(1:self%ntheta,2)
            self%t_spec = .true.

            If (averaging_weights(1,2) .gt. -1.0d-8) Then
                Allocate(self%theta_weights(1:self%ntheta))
                self%theta_weights(:) = averaging_weights(1:self%ntheta,2)
                self%sum_theta = .true.
            Endif

        Else
            self%ntheta = pfi%n2p
        Endif

        If (grid_pars(1,3) .ne. -1) Then
            self%nphi = grid_pars(3,5)
            Allocate(self%phi_local(1:self%nphi))
            self%phi_local(:) = grid_pars(1:self%nphi,3)
            self%p_spec = .true.
        Else
            self%nphi = pfi%n3p
            If (self%sum_phi) Then
                self%nphi=1
                self%phi_weight = 1.0d0/pfi%n3p
            Endif
        Endif

        If (self%rank .eq. 0) Write(6,*)'check: ', self%nr, self%ntheta, self%nphi

        If (self%spectral)  Then
            
            self%lmax = (2*self%ntheta-1)/3
            self%lmax_in = self%lmax
            self%nlm = (self%lmax+1)*(self%lmax+1) 
            self%nlm_out = (self%lmax+1)*(self%lmax+1)   

            If (lmax_in .gt. -1) Then
                self%lmax_in = lmax_in
            Endif
            self%nlm_in = (self%lmax_in+1)*(self%lmax_in+1)
            If (self%spec_comp) Then 
                self%nlm_out = ((self%lmax+1)**2 + self%lmax+1)/2
                self%nlm_in = ((self%lmax_in+1)**2 + self%lmax_in+1)/2
            Endif
            If (grid_pars(1,4) .ne. -1) Then
                self%l_spec = .true.
                self%n_l_samp = grid_pars(4,5)
                Allocate(self%l_values(1:self%n_l_samp))
                self%l_values(:) = grid_pars(1:self%n_l_samp,4)
                self%nlm_out = SUM(self%l_values)+self%n_l_samp
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

        If (self%r_general .and. self%spectral ) Then
            self%r_general_spectral = .true.
        Endif


    End Subroutine Initialize_Sampled_Grid

    Subroutine Load_Balance_IO(self)
        Implicit None
        Class(io_buffer) :: self
        Integer :: m, n, p,  nextra, shsize, nout_cols
        Integer :: my_min, my_max
        Integer :: r, rmin, rmax
        Integer :: t, tmax, tmin
        Integer :: nm, mp, mp_min, mp_max, m_value
        Integer, Allocatable :: tmp(:)
        Integer :: k, ntot, fcount(3,2), fcnt

        Allocate(self%ntheta_at_column(0:pfi%nprow-1))
        Allocate(self%npts_at_column(0:pfi%nprow-1))
        Allocate(self%nr_out_at_column(0:pfi%nprow-1))
        Allocate(self%nr_out_at_row(0:pfi%npcol-1))
        Allocate(self%nrecv_from_column(0:pfi%nprow-1))
        Allocate(self%nsend_to_column(0:pfi%nprow-1))

        self%nr_out_at_column(:)  = 0
        self%nr_out_at_row(:)     = 0
        self%npts_at_column(:)    = 0
        self%ntheta_at_column(:)  = 0
        self%nrecv_from_column(:) = 0
        self%nsend_to_column(:) = 0

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

            Allocate(tmp(1:(my_max-my_min+1)))

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
        
        If (self%spectral) Then
            Allocate(self%nlm_at_column(0:pfi%nprow-1))
            self%nlm_at_column(:) = 0
            Do p = 0, pfi%nprow-1
                nm = pfi%all_3s(p)%delta
                self%nlm_at_column(p) = (self%lmax+1)*nm
            Enddo
            self%nlm_local = self%nlm_at_column(self%row_rank)
            Do p = 0, pfi%nprow-1
                self%nrecv_from_column(p) = self%nr_out*self%nlm_at_column(p) 
                self%nsend_to_column(p) = self%nr_out_at_column(p)*self%nlm_local
            Enddo

            ! To save space, fields sampled for spectra are
            ! stored at different radii/quantity indices within
            ! a spectral buffer.  Logic for that mapping is taken
            ! care of here.
            ntot = self%ncache/2*self%nr_local
            fcnt = ntot/pfi%my_1p%delta
            k = Mod(ntot,pfi%my_1p%delta)
            If (k .gt. 0) fcnt = fcnt+1
            fcount(:,:) = fcnt
            Call self%spectral_buffer%init(field_count=fcount,config='p3b')  


        Else

            Do p = 0, pfi%nprow-1
                self%nrecv_from_column(p) = self%nr_out* &
                                            self%nphi*self%ntheta_at_column(p) 
                self%nsend_to_column(p) = self%nr_out_at_column(p)*self%ntheta_local*self%nphi
            Enddo
        Endif



    End Subroutine Load_Balance_IO

    Subroutine Set_Displacements(self)
        Implicit None
        Class(io_buffer) :: self
        Integer :: j,n, p
        Integer(kind=MPI_OFFSET_KIND) :: shsize
        ! This routine is only called by output ranks.
        ! Determine offsets (in bytes) for MPI-IO.
        shsize = self%ntheta*self%nphi*self%nbytes  ! size in bytes of a single shell

        If (self%spectral) shsize=self%nlm_out*self%nbytes

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
        If (self%spectral) Then
        self%base_disp_in = self%base_disp/self%nlm_out*self%nlm_in
        self%qdisp_in = self%qdisp/self%nlm_out*self%nlm_in
        Endif

        self%buffsize = self%nr_out*self%nphi*self%ntheta 
        if (self%spectral) self%buffsize=self%nr_out*self%nlm_out

        self%in_buffer_size = self%nlm_in*self%nr_out ! Only used in spectral mode, one field at a time

        If (self%write_mode .eq. 1) Then
            self%io_buffer_size = self%buffsize*self%nwrites
            If (self%l_spec .or. self%spec_comp) self%io_buffer_size = self%nr_out*self%nlm*self%nwrites
          
        Else
            self%io_buffer_size = self%buffsize
            
            If (self%l_spec) self%io_buffer_size = self%nr_out*self%nlm ! buffer still needs to hold a full spectrum temporarily
        Endif

        If (self%sum_theta) self%buffsize = self%buffsize/self%ntheta ! Do this AFTER setting io_buffer_size

        Allocate(self%file_disp(1:self%nwrites))
        Allocate(self%ind(1:self%nwrites))
        Allocate(self%buffer_disp(1:self%nwrites))
        Allocate(self%file_disp_in(1:self%nwrites))

        Do j = 1, self%nwrites
            self%file_disp(j) = self%base_disp + self%qdisp*(j-1) + self%rec_skip*((j-1)/self%ncache)
            self%file_disp_in(j) = self%base_disp_in + self%qdisp_in*(j-1) ! no timestamps on input
            self%ind(j) = j
            self%buffer_disp(j) = 1+(j-1)*self%buffsize
        Enddo
        If (self%write_mode .ne. 1) self%buffer_disp(:) = 1  
        
    End Subroutine Set_Displacements

    Subroutine Allocate_Receive_Buffers(self)
        Implicit None
        Class(io_buffer) :: self
        Integer :: p, np,nr,nt,lp1, mp_min,mp_max
        If (self%output_rank) Then
            If (self%spectral) Then
                nr = self%nr_out
                lp1 = self%lmax+1
                Do p = 0, pfi%nprow-1
                    mp_min = pfi%all_3s(p)%min
                    mp_max = pfi%all_3s(p)%max
                    If (self%write_mode .eq. 1) Then
                        Allocate(self%recv_buffers(p)%data(1:lp1,mp_min:mp_max,1:self%ncache,1:nr))
                    Else
                        Allocate(self%recv_buffers(p)%data(1:lp1,mp_min:mp_max,1:nr,1))
                    Endif
                    self%recv_buffers(p)%data(:,:,:,:) = 0.0d0
                Enddo
            Else
                np = self%nphi
                nr = self%nr_out
                Do p = 0, pfi%nprow-1

                    nt = self%ntheta_at_column(p)
                    
                    If (self%write_mode .eq. 1) Then
                        Allocate(self%recv_buffers(p)%data(1:np,1:nt,1:self%ncache,1:nr))
                    Else
                        Allocate(self%recv_buffers(p)%data(1:np,1:nt,1:nr,1))
                    Endif
                    self%recv_buffers(p)%data(:,:,:,:) = 0.0d0
                Enddo
            Endif
        Endif
    End Subroutine Allocate_Receive_Buffers

    Subroutine DeAllocate_Receive_Buffers(self)
        Implicit None
        Class(io_buffer) :: self
        Integer :: p
        If (self%output_rank) Then
            Do p = 0, pfi%nprow -1
                DeAllocate(self%recv_buffers(p)%data)
            Enddo
        Endif
    End Subroutine DeAllocate_Receive_Buffers

    Subroutine Initialize_IO_MPI(self)
        Implicit None
        Class(io_buffer) :: self
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
        Class(io_buffer) :: self

        If (.not. self%spectral) Then
            If (self%write_mode .eq. 1) Then
                Allocate(self%cache(1:self%nphi, 1:self%ntheta_local, &
                                1:self%ncache, 1:self%nr_local))
            Else
                Allocate(self%cache(1:self%nphi, 1:self%ntheta_local, &
                                1:self%nr_local,self%ncache))
            Endif
            self%cache(:,:,:,:) = 0.0d0

        Else
            If (self%cache_spectral) Then
                Call self%spectral_buffer%construct('s2b')
            Else    
                Call self%spectral_buffer%construct('p3b')
            Endif
        Endif
        If (self%write_timestamp) Then
            Allocate(self%iter(1:self%nrec))
            Allocate(self%time(1:self%nrec))
        Endif
    End Subroutine Allocate_Cache

    Subroutine Timestamp(self,ival,tval)
        Implicit None
        Class(io_buffer) :: self
        Integer, Intent(In) :: ival
        Real*8, Intent(In) :: tval
        If (self%write_timestamp) Then
            self%iter(self%time_index) = ival
            self%time(self%time_index) = tval
            self%time_index = MOD(self%time_index, self%nrec)+1
        Endif
    End Subroutine Timestamp

    Subroutine Cache_Data_Spectral(self, spec_vals,in_cache)
        Implicit None
        Class(io_buffer) :: self
        Type(rmcontainer4D), Intent(In) :: spec_vals(1:)
        Integer, Intent(In) :: in_cache
        Integer :: my_mp_min, my_mp_max, mp

        my_mp_min = pfi%all_3s(self%row_rank)%min
        my_mp_max = pfi%all_3s(self%row_rank)%max

        ! This bit was added for checkpoints
        ! It is assumed that no sampling takes place in this mode
        ! No caching either (1 at a time)
        Do mp = my_mp_min,my_mp_max                          
            self%spectral_buffer%s2b(mp)%data(:,:,:,1) = &
                spec_vals(mp-my_mp_min+1)%data(:,:,:,in_cache)
        Enddo
        
    End Subroutine Cache_Data_Spectral

    Subroutine Grab_Data_Spectral(self, spec_vals,in_cache)
        Implicit None
        Class(io_buffer) :: self
        Type(rmcontainer4D), Intent(InOut) :: spec_vals(1:)
        Integer, Intent(In) :: in_cache
        Integer :: my_mp_min, my_mp_max, mp

        my_mp_min = pfi%all_3s(self%row_rank)%min
        my_mp_max = pfi%all_3s(self%row_rank)%max

        ! This is also used for checkpoints.
        ! It is assumed that no sampling takes place in this mode.
        ! No caching either (1 at a time).
        If (self%nr_local .gt. 0) Then
            Do mp = my_mp_min,my_mp_max                          
                spec_vals(mp-my_mp_min+1)%data(:,:,:,in_cache) = &
                    self%spectral_buffer%s2b(mp)%data(:,:,:,1)
                    
            Enddo
        Endif
        
    End Subroutine Grab_Data_Spectral

    Subroutine Cache_Data(self,vals, spec_vals, in_cache)
        Implicit None
        Class(io_buffer) :: self
        Real*8, Intent(In) :: vals(1:,1:,1:)
        Type(rmcontainer4D), Intent(In), Optional :: spec_vals(1:)
        Integer, Intent(In), Optional :: in_cache
        Integer :: r, t, p,tind
        Integer :: counter, field_ind, rind
        Integer :: my_mp_min, my_mp_max, mp

        If (self%r_general_spectral) Then
   
            Do r = 1, self%nr_local
                counter = (self%cache_index-1)*self%nr_local +r-1
                field_ind = counter/pfi%my_1p%delta+1
                rind = MOD(counter, pfi%my_1p%delta)+pfi%my_1p%min

                Do t = 1, self%ntheta_local
                    tind = pfi%my_2p%min+t-1
                    Do p =1, self%nphi
                        self%spectral_buffer%p3b(p,rind,tind,field_ind) = vals(p,self%r_local(r),t)
                    Enddo
                Enddo

            Enddo

        Endif
        
        If (.not. self%spectral) Then
        If (self%write_mode .eq. 1) Then
            ! If write_mode is 1, we may be performing a weighted sum
            If (self%weighted_sum) Then

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
                Do t = 1, self%ntheta_local
                    Do r = 1, self%nr_local
                        Do p = 1, self%nphi
                            self%cache(p,t,r,self%cache_index) = vals(self%phi_local(p),r,t)
                        Enddo
                    Enddo
                Enddo
            Endif

        Endif
        Endif

        If (self%spectral) Then
            ! Adjust logic for real/imaginary
            self%cache_index = MOD(self%cache_index, self%ncache/2)+1    
        Else
            self%cache_index = MOD(self%cache_index, self%ncache)+1
        Endif
    End Subroutine Cache_Data

    Subroutine Spectral_Prep(self)
        Implicit None
        Class(io_buffer) :: self
        Integer :: p, mp, f, counter, field_ind, m, my_mp_min, my_mp_max, nf
        Integer :: nq, r, rind, lmax, i1, i2, my_nm, mstore
        ! This routine transforms phi-theta to ell-m space.
        ! It then puts the spectral data into the cache buffer
        ! in a format that is compatible with the collate and write routines.
        lmax = self%lmax
        my_mp_min = pfi%all_3s(self%row_rank)%min
        my_mp_max = pfi%all_3s(self%row_rank)%max
        my_nm     = my_mp_max-my_mp_min+1
        If (self%nr_local .gt. 0) Then
            If (.not. self%cache_spectral) Then
                !(1) Transform/transpose the phi-theta buffer to Ylm space
                Call FFT_To_Spectral(self%spectral_buffer%p3b, rsc = .true.)
                self%spectral_buffer%config ='p3b'
                Call self%spectral_buffer%reform()
                Call self%spectral_buffer%construct('s2b')
                Call Legendre_Transform(self%spectral_buffer%p2b,self%spectral_buffer%s2b)
                Call self%spectral_buffer%deconstruct('p2b')
            Endif
            !(2) Stripe the spectral data into the self%cache buffer
            !    While that buffer has nphi and ntheta dimensions, we 
            !    alias those to n_lm and nr
            If (self%write_mode .eq. 1) Then
                Allocate(self%cache(1:self%lmax+1, 1:my_nm, &
                    self%ncache, 1:self%nr_local))     
            Else
                Allocate(self%cache(1:self%lmax+1, 1:my_nm, &
                    1:self%nr_local, self%ncache))  
            Endif

            self%cache = 0.0d0

            nf = self%spectral_buffer%nf2b
            nq = self%ncache/2
            i2 = lmax+1

            Do p = 1, 2  ! Real and imaginary parts TODO:  I don't think this loop is here anymore...
                Do mp = my_mp_min,my_mp_max
                    m = pfi%inds_3s(mp)
                    counter = 0
                    i1 = m+1
                    mstore = mp-my_mp_min+1
                    Do f = 1, nq

                        field_ind = counter/pfi%my_1p%delta+1
                        Do r = 1, self%nr_local   
                         
                            rind = MOD(counter, pfi%my_1p%delta)+pfi%my_1p%min
                            If (self%write_mode .eq. 1) Then                            
                                self%cache(i1:i2,mstore,f,r) = &
                                    self%spectral_buffer%s2b(mp)%data(m:lmax,rind,1,field_ind)

                                self%cache(i1:i2,mstore,f+nq,r) = &
                                    self%spectral_buffer%s2b(mp)%data(m:lmax,rind,2,field_ind)
                            Else
                                self%cache(i1:i2,mstore,r,f) = &
                                    self%spectral_buffer%s2b(mp)%data(m:lmax,rind,1,field_ind)

                                self%cache(i1:i2,mstore,r, f+nq) = &
                                    self%spectral_buffer%s2b(mp)%data(m:lmax,rind,2,field_ind)

                            Endif
                            counter = counter+1
                        Enddo
                    Enddo
                Enddo
            Enddo

            If (.not. self%cache_spectral) Then
                Call self%spectral_buffer%deconstruct('s2b')
                self%spectral_buffer%config ='p3b'
                Call self%spectral_buffer%construct('p3b')  ! TODO:  Do we want p3b to effectively be persistent?
            Endif
            !~~~~~

        Endif
    End Subroutine Spectral_Prep

    Subroutine Gather_Data(self, cache_ind)
        Implicit None
        Class(io_buffer) :: self
        Integer, Intent(In) :: cache_ind
    
        If (self%write_mode .eq. 1) Then
            If (self%spectral) Call self%Spectral_Prep()
        Else
            If ( (cache_ind .eq. 1) .and. (self%spectral) ) Call self%spectral_prep()
        Endif
    
        Call self%cascade(cache_ind)
        
        If (self%spectral) Then
            Call self%collate_spectral(cache_ind)
        Else
            Call self%collate_physical(cache_ind)
        Endif
        
    End Subroutine Gather_Data

    Subroutine Distribute_Data(self, cache_ind)
        Implicit None
        Class(io_buffer) :: self
        Integer, Intent(In) :: cache_ind

        !NOTE:  Assuming we are spectral for distribution

        Call self%decollate_spectral()
        Call self%decascade()
        Call self%deSpectral_Prep(cache_ind)

    End Subroutine Distribute_Data

    Subroutine Collate_Spectral(self,cache_ind)
        Implicit None
        Class(io_buffer), Target :: self
        Integer, Intent(In) :: cache_ind
        Integer :: m, r, i, p, mp_min, mp_max
        Integer :: cend, ncache, mp, k, l, j, lval, lmax, ind, ind2
        Logical :: free_mem
        Real*8, Allocatable :: data_copy(:,:,:,:)

        free_mem = .false.
        If (self%write_mode .eq. 1) free_mem = .true.
        If (cache_ind .eq. self%nwrites) free_mem = .true.
        If (self%nr_local .eq. 0) free_mem=.false.

        If (free_mem) DeAllocate(self%cache)
        ncache =1
        If (self%write_mode .eq. 1) ncache = self%ncache
        
        If (self%output_rank) Then
            cend = 1
            If (self%write_mode .eq. 1) cend = self%ncache
            
            self%collated_data(1:self%lmax+1,1:self%lmax+1,1:self%nr_out,1:cend) => &
                self%buffer(1:self%io_buffer_size)

            Do p = 0, pfi%nprow-1    ! TODO:   we should really convert shell_spectra to compressed format
                mp_min = pfi%all_3s(p)%min
                mp_max = pfi%all_3s(p)%max
                If (self%nrecv_from_column(p) .gt. 0) Then
                    If (self%write_mode .eq. 1) Then
                        Do i = 1, self%ncache
                            Do r =1, self%nr_out
                                Do mp = mp_min, mp_max
                                    m = pfi%inds_3s(mp)
                                    self%collated_data(:,m+1,r,i) = self%recv_buffers(p)%data(:,mp,i,r)
                                Enddo           
                            Enddo
                        Enddo
                    Else
                        Do r =1, self%nr_out
                            Do mp = mp_min, mp_max
                                m = pfi%inds_3s(mp)
                                self%collated_data(:,m+1,r,1) = self%recv_buffers(p)%data(:,mp,r,1)
                            Enddo
                        Enddo                        
                    Endif

                Endif
            Enddo    
            If (free_mem) Call self%deallocate_receive_buffers()

            If (self%l_spec ) Then
                Allocate(data_copy(1:self%lmax+1,1:self%lmax+1,1:self%nr_out,1:cend))
                ! Restripe to m-l (i,j) from l-m
                Do k = 1, cend
                Do r = 1, self%nr_out
                Do j = 1, self%lmax+1
                Do i = 1, self%lmax+1
                    data_copy(i,j,r,k) = self%collated_data(j,i,r,k)
                Enddo
                Enddo
                Enddo
                Enddo

                ! next, repoint collated_data
                self%collated_data(1:self%nlm_out,1:1,1:self%nr_out,1:cend) => &
                    self%buffer(1:self%nlm_out*self%nr_out*cend)
                If (self%l_spec) Then   
                    ! Grab the desired l values
                    Do k = 1, cend
                    Do r = 1, self%nr_out
                        ind = 1
                        Do l = 1, self%n_l_samp
                            lval = self%l_values(l)
                            self%collated_data(ind:ind+lval,1,r,k) = &
                                data_copy(1:lval+1,lval+1,r,k)
                            ind = ind+lval+1
                        Enddo
                
                    Enddo
                    Enddo
                Endif
                DeAllocate(data_copy)
            Endif
            If (self%spec_comp) Then

                Allocate(data_copy(1:self%lmax+1,1:self%lmax+1,1:self%nr_out,1:cend))
                data_copy(:,:,:,:) = self%collated_data(:,:,:,:)

                self%collated_data(1:self%nlm_out,1:1,1:self%nr_out,1:cend) => &
                    self%buffer(1:self%nlm_out*self%nr_out*cend)

                lmax=self%lmax

                ! Grab the desired l values
                Do k = 1, cend
                Do r = 1, self%nr_out
                    ind = 1
                    Do m = 1, lmax+1
                        ind2 = ind + lmax+1-m
                        self%collated_data(ind:ind2,1,r,k) = &
                            data_copy(m:lmax+1,m,r,k)
                        ind = ind2+1
                    Enddo
            
                Enddo
                Enddo
                DeAllocate(data_copy)
            Endif

        Endif

    End Subroutine Collate_Spectral

    Subroutine DeCollate_Spectral(self)
        Implicit None
        Class(io_buffer), Target :: self
        Integer :: m, r, i, p, mp_min, mp_max,mm, lmax, maxl
        Integer :: mp, k, l, j,  ind1, ind2, loff, lmax_in, dind
        ! restripe i/o buffer into send/recv buffers for each rank
        lmax = self%lmax
        lmax_in = self%lmax_in
        maxl = lmax_in
        if (lmax_in .gt. lmax) maxl = lmax

        If (self%output_rank) Then
            Call self%Allocate_Receive_Buffers()
            self%collated_data(1:self%nlm_in,1:1,1:self%nr_out,1:1) => &
                self%buffer(1:self%nlm_in*self%nr_out)
            Do p = 0, pfi%nprow-1
                If (self%nrecv_from_column(p) .gt. 0) Then
                    mp_min = pfi%all_3s(p)%min
                    mp_max = pfi%all_3s(p)%max
                    Do r =1, self%nr_out
                        Do mp = mp_min, mp_max
                            m = pfi%inds_3s(mp)
                            loff = 0
                            Do mm = 1, m
                                loff = loff+(lmax_in-mm+2) 
                            Enddo

                            !ind1 = (r-1)*self%nlm_in+loff+1
                            ind1 = loff+1
                            ind2 = ind1+(maxl-m)
                            dind = ind2-ind1+1
                            !WRite(6,*)'m value: ', m
                            self%recv_buffers(p)%data(m+1:maxl+1,mp,r,1) = &
                                self%collated_data(ind1:ind2,1,r,1)
                        Enddo
                    Enddo                        
                Endif
            Enddo    

        Endif

    End Subroutine DeCollate_Spectral

    Subroutine DeCascade(self)
        Implicit None
        Class(io_buffer) :: self
        Integer :: p, n, nn, rstart, rend,  nrirq, nsirq
        Integer, Allocatable :: rirqs(:), sirqs(:)
        Integer :: inds(4)
        Integer :: i, r, t, ncache, my_nm

        ! Receive 1 cache ind at a time
        ! Data is presumed to be single field, real/imaginary
        ! So this routine will be called twice
        ! We Reuse variables but flip send/receive logic from cascade here.

        my_nm = pfi%all_3s(self%row_rank)%delta
        if (self%nr_local .gt. 0) Allocate(self%cache(1:self%lmax+1, 1:my_nm, &
            1:self%nr_local, 1))  

        If (self%output_rank) Then
            !Post SENDS
            Allocate(rirqs(1:pfi%nprow-1))
            nrirq = 0
            Do p = 0, pfi%nprow-1
                If (p .ne. self%row_rank) Then
                    n = self%nrecv_from_column(p)
                    If (n .gt. 0) Then
                        nrirq =nrirq+1
                        Call ISend(self%recv_buffers(p)%data, rirqs(nrirq),n_elements = n, &
                                &  dest= p,tag = self%tag, grp = pfi%rcomm)	            
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
        inds(4) = 1
        Do p = 0, self%nout_cols-1

            inds(3) = rstart

            If (p .ne. self%row_rank) Then

                n = self%nsend_to_column(p)
                If (n .gt. 0) Then

                    Call IReceive(self%cache, sirqs(nn),n_elements = n, source = p, tag = self%tag, & 
                        grp = pfi%rcomm, indstart = inds)
                    nn = nn+1
                Endif

            Else
                If (self%output_rank) Then
                    rend = rstart+self%nr_out-1

                    self%cache(:,:,rstart:rend,1) = self%recv_buffers(p)%data(:,:,:,1)

                Endif
            Endif
            rstart = rstart+self%nr_out_at_column(p)
        Enddo

        nsirq = nn-1  ! actual number of receives undertaken

        If (self%output_rank) Then
            if (nrirq .gt. 0) Call IWaitAll(nrirq, rirqs)
            DeAllocate(rirqs)
            Call self%Deallocate_receive_buffers()
        Endif
        If (nsirq .gt. 0) Call IWaitAll(nsirq,sirqs)
        DeAllocate(sirqs)

    End Subroutine DeCascade

    Subroutine DeSpectral_Prep(self, cache_ind)
        Implicit None
        Class(io_buffer) :: self
        Integer, Intent(In) :: cache_ind
        Integer :: p, mp, f,m, my_mp_min, my_mp_max, myrmin
        Integer :: r, lmax, i1, i2, my_nm, mstore
        ! This routine takes the cache and stripes it back into
        ! the s2b format.  Cache_ind is assumed to be 1 or 2 (real/imag)
        lmax = self%lmax
        my_mp_min = pfi%all_3s(self%row_rank)%min
        my_mp_max = pfi%all_3s(self%row_rank)%max
        my_nm     = my_mp_max-my_mp_min+1
        If (self%nr_local .gt. 0) Then

            i2 = lmax+1
            myrmin = pfi%my_1p%min

            Do mp = my_mp_min,my_mp_max
                m = pfi%inds_3s(mp)
                i1 = m+1
                mstore = mp-my_mp_min+1
                Do r = 1, self%nr_local   
                        self%spectral_buffer%s2b(mp)%data(m:lmax,myrmin+r-1,cache_ind,1) = &
                            self%cache(i1:i2,mstore,r,1) 
                Enddo
            Enddo
            DeAllocate(self%cache)
        Endif
    End Subroutine DeSpectral_Prep

    Subroutine Cascade(self,cache_ind)
        Implicit None
        Class(io_buffer) :: self
        Integer, Intent(In) :: cache_ind
        Integer :: p, n, nn, rstart, rend,  nrirq, nsirq
        Integer, Allocatable :: rirqs(:), sirqs(:)
        Integer :: inds(4)
        Integer :: i, r, t, ncache

        If ( (self%write_mode .eq. 1) .or. (cache_ind .eq. 1) ) Call self%Allocate_Receive_Buffers()

        ncache =1
        If (self%write_mode .eq.1) ncache = self%ncache

        If (self%output_rank) Then
            !Post receives
            Allocate(rirqs(1:pfi%nprow-1))
            nrirq = 0
            Do p = 0, pfi%nprow-1
                If (p .ne. self%row_rank) Then
                    n = self%nrecv_from_column(p)*ncache
                    If (n .gt. 0) Then
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
            If (self%write_mode .eq. 1) Then
                inds(3) = 1
                inds(4) = rstart
            Else
                inds(3) = rstart
                inds(4) = cache_ind
            Endif
            If (p .ne. self%row_rank) Then

                n = self%nsend_to_column(p)*ncache
                If (n .gt. 0) Then

                    Call ISend(self%cache, sirqs(nn),n_elements = n, dest = p, tag = self%tag, & 
                        grp = pfi%rcomm, indstart = inds)
                    nn = nn+1
                Endif

            Else
                If (self%output_rank) Then
                rend = rstart+self%nr_out-1
                If (self%write_mode .eq. 1) Then
                    self%recv_buffers(p)%data(:,:,:,:) = self%cache(:,:,:,rstart:rend)
                Else
                    self%recv_buffers(p)%data(:,:,:,1) =  self%cache(:,:,rstart:rend,cache_ind)
                Endif
                Endif
            Endif
            rstart = rstart+self%nr_out_at_column(p)
        Enddo

        nsirq = nn-1  ! actual number of sends undertaken

        If (self%output_rank) Then
            if (nrirq .gt. 0) Call IWaitAll(nrirq, rirqs)
            DeAllocate(rirqs)
        Endif

        If (nsirq .gt. 0) Call IWaitAll(nsirq,sirqs)

        DeAllocate(sirqs)

    End Subroutine Cascade

    Subroutine Collate_Physical(self,cache_ind)
        Implicit None
        Class(io_buffer), Target :: self
        Integer, Intent(In), Optional :: cache_ind
        Integer :: p, n, nn, rstart, rend,  nrirq, nsirq
        Integer, Allocatable :: rirqs(:), sirqs(:)
        Integer :: inds(4)
        Integer :: i, tstart,tend, r, t, ncache
        Real*8, Allocatable :: data_copy(:,:,:,:)
        Integer :: cend
        Logical :: free_mem
        ncache =1
        If (self%write_mode .eq. 1) ncache = self%ncache

        free_mem = .false.
        If (self%write_mode .eq. 1) free_mem = .true.
        If (cache_ind .eq. self%nwrites) free_mem = .true.
        
        If (self%output_rank) Then
            cend = 1
            If (self%write_mode .eq. 1) cend = self%ncache

            self%collated_data(1:self%nphi,1:self%ntheta,1:self%nr_out,1:cend) => &
                self%buffer(1:self%io_buffer_size)

            tstart = 1
            Do p = 0, pfi%nprow-1
                If (self%nrecv_from_column(p) .gt. 0) Then
                    tend = tstart+self%ntheta_at_column(p)-1
                    If (self%write_mode .eq. 1) Then
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
            If (free_mem) Call self%deallocate_receive_buffers()
        Endif

        If (self%output_rank .and. self%sum_theta) Then

            Allocate(data_copy(self%nphi,self%ntheta,self%nr_out,self%ncache))
            data_copy(:,:,:,:) = self%collated_data(:,:,:,:)
 
            self%collated_data(1:self%nphi,1:1,1:self%nr_out,1:self%ncache) => &
                self%buffer(1:self%io_buffer_size/self%ntheta)
            
            self%collated_data(:,:,:,:) = 0.0d0

            Do t = 1, self%ntheta
                self%collated_data(:,1,:,:) = self%collated_data(:,1,:,:) + &
                        data_copy(:,t,:,:)*self%theta_weights(t)
            Enddo   

            DeAllocate(data_copy)

        Endif

    End Subroutine Collate_Physical

    Subroutine Write_Data(self,disp,file_unit,filename)
        Implicit None
        Class(io_buffer) :: self
        Integer, Intent(In), Optional :: file_unit
        Character*120, Intent(In), Optional :: filename
        Integer :: funit, ierr, j
        Logical :: error
		Integer(kind=MPI_OFFSET_KIND), Intent(In), Optional :: disp
        Integer(kind=MPI_OFFSET_KIND) :: my_disp, hdisp, tdisp, fdisp, bdisp
		Integer :: mstatus(MPI_STATUS_SIZE)

        hdisp = 0
        If (present(disp)) hdisp = disp
      
        ! The file can be opened previously or opened by this routine
        error = .false.

        If (present(file_unit)) Then
            ! The file is already open. We are modifying it.
            funit = file_unit
        Else
            If (present(filename)) Then
                ! The file is not open.  We must create it.
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

            If (self%output_rank) Allocate(self%buffer(1:self%io_buffer_size))

            If (self%write_mode .eq. 1) Call self%gather_data(-1)
            
            ! First Write the Data
            Do j = 1, self%nwrites

                If (self%write_mode .gt. 1) Call self%gather_data(j)

                If (self%output_rank) Then
                    fdisp = self%file_disp(j)+hdisp
                    bdisp = self%buffer_disp(j)
                    Call MPI_File_Seek( funit, fdisp, MPI_SEEK_SET, ierr) 
                    Call MPI_FILE_WRITE(funit, self%buffer(bdisp), & 
                        self%buffsize, MPI_DOUBLE_PRECISION, mstatus, ierr)
                Endif
            Enddo

            ! Next, write timestamps as needed
            If (self%write_timestamp) Then  ! (output_rank 0)
                Do j = 1, self%nrec

                    tdisp = j*self%qdisp*self%ncache+hdisp  ! May need to account for real/imaginary here.
                    Call MPI_File_Seek(funit,tdisp,MPI_SEEK_SET,ierr)

                    Call MPI_FILE_WRITE(funit, self%time(j), 1, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)

                    Call MPI_FILE_WRITE(funit, self%iter(j), 1, & 
                           MPI_INTEGER, mstatus, ierr)
                Enddo
            Endif    

            If (self%output_rank) Then
                DeAllocate(self%buffer)
                If (present(filename)) Then
			        Call MPI_FILE_CLOSE(funit, ierr) 
                Endif
            Endif

        Endif
    End Subroutine Write_Data

    Subroutine Read_Data(self,disp,file_unit,filename)
        Implicit None
        Class(io_buffer) :: self
        Integer, Intent(In), Optional :: file_unit
        Character*120, Intent(In), Optional :: filename
        Integer :: funit, ierr, j
        Logical :: error
		Integer(kind=MPI_OFFSET_KIND), Intent(In), Optional :: disp
        Integer(kind=MPI_OFFSET_KIND) :: my_disp, hdisp, tdisp, fdisp, bdisp
		Integer :: mstatus(MPI_STATUS_SIZE)

        hdisp = 0
        If (present(disp)) hdisp = disp
      
        ! The file can be opened previously or opened by this routine

        error = .false.

        If (present(file_unit)) Then
            ! The file is already open. We are modifying it.
            funit = file_unit
        Else
            If (present(filename)) Then
                ! The file is not open.  We must create it.
                If (self%output_rank) Then
		            Call MPI_FILE_OPEN(self%ocomm%comm, filename, & 
                           MPI_MODE_RDONLY, & 
                           MPI_INFO_NULL, funit, ierr) 
                Endif
            Else
                error = .true.
            Endif
        Endif 

        If (.not. error) Then

            If (self%output_rank) Then
                Allocate(self%buffer(1:self%in_buffer_size))
                self%buffer(:) = 0.0d0
                !Write(6,*)'inbuffer: ', self%in_buffer_size, self%nlm_in, self%nr_out
            Endif
            ! Read the Data
            Do j = 1, self%nwrites
                If (self%output_rank) Then
                    fdisp = self%file_disp_in(j)+hdisp
                    !Write(6,*)'j is: ',j, fdisp, self%in_buffer_size
                    !bdisp = self%buffer_indisp(j)
                    Call MPI_File_Seek( funit, fdisp, MPI_SEEK_SET, ierr) 
                    Call MPI_FILE_READ(funit, self%buffer(1), & 
                        self%in_buffer_size, MPI_DOUBLE_PRECISION, mstatus, ierr)
                Endif
                Call self%distribute_data(j)
            Enddo

            If (self%output_rank) Then
                DeAllocate(self%buffer)
                If (present(filename)) Then
			        Call MPI_FILE_CLOSE(funit, ierr) 
                Endif
            Endif

        Endif
    End Subroutine Read_Data

End Module Parallel_IO
