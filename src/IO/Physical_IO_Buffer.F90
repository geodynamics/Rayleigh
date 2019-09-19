Module Physical_IO_Buffer
    Use RA_MPI_Base
    Use Parallel_Framework
    Type, Public :: IO_Buffer_Physical


        ! MPI-related variables
        Integer :: col_rank, row_rank
        
        !General variables, true for each column and row rank
        Integer :: nr            ! number of r-values in potentially subsampled array
        Integer :: ntheta        ! number of theta values in (potentially subsampled) array
        Integer :: nphi          ! number of phi values in the (potentially subsampled) array
        Integer :: ncache        ! number of fields (or fields/time indices) in the array
        Integer :: cache_index=1 ! current index within cache

        Integer :: nr_local     ! number of r-values subsampled locally (possibly 0)
        Integer :: ntheta_local     ! number of theta-values subsampled locally (possibly 0)
        
        Integer, Allocatable :: phi_ind(:)  ! phi indices in (potentially subsampled) array 

        !**Row-specific** variables related to the cascade operation
        !Everything here is true within a given row, but may different between
        !rows
        Integer, Allocatable :: ntheta_at_column(:)  ! number of theta points for each row rank
        Integer, Allocatable :: theta_inds_at_column(:,:) ! theta indices at each column
      
        !Variables only needed by output ranks
        Integer, Allocatable :: nr_out_at_column(:) ! number of radii output by row rank X in this row
        Integer, Allocatable :: nr_out_at_row(:)    ! number of radii output by column rank Y's row


        ! Variables that (may) help efficiency of subsampling
        Logical :: simple   = .false.   ! no subsampling (3-D output)
        Logical :: general  = .false.   ! subsampling in all 3 dimensions (point-probe output)
        Logical :: r_general = .false.   ! subsampling in radius only (Shell Slices)
        !Logical :: rp_simple = .false.  ! 

        ! Buffer-specific variables
        Real*8, Allocatable :: cache(:,:,:,:)   ! storage space
        Logical :: output_rank = .false.        ! True if this rank performs I/O
        Type(communicator) :: ocomm             ! communicator for parallel output
        Integer :: nout_cols = 0                ! Number of columns that participate in I/O for this subsample
    Contains
        Procedure :: init => Initialize_Physical_IO_Buffer
        Procedure :: cache_data
        Procedure :: allocate_cache
        Procedure :: cascade
        Procedure :: simple_cascade

    End Type IO_Buffer_Physical

Contains

    Subroutine Initialize_Physical_IO_Buffer(self,r_indices, theta_indices, &
                                             phi_indices,ncache)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Integer, Intent(In), Optional :: r_indices, theta_indices, phi_indices
        Integer, Intent(In), Optional :: ncache
        Integer :: p
        Integer :: nout_cols
        Write(6,*)'Initializing!'

        !Philosophy:  Initialization happens once.  It should be legible.

        !//////////////////////////////////////////////
        ! Establish my identity
        self%col_rank = pfi%ccomm%rank
        self%row_rank = pfi%rcomm%rank


        !///////////////////////////////////////////////////////////
        ! Establish how many fields will be output
        If (.not. present(ncache)) Then
            self%ncache = 1
        Else
            self%ncache = ncache
        Endif
        

        !///////////////////////////////////////////////////////////
        ! Establish how this buffer will be subsampled, if at all
        If (.not. present(r_indices)) Then
            If (.not. present(theta_indices)) Then
                If (.not. present(phi_indices)) Then
                    self%simple = .true.

                Endif
            Endif
        Else
            self%r_general=.true.
        Endif


        If (self%simple) Then
            Write(6,*)'ncache: ', self%ncache
            self%nr     = pfi%n1p
            self%ntheta = pfi%n2p
            self%nphi   = pfi%n3p

            self%ntheta_local = pfi%all_2p(self%row_rank)%delta
            self%nr_local     = pfi%all_1p(self%col_rank)%delta

            Allocate(self%ntheta_at_column(0:pfi%nprow-1))
            Do p = 0, pfi%nprow-1
                self%ntheta_at_column(p) = pfi%all_2p(p)%delta
            Enddo

            !Handle multi-column I/O
            nout_cols = pfi%output_columns
            If (nout_cols .gt. pfi%my_1p%delta) Then
                nout_cols = pfi%my_1p%delta
            Endif
            self%nout_cols = nout_cols
            If (self%col_rank .lt. nout_cols) self%output_rank = .true.
            Write(6,*)'output columns: ', self%nout_cols


        Endif



        Call self%Allocate_Cache()

    End Subroutine Initialize_Physical_IO_Buffer

    Subroutine Allocate_Cache(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Allocate(self%cache(1:self%nphi, 1:self%nr_local, &
                            1:self%ntheta_local, 1:self%ncache))
        self%cache(:,:,:,:) = 0.0d0
    End Subroutine Allocate_Cache

    Subroutine Cache_Data(self,vals)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Real*8, Intent(In) :: vals(1:,1:,1:)
        Integer :: r, t, p
        
        If (self%simple) Then
            Do t = 1, self%ntheta_local
                Do r = 1, self%nr_local
                    Do p = 1, self%nphi
                        self%cache(p,r,t,self%cache_index) = vals(p,r,t)
                    Enddo
                Enddo
            Enddo
        Endif

        self%cache_index = MOD(self%ncache,self%cache_index+1)
        Write(6,*)'Caching!'
    End Subroutine Cache_Data

    Subroutine Cascade(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self

        If (self%simple) call self%simple_cascade()

    End Subroutine Cascade

    Subroutine Simple_Cascade(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self
        Write(6,*)'Simple Cascade!'

    End Subroutine Simple_Cascade

    Subroutine Write_Data(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self

    End Subroutine Write_Data


End Module Physical_IO_Buffer
