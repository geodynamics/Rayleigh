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

!///////////////////////////////////////////////////////////////////////////////////////
!       "Let the yolk fall from your shoulders.  Don't carry it all, don't carry it all.
!           We are all our hands and holders beneath this bold and brilliant sun."
!                               - The Decemberists
!///////////////////////////////////////////////////////////////////////////////////////
Module Parallel_Framework
    Use MPI_LAYER
    Use General_MPI
    Use Load_Balance
    Use Structures
    Use BufferedOutput
    Implicit None
    Private
#ifdef usemkl
    include 'mkl_service.fi'
#endif
    !///////////////////////////////////////////////////////////////////////////
    !
    Integer, Parameter, Public :: Cartesian = 1, Cylindrical = 2, Spherical = 3
    Integer, Parameter, Public :: p1 =1 ,s1 = 2, p2 = 3,s2a =4, p3a=5,p3b=6
    Public :: Load_Config, Full_Barrier, Main_MPI_Init
    Character*6 :: ifmt = '(i4.4)' ! Integer format for indicating processor tag numbers in output
    Type, Public :: Parallel_Interface
        Type(Load_Config) :: my_1p, my_2p, my_3p    !  like my_r, my_theta in ASH
        Type(Load_Config) :: my_1s, my_2s, my_3s    !     like my_mp, my_n etc.
        Type(Load_Config), Allocatable :: all_1p(:), all_2p(:)
        Type(Load_Config), Allocatable :: all_3s(:)
        Integer :: Geometry
        ! Global dimensions in 'a' and 'b' configurations
        ! Allows primarily for dealiasing right now
        Integer :: n1p, n1s, n2p, n2s, n3p, n3s
        Integer :: npe, nprow, npcol, npio,npc
        Integer :: nthreads = 1
        Integer, Allocatable :: inds_3s(:)
        Type(communicator) :: rcomm ! row communicator
        Type(communicator) :: ccomm ! column communicator
        Type(communicator) :: gcomm ! global communicator for individual run.  Mirrors mpi_comm_world in standard mode
        Type(communicator) :: wcomm ! Mirror of MPI_COMM_WORLD
        !##################################################
        Type(Load_Config), Allocatable, Public :: lb_1p(:), lb_1s(:)


        Contains
        Procedure :: Init_World_Comm
        Procedure :: Init => Initialize_Parallel_Interface
        Procedure :: openmp_init
        Procedure :: Exit => Finalize_Framework
        Procedure :: Spherical_Init
        Procedure :: Init_Geometry
        Procedure :: Broadcast_Intarr
    End Type Parallel_Interface

    Type, Public :: mcontainer
        Complex*16, Allocatable :: data(:,:,:)
    End Type mcontainer
    !Type, Public :: rmcontainer
    !    Real*8, Allocatable :: data(:,:)
    !End Type rmcontainer
    Type, Public :: eqcontainer
        Real*8, Allocatable :: data(:,:,:)
    End Type eqcontainer




    Type(Parallel_Interface), Public :: pfi



Contains
    Subroutine Full_Barrier()
        Call Barrier(pfi%gcomm)
    End Subroutine Full_Barrier



    Subroutine Main_MPI_Init(ret_rank)
        Implicit None
        Integer, Intent(out) :: ret_rank
        Call pfi%init_world_comm()
        ret_rank = pfi%wcomm%rank
    End Subroutine Main_MPI_Init

    Subroutine Init_World_Comm(self)
        !This handles initialization of MPI_COMM_WORLD
        Implicit None
        Integer :: error
        Class(Parallel_Interface) :: self
        self%wcomm = init_main_group(error)
    End Subroutine Init_World_Comm


    Subroutine Initialize_Parallel_Interface(self, pars,ncpus,init_only)
        Implicit None
        Integer, Intent(In) ::  pars(1:)
        Integer, Intent(In) :: ncpus(1:)
        Integer :: pcheck, error
        Integer :: ierr
        Character*6 :: istr
        Logical, Intent(In) :: init_only
        Class(Parallel_Interface) :: self
        self%geometry = pars(1)
        self%n1p = pars(2)
        self%n1s = pars(3)
        self%n2p = pars(4)
        self%n2s = pars(5)
        self%n3p = pars(6)
        self%n3s = pars(7)
        pcheck = size(pars,1)
        If (pcheck .gt. 7) Then
            self%npe = pars(8)
            self%nprow = pars(9)
            self%npcol = pars(10)
            self%npio = pars(8)-pars(9)*pars(10)
            self%npc = pars(10)*pars(9)
        Else
            ! get the process grid info from the environment
        Endif
        If (size(ncpus) .gt. 1) Then
            !Multiple run Mode
            self%gcomm = Init_SubGroup(self%wcomm,ncpus,error)
        Else
            !Normal Mode -- gcomm is mirror of wcomm/MPI_COMM_WORLD
            self%gcomm%rank = self%wcomm%rank
            self%gcomm%np = self%wcomm%np
            self%gcomm%comm = self%wcomm%comm

        Endif
        If (self%gcomm%rank .eq. 0) Then

            Call stdout%print(" //////////////////////////////////////")
            Call stdout%print(" Initializating Rayleigh...")
            Write(istr,'(i6)')self%npe
            call stdout%print(" ")
            call stdout%print(" -- Initalizing MPI...")
            Call stdout%print(" ---- Specified parameters:")
            call stdout%print(" ---- NCPU  : "//trim(istr))
            Write(istr,'(i6)')self%nprow
            call stdout%print(" ---- NPROW : "//trim(istr))
            Write(istr,'(i6)')self%npcol
            call stdout%print(" ---- NPCOL : "//trim(istr))
        Endif
        If (self%nprow .le. 0) Then
            If (self%gcomm%rank .eq. 0) Then
                Call stdout%print('........................................................')
                Call stdout%print(' --- Error:  nprow must be at least 1.  Exiting...')
            Endif
            Call self%exit()
        Endif

        If (self%npcol .le. 0) Then
            If (self%gcomm%rank .eq. 0) Then
                Call stdout%print('........................................................')
                Call stdout%print(' --- Error:  npcol must be at least 1.  Exiting...')
            Endif
            Call self%exit()
        Endif
        if (self%gcomm%np .ne. self%npe) Then
            If (self%gcomm%rank .eq. 0) Then
                Call stdout%print('........................................................')
                Call stdout%print(' --- Error NCPU does not agree with number of MPI Ranks!')
                Write(istr,'(i6)')self%gcomm%np
                Call stdout%print(' --- NCPU from MPI   : '//trim(istr))
                Write(istr,'(i6)')self%npe
                Call stdout%print(' --- NCPU from Input : '//trim(istr))
                Call stdout%print('........................................................')
                Call stdout%print(' Exiting...')
            Endif
            Call self%exit()
        Endif


            Call rowcolsplit(self%gcomm,self%rcomm,self%ccomm,self%nprow,error)
        If (.not. init_only) Then
            !Load balancing does not HAVE to be initialized (though it usually is)
            Call self%Init_Geometry()
        Endif

        Call self%openmp_init()
        If (self%gcomm%rank .eq. 0) Then
            call stdout%print(" -- MPI initialized.")
            call stdout%print(" ")
        Endif

    End Subroutine Initialize_Parallel_Interface
    Subroutine Broadcast_Intarr(self,intarr,src,comm_option)
        Implicit None
        Class(Parallel_Interface) :: self
        Integer, Intent(In) :: intarr(:)
        Integer, Intent(In), Optional :: src, comm_option
        Integer :: arr_size, comm_handle, src_rank, ierr
        arr_size = size(intarr)
        If (present(src)) Then
            src_rank = src
        Else
            src_rank = 0
        Endif
        comm_handle = self%wcomm%comm
        If (present(comm_option)) Then
            Select Case(comm_option)
                Case(1)
                    comm_handle = self%gcomm%comm
                Case(2)
                    comm_handle = self%rcomm%comm
                Case(3)
                    comm_handle = self%ccomm%comm
                Case Default
                    comm_handle = self%wcomm%comm
            End Select
        Endif
        Call MPI_Bcast(intarr,arr_size, MPI_INTEGER, src_rank, comm_handle, ierr)
    End Subroutine Broadcast_Intarr

    Subroutine OpenMp_Init(self)
#ifdef useomp
        Use Omp_lib
#endif

        Class(Parallel_Interface) :: self
        Integer :: my_mpi_rank,my_thread
        my_mpi_rank = pfi%gcomm%rank
#ifdef useomp
        self%nthreads = omp_get_max_threads()
        !$OMP PARALLEL PRIVATE(my_thread)
        my_thread = omp_get_thread_num()
        write(6,*)'rank: ', my_mpi_Rank, 'thread: ', my_thread, 'nthread: ', omp_get_num_threads()
        !$OMP END PARALLEL
#endif
#ifdef usemkl
            my_thread = mkl_get_max_threads()
            write(6,*)"MKL MAX: ", my_thread
#endif
    End Subroutine OpenMp_Init

    Subroutine Init_Geometry(self)
        Class(Parallel_Interface) :: self
        If (self%geometry == Spherical) Then
            Call self%Spherical_Init()
        Endif
    End Subroutine Init_Geometry


    Subroutine Spherical_Init(self)
        Integer :: r, unit1
        Class(Parallel_Interface) :: self

        !Distribute radii amongst members of a column
        Allocate(self%all_1p(0:self%npcol-1))
        Call Standard_Balance(self%all_1p,self%n1p,self%ccomm)
        self%my_1p = self%all_1p(self%ccomm%rank)

        !Distribute theta amongst members of a row
        Allocate(self%all_2p(0:self%nprow-1))
        Call Standard_Balance(self%all_2p,self%n2p,self%rcomm)
        self%my_2p = self%all_2p(self%rcomm%rank)

        ! Determine number of m-values per process

        Allocate(self%all_3s(0:self%nprow-1))
        Call Standard_Balance(self%all_3s,self%n3s,self%rcomm)
        self%my_3s = self%all_3s(self%rcomm%rank)
        !Write(6,*), self%my_3s%min, self%my_3s%max,self%rcomm%rank, self%my_3s%delta
        ! We assume a high-low pairing of m-values.
        ! This means that we need to set up an index array
        ! for the m-values that tells how many each processor has.
        Allocate(self%inds_3s(1:self%n3s))
        Call m_balance(self%all_3s, self%inds_3s, self%rcomm)
        Call LM_Load_Balance(pfi%my_3s,self%inds_3s,self%ccomm)
        If (self%gcomm%rank .eq. -100) Then
            unit1 = 10
            Open(unit1, file = 'verification/m_check', status='replace', form = 'unformatted')
            Write(unit1) self%n3s
            Write(unit1) (self%inds_3s(r),r=1,self%n3s)

            Close(unit1)
        Endif
    End Subroutine Spherical_Init



    Subroutine Finalize_Framework(self)
        Class(Parallel_Interface) :: self
        Integer :: error
        self%n1p = 0    ! This line is here purely so that the intel compiler does not
        ! throw an unused variable warning when warn-all is used.
        Call Exit_Comm_Lib(error)
        !STOP
    End Subroutine Finalize_Framework

End Module Parallel_Framework
