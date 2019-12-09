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

Module Spherical_IO
    Use Spherical_Buffer
	Use Parallel_Framework
	Use SendReceive
    Use ISendReceive
    Use General_MPI
    Use Fourier_Transform
    Use Legendre_Transforms, Only : Legendre_Transform
    Use BufferedOutput
    Use Math_Constants
    Use MakeDir
#ifdef INTEL_COMPILER 
    USE IFPORT
#endif
    Use Parallel_IO
	Implicit None
	! This module contains routines for outputing spherical data as:
	! 1. Slices of sphere
	! 2. Phi-Averages over sphere (f(r,theta))
	! 3. Phi/theta Averages over sphere ((f(r))
	! 4. Volume averages over sphere
    ! 5. 3-D Output on Spherical grid
    ! 6. Spectra on slices of the sphere
    ! 7. Sampled via individual spherical harmonic modes
    ! 8. Meridional Slices
    ! 9. Equatorial Slices

	!////////////////////////////////////////////
    Integer, Parameter :: nqmax=4000, nshellmax=2048, nmeridmax=8192, nmodemax=2048
    Integer, Parameter :: nprobemax=4096
    Integer, Parameter :: endian_tag = 314      ! first 4 bits of each diagnostic file - used for assessing endianness on read-in
    Integer, Parameter :: reallybig = 90000000
    Integer, Parameter :: max_nheader=20
    ! Each diagnostic type has an associated version number that is written to the file
    ! following the endian tag.  Hopefully, reading routines can be backward compatible if
    ! significant changes are made to the output structure (and reflected in the version number)
    ! Version numbers are not assumed to be in sync.  Instead, they reflect how many times
    ! a particular output has been modified substantially.

    Integer, Parameter :: shellslice_version = 3
    Integer, Parameter :: azavg_version = 3
    Integer, Parameter :: shellavg_version = 4
    Integer, Parameter :: globalavg_version = 3
    Integer, Parameter :: shellspectra_version = 3
    Integer, Parameter :: equslice_version = 1  
    Integer, Parameter :: meridslice_version = 1
    Integer, Parameter :: sphmode_version =4
    INTEGER, PARAMETER :: probe_version = 1
    Integer, Parameter :: full3d_version = 3    !currently unused
    Type, Public :: DiagnosticInfo
        ! Need to see if we can make these allocatable, but for now..
        ! Each instance of this class has two static arrays used for reading in namelist input
        ! These arrays are large enough to hold nqmax and nshellmax values, but 
        ! typically only a small fraction of that amount will be specified at input.
        Integer :: values(1:nqmax) = -1 ! The list of values specified in an input namelist
        Integer :: levels(1:nshellmax) = -1 ! The radial indices output (shell slices, spectra only)
        Integer :: compute(1:nqmax)= -1 ! compute(i) = 1 if i was specified in values.  -1 otherwise
        
        Integer :: nq, nlevels ! Number of nonzero elements of values and levels
        Integer :: my_nlevels  ! Number of nonzero elements of levels that are in process

        Integer :: file_unit = 15
        Character*120 :: file_prefix = 'None'

        Integer, Allocatable :: oqvals(:)   ! Array of size nq used by I/O process to record output ordering of diagnostics

        Integer :: frequency = reallybig ! How often we write this diagnostic type
        Integer :: rec_per_file =1     ! How many of these records we write to a file

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! These variables control file opening and positioning
        Integer :: current_rec = 1        ! Which record we are on within a file
        Integer :: file_header_size =0    ! Size of file header in bytes
        Integer :: file_record_size = 0   ! Size of each record in bytes
        Integer :: file_position = 1      ! Keep track of where we are (byte-wise) within the file
        Logical :: file_open = .false.    ! Was the file opened successfully?
        Logical :: write_header = .false. ! Do we need to write header information
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        Integer :: avg_level = 0       ! global_averages = 3, shell_averages = 1, az_averages = 1, all others = 0

        Integer :: ind = 1              ! index of current diagnostic being stored (advances throughout an output calculation)
        Logical :: begin_output = .false.
        Integer :: mpi_tag = 1          ! For use when communicating before writing the file

        Integer :: output_version = -1  ! See note at top

        !This flag changes throughout the diagnostic evaluation
        !If the qty being output at the current iteration is supposed to be
        !output as this diagnostic type, the grab_this_q will be set to true.
        !Otherwise grab_this_q remains false.
        Logical :: grab_this_q = .false.

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Variables used for management of SPH_Mode_Samples output
        INTEGER, Allocatable :: l_values(:)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !  Variables used for cached output (currently only used for point probes)
        INTEGER :: cache_size = 1
        INTEGER :: cc = 0  ! cache counter
        INTEGER, ALLOCATABLE :: iter_save(:)
        REAL*8 , ALLOCATABLE :: time_save(:)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Legacy variables
        Integer :: nphi_global, nr_global, ntheta_global
        Integer, Allocatable :: phi_global(:), theta_global(:), r_global(:)
        Real*8, Allocatable :: r_values(:), theta_values(:), phi_values(:)
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Variables for subsampling the IO buffer
        Integer, Allocatable :: r_inds(:), theta_inds(:), phi_inds(:)
        Integer :: nr, ntheta, nphi, nell
        Real*8, Allocatable :: r_vals(:), theta_vals(:), phi_vals(:)

        !Communicatory Info for parallel writing (if used)
        Integer :: ocomm, orank, onp
        Logical :: master = .false.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! For header storage
        Type(rmcontainer1d) :: rvals(max_nheader)   ! array of real-arrays to write to header
        Type(imcontainer1d) :: ivals(max_nheader)   ! array of integer-array to write to header
        Integer :: buff_sizes(max_nheader)           ! size of each header line to write
        Integer :: buff_types(max_nheader)           ! buffer type (1 -> double; 2 -> integer)
        Integer :: nrbuffer = 0                     ! number of real and integer buffers
        Integer :: nibuffer = 0 
        Integer :: nheader = 0                      ! nibuffer+nrbuffer
        Integer(kind=MPI_OFFSET_KIND) :: hdisp = 12
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !I/O Buffer
        Type(io_buffer) :: buffer
        Integer :: write_mode =1

        ! Methods
        Contains
        Procedure :: Init => Initialize_Diagnostic_Info
        Procedure :: AdvanceInd
        Procedure :: AdvanceCC
        Procedure :: reset => diagnostic_output_reset
        Procedure :: OpenFile
        Procedure :: OpenFile_Par
        Procedure :: CloseFile
        Procedure :: CloseFile_Par
        Procedure :: Write_Header_Data
        Procedure :: Write_IO
        Procedure :: Add_DHeader
        Procedure :: Add_IHeader
        Procedure :: Set_IHeader_Line
        Procedure :: Store_Values
        Procedure :: update_position
        Procedure :: getq_now
        Procedure :: CleanUp
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    End Type DiagnosticInfo

    Type(DiagnosticInfo) :: Shell_Averages, Shell_Slices, Global_Averages, AZ_Averages, Full_3D, Shell_Spectra
    Type(DiagnosticInfo) :: Equatorial_Slices, Meridional_Slices, SPH_Mode_Samples, Point_Probes
    Type(DiagnosticInfo) :: temp_io

    Type(io_buffer) :: full_3d_buffer
    
    Integer :: current_averaging_level = 0
    Integer :: current_qval = 0

    Integer :: averaging_level(1:nqmax) = 0, compute_q(1:nqmax) = 0
    Integer :: shellavg_values(1:nqmax)=-1, globalavg_values(1:nqmax)=-1
    Integer :: shellslice_values(1:nqmax) =-1, shellslice_levels(1:nshellmax)=-1, azavg_values(1:nqmax)=-1
    Integer :: full3d_values(1:nqmax) = -1, shellspectra_values(1:nqmax)=-1, shellspectra_levels(1:nshellmax)=-1
    Integer :: equatorial_values(1:nqmax) = -1, meridional_values(1:nqmax) = -1, meridional_indices(1:nmeridmax) = -1
    Integer :: SPH_Mode_Values(1:nqmax) = -1, SPH_Mode_ell(1:nmodemax) = -1, SPH_Mode_Levels(1:nshellmax) = -1
    INTEGER :: point_probe_values(1:nqmax)

    INTEGER :: SPH_Mode_nell = 0, SPH_Mode_nmode = 0

    INTEGER :: point_probe_r(1:nprobemax) = -1, point_probe_theta(1:nprobemax) = -1, point_probe_phi(1:nprobemax) = -1

    Integer :: globalavg_nrec = 1, shellavg_nrec = 1, azavg_nrec = 1, shellslice_nrec =1, shellspectra_nrec =1
    Integer :: equatorial_nrec = 1, meridional_nrec=1, sph_mode_nrec, point_probe_nrec =1

    Integer :: globalavg_frequency = reallybig, shellavg_frequency = reallybig
    Integer :: azavg_frequency = reallybig, shellslice_frequency = reallybig
    Integer :: shellspectra_frequency=reallybig, equatorial_frequency = reallybig
    Integer :: meridional_frequency=reallybig, sph_mode_frequency = reallybig
    Integer :: point_probe_frequency = reallybig

    Integer :: full3d_frequency= reallybig
    Character*120 :: local_file_path=''
    Logical :: mem_friendly = .false.

    REAL*8 ::   shellslice_levels_nrm(1:nshellmax) = -3.0d0
    REAL*8 :: shellspectra_levels_nrm(1:nshellmax) = -3.0d0
    REAL*8 ::     sph_mode_levels_nrm(1:nshellmax) = -3.0d0
    REAL*8 ::  meridional_indices_nrm(1:nmeridmax) = -3.0d0

    REAL*8 ::     point_probe_r_nrm(1:nprobemax)  = -3.0d0
    REAL*8 :: point_probe_theta_nrm(1:nprobemax)  = -3.0d0
    REAL*8 ::   point_probe_phi_nrm(1:nprobemax)  = -3.0d0
    Integer :: point_probe_cache_size = 1

    Namelist /output_namelist/ mem_friendly, &
        ! All output types require that a minimal set of information is specified 
          meridional_values,   meridional_frequency,     meridional_nrec, &
          equatorial_values,   equatorial_frequency,     equatorial_nrec, &
          shellslice_values,   shellslice_frequency,     shellslice_nrec, &
        shellspectra_values, shellspectra_frequency,   shellspectra_nrec, &      
               azavg_values,        azavg_frequency,          azavg_nrec, &
            shellavg_values,     shellavg_frequency,       shellavg_nrec, &
           globalavg_values,    globalavg_frequency,      globalavg_nrec, &
            sph_mode_values,     sph_mode_frequency,       sph_mode_nrec, &
         point_probe_values,  point_probe_frequency,    point_probe_nrec, &
              full3d_values,       full3d_frequency,                      &
        !Some outputs have additional information that needs to be specified
         meridional_indices,      shellslice_levels, shellspectra_levels, &
              point_probe_r,      point_probe_theta,     point_probe_phi, &
               sph_mode_ell,                             sph_mode_levels, &

        shellslice_levels_nrm, shellspectra_levels_nrm, meridional_indices_nrm, &
        point_probe_r_nrm    ,     point_probe_phi_nrm,  point_probe_theta_nrm, &
        sph_mode_levels_nrm,        point_probe_cache_size

    Integer :: integer_zero = 0
    Real*8, Private, Allocatable :: qty(:,:,:), f_of_r_theta(:,:,:), f_of_r(:)
    Real*8, Private, Allocatable :: shellav_outputs(:,:,:), globav_outputs(:)

    Character*8, Public :: i_ofmt = '(i8.8)'

    Integer :: io_node = 0

    Integer, Private :: current_iteration
    !///////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ! We store some redundant information just to keep the IO level as independent of the Physics level as possible
    ! These variables are all private.
    Real*8, Allocatable,Private :: theta_integration_weights(:), r_integration_weights(:)
    Integer, Private :: my_theta_min, my_theta_max, my_rmin, my_rmax, my_mp_max, my_mp_min
    Integer, Private :: nr, nphi, ntheta, my_ntheta
    Integer, Private :: nproc1, nproc2, myid, nproc, my_row_rank, my_column_rank, my_nr
    Real*8, Allocatable, Private :: radius(:),sintheta(:),costheta(:) 
    Real*8 :: over_nphi_double

    !////////////////////////////////////////////////////////////
    ! And here are some variables used for storing averages so that moments may be taken
    ! These two arrays contain the ell0 and m0 representation of all quantities
    !   computed in shell- and (eventually) azimuthal- averaging.
    ! These two averages are used to compute moments
    Real*8, Private, Allocatable :: IOell0_values(:,:), IOm0_values(:,:,:)
    Integer, Private :: num_avg_store ! number of averages to store in the two arrays above
    Integer, Private :: IOavg_flag = -1
Contains

    Subroutine Begin_Outputting(iter)
		  Implicit None
		  Integer, Intent(In) :: iter
		  current_iteration = iter

        Call   Global_Averages%reset()
        Call    Shell_Averages%reset()
        Call       AZ_Averages%reset()
        Call      Shell_Slices%reset()
        Call     Shell_Spectra%reset()
        Call Equatorial_Slices%reset()
        Call Meridional_Slices%reset()
        Call  SPH_Mode_Samples%reset()
        Call      Point_Probes%reset()

        !Call Temp_IO%reset()

        Allocate(f_of_r_theta(1,my_rmin:my_rmax,my_theta_min:my_theta_max))
        Allocate(f_of_r(my_rmin:my_rmax))

        num_avg_store = shell_averages%nq
        Allocate(IOm0_values(my_rmin:my_rmax,my_theta_min:my_theta_max,1:num_avg_store))
        Allocate(IOell0_values(my_rmin:my_rmax,1:num_avg_store))
        IOm0_values(:,:,:) = 0.0d0
        IOell0_values(:,:) = 0.0d0
    End Subroutine Begin_Outputting

	Subroutine Initialize_Spherical_IO(rad_in,sintheta_in, rw_in, tw_in, costheta_in,file_path,digfmt)
		Implicit None
		Integer :: k, fcount(3,2), ntot, fcnt, master_rank, i
        Integer :: nn, cascade_type, wmode
		Real*8, Intent(In) :: rad_in(:), sintheta_in(:), rw_in(:), tw_in(:), costheta_in(:)
        Character*120 :: fdir
        Character*120, Intent(In) :: file_path
        Character*8, Intent(In) :: digfmt
        Integer, Allocatable :: gpars(:,:)

        ! Set format code for integer file names
        i_ofmt=digfmt

        local_file_path = file_path
		! Handles output bookkeeping

		my_theta_min = pfi%my_2p%min
		my_theta_max = pfi%my_2p%max
		my_rmin = pfi%my_1p%min
		my_rmax = pfi%my_1p%max
        
        my_mp_min = pfi%my_3s%min
        my_mp_max = pfi%my_3s%max

		my_nr = pfi%my_1p%delta
        my_ntheta = pfi%my_2p%delta

        nproc  = pfi%gcomm%np
		nproc1 = pfi%ccomm%np		! processor's per column (radius split among these)
		nproc2 = pfi%rcomm%np      ! processor's per row  (theta split among these)
		myid   = pfi%gcomm%rank
		my_row_rank = pfi%rcomm%rank
		my_column_rank = pfi%ccomm%rank
		nphi = pfi%n3p
		ntheta = pfi%n2p
		nr = pfi%n1p

        over_nphi_double = 1.0d0/nphi

		Allocate(sintheta(1:ntheta))
		Allocate(costheta(1:ntheta))

		sintheta(:) = sintheta_in(:)
		costheta(:) = costheta_in(:)
		Allocate(radius(1:nr))
		radius(:) = rad_in(:)

        Allocate(r_integration_weights(1:nr))
        Allocate(theta_integration_weights(1:ntheta))

        r_integration_weights(:) = rw_in(:)
        theta_integration_weights(:) = tw_in(:)		

        ! Before initializing the different output types, we check to see if
        ! the user has specified shell-slice levels etc. using normalized
        ! physical coordinates or negative-index shorthand

        Call Process_Coordinates()

        ! Map the various quantity lists etc. into their associated diagnostic structures
        ! Numbers here are the mpi_tags used in communication for each output
        ! In theory they can be the same, but it's probably a good idea to keep them unique.

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! First, initialize outputs that do not rely on the I/O Buffer data structure
        ! for caching etc.
        fdir = 'G_Avgs/'
        Call Global_Averages%Init(averaging_level,compute_q,myid, 55, fdir, &
                          globalavg_version, globalavg_nrec, globalavg_frequency, &
                          values = globalavg_values, avg_level = 3, nobuffer=.true.) 
       
        fdir = 'Shell_Avgs/'
        Call Shell_Averages%Init(averaging_level,compute_q,myid, 57, fdir, &
                          shellavg_version, shellavg_nrec, shellavg_frequency, &
                          values = shellavg_values, avg_level = 2, nobuffer=.true.) 
       
        fdir = 'Spherical_3D/'  ! Note: Full 3D only uses frequency, not nrec
        Call Full_3D%Init(averaging_level,compute_q,myid, 54, fdir, &
                          full3d_version, 1, full3d_frequency, &
                          values = full3d_values, nobuffer=.true.) 

        Allocate(gpars(1:5,1:5))
        gpars(:,:) = -1
        Call Full_3D_Buffer%Init(gpars, mpi_tag=54) ! Full 3-D uses a separate buffer (due to io semantics)
        DeAllocate(gpars)

        ! Now, initialize those outputs that use the new I/O for caching and parallel output.

        wmode = 1
        if (mem_friendly) wmode = 2

        fdir = 'Shell_Slices/'
        Call Shell_Slices%Init(averaging_level,compute_q,myid, 58, fdir, &
                          shellslice_version, shellslice_nrec, shellslice_frequency, &
                          values = shellslice_values, rinds = shellslice_levels, &
                          write_mode = wmode)

        fdir = 'Point_Probes/'
        Call Point_Probes%Init(averaging_level,compute_q,myid, 63, fdir, &
                          probe_version, point_probe_nrec, point_probe_frequency, &
                          values = point_probe_values, cache_size = point_probe_cache_size, &
                          rinds = point_probe_r, tinds = point_probe_theta, & 
                          pinds = point_probe_phi)


        fdir = 'Meridional_Slices/'
        Call Meridional_Slices%Init(averaging_level,compute_q,myid, 61, fdir, &
                          meridslice_version, meridional_nrec, meridional_frequency, &
                          values = meridional_values, pinds = meridional_indices)

        fdir = 'AZ_Avgs/'
        Call AZ_Averages%Init(averaging_level,compute_q,myid, 56, fdir, &
                          azavg_version, azavg_nrec, azavg_frequency, &
                          values = azavg_values, average_in_phi = .true. )

        fdir = 'Equatorial_Slices/'
        Call Equatorial_Slices%Init(averaging_level,compute_q,myid, 60, fdir, &
                          equslice_version, equatorial_nrec, equatorial_frequency, &
                          values = equatorial_values, tinds=(/ ntheta/2, ntheta/2+1 /), &
                          tweights = (/ 0.5d0, 0.5d0 /) ) 

        wmode =1
        If (mem_friendly) wmode =2
        fdir = 'Shell_Spectra/'
        Call Shell_Spectra%Init(averaging_level,compute_q,myid, 59, fdir, &
                          shellspectra_version, shellspectra_nrec, shellspectra_frequency, &
                          values = shellspectra_values, rinds=shellspectra_levels, &
                          is_spectral = .true. , write_mode = wmode) 

        fdir = 'SPH_Modes/'
        wmode =1
        If (mem_friendly) wmode =2
        Call SPH_Mode_Samples%Init(averaging_level,compute_q,myid, 62, fdir, &
                          sphmode_version, sph_mode_nrec, sph_mode_frequency, &
                          values = sph_mode_values, rinds=sph_mode_levels, &
                          is_spectral = .true. , write_mode = wmode, lvals=sph_mode_ell) 

        Call Initialize_Headers()

   End Subroutine Initialize_Spherical_IO

    Subroutine Initialize_Headers
        Implicit None
        Integer :: dims(1:4), lmax

        !///////////////////////////////////////////////////
        ! Point Probes Header
        ! Dimensions ------------------
        dims(1:4) = (/ Point_Probes%nr, Point_Probes%ntheta, &
                       Point_Probes%nphi, Point_Probes%nq /)
        Call Point_Probes%Add_IHeader(dims,4)
        Call Point_Probes%Add_IHeader(Point_Probes%oqvals*0,Point_Probes%nq)
        ! Radial grid -----------------------------------
        Call Point_Probes%Add_DHeader(Point_Probes%r_vals,Point_Probes%nr)
        Call Point_Probes%Add_IHeader(Point_Probes%r_inds, Point_Probes%nr)
        !Theta grid-----------------------------------------
        Call Point_Probes%Add_DHeader(Point_Probes%theta_vals,Point_Probes%ntheta)
        Call Point_Probes%Add_IHeader(Point_Probes%theta_inds, Point_Probes%ntheta)
        !Phi grid-----------------------------------------
        Call Point_Probes%Add_DHeader(Point_Probes%phi_vals,Point_Probes%nphi)
        Call Point_Probes%Add_IHeader(Point_Probes%phi_inds, Point_Probes%nphi)

        !////////////////////////////////////////////////////
        ! Shell Slices
        dims(1:3) = (/ ntheta, Shell_Slices%nr, Shell_Slices%nq /)
        Call Shell_Slices%Add_IHeader(dims,3)
        Call Shell_Slices%Add_IHeader(Shell_Slices%oqvals, Shell_Slices%nq)

        Call Shell_Slices%Add_DHeader(Shell_Slices%r_vals, Shell_Slices%nr)
        Call Shell_Slices%Add_IHeader(Shell_Slices%r_inds, Shell_Slices%nr)
        Call Shell_Slices%Add_DHeader(costheta,ntheta)  

        !//////////////////////////////////////////////////////
        ! Meridional Slices
        dims(1:4) = (/ nr, ntheta, Meridional_Slices%nphi, Meridional_Slices%nq /)       
        Call Meridional_Slices%Add_IHeader(dims,4)
        Call Meridional_Slices%Add_IHeader(Meridional_Slices%oqvals, Meridional_Slices%nq)
        Call Meridional_Slices%Add_DHeader(radius,nr)
        Call Meridional_Slices%Add_DHeader(costheta,ntheta)
        Call Meridional_Slices%Add_IHeader(Meridional_Slices%phi_inds, Meridional_Slices%nphi)

        !//////////////////////////////////////////////////////
        ! AZ_Averages
        dims(1:3) = (/ nr, ntheta,  AZ_Averages%nq /)       
        Call AZ_Averages%Add_IHeader(dims,3)
        Call AZ_Averages%Add_IHeader(AZ_Averages%oqvals, AZ_Averages%nq)
        Call AZ_Averages%Add_DHeader(radius,nr)
        Call AZ_Averages%Add_DHeader(costheta,ntheta)

        !//////////////////////////////////////////////////////
        ! Equatorial Slices
        dims(1:3) = (/ nphi, nr,  Equatorial_Slices%nq /)       
        Call Equatorial_Slices%Add_IHeader(dims,3)
        Call Equatorial_Slices%Add_IHeader(Equatorial_Slices%oqvals, Equatorial_Slices%nq)
        Call Equatorial_Slices%Add_DHeader(radius,nr)

        !/////////////////////////////////////////////////////
        ! Shell Spectra
        lmax = maxval(pfi%inds_3s)
        dims(1:3) = (/ lmax, Shell_Spectra%nr, Shell_Spectra%nq /)
        Call Shell_Spectra%Add_IHeader(dims,3)
        Call Shell_Spectra%Add_IHeader(Shell_Spectra%oqvals, Shell_Spectra%nq)
        Call Shell_Spectra%Add_DHeader(Shell_Spectra%r_vals, Shell_Spectra%nr)
        Call Shell_Spectra%Add_IHeader(Shell_Spectra%r_inds, Shell_Spectra%nr)

        !/////////////////////////////////////////////////////
        ! SPH Modes
        lmax = maxval(pfi%inds_3s)
        dims(1:3) = (/ SPH_Mode_Samples%nell, SPH_Mode_Samples%nr, SPH_Mode_Samples%nq /)
        Call SPH_Mode_Samples%Add_IHeader(dims,3)
        Call SPH_Mode_Samples%Add_IHeader(SPH_Mode_Samples%oqvals, SPH_Mode_Samples%nq)
        Call SPH_Mode_Samples%Add_DHeader(SPH_Mode_Samples%r_vals, SPH_Mode_Samples%nr)
        Call SPH_Mode_Samples%Add_IHeader(SPH_Mode_Samples%r_inds, SPH_Mode_Samples%nr)
        Call SPH_Mode_Samples%Add_IHeader(SPH_Mode_Samples%l_values, SPH_Mode_Samples%nell)

    End Subroutine Initialize_Headers

    Subroutine Store_Values(self,qty)
        Implicit None
        Class(DiagnosticInfo) :: self
        REAL*8, INTENT(IN) :: qty(:,my_rmin:,my_theta_min:)
 
        Call self%buffer%cache_data(qty)
        self%oqvals(self%ind) = current_qval
        Call self%AdvanceInd()

    End Subroutine Store_Values

	Subroutine Write_IO(self,this_iter,simtime)
        USE RA_MPI_BASE
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
        Class(DiagnosticInfo) :: self
		INTEGER(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_pdisp, new_disp, qdisp, full_disp, tdisp
        Logical :: responsible, output_rank
        Integer :: orank, funit, error, ncache, omode

	    If ((self%nq > 0) .and. (Mod(this_iter,self%frequency) .eq. 0 )) Then

            Call self%Buffer%timestamp(this_iter, simtime)

            If ((self%cache_size-1) .eq. self%cc) Then

                orank       = self%Buffer%ocomm%rank
                output_rank = self%Buffer%output_rank
                responsible = .false.
                ncache = self%cache_size

                If ((output_rank) .and. (orank .eq. 0) ) responsible = .true.
         
                If (output_rank) Call self%OpenFile_Par(this_iter, error)
                funit   = self%file_unit

                If ( responsible .and. self%file_open ) Then   
                    funit = self%file_unit
                    If (self%current_rec .eq. ncache) Then                                        
                        If ( responsible .and. (self%write_header) ) Then    !           
                            Call self%Set_IHeader_Line(2,self%oqvals) 
                            Call self%Write_Header_Data()
                        Endif
                    Endif
                Endif

                full_disp = self%buffer%qdisp*self%buffer%ncache+12  ! 12 is for the simtime+iteration at the end; use ncache instead of nq to account for real/imaginary in spectral i/o
                new_disp = self%hdisp+full_disp*(self%current_rec-ncache)

                !omode = 1
                !If (mem_friendly .and. (self%cascade_type .ne. 1)) omode = 3
                Call self%buffer%write_data(disp=new_disp,file_unit=funit) ! , mode = omode)
            
                If (output_rank) Call self%closefile_par()

            Endif            

            Call self%AdvanceCC() 

        Endif

	End Subroutine Write_IO

    Function Compute_Quantity(qval) result(yesno)
        integer, intent(in) :: qval 
        Logical             :: yesno 
        current_qval = qval
        yesno = .false.

        If (IOAvg_Flag .eq. 1) Then
            ! We check only the shell_averages so we can compute averages for moments
            If (Shell_Averages%compute(qval) .eq. 1) Then
                Call Shell_Averages%getq_now(yesno)
            Endif
        Else
            ! Otherwise, check everything - normal output
            If (compute_q(qval) .eq. 1) Then

                ! We check all diagnostic types and their respective frequencies
                ! While doing so, we modify the averaging level
                Call Full_3D%getq_now(yesno)
                Call Shell_Slices%getq_now(yesno)
                Call Equatorial_Slices%getq_now(yesno)

                Call Meridional_Slices%getq_now(yesno)

                Call SPH_Mode_Samples%getq_now(yesno)
                Call Point_Probes%getq_now(yesno)

                !Call Temp_IO%getq_now(yesno)

                Call Shell_Spectra%getq_now(yesno)
                Call AZ_Averages%getq_now(yesno)
                Call Shell_Averages%getq_now(yesno)
                Call Global_Averages%getq_now(yesno)

            Endif
        Endif
    End function Compute_quantity

    Function Sometimes_Compute(qval) result(yesno)
        integer, intent(in) :: qval 
        Logical             :: yesno 
        yesno = .false.
        ! If we ever compute qval, return true
        ! Return false otherwise
        If (compute_q(qval) .eq. 1) yesno = .true.
    End function Sometimes_Compute

    Function Time_To_Output(iter) result(yesno)
        integer, intent(in) :: iter 
        Logical             :: yesno 
        yesno = .false.

         If ((Global_Averages%nq > 0) .and. (Mod(iter,Global_Averages%Frequency) .eq. 0)) yesno = .true.
         If ((Shell_Averages%nq > 0) .and. (Mod(iter,Shell_Averages%Frequency) .eq. 0)) yesno = .true.        
         If ((Shell_Spectra%nq > 0) .and. (Mod(iter,Shell_Spectra%Frequency) .eq. 0)) yesno = .true.            
         If ((Equatorial_Slices%nq > 0) .and. (Mod(iter,Equatorial_Slices%Frequency) .eq. 0)) yesno = .true. 

         If ((Meridional_Slices%nq > 0) .and. (Mod(iter,Meridional_Slices%Frequency) .eq. 0)) yesno = .true. 

         If ((SPH_Mode_Samples%nq > 0) .and. (Mod(iter,SPH_Mode_Samples%Frequency) .eq. 0)) yesno = .true. 
         If ((Point_Probes%nq > 0) .and. (Mod(iter,Point_Probes%Frequency) .eq. 0)) yesno = .true. 
     

         If ((AZ_Averages%nq > 0) .and. (Mod(iter,AZ_Averages%Frequency) .eq. 0)) yesno = .true.
         If ((Shell_Slices%nq > 0) .and. (Mod(iter,Shell_Slices%Frequency) .eq. 0)) yesno = .true. 
         If ((Full_3D%nq > 0) .and. (Mod(iter,Full_3D%Frequency) .eq. 0)) yesno = .true.

         !If ((Temp_IO%nq > 0) .and. (Mod(iter,Temp_IO%Frequency) .eq. 0)) yesno = .true. 
       
    End function Time_To_Output

    Subroutine Finalize_Averages()
        ! Use the azimuthal averages to compute the ell=0 mean
        Call IOComputeEll0(IOm0_values,IOell0_values)
        Call Shell_Averages%reset() !need to reset the index counter
    End Subroutine Finalize_Averages

    Subroutine Set_Avg_Flag(flag_val)
        Integer, Intent(In) :: flag_val
        IOavg_flag = flag_val
    End Subroutine Set_Avg_Flag

	Subroutine Add_Quantity(qty)
		Implicit None
		Real*8, Intent(In) :: qty(:,:,:)

        If (IOavg_flag .eq. 1) Then
            !Compute and store the azimuthal average of qty
            If (shell_averages%grab_this_q) Then
                Call IOComputeM0(qty)
            Endif
        Else

            !If (Temp_IO%grab_this_q)           Call Temp_IO%Store_Values(qty)
            If (AZ_Averages%grab_this_q)       Call AZ_Averages%Store_Values(qty)  

            If (Equatorial_Slices%grab_this_q) Call Equatorial_Slices%Store_Values(qty) 
            If (Meridional_Slices%grab_this_q) Call Meridional_Slices%Store_Values(qty)
            If (SPH_Mode_Samples%grab_this_q)  Call SPH_Mode_Samples%Store_Values(qty)
            If (Point_Probes%grab_this_q)      Call Point_Probes%Store_Values(qty)

		    If (Shell_Slices%grab_this_q)      Call Shell_Slices%Store_Values(qty)
		    If (Shell_Spectra%grab_this_q)     Call Shell_Spectra%Store_Values(qty)

            Call Get_Averages(qty)  ! -> Global Averages, Shell Averages, Azimuthal Averages

		    If (full_3d%grab_this_q) Call write_full_3d(qty)
        Endif

	End Subroutine	Add_Quantity

	Subroutine Complete_Output(iter, sim_time)
		Integer, Intent(In) :: iter
		Real*8, Intent(In) :: sim_time

        Call Shell_Slices%write_io(iter,sim_time)    
        Call Equatorial_Slices%Write_IO(iter,sim_time)
        Call Meridional_Slices%write_io(iter,sim_time)
        Call Point_Probes%write_io(iter,sim_time)
        Call AZ_Averages%write_io(iter,sim_time)
        Call Shell_Spectra%write_io(iter,sim_time)
        !Call Temp_IO%write_io(iter,sim_time)
        Call SPH_Mode_Samples%write_io(iter, sim_time)

	    If ((Shell_Averages%nq > 0) .and. (Mod(iter,Shell_Averages%frequency) .eq. 0 )) Call Write_Shell_Average(iter,sim_time)
	    If ((Global_Averages%nq > 0) .and. (Mod(iter,Global_Averages%frequency) .eq. 0 )) Call Write_Global_Average(iter,sim_time)

        DeAllocate(f_of_r_theta)
        DeAllocate(f_of_r)
        DeAllocate(IOm0_values)
        DeAllocate(IOell0_values)

	End Subroutine Complete_Output

	Subroutine Get_Averages(qty)
        ! Takes azimuthal, shellular, and partial global averages of 3D array qty
		Implicit None
		Integer :: azav_ind,shellav_ind, globav_ind
        Integer :: i, r,t,p,m
        Real*8 :: this_average, wght
		Real*8, Intent(In) :: qty(1:,my_rmin:,my_theta_min:)

        !//////////////////////////////
        !First the azimuthal average
        If (current_averaging_level .ge. 1) Then

		    f_of_r_theta(1,:,:) = 0.0D0
		
		    Do t = my_theta_min, my_theta_max
			    Do r = my_rmin, my_rmax
				    f_of_r_theta(1,r,t) = sum(qty(:,r,t))
			    Enddo
		    Enddo

		    f_of_r_theta = f_of_r_theta*over_nphi_double     ! average in phi

        Endif

        !/////////////////////
        !Now Average over the partial sphere (shell_averages)
		If (current_averaging_level .ge. 2) Then
		    f_of_r(:) = 0.0D0

            Do t = my_theta_min, my_theta_max
                f_of_r(:) = f_of_r(:) + f_of_r_theta(1,:,t)*theta_integration_weights(t)
            Enddo

		    If (Shell_Averages%grab_this_q) Then

                If (.not. Allocated(shellav_outputs)) Then
			        Allocate(shellav_outputs(my_rmin:my_rmax,1:4,1:Shell_Averages%nq))  ! four moments
                    shellav_outputs(:,:,:) = 0.0d0
                Endif

                shellav_ind = Shell_Averages%ind

                !First, add the partially integrated spherically symmetric mean ( f_of_r )
                Do r = my_rmin, my_rmax
                    shellav_outputs(r,1,shellav_ind) = f_of_r(r)
                Enddo

                ! Now, the moments - rms, skewness, kurtosis
                Do m = 2, 4
                    Do t = my_theta_min, my_theta_max
                        wght = theta_integration_weights(t)*over_nphi_double
                        Do r = my_rmin, my_rmax
                            Do p = 1, nphi
                            shellav_outputs(r,m,shellav_ind) = shellav_outputs(r,m,shellav_ind) + &
                                & wght*(qty(p,r,t)-IOell0_values(r,shellav_ind))**m
                            Enddo
                            
                        Enddo
                    Enddo
                Enddo

                If (myid .eq. 0) Shell_Averages%oqvals(shellav_ind) = current_qval
                Call Shell_Averages%AdvanceInd()
		    Endif
        Endif

        !////////////////////////
        ! Finally, integrate in radius to get the partial global average
        If (current_averaging_level .ge. 3) Then

            this_average =0.0d0
            Do i = my_rmin, my_rmax
                this_average = this_average+f_of_r(i)*r_integration_weights(i)
            Enddo

		    If (Global_Averages%grab_this_q) Then
    		    If (.not. Allocated(globav_outputs)) Then			
    			    Allocate(globav_outputs(1:Global_Averages%nq))			
    		    Endif
                globav_ind = Global_Averages%ind
			    globav_outputs(globav_ind) = this_average
                If (myid .eq. 0) Global_Averages%oqvals(globav_ind) = current_qval
                Call Global_Averages%AdvanceInd()
		    Endif

        Endif

	End Subroutine Get_Averages

    Subroutine Write_Global_Average(this_iter,simtime)
        Implicit None
        Real*8, Allocatable :: buff(:), full_avg(:)

        Integer :: responsible
        Integer :: i,n, nq_globav, global_avg_tag, error, file_pos
        Integer :: funit
        Integer, Intent(In) :: this_iter
        Real*8, Intent(In) :: simtime

        nq_globav = Global_Averages%nq
        global_avg_tag = Global_Averages%mpi_tag
        funit = Global_Averages%file_unit

        !///////////////////////////////
        ! Sum across rows, and then across the first column

        Allocate(full_avg(nq_globav))
        Allocate(buff(nq_globav))

        full_avg(:) = 0.0d0
        buff(:) = 0.0d0
        Call dsum1d(globav_outputs,buff,pfi%rcomm)
        If (my_row_rank .eq. 0) Then
            Call dsum1d(buff,full_avg,pfi%ccomm)
        Endif

        !//////////////////////////////
        If (myid .eq. io_node) Then

            Call Global_Averages%OpenFile(this_iter, error)

            If (error .eq. 0) Then
                If (Global_Averages%write_header) Then
                    Write(funit)nq_globav
                    Write(funit)(Global_Averages%oqvals(i),i=1,nq_globav)
                    Call Global_Averages%update_position
                Endif

                Write(funit)(full_avg(i),i=1,nq_globav)
                Write(funit)simtime
                Write(funit)this_iter
                Call Global_Averages%CloseFile() 
            Endif

        Endif

        If (Allocated(globav_outputs)) DeAllocate(globav_outputs)
        DeAllocate(buff)
        DeAllocate(full_avg)

    End Subroutine Write_Global_Average

	Subroutine Write_Shell_Average(this_iter, simtime)
		Implicit None
        Integer, Intent(In) :: this_iter
		Real*8, Intent(In) :: simtime
		Integer :: i,j, k, n, nn, nq_shellav, shell_avg_tag, m 
        Real*8, Allocatable :: full_shellavg(:,:,:), buff(:,:,:), buff2(:,:,:)		
        Integer :: your_r_min, your_r_max, your_nr, your_id
        Integer :: funit, error, inds(3), ncount, sirq,nirq
        Integer, Allocatable :: rirqs(:)

        shell_avg_tag = Shell_Averages%mpi_tag
        nq_shellav = Shell_Averages%nq
        funit = Shell_Averages%file_unit

        !Sum across the row to complete integration in theta
        Allocate(buff(my_rmin:my_rmax,1:4,nq_shellav))
        buff(:,:,:) = 0.0d0
        Call dsum3d(shellav_outputs,buff,pfi%rcomm)

        If (my_row_rank .eq. 0) Then
            !now set up a series of isends/ireceives along the column
            If (my_column_rank .eq. 0) Then
                !Post ireceives before anything else is done
                nirq = nproc1-1
                Allocate(rirqs(1:nirq))
                Allocate(full_shellavg(1:nq_shellav,1:4,1:nr))
                Do n = 1, nproc1-1
                    your_r_min = pfi%all_1p(n)%min
                    your_nr = pfi%all_1p(n)%delta
                    ncount = your_nr*nq_shellav*4
                    inds(1) = 1
                    inds(2) = 1
                    inds(3) = your_r_min
                    Call IReceive(full_shellavg, rirqs(n),n_elements = ncount, &
                            &  source= n,tag = shell_avg_tag, grp = pfi%ccomm,indstart = inds)
                Enddo
                ! Load my buff into the full_shellavg array
                Do j = my_rmin,my_rmax
                    Do m = 1, 4
                        Do i = 1, nq_shellav
                            full_shellavg(i,m,j) = buff(j,m,i)
                        Enddo
                    Enddo
                Enddo
                DeAllocate(buff)
                Call IWaitAll(nirq, rirqs)
                DeAllocate(rirqs)
            Else
                Allocate(buff2(nq_shellav, 1:4,my_rmin:my_rmax)) ! Transpose for easier send logic
                Do j = my_rmin,my_rmax
                    Do m = 1, 4
                        Do i = 1, nq_shellav
                            buff2(i,m,j) = buff(j,m,i)
                        Enddo
                    Enddo
                Enddo
                ncount = my_nr*nq_shellav*4
                Call ISend(buff2, sirq,ncount, dest = 0, tag = shell_avg_tag, grp = pfi%ccomm)

                Call IWait(sirq)
                DeAllocate(buff2)
                DeAllocate(buff)
            Endif
        Endif

        If (myid .eq. io_node) Then ! Rank Zero writes the file
            Call Shell_Averages%OpenFile(this_iter, error)
            If (error .eq. 0) Then
                If (Shell_Averages%write_header) Then
                    Write(funit)nr,nq_shellav
                    Write(funit)(Shell_Averages%oqvals(i),i=1,nq_shellav)
                    Write(funit)(radius(i),i=1,nr)
                    Call Shell_Averages%update_position()
                Endif

                Write(funit)(((full_shellavg(k,m,i),i=1,nr),m=1,4),k=1,nq_shellav)
                Write(funit) simtime

                Write(funit)this_iter


                Call Shell_Averages%CloseFile()

                DeAllocate(full_shellavg)
            Endif
		Endif

        DeAllocate(shellav_outputs)

	End Subroutine Write_Shell_Average

	Subroutine Write_Full_3D(qty)
		Implicit None		
		Real*8, Intent(In) :: qty(:,my_rmin:,my_theta_min:)
		Integer :: i, funit
		Character*4 :: qstring
		Character*120 :: iterstring, data_file, grid_file

        ! Write the data file (Parallel I/O)
        Write(iterstring,i_ofmt) current_iteration
        Write(qstring,'(i4.4)') current_qval
        data_file = trim(local_file_path)//'Spherical_3D/'//trim(iterstring)//'_'//qstring
        Call full_3d_buffer%cache_data(qty)
        Call full_3d_buffer%write_data(filename=data_file)

        ! Rank 0 writes the grid file.
        If (myid .eq. 0) Then
            Write(iterstring,'(i8.8)') current_iteration
            grid_file = trim(local_file_path)//'Spherical_3D/'//trim(iterstring)//'_'//'grid'
            Open(newunit=funit,file=grid_file,form='unformatted', status='replace', access='stream')
            Write(funit)endian_tag
            Write(funit)nr
            Write(funit)ntheta
            Write(funit)nphi
            Write(funit)(radius(i),i=1,nr)
            Write(funit)(acos(costheta(i)),i = 1, ntheta)
            Close(funit)
        Endif

	End Subroutine Write_Full_3D

    !////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !       Diagnostic Class Methods

    Subroutine Write_Header_Data(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer :: mstatus(MPI_STATUS_SIZE)
        Integer :: ierr, bsize, btype, i, funit, nbuff, ii,di
        funit = self%file_unit
        nbuff = self%nheader
        di = 1
        ii = 1
        Do i = 1, nbuff
            bsize = self%buff_sizes(i)
            btype = self%buff_types(i)
            If (btype .eq. 1) Then
                CALL MPI_FILE_WRITE(funit, self%rvals(di)%data(1), bsize, &
                                    MPI_DOUBLE_PRECISION, mstatus, ierr) 
                di = di+1
            Else
                CALL MPI_FILE_WRITE(funit, self%ivals(ii)%data(1), bsize, &
                                    MPI_INTEGER, mstatus, ierr) 
                ii = ii+1
            Endif
        Enddo

    End Subroutine Write_Header_Data

    Subroutine Add_IHeader(self,invals,n)
        Class(DiagnosticInfo) :: self
        Integer, Intent(In) :: invals(1:), n
        self%nibuffer = self%nibuffer+1
        self%nheader = self%nheader+1
        self%buff_types(self%nheader) = 2
        self%buff_sizes(self%nheader) = n
        Allocate(self%ivals(self%nibuffer)%data(1:n))
        self%ivals(self%nibuffer)%data(1:n) = invals(1:n)
        self%hdisp = self%hdisp+4*n
    End Subroutine Add_IHeader

    Subroutine Set_IHeader_Line(self,ind,vals)
        Class(DiagnosticInfo) :: self
        Integer, Intent(In) :: vals(1:), ind
        self%ivals(ind)%data(:) = vals(:)
    End Subroutine Set_IHeader_Line

    Subroutine Add_DHeader(self,invals,n)
        Class(DiagnosticInfo) :: self
        Real*8, Intent(In) :: invals(1:)
        Integer, Intent(In) :: n
        self%nrbuffer = self%nrbuffer+1
        self%nheader = self%nheader+1
        self%buff_types(self%nheader) = 1
        self%buff_sizes(self%nheader) = n
        Allocate(self%rvals(self%nrbuffer)%data(1:n))
        self%rvals(self%nrbuffer)%data(1:n) = invals(1:n)
        self%hdisp = self%hdisp+8*n
    End Subroutine Add_DHeader

    Subroutine Initialize_Diagnostic_Info(self,avg_levels,computes,pid,mpi_tag, &
                 dir, version, nrec, frequency, avg_level,values, &
                 levels, phi_inds, cache_size, rinds, tinds, pinds, write_mode, &
                 tweights, is_spectral, lvals, nobuffer, &
                 average_in_phi, average_in_radius, average_in_theta)
        Implicit None
        Integer :: i,ind
        Integer, Intent(In) :: pid, mpi_tag, version, nrec, frequency
        Integer, Intent(In), Optional :: cache_size, write_mode
        Character*120, Intent(In) :: dir
        Integer, Optional, Intent(In) :: avg_level
        Integer, Optional, Intent(In) :: values(1:)
        Integer, Optional, Intent(In) :: levels(1:)
        Integer, Optional, Intent(In) :: phi_inds(1:)
        Integer, Intent(InOut) :: computes(1:), avg_levels(1:)
        Integer, Intent(In), Optional :: rinds(1:), tinds(1:), pinds(1:), lvals(1:)

        !Averaging variables
        Real*8, Intent(In), Optional :: tweights(1:)
        Logical, Intent(In), Optional :: is_spectral
        Logical, Intent(In), Optional :: nobuffer, average_in_phi, average_in_theta, average_in_radius
        !Real*8, Allocatable th_weights(:)
        Integer :: rcount, tcount, pcount, lcount, modcheck, nmax
        Logical ::  spectral_io, spectral_buffer
        Integer :: avg_axes(1:3), ecode
        Integer, Allocatable :: indices(:,:)
        Real*8, Allocatable :: avg_weights(:,:)
        Class(DiagnosticInfo) :: self 

        ! File info
        self%file_prefix = dir
        self%output_version = version
        self%rec_per_file = nrec
        self%frequency = frequency

        If (myid .eq. 0) Call Make_Directory(trim(local_file_path)//trim(self%file_prefix) ,ecode)
        ! Indices array
        nmax = Maxval((/ nr, ntheta, nphi /) )
        Allocate(indices(1:nmax,1:5))
        If (present(tweights)) Then
            Allocate(avg_weights(1:size(tweights),1:3))
            avg_weights(:,2) = tweights
        Else
            Allocate(avg_weights(1,1:3))
            avg_weights(:,:) = -1.0d0
        Endif
        indices(:,:) = -1


        !Check that the cache size is appropriate
        If (present(cache_size)) Then
            If (cache_size .ge. 1) Then
                self%cache_size = cache_size
            Endif
        Endif
        If (self%cache_size .lt. 1) Then
            If (myid .eq. 0) Then
                    Write(6,*)'////////////////////////////////////////////////////////////////////'
                    Write(6,*)'   Warning:  Incorrect cache_size specification for ',self%file_prefix
                    Write(6,*)'   Cache_size must be at least 1.'
                    Write(6,*)'   Specified cache_size: ', self%cache_size
                    Write(6,*)'   Caching has been deactivated for ', self%file_prefix
                    Write(6,*)'////////////////////////////////////////////////////////////////////'                    
            Endif
            self%cache_size = 1
        Endif
        If (self%cache_size .gt. self%rec_per_file) Then
            modcheck = MOD(self%rec_per_file,self%cache_size)
            If (modcheck .ne. 0) Then

                If (myid .eq. 0) Then
                    Write(6,*)'////////////////////////////////////////////////////////////////////'
                    Write(6,*)'   Warning:  Incorrect cache_size specification for ',self%file_prefix
                    Write(6,*)'   Cache_size cannot be larger than nrec.'
                    Write(6,*)'   Cache_size: ', self%cache_size
                    Write(6,*)'   nrec      : ', self%rec_per_file
                    Write(6,*)'   Cache_size has been set to nrec.'
                    Write(6,*)'////////////////////////////////////////////////////////////////////'                    
                Endif
                self%cache_size = self%rec_per_file
            Endif
        Endif
        
        If (self%cache_size .gt. 1) Then
            modcheck = MOD(self%rec_per_file,self%cache_size)
            If (modcheck .ne. 0) Then

                If (myid .eq. 0) Then
                    Write(6,*)'////////////////////////////////////////////////////////////////////'
                    Write(6,*)'   Warning:  Incorrect cache_size specification for ',self%file_prefix
                    Write(6,*)'   Cache_size must divide evenly into nrec.'
                    Write(6,*)'   Cache_size: ', self%cache_size
                    Write(6,*)'   nrec      : ', self%rec_per_file
                    Write(6,*)'   Caching has been deactivated for ', self%file_prefix
                    Write(6,*)'////////////////////////////////////////////////////////////////////'                    
                Endif
                self%cache_size = 1
            Endif
        Endif

        If (present(write_mode)) Then
            self%write_mode = write_mode
        Else
            self%write_mode = 1
        Endif

        If (present(avg_level)) Then
            if (avg_level .gt. 0) self%avg_level = avg_level
        Endif
        self%mpi_tag = mpi_tag
        self%nq = 0

        Allocate(self%iter_save(1:self%cache_size))
        Allocate(self%time_save(1:self%cache_size))

        If (present(values)) Then
            self%values(:) = values(:) 
            
            Do i = 1, nqmax
                if(self%values(i) .gt. 0) Then 
                    self%nq = self%nq+1
                    ind = self%values(i)
                    self%compute(ind) = 1
                    computes(ind) = 1
                    If (avg_levels(ind) .lt. self%avg_level) Then
                        avg_levels(ind) = self%avg_level
                    Endif
                endif 
            Enddo
        Endif

        Allocate(self%oqvals(1:self%nq))
        self%oqvals(:) = nqmax+100

        ! Buffer init
        If (present(rinds)) Then
            rcount = size(rinds)
            i = 1
            Do While ( i .le. rcount )
                If (rinds(i) .lt. 0) Then
                    rcount = i-1
                Endif
                i = i+1
            Enddo

            If (rcount .gt. 0) Then
                Allocate(self%r_inds(1:rcount))
                Allocate(self%r_vals(1:rcount))
                self%r_inds(1:rcount) = rinds(1:rcount)
                Do i = 1, rcount
                    self%r_vals(i) = radius(rinds(i))
                Enddo
                indices(1,5) = rcount
                indices(1:rcount,1) = self%r_vals(1:rcount)
            Endif

        Else
            Allocate(self%r_inds(1:1))
            Allocate(self%r_vals(1:1))
            self%r_inds(1) = -1
            self%r_vals(1) = -1
        Endif

        If (present(tinds)) Then
            tcount = size(tinds)
            i = 1
            DO While ( i .le. tcount )
                If (tinds(i) .lt. 0) Then
                    tcount = i-1
                Endif
                i = i+1
            Enddo

            If (tcount .gt. 0) Then
                Allocate(self%theta_inds(1:tcount))
                Allocate(self%theta_vals(1:tcount))
                self%theta_inds(1:tcount) = tinds(1:tcount)
                Do i = 1, tcount
                    self%theta_vals(i) = costheta(tinds(i))
                Enddo
                indices(2,5) = tcount
                indices(1:tcount,2) = self%theta_vals(1:tcount)
            Endif

        Else
            Allocate(self%theta_inds(1:1))
            Allocate(self%theta_vals(1:1))
            self%theta_inds(1) = -1
            self%theta_vals(1) = -1
        Endif

        If (present(pinds)) Then
            pcount = size(pinds)
            i = 1
            Do While ( i .le. pcount )
                If (pinds(i) .lt. 0) Then
                    pcount = i-1
                Endif
                i = i+1
            Enddo
            If (pcount .gt. 0) Then
                Allocate(self%phi_inds(1:pcount))
                Allocate(self%phi_vals(1:pcount))
                self%phi_inds(1:pcount) = pinds(1:pcount)
                Do i = 1, pcount
                    self%phi_vals(i) = (pinds(i)-1)*(two_pi/nphi)   
                Enddo
                indices(3,5) = pcount
                indices(1:pcount,3) = self%phi_vals(1:pcount)
            Endif
        Else
            Allocate(self%phi_inds(1:1))
            Allocate(self%phi_vals(1:1))
            self%phi_inds(1) = -1
            self%phi_vals(1) = -1
        Endif


        If (present(lvals)) Then
            lcount = size(lvals)
            i = 1
            Do While ( i .le. lcount )
                If (lvals(i) .lt. 0) Then
                    lcount = i-1
                Endif
                i = i+1
            Enddo

            If (lcount .gt. 0) Then
                Allocate(self%l_values(1:lcount))
                self%l_values(1:lcount) = lvals(1:lcount)
                self%nell = lcount
                indices(1:self%nell,4) = lvals(1:lcount)
                indices(4,5) = self%nell
            Endif
        Else
            Allocate(self%l_values(1:1))
            self%l_values(1) = -1
        Endif

        self%nr     = rcount
        self%ntheta = tcount
        self%nphi   = pcount

        ! Initialize the buffer (and then ocomm)
        avg_axes(:) = 0  ! phi, r, theta
        If (present(average_in_phi)) Then
            If (average_in_phi) avg_axes(1) = 1            
        Endif

        spectral_io = .false.
        If (present(is_spectral)) spectral_io = is_spectral

        Call self%buffer%init( indices, &
                              mode = self%write_mode, mpi_tag = self%mpi_tag, &
                              ncache  = self%nq*self%cache_Size, &
                              nrec = self%cache_size, skip = 12, &
                              write_timestamp = .true., averaging_axes = avg_axes, &
                              averaging_weights = avg_weights, spectral=spectral_io)

        DeAllocate(indices)
        DeAllocate(avg_weights)

        If (.not. present(nobuffer)) Then
            ! Once the buffer is initialized, set the ocomm info
            self%ocomm = self%buffer%ocomm%comm
            self%orank = self%buffer%ocomm%rank
            self%onp = self%buffer%ocomm%np
            If (self%orank .eq. 0) then
                self%master = .true.    ! This process handles file headers in parallel IO
            Endif
        Endif

    End Subroutine Initialize_Diagnostic_Info

    Subroutine AdvanceInd(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        self%ind = self%ind+1
        self%begin_output = .false.
    End Subroutine AdvanceInd

    Subroutine AdvanceCC(self)
        !Advances the cache counter
        Implicit None
        Class(DiagnosticInfo) :: self
        self%cc = self%cc+1
        self%cc = MOD(self%cc,self%cache_size)
    End Subroutine AdvanceCC

    Subroutine Diagnostic_Output_Reset(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        self%ind = 1
        self%begin_output = .true.
    End Subroutine Diagnostic_Output_Reset

    Subroutine OpenFile(self,iter,errcheck)
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer, Intent(In) :: iter
        Integer, Intent(InOut) :: errcheck
        Character*120 :: iterstring,istr
        Character*120 :: filename, omsg
        Integer :: modcheck, imod, file_iter, next_iter, ibelong
        Integer :: fstat
        Logical :: create_file
        Logical :: file_exists

        create_file = .false.
        file_exists = .false.

        self%file_open    = .false.
        self%write_header = .false.

        ! Determine which file # this data belongs to
        modcheck = self%frequency*self%rec_per_file
        imod = Mod(iter,modcheck) 
        If (imod .eq. 0) Then
            ibelong = iter
        Else
            ibelong = iter-imod+modcheck    ! Data belongs to file #ibelong
        Endif

        Write(iterstring,i_ofmt) ibelong
        filename = trim(local_file_path)//trim(self%file_prefix)//trim(iterstring)

        If ( (imod .eq. self%frequency) .or. (self%rec_per_file .eq. 1) ) Then   
            ! Time to begin a new file 
            create_file = .true.
        Else
            ! We are either:
            ! (a)  Appending to an existing file or
            ! (b)  The user has changed the output cadence in between checkpoints.
            !      In that case, we're "off-cycle" and need to create the file.
            Inquire(File=filename, Exist=file_exists)
            If (.not. file_exists) create_file = .true.
        Endif

        
        If (create_file) Then
            Call stdout%print(' Creating Filename: '//trim(filename))
            Open(unit=self%file_unit,file=filename,form='unformatted', &
                & status='replace',access='stream',iostat = errcheck)

            self%write_header = .true.            
        Else
            Open(unit=self%file_unit,file=filename,form='unformatted', status='old',access='stream', &
                & iostat = errcheck, POSITION = 'APPEND')    
            self%write_header = .false.
        Endif

        If (errcheck .eq. 0) Then
            self%file_open = .true.
        Else
            Call stdout%print(' Unable to open file: '//trim(filename))
        Endif
        
        If (self%file_open) Then

            If (self%write_header) Then

                Write(self%file_unit)endian_tag
                Write(self%file_unit)self%output_version
                ! We write zero initially - only update nrec after the data is actually written
                Write(self%file_unit)integer_zero 
                self%current_rec = 1            
            Else

                Call self%update_position                    ! save current position
                Read(self%file_unit,POS = 9)self%current_rec ! read previous record #
                self%current_rec = self%current_rec+1
#ifdef INTEL_COMPILER 
                fstat=fseek(self%file_unit, self%file_position, 0) ! return to end of file 
#else
                Call fseek(self%file_unit, self%file_position, 0, fstat) ! return to end of file 
#endif
            Endif
        Endif

    End Subroutine OpenFile

    Subroutine OpenFile_Par(self,iter,ierr)
        !Performs the same tasks as OpenFile, but uses MPI-IO
        ! Opens file, advances record count, writes header etc.
        Use RA_MPI_BASE
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer, Intent(In) :: iter
        Integer, Intent(InOut) :: ierr
        Character*120 :: iterstring
        Character*120 :: filename
        Integer :: modcheck, imod, file_iter, next_iter, ibelong, icomp
        Integer :: buffsize, funit
        Integer :: mstatus(MPI_STATUS_SIZE)
        integer(kind=MPI_OFFSET_KIND) :: disp
        Logical :: create_file
        Logical :: file_exists

        create_file = .false.
        file_exists = .false.

        self%file_open    = .false.
        self%write_header = .false.

        ! Determine which file # this data belongs to
        modcheck = self%frequency*self%rec_per_file
        icomp = iter - (self%cache_size-1)*self%frequency
        imod = Mod(icomp,modcheck) 
        If (imod .eq. 0) Then
            ibelong = iter
        Else
            ibelong = icomp-imod+modcheck    ! Data belongs to file #ibelong
        Endif

        Write(iterstring,i_ofmt) ibelong
        filename = trim(local_file_path)//trim(self%file_prefix)//trim(iterstring)

        If ( (imod .eq. self%frequency) .or. (self%rec_per_file .eq. 1) ) Then   
            ! Time to begin a new file 
            create_file = .true.
        Else
            ! We are either:
            ! (a)  Appending to an existing file or
            ! (b)  The user has changed the output cadence in between checkpoints.
            !      In that case, we're "off-cycle" and need to create the file.

            ! This sort of filesystem inquiry should be done by only one rank
            If (self%orank .eq. 0) Then
                Inquire(File=filename, Exist=file_exists)
                If (.not. file_exists) create_file = .true.
            Endif

            Call MPI_Bcast(create_file, 1, MPI_LOGICAL, 0, self%ocomm, ierr)

        Endif

        If (create_file) Then

            If (self%master) Call stdout%print(' Creating Filename: '//trim(filename))
    	    Call MPI_FILE_OPEN(self%ocomm, filename, & 
                 MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                 MPI_INFO_NULL, funit, ierr) 

            self%file_unit = funit
            self%write_header = .true.            

        Else

    		Call MPI_FILE_OPEN(self%ocomm, filename, & 
                 MPI_MODE_RDWR, & 
                 MPI_INFO_NULL, funit, ierr) 


            self%file_unit = funit
            self%write_header = .false.

        Endif

        If (ierr .eq. 0) Then
            self%file_open = .true.
        Else
            Call stdout%print(' Unable to open file: '//trim(filename))
        Endif


        If (self%file_open) Then

            If (self%write_header) Then

                If (self%master) Then     


                    buffsize = 1
                    Call MPI_FILE_WRITE(self%file_unit, endian_tag, buffsize, &
                        & MPI_INTEGER, mstatus, ierr)

                    buffsize = 1
                    Call MPI_FILE_WRITE(self%file_unit, self%Output_Version, & 
                        buffsize, MPI_INTEGER, mstatus, ierr)

                    ! We write zero initially --
                    ! update nrec only after the data is actually written.
                    buffsize = 1
                    Call MPI_FILE_WRITE(self%file_unit, integer_zero, &
                        & buffsize, MPI_INTEGER, mstatus, ierr) 
                    
                Endif
                
                self%current_rec = 1+self%cc

            Else

                !Read the previous record number / advance record
                If (self%orank .eq. 0) Then

                    disp = 8
                    Call MPI_File_Seek(self%file_unit,disp,MPI_SEEK_SET,ierr)
                    Call MPI_FILE_READ(self%file_unit, self%current_rec, 1, &
                        & MPI_INTEGER, mstatus, ierr)

                Endif

                Call MPI_Bcast(self%current_rec, 1, MPI_INTEGER, 0, self%ocomm, ierr)

                self%current_rec = self%current_rec+1+self%cc

            Endif

        Endif

    End Subroutine OpenFile_Par

    Subroutine CloseFile(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        Call self%update_position()
        Write(self%file_unit,POS = 9)self%current_rec
        Close(self%file_unit)
    End Subroutine CloseFile

    Subroutine CloseFile_Par(self)
        USE RA_MPI_BASE
        Implicit None
        integer :: ierr, buffsize
        Integer :: mstatus(MPI_STATUS_SIZE)
        integer(kind=MPI_OFFSET_KIND) :: disp
        !Parallel File Close
        !Peforms the same task as closefile, but using MPI-IO
        Class(DiagnosticInfo) :: self
        disp = 8
        Call MPI_File_Seek(self%file_unit,disp,MPI_SEEK_SET,ierr)
        If (ierr .ne. 0) Then
            Write(6,*)'Error rewinding to header.  Error code: ', ierr, myid, self%file_prefix
        Endif
        ! Without this barrier,  the record update below can occurs before some ranks
        ! have read the original record count (leading to oversized, corrupted files).

        If (self%master) Then  

            buffsize = 1

            call MPI_FILE_WRITE(self%file_unit,self%current_rec , buffsize, MPI_INTEGER, & 
                mstatus, ierr) 
            If (ierr .ne. 0) Write(6,*)'Error writing to header.  Error code: ', ierr, myid, self%file_prefix
        Endif
        Call MPI_FILE_CLOSE(self%file_unit, ierr)
        If (ierr .ne. 0) Write(6,*)'Error closing file.  Error code: ',ierr, myid, self%file_prefix
    End Subroutine CloseFile_Par

    Subroutine Update_Position(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        self%file_position=ftell(self%file_unit)
    End Subroutine Update_Position

    Subroutine getq_now(self,yesno)
        Implicit None
        Logical, Intent(InOut) :: yesno
        Integer :: modcheck
        Class(DiagnosticInfo) :: self
        self%grab_this_q = .false.
        If(self%compute(current_qval) .eq. 1) Then
            modcheck = Mod(current_iteration,self%frequency)
            If (modcheck .eq. 0) Then
                yesno = .true.
                current_averaging_level = Max(current_averaging_level,self%avg_level)
                self%grab_this_q = .true.
            Endif
        Endif
    End Subroutine getq_now

    !/////////////////////////////////////////
    ! (Presently) Random Utility routines
    Subroutine ComputeEll0(inbuff,outbuff)
        Real*8, Intent(In) :: inbuff(1:,my_rmin:,my_theta_min:,1:)
        Real*8, Intent(InOut) :: outbuff(my_rmin:,1:)
        Real*8, Allocatable :: tmp_buffer(:,:)
        Integer :: bdims(1:4)
        Integer :: q,nq,r,t,p
        !Averages over theta and phi to get the spherically symmetric mean of all
        ! fields in inbuff at each radii
        ! inbuff is expected to be dimensioned as (1:nphi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,1:nfields)
        ! ** Note that this routine should be used sparingly because it requires 
        ! **   a collective operation (allreduce) across process rows
        ! ** One way to do this is to aggregate several fields into the inbuff when calling this routine


        bdims = shape(inbuff)
        nq = bdims(4)

        Allocate(tmp_buffer(my_rmin:my_rmax,1:nq))
        tmp_buffer(:,:) = 0.0d0
        outbuff(:,:) = 0.0d0

        ! Perform phi-integration and partial averaging in theta
        Do q = 1, nq
            Do t = my_theta_min, my_theta_max
                Do r = my_rmin, my_rmax
                    Do p = 1, nphi
                        tmp_buffer(r,q) = tmp_buffer(r,q)+inbuff(p,r,t,q) &
                            & *theta_integration_weights(t)
                    Enddo
                Enddo
            Enddo
        Enddo

        ! Turn phi-integration into an average
        tmp_buffer(:,:) = tmp_buffer(:,:)*over_nphi_double

        ! Complete the averaging process in theta
        Call DALLSUM2D(tmp_buffer, outbuff, pfi%rcomm)

        DeAllocate(tmp_buffer)

    End Subroutine ComputeEll0

    Subroutine IOComputeEll0(inbuff,outbuff)
        !Works exactly like computeEll0, but inbuff has already been averaged in phi
        Real*8, Intent(In) :: inbuff(my_rmin:,my_theta_min:,1:)
        Real*8, Intent(InOut) :: outbuff(my_rmin:,1:)
        Real*8, Allocatable :: tmp_buffer(:,:)
        Integer :: bdims(1:3)
        Integer :: q,nq,r,t,p
        !Averages over theta and phi to get the spherically symmetric mean of all
        ! fields in inbuff at each radii
        ! inbuff is expected to be dimensioned as (my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,1:nfields)
        ! ** Note that this routine should be used sparingly because it requires 
        ! **   a collective operation (allreduce) across process rows
        ! ** One way to do this is to aggregate several fields into the inbuff when calling this routine


        bdims = shape(inbuff)
        nq = bdims(3)

        Allocate(tmp_buffer(my_rmin:my_rmax,1:nq))
        tmp_buffer(:,:) = 0.0d0
        outbuff(:,:) = 0.0d0

        ! Perform partial averaging in theta
        Do q = 1, nq
            Do t = my_theta_min, my_theta_max
                Do r = my_rmin, my_rmax
                    tmp_buffer(r,q) = tmp_buffer(r,q)+inbuff(r,t,q) &
                        & *theta_integration_weights(t)
                Enddo
            Enddo
        Enddo

        !This next line was erroneous.  inbuff has already been averaged in phi
        ! Turn phi-integration into an average
        !tmp_buffer(:,:) = tmp_buffer(:,:)*over_nphi_double

        ! Complete the averaging process in theta
        Call DALLSUM2D(tmp_buffer, outbuff, pfi%rcomm)

        DeAllocate(tmp_buffer)

    End Subroutine IOComputeEll0

    Subroutine Compute_Radial_Average(inbuff,outbuff)
        Real*8, Intent(In) :: inbuff(my_rmin:,1:)
        Real*8, Intent(InOut) :: outbuff(1:)
        Real*8, Allocatable :: tmp_buffer(:)
        Integer :: bdims(1:2)
        Integer :: q,nq,r,t,p
        !Averages over radius for all fields contained in inbuff
        ! fields in inbuff at each radii
        ! inbuff is expected to be dimensioned as (1:nphi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,1:nfields)
        ! ** Note that this routine should be used sparingly because it requires 
        ! **   a collective operation (allreduce) across process rows
        ! ** One way to do this is to aggregate several fields into the inbuff when calling this routine


        bdims = shape(inbuff)
        nq = bdims(2)

        Allocate(tmp_buffer(1:nq))
        tmp_buffer(:) = 0.0d0
        outbuff(:) = 0.0d0

        ! Perform partial averaging in r
        Do q = 1, nq
            Do r = my_rmin, my_rmax
                tmp_buffer(q) = tmp_buffer(q)+inbuff(r,q) &
                    & *r_integration_weights(r)
            Enddo
        Enddo

        ! Complete the averaging process in theta
        Call DALLSUM1D(tmp_buffer, outbuff, pfi%ccomm)

        DeAllocate(tmp_buffer)

    End Subroutine Compute_Radial_Average

    Subroutine IOComputeM0(qty)
        Real*8, Intent(In) :: qty(1:,my_rmin:,my_theta_min:)
        Integer :: q,nq,r,t,p, ind
        !Averages over phi to get the azimuthally symmetric mean of all
        ! fields in inbuff at each radii and theta value.
        ! inbuff is expected to be dimensioned as (1:nphi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,my_theta_min:my_theta_max,1:nfields)

        ind = Shell_Averages%ind

        ! Perform phi-integration

        Do t = my_theta_min, my_theta_max
            Do r = my_rmin, my_rmax
                Do p = 1, nphi
                    IOm0_values(r,t,ind) = IOm0_values(r,t,ind)+qty(p,r,t)
                Enddo
                IOm0_values(r,t,ind) = IOm0_values(r,t,ind)*over_nphi_double ! Turn integration into an average
            Enddo
        Enddo

        Call Shell_Averages%AdvanceInd()

    End Subroutine IOComputeM0

    Subroutine ComputeM0(inbuff,outbuff)
        Real*8, Intent(In) :: inbuff(1:,my_rmin:,my_theta_min:,1:)
        Real*8, Intent(InOut) :: outbuff(my_rmin:,my_theta_min:,1:)

        Integer :: bdims(1:4)
        Integer :: q,nq,r,t,p
        !Averages over phi to get the azimuthally symmetric mean of all
        ! fields in inbuff at each radii and theta value.
        ! inbuff is expected to be dimensioned as (1:nphi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:nfields)
        ! outbuff is dimensioned as outbuff(my_r%min:my_rmax,my_theta_min:my_theta_max,1:nfields)

        bdims = shape(inbuff)
        nq = bdims(4)

        outbuff(:,:,:) = 0.0d0

        ! Perform phi-integration and partial averaging in theta
        Do q = 1, nq
            Do t = my_theta_min, my_theta_max
                Do r = my_rmin, my_rmax
                    Do p = 1, nphi
                        outbuff(r,t,q) = outbuff(r,t,q)+inbuff(p,r,t,q)
                    Enddo
                Enddo
            Enddo
        Enddo

        ! Turn phi-integration into an average
        outbuff = outbuff*over_nphi_double

    End Subroutine ComputeM0

    !/////////////////////////////////////////
    ! Cleanup Routines
    Subroutine CleanUp(self)
        Implicit None
        Class(DiagnosticInfo) :: self
        ! DeAllocates all allocated components of the info structure.
        ! Reinitializes all static components to their initial values (if they had one).
        self%values(1:nqmax) = -1 
        self%levels(1:nshellmax) = -1 
        self%compute(1:nqmax)= -1
        self%nq = 0
        self%nlevels = 0
        self%my_nlevels = 0
        self%file_unit = 15
        self%file_prefix = 'None'
        If (Allocated(self%oqvals))  DeAllocate(self%oqvals)

        self%frequency = reallybig 
        self%rec_per_file =1     
        self%current_rec = 1       
        self%file_header_size =0 
        self%file_record_size = 0 
        self%file_position = 1   
        self%avg_level = 0       

        self%ind = 1              
        self%begin_output = .false.
        self%mpi_tag = 1          

        self%output_version = -1 
        self%grab_this_q = .false.       

        self%master = .false.
    End Subroutine CleanUP

    Subroutine CleanUP_Spherical_IO
        Implicit None
        !DeAllocates all allocated module arrays
        !Re-initializes module variables to their default values (if any)
        current_averaging_level = 0
        current_qval = 0

        averaging_level(1:nqmax) = 0
        compute_q(1:nqmax) = 0
        shellavg_values(1:nqmax)=-1 
        globalavg_values(1:nqmax)=-1
        shellslice_values(1:nqmax) =-1
        shellslice_levels(1:nshellmax)=-1
        azavg_values(1:nqmax)=-1
        full3d_values(1:nqmax) = -1 
        shellspectra_values(1:nqmax)=-1
        shellspectra_levels(1:nshellmax)=-1

        globalavg_nrec = 1
        shellavg_nrec = 1
        azavg_nrec = 1
        shellslice_nrec =1
        shellspectra_nrec =1

        globalavg_frequency = reallybig
        shellavg_frequency = reallybig
        azavg_frequency = reallybig
        shellslice_frequency = reallybig
        shellspectra_frequency=reallybig

        full3d_frequency= reallybig
        local_file_path=''

        integer_zero = 0
        If (Allocated(qty)) DeAllocate(qty)
        If (Allocated(f_of_r_theta)) DeAllocate(f_of_r_theta)
        If (Allocated(f_of_r)) DeAllocate(f_of_r)
        If (Allocated(shellav_outputs)) DeAllocate(shellav_outputs)
        If (Allocated(globav_outputs)) DeAllocate(globav_outputs)
    
        i_ofmt = '(i8.8)'  ! These should never change during a run, but just in case...

        io_node = 0
        If (Allocated(theta_integration_weights)) DeAllocate(theta_integration_weights)
        If (Allocated(r_integration_weights))   DeAllocate(r_integration_weights)
        If (Allocated(radius)) DeAllocate(radius)
        If (Allocated(sintheta)) DeAllocate(sintheta)
        If (Allocated(costheta)) DeAllocate(costheta)


        Call Full_3D%cleanup
        Call Global_Averages%cleanup
        Call AZ_Averages%cleanup
        Call Shell_Averages%cleanup
        Call Shell_Slices%cleanup
        Call Shell_Spectra%cleanup
        Call Point_Probes%cleanup

        !Call Temp_IO%cleanup

        !Call spectra_buffer%cleanup

    End Subroutine CleanUP_Spherical_IO

    SUBROUTINE nrm_to_index(nrm_coords, coord_arr, indices, rev_inds)
        ! Converts a list of normalized coordinations (nrm_coords)
        ! into grid indices (indices) based on the indicated
        ! coordinate axes (coord_arr; e.g. radius)
        !
        ! nrm_coords and coord_arr are assumed to be in ascending order
        ! If necessary, the indices can be reversed using rev_inds.  
        ! This is mostly here to deal with the "reversed" radius array.
        IMPLICIT NONE
        REAL*8, Intent(In)  :: nrm_coords(:)
        REAL*8, Intent(In)  :: coord_arr(:)
        INTEGER, Intent(Out) :: indices(:)
        INTEGER, ALLOCATABLE :: ind_copy(:)
        LOGICAL, INTENT(In), Optional :: rev_inds
        INTEGER :: nnorm_coord, nbase_coord, numi
        INTEGER :: i,j,k, ind, last_index
        REAL*8 :: cdel, del1,del2, ccheck
        REAL*8 :: nrmc, cmin, tolchk = 2.0d0

        ! User-specified normalized coordinates are assumed to
        ! be in ascending order, but rayleigh's grids are striped in
        ! in descending order.  We need to reverse the user-defined
        ! normalized coordinates in order for the coordinate-to-index
        ! mapping to work properly.
        !
        ! This also means we need to account for the use of negative
        ! number to indicate a range in coordinates (hence the if/else below)

        cmin = MINVAL(coord_arr)
        cdel = MAXVAL(coord_arr)-cmin 

        indices(:) = -1   ! Set to the default index value
        nnorm_coord = SIZE(nrm_coords)
        nbase_coord = SIZE(coord_arr)

        last_index = 0  ! The user might inadvertantly specify redundant coordinates.
                        ! This can happen if resolution is low.
                        ! We make sure that no redundant coordinates are output.
        numi = 0 ! number of valid coordinates specified
        DO i = 1, nnorm_coord

            IF (i .ne. 1) last_index = abs(indices(i-1))

            ccheck = nrm_coords(i)  

            IF (ABS(ccheck) .lt. tolchk) THEN
                numi = numi+1
                ! This is a valid normalized-coordinate value.
                ! Translate it into a grid-based index.
                ! Find the index within coord_arr that
                ! corresponds to the
                indices(i) = 1
                IF (ccheck .lt. 0) indices(i) = -1
                nrmc = (coord_arr(1)-cmin)/cdel
                del1 = ABS(nrmc-abs(ccheck))

                DO j = 2, nbase_coord
                    nrmc = (coord_arr(j)-cmin)/cdel
                    del2 = ABS(nrmc-abs(ccheck))
                    If ( (del2 .lt. del1) .and. (j .ne. last_index) ) THEN
                        indices(i) = j
                        If (ccheck .lt. 0) indices(i) = -j
                        del1 = del2
                    ENDIF
                ENDDO

            ENDIF


        ENDDO

        IF (PRESENT(rev_inds)) THEN
           IF (rev_inds) THEN ! in case someone specified rev_inds = .false. ...
                i = 1
                DO WHILE(i .le. numi)
                    IF (indices(i) .lt. 0) THEN
                        indices(i) = -indices(i)
                        indices(i-1) = -indices(i-1)
                        !i = i+1
                    ENDIF
                    i = i+1
                ENDDO
                ALLOCATE(ind_copy(1:numi))
                ind_copy(1:numi) = indices(1:numi)

                ! Second pass:  rearrange indices in ascending order
                DO i = 1, numi
                    indices(i) = ind_copy(numi-i+1)
                ENDDO

                DEALLOCATE(ind_copy)
           ENDIF 
        ENDIF

    END SUBROUTINE nrm_to_index

    Subroutine Parse_Inds(indices_inout, indcount)
        IMPLICIT NONE
        INTEGER, INTENT(InOut) :: indices_inout(:)
        INTEGER, INTENT(InOut), Optional :: indcount
        INTEGER, ALLOCATABLE :: indices_out(:)
        INTEGER :: i, j,ind,ni, jmin, jmax
        !Converts index lists of the form [1,-4] to [1,2,3,4]
        i = 1
        ind = 1
        ni = size(indices_inout)
        ALLOCATE(indices_out(1:ni))
        indices_out(:) = -1
        DO WHILE( (indices_inout(i) .gt. -1) .and. (i .le. ni))
            indices_out(ind) = indices_inout(i)
            ind = ind+1
            IF ( indices_inout(i+1) .lt. -1) THEN
                ! User has specified a sub-range in this coordinate
                jmax = -indices_inout(i+1) 
                jmin =  indices_inout(i)+1
                Do j = jmin,jmax
                    indices_out(ind) = j
                    ind = ind+1
                ENDDO
                i = i + 1  ! We increment an extra time here to skip the next negative number
            ENDIF
            i = i + 1
        ENDDO
        indices_inout(1:ni) = indices_out(1:ni)
        DEALLOCATE(indices_out)
        IF (present(indcount)) indcount = ind-1        
    End Subroutine Parse_Inds

    SUBROUTINE PROCESS_COORDINATES()
        IMPLICIT NONE
        INTEGER :: i
        Real*8, Allocatable :: tmp_theta(:), tmp_phi(:)
        ALLOCATE(tmp_theta(1:ntheta), tmp_phi(1:nphi))


        ! We want theta instead of costheta
        DO i = 1, ntheta
            tmp_theta(i) = ACOS(costheta(i))
        ENDDO

        DO i = 1, nphi
            tmp_phi(i) = (i-1)*two_pi/DBLE(nphi)
        ENDDO

        CALL INTERPRET_INDICES(      point_probe_r_nrm, radius   , point_probe_r,      revg =.true.)
        CALL INTERPRET_INDICES(  point_probe_theta_nrm, tmp_theta, point_probe_theta,  revg=.true.)
        CALL INTERPRET_INDICES(    point_probe_phi_nrm, tmp_phi  , point_probe_phi)
        CALL INTERPRET_INDICES(  shellslice_levels_nrm, radius   , shellslice_levels,  revg=.true.)
        CALL INTERPRET_INDICES(shellspectra_levels_nrm, radius   , shellspectra_levels,revg=.true.)
        CALL INTERPRET_INDICES(sph_mode_levels_nrm, radius   , sph_mode_levels,revg=.true.)
        CALL INTERPRET_INDICES( meridional_indices_nrm, tmp_phi  , meridional_indices, revg=.true.)

        DeALLOCATE(tmp_theta,tmp_phi)
        
        ! Parse the SPH_MODE_ELL list:
        Call Parse_inds(sph_mode_ell, SPH_MODE_Nell)
        If (SPH_Mode_nell .gt. 0) Then
            Do i =1, sph_mode_nell
                sph_mode_nmode = sph_mode_nmode+(sph_mode_ell(i)+1)
            Enddo
        ENDIF

    END SUBROUTINE PROCESS_COORDINATES

    SUBROUTINE Interpret_Indices(indices_nrm, coord_grid, indices, revg)
        IMPLICIT NONE
        INTEGER, INTENT(InOut) :: indices(:)
        REAL*8, INTENT(InOut) :: coord_grid(:), indices_nrm(:)
        LOGICAL, INTENT(In), Optional :: revg
        LOGICAL :: reverse_grid 
        IF (present(revg)) THEN
            IF (revg) reverse_grid = .true.
        ELSE
            reverse_grid = .false.
        ENDIF

        ! If needed, convert normalized coordinates to indices
        ! e.g., map [0.1, 0.25, 0.5] to [4, 15, 32]  
        IF (maxval(indices_nrm) .gt. -2.9d0) THEN
            
            CALL nrm_to_index(indices_nrm, coord_grid, indices, rev_inds =reverse_grid)
        ENDIF

        ! Next, interpret any range shorthand used.
        ! e.g., convert [1,-4,8] to [1,2,3,4,8]
        Call Parse_Inds(indices)
    END SUBROUTINE Interpret_Indices

End Module Spherical_IO
