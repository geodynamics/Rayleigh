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
#ifdef INTEL_COMPILER 
    USE IFPORT
#endif
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
    Integer, Parameter :: sphmode_version =3
    INTEGER, PARAMETER :: probe_version = 1
    Integer, Parameter :: full3d_version = 3    !currently unused
    Type, Public :: DiagnosticInfo
        ! Need to see if we can make these allocatable, but for now..
        ! Each instance of this class has two static arrays used for reading in namelist input
        ! These arrays are large enough to hold nqmax and nshellmax values, but 
        ! typically only a small fraction of that amount will be specified at input.
        Integer :: values(1:nqmax) = -1 ! The list of values specified in an input namelist
        Integer :: levels(1:nshellmax) = -1 ! The radial indices output (shell slices, spectra, and histograms only)
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
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Additional Organizational Arrays are required for managing shell-slice like outputs
        Integer, Allocatable :: my_shell_levs(:), have_shell(:), my_shell_ind(:)
        Integer, Allocatable :: shell_r_ids(:), nshells_at_rid(:)
        Integer :: nshell_r_ids

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Variables used for management of SPH_Mode_Samples output
        INTEGER, Allocatable :: l_values(:), m_values(:)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Variables used for management of Meridional Slices
        INTEGER, Allocatable :: phi_indices(:)
        INTEGER :: nphi_indices

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !  Variables used for cached output (currently only used for point probes)
        INTEGER :: cache_size = 1
        INTEGER :: cc = 0  ! cache counter
        INTEGER, ALLOCATABLE :: iter_save(:)
        REAL*8 , ALLOCATABLE :: time_save(:)

        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Point-probe Variables
        INTEGER, ALLOCATABLE :: probe_p_global(:)  !phi-indces (phi is in-process at I/O time) 
        INTEGER, ALLOCATABLE :: probe_r_global(:), probe_r_local(:) ! global/local r-indices
        INTEGER, ALLOCATABLE :: probe_t_global(:), probe_t_local(:) ! " " theta-indices
        INTEGER :: probe_nr_local, probe_nr_global  !number of local/global r-indices
        INTEGER :: probe_nt_local, probe_nt_global, probe_np_global !" " theta-indices & phi-indices
        INTEGER, ALLOCATABLE :: probe_nr_atrank(:)  ! Number of radial indices held by each column rank
        INTEGER, ALLOCATABLE :: probe_nt_atrank(:)  ! Number of theta indices held by each theta rank
        INTEGER, ALLOCATABLE :: npts_at_colrank(:)  ! Number of points contained on a given row
        INTEGER, ALLOCATABLE :: npts_at_rowrank(:)  ! Number of points contained at a given rank within a row
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        !Communicatory Info for parallel writing (if used)
        Integer :: ocomm, orank, onp
        Logical :: master = .false.
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! Methods
        Contains
        Procedure :: Init => Initialize_Diagnostic_Info
        Procedure :: AdvanceInd
        Procedure :: AdvanceCC
        Procedure :: reset => diagnostic_output_reset
        Procedure :: Shell_Balance
        Procedure :: Scattered_Balance
        Procedure :: OpenFile
        Procedure :: OpenFile_Par
        Procedure :: Set_File_Info
        Procedure :: CloseFile
        Procedure :: CloseFile_Par
        Procedure :: update_position
        Procedure :: getq_now
        Procedure :: init_ocomm
        Procedure :: CleanUp
        !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ! 


    End Type DiagnosticInfo

    Type(DiagnosticInfo) :: Shell_Averages, Shell_Slices, Global_Averages, AZ_Averages, Full_3D, Shell_Spectra
    Type(DiagnosticInfo) :: Equatorial_Slices, Meridional_Slices, SPH_Mode_Samples, Point_Probes


    Integer :: current_averaging_level = 0
    Integer :: current_qval = 0

    Integer :: averaging_level(1:nqmax) = 0, compute_q(1:nqmax) = 0
    Integer :: shellavg_values(1:nqmax)=-1, globalavg_values(1:nqmax)=-1
    Integer :: shellslice_values(1:nqmax) =-1, shellslice_levels(1:nshellmax)=-1, azavg_values(1:nqmax)=-1
    Integer :: full3d_values(1:nqmax) = -1, shellspectra_values(1:nqmax)=-1, shellspectra_levels(1:nshellmax)=-1
    Integer :: histo_values(1:nqmax) = -1, histo_levels(1:nshellmax)=-1
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
    Real*8, Private, Allocatable :: circumference(:,:), qty(:,:,:), f_of_r_theta(:,:)
    Real*8, Private, Allocatable :: azav_outputs(:,:,:), f_of_r(:), rdtheta_total(:)
    Real*8, Private, Allocatable :: shellav_outputs(:,:,:), globav_outputs(:), shell_slice_outputs(:,:,:,:)
    Real*8, Private, Allocatable :: meridional_outputs(:,:,:,:), probe_outputs(:,:,:)

    Type(SphericalBuffer) :: spectra_buffer, sph_sample_buffer
    Real*8, Private :: da_total, int_vol, int_dphi, int_rsquared_dr, int_sintheta_dtheta
    Real*8, Private, Allocatable :: sintheta_dtheta(:), rsquared_dr(:)
    Character*6, Public :: i_ofmt = '(i8.8)', i_pfmt = '(i5.5)'
    integer :: output_ireq(1), output_status(1)

    !//////////////////////////////////////////////////////////////////////
    ! *****           Equatorial Slice Output Variables
    REAL*8, ALLOCATABLE :: equslice_outputs(:,:,:)
    INTEGER :: th1_owner = -1, th2_owner = -1 ! Angular ranks of the owners of ntheta/2 and ntheta/2+1 
    INTEGER :: my_nth_owned = -1


    !////////////////////////////////////////////////////////////////////
    Logical :: start_az_average, start_shell_average  !, start_shell_slice
    Logical :: start_pdf, start_global_average, start_spectra

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

        Allocate(f_of_r_theta(my_rmin:my_rmax,my_theta_min:my_theta_max))
        Allocate(f_of_r(my_rmin:my_rmax))

        num_avg_store = shell_averages%nq
        Allocate(IOm0_values(my_rmin:my_rmax,my_theta_min:my_theta_max,1:num_avg_store))
        Allocate(IOell0_values(my_rmin:my_rmax,1:num_avg_store))
        IOm0_values(:,:,:) = 0.0d0
        IOell0_values(:,:) = 0.0d0
    End Subroutine Begin_Outputting



	Subroutine Initialize_Spherical_IO(rad_in,sintheta_in, rw_in, tw_in, costheta_in,file_path)
		Implicit None
		Integer :: k, fcount(3,2), ntot, fcnt, master_rank, i
        Integer :: thmax, thmin, nth1, nth2, nn
		Real*8, Intent(In) :: rad_in(:), sintheta_in(:), rw_in(:), tw_in(:), costheta_in(:)
        Character*120 :: fdir
        Character*120, Intent(In) :: file_path
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

        !Numbers here are the mpi_tags used in communication for each output
        !In theory they can be the same, but it's probably a good idea to keep them unique
        Call            Full_3D%Init(averaging_level,compute_q,myid, &
            & 54,values = full3d_values)

        Call    Global_Averages%Init(averaging_level,compute_q,myid, &
            & 55,avg_level = 3,values = globalavg_values)

		Call        AZ_Averages%Init(averaging_level,compute_q,myid, &
            & 56,avg_level = 1,values = azavg_values)

        Call     Shell_Averages%Init(averaging_level,compute_q,myid, &
            & 57,avg_level = 2,values = shellavg_values)

        Call       Shell_Slices%Init(averaging_level,compute_q,myid, &
            & 58,values = shellslice_values, levels = shellslice_levels)

        Call       Shell_Spectra%Init(averaging_level,compute_q,myid, &
            & 59,values = shellspectra_values, levels = shellspectra_levels)

        Call       Equatorial_Slices%Init(averaging_level,compute_q,myid, &
            & 60,values = equatorial_values)

        Call       Meridional_Slices%Init(averaging_level,compute_q,myid, &
            & 61,values = meridional_values, phi_inds = meridional_indices)

        Call       SPH_Mode_Samples%Init(averaging_level,compute_q,myid, &
            & 62,values = sph_mode_values, levels = sph_mode_levels)


        Call       Point_Probes%Init(averaging_level,compute_q,myid, &
            & 63,values = point_probe_values, cache_size = point_probe_cache_size)


        !Outputs involve saving and communicating partial shell slices (e.g. Shell_Slices or spectra)
        !require an additional initialization step to load-balance the shells
        Call Shell_Slices%Shell_Balance()
        Call Shell_Spectra%Shell_Balance()
        Call SPH_Mode_Samples%Shell_Balance()
        Call Point_Probes%Scattered_Balance(point_probe_r, point_probe_theta, &
            & point_probe_phi)
        if (my_row_rank .eq. 0) Then

            If (Shell_Slices%nshell_r_ids .gt. 0) Then
                master_rank = shell_slices%shell_r_ids(1)
                Call Shell_Slices%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,master_rank) ! For parallel IO
            Endif
            Call AZ_Averages%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,0) ! 0 handles file headers etc. for AZ Average output
            Call Equatorial_Slices%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,0)
            Call Meridional_Slices%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,0)

            If (Shell_Spectra%nshell_r_ids .gt. 0) Then
                master_rank = shell_spectra%shell_r_ids(1)
                Call Shell_Spectra%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,master_rank) 
            Endif

            If (SPH_Mode_Samples%nshell_r_ids .gt. 0) Then
                master_rank = SPH_Mode_Samples%shell_r_ids(1)
                Call SPH_Mode_Samples%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,master_rank) 
            Endif
            If (maxval(point_probes%npts_at_colrank) .gt. 0) THEN

                i = 0    
                master_rank = -1
                Do While (master_rank .eq. -1) 
                    if (point_probes%npts_at_colrank(i) .gt. 0) master_rank = i
                    i = i+1
                Enddo

                Call Point_Probes%init_ocomm(pfi%ccomm%comm,nproc1,my_column_rank,master_rank)
            ENDIF
        Endif
        
        If (Shell_Spectra%nlevels .gt. 0) Then
            !Shell spectra require an additional step
            !Initialize the buffer object that we use for transposing spectra
            !Some row ranks might have no buffer initialized.
            ntot = Shell_Spectra%nq*Shell_Spectra%my_nlevels
            fcnt = ntot/my_nr
            k = Mod(ntot,my_nr)
            if (k .gt. 0) fcnt = fcnt+1
            If (fcnt .gt. 0) Then
                fcount(:,:) = fcnt
           		Call spectra_buffer%init(field_count = fcount, config = 'p3b')		
            Endif	
        Endif 

        If (SPH_Mode_Samples%nlevels .gt. 0) Then
            !Similarly so for SPH_Mode_Samples
            ntot = SPH_Mode_Samples%nq*SPH_Mode_Samples%my_nlevels
            fcnt = ntot/my_nr
            k = Mod(ntot,my_nr)
            if (k .gt. 0) fcnt = fcnt+1
            If (fcnt .gt. 0) Then
                fcount(:,:) = fcnt
           		Call sph_sample_buffer%init(field_count = fcount, config = 'p3b')		
            Endif	
        Endif 



        !/////////////////////////////////////////////////
        !Some BookKeeping for Equatorial Slices
        nth1 = ntheta/2   ! These two points bracket the equator
        nth2 = ntheta/2+1
        Do nn = 0, nproc2-1

			thmax = pfi%all_2p(nn)%max
			thmin = pfi%all_2p(nn)%min

            If ( (nth1 .ge. thmin) .and. (nth1 .le. thmax) ) Then
                th1_owner = nn
            Endif
            If ( (nth2 .ge. thmin) .and. (nth2 .le. thmax) ) Then
                th2_owner = nn
            Endif

        Enddo
        my_nth_owned = 0
        If (my_row_rank .eq. th1_owner) my_nth_owned = 1
        If (my_row_rank .eq. th2_owner) my_nth_owned = my_nth_owned+1
        !//////////////////////////////////////////////////

        !Next, provide the file directory, frequency info, output versions etc.

        ! Global Averages
        fdir = 'G_Avgs/'
        Call Global_Averages%set_file_info(globalavg_version,globalavg_nrec,globalavg_frequency,fdir)    

        ! Shell Averages
        fdir = 'Shell_Avgs/'
        Call Shell_Averages%set_file_info(shellavg_version,shellavg_nrec,shellavg_frequency,fdir)    

        ! Azimuthal Averages
        fdir = 'AZ_Avgs/'
        Call AZ_Averages%set_file_info(azavg_version,azavg_nrec,azavg_frequency,fdir)    

        ! Shell Slices
        fdir = 'Shell_Slices/'
        Call Shell_Slices%set_file_info(shellslice_version,shellslice_nrec,shellslice_frequency,fdir)   

        ! Shell Spectra
        fdir = 'Shell_Spectra/'
        Call Shell_Spectra%set_file_info(shellspectra_version,shellspectra_nrec,shellspectra_frequency,fdir) 

        ! Equatorial Slices
        fdir = 'Equatorial_Slices/'
        Call Equatorial_Slices%set_file_info(equslice_version,equatorial_nrec,equatorial_frequency,fdir) 

        ! Meridional Slices
        fdir = 'Meridional_Slices/'
        Call Meridional_Slices%set_file_info(meridslice_version,meridional_nrec,meridional_frequency,fdir) 

        fdir = 'SPH_Modes/'
        Call SPH_Mode_Samples%set_file_info(sphmode_version,sph_mode_nrec,sph_mode_frequency,fdir) 

        !Point-wise Probes
        fdir = 'Point_Probes/'
        Call Point_Probes%set_file_info(probe_version,point_probe_nrec, &
                point_probe_frequency,fdir) 

        ! Full 3D (special because it only cares about the frequency, not nrec)
        fdir = 'Spherical_3D/'
        Call Full_3D%set_file_info(full3d_version,shellslice_nrec,full3d_frequency,fdir)    

   End Subroutine Initialize_Spherical_IO

    Subroutine Get_Meridional_Slice(qty)
        Implicit None
        REAL*8, INTENT(IN) :: qty(:,my_rmin:,my_theta_min:)
        INTEGER :: nphi_grab, mer_ind
        INTEGER :: tind, pind, rind
        nphi_grab = Meridional_Slices%nphi_indices
        IF (Meridional_Slices%begin_output) THEN
            ! Stripe with theta-index slowest, so data is ordered for row-communication
			Allocate(meridional_outputs(my_rmin:my_rmax,1:nphi_grab,&
                1:Meridional_Slices%nq, my_theta_min:my_theta_max))
        ENDIF


        mer_ind = Meridional_Slices%ind
        Do tind = my_theta_min, my_theta_max
            DO pind = 1,nphi_grab
                DO rind = my_rmin, my_rmax
                    meridional_outputs(rind,pind,mer_ind,tind) = &
                        qty(meridional_slices%phi_indices(pind),rind,tind)
                ENDDO
            ENDDO
        ENDDO    
        Meridional_Slices%oqvals(mer_ind) = current_qval
        Call Meridional_Slices%AdvanceInd()

    End Subroutine Get_Meridional_Slice

	Subroutine Write_Meridional_Slices(this_iter,simtime)
        USE RA_MPI_BASE
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:), all_slices(:,:,:,:)
		Integer :: responsible, current_rec, s_start, s_end, this_rid
		Integer :: i, j, k,qq, p, sizecheck, ii
		Integer :: n, nn, this_nshell, nq_merid, merid_tag,nphi_grab
		Integer :: your_theta_min, your_theta_max, your_ntheta
		Integer :: nelem, buffsize
        Integer :: file_pos, funit, error, dims(1:4)
        Integer :: inds(4), nirq,sirq
        Integer, Allocatable :: rirqs(:)
        
        integer :: ierr, rcount
		integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp, qdisp, full_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
        sizecheck = sizeof(disp)
        if (sizecheck .lt. 8) Then
            if (myid .eq. 0) Then
            Write(6,*)"Warning, MPI_OFFSET_KIND is less than 8 bytes on your system."
            Write(6,*)"Your size (in bytes) is: ", sizecheck
            Write(6,*)"A size of 4 bytes means that shell slices files are effectively limited to 2 GB in size."
            Endif
        endif 


		responsible = 0
        nq_merid     = Meridional_Slices%nq
        merid_tag    = Meridional_Slices%mpi_tag
        nphi_grab    = Meridional_Slices%nphi_indices
        funit        = Meridional_Slices%file_unit
        

		If (my_row_rank .eq. 0) Then
                responsible = 1
        Endif

        
		If (responsible .eq. 1) Then
			! Rank 0 in reach row receives from all other row members
		
			Allocate(all_slices(my_rmin:my_rmax,1:nphi_grab,1:nq_merid, 1:ntheta))
			all_slices(:,:,:,:) = 0.0d0

            nirq = nproc2-1
            Allocate(rirqs(1:nirq))

            Do nn = 1, nproc2-1

				your_ntheta    = pfi%all_2p(nn)%delta
				your_theta_min = pfi%all_2p(nn)%min
				your_theta_max = pfi%all_2p(nn)%max
                inds(1) = 1
                inds(2) = 1
                inds(3) = 1
                inds(4) = your_theta_min
                nelem = your_ntheta*my_nr*nq_merid*nphi_grab

                Call IReceive(all_slices, rirqs(nn),n_elements = nelem, &
                            &  source= nn,tag = merid_tag, grp = pfi%rcomm,indstart = inds)				
			Enddo

            all_slices(my_rmin:my_rmax, 1:nphi_grab, 1:nq_merid,my_theta_min:my_theta_max) = &
               meridional_outputs(my_rmin:my_rmax, 1:nphi_grab, 1:nq_merid,my_theta_min:my_theta_max)

            Call IWaitAll(nirq, rirqs)
            DeAllocate(rirqs)
		Else
			!  Rest of the row sends to process 0 within the row
            inds(1) = 1 !my_rmin
            inds(2) = 1
            inds(3) = 1
            inds(4) = 1 
            nelem = my_nr*my_ntheta*nq_merid*nphi_grab
            Call ISend(meridional_outputs, sirq,n_elements = nelem, dest = 0, tag = merid_tag, & 
                grp = pfi%rcomm, indstart = inds)
            Call IWait(sirq)
		Endif
        If (allocated(meridional_outputs)) DeAllocate(meridional_outputs)

        ! Communication is complete.  Now we open the file using MPI-IO
        

      

        If (responsible .eq. 1) Then   
            Call Meridional_Slices%OpenFile_Par(this_iter, error)
            current_rec = Meridional_Slices%current_rec
            funit = Meridional_Slices%file_unit
            ! before we do anything else, we need to restripe the data so
            ! that q-index is the slowest

            Allocate(buff(1:nphi_grab, 1:ntheta,my_rmin:my_rmax, 1:nq_merid))
            Do k = 1, ntheta
                Do j = 1, nq_merid
                    Do ii = 1, nphi_grab
                    Do i = my_rmin, my_rmax    
                        buff(ii,k,i,j) = all_slices(i,ii,j,k)
                    Enddo
                    Enddo
                Enddo
            Enddo        
            DeAllocate(all_slices)
            

            If ( (my_column_rank .eq. 0) .and. (Meridional_Slices%write_header) ) Then
                If (Meridional_Slices%file_open) Then
                    ! Rank 0 in column and row writes the header
                    dims(1) =  nr
                    dims(2) =  ntheta
                    dims(3) =  nphi_grab
                    dims(4) =  nq_merid
                    buffsize = 4
                    Call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 

                    buffsize = nq_merid
                    Call MPI_FILE_WRITE(funit,Meridional_Slices%oqvals, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 

                    buffsize = nr
                    Call MPI_FILE_WRITE(funit, radius, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 

                    buffsize = ntheta
                    Call MPI_FILE_WRITE(funit, costheta, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 

                    buffsize = nphi_grab
                    Call MPI_FILE_WRITE(funit, Meridional_Slices%phi_indices, buffsize, &
                        MPI_INTEGER, mstatus, ierr) 
                Endif

            Endif

            hdisp = 28 ! dimensions+endian+version+record count
            hdisp = hdisp+nq_merid*4 ! nq
            hdisp = hdisp+nr*8  ! The radius array
            hdisp = hdisp+ ntheta*8  ! costheta
            hdisp = hdisp+ nphi_grab*4  ! phi indices

            qdisp = ntheta*nr*nphi_grab*8
            full_disp = qdisp*nq_merid+12  ! 12 is for the simtime+iteration at the end
            disp = hdisp+full_disp*(current_rec-1)
            
            buffsize = my_nr*ntheta*nphi_grab
            ! The file is striped with time step slowest, followed by q

            my_rdisp = (my_rmin-1)*ntheta*nphi_grab*8

            If (Meridional_Slices%file_open) Then

                Do i = 1, nq_merid
                    new_disp = disp+qdisp*(i-1)+my_rdisp                
                    Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                    
                    Call MPI_FILE_WRITE(funit, buff(1,1,my_rmin,i), buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                Enddo

                disp = hdisp+full_disp*current_rec
                disp = disp-12
                Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)

                If (my_column_rank .eq. 0) Then
                    buffsize = 1
                    Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                    Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                           MPI_INTEGER, mstatus, ierr)
                Endif

            Endif

			DeAllocate(buff)

            Call Meridional_Slices%Closefile_Par()

        Endif  ! Responsible

	End Subroutine Write_Meridional_slices


    Subroutine Get_Point_Probes(qty)
        Implicit None
        REAL*8, INTENT(IN) :: qty(:,my_rmin:,my_theta_min:)
        INTEGER :: q_ind, npts, ind, cache_ind
        INTEGER :: i,j,k
        INTEGER :: rind, tind
 
        q_ind = Point_Probes%ind
        npts = Point_Probes%npts_at_rowrank(my_row_rank)
        IF (Point_Probes%begin_output) THEN
            If ((npts .gt. 0) .and. (Point_Probes%cc .eq. 0)) THEN
                ALLOCATE( probe_outputs( &
                     1:Point_Probes%nq, 1:Point_Probes%cache_size,1:npts) )
                probe_outputs(:,:,:) = 0.0d0
            ENDIF
        ENDIF

        IF (npts .gt. 0) THEN
            cache_ind = Point_Probes%cc+1
            ind = 1
            ! We stripe phi,r,theta (just as qty's indexing)

            DO k = 1, Point_Probes%probe_nt_local
                tind = Point_Probes%probe_t_local(k)
                DO j = 1, Point_Probes%probe_nr_local
                    rind = Point_Probes%probe_r_local(j)
                    DO i = 1, Point_Probes%probe_np_global
                        probe_outputs(q_ind,cache_ind,ind) = &
                            qty(Point_Probes%probe_p_global(i),rind,tind)
                        ind = ind+1
                    ENDDO
                ENDDO
            ENDDO
        ENDIF

        Point_Probes%oqvals(q_ind) = current_qval
        Call Point_Probes%AdvanceInd()

    End Subroutine Get_Point_Probes


	Subroutine Write_Point_Probes(this_iter,simtime)
        USE RA_MPI_BASE
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter

		Real*8, Allocatable :: buff(:,:,:), row_probes(:,:,:), slice(:,:,:)

		INTEGER(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_pdisp, new_disp, qdisp, full_disp, tdisp

		INTEGER :: responsible, current_shell, s_start, s_end, this_rid
		INTEGER :: i, j, k,qq, p, sizecheck, t
		INTEGER :: n, nn, this_nshell, nq, probe_tag
		INTEGER :: your_theta_min, your_theta_max, your_ntheta, your_id
		INTEGER :: nelem, buffsize, sirq, nrirqs, inds(1:3), buffsize2
        INTEGER :: file_pos, funit, error, dims(1:4), first_shell_rank        
        INTEGER :: ierr, pcount, ichunk, ii, jj, ind, nt_row, rind, rslab_size
		INTEGER :: mstatus(MPI_STATUS_SIZE)
        INTEGER :: npts_this_row, this_nr, nphi_probe, irqc, ncache, npts
        INTEGER :: probe_nr, probe_nt, probe_np, nvals, this_nt, tind
        Real*8, Allocatable :: probe_vals(:)
        INTEGER, Allocatable :: level_inds(:), rirqs(:)


        nq      = Point_Probes%nq
        ncache  = Point_Probes%cache_size
        probe_tag = Point_Probes%mpi_tag
        funit   = Point_Probes%file_unit

        probe_nr = Point_Probes%probe_nr_global
        probe_nt = Point_Probes%probe_nt_global
        probe_np = Point_Probes%probe_np_global

        ! Initially, we carry out communication down-row, but only if our row has points
		responsible = 0
        npts_this_row = Point_Probes%npts_at_colrank(my_column_rank)
		If (my_row_rank .eq. 0) Then
            If (npts_this_row .gt. 0) Then
                responsible = 1
            Endif
        Endif

        !////////////////////////////////




		If (responsible .eq. 1) Then
            this_nr = Point_Probes%probe_nr_atrank(my_column_rank)
            nphi_probe = Point_Probes%probe_np_global

            ALLOCATE( row_probes(1:npts_this_row, 1:nq, ncache) )

            ALLOCATE(buff(1:nq, 1:ncache, 1:npts_this_row ))            

            row_probes(:,:,:) = 0.0d0
                  buff(:,:,:) = 0.0d0


            ! Post Ireceives
            nrirqs = nproc2-1
            Allocate(rirqs(1:nrirqs))          
            irqc = 0  !irq counter -- everyone in the row does not necessarily have points
            inds(3) = 1
            Do nn = 1, nproc2-1
                npts = Point_Probes%npts_at_rowrank(nn)
                inds(3) = inds(3)+Point_Probes%npts_at_rowrank(nn-1) ! Keep this out of the conditional

                IF (npts .gt. 0) THEN
                    irqc = irqc+1
				    your_id = nn

                    inds(1) = 1
                    inds(2) = 1

				    nelem = npts*ncache*nq

             		Call IReceive(buff, rirqs(irqc),n_elements = nelem, & 
                        source= your_id, tag=probe_tag, grp = pfi%rcomm, indstart = inds)
                ENDIF
            Enddo  
            ! Stripe my own data
            If (Point_Probes%npts_at_rowrank(my_row_rank) .gt. 0) THEN
                buff(:,:,1:Point_Probes%npts_at_rowrank(my_row_rank)) = probe_outputs(:,:,:)
            ENDIF
            Call IWaitAll(irqc,rirqs(1:irqc))  ! wait on the sends to come through



            !At this point, buff is striped phi,r,theta.  We need to stripe it phi,theta,r
            Allocate(slice(1:nq,1:ncache,1:nphi_probe))
            nt_row = SUM(Point_Probes%probe_nt_atrank)
            rslab_size = nt_row*nphi_probe
            !write(6,*)'Check: ', nt_row, rslab_size, this_nr
            ichunk = 1
            tind = 1
            DO nn = 0, nproc2-1
                this_nt = Point_Probes%probe_nt_atrank(nn)
                IF (this_nt .gt. 0) THEN
                    !Iterate through the receive buffer, grabbing one r-theta combination at a time
                    Do j = 1, this_nt
                        rind = 0
                        Do i = 1, this_nr
                            slice(:,:,1:nphi_probe) = buff(:,:,ichunk:ichunk+nphi_probe-1)
                            ichunk = ichunk+nphi_probe

                            ind = rind+tind
                            Do jj = 1, ncache
                                Do ii = 1, nq
                                    row_probes(ind:ind+nphi_probe-1,ii,jj) = slice(ii,jj,1:nphi_probe)
                                Enddo
                            Enddo

                            rind = rind+rslab_size
                        Enddo
                        tind = tind+nphi_probe
                    ENDDO
                ENDIF
            ENDDO

            If (allocated(probe_outputs)) DeAllocate(probe_outputs)
            DeAllocate(slice)
            DeAllocate(buff)
		Else
			!  Non responsible nodes send their info

            nelem = Point_Probes%npts_at_rowrank(my_row_rank)*ncache*nq
            IF (nelem .gt. 0) THEN
                inds(:) = 1
			    CALL Isend(probe_outputs,sirq,n_elements = nelem,dest = 0, &
                    tag=probe_tag, grp = pfi%rcomm, indstart = inds)            
                CALL IWait(sirq)
                If (allocated(probe_outputs)) Deallocate(probe_outputs)
            ENDIF

		Endif




        ! Communication is complete.  Now we open the file using MPI-IO
        

        ! For the moment, every process in column 0 participates in the mpi file-open operation
 
        If (my_row_rank .eq. 0) Call Point_Probes%OpenFile_Par(this_iter, error)

        If ( (responsible .eq. 1) .and. (Point_Probes%file_open) ) Then   
            funit = Point_Probes%file_unit
            If (Point_Probes%current_rec .eq. ncache) Then                
                
                If ( (Point_Probes%master) .and. (Point_Probes%write_header) ) Then    !           
                    ! The master rank (whoever owns the first radius output) writes the header
                    dims(1) = probe_nr
                    dims(2) = probe_nt
                    dims(3) = probe_np
                    dims(4) = nq

                    buffsize = 4
                    Call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 

                    buffsize = nq
                    Call MPI_FILE_WRITE(funit,Point_Probes%oqvals, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 

                    ! Radial grid -----------------------------------
                    nvals = max(probe_nr,probe_nt)
                    nvals = max(nvals,probe_np)
                    Allocate(probe_vals(1:nvals))
                    probe_vals(:) = 0
                    Do i = 1, probe_nr
                        probe_vals(i) = radius(Point_Probes%probe_r_global(i))                       
                    Enddo
                    buffsize = probe_nr

                    Call MPI_FILE_WRITE(funit, probe_vals, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 
                    Call MPI_FILE_WRITE(funit, Point_Probes%probe_r_global, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 


                    !Theta grid-----------------------------------------
                    DO i = 1, probe_nt
                        probe_vals(i) = costheta(Point_Probes%probe_t_global(i))                       
                    ENDDO
                    buffsize = probe_nt
                    CALL MPI_FILE_WRITE(funit, probe_vals, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 
                    CALL MPI_FILE_WRITE(funit, Point_Probes%probe_t_global, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 


                    !Phi grid----------------------------------

                    DO i = 1, probe_np
                        probe_vals(i) = (Point_Probes%probe_p_global(i)-1)*(two_pi/nphi)                        
                    ENDDO
                    buffsize = probe_np
                    CALL MPI_FILE_WRITE(funit, probe_vals, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 
                    CALL MPI_FILE_WRITE(funit, Point_Probes%probe_p_global, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 
                    DEALLOCATE(probe_vals)


                Endif
            Endif


            hdisp = 28 ! dimensions+endian+version+record count
            hdisp = hdisp+nq*4 ! nq
            hdisp = hdisp+probe_nr*12  ! radial indices and values
            hdisp = hdisp+probe_nt*12  ! theta  indices and values
            hdisp = hdisp+probe_np*12  ! phi indices and values

            qdisp = probe_nr*probe_nt*probe_np*8
            full_disp = qdisp*nq+12  ! 12 is for the simtime+iteration at the end
            disp = hdisp+full_disp*(Point_Probes%current_rec-ncache)

            buffsize = Point_Probes%npts_at_colrank(my_column_rank)
            ! The file is striped with time step slowest, followed by q


            pcount = 0
            Do p = 0, nproc1
                if (p .lt. my_column_rank) Then
                    pcount = pcount+ Point_Probes%npts_at_colrank(p)
                Endif
            Enddo

            my_pdisp = pcount*8

            Do j = 1, ncache


                Do i = 1, nq
                    new_disp = disp+qdisp*(i-1)+my_pdisp                
                    Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                    
                    Call MPI_FILE_WRITE(funit, row_probes(1,i,j), buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                Enddo

                tdisp = disp+full_disp-12

                Call MPI_File_Seek(funit,tdisp,MPI_SEEK_SET,ierr)


                If (Point_Probes%master) Then
                    buffsize2 = 1

                    Call MPI_FILE_WRITE(funit, Point_Probes%time_save(j), buffsize2, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                    Call MPI_FILE_WRITE(funit, Point_Probes%iter_save(j), buffsize2, & 
                           MPI_INTEGER, mstatus, ierr)
                Endif


                disp = disp+full_disp 
            Enddo

        Endif  ! Responsible

        If (responsible .eq. 1) DeAllocate(row_probes)
        If (my_row_rank .eq. 0) Call Point_Probes%closefile_par()

	End Subroutine Write_Point_Probes



	Subroutine Get_SPH_Modes(qty)
		Implicit None
		Integer :: j, ilocal, shell_ind, field_ind, rind, counter
        Integer :: k, jj
		Real*8, Intent(In) :: qty(1:,1:,my_theta_min:)

        If (SPH_Mode_Samples%nlevels .gt. 0) Then
            shell_ind = SPH_Mode_Samples%ind

            SPH_Mode_Samples%oqvals(shell_ind) = current_qval

        

		    If (SPH_Mode_Samples%my_nlevels .gt. 0) Then

		        If (SPH_Mode_Samples%begin_output) Then
			        Call sph_sample_buffer%construct('p3b')
                    sph_sample_buffer%p3b(:,:,:,:) = 0.0d0
		        Endif

                            

		        Do j = 1, SPH_Mode_Samples%my_nlevels



		            ilocal = SPH_Mode_Samples%my_shell_levs(j)-my_rmin+1

                    counter = (shell_ind-1)*SPH_Mode_Samples%my_nlevels+ j-1 

                    field_ind = counter/my_nr+1
                    rind = MOD(counter,my_nr)+my_rmin

                    Do k = 1, nphi
                    Do jj = my_theta_min, my_theta_max
				        sph_sample_buffer%p3b(k,rind,jj,field_ind) = &
                        & qty(k, ilocal, jj)
                    Enddo
                    Enddo

		        Enddo
            Endif
            Call SPH_Mode_Samples%AdvanceInd()
		Endif

	End Subroutine Get_SPH_Modes

	Subroutine Write_SPH_Modes(this_iter,simtime)
        ! This version mirrors the mem-friendly shell_spectra writing routine
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:,:), all_spectra(:,:,:,:,:)
        Real*8, Allocatable :: sendbuffer(:,:,:,:,:), out_radii(:)
        Real*8, Allocatable :: bsendbuffer(:,:,:,:,:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,qq, m, mp, lmax,rind,field_ind,f,r
        Integer :: rone,  p,  counter, nf
		Integer :: n, nn, this_nshell, nq_shell, sph_samples_tag, nmodes
        Integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp
        Integer(kind=MPI_OFFSET_KIND)  :: qsize, qdisp, rec_size
		Integer :: your_mp_min, your_mp_max, your_nm, your_id
		Integer :: nelem, m_ind, m_val, current_rec
        Integer :: funit, error, sirq, inds(5), dims(3)
        Integer :: my_nlevels, nlevels, qindex
        Integer :: lp1, nrirqs, ind5
        Integer :: ierr, rcount, buffsize, lv,lval
        Integer, Allocatable :: rirqs(:)
        Integer :: mstatus(MPI_STATUS_SIZE)       

        nlevels = SPH_Mode_Samples%nlevels             ! The total number of spectra levels that needs to be output
        my_nlevels = SPH_Mode_Samples%my_nlevels       ! The number of radial levels that this rank needs to write out
        nq_shell = SPH_Mode_Samples%nq                 ! The number of quantities 
        sph_samples_tag = SPH_Mode_Samples%mpi_tag
        funit = SPH_Mode_Samples%file_unit
        lmax = maxval(pfi%inds_3s)
        lp1 = lmax+1
        nmodes = sph_mode_nmode
		responsible = 0
		If ( (my_row_rank .eq. 0) .and. (my_nlevels .gt. 0) )  Then
            responsible = 1
            Allocate(all_spectra(0:lmax,0:lmax, my_nlevels,1, 1:2))
            Allocate(buff(0:lmax,my_nlevels,1,1:2,1:lp1))  !note - indexing starts at 1 not zero for mp_min etc.
            nrirqs = nproc2-1
            Allocate(rirqs(1:nrirqs))
        Endif


        !Before we start the main communication, all processes that contribute to the
        ! spectral output must get their buffers in the correct form
        If (my_nlevels .gt. 0) Then
            !//////////////////////
            ! First thing we do is FFT/reform the buffer/Legendre Transform
            !
            Call FFT_To_Spectral(sph_sample_buffer%p3b, rsc = .true.)
            sph_sample_buffer%config ='p3b'
            Call sph_sample_buffer%reform()
            Call sph_sample_buffer%construct('s2b')
            Call Legendre_Transform(sph_sample_buffer%p2b,sph_sample_buffer%s2b)
            Call sph_sample_buffer%deconstruct('p2b')

            Allocate(bsendbuffer(0:lmax,my_nlevels,nq_shell,2, my_mp_min:my_mp_max ))
            Allocate(sendbuffer(0:lmax,my_nlevels,1,2, my_mp_min:my_mp_max ))
            bsendbuffer = 0.0d0 
            sendbuffer = 0.0d0
            nf = sph_sample_buffer%nf2b
            Do p = 1, 2  ! Real and imaginary parts
            Do mp = my_mp_min,my_mp_max
                m = pfi%inds_3s(mp)
                    counter = 0
                    Do f = 1, nq_shell

                        field_ind = counter/my_nr+1
                        Do r = 1, SPH_Mode_Samples%my_nlevels   
                                
                            rind = MOD(counter,my_nr)+my_rmin
                            bsendbuffer(m:lmax,r,f,p,mp) = &
                                & sph_sample_buffer%s2b(mp)%data(m:lmax,rind,p,field_ind)
                            counter = counter+1
                        Enddo
                    Enddo
                Enddo

            Enddo
            call sph_sample_buffer%deconstruct('s2b')

        Endif


        If (my_row_rank .eq. 0) Call SPH_Mode_Samples%OpenFile_Par(this_iter, error)

        If ( (responsible .eq. 1) .and. (SPH_Mode_Samples%file_open) ) Then
            ! Processes that take part in the write have some extra work to do
            funit = SPH_Mode_Samples%file_unit
            current_rec = SPH_Mode_Samples%current_rec  ! Note that we have to do this after the file is opened
            If  ( (SPH_Mode_Samples%write_header) .and. (SPH_Mode_Samples%master) ) Then                

                dims(1) =  sph_mode_nell
                dims(2) =  nlevels
                dims(3) =  nq_shell
                buffsize = 3
                call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nq_shell
                call MPI_FILE_WRITE(funit,SPH_Mode_Samples%oqvals, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                allocate(out_radii(1:nlevels))
                Do i = 1, nlevels
                    out_radii(i) = radius(SPH_Mode_Samples%levels(i))
                Enddo
                buffsize = nlevels
	            call MPI_FILE_WRITE(funit, out_radii, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 
                DeAllocate(out_radii)
                
	            call MPI_FILE_WRITE(funit, SPH_Mode_Samples%levels, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = sph_mode_nell
	            call MPI_FILE_WRITE(funit, SPH_Mode_Ell, buffsize, MPI_INTEGER, & 
                    mstatus, ierr)                 

            Endif

            hdisp = 24 ! dimensions+endian+version+record count
            hdisp = hdisp+nq_shell*4 ! nq
            hdisp = hdisp+nlevels*12  ! level indices and level values
            hdisp = hdisp+sph_mode_nell*4  !ell-values
            rcount = 0
            Do p = 1, SPH_Mode_Samples%nshell_r_ids
                if (SPH_Mode_Samples%shell_r_ids(p) .lt. my_column_rank) Then
                    rcount = rcount+ SPH_Mode_Samples%nshells_at_rid(p)
                Endif
            Enddo
            my_rdisp = rcount*nmodes*8

                

            ! This is the LOCAL number ELEMENTS in the real or imaginary component of
            ! of a single quantity  (This is not in bytes)
            !buffsize = my_nlevels*nmodes 

            !This is the half-size (bytes) of a single quantity's information
            !Each quantity has real/imaginary components, and
            ! so the full size is twice this value.  THIS IS GLOBAL
            qsize = nlevels*nmodes*8

            !This is the size (bytes) of a single iteration's record
            rec_size = qsize*2*nq_shell+12  ! 12 is for the simtime+iteration at the end

            disp = hdisp+rec_size*(current_rec-1)

        Endif


        Do qindex = 1, nq_shell  ! Q LOOP starts here!

            !Load the current quantity into the sendbuffer
            If (my_nlevels .gt. 0) Then
                sendbuffer(:,:,1,:,:) = & 
                    & bsendbuffer(:,:,qindex,:,:)
            Endif



            If (responsible .eq. 1) Then
                ! Rank 0 in reach row receives  all pieces of the shell spectra from the other nodes

                all_spectra(:,:,:,:,:) = 0.0d0
                buff(:,:,:,:,:) = 0.0d0


                rirqs(:) = 0
                ind5 = pfi%all_3s(0)%delta+1
                Do nn = 1, nrirqs
                    !Write(6,*)'Ind5: ', ind5
                    your_id = nn

                    your_nm     = pfi%all_3s(nn)%delta
                    your_mp_min = pfi%all_3s(nn)%min
                    your_mp_max = pfi%all_3s(nn)%max


                    nelem = your_nm*my_nlevels*2*lp1

                    inds(:) = 1
                    inds(5) = ind5  !This is the mp_index here.

                    Call Ireceive(buff, rirqs(nn), n_elements = nelem,source= your_id, &
                        &  tag=sph_samples_tag,grp = pfi%rcomm, indstart = inds)
                    ind5 = ind5+your_nm
                Enddo

                ! Stripe my own data into the receive buffer

                Do mp = my_mp_min,  my_mp_max
                    m = pfi%inds_3s(mp)
                    Do p = 1,2
                        Do r = 1, my_nlevels   
                            buff(m:lmax,r,1,p,mp) = sendbuffer(m:lmax,r,1,p,mp) 
                        Enddo
                    Enddo
                Enddo

                Call IWaitAll(nrirqs,rirqs)

                !Stripe the receiver buffer into the spectra buffer

                !Modified stripe (we stripe m, ell  vs. ell, m as in shell_spectra)

                Do mp = 1,lp1
                    m = pfi%inds_3s(mp)
                    Do p = 1, 2  ! Real and imaginary parts
                        Do r = 1, my_nlevels   
                            all_spectra(m,m:lmax,r,1,p) = buff(m:lmax,r,1,p,mp)
                        Enddo
                    Enddo
                Enddo





                If (SPH_Mode_Samples%file_open) Then
                    Do p = 1, 2
                        new_disp = disp+  (qindex-1)*qsize*2 +(p-1)*qsize +my_rdisp  
                        Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                        Do r = 1, my_nlevels
                        Do lv = 1, SPH_MODE_NELL                    
                            lval = SPH_MODE_ELL(lv)
                            buffsize = lval+1
                            !if ((lval .eq. 0)) Then
                            !    Write(6,*)p, myid, my_rdisp, all_spectra(0,lval,r,1,p)
                            !Endif
                            Call MPI_FILE_WRITE(funit, all_spectra(0,lval,r,1,p), buffsize, & 
                                   MPI_DOUBLE_PRECISION, mstatus, ierr)
                        Enddo
                        Enddo
                    Enddo
                Endif

            Else
			    !  Non-responsible nodes send their info
			    If (my_nlevels .gt. 0) Then
                    inds(:) = 1
				    Call Isend(sendbuffer,sirq, dest = 0,tag=sph_samples_tag, grp = pfi%rcomm, indstart = inds)
                    Call IWait(sirq)
			    Endif
		    Endif

        Enddo  ! Q-LOOP

        If ( (responsible .eq. 1) ) Then
            disp = hdisp+rec_size*current_rec
            disp = disp-12

            If (SPH_Mode_Samples%file_open) Then
                Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)

                If (SPH_Mode_Samples%master) Then

                    buffsize = 1
                    Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                    Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                           MPI_INTEGER, mstatus, ierr)
                Endif

            Endif

            DeAllocate(all_spectra)
            DeAllocate(buff)
            DeAllocate(rirqs)
        Endif


        If (my_row_rank .eq. 0) Call SPH_Mode_Samples%closefile_par()
        If (my_nlevels .gt. 0) Then 
            DeAllocate(sendbuffer, bsendbuffer)
        Endif



	End Subroutine Write_SPH_Modes


    Subroutine Get_Equatorial_Slice(qty)
        Implicit None
        REAL*8, INTENT(IN) :: qty(:,:,my_theta_min:)
        INTEGER :: eq_ind 
        eq_ind = Equatorial_Slices%ind
        IF (my_nth_owned .gt. 0) THEN
            IF (Equatorial_Slices%begin_output) THEN
                ALLOCATE( equslice_outputs(1:nphi, my_rmin:my_rmax, 1:Equatorial_Slices%nq) )
                equslice_outputs(:,:,:) = 0.0d0
            ENDIF
            IF (th1_owner .eq. my_row_rank) THEN
                equslice_outputs(:,:,eq_ind) = 0.5d0*qty(:,:,ntheta/2)
            ENDIF
            IF (th2_owner .eq. my_row_rank) THEN
                equslice_outputs(:,:,eq_ind) = equslice_outputs(:,:,eq_ind)+0.5d0*qty(:,:,ntheta/2+1)
            ENDIF

            Equatorial_Slices%oqvals(eq_ind) = current_qval
        ENDIF
        Equatorial_Slices%oqvals(eq_ind) = current_qval
        Call Equatorial_Slices%AdvanceInd()
    End Subroutine Get_Equatorial_Slice

    Subroutine Write_Equatorial_Slices(this_iter, simtime)
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter

        integer :: sizecheck, responsible, current_rec, buffsize
        INTEGER :: nq_eqs, eqs_tag, funit, error, dims(3), i, ierr
		integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, qdisp, new_disp, my_rdisp, full_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
        Real*8, Allocatable :: buff(:,:,:), tbuff(:,:,:)
        INTEGER :: rirq1, rirq2, sirq, nelem

        sizecheck = sizeof(disp)
        IF (sizecheck .lt. 8) Then
            if (myid .eq. 0) Then
            Write(6,*)"Warning, MPI_OFFSET_KIND is less than 8 bytes on your system."
            Write(6,*)"Your size (in bytes) is: ", sizecheck
            Write(6,*)"A size of 4 bytes means that files are effectively limited to 2 GB in size."
            Endif
        ENDIF 
		responsible = 0
        nq_eqs   = Equatorial_Slices%nq
        eqs_tag  = Equatorial_Slices%mpi_tag
        funit    = Equatorial_Slices%file_unit

        If (my_row_rank .eq. 0) Then
            responsible = 1
        Endif

        !////////////////////////////////////////
        ! Communication down-row
        nelem = my_nr*nphi*nq_eqs
        If (responsible .eq. 1) Then
            ! We average over latitudes that bracket the equator.
            ! The 0.5 factor was already taken into account when
            ! the equatorial slices were sampled.
            ALLOCATE(buff(1:nphi,my_rmin:my_rmax,1:nq_eqs))
            
            IF (th1_owner .eq. th2_owner) THEN
                IF (th1_owner .eq. my_row_rank) THEN
                    buff(:,:,:) = equslice_outputs(:,:,:)
                ELSE
                    Call IReceive(buff, rirq1,n_elements = nelem, &
                        &  source= th1_owner,tag = eqs_tag, grp = pfi%rcomm)	
                    Call IWait(rirq1) 
                ENDIF

            ELSE
                Allocate(tbuff(1:nphi,my_rmin:my_rmax,1:nq_eqs))
                IF (th1_owner .eq. my_row_rank) THEN
                    buff(:,:,:) = equslice_outputs(:,:,:)
                ELSE
                    Call IReceive(buff, rirq1,n_elements = nelem, &
                        &  source= th1_owner,tag = eqs_tag, grp = pfi%rcomm)	
                    Call IWait(rirq1) 
                ENDIF
                IF (th2_owner .eq. my_row_rank) THEN
                    tbuff(:,:,:) = equslice_outputs(:,:,:)
                ELSE
                    Call IReceive(tbuff, rirq2,n_elements = nelem, &
                        &  source= th2_owner,tag = eqs_tag, grp = pfi%rcomm)	
                    Call IWait(rirq2) 
                ENDIF
                buff(:,:,:) = buff(:,:,:)+tbuff(:,:,:)
                DEALLOCATE(tbuff)
            ENDIF
            
            
            
            !This is where we would receive
        Else    
            IF ((th1_owner .eq. my_row_rank) .or. &
                & (th2_owner .eq. my_row_rank)) THEN
			    !  Rest of the row sends to process 0 within the row

                Call ISend(equslice_outputs, sirq,n_elements = nelem, dest = 0, tag = eqs_tag, & 
                    grp = pfi%rcomm)
                Call IWait(sirq)                
            
            ENDIF
        Endif



        !/////////////////////////////////////////
        ! Parallel Write amongst all processes in column 0
        If (responsible .eq. 1) Then

            Call Equatorial_Slices%OpenFile_Par(this_iter, error)

            current_rec = Equatorial_Slices%current_rec
            funit = Equatorial_Slices%file_unit

            If (Equatorial_Slices%file_open) Then

                If ( (my_column_rank .eq. 0) .and. (Equatorial_Slices%write_header) ) Then
                    dims(1) =  nphi
                    dims(2) =  nr
                    dims(3) =  nq_eqs
                    buffsize = 3
                    CALL MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 

                    buffsize = nq_eqs
                    CALL MPI_FILE_WRITE(funit,Equatorial_Slices%oqvals, buffsize, &
                        & MPI_INTEGER, mstatus, ierr) 

                    buffsize = nr
                    CALL MPI_FILE_WRITE(funit, radius, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 

                Endif

                hdisp = 24 ! dimensions+endian+version+record count
                hdisp = hdisp+nq_eqs*4 ! nq
                hdisp = hdisp+nr*8  ! The radius array

                qdisp = nphi*nr*8
                full_disp = qdisp*nq_eqs+12  ! 12 is for the simtime+iteration at the end
                disp = hdisp+full_disp*(current_rec-1)
                
                buffsize = my_nr*nphi
                ! The file is striped with time step slowest, followed by q

                my_rdisp = (my_rmin-1)*nphi*8
                Do i = 1, nq_eqs
                    new_disp = disp+qdisp*(i-1)+my_rdisp                
                    Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                    
                    Call MPI_FILE_WRITE(funit, buff(1,my_rmin,i), buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                Enddo

                disp = hdisp+full_disp*current_rec
                disp = disp-12
                Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)

                If (my_column_rank .eq. 0) Then
                    buffsize = 1
                    Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                    Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                           MPI_INTEGER, mstatus, ierr)
                Endif

            Endif

            Call Equatorial_Slices%closefile_par()
            Deallocate(buff)

        Endif
        

        !//////////////////////////////////////////////
        ! DeAllocation
        If (my_nth_owned .gt. 0) Then
            If (allocated(equslice_outputs)) DEALLOCATE( equslice_outputs)
        Endif


    End Subroutine Write_Equatorial_Slices


	Subroutine Get_Shell_Slice(qty)
		Implicit None
		Integer :: j, ilocal, shell_lev_ind, shell_ind
		Real*8, Intent(In) :: qty(:,:,my_theta_min:)
        If (Shell_Slices%nlevels .gt. 0) Then
            shell_ind = Shell_Slices%ind
            !If (myid .eq. 0) THen
                !NOTE:  Same as other remark when allocating.  Really should use master node here
                Shell_Slices%oqvals(shell_ind) = current_qval
            !Endif
        

		    If (Shell_Slices%my_nlevels .gt. 0) Then

		      If (Shell_Slices%begin_output) Then
			    Allocate(shell_slice_outputs(1:nphi,my_theta_min:my_theta_max,1:Shell_Slices%my_nlevels,1:Shell_Slices%nq))
		      Endif

            
		      shell_lev_ind =1
		      Do j = 1, Shell_Slices%nlevels
		        ilocal = Shell_Slices%levels(j)-my_rmin+1
		        If (Shell_Slices%have_shell(j) .eq. 1) Then ! my processor has this radius
				    shell_slice_outputs(:,my_theta_min:my_theta_max,shell_lev_ind,shell_ind) = &
                        &  qty(:,ilocal,my_theta_min:my_theta_max)
				    shell_lev_ind = shell_lev_ind +1
		        Endif
		      Enddo


		    Endif
            ! advance counter for next quantity to store (if any)
            Call Shell_Slices%AdvanceInd()
        Endif
	End Subroutine Get_Shell_Slice

	Subroutine Get_Shell_Spectra(qty)
		Implicit None
		Integer :: j, ilocal, shell_ind, field_ind, rind, counter
        Integer :: k, jj
		Real*8, Intent(In) :: qty(1:,1:,my_theta_min:)

        If (Shell_Spectra%nlevels .gt. 0) Then
            shell_ind = Shell_Spectra%ind
            !If (myid .eq. 0) Then
                Shell_Spectra%oqvals(shell_ind) = current_qval
            !Endif
        

		    If (Shell_Spectra%my_nlevels .gt. 0) Then

		        If (Shell_Spectra%begin_output) Then
			        Call spectra_buffer%construct('p3b')
                    spectra_buffer%p3b(:,:,:,:) = 0.0d0
		        Endif

                            

		        Do j = 1, Shell_Spectra%my_nlevels



		            ilocal = Shell_Spectra%my_shell_levs(j)-my_rmin+1

                    counter = (shell_ind-1)*Shell_Spectra%my_nlevels+ j-1 

                    field_ind = counter/my_nr+1
                    rind = MOD(counter,my_nr)+my_rmin

                    Do k = 1, nphi
                    Do jj = my_theta_min, my_theta_max
				        spectra_buffer%p3b(k,rind,jj,field_ind) = &
                        & qty(k, ilocal, jj)
                    Enddo
                    Enddo

		        Enddo
            Endif
            Call Shell_Spectra%AdvanceInd()
		Endif

	End Subroutine Get_Shell_Spectra

	Subroutine Write_Shell_Spectra(this_iter,simtime)
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:,:), all_spectra(:,:,:,:,:)
        Real*8, Allocatable :: sendbuffer(:,:,:,:,:), out_radii(:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,qq, m, mp, lmax,rind,field_ind,f,r
        Integer :: rone,  p,  counter, nf
		Integer :: n, nn, this_nshell, nq_shell, shell_spectra_tag, nmodes
        Integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp
        Integer(kind=MPI_OFFSET_KIND)  :: qsize, qdisp, rec_size
		Integer :: your_mp_min, your_mp_max, your_nm, your_id
		Integer :: nelem, m_ind, m_val, current_rec
        Integer :: funit, error, sirq, inds(5), dims(3)
        Integer :: my_nlevels, nlevels
        Integer :: lp1, nrirqs, ind5
        Integer :: ierr, rcount, buffsize
        Integer, Allocatable :: rirqs(:)
        Integer :: mstatus(MPI_STATUS_SIZE)       

        nlevels = Shell_Spectra%nlevels             ! The total number of spectra levels that needs to be output
        my_nlevels = Shell_Spectra%my_nlevels       ! The number of radial levels that this rank needs to write out
        nq_shell = Shell_Spectra%nq                 ! The number of quantities 
        shell_spectra_tag = Shell_Spectra%mpi_tag
        funit = Shell_Spectra%file_unit
        lmax = maxval(pfi%inds_3s)
        lp1 = lmax+1
        nmodes = lp1*lp1
		responsible = 0
		If ( (my_row_rank .eq. 0) .and. (my_nlevels .gt. 0) )  responsible = 1


        
        !/////////////
        If (my_nlevels .gt. 0) Then
            !//////////////////////
            ! First thing we do is FFT/reform the buffer/Legendre Transform
            !
            Call FFT_To_Spectral(spectra_buffer%p3b, rsc = .true.)
            spectra_buffer%config ='p3b'
            Call spectra_buffer%reform()
            Call spectra_buffer%construct('s2b')
            Call Legendre_Transform(spectra_buffer%p2b,spectra_buffer%s2b)
            Call spectra_buffer%deconstruct('p2b')

            Allocate(sendbuffer(0:lmax,my_nlevels,nq_shell,2, my_mp_min:my_mp_max ))
            sendbuffer = 0.0d0 


            
            nf = spectra_buffer%nf2b
            Do p = 1, 2  ! Real and imaginary parts
            Do mp = my_mp_min,my_mp_max
                m = pfi%inds_3s(mp)
                    counter = 0
                    Do f = 1, nq_shell

                        field_ind = counter/my_nr+1
                        Do r = 1, shell_spectra%my_nlevels   
                                
                            rind = MOD(counter,my_nr)+my_rmin
                            sendbuffer(m:lmax,r,f,p,mp) = &
                                & spectra_buffer%s2b(mp)%data(m:lmax,rind,p,field_ind)
                            counter = counter+1
                        Enddo
                    Enddo
                Enddo

            Enddo
            Call spectra_buffer%deconstruct('s2b')


        Endif

        If (responsible .eq. 1) Then
            ! Rank 0 in reach row receives  all pieces of the shell spectra from the other nodes

            Allocate(all_spectra(0:lmax,0:lmax, my_nlevels,nq_shell, 1:2))
            Allocate(buff(0:lmax,my_nlevels,nq_shell,1:2,1:lp1))  !note - indexing starts at 1 not zero for mp_min etc.
            all_spectra(:,:,:,:,:) = 0.0d0
            buff(:,:,:,:,:) = 0.0d0

            nrirqs = nproc2-1
            Allocate(rirqs(1:nrirqs))
            rirqs(:) = 0
            ind5 = pfi%all_3s(0)%delta+1
            Do nn = 1, nrirqs
                !Write(6,*)'Ind5: ', ind5
                your_id = nn

                your_nm     = pfi%all_3s(nn)%delta
                your_mp_min = pfi%all_3s(nn)%min
                your_mp_max = pfi%all_3s(nn)%max


                nelem = your_nm*my_nlevels*2*lp1*nq_shell

                inds(:) = 1
                inds(5) = ind5  !This is the mp_index here.

                Call Ireceive(buff, rirqs(nn), n_elements = nelem,source= your_id, &
                    &  tag=shell_spectra_tag,grp = pfi%rcomm, indstart = inds)
                ind5 = ind5+your_nm
            Enddo

            ! Stipe my own data into the receive buffer


            Do mp = my_mp_min,  my_mp_max
                m = pfi%inds_3s(mp)
                Do p = 1,2
                    Do f = 1, nq_shell
                        Do r = 1, my_nlevels   

                            buff(m:lmax,r,f,p,mp) = sendbuffer(m:lmax,r,f,p,mp) 

                        Enddo
                    Enddo
                Enddo

            Enddo
            !DeAllocate(sendbuffer)

            Call IWaitAll(nrirqs,rirqs)

            !Stripe the receiver buffer into the spectra buffer
           
            Do mp = 1,lp1
                m = pfi%inds_3s(mp)
                Do p = 1, 2  ! Real and imaginary parts
                    Do f = 1, nq_shell
                        Do r = 1, my_nlevels   
                            all_spectra(m:lmax,m,r,f,p) = buff(m:lmax,r,f,p,mp)  
                        Enddo
                    Enddo
                Enddo

            Enddo

            DeAllocate(sendbuffer)
            DeAllocate(buff)
            DeAllocate(rirqs)
        Else
			!  Non responsible nodes send their info
			If (my_nlevels .gt. 0) Then
                inds(:) = 1

				Call Isend(sendbuffer,sirq, dest = 0,tag=shell_spectra_tag, grp = pfi%rcomm, indstart = inds)
                Call IWait(sirq)
                DeAllocate(sendbuffer)
			Endif
		Endif


        If (my_row_rank .eq. 0) Call Shell_Spectra%OpenFile_Par(this_iter, error)



        If ( (responsible .eq. 1) .and. (Shell_Spectra%file_open) ) Then   
            !Write(6,*)'I am responsible: ', my_column_rank
            funit = shell_spectra%file_unit
            current_rec = Shell_Spectra%current_rec  ! Note that we have to do this after the file is opened
            If  ( (Shell_Spectra%write_header) .and. (Shell_Spectra%master) ) Then                

                dims(1) =  lmax
                dims(2) =  nlevels
                dims(3) =  nq_shell
                buffsize = 3
                call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nq_shell
                call MPI_FILE_WRITE(funit,Shell_Spectra%oqvals, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                allocate(out_radii(1:nlevels))
                Do i = 1, nlevels
                    out_radii(i) = radius(Shell_Spectra%levels(i))
                Enddo
                buffsize = nlevels
	            call MPI_FILE_WRITE(funit, out_radii, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 
                DeAllocate(out_radii)
                
	            call MPI_FILE_WRITE(funit, Shell_Spectra%levels, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

            Endif


            ! Depending on the offset size mpi type, disp will crap out past 2GB
            hdisp = 24 ! dimensions+endian+version+record count
            hdisp = hdisp+nq_shell*4 ! nq
            hdisp = hdisp+nlevels*12  ! level indices and level values
            


            ! The file is striped with time step slowest, followed by q


            rcount = 0
            Do p = 1, Shell_Spectra%nshell_r_ids
                if (Shell_Spectra%shell_r_ids(p) .lt. my_column_rank) Then
                    rcount = rcount+ Shell_Spectra%nshells_at_rid(p)
                Endif
            Enddo
            my_rdisp = rcount*nmodes*8

                

            ! This is the LOCAL number ELEMENTS in the real or imaginary component of
            ! of a single quantity  (This is not in bytes)
            buffsize = my_nlevels*nmodes 

            !This is the half-size (bytes) of a single quantity's information
            !Each quantity has real/imaginary components, and
            ! so the full size is twice this value.  THIS IS GLOBAL
            qsize = nlevels*nmodes*8

            !This is the size (bytes) of a single iteration's record
            rec_size = qsize*2*nq_shell+12  ! 12 is for the simtime+iteration at the end

            disp = hdisp+rec_size*(current_rec-1)


            Do p = 1, 2
                new_disp = disp+my_rdisp +(p-1)*qsize*nq_shell

                Do i = 1, nq_shell
                         
          
                    Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                    
                    Call MPI_FILE_WRITE(funit, all_spectra(0,0,1,i,p), buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)

                    new_disp = new_disp+qsize
                Enddo
            Enddo
            disp = hdisp+rec_size*current_rec
            disp = disp-12
            Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)


            If (shell_spectra%master) Then

                buffsize = 1
                Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
                Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                       MPI_INTEGER, mstatus, ierr)
            Endif



        Endif  ! Responsible & File Open

        If (responsible .eq. 1) DeAllocate(all_spectra)

        If (my_row_rank .eq. 0) Call Shell_Spectra%Closefile_Par()

	End Subroutine Write_Shell_Spectra


	Subroutine Write_Shell_Spectra_MEM(this_iter,simtime)
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:,:), all_spectra(:,:,:,:,:)
        Real*8, Allocatable :: sendbuffer(:,:,:,:,:), out_radii(:)
        Real*8, Allocatable :: bsendbuffer(:,:,:,:,:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,qq, m, mp, lmax,rind,field_ind,f,r
        Integer :: rone,  p,  counter, nf
		Integer :: n, nn, this_nshell, nq_shell, shell_spectra_tag, nmodes
        Integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp
        Integer(kind=MPI_OFFSET_KIND)  :: qsize, qdisp, rec_size
		Integer :: your_mp_min, your_mp_max, your_nm, your_id
		Integer :: nelem, m_ind, m_val, current_rec
        Integer :: funit, error, sirq, inds(5), dims(3)
        Integer :: my_nlevels, nlevels, qindex
        Integer :: lp1, nrirqs, ind5
        Integer :: ierr, rcount, buffsize
        Integer, Allocatable :: rirqs(:)
        Integer :: mstatus(MPI_STATUS_SIZE)       

        nlevels = Shell_Spectra%nlevels             ! The total number of spectra levels that needs to be output
        my_nlevels = Shell_Spectra%my_nlevels       ! The number of radial levels that this rank needs to write out
        nq_shell = Shell_Spectra%nq                 ! The number of quantities 
        shell_spectra_tag = Shell_Spectra%mpi_tag
        funit = Shell_Spectra%file_unit
        lmax = maxval(pfi%inds_3s)
        lp1 = lmax+1
        nmodes = lp1*lp1
		responsible = 0
		If ( (my_row_rank .eq. 0) .and. (my_nlevels .gt. 0) )  Then
            responsible = 1
            Allocate(all_spectra(0:lmax,0:lmax, my_nlevels,1, 1:2))
            Allocate(buff(0:lmax,my_nlevels,1,1:2,1:lp1))  !note - indexing starts at 1 not zero for mp_min etc.
            nrirqs = nproc2-1
            Allocate(rirqs(1:nrirqs))
        Endif


        !Before we start the main communication, all processes that contribute to the
        ! spectral output must get their buffers in the correct form
        If (my_nlevels .gt. 0) Then
            !//////////////////////
            ! First thing we do is FFT/reform the buffer/Legendre Transform
            !
            Call FFT_To_Spectral(spectra_buffer%p3b, rsc = .true.)
            spectra_buffer%config ='p3b'
            Call spectra_buffer%reform()
            Call spectra_buffer%construct('s2b')
            Call Legendre_Transform(spectra_buffer%p2b,spectra_buffer%s2b)
            Call spectra_buffer%deconstruct('p2b')

            Allocate(bsendbuffer(0:lmax,my_nlevels,nq_shell,2, my_mp_min:my_mp_max ))
            Allocate(sendbuffer(0:lmax,my_nlevels,1,2, my_mp_min:my_mp_max ))
            bsendbuffer = 0.0d0 
            sendbuffer = 0.0d0
            nf = spectra_buffer%nf2b
            Do p = 1, 2  ! Real and imaginary parts
            Do mp = my_mp_min,my_mp_max
                m = pfi%inds_3s(mp)
                    counter = 0
                    Do f = 1, nq_shell

                        field_ind = counter/my_nr+1
                        Do r = 1, shell_spectra%my_nlevels   
                                
                            rind = MOD(counter,my_nr)+my_rmin
                            bsendbuffer(m:lmax,r,f,p,mp) = &
                                & spectra_buffer%s2b(mp)%data(m:lmax,rind,p,field_ind)
                            counter = counter+1
                        Enddo
                    Enddo
                Enddo

            Enddo
            call spectra_buffer%deconstruct('s2b')

        Endif


        If (my_row_rank .eq. 0) Call Shell_Spectra%OpenFile_Par(this_iter, error)

        If ( (responsible .eq. 1) .and. (Shell_Spectra%file_open) ) Then
            ! Processes that take part in the write have some extra work to do
            funit = Shell_Spectra%file_unit
            current_rec = Shell_Spectra%current_rec  ! Note that we have to do this after the file is opened
            If  ( ( Shell_Spectra%write_header ) .and. ( shell_spectra%master ) ) Then                

                dims(1) =  lmax
                dims(2) =  nlevels
                dims(3) =  nq_shell
                buffsize = 3
                Call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nq_shell
                Call MPI_FILE_WRITE(funit,Shell_Spectra%oqvals, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                Allocate(out_radii(1:nlevels))
                Do i = 1, nlevels
                    out_radii(i) = radius(Shell_Spectra%levels(i))
                Enddo
                buffsize = nlevels
	            Call MPI_FILE_WRITE(funit, out_radii, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 
                DeAllocate(out_radii)
                
	            Call MPI_FILE_WRITE(funit, Shell_Spectra%levels, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

            Endif

            hdisp = 24 ! dimensions+endian+version+record count
            hdisp = hdisp+nq_shell*4 ! nq
            hdisp = hdisp+nlevels*12  ! level indices and level values
            
            rcount = 0
            Do p = 1, Shell_Spectra%nshell_r_ids

                If (Shell_Spectra%shell_r_ids(p) .lt. my_column_rank) Then

                    rcount = rcount + Shell_Spectra%nshells_at_rid(p)

                Endif

            Enddo
            my_rdisp = rcount*nmodes*8

                

            ! This is the LOCAL number ELEMENTS in the real or imaginary component of
            ! of a single quantity  (This is not in bytes)
            buffsize = my_nlevels*nmodes 

            !This is the half-size (bytes) of a single quantity's information
            !Each quantity has real/imaginary components, and
            ! so the full size is twice this value.  THIS IS GLOBAL
            qsize = nlevels*nmodes*8

            !This is the size (bytes) of a single iteration's record
            rec_size = qsize*2*nq_shell+12  ! 12 is for the simtime+iteration at the end

            disp = hdisp+rec_size*(current_rec-1)

        Endif


        Do qindex = 1, nq_shell  ! Q LOOP starts here!

            !Load the current quantity into the sendbuffer
            If (my_nlevels .gt. 0) Then
                sendbuffer(:,:,1,:,:) = & 
                    & bsendbuffer(:,:,qindex,:,:)
            Endif



            If (responsible .eq. 1) Then
                ! Rank 0 in reach row receives  all pieces of the shell spectra from the other nodes

                all_spectra(:,:,:,:,:) = 0.0d0
                buff(:,:,:,:,:) = 0.0d0


                rirqs(:) = 0
                ind5 = pfi%all_3s(0)%delta+1
                Do nn = 1, nrirqs
                    !Write(6,*)'Ind5: ', ind5
                    your_id = nn

                    your_nm     = pfi%all_3s(nn)%delta
                    your_mp_min = pfi%all_3s(nn)%min
                    your_mp_max = pfi%all_3s(nn)%max


                    nelem = your_nm*my_nlevels*2*lp1

                    inds(:) = 1
                    inds(5) = ind5  !This is the mp_index here.

                    Call Ireceive(buff, rirqs(nn), n_elements = nelem,source= your_id, &
                        &  tag=shell_spectra_tag,grp = pfi%rcomm, indstart = inds)
                    ind5 = ind5+your_nm
                Enddo

                ! Stripe my own data into the receive buffer

                Do mp = my_mp_min,  my_mp_max
                    m = pfi%inds_3s(mp)
                    Do p = 1,2
                        Do r = 1, my_nlevels   
                            buff(m:lmax,r,1,p,mp) = sendbuffer(m:lmax,r,1,p,mp) 
                        Enddo
                    Enddo
                Enddo

                Call IWaitAll(nrirqs,rirqs)

                !Stripe the receiver buffer into the spectra buffer
               
                Do mp = 1,lp1
                    m = pfi%inds_3s(mp)
                    Do p = 1, 2  ! Real and imaginary parts
                        Do r = 1, my_nlevels   
                            all_spectra(m:lmax,m,r,1,p) = buff(m:lmax,r,1,p,mp)  
                        Enddo
                    Enddo
                Enddo

                !Write the slice we just received
                If (Shell_Spectra%file_open) Then

                    Do p = 1, 2

                        new_disp = disp+my_rdisp +(p-1)*qsize*nq_shell +(qindex-1)*qsize        
                        Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                        
                        Call MPI_FILE_WRITE(funit, all_spectra(0,0,1,1,p), buffsize, & 
                               MPI_DOUBLE_PRECISION, mstatus, ierr)
                    Enddo

                Endif

            Else
			    !  Non-responsible nodes send their info
			    If (my_nlevels .gt. 0) Then
                    inds(:) = 1
				    Call Isend(sendbuffer,sirq, dest = 0,tag=shell_spectra_tag, grp = pfi%rcomm, indstart = inds)
                    Call IWait(sirq)
			    Endif
		    Endif

        Enddo  ! Q-LOOP

        If (responsible .eq. 1) Then
            disp = hdisp+rec_size*current_rec
            disp = disp-12

            If (Shell_Spectra%file_open) Then

                Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)

                If (Shell_Spectra%master) Then

                    buffsize = 1
                    Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                    Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                           MPI_INTEGER, mstatus, ierr)
                Endif

            Endif

            DeAllocate(all_spectra)
            DeAllocate(buff)
            DeAllocate(rirqs)

        Endif


        If (my_row_rank .eq. 0) Call Shell_Spectra%Closefile_Par()
        If (my_nlevels .gt. 0) Then 
            DeAllocate(sendbuffer, bsendbuffer)
        Endif



	End Subroutine Write_Shell_Spectra_MEM




	Subroutine Write_Shell_Slices(this_iter,simtime)
        USE RA_MPI_BASE
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:), all_shell_slices(:,:,:,:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,qq, p, sizecheck, t
		Integer :: n, nn, this_nshell, nq_shell, shell_slice_tag
		Integer :: your_theta_min, your_theta_max, your_ntheta, your_id
		Integer :: nelem, buffsize, sirq, nrirqs, inds(1:4)
        Integer :: file_pos, funit, error, dims(1:3), first_shell_rank
        Real*8, Allocatable :: out_radii(:)
        Integer, Allocatable :: level_inds(:), rirqs(:)
        
        integer :: ierr, rcount
		integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp, qdisp, full_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
        sizecheck = sizeof(disp)
        if (sizecheck .lt. 8) Then
            if (myid .eq. 0) Then
            Write(6,*)"Warning, MPI_OFFSET_KIND is less than 8 bytes on your system."
            Write(6,*)"Your size (in bytes) is: ", sizecheck
            Write(6,*)"A size of 4 bytes means that shell slices files are effectively limited to 2 GB in size."
            Endif
        endif 

        nq_shell = Shell_Slices%nq
        shell_slice_tag = Shell_Slices%mpi_tag
        funit = Shell_Slices%file_unit

		responsible = 0
		If (my_row_rank .eq. 0) Then
            If (Shell_Slices%my_nlevels .gt. 0) Then
                responsible = 1
            Endif
        Endif




        this_nshell = Shell_Slices%my_nlevels
		If (responsible .eq. 1) Then
			! Responsible node receives  all the pieces of the shell slices from the other nodes
			Allocate(all_shell_slices(1:nphi,1:ntheta,1:Shell_Slices%my_nlevels,nq_shell))

            Allocate(buff(1:nphi,1:this_nshell,1:nq_shell, 1:ntheta))            
			all_shell_slices(:,:,:,:) = 0.0d0
            buff(:,:,:,:) = 0.0d0

            ! Post Ireceives
            nrirqs = nproc2-1
            Allocate(rirqs(1:nrirqs))          

            Do nn = 1, nproc2-1
				your_id = nn

				your_ntheta    = pfi%all_2p(nn)%delta
				your_theta_min = pfi%all_2p(nn)%min

                inds(1) = 1
                inds(2) = 1
                inds(3) = 1
                inds(4) = your_theta_min

				nelem = nphi*your_ntheta*this_nshell*nq_shell

         		Call IReceive(buff, rirqs(nn),n_elements = nelem, source= your_id,tag=shell_slice_tag,grp = pfi%rcomm, indstart = inds)

            Enddo

            ! Stripe my own data into buff

            Do k = 1, nq_shell
                Do j = 1, this_nshell
                    Do t = my_theta_min, my_theta_max
                        Do i = 1, nphi
                            buff(i,j,k,t) = shell_slice_outputs(i,t,j,k)
                        Enddo
                    Enddo
                Enddo
            Enddo
            

            Call IWaitAll(nrirqs,rirqs)


            Do k = 1, nq_shell
                Do j = 1, this_nshell
                    Do t = 1, ntheta
                        Do i = 1, nphi
                            all_shell_slices(i,t,j,k) = buff(i,j,k,t)
                        Enddo
                    Enddo
                Enddo
            Enddo

            DeAllocate(buff)
            DeAllocate(rirqs)
            DeAllocate(shell_slice_outputs)
		Else
			!  Non responsible nodes send their info
			If (Shell_Slices%my_nlevels .gt. 0) Then
                !Everyone needs to restripe their data before sending it down the row
                !Stripe so that theta is slowest 
                Allocate(buff(1:nphi,1:this_nshell,1:nq_shell, my_theta_min:my_theta_max))
                Do t = my_theta_min, my_theta_max
                    Do k = 1, nq_shell
                        Do j = 1, this_nshell
                            Do i = 1, nphi
                                buff(i,j,k,t) = shell_slice_outputs(i,t,j,k)
                            Enddo
                        Enddo
                    Enddo
                Enddo
                nelem = nphi*my_ntheta*this_nshell*nq_shell
                inds(:) = 1
				Call Isend(buff,sirq,n_elements = nelem,dest = 0,tag=shell_slice_tag, grp = pfi%rcomm, indstart = inds)
            
                Call IWait(sirq)
                DeAllocate(shell_slice_outputs)
                DeAllocate(buff)
			Endif
		Endif


        ! Communication is complete.  Now we open the file using MPI-IO
        

        ! For the moment, every process in column 0 participates in the mpi operation
        ! The plan is to tune this later so that 
        If (my_row_rank .eq. 0) Call Shell_Slices%OpenFile_Par(this_iter, error)

        If ( (responsible .eq. 1) .and. (shell_slices%file_open) ) Then   

            funit = shell_slices%file_unit

            If (Shell_Slices%write_header) Then                
                
                If (Shell_Slices%master) Then            
                    ! The master rank (whoever owns the first output shell level) writes the header
                    dims(1) = ntheta
                    dims(2) = Shell_Slices%nlevels
                    dims(3) =  nq_shell
                    buffsize = 3
                    Call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 

                    buffsize = nq_shell
                    Call MPI_FILE_WRITE(funit,Shell_Slices%oqvals, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 

                    Allocate(out_radii(1:Shell_Slices%nlevels))
                    Do i = 1, Shell_Slices%nlevels
                        out_radii(i) = radius(Shell_Slices%levels(i))
                    Enddo
                    buffsize = Shell_Slices%nlevels
	                Call MPI_FILE_WRITE(funit, out_radii, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 
                    DeAllocate(out_radii)
                    

                    Allocate(level_inds(1:Shell_Slices%nlevels))
                    Do i = 1, Shell_Slices%nlevels
                        level_inds(i) = Shell_Slices%levels(i)
                    Enddo

	                Call MPI_FILE_WRITE(funit, Shell_Slices%levels, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 
                    DeAllocate(level_inds)
                    buffsize = ntheta
	                Call MPI_FILE_WRITE(funit, costheta, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 

                Endif
            Endif

            ! Maybe look into file views later.
            ! Depending on the offset size mpi type, disp will crap out past 2GB
            hdisp = 24 ! dimensions+endian+version+record count
            hdisp = hdisp+nq_shell*4 ! nq
            hdisp = hdisp+Shell_Slices%nlevels*12  ! level indices and level values
            hdisp = hdisp+ ntheta*8  ! costheta

            qdisp = ntheta*Shell_Slices%nlevels*nphi*8
            full_disp = qdisp*nq_shell+12  ! 12 is for the simtime+iteration at the end
            disp = hdisp+full_disp*(Shell_Slices%current_rec-1)
            
            buffsize = Shell_Slices%my_nlevels*ntheta*nphi
            ! The file is striped with time step slowest, followed by q

            rcount = 0
            Do p = 1, Shell_Slices%nshell_r_ids
                if (Shell_Slices%shell_r_ids(p) .lt. my_column_rank) Then
                    rcount = rcount+ Shell_Slices%nshells_at_rid(p)
                Endif
            Enddo
            my_rdisp = rcount*ntheta*nphi*8
            Do i = 1, nq_shell
                new_disp = disp+qdisp*(i-1)+my_rdisp                
                Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                
                Call MPI_FILE_WRITE(funit, all_shell_slices(1,1,1,i), buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
            Enddo
            disp = hdisp+full_disp*Shell_Slices%current_rec
            disp = disp-12
            Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)


            If (shell_slices%master) Then
                buffsize = 1
                Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)
                Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                       MPI_INTEGER, mstatus, ierr)
            Endif


			DeAllocate(all_shell_slices)
        Endif  ! Responsible

        If (my_row_rank .eq. 0) Call Shell_Slices%Closefile_Par()


	End Subroutine Write_Shell_Slices


	Subroutine Write_Shell_Slices_MEM(this_iter,simtime)
        ! A "more" memory friendly version of write_shell_Slices.  
        ! Writes one quantity at a time
        USE RA_MPI_BASE
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:,:), all_shell_slices(:,:,:,:)
		Integer :: responsible, current_shell, s_start, s_end, this_rid
		Integer :: i, j, k,qq, p, sizecheck, t
		Integer :: n, nn, this_nshell, nq_shell, shell_slice_tag
		Integer :: your_theta_min, your_theta_max, your_ntheta, your_id
		Integer :: nelem, buffsize, sirq, nrirqs, inds(1:4),qbuffsize
        Integer :: file_pos, funit, error, dims(1:3), first_shell_rank
        Real*8, Allocatable :: out_radii(:)
        Integer, Allocatable :: level_inds(:), rirqs(:)
        
        integer :: ierr, rcount, qindex
		integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp, qdisp, full_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
        sizecheck = sizeof(disp)
        If (sizecheck .lt. 8) Then
            if (myid .eq. 0) Then
            Write(6,*)"Warning, MPI_OFFSET_KIND is less than 8 bytes on your system."
            Write(6,*)"Your size (in bytes) is: ", sizecheck
            Write(6,*)"A size of 4 bytes means that shell slices files are effectively limited to 2 GB in size."
            Endif
        Endif 

        nq_shell = Shell_Slices%nq
        shell_slice_tag = Shell_Slices%mpi_tag
        

		responsible = 0
        this_nshell = Shell_Slices%my_nlevels
		If (my_row_rank .eq. 0) Then
            Call Shell_Slices%OpenFile_Par(this_iter, error)
            If (Shell_Slices%my_nlevels .gt. 0) Then
                responsible = 1
			    Allocate(all_shell_slices(1:nphi,1:ntheta,1:Shell_Slices%my_nlevels,1))
                Allocate(buff(1:nphi,1:this_nshell,1:1, 1:ntheta))  
                nrirqs = nproc2-1
                Allocate(rirqs(1:nrirqs))  

                ! Some displacements for accessing the file
                hdisp = 24 ! dimensions+endian+version+record count
                hdisp = hdisp+nq_shell*4 ! nq
                hdisp = hdisp+Shell_Slices%nlevels*12  ! level indices and level values
                hdisp = hdisp+ ntheta*8  ! costheta

                qdisp = ntheta*Shell_Slices%nlevels*nphi*8
                full_disp = qdisp*nq_shell+12  ! 12 is for the simtime+iteration at the end
                disp = hdisp+full_disp*(Shell_Slices%current_rec-1)
                ! The file is striped with time step slowest, followed by q
                rcount = 0
                Do p = 1, Shell_Slices%nshell_r_ids
                    If (Shell_Slices%shell_r_ids(p) .lt. my_column_rank) Then
                        rcount = rcount+ Shell_Slices%nshells_at_rid(p)
                    Endif
                Enddo
                my_rdisp = rcount*ntheta*nphi*8

                qbuffsize = Shell_Slices%my_nlevels*ntheta*nphi ! Number of elements in one q's worth of shells

            Endif
        Endif
        If (responsible .eq. 0) Then
            If (Shell_Slices%my_nlevels .gt. 0) Then
                Allocate(buff(1:nphi,1:this_nshell,1:1, my_theta_min:my_theta_max))
            Endif
        Endif
        funit = Shell_Slices%file_unit
        !////////////////////////////
        !Write a header
        If (Shell_Slices%file_open) Then                
            
            If (shell_slices%master .and. shell_slices%write_header) Then            
                ! The master rank (whoever owns the first output shell level) writes the header
                dims(1) = ntheta
                dims(2) = Shell_Slices%nlevels
                dims(3) =  nq_shell
                buffsize = 3
                Call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                buffsize = nq_shell
                Call MPI_FILE_WRITE(funit,Shell_Slices%oqvals, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 

                Allocate(out_radii(1:Shell_Slices%nlevels))
                Do i = 1, Shell_Slices%nlevels
                    out_radii(i) = radius(Shell_Slices%levels(i))
                Enddo
                buffsize = Shell_Slices%nlevels
	            Call MPI_FILE_WRITE(funit, out_radii, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 
                DeAllocate(out_radii)
                

                Allocate(level_inds(1:Shell_Slices%nlevels))
                Do i = 1, Shell_Slices%nlevels
                    level_inds(i) = Shell_Slices%levels(i)
                Enddo

	            Call MPI_FILE_WRITE(funit, Shell_Slices%levels, buffsize, MPI_INTEGER, & 
                    mstatus, ierr) 
                DeAllocate(level_inds)
                buffsize = ntheta
	            call MPI_FILE_WRITE(funit, costheta, buffsize, MPI_DOUBLE_PRECISION, & 
                    mstatus, ierr) 

            Endif
        Endif


        !///////////////////////////////////////////////////////////////////////////
        ! In this revised version, we send one quantity's worth of shells at a time
        Do qindex = 1, nq_shell

		    If (responsible .eq. 1) Then
			    ! Responsible node receives  all the pieces of the shell slices from the other nodes
          
			    all_shell_slices(:,:,:,:) = 0.0d0
                buff(:,:,:,:) = 0.0d0
            
                Do nn = 1, nproc2-1
				    your_id = nn

				    your_ntheta    = pfi%all_2p(nn)%delta
				    your_theta_min = pfi%all_2p(nn)%min

                    inds(1) = 1
                    inds(2) = 1
                    inds(3) = 1
                    inds(4) = your_theta_min

				    nelem = nphi*your_ntheta*this_nshell

             		Call IReceive(buff, rirqs(nn),n_elements = nelem, source= your_id, &
                        & tag=shell_slice_tag,grp = pfi%rcomm, indstart = inds)

                Enddo

                ! Stripe my own data into the receive buffer
                Do j = 1, this_nshell
                    Do t = my_theta_min, my_theta_max
                        Do i = 1, nphi
                            buff(i,j,1,t) = shell_slice_outputs(i,t,j,qindex)
                        Enddo
                    Enddo
                Enddo
                
                Call IWaitAll(nrirqs,rirqs)

                !Re-organize the buffer
                Do j = 1, this_nshell
                    Do t = 1, ntheta
                        Do i = 1, nphi
                            all_shell_slices(i,t,j,1) = buff(i,j,1,t)
                        Enddo
                    Enddo
                Enddo


		    Else
			    !  Non responsible nodes send their info
			    If (Shell_Slices%my_nlevels .gt. 0) Then
                    !Everyone needs to restripe their data before sending it down the row
                    !Stripe so that theta is slowest 

                    Do t = my_theta_min, my_theta_max
                        Do j = 1, this_nshell
                            Do i = 1, nphi
                                buff(i,j,1,t) = shell_slice_outputs(i,t,j,qindex)
                            Enddo
                        Enddo
                    Enddo
                    nelem = nphi*my_ntheta*this_nshell
                    inds(:) = 1
				    Call Isend(buff,sirq,n_elements = nelem,dest = 0,tag=shell_slice_tag, &
                        & grp = pfi%rcomm, indstart = inds)
                
                    Call IWait(sirq)

			    Endif
		    Endif


            ! Communication is complete.  Write this q-value using MPI-IO

            If ( (responsible .eq. 1) .and. (Shell_Slices%file_open) ) Then   

                new_disp = disp+qdisp*(qindex-1)+my_rdisp                
                Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                
                Call MPI_FILE_WRITE(funit, all_shell_slices(1,1,1,1), qbuffsize, & 
                       MPI_DOUBLE_PRECISION, mstatus, ierr)

            Endif  ! Responsible
        Enddo

        If (responsible .eq. 1) Then
			DeAllocate(all_shell_slices)
            DeAllocate(rirqs)
            disp = hdisp+full_disp*Shell_Slices%current_rec
            disp = disp-12

            If (Shell_Slices%file_open) Then

                Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)

                If (Shell_Slices%master) Then
                    buffsize = 1
                    Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                    Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                           MPI_INTEGER, mstatus, ierr)
                Endif

            Endif

        Endif
        If (my_row_rank .eq. 0) Call Shell_Slices%CloseFile_Par()
        If (Shell_Slices%my_nlevels .gt. 0) Then
            DeAllocate(shell_slice_outputs)
            DeAllocate(buff)
        Endif
	End Subroutine Write_Shell_Slices_MEM



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

            If (Equatorial_Slices%grab_this_q) Call Get_Equatorial_Slice(qty)
            If (Meridional_Slices%grab_this_q) Call Get_Meridional_Slice(qty)
            If (SPH_Mode_Samples%grab_this_q)  Call Get_SPH_Modes(qty)
            If (Point_Probes%grab_this_q)      Call Get_Point_Probes(qty)
		    If (Shell_Slices%grab_this_q)      Call Get_shell_slice(qty)
		    If (Shell_Spectra%grab_this_q)     Call Get_shell_spectra(qty)

            Call Get_Averages(qty)  ! -> Global Averages, Shell Averages, Azimuthal Averages

		    If (full_3d%grab_this_q) Call write_full_3d(qty)
        Endif


	End Subroutine	Add_Quantity

	Subroutine Complete_Output(iter, sim_time)
		Integer, Intent(In) :: iter
		Real*8, Intent(In) :: sim_time

	    If ((Shell_Slices%nq > 0) .and. (Mod(iter,Shell_Slices%frequency) .eq. 0 )) Then
            If (mem_friendly) Then
                Call Write_Shell_Slices_MEM(iter,sim_time)
            Else
                Call Write_Shell_Slices(iter,sim_time)
            Endif
        Endif
	    If ((Shell_Spectra%nq > 0) .and. (Mod(iter,Shell_Spectra%frequency) .eq. 0 )) Then
            If (mem_friendly) Then
                Call Write_Shell_Spectra_MEM(iter,sim_time)
            else
                Call Write_Shell_Spectra(iter,sim_time)
            Endif
        Endif
	    If ((Equatorial_Slices%nq > 0) .and. (Mod(iter,Equatorial_Slices%frequency) .eq. 0 )) Then

            Call Write_Equatorial_Slices(iter,sim_time)

        Endif
	    If ((Meridional_Slices%nq > 0) .and. (Mod(iter,Meridional_Slices%frequency) .eq. 0 )) Then
            Call Write_Meridional_Slices(iter,sim_time)
        Endif
	    If ((SPH_Mode_Samples%nq > 0) .and. (Mod(iter,SPH_Mode_Samples%frequency) .eq. 0 )) Then
            Call Write_SPH_Modes(iter,sim_time)
        Endif
	    If ((Point_Probes%nq > 0) .and. (Mod(iter,Point_Probes%frequency) .eq. 0 )) Then
            Point_Probes%time_save(Point_Probes%cc+1) = sim_time
            Point_Probes%iter_save(Point_Probes%cc+1) = iter
            If ((Point_Probes%cache_size-1) .eq. Point_Probes%cc) Then
                Call Write_Point_Probes(iter,sim_time)
            Endif            
            Call Point_Probes%AdvanceCC() 
        Endif
	    If ((AZ_Averages%nq > 0) .and. (Mod(iter,AZ_Averages%frequency) .eq. 0 )) Call Write_Azimuthal_Average(iter,sim_time)
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

		    f_of_r_theta(:,:) = 0.0D0
		
		    Do t = my_theta_min, my_theta_max
			    Do r = my_rmin, my_rmax
				    f_of_r_theta(r,t) = sum(qty(:,r,t))
			    Enddo
		    Enddo

		    f_of_r_theta = f_of_r_theta*over_nphi_double     ! average in phi

		    If (AZ_Averages%grab_this_q) Then
                If (.not. Allocated(azav_outputs)) Then
			        Allocate(azav_outputs(my_rmin:my_rmax,my_theta_min:my_theta_max,1:AZ_Averages%nq))
                Endif
    		    azav_ind = AZ_Averages%ind
			    azav_outputs(:,:,azav_ind) = f_of_r_theta
                If (myid .eq. 0) AZ_Averages%oqvals(azav_ind) = current_qval
			    Call AZ_Averages%AdvanceInd()
		    Endif
        Endif

        !/////////////////////
        !Now Average over the partial sphere (shell_averages)
		If (current_averaging_level .ge. 2) Then
		    f_of_r(:) = 0.0D0

            Do t = my_theta_min, my_theta_max
                f_of_r(:) = f_of_r(:) + f_of_r_theta(:,t)*theta_integration_weights(t)
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
            do i = my_rmin, my_rmax
                this_average = this_average+f_of_r(i)*r_integration_weights(i)
            enddo

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

	END Subroutine Get_Averages





	Subroutine Write_Azimuthal_Average(this_iter,simtime)
        USE RA_MPI_BASE
		Implicit None
		Real*8, Intent(in) :: simtime
		Integer, Intent(in) :: this_iter
		Real*8, Allocatable :: buff(:,:,:), all_azavgs(:,:,:)
		Integer :: responsible, current_rec, s_start, s_end, this_rid
		Integer :: i, j, k,qq, p, sizecheck
		Integer :: n, nn, this_nshell, nq_azav, az_avg_tag
		Integer :: your_theta_min, your_theta_max, your_ntheta
		Integer :: nelem, buffsize
        Integer :: file_pos, funit, error, dims(1:3)
        Integer :: inds(3), nirq,sirq
        Integer, Allocatable :: rirqs(:)
        
        integer :: ierr, rcount
		integer(kind=MPI_OFFSET_KIND) :: disp, hdisp, my_rdisp, new_disp, qdisp, full_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
        sizecheck = sizeof(disp)
        if (sizecheck .lt. 8) Then
            if (myid .eq. 0) Then
            Write(6,*)"Warning, MPI_OFFSET_KIND is less than 8 bytes on your system."
            Write(6,*)"Your size (in bytes) is: ", sizecheck
            Write(6,*)"A size of 4 bytes means that shell slices files are effectively limited to 2 GB in size."
            Endif
        endif 


		responsible = 0
        nq_azav     = AZ_Averages%nq
        az_avg_tag  = AZ_Averages%mpi_tag
        funit       = AZ_Averages%file_unit
        

		If (my_row_rank .eq. 0) Then
                responsible = 1
        Endif

        !Everyone needs to restripe their data for the sends so that theta is slowest
        Allocate(buff(my_rmin:my_rmax,1:nq_azav, my_theta_min:my_theta_max))
        Do k = my_theta_min, my_theta_max
            Do j = 1, nq_azav
                Do i = my_rmin, my_rmax
                    buff(i,j,k) = azav_outputs(i,k,j)  ! restripe
                Enddo
            Enddo
        Enddo
        DeAllocate(azav_outputs)


		If (responsible .eq. 1) Then
			! Rank 0 in reach row receives from all other row members
		
			Allocate(all_azavgs(my_rmin:my_rmax,1:nq_azav, 1:ntheta))
			all_azavgs(:,:,:) = 0.0d0

            nirq = nproc2-1
            Allocate(rirqs(1:nirq))

            Do nn = 1, nproc2-1

				your_ntheta    = pfi%all_2p(nn)%delta
				your_theta_min = pfi%all_2p(nn)%min
				your_theta_max = pfi%all_2p(nn)%max
                inds(1) = 1
                inds(2) = 1
                inds(3) = your_theta_min
                nelem = your_ntheta*my_nr*nq_azav

                Call IReceive(all_azavgs, rirqs(nn),n_elements = nelem, &
                            &  source= nn,tag = az_avg_tag, grp = pfi%rcomm,indstart = inds)				
			Enddo

            all_azavgs(my_rmin:my_rmax,1:nq_azav,my_theta_min:my_theta_max) = &
                & buff(my_rmin:my_rmax,1:nq_azav,my_theta_min:my_theta_max)

            Call IWaitAll(nirq, rirqs)
            DeAllocate(rirqs)
		Else
			!  Rest of the row sends to process 0 within the row
            inds(1) = 1 !my_rmin
            inds(2) = 1
            inds(3) = 1 ! my_theta_min
            nelem = my_nr*my_ntheta*nq_azav
            Call ISend(buff, sirq,n_elements = nelem, dest = 0, tag = az_avg_tag, & 
                grp = pfi%rcomm, indstart = inds)
            Call IWait(sirq)
		Endif
        DeAllocate(buff)

        ! Communication is complete.  Now we open the file using MPI-IO
        

      

        If (responsible .eq. 1) Then   
            Call AZ_Averages%OpenFile_Par(this_iter, error)


            current_rec = AZ_Averages%current_rec
            funit = AZ_Averages%file_unit
            !before we do anything else, we need to restripe the data yet again (might be able to work around this later)

            Allocate(buff( 1:ntheta,my_rmin:my_rmax, 1:nq_azav))
            Do k = 1, ntheta
                Do j = 1, nq_azav
                    Do i = my_rmin, my_rmax    
                        buff(k,i,j) = all_azavgs(i,j,k)
                    Enddo
                Enddo
            Enddo        
            DeAllocate(all_azavgs)
            
            If (AZ_Averages%file_open) Then
                If ((my_column_rank .eq. 0) .and. (AZ_Averages%write_header) ) Then            
                    ! Rank 0 in column and row writes the header
                    dims(1) =  nr
                    dims(2) =  ntheta
                    dims(3) =  nq_azav
                    buffsize = 3
                    call MPI_FILE_WRITE(funit, dims, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 

                    buffsize = nq_azav
                    call MPI_FILE_WRITE(funit,AZ_Averages%oqvals, buffsize, MPI_INTEGER, & 
                        mstatus, ierr) 

                    buffsize = nr
                    call MPI_FILE_WRITE(funit, radius, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 

                    buffsize = ntheta
                    call MPI_FILE_WRITE(funit, costheta, buffsize, MPI_DOUBLE_PRECISION, & 
                        mstatus, ierr) 
                Endif


                hdisp = 24 ! dimensions+endian+version+record count
                hdisp = hdisp+nq_azav*4 ! nq
                hdisp = hdisp+nr*8  ! The radius array
                hdisp = hdisp+ ntheta*8  ! costheta

                qdisp = ntheta*nr*8
                full_disp = qdisp*nq_azav+12  ! 12 is for the simtime+iteration at the end
                disp = hdisp+full_disp*(current_rec-1)
                
                buffsize = my_nr*ntheta
                ! The file is striped with time step slowest, followed by q

                my_rdisp = (my_rmin-1)*ntheta*8

                Do i = 1, nq_azav
                    new_disp = disp+qdisp*(i-1)+my_rdisp                
                    Call MPI_File_Seek(funit,new_disp,MPI_SEEK_SET,ierr)
                    
                    Call MPI_FILE_WRITE(funit, buff(1,my_rmin,i), buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                Enddo

                disp = hdisp+full_disp*current_rec
                disp = disp-12
                Call MPI_File_Seek(funit,disp,MPI_SEEK_SET,ierr)


                If (my_column_rank .eq. 0) Then
                    buffsize = 1
                    Call MPI_FILE_WRITE(funit, simtime, buffsize, & 
                           MPI_DOUBLE_PRECISION, mstatus, ierr)
                    Call MPI_FILE_WRITE(funit, this_iter, buffsize, & 
                           MPI_INTEGER, mstatus, ierr)
                Endif

            Endif
			DeAllocate(buff)
            Call AZ_Averages%closefile_par()
        Endif  ! Responsible

	End Subroutine Write_Azimuthal_Average


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
		Use RA_MPI_BASE ! Doing this here for now.  No other routine above sees MPI_Base, and I may want to keep it that way.
		Implicit None		
		Real*8, Intent(In) :: qty(:,my_rmin:,my_theta_min:)
		Real*8, Allocatable :: my_shells(:,:,:), buff(:,:,:)
		Integer :: i, j
		Character*4 :: qstring
		Character*8 :: iterstring
		Character*120 :: cfile

		Integer :: your_theta_min, your_theta_max, your_ntheta
		Integer :: np, buffsize, p
		Integer(kind=MPI_OFFSET_KIND) :: my_disp
		Integer :: mstatus(MPI_STATUS_SIZE)
		Integer :: funit, ierr, full_3d_tag

		! qty is dimensioned 1:n_phi, my_rmin:my_rmax, my_theta_min:my_theta_max

        full_3d_tag = Full_3D%mpi_tag
		If (my_row_rank .eq. 0) Then
			! Everyone in the row communicates to row-rank zero.
			! Each row owns a specific rank of radii
			Allocate(my_shells(1:nphi, 1:ntheta, my_rmin:my_rmax))
			my_shells(:,:,:) = 0.0d0					

			! First each rank stripes its own data into the new array
			! Not that striping is not in the "native" order
			Do j = my_theta_min, my_theta_max
				Do i = my_rmin, my_rmax
					my_shells(:,j,i) = qty(:,i,j)
				Enddo
			Enddo
			np = pfi%rcomm%np
			Do p = 1, np-1	
				your_theta_min = pfi%all_2p(p)%min
				your_theta_max = pfi%all_2p(p)%max
				your_ntheta    = pfi%all_2p(p)%delta
				Allocate(buff(1:nphi,my_rmin:my_rmax, your_theta_min:your_theta_max))
				Call receive(buff, source= p,tag=full_3d_tag,grp = pfi%rcomm)
				Do j = your_theta_min, your_theta_max
					Do i = my_rmin, my_rmax
						my_shells(:,j,i) = buff(:,i,j)
					Enddo
				Enddo
				DeAllocate(buff)
			Enddo
			! Now do the MPI write
			np = pfi%ccomm%np
			my_disp = 0
			Do p = 1, my_column_rank
				my_disp = my_disp+pfi%all_1p(p-1)%delta
			Enddo
			my_disp = my_disp*ntheta*nphi*8	! Displacment of this rank in the MPI File
			buffsize = my_nr*nphi*ntheta		! Number of elements to write


			write(iterstring,i_ofmt) current_iteration
			write(qstring,'(i4.4)') current_qval
         cfile = trim(local_file_path)//'Spherical_3D/'//trim(iterstring)//'_'//qstring
			!Write(6,*)cfile
			call MPI_FILE_OPEN(pfi%ccomm%comm, cfile, & 
                   MPI_MODE_WRONLY + MPI_MODE_CREATE, & 
                   MPI_INFO_NULL, funit, ierr) 

			call MPI_FILE_SET_VIEW(funit, my_disp, MPI_DOUBLE_PRECISION, & 
                   MPI_DOUBLE_PRECISION, 'native', & 
                   MPI_INFO_NULL, ierr) 
			call MPI_FILE_WRITE(funit, my_shells, buffsize, MPI_DOUBLE_PRECISION, & 
                   mstatus, ierr) 

			call MPI_FILE_CLOSE(funit, ierr) 
			If (my_column_rank .eq. 0) Then
				! row/column 0 writes out a file with the grid, etc.
				! This file should contain everything that needs to be known for processing later
	         write(iterstring,'(i8.8)') current_iteration
            cfile = trim(local_file_path)//'Spherical_3D/'//trim(iterstring)//'_'//'grid'
	         open(unit=15,file=cfile,form='unformatted', status='replace', access='stream')
           Write(15)endian_tag
	         Write(15)nr
				Write(15)ntheta
				Write(15)nphi
	         Write(15)(radius(i),i=1,nr)
				Write(15)(acos(costheta(i)),i = 1, ntheta)
	         Close(15)
			Endif


		Else
			! Send an array that's indexed starting at 1.  Shouldn't be necessary, but just in case.
			!Allocate(buff(1:nphi,1:my_ntheta,1:my_nr))
			
			!Do j = my_theta_min, my_theta_max
			!	jind = j-my_theta_min+1
			!	Do i = my_rmin, my_rmax
			!		buff(:,jind,i) = qty(:,i,j)
			!	Enddo
			!Enddo
			Call send(qty, dest= 0,tag=full_3d_tag,grp = pfi%rcomm)
			!DeAllocate(buff)

		Endif	


	End Subroutine Write_Full_3D


    !////////////////////////////////////////////////////////////////////////////////////////////////////////////
    !       Diagnostic Class Methods

    Subroutine Initialize_Diagnostic_Info(self,avg_levels,computes,pid,mpi_tag,avg_level,values, &
            levels, phi_inds, cache_size)
        Implicit None
        Integer :: i,ind
        Integer, Intent(In) :: pid, mpi_tag
        Integer, Intent(In), Optional :: cache_size
        Integer, Optional, Intent(In) :: avg_level
        Integer, Optional, Intent(In) :: values(1:)
        Integer, Optional, Intent(In) :: levels(1:)
        Integer, Optional, Intent(In) :: phi_inds(1:)
        Integer, Intent(InOut) :: computes(1:), avg_levels(1:)
        Class(DiagnosticInfo) :: self 
        If (present(avg_level)) Then
            if (avg_level .gt. 0) self%avg_level = avg_level
        Endif
        self%mpi_tag = mpi_tag
        self%nq = 0
        IF (present(cache_size)) THEN
            IF (cache_size .ge. 1) THEN
                self%cache_size = cache_size
            ENDIF
        ENDIF
        Allocate(self%iter_save(1:self%cache_size))
        Allocate(self%time_save(1:self%cache_size))


        If (present(values)) Then
            self%values(:) = values(:)  ! This is clunky - will look into getting the object attributes directly into a namelist later
            
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
        self%my_nlevels = 0
        If (present(levels)) Then
            self%levels(:) = levels(:)
            Do i = 1, nshellmax
                if( (self%levels(i) .gt. 0) .and. (self%levels(i) .le. nr) ) Then
                    self%nlevels = self%nlevels+1
                Endif
            Enddo
        Endif

        If (present(phi_inds)) Then
            self%nphi_indices = 0
            Do i=1,nmeridmax
                IF (meridional_indices(i) .gt. 0) THEN
                    IF (meridional_indices(i) .le. nphi) THEN
                        self%nphi_indices = self%nphi_indices+1
                    ENDIf
                ENDIF
            Enddo
            IF (self%nphi_indices .gt. 0) THEN
                ALLOCATE(self%phi_indices(1:self%nphi_indices))
                DO i = 1, self%nphi_indices
                    IF (meridional_indices(i) .le. nphi) THEN
                        self%phi_indices(i) = meridional_indices(i)
                    ENDIf
                ENDDO

            ENDIF
        Endif

        !if (pid .eq. 0) Then
            !NOTE:  Later, we may want to do this only on the master node later, but this isn't a huge memory issue
            Allocate(self%oqvals(1:self%nq))
            self%oqvals(:) = nqmax+100
        !Endif
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

    Subroutine Shell_Balance(self)
        ! I'm being a little sloppy here.  This method of the diagnostic class still uses
        ! a number of module-wide variables.
        ! Eventually, the code might be cleaner if this class were moved to its own module.
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer :: j, ilocal, i, pcount, your_rmax, your_rmin
        Integer, Allocatable :: ptemp(:), ptemp2(:), ptemp3(:)
		self%my_nlevels = 0
		Allocate(self%my_shell_levs(1:self%nlevels))
		Allocate(self%have_shell(1:self%nlevels))
		Allocate(self%my_shell_ind(1:self%nlevels))
		self%have_shell(:) = 0		
		Do j = 1, self%nlevels
			ilocal = self%levels(j)
			If ((ilocal .ge. my_rmin) .and. (ilocal .le. my_rmax)) Then ! my processor has this radius
			   self%have_shell(j) = 1
			   self%my_nlevels = self%my_nlevels+1
			   self%my_shell_levs(self%my_nlevels) = self%levels(j)
			   self%my_shell_ind(self%my_nlevels) = j
			Endif
		Enddo

        !/// ID 0 has a little more work to do
		If (my_row_rank .eq. 0) Then
			!Use process zero to figure out which radial processors will
			! actually output shells.  This will make it easier to read in the input later.
			Allocate(ptemp(1:nproc1))
			Allocate(ptemp2(1:nproc1))
			Allocate(ptemp3(1:nproc1))
			ptemp(:) = 0
			ptemp2(:) = 0
			do i = 0, nproc1-1
				your_rmin = pfi%all_1p(i)%min
				your_rmax = pfi%all_1p(i)%max
				Do j = 1, self%nlevels
					ilocal = self%levels(j)
					If ( (ilocal .ge. your_rmin) .and. (ilocal .le. your_rmax) ) Then
						ptemp(i+1) = ptemp(i+1)+1	!  radial id "i" has another shell that we want to output
					Endif	
				Enddo
			enddo
			pcount =0
			do i = 0, nproc1-1
				if (ptemp(i+1) .ge. 1) then 
					pcount = pcount+1					! pcount is the number of unique radial id's that have shells
					ptemp2(pcount) = i				! ptemp2 is for storing the radial id's of that have shells to output
					ptemp3(pcount) = ptemp(i+1)	! ptemp3 is the number of shells this radial id has
				endif
			enddo
            
			!   Resize the temporary arrays
			Allocate(self%shell_r_ids(1:pcount))
			Allocate(self%nshells_at_rid(1:pcount))
			self%shell_r_ids(:) = ptemp2(1:pcount)	! These are the radial ids of the processors that have shells
			self%nshells_at_rid(:) = ptemp3(1:pcount) ! How many shells this rid has
			self%nshell_r_ids = pcount
			DeAllocate(ptemp3)
			DeAllocate(ptemp2)
			DeAllocate(ptemp)
		Endif

    End Subroutine Shell_Balance

    Subroutine Init_OComm(self,pcomm,pnp,prank,mrank)
        Implicit None
        Integer, Intent(In) :: pcomm, pnp, prank, mrank
        Class(DiagnosticInfo) :: self
        self%ocomm = pcomm
        self%orank = prank
        self%onp = pnp
        If (prank .eq. mrank) then
            self%master = .true.    ! This process handles file headers in parallel IO
        Endif
    End Subroutine Init_OComm

    Subroutine set_file_info(self,oversion,rcount, freq, fpref,funit)
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer :: oversion, rcount, freq, modcheck
        Integer, Optional, Intent(In) :: funit
        Character*120 :: fpref
        If (present(funit)) Then
            self%file_unit = funit
        Endif
        self%file_prefix = fpref
        self%output_version = oversion
        self%rec_per_file = rcount
        self%frequency = freq

        !Now that we have set the file info, let's make sure the
        !that the cache size is appropriate

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
            modcheck = MOD(rcount,self%cache_size)
            IF (modcheck .ne. 0) THEN

                If (myid .eq. 0) THEN
                    Write(6,*)'////////////////////////////////////////////////////////////////////'
                    Write(6,*)'   Warning:  Incorrect cache_size specification for ',self%file_prefix
                    Write(6,*)'   Cache_size cannot be larger than nrec.'
                    Write(6,*)'   Cache_size: ', self%cache_size
                    Write(6,*)'   nrec      : ', self%rec_per_file
                    Write(6,*)'   Cache_size has been set to nrec.'
                    Write(6,*)'////////////////////////////////////////////////////////////////////'                    
                ENDIF
                self%cache_size = self%rec_per_file
            ENDIF
        Endif
        
        If (self%cache_size .gt. 1) Then
            modcheck = MOD(rcount,self%cache_size)
            IF (modcheck .ne. 0) THEN

                If (myid .eq. 0) THEN
                    Write(6,*)'////////////////////////////////////////////////////////////////////'
                    Write(6,*)'   Warning:  Incorrect cache_size specification for ',self%file_prefix
                    Write(6,*)'   Cache_size must divide evenly into nrec.'
                    Write(6,*)'   Cache_size: ', self%cache_size
                    Write(6,*)'   nrec      : ', self%rec_per_file
                    Write(6,*)'   Caching has been deactivated for ', self%file_prefix
                    Write(6,*)'////////////////////////////////////////////////////////////////////'                    
                ENDIF
                self%cache_size = 1
            ENDIF
        Endif

    End Subroutine set_file_info

    Subroutine OpenFile(self,iter,errcheck)
        Implicit None
        Class(DiagnosticInfo) :: self
        Integer, Intent(In) :: iter
        Integer, Intent(InOut) :: errcheck
        Character*8 :: iterstring,istr
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
        Character*8 :: iterstring
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
                disp = 8
                Call MPI_File_Seek(self%file_unit,disp,MPI_SEEK_SET,ierr)
                Call MPI_FILE_READ(self%file_unit, self%current_rec, 1, &
                    & MPI_INTEGER, mstatus, ierr)

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

        If (Allocated(self%my_shell_levs)) DeAllocate(self%my_shell_levs)
        If (Allocated(self%have_shell)) DeAllocate(self%have_shell)
        If (Allocated(self%my_shell_ind)) DeAllocate(self%my_shell_ind)
        If (Allocated(self%shell_r_ids)) DeAllocate(self%shell_r_ids)
        If (Allocated(self%nshells_at_rid)) DeAllocate(self%nshells_at_rid)
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
        histo_values(1:nqmax) = -1
        histo_levels(1:nshellmax)=-1

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
        If (Allocated(circumference)) DeAllocate(circumference)
        If (Allocated(qty)) DeAllocate(qty)
        If (Allocated(f_of_r_theta)) DeAllocate(f_of_r_theta)
        If (Allocated(azav_outputs)) DeAllocate(azav_outputs)
        If (Allocated(f_of_r)) DeAllocate(f_of_r)
        If (Allocated(rdtheta_total)) DeAllocate(rdtheta_total)
        If (Allocated(shellav_outputs)) DeAllocate(shellav_outputs)
        If (Allocated(globav_outputs)) DeAllocate(globav_outputs)
        If (Allocated(shell_slice_outputs)) DeAllocate(shell_slice_outputs)
    
        
        If (Allocated(sintheta_dtheta)) DeAllocate(sintheta_dtheta)
        If (Allocated(rsquared_dr)) DeAllocate(rsquared_dr)
        i_ofmt = '(i8.8)'  ! These should never change during a run, but just in case...
        i_pfmt = '(i5.5)'
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


    Subroutine Scattered_Balance(self, rinds, tinds, pinds)
        IMPLICIT NONE
        CLASS(DiagnosticInfo) :: self
        INTEGER, INTENT(INOUT) :: rinds(:), tinds(:), pinds(:)
        INTEGER, ALLOCATABLE :: tmp(:)
        INTEGER :: i, j,k, ind, this_min, this_max
        INTEGER :: rcount, tcount, pcount, ntmp
        INTEGER :: error_code = 0
        ! Few steps here:

        ! First, find out how many indices we need to output...

        rcount = size(rinds)
        i = 1
        DO WHILE ( i .le. rcount )
            IF (rinds(i) .lt. 0) THEN
                rcount = i-1
            ENDIF
            i = i+1
        ENDDO

        tcount = size(tinds)
        i = 1
        DO WHILE ( i .le. tcount )
            IF (tinds(i) .lt. 0) THEN
                tcount = i-1
            ENDIF
            i = i+1
        ENDDO

        pcount = size(pinds)
        i = 1
        DO WHILE ( i .le. pcount )
            IF (pinds(i) .lt. 0) THEN
                pcount = i-1
            ENDIF
            i = i+1
        ENDDO

        IF (rcount .gt. 0) THEN
            Allocate(self%probe_r_global(1:rcount))
            self%probe_r_global(1:rcount) = rinds(1:rcount)
        ENDIF

        IF (tcount .gt. 0) THEN
            Allocate(self%probe_t_global(1:tcount))
            self%probe_t_global(1:tcount) = tinds(1:tcount)
        ENDIF
        IF (pcount .gt. 0) THEN
            Allocate(self%probe_p_global(1:pcount))
            self%probe_p_global(1:pcount) = pinds(1:pcount)
        ENDIF

        self%probe_nr_global = rcount
        self%probe_nt_global = tcount
        self%probe_np_global = pcount
        !probe_nr_atrank(:)
		!my_theta_min = pfi%my_2p%min
		!my_theta_max = pfi%my_2p%max
		!my_rmin = pfi%my_1p%min
		!my_rmax = pfi%my_1p%max
        !nproc1 => processes per column
        !nproc2 => processes per row
        Allocate(self%probe_nr_atrank(0:nproc1-1))
        Allocate(self%probe_nt_atrank(0:nproc2-1))
        ALLOCATE(self%npts_at_colrank(0:nproc1-1))
        ALLOCATE(self%npts_at_rowrank(0:nproc2-1))

        self%npts_at_colrank(:) = 0
        self%npts_at_rowrank(:) = 0
        self%probe_nr_atrank(:) = 0
        self%probe_nr_local     = 0
        self%probe_nt_atrank(:) = 0
        self%probe_nt_local     = 0

        Do i = 0, nproc1-1 
            this_min = pfi%all_1p(i)%min
            this_max = pfi%all_1p(i)%max
            Do j = 1, rcount
                IF ( (rinds(j) .le. this_max) .and. (rinds(j) .ge. this_min) ) THEN
                    self%probe_nr_atrank(i) = self%probe_nr_atrank(i)+1
                ENDIF
            Enddo
        Enddo

        self%probe_nr_local = self%probe_nr_atrank(my_column_rank)
        Allocate(self%probe_r_local(1:self%probe_nr_local))
        ind = 1
        Do j = 1, rcount
            IF ( (rinds(j) .le. my_rmax) .and. (rinds(j) .ge. my_rmin) ) THEN
               self%probe_r_local(ind) = rinds(j) 
               ind = ind+1
            ENDIF
        Enddo       


        Do j = 1, tcount
            Do i = 0, nproc2-1 
                this_min = pfi%all_2p(i)%min
                this_max = pfi%all_2p(i)%max
                IF ( (tinds(j) .le. this_max) .and. (tinds(j) .ge. this_min) ) THEN
                    self%probe_nt_atrank(i) = self%probe_nt_atrank(i)+1
                ENDIF
            Enddo
        Enddo

        self%probe_nt_local = self%probe_nt_atrank(my_row_rank)
        Allocate(self%probe_t_local(1:self%probe_nt_local))
        ind = 1
        Do j = 1, tcount
            IF ( (tinds(j) .le. my_theta_max) .and. (tinds(j) .ge. my_theta_min) ) THEN
               self%probe_t_local(ind) = tinds(j) 
               ind = ind+1
            ENDIF
        Enddo 
        !WRite(6,*)'Check: ', my_theta_min, my_theta_max, ':', self%probe_t_local

        ! Calculate how many points are located on rank 0 of each row (just prior to parallel write)
        Do i = 0, nproc1-1
            self%npts_at_colrank(i) = self%probe_nr_atrank(i)*self%probe_nt_global*&
                self%probe_np_global
        Enddo


        ! Calculate how many points are located at each rank within a row
        Do i = 0, nproc2-1
            self%npts_at_rowrank(i) = self%probe_nr_atrank(my_column_rank)* &
                self%probe_nt_atrank(i)*self%probe_np_global
        Enddo



        !Write(6,*)'INITIALIZED!'

        ! The indexing seems to be working, but let's leave this code around for a bit
        ! just in case.
        !IF ( ( rcount .gt. 0) .and. (tcount .gt. 0) .and. (pcount .gt. 0) ) THEN
        !    IF (myid .eq. 0) THEN0
        !        WRITE(6,*)'Probes specified (phi,theta,r) :  '
        !        DO k = 1, rcount
        !            DO j = 1, tcount
        !                DO i = 1, pcount
        !                    Write(6,*)self%probe_p_global(i), self%probe_t_global(j), self%probe_r_global(k)
        !                ENDDO
        !            ENDDO
        !        ENDDO
        !    ENDIF
        !ELSE
        !    WRITE(6,*)'No probes specified.'
        !    error_code = 1
        !ENDIF
    End Subroutine Scattered_Balance

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
