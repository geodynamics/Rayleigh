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

Module ProblemSize
    Use Parallel_Framework, Only : pfi, Load_Config, Spherical
    Use Legendre_Polynomials, Only : Initialize_Legendre,coloc
    Use Spectral_Derivatives, Only : Initialize_Angular_Derivatives
    Use Controls, Only : Chebyshev, use_parity, multi_run_mode, run_cpus, my_path, outputs_per_row, m_balance_contiguous
    Use Chebyshev_Polynomials, Only : Cheby_Grid
    Use Finite_Difference, Only  : Initialize_Derivatives, Rescale_Grid_FD
    Use Math_Constants
    Use BufferedOutput
    Use Timers

    Implicit None

    !//////////////////////////////////////////////////////////////
    ! Processor Configuration
    Integer :: ncpu = 1, nprow = -1, npcol =-1 , npout = 1
    Integer :: my_rank      ! This is the rank within a run communicator
    Integer :: global_rank  ! This differs from my_rank only when multi-run mode is active
    Integer :: ncpu_global  ! Same as ncpu unless multi-run mode is active
    Integer :: my_row_rank, my_column_rank ! rank *within* row and rank *within* column

    !//////////////////////////////////////////////////////////////
    ! Horizontal Grid Variables
    Integer              :: n_theta = -1, n_phi
    Integer              :: l_max = -1, n_l=-1
    Integer              :: m_max, n_m
    Logical              :: dealias = .True.
    Integer, Allocatable :: m_values(:)
    Real*8, Allocatable  :: l_l_plus1(:), over_l_l_plus1(:)
    Real*8, Allocatable  :: costheta(:), sintheta(:), cos2theta(:), sin2theta(:), cottheta(:), csctheta(:)
    Real*8, Allocatable  :: phivals(:), cosphi(:), sinphi(:)
    Real*8               :: delta_phi
    Integer              :: k
    Type(Load_Config)    :: my_mp,  my_theta

    !//////////////////////////////////////////////////////////////
    !  Radial Grid Variables
    Integer             :: n_r, tnr
    Integer             :: grid_type = 1
    Real*8              :: rmin = -1.0d0, rmax = -1.0d0, r_inner, r_outer
    Real*8              :: aspect_ratio = -1.0d0
    Real*8              :: shell_depth = -1.0d0
    Real*8              :: shell_volume
    Real*8              :: stretch_factor = 0.0d0
    Real*8, Allocatable :: Radius(:), R_squared(:), One_Over_R(:)
    Real*8, Allocatable :: Two_Over_R(:), OneOverRSquared(:), Delta_R(:)
    Real*8, Allocatable :: radial_integral_weights(:)
    Integer :: precise_bounds = -1
    Type(Load_Config)   :: my_r

    !///////////////////////////////////////////////////////////////
    ! Radial Grid Variables Related to the Finite-Element Approach

    Integer, Parameter :: nsubmax = 256
    Integer :: ncheby(1:nsubmax)=-1
    Integer :: dealias_by(1:nsubmax) = -1
    Integer :: ndomains = 1
    Integer :: n_uniform_domains =1
    Real*8  :: domain_bounds(1:nsubmax+1)=-1.0d0
    Real*8  :: dr_input(4096) = 0.0d0
    Logical :: uniform_bounds = .false.
    Type(Cheby_Grid), Target :: gridcp

    !///////////////////////////////////////////////////////////////
    ! Error-handling variables
    Integer, Parameter :: maxerr = 32
    Integer :: perr(1:maxerr)
    Integer :: eindex = 1
    Logical :: grid_error = .false.


    Namelist /ProblemSize_Namelist/ n_r,n_theta, nprow, npcol,rmin,rmax,npout, dr_input, &
            &  precise_bounds,grid_type, l_max, n_l, &
            &  aspect_ratio, shell_depth, ncheby, domain_bounds, dealias_by, &
            &  n_uniform_domains, uniform_bounds, stretch_factor
Contains

    Subroutine Init_ProblemSize()
        Implicit None


        perr(:) = 0                       ! initialize the error-checking array

        Call Establish_Grid_Parameters()  ! Discern rmax,cheby-domain bounds, l_max, etc.
        Call Init_Comm()                  ! Initialize the MPI and the parallel framework

        If (my_rank .eq. 0) Then
            call stdout%print(" ")
            If (chebyshev) Then
                call stdout%print(" -- Initializing Chebyshev Grid...")
            Else
                call stdout%print(" -- Initializing Finite-Difference Grid...")
            EndIf
        Endif

        Call Report_Grid_Parameters()     ! Print some grid-related info to the screen
        Call Initialize_Horizontal_Grid() ! Init theta-grid and Legendre transforms
        Call Initialize_Radial_Grid()     ! Init radial grid and Chebyshev transforms

        If (my_rank .eq. 0) Then
            call stdout%print(" -- Grid initialized.")
            call stdout%print(" ")
        Endif

    End Subroutine Init_ProblemSize

    Subroutine Write_Grid()
        Implicit None
        Integer :: r
        Character*120 :: grid_file
        If (my_rank .eq. 0) Then
            grid_file = Trim(my_path)//'grid'
            Open(101, file = grid_file, status='replace', form = 'unformatted')
            Write(101) n_r
            Write(101) (radius(r),r=1,n_r)
            Write(101) n_theta
            Write(101) (costheta(r),r=1,n_theta)
            Close(101)
        Endif
    End Subroutine Write_Grid

    Subroutine Init_Comm()
        Implicit None
        Integer :: cpu_tmp(1)
        Integer :: ppars(1:12)
        !/////////////////////////////////////////////////////////
        ! Initialize MPI and load balancing (if no errors detected)
        ncpu = nprow*npcol
        ppars(1) = Spherical
        ppars(2) = n_r
        ppars(3) = n_r
        ppars(4) = n_theta
        ppars(5) = n_l
        ppars(6) = n_phi
        ppars(7) = n_m
        ppars(8) = ncpu
        ppars(9) = nprow
        ppars(10) = npcol
        ppars(11) = outputs_per_row
        If (m_balance_contiguous) Then
            ppars(12) = 1 ! use version 1
        Else
            ppars(12) = 0 ! use version 0, i.e., the original
        Endif
        If (multi_run_mode) Then
            Call pfi%init(ppars,run_cpus, grid_error)
        Else
            cpu_tmp(1) = ncpu
            Call pfi%init(ppars,cpu_tmp,grid_error)
        Endif
        my_rank = pfi%gcomm%rank
        my_row_rank = pfi%rcomm%rank
        my_column_rank = pfi%ccomm%rank


        my_mp    = pfi%my_3s
        my_r     = pfi%my_1p
        my_theta = pfi%my_2p

        tnr = 2*my_r%delta
        ! Once MPI has been initialized, we can start timing
        Call Initialize_Timers(my_path)
        Call StopWatch(init_time)%startclock()
        Call StopWatch(walltime)%startclock()
    End Subroutine Init_Comm

    Subroutine Auto_Assign_Domain_Decomp()
        Implicit None
        Integer :: f1, f2, imax, i, j, npairs, i3, i4, ii, jj
        Integer :: fcount
        Integer :: nprow_max, npcol_max, jinds(1:3), iinds(1:4)
        Integer, Allocatable, Dimension(:,:) :: factors_of_ncpu
        Integer, Allocatable :: factor_type(:), suitable_factors(:,:,:)
        Logical, Allocatable :: factor_balanced(:), have_pair(:,:)
        Real*8 :: ncpu_sqrt, rval, min_ratio, rdiff, rtol
        Real*8, Allocatable :: ratio_measure(:,:)
        Logical :: fewer_npcol, hbal, rbal, need_pair
        ncpu_sqrt = sqrt(dble(ncpu))
        npairs = n_l/2   ! Number of high-m, low-m pairs
        imax = int(ncpu_sqrt)+1
        npcol_max = n_r
        nprow_max = n_l/2 + Mod(n_l,2)

        Allocate(factors_of_ncpu(2,imax))
        Allocate(ratio_measure(1:4,1:3), have_pair(1:4,1:3))  
        Allocate(factor_type(1:imax*2))
        Allocate(factor_balanced(1:imax*2))
        Allocate(suitable_factors(1:2,1:4,1:3)) ! 1=> (f1,f2) -> (nprow, npcol), 2=> (f1,f2) -> (npcol,nprow)

        fewer_npcol = .true.   ! favor configurations with npcol < nprow when nprow == npcol is not possible

        factor_balanced(:) = .false.
        factor_type(:) = 0

        ratio_measure(:,:) = 2*ncpu  ! unachievably large, undesired ratio
        have_pair(:,:) = .false.
        suitable_factors(:,:,:) = 0

        ! Identify all factors of ncpu
        fcount = 0
        Do f1 = 1, imax
            If (Mod(ncpu,f1) .eq. 0) Then
                f2 = ncpu/f1
                If (Min(f1,f2) .ge. 2) Then  ! nprow and npcol must be at least 2
                    fcount = fcount+1
                    factors_of_ncpu(1:2,fcount) = (/ f1, f2 /)
                Endif
            Endif
        Enddo

        ! Test whether each (f1,f2) pair is suitable,
        ! but not necessarily optimal, as a choice for
        ! (nprow, npcol), (npcol, nprow), or both.
        ! The ideal combination is load balanced in radius and in the horizontal,
        ! with f1 = f2.  That combination may not be achieved, and so we store
        ! each combo in suitable_pairs(1:2, i, j)

        ! i = 1 => load balanced in both directions
        ! i = 2 => load balanced in the horizontal only
        ! i = 3 => load balanced in radius only
        ! i = 4 => not load balanced at all

        ! j = 1 => nprow == npcol
        ! j = 2 nprow > npcol
        ! j = 3 npcol < nprow 

        ! For each i,j combination, we store the pair with the smallest ratio
        ! of nprow/npcol or npcol/nprow (smaller of the two is used via MIN function)
        Do i = 1, fcount
            Do j = 1, 2
                ! Check both possible combinations (f1 = nprow, f2 = npcol) and (f2 = nprow, f1 = npcol)
                If (j .eq. 1) Then
                    f1 = factors_of_ncpu(1,i)
                    f2 = factors_of_ncpu(2,i)
                Else
                    f2 = factors_of_ncpu(1,i)
                    f1 = factors_of_ncpu(2,i)
                Endif
                If ((f1 .le. nprow_max) .and. (f2 .le. npcol_max) ) Then
                    ! Determine whether row > col or col > row or row == col
                    i4 = 1              ! assume equal
                    If (f1 > f2) i4 =2  ! row > col
                    If (f2 > f1) i4 = 3 ! col > row

                    ! Determine how these factors load balance
                    rbal = .false.
                    hbal = .false.
                    If ( Mod(N_R,f2) .eq. 0) rbal = .true.
                    If ( Mod(npairs,f1) .eq. 0) hbal = .true.
                    i3 = 4                              ! Not balanced in either direction.
                    If (rbal .and. (.not. hbal)) i3 = 3 ! Balanced in radial direction.
                    If (hbal .and. (.not. rbal)) i3 = 2 ! Balanced in theta direction.
                    If (rbal .and. hbal)         i3 = 1 ! Balanced in both directions.
                    
                    rval = ABS(1.0d0-MAX(DBLE(f1)/DBLE(f2),DBLE(f2)/DBLE(f1)))
    
                    If (rval .lt. ratio_measure(i3,i4) ) Then
                        ! This is a better combination for these parameters than we've encountered yet
                        suitable_factors(1:2,i3,i4) = (/ f1, f2 /)
                        ratio_measure(i3,i4) = rval
                        have_pair(i3,i4) = .true.
                    Endif                
                Endif
            Enddo
        Enddo


        ! Next, identify the smallest possible ratio
        min_ratio = MINVAL(ratio_measure)
        rtol = 2.0d0
        need_pair = .true.

        ! Now, parse the array of suitable pairs in such an order
        ! that our preferences in terms of nprow >,<,== npcol, and
        ! in terms of load balancing are respected
        iinds(1:4) = (/ 1, 2, 3, 4 /)  ! This really is probably the best order
        jinds(1:3) = (/ 1, 2, 3 /)
        If (.not. fewer_npcol) jinds = (/ 1, 3, 2 /)
        jj = 1
        Do While( need_pair .and. (jj .le. 3) )
            j = jinds(jj)
            ii = 1
            Do While (need_pair .and. (ii .le. 4) )
                i = iinds(ii)
                If (have_pair(i,j)) Then
                    rdiff = ratio_measure(i,j) - min_ratio
                    If (rdiff .le. rtol) Then
                        need_pair = .false.
                        nprow = suitable_factors(1,i,j)
                        npcol = suitable_factors(2,i,j)
                    Endif
                Endif

                ii = ii + 1
            Enddo
            jj = jj+1
        Enddo
        If (need_pair .and. (global_rank .eq. 0) ) Then
             Call stdout%print(' ')
             Call stdout%print('  /////////////////////////////////////////////////////////////////////////////// ')
             Call stdout%print('  ERROR:  Could not suitably factor specified process count for this problem size.')
             Call stdout%print('          Setting nprow and npcol to -1.')
             Call stdout%print('          Make sure not to run with a prime number of MPI ranks.')
             Call stdout%print(' ')
        Endif
    End Subroutine Auto_Assign_Domain_Decomp

    Subroutine Establish_Grid_Parameters()
        Implicit None
        Integer :: cheby_count, bounds_count, i
        Real*8 :: rdelta
        !Initialize everything related to grid resolution and domain bounds.

        If ((aspect_ratio .gt. 0.0d0) .and. (shell_depth .gt. 0.0d0) ) Then
            ! Set the bounds based on the aspect ratio and shell depth
            rmax = shell_depth/(1-aspect_ratio)
            rmin = rmax*aspect_ratio
        Endif

        !////////////////////////////////////////////////////////////////////
        ! Decide how many chebyshev domains we have
        !   - Users may specify:
        !       a.) a single Chebyshev domain (normal mode)
        !       b.) a uniform set of N Chebyshev domains (original SFE mode)
        !       c.) N Chebyshev domains with differing number of polynomials
        If ((ncheby(1) .le. 0 ) .and. (n_uniform_domains .lt. 2) ) Then
            ! Case (a)
            ncheby(1) = n_r
            domain_bounds(1) = rmin
            domain_bounds(2) = rmax
        Endif


        cheby_count = 0
        bounds_count = 1
        Do i = 1, nsubmax
            If (ncheby(i) .gt. 0) Then
                cheby_count = cheby_count+1
            Endif
            If (domain_bounds(i+1) .ge. 0) Then
                bounds_count = bounds_count+1
            Endif
        Enddo

        If (n_uniform_domains .gt. 1) Then
            ! Case (b)
            If (ncheby(1) .le. 0) Then
                ncheby(1) = n_r/n_uniform_domains
            Endif
            n_r =n_uniform_domains*ncheby(1)
            ncheby(1:n_uniform_domains) = ncheby(1)
            dealias_by(1:n_uniform_domains) = dealias_by(1)
            cheby_count = n_uniform_domains
            domain_bounds(1) = rmin
            rdelta = (rmax-rmin)/DBLE(cheby_count)
            bounds_count = cheby_count+1
            Do i = 2,bounds_count
                domain_bounds(i) = domain_bounds(i-1)+rdelta
            Enddo
        Else
            ! Case (c)
            n_r = 0
            Do i = 1, cheby_count
                n_r = n_r+ncheby(i)
            Enddo

            ! Here we check that a sufficient # of domain bounds were specified.
            ! If not, we check the value of uniform_bounds
            !   case i  : True => all domains are given the same width in radius
            !   case ii : False => code exits due to rmax < 0
            If (domain_bounds(cheby_count+1) .lt. 0) Then
                If (uniform_bounds) Then
                    n_uniform_domains = cheby_count
                    domain_bounds(1) = rmin
                    rdelta = (rmax-rmin)/DBLE(cheby_count)
                    bounds_count = cheby_count+1
                    Do i = 2,bounds_count
                        domain_bounds(i) = domain_bounds(i-1)+rdelta
                    Enddo
                Endif
            Endif
        Endif
        ndomains = cheby_count

        rmin = domain_bounds(1)
        rmax = domain_bounds(bounds_count)

        shell_volume = four_pi*one_third*(rmax**3-rmin**3)

        If ( (l_max .le. 0) .and. (n_l .le. 0) ) Then

            If (dealias) Then
                l_max = (2*n_theta/3) - 1
            Else
                l_max = n_theta-1
            Endif

        Else 
            If (l_max .le. 0) l_max = n_l-1
            !base n_theta on l_max
            If (dealias) Then
                n_theta = (3*(l_max+1))/2
            Else
                n_theta = l_max+1
            Endif
        Endif

        n_phi = 2*n_theta
        m_max = l_max
        n_l = l_max+1
        n_m = m_max+1
        If (.not. multi_run_mode) Then
            ! For single runs, Rayleigh will attempt to identify a
            ! sensible value of nprow and npcol if they were not
            ! specified by the user.
            If ((nprow .lt. 0) .and. (npcol .lt. 0) ) Then
                ncpu = ncpu_global
                Call Auto_Assign_Domain_Decomp()
            Endif
        Endif
        Call Consistency_Check()  ! Identify issues related to the grid setup

    End Subroutine Establish_Grid_Parameters

    Subroutine Initialize_Horizontal_Grid()
        Implicit None
        Integer :: tmp, l
        Integer, Allocatable :: m_vals(:)
        Real*8 :: ell
        !//////////////////////////////////////////////////
        ! Intialize Legendre Transforms & Horizontal Grid
        Allocate(m_values(1:n_m))
        m_values = pfi%inds_3s

        tmp = my_mp%delta
        allocate(m_vals(1:tmp))
        m_vals(:) = m_values(my_mp%min:my_mp%max)

        Call Initialize_Legendre(n_theta,l_max,m_vals,use_parity)
        tmp = my_r%delta
        Call Initialize_Angular_Derivatives(m_vals,l_max,tmp)
        DeAllocate(m_vals)


        Allocate(costheta(1:n_theta),cos2theta(1:n_theta))
        Allocate(sintheta(1:n_theta),sin2theta(1:n_theta))
        Allocate(csctheta(1:n_theta), cottheta(1:n_theta))
        Allocate(phivals(1:n_phi), cosphi(1:n_phi), sinphi(1:n_phi))
        costheta  = coloc    ! coloc computed in init_legendre
        cos2theta = costheta*costheta
        sin2theta = 1-cos2theta
        sintheta  = sqrt(sin2theta)
        csctheta = 1/sintheta
        cottheta = costheta/sintheta
        ! Calculate spacing of equally distributed phi points, then the phi grid
        ! The range of 0 to just below 2*pi (increasing) agrees with the Meridional
        ! and Equatorial Slices
        delta_phi = two_pi/n_phi
        Do k = 1, n_phi
            phivals(k) = (k-1)*delta_phi
            cosphi(k) = cos(phivals(k))
            sinphi(k) = sin(phivals(k))
        Enddo

        Allocate(l_l_plus1(0:l_max))
        Allocate(over_l_l_plus1(0:l_max))
        over_l_l_plus1(0) = 1.0d0 ! Prevent compiler warnings related to dividing by zero
        Do l = 0, l_max
            ell = l*1.0d0
            l_l_plus1(l) = ell*(ell+1)
            if (l .ne. 0) over_l_l_plus1(l) = 1.0d0/l_l_plus1(l)
        Enddo

    End Subroutine Initialize_Horizontal_Grid

    Subroutine Initialize_Radial_Grid()
        Implicit None
        Integer :: r, nthr,i ,n

        real*8 :: uniform_dr, arg, pi_over_N, rmn, rmx, delta, scaling
        real*8 :: delr0

        nthr = pfi%nthreads
        Allocate(Delta_r(1:N_R))
        Allocate( Radius(1:N_R))
        Allocate(Radial_Integral_Weights(1:N_R))



        If (chebyshev) Then ! No other choice at this time
            grid_type = 2

            Call gridcp%Init(radius, radial_integral_weights, &
                & ndomains,ncheby,domain_bounds,nthread = nthr, &
                & dealias_by = dealias_by)
            r = 1
            Do n = 1, ndomains
                Delta_r(r) = radius(r)-radius(r+1)
                r = r+1
                Do i = 2, gridcp%npoly(n)
                    Delta_r(r) = radius(r-1)-radius(r)
                    r = r+1
                Enddo
            Enddo
        Else
            grid_type = 1
            Radius(N_R) = rmin ! Follow ASH convention of reversed radius
            uniform_dr = 1.0d0/(N_R-1.0d0)*(rmax-rmin)
            If (dr_input(N_R) .eq. 0.0) Then
                Delta_r(N_R) = uniform_dr
            Else If (dr_input(N_R) .gt. 0.0) Then
                Delta_r(N_R) = dr_input(N_R)
            Endif

            Do r=N_R-1,1,-1
                    If (dr_input(r) .eq. 0.0) Then
                        Delta_r(r) = uniform_dr
                    Else If (dr_input(r) .gt. 0.0) Then
                        Delta_r(r) = dr_input(r)
                    Endif
                    Radius(r) = Delta_r(r) + Radius(r+1)
            Enddo
        Endif

        Allocate(OneOverRSquared(1:N_R),r_squared(1:N_R),One_Over_r(1:N_R),Two_Over_r(1:N_R))
        R_squared       = Radius**2
        One_Over_R      = (1.0d0)/Radius
        Two_Over_R      = (2.0d0)/Radius
        OneOverRSquared = (1.0d0)/r_Squared
        r_inner = rmin
        r_outer = rmax
        If (.not. chebyshev) Call Initialize_Derivatives(Radius,radial_integral_weights)
    End Subroutine Initialize_Radial_Grid

    Subroutine Report_Grid_Parameters()
        Implicit None
        Integer :: i
        Character*6 :: istr
        Character*12 :: dstring
        Character*8 :: dofmt = '(ES12.5)'

        Call Halt_On_Error()  ! Halt the program if errors were found in the consistency check
        If (my_rank .eq. 0) Then
            Write(istr,'(i6)')n_r
            Call stdout%print(" ---- Specified parameters ----")
            call stdout%print(" ---- N_R                 : "//trim(istr))
            Write(istr,'(i6)')n_theta
            call stdout%print(" ---- N_THETA             : "//trim(istr))
            Write(istr,'(i6)')l_max
            call stdout%print(" ---- Ell_MAX             : "//trim(istr))
            Write(dstring,dofmt)rmin
            call stdout%print(" ---- R_MIN               : "//trim(dstring))
            Write(dstring,dofmt)rmax
            Call stdout%print(" ---- R_MAX               : "//trim(dstring))
            Write(istr,'(i6)')ndomains
            If (chebyshev) Then
                call stdout%print(" ---- Chebyshev Domains   : "//trim(istr))
            Else
                call stdout%print(" ---- Finite Difference Domains   : "//trim(istr))
            Endif

            Do i = 1, ndomains
                call stdout%print(" ")
                Write(istr,'(i6)')i
                call stdout%print("        Domain "//&
                     & trim(adjustl(istr)))
                Write(istr,'(i6)')ncheby(i)
                call stdout%print("          Grid points           :  "//trim(adjustl(istr)))

                If (dealias_by(i) .eq. -1) Then
                    Write(istr,'(i6)')(ncheby(i)*2)/3
                Else
                    Write(istr,'(i6)')ncheby(i)-dealias_by(i)
                Endif
                call stdout%print("          Dealiased Polynomials :  "//trim(adjustl(istr)))
                Write(dstring,dofmt)domain_bounds(i)
                call stdout%print("          Domain Lower Bound    : "//trim(dstring))
                Write(dstring,dofmt)domain_bounds(i+1)
                call stdout%print("          Domain Upper Bound    : "//trim(dstring))
            Enddo
            call stdout%print(" ")
        Endif
    End Subroutine Report_Grid_Parameters



    Subroutine Add_Ecode(ecode)
        Integer, Intent(In) :: ecode
        perr(eindex) = ecode
        eindex = eindex+1
        grid_error = .true.
    End Subroutine Add_Ecode

    Subroutine Consistency_Check()
        Implicit None
        Integer :: tmp
        If (rmin .le. 0)    Call Add_Ecode(1)
        If (rmax .le. 0)    Call Add_Ecode(2)
        If (rmax .le. rmin) Call Add_Ecode(3)
        If (n_r .le. 0)     Call Add_Ecode(4)
        If (n_theta .le. 0) Call Add_Ecode(5)
        ! Check horizontal and radial load balancing
        tmp = (l_max+1)/2+MOD(l_max+1,2)  ! max nprow based on number of spherical harmonic pairs
        If (tmp .lt. nprow) call add_ecode(6)
        If (npcol .gt. n_r) call add_ecode(7)
        !Check parity
        If (MOD(n_theta,2) .ne. 0) Call Add_Ecode(8)
        tmp = 1
        Do While (tmp .le. ndomains)
            If (MOD(ncheby(tmp),2) .ne. 0) Then
                Call Add_Ecode(9)
                tmp = ndomains+1
            Endif
            tmp = tmp+1
        Enddo
        If ( (npcol .eq. 1) .or. (nprow .eq. 1) ) Call Add_Ecode(10)
    End Subroutine Consistency_Check

    Subroutine Halt_On_Error()

        Integer :: i,j,tmp, ecode
        Character*6 :: istr, istr2
        If (maxval(perr) .gt. 0) Then
            ecode = maxval(perr)
            If (my_rank .eq. 0) Then
                Call stdout%print(' /////////////////////////////////////////////////////////////////')
                Call stdout%print(' The following errors(s) were detected during grid initialization: ')
                Call stdout%print(' ')
                Do i = 1, maxerr
                    Select Case(perr(i))
                    Case(1)
                        Call stdout%print('  ERROR:  rmin must be greater than zero.')
                    Case(2)
                        Call stdout%print('  ERROR:  rmax must be greater than zero.')
                    Case(3)
                        Call stdout%print('  ERROR:  rmax must be greater than rmin.')
                    Case(4)
                        Call stdout%print('  ERROR:  n_r must be greater than zero.')
                        Write(istr,'(i6)')n_r
                        Call stdout%print('          current n_r   :'//trim(istr))

                    Case(5)
                        Call stdout%print('  ERROR:  n_theta must be greater than zero.')
                        Write(istr,'(i6)')n_theta
                        Call stdout%print('          specified n_theta   :'//trim(istr))
                    Case(6)
                        Call stdout%print('  ERROR:  NPROW is too large for this resolution.')
                        Call stdout%print('          NPROW may not be larger than:          ')
                        Call stdout%print('          (L_MAX+1)/2 + MOD(l_max+1,2)  ')
                        Write(istr,'(i6)')n_theta
                        Call stdout%print('          N_THETA           :'//trim(istr))
                        Write(istr,'(i6)')l_max
                        Call stdout%print('          L_MAX             :'//trim(istr))
                        tmp = (l_max+1)/2+MOD(l_max+1,2)


                        Write(istr,'(i6)')tmp
                        Call stdout%print('          Max allowed NPROW :'//trim(istr))
                        Write(istr,'(i6)')nprow
                        Call stdout%print('          Specified NPROW   :'//trim(istr))
                    Case(7)
                        Call stdout%print('  ERROR:  NPCOL is too large for this resolution.')
                        Call stdout%print('          NPCOL may not be larger than N_R.')
                        Write(istr,'(i6)')npcol
                        Call stdout%print('          current NPCOL :'//trim(istr))
                        Write(istr,'(i6)')n_r
                        Call stdout%print('          current N_R   :'//trim(istr))
                    Case(8)
                        Call stdout%print('  ERROR:  N_THETA must be Even.')
                        Call stdout%print('          (Legendre transforms exploit parity).')
                        Write(istr,'(i6)')n_theta
                        Call stdout%print('          current N_THETA :'//trim(istr))
                    Case(9)
                        Call stdout%print('  ERROR:  Each Chebyshev domains must contain an')
                        Call stdout%print('          even number of collocation points.')
                        Call stdout%print(" ")
                        Do j = 1, ndomains
                            Write(istr,'(i6)')j
                            Write(istr2,'(i6)')ncheby(j)
                            call stdout%print("          Grid points (domain "// &
                                &trim(adjustl(istr))// &
                                & ")          :  "//trim(adjustl(istr2)))
                        Enddo
                    Case(10)
                        Call stdout%print('  ERROR:  Neither nprow nor npcol may equal 1.')
                        Write(istr,'(i6)')npcol
                        Call stdout%print('          current NPCOL :'//trim(istr))
                        Write(istr,'(i6)')n_r
                        Call stdout%print('          current N_R   :'//trim(istr))
                    End Select
                    If (perr(i) .gt. 0) Call stdout%print(' ')
                Enddo
                Call stdout%print(' Exiting...')
                Call stdout%print(' /////////////////////////////////////////////////////////////////')

                Call stdout%finalize()
            Endif
            Call pfi%exit(ecode)
        Endif
    End Subroutine Halt_On_Error
End Module ProblemSize
