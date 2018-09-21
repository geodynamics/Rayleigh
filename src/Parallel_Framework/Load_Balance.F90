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

Module Load_Balance
    USE MPI_LAYER

    Implicit None

    Type, Public :: Load_Config
        Integer :: min, max !min and max index in global array.
        Integer :: delta ! The number of LOCAL indices stored on this rank.  delta = max-min+1
        Integer :: np  ! The number of processors the load is distributed across
        Integer :: nx  ! The number of GLOBAL indices distributed across np MPI ranks
    End Type Load_Config

    !///////////////////////////////////////////////////////////////////////////
    ! Spherical Geometry Load Balancing
    Type(Load_Config) :: my_mp, my_theta, my_r

    !///////////////////////////////////////////////////////////////////////////
    ! Cartesian Geometry Load Balancing
    Type(Load_Config), Allocatable, Public :: lb_z(:), lb_y(:), lb_x(:)


    !///////////////////////////////////////////////////////////////////////////////////////////////////////
    !            This is Absolutely Ugly and needs to be re-organized
    !             **FOR NOW** make these variables part of this module - rather than attributes of an object.
    !===========================================
    ! L-M load balancing / Parallelization
    !  --this information needs to get out to the rest of the code ultimately
    Integer, Allocatable :: l_lm_values(:)    ! The l-values of each l-m combination in this radial group
    Integer, Allocatable :: m_lm_values(:)    ! The m-values of each l-m combination in this radial group
    Integer, Allocatable :: mp_lm_values(:)   ! The mp index in the IDX2 configuration corresponding to each l-m combo
    Integer, Allocatable :: lm_owner(:)       ! The radial group rank of the processor responsible for solving for each l-m RHS
    Integer, Allocatable :: num_lm(:)         ! The number of l-m combinations assigned to each process in the radial group
    Integer, Allocatable :: my_lm_lval(:)     ! The l-values of the matrices I am responsible for
    Integer, Allocatable :: my_nm_lm(:)       ! The number of m-values associated with each l that I own
    Integer :: my_num_lm                      ! The number of l-lm combinations I am responsible for
    Integer :: my_nl_lm                       ! The number of l-values I am responsible for
    Integer :: my_lm_min, my_lm_max           ! The min and max of the indices I own within the l/m_lm_values arrays
    Integer :: lm_count                       ! The number of l-m combinations in this radial group
    !////////////////////////////////////////////////////////////////////////////////////////////////////////

Contains

    !////////////////////////////////////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  Standard_Balance
    !
    ! DESCRIPTION:  Handles load balancing of grid indices across the members of
    !                 of the MPI communicator group comm_group.  Grid indices are
    !                 assumed to run from 1 through nx.  Balancing is accomplished
    !                 in the standard, contiguous fashion.
    !
    ! INPUTS:
    !            nx         - the number of indices (i.e. gridpoints) to load balance across comm_group
    !            comm_group - MPI communicator across which load balancing is to be done
    !
    ! OUTPUTS:
    !
    !            lb_in      - array of load_config objects containing the load balancing
    !                           information for each member of comm_group
    !//////////////////////////////////////////////////////////////////////////////////////////////////////
    Subroutine Standard_Balance(lb_in,nx,comm_group)
        Implicit None
        Type(Load_Config), Intent(InOut) :: lb_in(0:)
        Type(communicator), Intent(In) :: comm_group
        Integer, Intent(In) :: nx        ! the size of the dimension we want to load balance (e.g. nr_global)

        Integer :: p,np, mcheck, nx_local

        np = comm_group%np ! Number of ranks in comm_group
        nx_local = nx/np   ! All ranks get at least nx_local points

        mcheck = mod(nx,np)-1 ! Processors with rank less than mcheck get nx_local+1 points
        lb_in(0)%min = 1      ! First grid index is assumed to be 1 and is assigned to process 0

        Do p = 0, np-1
            If (p .gt. 0) lb_in(p)%min = lb_in(p-1)%max+1
            lb_in(p)%delta = nx_local
            lb_in(p)%nx = nx
            If (p .le. mcheck) then
                !This rank gets 1 additional point
                lb_in(p)%delta = nx_local+1
            Endif
            lb_in(p)%max = lb_in(p)%min+lb_in(p)%delta-1
            lb_in(p)%nx = nx
        Enddo
    End Subroutine Standard_Balance

    !////////////////////////////////////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  M_Balance
    !
    ! DESCRIPTION:  Handles load balancing of grid indices across the members of
    !                 of the MPI communicator group 'comm'.  Grid indices are
    !                 assumed to run from 1 through nx.  Balancing is accomplished
    !                 in the "m"-fashion such that pairs of high- and low-valued
    !                 grid indices are distributed amongst the processes in comm_group.
    !
    ! INPUTS:
    !            comm   - MPI communicator across which load balancing is to be done
    !
    ! OUTPUTS:
    !
    !            m_values - integer array which, upon exit of the subroutine, contains a list of
    !                       grid indices ordered as they have been distributed across processes.
    !            lb_in    - array of load_config objects containing the load balancing
    !                           information for each member of comm_group
    !//////////////////////////////////////////////////////////////////////////////////////////////////////
    Subroutine M_Balance(lb_in, m_values,comm)
        ! Pair up high and low m's
        ! Does not have to be m's.  Can be any index
        Implicit None
        Type(Load_Config), Intent(InOut) :: lb_in(0:)
        Type(Communicator), Intent(In) :: comm
        Logical :: found
        Integer, Intent(InOut) :: m_values(1:)
        Integer :: n_m, ind, dpair,k
        Integer :: npairs , np, i, mcheck, mextra
        Integer :: offset, current_pair,p,my_npairs
        Integer, Allocatable :: pairs(:,:), unpaired(:), i_am_holding(:)
        n_m = lb_in(0)%nx

        npairs = n_m/2        ! Total number of high/low m-pairs
        mcheck = Mod(n_m,2)   ! possibly have a lone-m value (a sign of poor user-specified grid dimensions)
        Allocate(pairs(1:2,1:npairs))
        np = comm%np

        !Generate a list of high/low pairs
        Do i = 0, npairs-1
            pairs(1,i+1) = i
            pairs(2,i+1) = n_m-i-1
        Enddo

        Allocate(i_am_holding(0:np-1))  ! Keep track of how many m-values each process has
        i_am_holding(:) = 0

        !  Distribute complete pairs amongst the different cpus
        !  Do not break up any pairs at this stage
        offset = 0
        current_pair = 1
        do p = 0, np-1
            my_npairs = lb_in(p)%delta/2
            ind = lb_in(p)%min
            do i = 1, my_npairs
                m_values(ind) = pairs(1,current_pair)
                m_values(ind+1) = pairs(2,current_pair)
                current_pair = current_pair+1
                ind = ind+2
                i_am_holding(p) = i_am_holding(p)+2
            enddo
        enddo




        !////////////
        ! find the first process that has a 'hole' to stuff
        ! a broken-up pair into.
        p = 0
        found = .false.

        do while(.not. found)
            if (i_am_holding(p) .eq. lb_in(p)%delta) then
                p = p+1
            else
                found = .true.
            endif
            if (p .eq. np) then
                found = .true.    ! Everyone was full on m-values
            endif
        enddo

        !///////////
        !Break up any pairs that are left
        If (current_pair .ne. (npairs+1)) then

            dpair = npairs+1-current_pair

            Allocate(unpaired(1:dpair*2))
            ind = 1
            do i = 1, dpair
                unpaired(ind) = pairs(1,current_pair)
                unpaired(ind+1) = pairs(2,current_pair)
                current_pair = current_pair+1
                ind = ind+2
            enddo

            do i = 1, dpair
                do k = 1, 2
                    my_npairs = lb_in(p)%delta/2
                    ind = lb_in(p)%min+my_npairs*2
                    m_values(ind) = unpaired((i-1)*2+k)
                    p = p+1
                enddo
            enddo
            DeAllocate(unpaired)
        Endif

        ! Account for middle m_value if it was left unpaired
        if (mcheck .eq. 1) then
            mextra = n_m/2
            my_npairs = lb_in(p)%delta/2
            ind = lb_in(p)%min+my_npairs*2
            m_values(ind) = mextra
        endif

        DeAllocate(i_am_holding)
        DeAllocate(pairs)
    End Subroutine M_Balance


    !////////////////////////////////////////////////////////////////////////////////////////////////////
    ! SUBROUTINE:  LM_Load_Balance
    !
    ! DESCRIPTION:  Handles load balancing of spherical harmonic modes (l-m pairs) across the members of
    !                 of the MPI communicator group 'comm'.  This load balancing is appropriate for
    !                 a triangular truncation.
    !
    !
    ! INPUTS:
    !            comm   - MPI communicator across which load balancing is to be done
    !            mlb    - load_config object returned from M_Balance
    !            inds   - array of m_values as returned from M_Balance
    ! OUTPUTS:
    !           NONE
    !
    ! Side Effects:
    !                All variables declared in the L-M load balancing / Parallelization section at the
    !                top of this module are initialized.  These will be wrapped into their own load balancing
    !                object at some point in the future.
    ! Notes:
    !                Editing this routine is not a good idea.
    !//////////////////////////////////////////////////////////////////////////////////////////////////////

    Subroutine LM_Load_Balance(mlb,inds,comm)
        ! Handles the load balancing of modes for the implicit transpose
        Implicit None
        Type(Load_Config), Intent(In) :: mlb
        Integer :: nmv, nlv, l_min, lm_per, lm_remainder, l, m
        Integer :: l_lm_min, l_lm_max, mp, m_count, i, ind
        Integer, Allocatable :: mv_temp(:), lv_temp(:)
        Integer ::  mpmin,mpmax

        Integer :: l_max,n_radial_cpus,my_radial_rank
        Integer, Intent(In) :: inds(:)
        Type(Communicator), Intent(In) :: comm

        my_radial_rank = comm%rank

        l_max = maxval(inds)
        !======================================
        !  Count the l-m pairs that I (and others in my radial group) have.
        !  I.e., count the number of spherical harmonic modes owned by my group.
        !  Create an array of l's in processor and m's in processor
        lm_count = 0


        mpmin = mlb%min
        mpmax = mlb%max

        nmv = mpmax-mpmin+1
        Allocate(mv_temp(1:nmv))
        m_count = 1
        Do mp = mpmin, mpmax
            m = inds(mp)
            mv_temp(m_count) = m
            m_count = m_count+1
            Do l = m, l_max
                lm_count = lm_count+1
            Enddo
        Enddo

        l_min = minval(mv_temp)
        nlv = l_max-l_min+1
        Allocate(lv_temp(1:nlv))
        lv_temp(1) = l_min
        Do i = 2, nlv
            lv_temp(i) = lv_temp(i-1)+1
        Enddo

        !=============================================
        ! Order the l-m pairs by l-value (so that l varies the slowest and m the fastest)
        ! Ordering m within each l group is not crucial
        Allocate(l_lm_values(1:lm_count))
        Allocate(m_lm_values(1:lm_count))
        ind = 1
        Do i = 1, nlv
            l = lv_temp(i)
            Do mp = 1, nmv
                m = mv_temp(mp)
                If ( m .le. l) Then    ! This is one of our pairs
                    l_lm_values(ind) = l
                    m_lm_values(ind) = m
                    ind = ind+1
                Endif
            Enddo
        Enddo
        DeAllocate(lv_temp,mv_temp)


        !===========================================
        ! Build a reverse lookup table of mp values corresponding to each l-m combo
        Allocate(mp_lm_values(1:lm_count))
        mp_lm_values(:) = -1
        Do i = 1, lm_count
            m = m_lm_values(i)
            Do mp = mpmin, mpmax
                If (m .eq. inds(mp)) Then
                    mp_lm_values(i) = mp
                Endif
            Enddo
        Enddo


        !===========================================
        ! Decide how many pairs each processor gets.
        ! Assign pairs to each processor.
        ! Distributing pairs this way shouldn't be too bad,
        ! but it might be possible to better load-balance.
        N_radial_Cpus = comm%np
        lm_per = lm_count/N_Radial_Cpus
        lm_remainder = MOD(lm_count, N_Radial_Cpus)
        Allocate(lm_owner(1:lm_count))
        Allocate(num_lm(0:N_Radial_Cpus-1))
        num_lm(:) = lm_per                    ! Number of lm_pairs assigned to each processor
        ind = 1
        Do i = 1, N_Radial_Cpus
            lm_owner(ind:ind+lm_per-1) = i-1
            ind = ind+lm_per
            If ( i .le. lm_remainder) Then
                lm_owner(ind) = i-1
                ind = ind+1
                num_lm(i-1) = num_lm(i-1)+1
            Endif
        Enddo
        my_num_lm = num_lm(my_radial_rank)




        !===============================
        !  Identify the index range within the lm values that I own
        my_lm_min = -1
        my_lm_max = -1
        If (my_num_lm .gt. 0) Then
            If (my_radial_rank .eq. 0) Then
                my_lm_min = 1
            Else
                my_lm_min = SUM(num_lm(0:my_radial_rank-1))+1
            Endif
            If (num_lm(my_radial_rank) .gt. 0) Then
                my_lm_max = my_lm_min+num_lm(my_radial_rank)-1
            Endif
        Endif


        !=====================================
        ! Discern the range of l-values I have been assigned.
        ! Create an array of unique l-values

        If (my_num_lm .gt. 0) Then

            ind = 0
            l_lm_min = l_lm_values(my_lm_min)
            l_lm_max = l_lm_values(my_lm_max)
            my_nl_lm = l_lm_max-l_lm_min+1

            Allocate(my_lm_lval(1:my_nl_lm))
            Do i = 1, my_nl_lm
                my_lm_lval(i) = l_lm_min+(i-1)
            Enddo


            !==================================
            ! Count the m-values I own for each l-value
            Allocate(my_nm_lm(1:my_nl_lm))
            my_nm_lm(:) = 0
            Do i = 1, my_nl_lm
                l = my_lm_lval(i)
                Do ind = my_lm_min, my_lm_max
                    If (l_lm_values(ind) .eq. l) Then
                        my_nm_lm(i) = my_nm_lm(i)+1
                    Endif
                Enddo
            Enddo
        Else
            ! I own nothing - initialize arrays to reflect
            Allocate(my_lm_lval(1:1))
            my_lm_lval(1) = -1
            Allocate(my_nm_lm(1:1))
            my_nm_lm(1) = 0
        Endif

    End Subroutine LM_load_balance

End Module Load_Balance
