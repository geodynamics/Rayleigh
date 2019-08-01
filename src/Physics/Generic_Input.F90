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

module Generic_Input

  use RA_MPI_BASE
  use Parallel_Framework
  use Spherical_Buffer
  use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values, &
                           my_lm_min, my_nl_lm, my_nm_lm, my_lm_lval, my_lm_max, &
                           lm_count, lm_owner

  use SendReceive
  use ProblemSize

  implicit none

  integer,private :: genericinput_tag = 426

contains

  subroutine read_input(filename, field_ind, field)
    ! read the input from the given filename
    implicit none
    integer, intent(in) :: field_ind
    character*120, intent(in) :: filename
    type(SphericalBuffer), intent(inout) :: field
    
    integer :: l_endian_tag, version, fmode, n_lmn, l_l_max, l_n_max, k_lm_owner, l_max_n
    integer :: i, j, k, l, m, n, lm, mp, j1, jc, jf, jn, jo, ji, p, col_np, col_mp_min
    integer :: start_count, end_count, my_lmn_count, total_lmn_count, col_lmn_n, col_lmn_i, l_lm_count
    integer :: int_size, real_size, offset, ierr, funit
    integer(kind=MPI_OFFSET_KIND) :: disp1, disp2
    integer, dimension(3) :: pars
    integer, dimension(MPI_STATUS_SIZE) :: mstatus
    integer, allocatable, dimension(:) :: ls, ms, ns
    integer, allocatable, dimension(:) :: proc_lmn_count, col_lmn_count, col_lmn_ind1, col_coeffs_ind1, sendarri
    integer, allocatable, dimension(:,:) :: my_lmn_inds, lmn_inds, sendarri2
    real*8, allocatable, dimension(:,:) :: my_lmn_coeffs, lmn_coeffs, col_coeffs, sendarr2

    ! set up some sizes of mpi types
    call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, real_size, ierr)
    call MPI_TYPE_EXTENT(MPI_INTEGER, int_size, ierr)

    ! process 0 reads the first few parameters describing the file
    pars(:) = 0
    if (my_rank .eq. 0) then
      open(unit=15,file=trim(filename),form='unformatted',status='old',access='stream')
      read(15)l_endian_tag   ! endian tag - FIXME: not currently used
      read(15)version        ! version - not needed (yet)
      read(15)fmode          ! fmode - tells us whether this file uses sparse (0) or dense (1) storage
      pars(1) = fmode
      if (fmode == 0) then
        ! sparse storage, process 0 needs to read in the indices
        read(15)n_lmn        ! number of coefficients
        pars(2) = n_lmn
        allocate(ls(n_lmn), ms(n_lmn), ns(n_lmn))
        read(15)(ns(i),i=1,n_lmn)
        read(15)(ls(i),i=1,n_lmn)
        read(15)(ms(i),i=1,n_lmn)
      else if (fmode == 1) then
        ! dense storage, no indices to read but need to know the number of coefficients
        read(15)l_n_max      ! n_max
        pars(2) = l_n_max
        read(15)l_l_max      ! l_max
        pars(3) = l_l_max
      else
        n_lmn = 0
        write(6,*)'Unknown generic input file mode: ', fmode
        write(6,*)'Acceptable modes are 0 (sparse) or 1 (full).'
      endif
      close(15)
    end if

    ! process 0 of each column receives file parameters from global process 0
    if (my_column_rank .eq. 0) then
      call MPI_Bcast(pars, 3, MPI_INTEGER, 0, pfi%rcomm%comm, ierr)
    end if
    call MPI_Bcast(pars, 3, MPI_INTEGER, 0, pfi%ccomm%comm, ierr)
    fmode = pars(1)

    ! get min mp indices that this column owns
    col_mp_min = pfi%my_3s%min

    if (fmode == 0) then
      n_lmn = pars(2)       ! expected total number of coefficients

      if (my_column_rank .eq. 0) then
        ! sparse file storage... need to only read coefficients that exist

        if (my_rank .ne. 0) then
          allocate(ls(n_lmn), ms(n_lmn), ns(n_lmn))
        end if
        ! broadcast the ms, ls and ms to all processes with column rank 0
        call MPI_Bcast(ns, n_lmn, MPI_INTEGER, 0, pfi%rcomm%comm, ierr)
        call MPI_Bcast(ls, n_lmn, MPI_INTEGER, 0, pfi%rcomm%comm, ierr)
        call MPI_Bcast(ms, n_lmn, MPI_INTEGER, 0, pfi%rcomm%comm, ierr)
        ! now all column rank 0 processes know what indices are in the file

        ! next we need to identify which ms this column owns and how many coefficients there are of each m
        ! here we assume that the generic input file is indexed (n,l,m) and ordered using column-major ordering
        ! (i.e. m varies last, n varies fastest)

        allocate(col_lmn_count(pfi%my_3s%delta))
        allocate(col_lmn_ind1(pfi%my_3s%delta))
        col_lmn_count(:) = 0  ! how many coeffs there are at each m we find
        col_lmn_ind1(:) = -1 ! starting index of this m value in file (assuming ordering here!)
        do i = 1, pfi%my_3s%delta
          mp = col_mp_min + i - 1
          m = pfi%inds_3s(mp)
          do j = 1, n_lmn
            if (ms(j) .eq. m) then
              ! NOTE: here we'll automatically drop ms(j) > l_max and ms(j) < 0 but can still pick up arbitrary ls(j) and ns(j)
              ! because we want reading in (below) to be contiguous so dropping ls and ns isn't easy at this stage
              col_lmn_count(i) = col_lmn_count(i) + 1
              if (col_lmn_ind1(i) .eq. -1) then
                col_lmn_ind1(i) = j
              end if
            end if
          end do
        end do
        total_lmn_count = sum(col_lmn_count)   ! the total number of coefficients 
                                               ! we're expecting to read in is sum(col_lmn_count)


        ! now we know which coefficients we want for our column we need to open the file and read them (and only them) in
        allocate(col_coeffs(total_lmn_count,2))
        allocate(col_coeffs_ind1(pfi%my_3s%delta))  ! just a cumulative sum of col_lmn_count
        col_coeffs(:,:) = 0.0
        col_coeffs_ind1(:) = -1
        ! open the file
        call MPI_FILE_OPEN(pfi%rcomm%comm, trim(filename), &
                           MPI_MODE_RDONLY, MPI_INFO_NULL, &
                           funit, ierr)
        if (ierr .eq. 0) then
          ! calculate the offset of the integer header
          disp1 = (4 + 3*n_lmn)*int_size  ! 4 parameters then 3*n_lmn indices
          offset = 1
          ! loop over all the values of m that we know this column owns
          ! and read in all coeffs associated with that m value
          ! NOTE: assuming that real and imaginary parts are split into separate sections of the file
          do i = 1, pfi%my_3s%delta
            if (col_lmn_count(i) .gt. 0) then
              col_coeffs_ind1(i) = offset
              disp2 = disp1 + (col_lmn_ind1(i)-1)*real_size
              call MPI_FILE_SEEK(funit, disp2, MPI_SEEK_SET, ierr)
              call MPI_FILE_READ(funit, col_coeffs(offset,1), col_lmn_count(i), &
                                 MPI_DOUBLE_PRECISION, &
                                 mstatus, ierr)  ! real
              disp2 = disp1 + (n_lmn+col_lmn_ind1(i)-1)*real_size
              call MPI_FILE_SEEK(funit, disp2, MPI_SEEK_SET, ierr)
              call MPI_FILE_READ(funit, col_coeffs(offset,2), col_lmn_count(i), &
                                 MPI_DOUBLE_PRECISION, &
                                 mstatus, ierr) ! imaginary
              offset = offset + col_lmn_count(i)
            end if
          end do
        else
          write(6,*)'Error opening generic input file: ', pfi%rcomm%rank
        end if
        call MPI_FILE_CLOSE(funit, ierr)

        ! We've now read in all the coefficients for this column but they're currently ordered in (n,l,m) sets such that n varies
        ! fastest.  This isn't very useful for distributing them to the processes of the column as they are subdivided in lm pairs
        ! with m varying fastest non-sequentially.  Let's reorder all the coefficients so they can be packaged for the processes.

        col_np = pfi%ccomm%np         ! number of processes in this column
        allocate(proc_lmn_count(0:col_np-1))
        allocate(lmn_inds(total_lmn_count,2))   ! each mode gets its lm index (:,1) and n index (:,2) stored as a pair
        allocate(lmn_coeffs(total_lmn_count,2)) ! size may be an overestimate as we drop ls and ns that are not in range
        proc_lmn_count(:) = 0     ! how many total (n,l,m) combinations rank p of a col owns
        lmn_inds(:,:) = -1        ! indices of file coeff values in lm ordering
        lmn_coeffs(:,:) = 0.0     ! re-ordered coefficients, now m is ordered according to mp ordering (not necessarily sequential)
        ! loop over the lm modes in the order they are stored and partitioned and reorder those we've just read from the file to
        ! match
        jn = 1
        do k = 1, lm_count
          k_lm_owner = lm_owner(k) ! because of how each column is partitioned this will be sequential with the order of lms
          l = l_lm_values(k)
          mp = mp_lm_values(k)
          i = mp - col_mp_min + 1
          j1 = col_lmn_ind1(i)     ! first index of this m value in the file
          if (j1 .gt. 0) then
            jc = col_lmn_count(i)  ! number of coefficients with this m value in the file 
            do ji = 1, jc
              j = j1 + ji -1
              if ((ls(j) .eq. l) .and. (ns(j) .le. n_r)) then
                ! filter out extraneous ls and ns
                proc_lmn_count(k_lm_owner) = proc_lmn_count(k_lm_owner) + 1  ! increment count for this process (NOTE: assumed sequential!)
                lmn_inds(jn,1) = k                    ! record lm index of coeff
                lmn_inds(jn,2) = ns(j)+1              ! record n index of coeff
                jo = col_coeffs_ind1(i) + ji - 1
                lmn_coeffs(jn,1) = col_coeffs(jo,1)    ! record coeffs
                lmn_coeffs(jn,2) = col_coeffs(jo,2)
                jn = jn + 1
              end if
            end do
          end if
        end do
        ! because we've now gone through lm pairs in the order they're partitioned, lmn_inds and lmn_coeffs should be contiguous 
        ! by process

        deallocate(ls, ms, ns)
        deallocate(col_coeffs)
        deallocate(col_coeffs_ind1)
        deallocate(col_lmn_count)
        deallocate(col_lmn_ind1)

      end if

      ! communicate the number of coeffs
      if (my_column_rank .eq. 0) then
        my_lmn_count = proc_lmn_count(0)
        ! send sizes
        allocate(sendarri(1))
        do p = 1, col_np-1
          sendarri = proc_lmn_count(p)
          call send(sendarri, dest=p, tag=genericinput_tag, grp=pfi%ccomm)
        end do
        deallocate(sendarri)

      else
        allocate(sendarri(1))
        call receive(sendarri, source=0, tag=genericinput_tag, grp=pfi%ccomm)
        my_lmn_count = sendarri(1)
        deallocate(sendarri)
      end if

      ! allocate arrays to temporarily store the indices and coeffs
      allocate(my_lmn_inds(my_lmn_count,2))
      allocate(my_lmn_coeffs(my_lmn_count,2))

      ! communicate the indices
      if (my_column_rank .eq. 0) then

        my_lmn_inds(:,:) = lmn_inds(1:proc_lmn_count(0),:)
        start_count = proc_lmn_count(0)
        do p = 1, col_np-1
          if (proc_lmn_count(p) .gt. 0) then
            allocate(sendarri2(proc_lmn_count(p),2))
            end_count = start_count + proc_lmn_count(p)
            sendarri2(1:proc_lmn_count(p),:) = lmn_inds(start_count+1:end_count,:)
            call send(sendarri2, dest=p, tag=genericinput_tag, grp = pfi%ccomm)
            deallocate(sendarri2)
            start_count = end_count
          end if
        end do
        deallocate(lmn_inds)
      else
        if (my_lmn_count .gt. 0) then
          call receive(my_lmn_inds, source=0, tag=genericinput_tag, grp=pfi%ccomm)
        end if
      end if

      ! communicate the coeffs
      if (my_column_rank .eq. 0) then

        my_lmn_coeffs(:,:) = lmn_coeffs(1:proc_lmn_count(0),:)
        start_count = proc_lmn_count(0)
        do p = 1, col_np-1
          if (proc_lmn_count(p) .gt. 0) then
            allocate(sendarr2(proc_lmn_count(p),2))
            end_count = start_count + proc_lmn_count(p)
            sendarr2(1:proc_lmn_count(p),:) = lmn_coeffs(start_count+1:end_count,:)
            call send(sendarr2, dest=p, tag=genericinput_tag, grp = pfi%ccomm)
            deallocate(sendarr2)
            start_count = end_count
          end if
        end do
        deallocate(lmn_coeffs)
        deallocate(proc_lmn_count)
      else
        if (my_lmn_count .gt. 0) then
          call receive(my_lmn_coeffs, source=0, tag=genericinput_tag, grp=pfi%ccomm)
        end if
      end if
      
      if (allocated(field%p1b)) then
        ! place the coeffs into a spherical buffer
        field%p1b(:,:,:,field_ind) = 0.0
        do i = 1, my_lmn_count
          lm = my_lmn_inds(i,1) - my_lm_min + 1
          n = my_lmn_inds(i,2)
          field%p1b(n, 1, lm, field_ind) = my_lmn_coeffs(i, 1)
          field%p1b(n, 2, lm, field_ind) = my_lmn_coeffs(i, 2)
        end do
      else
        write(6,*)'field%p1b not allocated in genericinput read_input!'
      end if
      deallocate(my_lmn_inds)
      deallocate(my_lmn_coeffs)

    else if (fmode == 1) then
      l_n_max = pars(2)
      l_l_max = pars(3)
      n_lmn = (l_n_max+1)*(l_l_max+1)*(l_l_max+2)/2
      l_max_n = min(l_n_max, n_r - 1)

      if (my_column_rank .eq. 0) then

        total_lmn_count = 0
        do i = 1, pfi%my_3s%delta
          mp = col_mp_min + i - 1
          m = pfi%inds_3s(mp)
          total_lmn_count = total_lmn_count + (l_n_max+1)*max(l_l_max - m + 1, 0)
        end do

        ! now we know how many coefficients we want for our column we need to open the file and read them (and only them) in
        allocate(col_coeffs(total_lmn_count,2))
        allocate(col_coeffs_ind1(pfi%my_3s%delta))  ! just a cumulative sum of col_lmn_ns
        col_coeffs(:,:) = 0.0
        col_coeffs_ind1(:) = -1
        ! open the file
        call MPI_FILE_OPEN(pfi%rcomm%comm, trim(filename), &
                           MPI_MODE_RDONLY, MPI_INFO_NULL, &
                           funit, ierr)
        if (ierr .eq. 0) then
          ! calculate the offset of the integer header
          disp1 = 5*int_size  ! just 5 parameters
          offset = 1
          ! loop over all the values of m that we know this column owns
          ! and read in all coeffs associated with that m value
          ! NOTE: assuming that real and imaginary parts are split into separate sections of the file
          do i = 1, pfi%my_3s%delta
            mp = col_mp_min + i - 1
            m = pfi%inds_3s(mp)
            col_lmn_n = (l_n_max+1)*max(l_l_max - m + 1, 0)
            if (col_lmn_n .gt. 0) then
              col_coeffs_ind1(i) = offset
              col_lmn_i = (l_n_max+1)*m*(2*l_l_max + 3 - m)/2
              disp2 = disp1 + col_lmn_i*real_size
              call MPI_FILE_SEEK(funit, disp2, MPI_SEEK_SET, ierr)
              call MPI_FILE_READ(funit, col_coeffs(offset,1), col_lmn_n, &
                                 MPI_DOUBLE_PRECISION, &
                                 mstatus, ierr)  ! real
              disp2 = disp1 + (n_lmn+col_lmn_i)*real_size
              call MPI_FILE_SEEK(funit, disp2, MPI_SEEK_SET, ierr)
              call MPI_FILE_READ(funit, col_coeffs(offset,2), col_lmn_n, &
                                 MPI_DOUBLE_PRECISION, &
                                 mstatus, ierr) ! imaginary
              offset = offset + col_lmn_n
            end if
          end do
        else
          write(6,*)'Error opening generic input file: ', pfi%rcomm%rank
        end if
        call MPI_FILE_CLOSE(funit, ierr)

        ! We've now read in all the coefficients for this column but they're currently ordered in (n,l,m) sets such that n varies
        ! fastest.  This isn't very useful for distributing them to the processes of the column as they are subdivided in lm pairs
        ! with m varying fastest non-sequentially.  Let's reorder all the coefficients so they can be packaged for the processes.

        col_np = pfi%ccomm%np     ! number of processes in this column
        allocate(proc_lmn_count(0:col_np-1))
        allocate(lmn_coeffs(total_lmn_count,2)) ! size may be an overestimate as we drop ls and ns that are not in range
        proc_lmn_count(:) = 0     ! how many total (n,l,m) combinations rank p of a col owns
        lmn_coeffs(:,:) = 0.0     ! re-ordered coefficients, now m is ordered according to mp ordering (not necessarily sequential)
        ! loop over the lm modes in the order they are stored and partitioned and reorder those we've just read from the file to
        ! match
        jn = 1
        do k = 1, lm_count
          k_lm_owner = lm_owner(k) ! because of how each column is partitioned this will be sequential with the order of lms
          m = m_lm_values(k)
          l = l_lm_values(k)
          if ((m .lt. l_l_max + 1) .and. (l .lt. l_l_max + 1)) then
            mp = mp_lm_values(k)
            i = mp - col_mp_min + 1
            j1 = col_coeffs_ind1(i) + (l - m)*(l_n_max + 1)
            proc_lmn_count(k_lm_owner) = proc_lmn_count(k_lm_owner) + l_max_n + 1
            lmn_coeffs(jn:jn+l_max_n,1) = col_coeffs(j1:j1+l_max_n,1)
            lmn_coeffs(jn:jn+l_max_n,2) = col_coeffs(j1:j1+l_max_n,2)
            jn = jn + l_max_n + 1
          end if
        end do
        ! because we've now gone through lm pairs in the order they're partitioned, lmn_coeffs should be contiguous 
        ! by process

        deallocate(col_coeffs)
        deallocate(col_coeffs_ind1)
      end if
        
      ! communicate the number of coeffs
      if (my_column_rank .eq. 0) then
        my_lmn_count = proc_lmn_count(0)
        ! send sizes
        allocate(sendarri(1))
        do p = 1, col_np-1
          sendarri = proc_lmn_count(p)
          call send(sendarri, dest=p, tag=genericinput_tag, grp=pfi%ccomm)
        end do
        deallocate(sendarri)

      else
        allocate(sendarri(1))
        call receive(sendarri, source=0, tag=genericinput_tag, grp=pfi%ccomm)
        my_lmn_count = sendarri(1)
        deallocate(sendarri)
      end if

      ! allocate an array to temporarily store the coeffs
      allocate(my_lmn_coeffs(my_lmn_count,2))

      ! communicate the coeffs
      if (my_column_rank .eq. 0) then

        my_lmn_coeffs(:,:) = lmn_coeffs(1:proc_lmn_count(0),:)
        start_count = proc_lmn_count(0)
        do p = 1, col_np-1
          if (proc_lmn_count(p) .gt. 0) then
            allocate(sendarr2(proc_lmn_count(p),2))
            end_count = start_count + proc_lmn_count(p)
            sendarr2(1:proc_lmn_count(p),:) = lmn_coeffs(start_count+1:end_count,:)
            call send(sendarr2, dest=p, tag=genericinput_tag, grp = pfi%ccomm)
            deallocate(sendarr2)
            start_count = end_count
          end if
        end do
        deallocate(lmn_coeffs)
        deallocate(proc_lmn_count)
      else
        if (my_lmn_count .gt. 0) then
          call receive(my_lmn_coeffs, source=0, tag=genericinput_tag, grp=pfi%ccomm)
        end if
      end if
      
      if (allocated(field%p1b)) then
        ! place the coeffs into a spherical buffer
        field%p1b(:,:,:,field_ind) = 0.0
        l_lm_count = my_lmn_count/(l_max_n+1)
        do k = 1, l_lm_count
          field%p1b(1:l_max_n+1, 1, k, field_ind) = my_lmn_coeffs((k-1)*(l_max_n+1)+1:k*(l_max_n+1), 1)
          field%p1b(1:l_max_n+1, 2, k, field_ind) = my_lmn_coeffs((k-1)*(l_max_n+1)+1:k*(l_max_n+1), 2)
        end do
      else
        write(6,*)'field%p1b not allocated in genericinput read_input!'
      end if
      deallocate(my_lmn_coeffs)


    end if

  end subroutine read_input


end module Generic_Input

