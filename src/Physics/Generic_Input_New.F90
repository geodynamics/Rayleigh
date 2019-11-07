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

Module Generic_Input

    Use RA_MPI_BASE
    Use Parallel_Framework
    Use Spherical_Buffer
    Use ProblemSize

    Implicit None
    Integer,Private :: genericinput_tag = 426

Contains

    Subroutine Read_Input(infile, field_ind, field)
        ! read the input from the given filename
        implicit none
        Integer, Intent(In) :: field_ind
        Character*120, Intent(In) :: infile
        Type(SphericalBuffer), Intent(InOut) :: field
        Type(SphericalBuffer) :: temp_field
        Type(IO_Buffer) :: inbuffer
        Integer :: fmode, funit, i, l_endian_tag, n_n,n_n_in version, pars(4)
        Integer, Allocatable :: nvals_in(:), nvals(:)
        Integer(kind=MPI_OFFSET_KIND) :: hdisp

        ! process 0 reads the first few parameters describing the file
        pars(:) = 0
        If (my_rank .eq. 0) Then
            Open(unit=15,file=trim(infile),form='unformatted',status='old',access='stream')
            Read(15)l_endian_tag   ! endian tag - FIXME: not currently used
            Read(15)version        ! version - not needed (yet)
            Read(15)fmode          ! fmode - tells us whether this file uses sparse (0) or dense (1) storage
            Read(15)n_n_in        ! number of coefficients
            Read(15)l_max_in

            If (fmode == 0) Then
                ! sparse storage in n, process 0 needs to read in the indices
                pars(2) = n_n
                Allocate(nvals_in(n_n))
                Read(15)(nvals_in(i),i=1,n_n)
            Else If (fmode != 1) then
                n_lmn = 0
                write(6,*)'Unknown generic input file mode: ', fmode
                write(6,*)'Acceptable modes are 0 (sparse) or 1 (full).'
            Endif
            Close(15)

            ! need to decide on actual n-values  TODO!
            ! nvals = SOMETHING
            ! DeAllocate(nvals_in)

            pars(1) = fmode
            pars(2) = n_n
            pars(3) = l_max_in
        Endif

        ! process 0 of each column receives file parameters from global process 0
        If (my_column_rank .eq. 0) Then
            Call MPI_Bcast(pars, 3, MPI_INTEGER, 0, pfi%rcomm%comm, ierr)
        Endif
        Call MPI_Bcast(pars, 3, MPI_INTEGER, 0, pfi%ccomm%comm, ierr)
        fmode    = pars(1)
        n_n      = pars(2)
        l_max_in = pars(3)

        hdisp = 20 
        if (fmode .ne. 0) hdisp = hdisp+4*n_n_in

        !~~~~~~~~~~~~~~~~>>>>   assume rank 0 takes care of the ns one read (nothing more than n_max)
        If (my_rank .ne. 0) Allocate(ns(1:n_n))
        If (my_column_rank .eq. 0) Then
            Call MPI_Bcast(nvals, n_n, MPI_INTEGER, 0, pfi%rcomm%comm, ierr)
        Endif
        Call MPI_Bcast(nvals, n_n, MPI_INTEGER, 0, pfi%ccomm%comm, ierr)

        Call inbuffer%Init(mpi_tag=checkpoint_tag,spectral=.true., & 
               cache_spectral = .true., spec_comp = .true., lmax_in = l_max_in, &
               r_indices = nvals, mode = 2)  ! TODO:  Need to add small bit of logic to io_buffer for nmax_in > nmax (qdisp)

        Call field%construct('s2b')
        temp_field%config = 's2b'
        Do mp = my_mp%min, my_mp%max
            temp_field%s2b(mp)%data(:,:,:,field_ind) = 0.0d0
        Enddo

        Call inbuffer%Read_Data(filename=infile, disp = hdisp)
        Call inbuffer%Grab_Data_spectral(temp_field%s2b,1)

        Call temp_field%reform()  ! move to p1b
          
        If (allocated(field%p1b)) Then  ! place the coeffs into a spherical buffer
            field%p1b(:,:,:,field_ind) = temp_field(:,:,:,1) 
        Else
            Write(6,*)'field%p1b not allocated in genericinput read_input!'
        Endif

    End Subroutine read_input

end module Generic_Input

