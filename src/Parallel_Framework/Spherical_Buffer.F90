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

Module Spherical_Buffer
    Use MPI_LAYER
    Use Parallel_Framework
    Use Structures
    Use Load_Balance
    Use General_MPI
    Implicit None
    Private

    Character*6 :: ifmt = '(i4.4)' ! Integer format for indicating processor tag numbers in output
    Type, Public :: SphericalBuffer
        ! The buffer object for buffer moving between spaces
        ! The buffer can advance configurations

        Type(rmcontainer4D), Allocatable :: s2b(:),s2a(:)

        Real*8, Allocatable :: p2b(:,:,:),p2a(:,:,:)        ! for dgemm
        Real*8, Allocatable :: p3a(:,:,:,:)  ! m/r/delta_theta/field
        Real*8, Allocatable :: p3b(:,:,:,:)  ! something
        Real*8, Allocatable :: p1b(:,:,:,:),p1a(:,:,:,:)  ! something

        Integer :: nf1a = 1
        Integer :: nf2a = 1
        Integer :: nf3a = 1
        Integer :: nf3b = 1
        Integer :: nf2b = 1
        Integer :: nf1b = 1


        Integer, Allocatable :: sdisp12(:), rdisp12(:)
        Integer, Allocatable :: scount12(:), rcount12(:)

        Integer, Allocatable :: sdisp21(:), rdisp21(:)
        Integer, Allocatable :: scount21(:), rcount21(:)

        Integer, Allocatable :: sdisp23(:), rdisp23(:)
        Integer, Allocatable :: scount23(:), rcount23(:)

        Integer, Allocatable :: sdisp32(:), rdisp32(:)
        Integer, Allocatable :: scount32(:), rcount32(:)

        Integer :: send_size12, recv_size12, send_size21, recv_size21
        Integer :: send_size23, recv_size23, send_size32, recv_size32


        Character*3 :: config
        ! When memory is not a consideration, it may be advantageous
        ! to keep the send/receive and/or the configuration buffers
        ! static in memory.  The default is to allocate/deallocate them
        ! throughout each iteration
        Logical :: dynamic_transpose_buffers = .true.
        Logical :: dynamic_config_buffers = .true.
        Logical :: pad_buffer = .false.

        Real*8, Allocatable :: recv_buff(:), send_buff(:)
        Integer :: max_recv, max_send

        !Additional information can be loaded into the transpose buffers.
        !This can timestep information, kill signals, or global walltime etc.
        !Following the cargo transposes 3b2b and 2b1b, each element of the cargo
        ! is replaced by the row-max (3b2b) or column-max (2b1b) value of that element
        Logical :: have_cargo = .false.
        Real*8, Allocatable :: cargo(:), tmp_cargo(:)
        Integer, Allocatable :: istop(:)
        Integer :: ncargo


        !scount12(0) => number I send to rank 0 when going FROM 1 TO 2
        !scount21(0) => number I send to rank 0 when going FROM 2 TO 1
        Contains

        Procedure :: Init  => Initialize_Spherical_Buffer
        Procedure :: construct => Allocate_Spherical_Buffer
        Procedure :: deconstruct => DeAllocate_Spherical_Buffer
        Procedure :: reform => advance_configuration
        Procedure :: set_buffer_sizes
        Procedure :: transpose_1a2a
        Procedure :: Compute_Packet_Sizes12
        Procedure :: transpose_2a3a  ! Move from 2a to 3a
        Procedure :: Compute_Packet_Sizes23
        Procedure :: transpose_3b2b
        Procedure :: transpose_2b1b
        Procedure :: write_space
        Procedure :: load_cargo
        Procedure :: unload_cargo
    End Type SphericalBuffer

Contains

    Subroutine load_cargo(self,vals)
        Implicit None
        Real*8, Intent(In) :: vals(1:)
        Class(SphericalBuffer) :: self
        self%cargo(1:self%ncargo) = vals(1:self%ncargo)
    End Subroutine load_cargo

    Subroutine unload_cargo(self,vals)
        Implicit None
        Real*8, Intent(InOut) :: vals(1:)
        Class(SphericalBuffer) :: self
        vals(1:self%ncargo) = self%cargo(1:self%ncargo)
    End Subroutine unload_cargo

    Subroutine Compute_Packet_Sizes23(self,numfields)
        Implicit None
        Class(SphericalBuffer) :: self
        Integer, Intent(In), Optional :: numfields
        Integer :: stemp1, rtemp1, my_max,gmax, np
        Integer :: nfs, p
        If (present(numfields)) Then
            nfs = numfields
        Else
            nfs = self%nf2a
        Endif
        np = pfi%rcomm%np
        Do p = 0, pfi%rcomm%np-1
            self%scount23(p) = (pfi%my_1p%delta) * (pfi%all_2p(p)%delta) * (pfi%my_3s%delta) * nfs
            self%rcount23(p) = (pfi%my_1p%delta) * (pfi%my_2p%delta) * (pfi%all_3s(p)%delta) * nfs
            self%scount23(p) = self%scount23(p)*2    ! times 2 for real/complex part
            self%rcount23(p) = self%rcount23(p)*2
        Enddo
        If (self%pad_buffer) Then
            ! adjustment for using alltoall vs. alltoallv
            ! Everyone in the row or column needs to have the same buffer size
            ! (hence the call to allreduce)
            stemp1 = maxval(self%scount23)
            rtemp1 = maxval(self%rcount23)
            my_max = max(stemp1,rtemp1)
            Call Global_Imax(my_max,gmax,pfi%rcomm)
            self%scount23(:) = gmax
            self%rcount23(:) = gmax
        Endif
        self%sdisp23(0) = 0
        self%rdisp23(0) = 0
        If (np .gt. 1) Then
            Do p = 1, np -1
                self%sdisp23(p) = self%sdisp23(p-1)+self%scount23(p-1)
                self%rdisp23(p) = self%rdisp23(p-1)+self%rcount23(p-1)
            Enddo
        Endif

        self%send_size23 = sum(self%scount23)
        self%recv_size23 = sum(self%rcount23)
    End Subroutine Compute_Packet_Sizes23


    Subroutine Set_Buffer_Sizes(self)
            Implicit None
            Class(SphericalBuffer) :: self
            Integer :: new_recv, new_send

            ! Start with 2a3a sizes
            self%max_send = sum(self%scount23)
            self%max_recv = sum(self%rcount23)

            ! Compare against 3b2b
            new_send = sum(self%scount32)
            new_recv = sum(self%rcount32)
            if (new_send .gt. self%max_send) then
                self%max_send = new_send
            endif
            if (new_recv .gt. self%max_recv) then
                self%max_recv = new_recv
            endif


            Allocate(self%recv_buff(1:self%max_recv))
            Allocate(self%send_buff(1:self%max_send))
    End Subroutine Set_Buffer_Sizes

    Subroutine Transpose_2a3a(self,extra_recv)
        Class(SphericalBuffer) :: self
        !Real*8, Allocatable :: send_buff(:), recv_buff(:)
        Integer :: np
        Integer :: imin, imax, jmin, jmax, kmin,kmax,ii,nf
        Integer :: i,f,j,p,k,k_ind,delf,delj
        Integer, Intent(In), Optional :: extra_recv
        Integer :: numalloc
        ! This is where we we move from theta, delta_r, delta_m
        !  to m, delta_r, delta_theta
        If (self%dynamic_transpose_buffers) Then


            Allocate(self%send_buff(1:self%send_size23))
        Endif
        !--- Not sure if this is good or bad, but copy out the bounds of the loop for now
        nf = self%nf2a
        kmin = pfi%my_3s%min
        kmax = pfi%my_3s%max
        !jmin = pfi%my_1p%min
        !jmax = pfi%my_1p%max

        jmin = 1
        jmax = pfi%my_1p%delta

        np = pfi%rcomm%np
        ! Possibly better ways to stripe, but for now, we will stripe each processors data all at once
        ! This means the the send buffer is accessed "naturally," but that the p2a array
        ! gets jumped around in
        delj = pfi%my_1p%delta


        Do p = 0, np -1
            ii = self%sdisp23(p)+1
            imin = pfi%all_2p(p)%min
            imax = pfi%all_2p(p)%max
            Do k = kmin, kmax
            Do f = 1,nf
            delf = (f-1)*delj*2
            Do j = jmin,jmax
            Do i = imin,imax

                self%send_buff(ii) = self%p2a(i,j+delf,k)
                self%send_buff(ii+1) = self%p2a(i,j+delf+delj,k)
                ii = ii+2
            Enddo
            Enddo
            Enddo
            Enddo
        Enddo

        Call self%deconstruct('p2a')
        If (self%dynamic_transpose_buffers) Allocate(self%recv_buff(1:self%recv_size23))

        self%recv_buff(:) = 0.0d0


        Call Standard_Transpose(self%send_buff, self%recv_buff, self%scount23, &
                self%sdisp23, self%rcount23, self%rdisp23, pfi%rcomm, self%pad_buffer)
        !--------------------------------------------------
        If (self%dynamic_transpose_buffers) DeAllocate(self%send_buff)

        !Here, we need self%construct('p3a',extra =nfextra or another number)

        If (present(extra_recv)) Then
            !Allocate an s2a buffer that is larger than normal
            numalloc = extra_recv+self%nf3a
            Call self%construct('p3a',numfields = numalloc)
        Else
            Call self%construct('p3a')
        Endif


        self%p3a(:,:,:,:) = 0.0d0    ! This is important because we are going to take an fft later (De-aliasing is implicit here because
        ! we only stripe in data of the de-aliased m's, but we need to make sure the higher m's are zero!
        !Stripe from the receive buff
        !Here we access the receive buffer in the 'natural' order
        imin = pfi%my_2p%min
        imax = pfi%my_2p%max
        jmin = pfi%my_1p%min
        jmax = pfi%my_1p%max
        Do p = 0, np -1
            ii = self%rdisp23(p)+1
            kmin = pfi%all_3s(p)%min
            kmax = pfi%all_3s(p)%max
            Do k = kmin, kmax
            Do f = 1,nf
            Do j = jmin,jmax
            Do i = imin,imax

                k_ind = pfi%inds_3s(k)*2+1  ! (real) m=0 stored in p3b(1,:,:,:) (img in p3b(2,:,:,:))

                self%p3a(k_ind,j,i,f)=self%recv_buff(ii)
                self%p3a(k_ind+1,j,i,f)=self%recv_buff(ii+1)

                ii = ii+2

            Enddo
            Enddo
            Enddo
            ! <----------- Here we want another loop if additional fields were added
            !              send counts, recv counts, send displ, and recv displ will need to be adjusted
            !   foff = self%nf3a
            !   nfextra = number of additional fields being sent
            !   p3a also needs to be bigger (see pseudo code above)

            !Do f = 1,nfextra
            !Do j = jmin,jmax
            !Do i = imin,imax

            !    k_ind = pfi%inds_3s(k)*2+1  ! (real) m=0 stored in p3b(1,:,:,:) (img in p3b(2,:,:,:))

            !    self%p3a(k_ind,j,i,foff+1)=self%recv_buff(ii)
            !    self%p3a(k_ind+1,j,i,foff+f)=self%recv_buff(ii+1)

            !    ii = ii+2

            !Enddo
            !Enddo
            !Enddo



            Enddo ! K



        Enddo
        self%config = 'p3a'
        If (self%dynamic_transpose_buffers) DeAllocate(self%recv_buff)


    End Subroutine Transpose_2a3a

    Subroutine Transpose_3b2b(self)
        Class(SphericalBuffer) :: self
        Integer :: np
        Integer :: imin, imax, jmin, jmax, kmin,kmax,ii,nf
        Integer :: i,f,j,p,k,k_ind
        Integer :: delf, delj


        ! This is where we we move from theta, delta_r, delta_m
        !  to m, delta_r, delta_theta
        If (self%dynamic_transpose_buffers) Then
            Allocate(self%send_buff(1:self%send_size32))
        Endif
        !write(6,*)'executing new transpose'
        !--- Not sure if this is good or bad, but copy out the bounds of the loop for now
        nf = self%nf3b

        jmin = pfi%my_1p%min
        jmax = pfi%my_1p%max

        np = pfi%rcomm%np


        !///////////////////
        !  Again stripe in the natural order of the send buffer
        imin = pfi%my_2p%min
        imax = pfi%my_2p%max
        Do p = 0, np -1
            ii = self%sdisp32(p)+1
            kmin = pfi%all_3s(p)%min
            kmax = pfi%all_3s(p)%max
            Do k = kmin, kmax
            Do f = 1,nf
            Do j = jmin,jmax
            Do i = imin,imax        ! interleave real and imaginary parts
                k_ind = pfi%inds_3s(k)*2+1  ! (real) m=0 stored in p3b(1,:,:,:) (img in p3b(2,:,:,:))

                self%send_buff(ii) = self%p3b(k_ind,j,i,f)! real part
                self%send_buff(ii+1) = self%p3b(k_ind+1,j,i,f)! complex part

                ii = ii+2
            Enddo
            Enddo
            Enddo
            Enddo
         if (self%have_cargo) Then
             self%send_buff(ii:ii+self%ncargo-1) = self%cargo(1:self%ncargo) !self%mrv
             ii = ii+self%ncargo ! 1
         Endif
        Enddo

        !/////////////////////////////////////
        Call self%deconstruct('p3b')
        If (self%dynamic_transpose_buffers) Allocate(self%recv_buff(1:self%recv_size32))
        self%recv_buff(:) = 0.0d0
        !----- This is where alltoall will be called

        Call Standard_Transpose(self%send_buff, self%recv_buff, self%scount32, &
                self%sdisp32, self%rcount32, self%rdisp32, pfi%rcomm, self%pad_buffer)

        !--------------------------------------------------
        If (self%dynamic_transpose_buffers) DeAllocate(self%send_buff)

        Call self%construct('p2b')        ! p2a and p2b can share the same buffer space... maybe just call this p2...



        delj = pfi%my_1p%delta


        !///////////////////////////////////////
        ! Read from the receive buffer in its natural order
        kmin = pfi%my_3s%min
        kmax = pfi%my_3s%max

        jmin = 1
        jmax = pfi%my_1p%delta

        Do p = 0, np -1
            ii = self%rdisp32(p)+1
            imin = pfi%all_2p(p)%min
            imax = pfi%all_2p(p)%max
            Do k = kmin, kmax
            Do f = 1,nf
                delf = (f-1)*delj*2
            Do j = jmin,jmax
            Do i = imin,imax
                !p2b needs to be reshaped
                !   self%p2b(i,j*2*f,k) -- since dgemm needs a 2D array

                self%p2b(i,j+delf ,k)   = self%recv_buff(ii)  ! real
                self%p2b(i,j+delf+delj,k) = self%recv_buff(ii+1)    ! complex
                ii = ii+2 ! added +2
            Enddo
            Enddo
            Enddo
            Enddo
            if (self%have_cargo) then
                self%tmp_cargo(1:self%ncargo) = self%recv_buff(ii:ii+self%ncargo-1)
                Do i = 1, self%ncargo
                    If ( self%tmp_cargo(i) .gt. self%cargo(i) ) self%cargo(i) = self%tmp_cargo(i)
                Enddo


                !self%mrv = max(self%mrv,self%cargo(1))  ! not quite sure how to handle mrv yet
                ii = ii+self%ncargo ! 1
            endif
        Enddo



        self%config = 'p2b'
        If (self%dynamic_transpose_buffers) DeAllocate(self%recv_buff)
    End Subroutine Transpose_3b2b

    Subroutine Transpose_2b1b(self)
        ! Go from Explicit_Part(IDX) configuration to the new Implicit%RHS configuration
        ! Communication is now done entirely within the radial group (processors that
        ! share a common set of ell-m values).
        Implicit None

        Integer :: r,l, mp, lp, indx, r_min, r_max, dr,  cnt,i, imi, rind
        Integer :: n1, n, nfields, offset, delta_r, rmin, rmax, np,p

        Real*8, Allocatable :: send_buff(:),recv_buff(:)
        Integer :: tnr, send_offset

        ! cargo information
        Integer :: inext, pcurrent

        Class(SphericalBuffer) :: self

        n1 = pfi%n1p
        np = pfi%ccomm%np

        !nfields = self%nf2b
        ! In some situations, we may not want to transmit all fields
        ! and the next buffer may have less space.
        ! Use the smaller number of fields from this config and the next
        ! to determine how many fields are sent
        nfields = MIN(self%nf2b,self%nf1b)


        !/////
        ! First loading, should be relatively straight forward
        ! We have nfields (nf2b)  and the data is dimensioned x%(l,r,field_num)

        ! Load the send array in the ordering of the l-m values




        Allocate(send_buff(1:self%send_size21))




        r_min = pfi%my_1p%min
        r_max = pfi%my_1p%max

        send_buff = 0.0d0
        pcurrent = 0
        inext = num_lm(0)
        send_offset = 1
        Do i = 1, lm_count
            mp = mp_lm_values(i)
            l = l_lm_values(i)
            !offset = 0

            Do n = 1, nfields
                Do imi = 1, 2
                Do r = r_min,r_max
                send_buff(send_offset) = self%s2b(mp)%data(l,r,imi,n)
                send_offset = send_offset+1
                EndDo
                Enddo
            Enddo
            If (i .eq. inext) Then
                If (self%have_cargo) Then
                    ! load the cargo values into the next buffer slot.
                    send_buff(send_offset:send_offset+self%ncargo-1) = self%cargo(1:self%ncargo)
                    send_offset = send_offset+self%ncargo

                    !send_buff(send_offset) = self%mrv
                    !send_offset = send_offset+1
                Endif
                !I think the next 3 lines work with the new cargo structure
                pcurrent = pcurrent+1
                inext = self%istop(pcurrent)
                If (i .lt. lm_count) send_offset = self%sdisp21(pcurrent) + 1
            Endif
        Enddo

        !///////////////////////////////////////


        Call self%deconstruct('s2b')


        Allocate(recv_buff(1:self%recv_size21))

        recv_buff = 0.0d0
        Call Standard_Transpose(send_buff, recv_buff, self%scount21, self%sdisp21, self%rcount21, &
            self%rdisp21, pfi%ccomm, self%pad_buffer)
        DeAllocate(send_buff)

        Call self%construct('p1b')
        ! Now, the receive striping needs a little mapping
        ! WPS are coupled, but Z,Btor, and Bpol are not
        ! Let's assume that those buffers are dimensioned: (r,real/imag,mode,field)
        ! We may want to modify this later on to mesh with linear equation structure

        self%p1b = 0.0d0
        Do p = 0, np - 1
            indx = self%rdisp21(p)+1

            r_min = pfi%all_1p(p)%min
            r_max = pfi%all_1p(p)%max

            dr = r_max - r_min
            ! Each processor in the radial group will have given me the same number of
            ! l-m combos in the same (correct) order
            cnt = 1
            Do lp = 1, my_num_lm
                Do n = 1, nfields
                self%p1b(r_min:r_max,1,cnt,n) = recv_buff(indx:indx+dr); indx = indx + dr + 1
                self%p1b(r_min:r_max,2,cnt,n) = recv_buff(indx:indx+dr); indx = indx + dr + 1
                Enddo
                cnt = cnt+1
            Enddo
            if (self%have_cargo) then
                self%tmp_cargo(1:self%ncargo) = recv_buff(indx:indx+self%ncargo-1)
                Do i = 1, self%ncargo
                    If ( self%tmp_cargo(i) .gt. self%cargo(i) ) self%cargo(i) = self%tmp_cargo(i)
                Enddo
                indx = indx+self%ncargo

                !self%mrv = max(self%mrv,recv_buff(indx))
                !indx = indx+1
            endif
        EndDo

        self%config='p1b'
      Deallocate(recv_buff)

    End Subroutine Transpose_2b1b

    Subroutine Transpose_1a2a(self,extra_recv)
        ! Go from implicit configuration (1 physical) to configuration 2 (spectral)
        Implicit None

        Integer :: r,l, mp, lp, indx, r_min, r_max, dr, cnt,i
        Integer :: n, nfields, offset, delta_r, rmin, rmax, np,p,rind
        Integer :: recv_offset, tnr
        Integer :: imi, numalloc
        Integer, Intent(In), Optional :: extra_recv
        Real*8, Allocatable :: send_buff(:),recv_buff(:)
        Integer :: inext, pcurrent

        Class(SphericalBuffer) :: self

        nfields = self%nf1a

        Allocate(send_buff(1:self%send_size12))


        ! Now, the receive striping needs a little mapping
        ! WPS are coupled, but Z,Btor, and Bpol are not
        ! Let's assume that those buffers are dimensioned: (r,real/imag,mode,field)
        ! We may want to modify this later on to mesh with linear equation structure

        np = pfi%ccomm%np
        Do p = 0, np - 1
            indx = self%sdisp12(p)+1
            r_min = pfi%all_1p(p)%min
            r_max = pfi%all_1p(p)%max
            dr = r_max - r_min
            ! Each processor in the radial group will have given me the same number of
            ! l-m combos in the same (correct) order
            cnt = 1
            Do lp = 1, my_num_lm
                Do n = 1, nfields
                send_buff(indx:indx+dr) = self%p1a(r_min:r_max,1,cnt,n) ; indx = indx + dr + 1
                send_buff(indx:indx+dr) = self%p1a(r_min:r_max,2,cnt,n) ; indx = indx + dr + 1
                Enddo
                cnt = cnt+1
            Enddo
        End Do
        Call self%deconstruct('p1a')

        ! Load the send array in the ordering of the l-m values
        indx = 1

        Allocate(recv_buff(1:self%recv_size12))


        Call Standard_Transpose(send_buff, recv_buff, self%scount12, self%sdisp12, &
             self%rcount12, self%rdisp12, pfi%ccomm, self%pad_buffer)
        DeAllocate(send_buff)

        If (present(extra_recv)) Then
            !Allocate an s2a buffer that is larger than normal
            numalloc = extra_recv+self%nf2a
            Call self%construct('s2a',numfields = numalloc)
        Else
            Call self%construct('s2a')
        Endif
        pcurrent = 0
        inext = num_lm(0)
        recv_offset = 1
        r_min = pfi%my_1p%min
        r_max = pfi%my_1p%max
        Do i = 1, lm_count
            mp = mp_lm_values(i)
            l = l_lm_values(i)


            Do n = 1, nfields
            Do imi = 1, 2
                Do r = r_min,r_max
                    self%s2a(mp)%data(l,r,imi,n) = recv_buff(recv_offset) ! +rind)
                    recv_offset = recv_offset+1
                EndDo
            Enddo
            Enddo
            If (i .eq. inext) Then
                pcurrent = pcurrent+1
                inext = self%istop(pcurrent)
                if (i .lt. lm_count) recv_offset = self%rdisp12(pcurrent)+1  ! <---- Added +1
            Endif
        Enddo


        self%config='s2a'
        Deallocate(recv_buff)

    End Subroutine Transpose_1a2a

    Subroutine Compute_Packet_Sizes12(self,numfields)
        Implicit None
        Class(SphericalBuffer) :: self
        Integer, Intent(In), Optional :: numfields
        Integer :: stemp1, rtemp1, my_max,gmax, np
        Integer :: nfs,p
        If (present(numfields)) Then
            nfs = numfields
        Else
            nfs = self%nf1a
        Endif
        np = pfi%ccomm%np
        Do p = 0, pfi%ccomm%np-1
            self%scount12(p) = (pfi%all_1p(p)%delta) * (my_num_lm)  * (nfs)*2
            self%rcount12(p) = (pfi%my_1p%delta)     * (num_lm(p)) * (nfs)*2
        Enddo
        If (self%pad_buffer) Then
            ! adjustment for using alltoall vs. alltoallv
            stemp1 = maxval(self%scount12)
            rtemp1 = maxval(self%rcount12)
            my_max = max(stemp1,rtemp1)
            Call Global_Imax(my_max,gmax,pfi%ccomm)
            self%scount12(:) = gmax
            self%rcount12(:) = gmax
        Endif

        self%sdisp12(0) = 0
        self%rdisp12(0) = 0
        If (np .gt. 1) Then
            Do p = 1, np -1
                self%sdisp12(p) = self%sdisp12(p-1)+self%scount12(p-1)
                self%rdisp12(p) = self%rdisp12(p-1)+self%rcount12(p-1)
            Enddo
        Endif
        self%send_size12 = sum(self%scount12)
        self%recv_size12 = sum(self%rcount12)


    End Subroutine Compute_Packet_Sizes12

    Subroutine Write_Space(self,extra_tag)
        Class(SphericalBuffer) :: self
        character*120, Optional, Intent(In) :: extra_tag
        Character*120 :: report_file,report_tag
        Character*10 :: gtag, rtag, ctag
        Integer :: report_unit = 500
        Integer :: i,j,k,f
        Integer :: imin,imax,jmin,jmax,kmin,kmax,nf
        Write(gtag,ifmt)pfi%gcomm%rank
        Write(rtag,ifmt)pfi%rcomm%rank
        Write(ctag,ifmt)pfi%ccomm%rank
        report_tag = 'g'//Trim(gtag)//'_r'//Trim(rtag)//'_c'//Trim(ctag)
        Select Case(self%config)

            Case ('p3a')
                imin = 1
                imax = pfi%n3p
                jmin = pfi%my_1p%min
                jmax = pfi%my_1p%max
                nf = self%nf3a
                kmin = pfi%my_2p%min
                kmax = pfi%my_2p%max
                if (present(extra_tag)) Then

                    report_file = 'parallel_framework/reports/workspace/p3a_'//Trim(report_tag)//'_'//Trim(extra_tag)
                else
                    report_file = 'parallel_framework/reports/workspace/p3a_'//report_tag
                endif
                Open(unit=report_unit,file = report_file,form='unformatted',status='replace')
                Write(report_unit)pfi%n3p
                Write(report_unit)pfi%my_1p%delta
                Write(report_unit)pfi%my_2p%delta
                Write(report_unit)nf
                Write(report_unit)((((self%p3a(i,j,k,f),i=imin,imax),j=jmin,jmax),k = kmin,kmax),f = 1, nf)
                Write(report_unit)pfi%my_1p%min
                Write(report_unit)pfi%my_1p%max
                Write(report_unit)pfi%my_2p%min
                Write(report_unit)pfi%my_2p%max
                Close(report_unit)
            Case ('p3b')
                imin = 1
                imax = pfi%n3p
                jmin = pfi%my_1p%min
                jmax = pfi%my_1p%max
                nf = self%nf3a
                kmin = pfi%my_2p%min
                kmax = pfi%my_2p%max
                if (present(extra_tag)) Then

                    report_file = 'parallel_framework/reports/workspace/p3b_'//Trim(report_tag)//'_'//Trim(extra_tag)
                else
                    report_file = 'parallel_framework/reports/workspace/p3b_'//report_tag
                endif
                Open(unit=report_unit,file = report_file,form='unformatted',status='replace')
                Write(report_unit)pfi%n3p
                Write(report_unit)pfi%my_1p%delta
                Write(report_unit)pfi%my_2p%delta
                Write(report_unit)nf
                Write(report_unit)((((self%p3b(i,j,k,f),i=imin,imax),j=jmin,jmax),k = kmin,kmax),f = 1, nf)
                Write(report_unit)pfi%my_1p%min
                Write(report_unit)pfi%my_1p%max
                Write(report_unit)pfi%my_2p%min
                Write(report_unit)pfi%my_2p%max
                Close(report_unit)
        End Select
    End Subroutine Write_Space

    Subroutine Initialize_Spherical_Buffer(self,report, field_count, &
                                                        config,dynamic_transpose, dynamic_config, &
                                            hold_cargo,padding,num_cargo)
        ! Buffer initialization
        ! Handles send/receive disp/counts
        Implicit None
        Integer :: np, p
        Logical, Intent(In), Optional :: report, dynamic_transpose,dynamic_config, hold_cargo, padding
        Integer, Intent(In), Optional :: field_count(3,2), num_cargo
        Character*120 :: report_tag
        Character*10 :: gtag, rtag, ctag
        Character*3, Intent(In), Optional :: config
        Integer :: stemp1, rtemp1

        Integer*4 :: gmax, my_max
        Class(SphericalBuffer) :: self
        If (present(report)) Then
            Write(gtag,ifmt)pfi%gcomm%rank
            Write(rtag,ifmt)pfi%rcomm%rank
            Write(ctag,ifmt)pfi%ccomm%rank
            report_tag = 'g'//Trim(gtag)//'_r'//Trim(rtag)//'_c'//Trim(ctag)
        Endif
        If (present(padding)) Then
            If (padding) Then
                self%pad_buffer = .true.
            Endif
        Endif
        If (present(dynamic_config)) Then
            self%dynamic_config_buffers = dynamic_config
        Endif
        If (present(dynamic_transpose)) Then
            self%dynamic_transpose_buffers = dynamic_transpose
        Endif

        If (present(config)) Then
            self%config = config
        Else
            self%config = 'p1a'    ! Physical space by default, configuration 1 by default
        Endif
        self%nf1a = field_count(1,1)
        self%nf2a = field_count(2,1)
        self%nf3a = field_count(3,1)
        self%nf3b = field_count(3,2)
        self%nf2b = field_count(2,2)
        self%nf1b = field_count(1,2)
        ! Going from configuration 1 to configuration 2
        ! Suppose we are going from r-in processor, modes distributed
        ! --- immediately following the implicit solve
        ! -- to r-distributed, delta_lm(column) in processor
        ! The message size I would send to rank r is
        !scount(r,1) = my_nl_lm*delta_r(r)*nfsend(1,2)
        ! where my_nl_lm is the number of l-m modes I hold during the solve
        ! and nfsend is fields sent+radial derivatives sent
        ! i.e. WSZ, and radial derivatives (P + radial derivatives do not need to be transposed)

        !Configuration 2 to 3
        !Next, suppose we are going from theta in processor, r distributed, m distributed
        ! --- to --- phi in-processor, r-distributed, theta-distributed
        ! The message size i would send to rank r is
        !scount23(r,1) = delta_r(my_row_rank)*delta_theta(column(r))*my_nm(local)*nfsend(2,3)

        np = pfi%rcomm%np
        Allocate(self%scount23(0:np-1))
        Allocate(self%rcount23(0:np-1))
        Allocate(self%sdisp23(0:np-1))
        Allocate(self%rdisp23(0:np-1))

        Call self%Compute_Packet_Sizes23()


        !////////////////////////////////////////////////////////////
        ! ------------ Now the reverse-----------------
        !Configuration 3 to 2
        ! if nf2a and nf3b were the same, rcount32 would equal scount23, and scount32 would equal rcount23
        Allocate(self%scount32(0:np-1))
        Allocate(self%rcount32(0:np-1))
        Allocate(self%sdisp32(0:np-1))
        Allocate(self%rdisp32(0:np-1))
        Do p = 0, pfi%rcomm%np-1
            self%rcount32(p) = (pfi%my_1p%delta) * (pfi%all_2p(p)%delta) * (pfi%my_3s%delta) * (self%nf3b)
            self%scount32(p) = (pfi%my_1p%delta) * (pfi%my_2p%delta) * (pfi%all_3s(p)%delta) * (self%nf3b)
            self%rcount32(p) = self%rcount32(p)*2    ! for real/complex parts
            self%scount32(p) = self%scount32(p)*2
        Enddo
        If (self%pad_buffer) Then
            ! adjustment for using alltoall vs. alltoallv
            stemp1 = maxval(self%scount32)
            rtemp1 = maxval(self%rcount32)
            my_max = max(stemp1,rtemp1)
            Call Global_Imax(my_max,gmax,pfi%rcomm)
            self%scount32(:) = gmax
            self%rcount32(:) = gmax
        Endif
        self%sdisp32(0) = 0
        self%rdisp32(0) = 0
        If (np .gt. 1) Then
            Do p = 1, np -1
                self%sdisp32(p) = self%sdisp32(p-1)+self%scount32(p-1)
                self%rdisp32(p) = self%rdisp32(p-1)+self%rcount32(p-1)
            Enddo
        Endif
        self%send_size32 = sum(self%scount32)
        self%recv_size32 = sum(self%rcount32)

        ! Suppose we are going from 3 to 2 (reverse of above)
        !  The message size I would send to rank r is
        ! scount(r,2) =  delta_r(my_row_rank)*delta_theta(column(my_rank))*nm(r)*nfsend(3,2)

        ! Configuration 1 to 2
        np = pfi%ccomm%np
        Allocate(self%scount12(0:np-1))
        Allocate(self%rcount12(0:np-1))
        Allocate(self%sdisp12(0:np-1))
        Allocate(self%rdisp12(0:np-1))

        Call self%Compute_Packet_Sizes12()

        ! Configuration 2 to 1
        !  Suppose we are going from 2 to 1 (reverse direction)
        ! The message size I would send to rank r is
        !scount(r,1) = nl_lm(r)*delta_r(my_row_rank)*nfsend(2,1)
        np = pfi%ccomm%np
        Allocate(self%scount21(0:np-1))
        Allocate(self%rcount21(0:np-1))
        Allocate(self%sdisp21(0:np-1))
        Allocate(self%rdisp21(0:np-1))
        Do p = 0, np-1
            self%rcount21(p) = (pfi%all_1p(p)%delta) * (my_num_lm)  * (self%nf2b)*2    ! 2 is because we split the complex into two real pieces
            self%scount21(p) = (pfi%my_1p%delta)     * (num_lm(p)) * (self%nf2b)*2
        Enddo
        If (self%pad_buffer) Then
            ! adjustment for using alltoall vs. alltoallv
            stemp1 = maxval(self%scount21)
            rtemp1 = maxval(self%rcount21)
            my_max = max(stemp1,rtemp1)
            Call Global_Imax(my_max,gmax,pfi%ccomm)
            self%scount21(:) = gmax
            self%rcount21(:) = gmax
        Endif

        self%sdisp21(0) = 0
        self%rdisp21(0) = 0
        If (np .gt. 1) Then
            Do p = 1, np -1
                self%sdisp21(p) = self%sdisp21(p-1)+self%scount21(p-1)
                self%rdisp21(p) = self%rdisp21(p-1)+self%rcount21(p-1)
            Enddo
        Endif
        self%send_size21 = sum(self%scount21)
        self%recv_size21 = sum(self%rcount21)


        If (present(hold_cargo)) Then
            If (hold_cargo) Then

                self%have_cargo = .true.
                Allocate(self%cargo(1:num_cargo))
                Allocate(self%tmp_cargo(1:num_cargo))
                self%ncargo = num_cargo
                ! modify the send and receive sizes for 3b2b
                np = pfi%rcomm%np
                self%scount32 = self%scount32+num_cargo
                self%rcount32 = self%rcount32+num_cargo
                ! Adjust the send and receive displacements
                Do p = 1, np-1
                    self%sdisp32(p) = self%sdisp32(p-1)+self%scount32(p-1)
                    self%rdisp32(p) = self%rdisp32(p-1)+self%rcount32(p-1)
                Enddo
                     self%send_size32 = self%send_size32+np*num_cargo
                     self%recv_size32 = self%recv_size32+np*num_cargo

                !Now do the same for 2b1b transposes
                np = pfi%ccomm%np
                self%scount21 = self%scount21+num_cargo
                self%rcount21 = self%rcount21+num_cargo
                Do p = 1, np -1
                    self%sdisp21(p) = self%sdisp21(p-1)+self%scount21(p-1)
                    self%rdisp21(p) = self%rdisp21(p-1)+self%rcount21(p-1)
                enddo
                     self%send_size21 = self%send_size21+np*num_cargo
                     self%recv_size21 = self%recv_size21+np*num_cargo


            Endif
        Endif



            ! Allocate and initialize the inext array
          Allocate(self%istop(0:np))    ! goes to np, but only np-1 used
             self%istop(np) = -1
          self%istop(0)  = num_lm(0)
          Do p = 1, np -1
              self%istop(p) = self%istop(p-1)+num_lm(p)
          Enddo


        !//  Allocate the static send/receive buffers if desired
        If (.not. self%dynamic_transpose_buffers) then
            Call self%set_buffer_sizes()
        Endif

    End Subroutine Initialize_Spherical_Buffer

    Subroutine DeAllocate_Spherical_Buffer(self,config,override)
        Class(SphericalBuffer) :: self
        Logical, Optional, Intent(In) :: override
        Logical :: free_config_memory
        Character*3, Intent(In) :: config
        Integer :: mn1
        Integer :: mx1,i

        free_config_memory = .true.
        If (.not. self%dynamic_config_buffers ) Then
            free_config_memory = .false.
        Endif
        if (present(override)) then
            free_config_memory = override
        endif
        If (free_config_memory) Then
        Select Case(config)
            Case('p1a')
                If (allocated(self%p1a)) Then
                    DeAllocate(self%p1a)
                Else
                    Write(6,*)'p1a does not appear to be allocated'
                Endif
            Case('p1b')
                If (allocated(self%p1b)) Then
                    DeAllocate(self%p1b)
                Else
                    Write(6,*)'p1b does not appear to be allocated'
                Endif
            Case('p2a')
                If (allocated(self%p2a)) Then
                    DeAllocate(self%p2a)
                Else
                    Write(6,*)'p2a does not appear to be allocated'
                Endif
            Case('p2b')
                If (allocated(self%p2b)) Then
                    DeAllocate(self%p2b)
                Else
                    Write(6,*)'p2b does not appear to be allocated'
                Endif

            Case('p3a')
                If (allocated(self%p3a)) Then
                    DeAllocate(self%p3a)
                Else
                    Write(6,*)'p3a does not appear to be allocated'
                Endif
            Case('p3b')
                If (allocated(self%p3b)) Then
                    DeAllocate(self%p3b)
                Else
                    Write(6,*)'p3b does not appear to be allocated'
                Endif
            Case('s2a')
                ! Appropriate for a triangular truncation
                If (allocated(self%s2a)) Then
                    mn1 = pfi%my_3s%min
                    mx1 = pfi%my_3s%max
                    Do i = mn1, mx1
                        If (Allocated(self%s2a(i)%data)) Then
                            DeAllocate(self%s2a(i)%data)
                        Endif
                    Enddo
                    DeAllocate(self%s2a)
                Else
                    Write(6,*)'s2a does not appear to be allocated.'
                Endif

            Case('s2b')
                ! Appropriate for a triangular truncation
                If (allocated(self%s2b)) Then
                    mn1 = pfi%my_3s%min
                    mx1 = pfi%my_3s%max
                    Do i = mn1, mx1
                        If (Allocated(self%s2b(i)%data)) Then
                            DeAllocate(self%s2b(i)%data)
                        Endif
                    Enddo
                    DeAllocate(self%s2b)
                Else
                    Write(6,*)'s2b does not appear to be allocated.'
                Endif

        End Select
        Endif
    End Subroutine DeAllocate_Spherical_Buffer

    Subroutine Allocate_Spherical_Buffer(self,config,numfields)
        Class(SphericalBuffer) :: self
        Character*3, Intent(In) :: config
        Integer, Intent(In), Optional :: numfields
        Integer :: mn1, mn2, mn3, mn4
        Integer :: mx1,mx2,mx3,mx4
        Integer :: i
        Select Case(config)
            Case('p1a')
                If (.not. Allocated(self%p1a)) Then
                    mn1 = 1
                    mx1 = pfi%n1p
                    mn2 = 1
                    mx2 = 2
                    mn3 = 1
                    mx3 = my_num_lm
                    mn4 = 1
                    If (present(numfields)) Then
                        mx4 = numfields
                    Else
                        mx4 = self%nf1a
                    Endif
                    Allocate(self%p1a(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))
                Endif
            Case('p1b')
                If (.not. Allocated(self%p1b)) Then
                    mn1 = 1
                    mx1 = pfi%n1p
                    mn2 = 1
                    mx2 = 2
                    mn3 = 1
                    mx3 = my_num_lm
                    mn4 = 1
                    If (present(numfields)) Then
                        mx4 = numfields
                    Else
                        mx4 = self%nf1b
                    Endif
                    Allocate(self%p1b(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))
                Endif

            Case('p2a')
                If (.not. Allocated(self%p2a)) Then
                    mn1 = 1
                    mx1 = pfi%n2p
                    mn2 = pfi%my_1p%min
                    mx2 = pfi%my_1p%max
                    mn4 = pfi%my_3s%min
                    mx4 = pfi%my_3s%max
                    If (present(numfields)) Then
                        mx3 = numfields
                    Else
                        mx3 = self%nf2a
                    Endif
                    mx3 = mx3*2*(mx2-mn2+1)
                    Allocate(self%p2a(mn1:mx1, 1:mx3, mn4:mx4))
                Endif
            Case('p2b')
                If (.not. Allocated(self%p2b)) Then
                    mn1 = 1
                    mx1 = pfi%n2p
                    mn2 = pfi%my_1p%min
                    mx2 = pfi%my_1p%max
                    mn4 = pfi%my_3s%min
                    mx4 = pfi%my_3s%max

                    If (present(numfields)) Then
                        mx3 = numfields
                    Else
                        mx3 = self%nf2b
                    Endif
                    mx3 = mx3*2*(mx2-mn2+1)
                    Allocate(self%p2b(mn1:mx1, 1:mx3, mn4:mx4))
                Endif
            Case('p3a')
                If (.not. Allocated(self%p3a)) Then
                    mn1 = 1
                    mx1 = pfi%n3p+2
                    mn2 = pfi%my_1p%min
                    mx2 = pfi%my_1p%max
                    mn3 = pfi%my_2p%min
                    mx3 = pfi%my_2p%max
                    mn4 = 1
                    If (present(numfields)) Then
                        mx4 = numfields
                    Else
                        mx4 = self%nf3a
                    Endif
                    !Write(6,*)'p3a -- mx4 is: ', mx4, numfields, self%nf3a
                    Allocate(self%p3a(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))
                Endif
            Case('p3b')
                If (.not. Allocated(self%p3b)) Then
                    mn1 = 1
                    mx1 = pfi%n3p+2 ! necessary for an in place transform
                    mn2 = pfi%my_1p%min
                    mx2 = pfi%my_1p%max
                    mn3 = pfi%my_2p%min
                    mx3 = pfi%my_2p%max
                    mn4 = 1
                    If (present(numfields)) Then
                        mx4 = numfields
                    Else
                        mx4 = self%nf3b
                    Endif
                    !Write(6,*)'Allocating p3b'
                    Allocate(self%p3b(mn1:mx1, mn2:mx2, mn3:mx3, mn4:mx4))
                Endif

            Case('s2a')
                If (.not. Allocated(self%s2a)) Then
                    ! We use a real array here instead of a complex array
                    ! Fields are striped rmin-rmax real then rmax+1-2*rmax imaginary in second index
                    ! Appropriate for a triangular truncation
                    ! --- partly why spherical and cartesian buffers will be separate
                    mn1 = pfi%my_3s%min
                    mx1 = pfi%my_3s%max

                    If (present(numfields)) Then
                        mx4 = numfields
                    Else
                        mx4 = self%nf2a
                    Endif

                    mn3 = pfi%my_1p%min
                    mx3 = pfi%my_1p%max

                    Allocate(self%s2a(mn1:mx1))
                    mx2 = maxval(pfi%inds_3s)    ! l_max = m_max
                    !Write(6,*)'mx4 is: ', mx4, numfields, self%nf2a
                    Do i = mn1, mx1
                        mn2 = pfi%inds_3s(i)        !l_min = m
                        Allocate(self%s2a(i)%data(mn2:mx2,mn3:mx3,1:2,1:mx4))
                    Enddo
                Endif
            Case('s2b')
                If (.not. Allocated(self%s2b)) Then
                    ! This is just s2b, but we are no longer complex, and
                    ! have fields, imaginary/real, and radius all in second index
                    ! Appropriate for a triangular truncation
                    ! --- partly why spherical and cartesian buffers will be separate
                    mn1 = pfi%my_3s%min
                    mx1 = pfi%my_3s%max

                    If (present(numfields)) Then
                        mx4 = numfields
                    Else
                        mx4 = self%nf2b
                    Endif


                    mn3 = pfi%my_1p%min
                    mx3 = pfi%my_1p%max

                    Allocate(self%s2b(mn1:mx1))
                    mx2 = maxval(pfi%inds_3s)    ! l_max = m_max

                    Do i = mn1, mx1
                        mn2 = pfi%inds_3s(i)        !l_min = m
                        Allocate(self%s2b(i)%data(mn2:mx2,mn3:mx3,1:2,1:mx4))
                    Enddo

                Endif

        End Select

    End Subroutine Allocate_Spherical_Buffer
    Subroutine Advance_Configuration(self,nextra_recv)
        Integer, Intent(In), Optional :: nextra_recv
        Class(SphericalBuffer) :: self
        Select Case(self%config)
            Case ('p2a')
                If (present(nextra_recv)) Then
                    Call self%transpose_2a3a(extra_recv = nextra_recv)
                Else
                    Call self%transpose_2a3a()
                Endif
            Case ('p3b')

                Call self%transpose_3b2b()

            Case ('s2b')

                Call self%transpose_2b1b()

            Case ('p1a')
                If (present(nextra_recv)) Then
                    Call self%transpose_1a2a(extra_recv = nextra_recv)
                Else
                    Call self%transpose_1a2a()
                Endif
        End Select
    End Subroutine Advance_Configuration
End Module Spherical_Buffer
