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

Module Fourier_Transform
    Implicit None
    Include 'fftw3.f'

    Interface FFT_To_Spectral
        Module Procedure r2c_ip4d_fftw
    End Interface

    Interface FFT_To_Spectral_RSC
        Module Procedure r2c_ip4d_fftw_rsc
    End Interface

    Interface FFT_To_Physical
        Module Procedure c2r_ip4d_fftw
    End Interface
Contains
    Subroutine Initialize_FFTs()
        Use Parallel_Framework, Only : pfi
        Implicit None
        ! Really just here for openmp init
        integer :: iret, nthread

        nthread = pfi%nthreads
#ifdef useomp
        if (nthread .gt. 1) Then
            call dfftw_init_threads(iret)
                !Note that when using MKL, iret will be 0 always.
                !iret = 0 means an error when using normal FFTW, but for MKL
                ! This routine is a wrapper that does nothing.
            !write(6,*)"iret is: ", iret
            ! comment
            !write(6,*)"FFTW planning with nthreads: ", nthread
            call dfftw_plan_with_nthreads(nthread)
        Endif
#endif
    End Subroutine Initialize_FFts

    Subroutine r2c_ip4D_fftw_rsc(x,plan)
        ! r2c: Real to Complex (forward transform)
        ! ip4D: In-place transform of 4D array
        ! fftw: FFTW
        ! rsc:  carries out rescaling (for use when dropping into spectral space)
        Real*8, Intent(InOut) :: x(:,:,:,:)
        Integer*8, Intent(In), Optional :: plan
        Integer*8 :: fresh_plan
        Integer :: n, howmany
        Integer :: inembed, istride, idist
        Integer :: onembed, ostride, odist
        Integer :: xshape(4)

        If (present(plan)) Then
            ! Plan exists - x assumed to keep memory location
            call dfftw_execute(plan)
        Else
            ! X will typically be deallocated and reallocated (unfortunately)
            ! So we will need to recreate the plan (because x's memory location may change)
            xshape = shape(x)
            n = xshape(1) -2 ! We assume nphi is even that the arrays has been padded by 2
            howmany = xshape(2)*xshape(3)*xshape(4)
            inembed = 0
            onembed = 0
            istride = 1
            ostride = 1
            ! idist is the distance between successive arrays to be transformed are stored
            ! odist is the distance between where successive results stored
            !    -- we assume the stride is 1 here, so this is just the length of the first dimension of the array
            idist = 2*(n/2+1)        ! In place transforms require extra padding - 2 extra for even n.  1 extra for odd n.
            odist = n/2+1            ! This is the size of the corresponding complex arrays
                                        ! At least it would be if we weren't doing this in place

            call dfftw_plan_many_dft_r2c(fresh_plan,1,n, howmany, &
                                                & x, inembed, istride, idist, &
                                                & x, onembed, ostride, odist, &
                                                & FFTW_ESTIMATE)
            call dfftw_execute(fresh_plan)
            call dfftw_destroy_plan(fresh_plan)
        Endif
        x(3:,:,:,:) = x(3:,:,:,:)/(n)
        x(1,:,:,:) = x(1,:,:,:)/(2*n)
    End Subroutine r2c_ip4D_fftw_rsc

    Subroutine r2c_ip4D_fftw(x,plan,rsc)
        ! r2c: Real to Complex (forward transform)
        ! ip4D: In-place transform of 4D array
        ! fftw : FFTW
        Real*8, Intent(InOut) :: x(:,:,:,:)
        Integer*8, Intent(In), Optional :: plan
        Integer*8 :: fresh_plan
        Integer :: n, howmany
        Integer :: inembed, istride, idist
        Integer :: onembed, ostride, odist
        Integer :: xshape(4)
        Logical, Intent(In), Optional :: rsc
        If (present(plan)) Then
            ! Plan exists - x assumed to keep memory location
            call dfftw_execute(plan)
        Else
            ! X will typically be deallocated and reallocated (unfortunately)
            ! So we will need to recreate the plan (because x's memory location may change)
            xshape = shape(x)
            n = xshape(1) -2 ! We assume nphi is even that the arrays has been padded by 2
            howmany = xshape(2)*xshape(3)*xshape(4)
            inembed = 0
            onembed = 0
            istride = 1
            ostride = 1
            ! idist is the distance between successive arrays to be transformed are stored
            ! odist is the distance between where successive results stored
            !    -- we assume the stride is 1 here, so this is just the length of the first dimension of the array
            idist = 2*(n/2+1)        ! In place transforms require extra padding - 2 extra for even n.  1 extra for odd n.
            odist = n/2+1            ! This is the size of the corresponding complex arrays
                                        ! At least it would be if we weren't doing this in place

            call dfftw_plan_many_dft_r2c(fresh_plan,1,n, howmany, &
                                                & x, inembed, istride, idist, &
                                                & x, onembed, ostride, odist, &
                                                & FFTW_ESTIMATE)
            call dfftw_execute(fresh_plan)
            call dfftw_destroy_plan(fresh_plan)
        Endif

    End Subroutine r2c_ip4D_fftw

    Subroutine c2r_ip4D_fftw(x,plan,rsc)
        ! r2c: Complex to Real (inverse transform)
        ! ip4D: In-place transform of 4D array
        ! fftw : FFTW
        Real*8, Intent(InOut) :: x(:,:,:,:)
        Logical, Intent(In), optional :: rsc
        Integer*8, Intent(In), Optional :: plan
        Integer*8 :: fresh_plan
        Integer :: n, howmany
        Integer :: inembed, istride, idist
        Integer :: onembed, ostride, odist
        Integer :: xshape(4)

        If (present(plan)) Then
            ! Plan exists - x assumed to keep memory location
            call dfftw_execute(plan)
        Else
            ! X will typically be deallocated and reallocated (unfortunately)
            ! So we will need to recreate the plan (because x's memory location may change)
            xshape = shape(x)
            n = xshape(1) -2 ! We assume nphi is even that the arrays has been padded by 2
            howmany = xshape(2)*xshape(3)*xshape(4)
            inembed = 0
            onembed = 0
            istride = 1
            ostride = 1
            ! idist is the distance between successive arrays to be transformed are stored
            ! odist is the distance between where successive results stored
            !    -- we assume the stride is 1 here, so this is just the length of the first dimension of the array
            odist = 2*(n/2+1)        ! In place transforms require extra padding - 2 extra for even n.  1 extra for odd n.
            idist = n/2+1            ! This is the size of the corresponding complex arrays
                                        ! At least it would be if we weren't doing this in place

            call dfftw_plan_many_dft_c2r(fresh_plan,1,n, howmany, &
                                                & x, inembed, istride, idist, &
                                                & x, onembed, ostride, odist, &
                                                & FFTW_ESTIMATE)
            call dfftw_execute(fresh_plan)
            call dfftw_destroy_plan(fresh_plan)
        Endif

    End Subroutine c2r_ip4D_fftw
End Module Fourier_Transform
