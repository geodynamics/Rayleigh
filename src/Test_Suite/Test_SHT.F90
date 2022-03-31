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


Module Test_SHT
    !/////////////////////////////////////////////////
    !  Spherical Harmonic Transform Testing Module
    Use ProblemSize
    Use Fields
    Use Parallel_Framework
    Use General_MPI
    Use Fourier_Transform
    Use Legendre_Transforms, Only : Legendre_Transform
    Use SendReceive
    Implicit None
    Integer :: ntest_transform_loops = 1
    Real*8 :: passing_tolerance = 1.0d-10
Contains

    Subroutine Test_Spherical_Transforms()
        Call test_LT_FFT()
    End Subroutine Test_Spherical_Transforms

    Subroutine Test_Legendre_Transforms()
        Call single_mode_LT()
        Call single_mode_loops_LT()
        Call white_noise_LT()
    End Subroutine Test_Legendre_Transforms

    Subroutine single_mode_LT()
        ! -intialize single (l,m) mode in spectral space
        ! -build analytic solution in physical space
        ! -compute spectral -> physical and compute the error
        ! -compute spectral -> physical -> spectral and compute the error
        Implicit None
        Integer :: colrank, rowrank, nf=3, lval, mval
        Integer :: mp, l, m, r, f, ind, my_Nr
        Integer :: fcount(3,2)
        Real*8 :: mxdiff_spec_to_phys, mxdiff_spec_to_spec, diff, diff1, diff2, ans, norm, amp
        Real*8 :: mxdiff_spec_to_phys_reduced, mxdiff_spec_to_spec_reduced
        Real*8, allocatable :: to_phys_norm(:), to_spec_norm(:), true_soln(:)
        type(SphericalBuffer) :: mytest

        if ((pfi%rcomm%rank .eq. 0) .and. (pfi%ccomm%rank .eq. 0)) then
            write(6,*) '////////////////////////////////////////////////////////////'
            write(6,*) 'Test Legendre Transform: single mode, forward/inverse LT'
        endif

        fcount(:,:) = nf

        ! the FFT normalizations are built into the normalizations used in the Legendre
        ! transform. build the FFT normalizations here so they can be accounted for later
        !     these are set/defined in Legendre_Polynomials.F90 (search PTS or STP)
        allocate(to_phys_norm(0:l_max), to_spec_norm(0:l_max), true_soln(1:n_theta))
        to_phys_norm(0) = 1.0d0  ! m=0
        to_phys_norm(1:) = 0.5d0 ! m/=0
        to_spec_norm(0) = 1.0d0/(n_phi)    ! m=0
        to_spec_norm(1:) = 1.0d0/(n_theta) ! m/=0

        ! allocate fields in spectral space
        Call mytest%init(field_count=fcount, config='s2a')
        Call mytest%construct('s2a')

        ! build true solution for (l,m) = (2,1)
        lval = 2; mval = 1
        true_soln(:) = sqrt(5.0d0/(4.0d0*pi*6.0d0))*(-3.0d0*costheta(:)*sintheta(:))

        ! initialize specific (l,m) values
        do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            norm = 1.0d0/to_phys_norm(m) ! "undo" FFT normalization of to-physical direction
            do l=m,l_max
                mytest%s2a(mp)%data(l,:,:,:) = 0.0d0
                if ((l .eq. lval) .and. (m .eq. mval)) then
                    do r=my_r%min,my_r%max
                        do f=1,nf
                            mytest%s2a(mp)%data(l,r,1,f) = (2.0d0*f+1.0d0)*radius(r)*norm
                            mytest%s2a(mp)%data(l,r,2,f) = -(3.0d0*f*f+2.0d0)*radius(r)**2*norm
                        enddo
                    enddo
                endif
            enddo
        enddo

        ! build physical space
        Call mytest%construct('p2a')

        Call Legendre_Transform(mytest%s2a,mytest%p2a) ! LT to physical

        Call Legendre_Transform(mytest%p2a,mytest%s2a) ! LT back to spectral

        ! zero out l_max mode (done in production runs too), undo FFT norm of to-spectral
        do mp=my_mp%min,my_mp%max
            m = m_values(mp)
            norm = 1.0d0/to_spec_norm(m)
            mytest%s2a(mp)%data(l_max,:,:,:) = 0.0d0
            mytest%s2a(mp)%data(:,:,:,:) = norm*mytest%s2a(mp)%data(:,:,:,:)
        enddo

        ! compute Spec->Phys error using analytic solution: LT(spectral) vs. true_physical
        my_Nr = my_r%delta
        mxdiff_spec_to_phys = -1.0d0
        do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            if (m .eq. mval) then
                diff = -1.0d0
                do r=my_r%min,my_r%max
                    do f=1,nf
                        ! index = (r-rlo+1) + (imi-1)*Nr + (f-1)*Nr*2
                        ind = (r-my_r%min+1) + (f-1)*my_Nr*2 ! real part for this r/field
                        amp = (2.0d0*f+1.0d0)*radius(r)
                        diff1 = maxval(abs(amp*true_soln(:) - mytest%p2a(:,ind,mp)))

                        ind = (r-my_r%min+1) + my_Nr + (f-1)*my_Nr*2 ! imag part
                        amp = -(3.0d0*f*f+2.0d0)*radius(r)**2
                        diff2 = maxval(abs(amp*true_soln(:) - mytest%p2a(:,ind,mp)))

                        diff = max(diff1, diff2, diff)
                    enddo
                enddo
            else
                diff = maxval(abs(mytest%p2a(:,:,mp))) ! max th/r/real/imag/field at this m
            endif
            if (diff .gt. mxdiff_spec_to_phys) mxdiff_spec_to_phys = diff
        enddo

        ! compute Spec->Phys->Spec error: LT^{-1}[ LT(spectral) ] vs. spectral
        diff = -1.0d0
        mxdiff_spec_to_spec = -1.0d0
        do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            do l=m,l_max
                if ((l .eq. lval) .and. (m .eq. mval)) then
                    do r=my_r%min,my_r%max
                        do f=1,nf
                           ans = (2.0d0*f + 1.0d0)*radius(r)
                           diff1 = abs(ans - mytest%s2a(mp)%data(l,r,1,f))

                           ans = -(3.0d0*f*f + 2.0d0)*radius(r)**2
                           diff2 = abs(ans - mytest%s2a(mp)%data(l,r,2,f))

                           diff = max(diff1, diff2)
                        enddo
                    enddo
                else
                    diff = maxval(abs(mytest%s2a(mp)%data(l,:,:,:)))
                endif
                if (diff .gt. mxdiff_spec_to_spec) mxdiff_spec_to_spec = diff
            enddo
        enddo

        ! gather max data to root and report results to stdout
        rowrank = pfi%rcomm%rank
        colrank = pfi%ccomm%rank
        Call Global_Max(mxdiff_spec_to_phys, mxdiff_spec_to_phys_reduced)
        Call Global_Max(mxdiff_spec_to_spec, mxdiff_spec_to_spec_reduced)

        if ((rowrank .eq. 0) .and. (colrank .eq. 0)) then
            write(6,*)
            write(6,*) 'Using tolerance = ', passing_tolerance
            write(6,*)
            if (mxdiff_spec_to_phys_reduced .lt. passing_tolerance) then
                write(6,*) 'Spec -> Phys result = PASSING'
            else
                write(6,*) 'Spec -> Phys result = FAIL, max error: ', &
                                                        mxdiff_spec_to_phys_reduced
            endif
            if (mxdiff_spec_to_spec_reduced .lt. passing_tolerance) then
                write(6,*) 'Spec -> Phys -> Spec result = PASSING'
            else
                write(6,*) 'Spec -> Phys -> Spec result = FAIL, max error: ', &
                                                                mxdiff_spec_to_spec_reduced
            endif
            write(6,*)
        endif

        ! cleanup
        Call mytest%deconstruct('p2a')
        Call mytest%deconstruct('s2a')

        deallocate(to_phys_norm, to_spec_norm, true_soln)

    End Subroutine single_mode_LT

    Subroutine single_mode_loops_LT()
        ! -intialize single (l,m) mode in spectral space
        ! -loop over spectral -> physical -> spectral
        ! -compute the error
        Implicit None
        Integer :: colrank, rowrank, nf=3, lval, mval
        Integer :: mp, l, m, r, f, i
        Integer :: fcount(3,2)
        Real*8 :: mxdiff, diff, diff1, diff2, ans, norm, mxdiff_reduced
        Real*8, allocatable :: to_phys_norm(:), to_spec_norm(:)
        type(SphericalBuffer) :: mytest

        if ((pfi%rcomm%rank .eq. 0) .and. (pfi%ccomm%rank .eq. 0)) then
            write(6,*) '////////////////////////////////////////////////////////////'
            write(6,*) 'Test Legendre Transform: single mode, loop over forward/inverse LT'
        endif

        fcount(:,:) = nf

        ! the FFT normalizations are built into the normalizations used in the Legendre
        ! transform. build the FFT normalizations here so they can be accounted for later
        !     these are set/defined in Legendre_Polynomials.F90 (search PTS or STP)
        allocate(to_phys_norm(0:l_max), to_spec_norm(0:l_max))
        to_phys_norm(0) = 1.0d0  ! m=0
        to_phys_norm(1:) = 0.5d0 ! m/=0
        to_spec_norm(0) = 1.0d0/(n_phi)    ! m=0
        to_spec_norm(1:) = 1.0d0/(n_theta) ! m/=0

        ! allocate fields in spectral space
        Call mytest%init(field_count=fcount, config='s2a')
        Call mytest%construct('s2a')

        ! initialize specific (l,m) values
        lval = 2; mval = 1
        do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            do l=m,l_max
                mytest%s2a(mp)%data(l,:,:,:) = 0.0d0
                if ((l .eq. lval) .and. (m .eq. mval)) then
                    do r=my_r%min,my_r%max
                        do f=1,nf
                            mytest%s2a(mp)%data(l,r,1,f) = (2.0d0*f-1.0d0)*radius(r)
                            mytest%s2a(mp)%data(l,r,2,f) = -(3.0d0*f*f-5.0d0)*radius(r)**2
                        enddo
                    enddo
                endif
            enddo
        enddo

        ! build physical space
        Call mytest%construct('p2a')

        ! compute spectral -> physical -> spectral -> physical -> ...
        do i=1, ntest_transform_loops

            do mp=my_mp%min, my_mp%max ! "undo" FFT norm of to-physical direction
               m = m_values(mp)
               norm = 1.0d0/to_phys_norm(m)
               mytest%s2a(mp)%data(:,:,:,:) = norm*mytest%s2a(mp)%data(:,:,:,:)
            enddo

            Call Legendre_Transform(mytest%s2a,mytest%p2a) ! to physical

            Call Legendre_Transform(mytest%p2a,mytest%s2a) ! back to spectral

            ! zero out l_max mode (done in production runs too), undo FFT norm of to-spectral
            do mp=my_mp%min,my_mp%max
                m = m_values(mp)
                norm = 1.0d0/to_spec_norm(m)
                mytest%s2a(mp)%data(l_max,:,:,:) = 0.0d0
                mytest%s2a(mp)%data(:,:,:,:) = norm*mytest%s2a(mp)%data(:,:,:,:)
            enddo
        enddo

        ! compute error: LT^{-1}[ LT(spectral) ] vs. initial_spectral
        diff = -1.0d0
        mxdiff = -1.0d0
        do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            do l=m,l_max
                if ((l .eq. lval) .and. (m .eq. mval)) then
                    do r=my_r%min,my_r%max
                        do f=1,nf
                           ans = (2.0d0*f - 1.0d0)*radius(r)
                           diff1 = abs(ans - mytest%s2a(mp)%data(l,r,1,f))

                           ans = -(3.0d0*f*f - 5.0d0)*radius(r)**2
                           diff2 = abs(ans - mytest%s2a(mp)%data(l,r,2,f))

                           diff = max(diff1, diff2)
                        enddo
                    enddo
                else
                    diff = maxval(abs(mytest%s2a(mp)%data(l,:,:,:)))
                endif
                if (diff .gt. mxdiff) mxdiff = diff
            enddo
        enddo

        ! gather max data to root and report results to stdout
        rowrank = pfi%rcomm%rank
        colrank = pfi%ccomm%rank
        Call Global_Max(mxdiff, mxdiff_reduced)

        if ((rowrank .eq. 0) .and. (colrank .eq. 0)) then
            write(6,*)
            write(6,*) 'Using tolerance = ', passing_tolerance
            write(6,*) 'Completed number of loops = ', ntest_transform_loops
            write(6,*)
            if ((mxdiff_reduced/ntest_transform_loops) .lt. passing_tolerance) then
                write(6,*) 'Spec -> Phys -> Spec -> ... result = PASSING'
            else
                write(6,*) 'Spec -> Phys -> Spec -> ... result = FAIL, error per loop: ', &
                                                         mxdiff_reduced/ntest_transform_loops
            endif
            write(6,*)
        endif

        ! cleanup
        Call mytest%deconstruct('p2a')
        Call mytest%deconstruct('s2a')

        deallocate(to_phys_norm, to_spec_norm)

    End Subroutine single_mode_loops_LT

    Subroutine white_noise_LT()
        ! -intialize all modes to have unit amplitude
        ! -loop over spectral -> physical -> spectral
        ! -compute the error
        Implicit None
        Integer :: i, m, mp, l, nf=3
        Integer :: fcount(3,2)
        Integer :: colrank, rowrank
        Real*8  :: mxdiff, diff, mxdiff_reduced, norm
        Real*8, allocatable :: to_phys_norm(:), to_spec_norm(:)
        type(SphericalBuffer) :: test

        if ((pfi%rcomm%rank .eq. 0) .and. (pfi%ccomm%rank .eq. 0)) then
            write(6,*) '////////////////////////////////////////////////////////////'
            write(6,*) 'Test Legendre Transform: all modes unit amplitude, loop over forward/inverse LT'
        endif

        ! the FFT normalizations are built into the normalizations used in the Legendre
        ! transform. build the FFT normalizations here so they can be accounted for later
        !     these are set/defined in Legendre_Polynomials.F90 (search PTS or STP)
        allocate(to_phys_norm(0:l_max), to_spec_norm(0:l_max))
        to_phys_norm(0) = 1.0d0  ! m=0
        to_phys_norm(1:) = 0.5d0 ! m/=0
        to_spec_norm(0) = 1.0d0/(n_phi)    ! m=0
        to_spec_norm(1:) = 1.0d0/(n_theta) ! m/=0

        fcount(:,:) = nf ! specify number of fields

        ! allocate fields in spectral space
        Call test%init(field_count=fcount, config='s2a')
        Call test%construct('s2a')

        ! initialize modes
        Do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            if (m .gt. 0) then
                Do l=m, l_max-1 ! exclude l_max mode
                    test%s2a(mp)%data(l,:,:,:) = 1.0d0 ! unit amplitude for all modes
                Enddo
            else
                ! exclude l_max mode
                test%s2a(mp)%data(:l_max-1,:,1,:) = 1.0d0 ! real part of m=0 set to unit amp
                test%s2a(mp)%data(:l_max-1,:,2,:) = 0.0d0 ! imag part of m=0 must be zero
            endif

            ! include l_max mode
            test%s2a(mp)%data(l_max,:,:,:) = 0.0d0 ! zero out the l_max mode
        Enddo

        ! build physical space
        Call test%construct('p2a')

        ! compute spectral -> physical -> spectral -> physical -> ...
        Do i=1, ntest_transform_loops

            do mp=my_mp%min, my_mp%max ! "undo" FFT norm of to-physical direction
               m = m_values(mp)
               norm = 1.0d0/to_phys_norm(m)
               test%s2a(mp)%data(:,:,:,:) = norm*test%s2a(mp)%data(:,:,:,:)
            enddo

            Call Legendre_Transform(test%s2a,test%p2a) ! to physical

            Call Legendre_Transform(test%p2a,test%s2a) ! back to spectral

            ! zero out l_max mode (done in production runs too), undo FFT norm of to-spectral
            do mp=my_mp%min,my_mp%max
                m = m_values(mp)
                norm = 1.0d0/to_spec_norm(m)
                test%s2a(mp)%data(l_max,:,:,:) = 0.0d0
                test%s2a(mp)%data(:,:,:,:) = norm*test%s2a(mp)%data(:,:,:,:)
            enddo
        Enddo

        diff = -1.0d0
        mxdiff = -1.0d0
        Do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            if (m .gt. 0) then
                Do l=m, l_max-1 ! exclude the l_max mode
                    diff = maxval(abs(1.0d0-test%s2a(mp)%data(l,:,:,:))) ! all modes = 1
                Enddo
            else
                ! exclude the l_max mode
                diff = maxval(abs(1.0d0 - test%s2a(mp)%data(:l_max-1,:,1,:)))    ! real part=1
                diff = max(diff, maxval(abs(test%s2a(mp)%data(:l_max-1,:,2,:)))) ! imag part=0
            endif

            ! include the l_max mode
            diff = max(diff, maxval(abs(test%s2a(mp)%data(l_max,:,:,:)))) ! l_max mode = 0

            if (diff .gt. mxdiff) mxdiff = diff
        Enddo

        ! gather max data to root and report results to stdout
        rowrank = pfi%rcomm%rank
        colrank = pfi%ccomm%rank
        Call Global_Max(mxdiff, mxdiff_reduced)

        if ((rowrank .eq. 0) .and. (colrank .eq. 0)) then
            write(6,*)
            write(6,*) 'Using tolerance = ', passing_tolerance
            write(6,*) 'Completed number of loops = ', ntest_transform_loops
            write(6,*)
            if (mxdiff_reduced/ntest_transform_loops .lt. passing_tolerance) then
                write(6,*) 'Spec -> Phys -> Spec -> ... result = PASSING'
            else
                write(6,*) 'Spec -> Phys -> Spec -> ... result = FAIL, error per loop: ', &
                                                         mxdiff_reduced/ntest_transform_loops
            endif
            write(6,*)
        endif

        ! cleanup
        Call test%deconstruct('p2a')
        Call test%deconstruct('s2a')

        deallocate(to_phys_norm, to_spec_norm)

    End Subroutine white_noise_LT

    Subroutine test_LT_FFT()
        ! -intialize all modes to "random" values
        ! -loop that includes LegendreTransform, FFT, and transposes, i.e., full SHT
        ! -compute the error
        Implicit None
        Integer :: i, m, mp, l, f, r
        Integer :: fcount(3,2)
        Integer :: nf = 5
        Integer :: colrank, rowrank
        Real*8  :: mxdiff, ans, diff, mxdiff_reduced
        Real*8, allocatable :: norm(:,:)
        type(SphericalBuffer) :: test

        if ((pfi%rcomm%rank .eq. 0) .and. (pfi%ccomm%rank .eq. 0)) then
            write(6,*) '////////////////////////////////////////////////////////////'
            write(6,*) 'Test Legendre Transforms, FFT, and transpose'
        endif

        fcount(:,:) = nf

        ! initialize fields in spectral space
        Call test%init(field_count=fcount, config='s2a')
        Call test%construct('s2a')

        allocate(norm(my_mp%min:my_mp%max,1:2))
        Do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            Do l=m, l_max-1 ! exclude l_max mode
                Do r=my_r%min,my_r%max
                    Do f=1, nf
                        test%s2a(mp)%data(l,r,1,f) = (l+1)*1.0d0*(2*f+3)*r   ! real part
                        test%s2a(mp)%data(l,r,2,f) = (m+1)*1.0d0*(4*f+5)*r*r ! imag part
                    Enddo
                Enddo
            Enddo
            ! normalize the entries, makes for more reliable error computation
            norm(mp,1) = maxval(abs(test%s2a(mp)%data(:,:,1,:)))
            norm(mp,2) = maxval(abs(test%s2a(mp)%data(:,:,2,:)))
            test%s2a(mp)%data(:,:,1,:) = test%s2a(mp)%data(:,:,1,:) / norm(mp,1)
            test%s2a(mp)%data(:,:,2,:) = test%s2a(mp)%data(:,:,2,:) / norm(mp,2)

            test%s2a(mp)%data(l_max,:,:,:) = 0.0d0 ! set l_max mode to zero

            if (m .eq. 0) then
                test%s2a(mp)%data(:,:,2,:) = 0.0d0 ! m=0 has no imag part
            Endif
        Enddo

        ! build physical space
        Call test%construct('p2a')

        ! model the full simulation loop
        Do i=1, ntest_transform_loops

            do mp=my_mp%min,my_mp%max ! zero out l_max mode (done in production runs too)
                test%s2a(mp)%data(l_max,:,:,:) = 0.0d0
            enddo

            Call Legendre_Transform(test%s2a,test%p2a) ! LT to physical: now in (th,r,m)
            test%config = 'p2a'
            Call test%reform()                         ! transpose
            Call fft_to_physical(test%p3a,rsc=.true.)  ! FFT to physical: now in (th,r,phi)

            Call test%construct('p3b')
            test%p3b = test%p3a
            test%config = 'p3b'
            Call test%deconstruct('p3a')
            Call fft_to_spectral(test%p3b, rsc=.true.) ! FFT to spectral: now in (th,r,m)
            Call test%reform()                         ! transpose

            Call Legendre_Transform(test%p2b,test%s2a) ! LT to spectral: now in (l,r,m)
            Call test%deconstruct('p2b')
            Call test%construct('p2a')

            do mp=my_mp%min,my_mp%max ! zero out l_max mode (done in production runs too)
                test%s2a(mp)%data(l_max,:,:,:) = 0.0d0
            enddo
        Enddo

        ! compute error
        diff = -1.0d0
        mxdiff = -1.0d0
        Do mp=my_mp%min, my_mp%max
            m = m_values(mp)
            Do l=m, l_max-1 ! exclude l_max mode
                Do r=my_r%min,my_r%max
                    Do f=1, nf
                        ans = (l+1)*1.0d0*(2*f+3)*r/norm(mp,1)   ! real part
                        diff = abs(ans - test%s2a(mp)%data(l,r,1,f))

                        if (m .gt. 0) then
                            ans = (m+1)*1.0d0*(4*f+5)*r*r/norm(mp,2) ! imag part
                        else
                            ans = 0.0d0 ! imag part
                        endif
                        diff = max(diff, abs(ans - test%s2a(mp)%data(l,r,2,f)))
                    Enddo
                Enddo
            Enddo
            diff = max(diff, maxval(abs(test%s2a(mp)%data(l_max,:,:,:)))) ! l_max mode
            if (diff .gt. mxdiff) mxdiff = diff
        Enddo

        ! gather max data to root and report results to stdout
        rowrank = pfi%rcomm%rank
        colrank = pfi%ccomm%rank
        Call Global_Max(mxdiff, mxdiff_reduced)

        if ((rowrank .eq. 0) .and. (colrank .eq. 0)) then
            write(6,*)
            write(6,*) 'Using tolerance = ', passing_tolerance
            write(6,*) 'Completed number of loops = ', ntest_transform_loops
            write(6,*)
            if ((mxdiff_reduced/ntest_transform_loops) .lt. passing_tolerance) then
                write(6,*) 'LT/FFT/transpose result = PASSING'
            else
                write(6,*) 'LT/FFT/transpose result = FAIL, error per loop: ', &
                                                         mxdiff_reduced/ntest_transform_loops
            endif
            write(6,*)
        endif

        ! cleanup
        Call test%deconstruct('p2a')
        Call test%deconstruct('s2a')

    End Subroutine test_LT_FFT

End Module Test_SHT
