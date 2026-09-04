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

Module Legendre_Polynomials
#ifdef NVHPC_COMPILER
    use Iso_Fortran_Env, Only : REAL64
    integer, parameter :: qp = REAL64 ! switch to double for nvfortran
#else
    use Iso_Fortran_Env, Only : REAL128
    integer, parameter :: qp = REAL128
#endif
    ! NOTE - need to convert everything here except for the last step to quad precision eventually

    ! This module contains all code necessary to generate and store the Legendre polynomials,
    ! their colocation points, and their associated integration weights.
    ! The polynomials computed are actually the renormalized associated legendre polynomials -
    !  - meaning that they carry the spherical harmonic normalization.
    Real(kind=qp), Allocatable :: coloc(:), gl_weights(:)
    Integer :: n_theta
    Integer :: l_max, n_m
    Integer :: m_mod = 1    ! Only calculate p_lms for every m_mod'th m
    Integer, Allocatable :: m_values(:),n_l(:),n_l_even(:),n_l_odd(:)
    Logical :: parity = .true.
    Real(kind=qp), parameter ::    PiQuad  = 3.1415926535897932384626433832795028841972_qp
    Type, Public :: even_odd_sep
        Real*8, Allocatable :: even(:)
        Real*8, Allocatable :: odd(:)
    End Type even_odd_sep
    Type, Public :: even_odd_sepi
        Integer, Allocatable :: even(:)
        Integer, Allocatable :: odd(:)
    End Type even_odd_sepi
    Type, Public :: p_lm_array
        Real*8, Allocatable :: data(:,:)
    End Type p_lm_array

    Type, Public :: p_lm_array_quad
        Real(kind=qp), Allocatable :: data(:,:)
    End Type p_lm_array_quad

    Type(p_lm_array_quad), Allocatable :: p_lmq(:)
    Type(p_lm_array), Allocatable :: p_lm(:), ip_lm(:)
    Type(p_lm_array), Allocatable :: p_lm_odd(:), p_lm_even(:)
    Type(p_lm_array), Allocatable :: ip_lm_odd(:), ip_lm_even(:) ! i means 'integration weights included'
    ! Spin-weighted (s = +/-1) tables for the compressible horizontal pair
    ! (Vasil et al. 2019 eq. 26-28 construction; same 1/(2pi) dx-normalization
    ! family as p_lm, so the 2pi*gl_weights analysis convention is shared).
    Logical :: build_spin_tables = .false.
    Type(p_lm_array), Allocatable :: wp_lm(:), wm_lm(:), iwp_lm(:), iwm_lm(:)
    ! Theta-derivative synthesis tables d(W+-)/dtheta for the advection supply
    Type(p_lm_array), Allocatable :: wdp_lm(:), wdm_lm(:)
    ! s=0 (scalar) V19 tables at FULL theta, parity-unsorted: the slot<->spin
    ! converters need full-theta scalar synthesis/analysis, which p_lm/ip_lm
    ! cannot provide under use_parity (those are half-theta, resorted).
    Type(p_lm_array), Allocatable :: ws0_lm(:), iws0_lm(:)
    Type(even_odd_sep), Allocatable :: lvals(:)
    Type(even_odd_sepi), Allocatable :: lvalsi(:)
Contains

Subroutine Compute_Spin_Wlms()
    ! Build the s=+1 (wp) and s=-1 (wm) theta-tables from the V19 eq-(26)
    ! Jacobi construction, quad precision, on the same coloc/gl_weights grid.
    ! Certified against the standalone reference (build_spin_tables.f90) at 1.6e-13.
    ! NOTE: use_parity resorting is NOT applied to the spin tables (parity maps
    ! s -> -s and mixes the streams); the spin transform pass runs unparitied.
    Implicit None
    Real(kind=qp) :: shq, chq, nrm, fr, a2b2, c1q, c2q, c3q, xq
    Real(kind=qp) :: pts_facq, stp_facq, wrawq, wdrawq
    Real(kind=qp), Allocatable :: jac0(:), jac1(:), jacn(:), jacd(:)
    Integer :: i,l,m,mv,s,is,a,b,l0,l1,n,k

    Allocate(wp_lm(1:n_m), wm_lm(1:n_m), iwp_lm(1:n_m), iwm_lm(1:n_m))
    Allocate(wdp_lm(1:n_m), wdm_lm(1:n_m))
    Allocate(ws0_lm(1:n_m), iws0_lm(1:n_m))
    Allocate(jac0(1:n_theta), jac1(1:n_theta), jacn(1:n_theta), jacd(1:n_theta))
    Do m = 1, n_m
        mv = m_values(m)
        ! Match the scalar tables' FFT-normalization convention (see the
        ! PTS_normalization/STP_normalization wrap for ip_lm/p_lm): the
        ! synthesis tables carry STP (1 at m=0, 1/2 otherwise) and the
        ! analysis tables carry PTS (1/(2*n_theta) at m=0, 1/n_theta otherwise).
        If (mv .eq. 0) Then
            pts_facq = 1.0_qp/(2.0_qp*n_theta)
            stp_facq = 1.0_qp
        Else
            pts_facq = 1.0_qp/real(n_theta,qp)
            stp_facq = 0.5_qp
        Endif
        Allocate( wp_lm(m)%data(mv:l_max,1:n_theta),  wm_lm(m)%data(mv:l_max,1:n_theta))
        Allocate(iwp_lm(m)%data(1:n_theta,mv:l_max), iwm_lm(m)%data(1:n_theta,mv:l_max))
        Allocate(wdp_lm(m)%data(mv:l_max,1:n_theta), wdm_lm(m)%data(mv:l_max,1:n_theta))
        Allocate(ws0_lm(m)%data(mv:l_max,1:n_theta), iws0_lm(m)%data(1:n_theta,mv:l_max))
        wp_lm(m)%data = 0.0d0;  wm_lm(m)%data = 0.0d0
        iwp_lm(m)%data = 0.0d0; iwm_lm(m)%data = 0.0d0
        wdp_lm(m)%data = 0.0d0; wdm_lm(m)%data = 0.0d0
        ws0_lm(m)%data = 0.0d0; iws0_lm(m)%data = 0.0d0
        Do is = 1, 3
            If (is .le. 2) Then
                s = 3-2*is
            Else
                s = 0
            Endif
            a = abs(mv+s); b = abs(mv-s); l0 = max(abs(mv),abs(s)); l1 = min(abs(mv),abs(s))
            Do l = max(l0,mv), l_max
                n = l-l0
                jac0(:) = 1.0_qp
                If (n .ge. 1) Then
                    Do i = 1, n_theta
                        jac1(i) = 0.5_qp*(a-b)+0.5_qp*(a+b+2.0_qp)*coloc(i)
                    Enddo
                Endif
                If (n .eq. 0) Then
                    jacn = jac0
                Else If (n .eq. 1) Then
                    jacn = jac1
                Else
                    Do k = 2, n
                        c1q = 2.0_qp*k*(k+a+b)*(2.0_qp*k+a+b-2.0_qp)
                        a2b2 = real(a*a-b*b,qp)
                        Do i = 1, n_theta
                            c2q = (2.0_qp*k+a+b-1.0_qp)*((2.0_qp*k+a+b)*(2.0_qp*k+a+b-2.0_qp)*coloc(i)+a2b2)
                            c3q = 2.0_qp*(k+a-1.0_qp)*(k+b-1.0_qp)*(2.0_qp*k+a+b)
                            jacn(i) = (c2q*jac1(i)-c3q*jac0(i))/c1q
                        Enddo
                        jac0 = jac1; jac1 = jacn
                    Enddo
                Endif
                fr = 1.0_qp
                Do k = l+l1+1, l+l0
                    fr = fr*real(k,qp)
                Enddo
                Do k = l-l0+1, l-l1
                    fr = fr/real(k,qp)
                Enddo
                nrm = sqrt((2.0_qp*l+1.0_qp)*fr/(4.0_qp*PiQuad))
                If (mod(max(mv,-s),2) .eq. 1) nrm = -nrm
                ! Jacobi derivative for the theta-derivative table:
                ! J'(x) = (n+a+b+1)/2 * P^(a+1,b+1)_{n-1}(x)
                If (n .ge. 1) Then
                    jac0(:) = 1.0_qp
                    If (n-1 .ge. 1) Then
                        Do i = 1, n_theta
                            jac1(i) = 0.5_qp*(a-b)+0.5_qp*(a+b+4.0_qp)*coloc(i)
                        Enddo
                    Endif
                    If (n-1 .eq. 0) Then
                        jacd = jac0
                    Else If (n-1 .eq. 1) Then
                        jacd = jac1
                    Else
                        Do k = 2, n-1
                            c1q = 2.0_qp*k*(k+a+b+2.0_qp)*(2.0_qp*k+a+b)
                            a2b2 = real((a-b),qp)*real((a+b+2),qp)
                            Do i = 1, n_theta
                                c2q = (2.0_qp*k+a+b+1.0_qp)*((2.0_qp*k+a+b+2.0_qp)*(2.0_qp*k+a+b)*coloc(i)+a2b2)
                                c3q = 2.0_qp*(k+a)*(k+b)*(2.0_qp*k+a+b+2.0_qp)
                                jacd(i) = (c2q*jac1(i)-c3q*jac0(i))/c1q
                            Enddo
                            jac0 = jac1; jac1 = jacd
                        Enddo
                    Endif
                    jacd = jacd*0.5_qp*(n+a+b+1.0_qp)
                Else
                    jacd(:) = 0.0_qp
                Endif
                Do i = 1, n_theta
                    xq  = coloc(i)
                    shq = sqrt((1.0_qp-xq)/2.0_qp)
                    chq = sqrt((1.0_qp+xq)/2.0_qp)
                    wrawq  = nrm*shq**a*chq**b*jacn(i)
                    wdrawq = nrm*( 0.5_qp*a*shq**max(a-1,0)*chq**(b+1)*jacn(i) &
                            - 0.5_qp*b*shq**(a+1)*chq**max(b-1,0)*jacn(i) &
                            - 2.0_qp*shq**(a+1)*chq**(b+1)*jacd(i) )
                    If (s .eq. 1) Then
                        wp_lm(m)%data(l,i)  = wrawq*stp_facq
                        iwp_lm(m)%data(i,l) = wrawq*2.0_qp*PiQuad*gl_weights(i)*pts_facq
                        wdp_lm(m)%data(l,i) = wdrawq*stp_facq
                    Else If (s .eq. -1) Then
                        wm_lm(m)%data(l,i)  = wrawq*stp_facq
                        iwm_lm(m)%data(i,l) = wrawq*2.0_qp*PiQuad*gl_weights(i)*pts_facq
                        wdm_lm(m)%data(l,i) = wdrawq*stp_facq
                    Else
                        ws0_lm(m)%data(l,i)  = wrawq*stp_facq
                        iws0_lm(m)%data(i,l) = wrawq*2.0_qp*PiQuad*gl_weights(i)*pts_facq
                    Endif
                Enddo
            Enddo
        Enddo
    Enddo
    DeAllocate(jac0,jac1,jacn,jacd)
End Subroutine Compute_Spin_Wlms

Subroutine Finalize_Legendre()
    Implicit None
    DeAllocate(coloc, gl_weights)
    DeAllocate(m_values,n_l)
    Call DeAllocate_Plms(depar = .true.)
End Subroutine Finalize_Legendre

Subroutine DeAllocate_Plms(depar)
    Implicit None
    Integer :: i
    Logical, Optional, Intent(In) :: depar
    If (allocated(p_lm)) Then
        Do i = 1, n_m
            If (allocated(p_lm(i)%data)) Then
                DeAllocate(p_lm(i)%data)
            Endif
        Enddo
        DeAllocate(p_lm)
    Endif
    If (present(depar)) Then
        If (depar .and. parity) Then
            Call DeAllocate_Parity_Plms()
        Endif
    Endif
End Subroutine DeAllocate_Plms

Subroutine DeAllocate_Parity_Plms()
    Implicit None
    Integer :: i,m
    If (allocated(p_lm_odd)) Then
        Do i = 1, n_m
            If (allocated(p_lm_odd(i)%data)) Then
                DeAllocate(p_lm_odd(i)%data)
            Endif
        Enddo
        DeAllocate(p_lm_odd)
    Endif
    If (allocated(p_lm_even)) Then
        Do i = 1, n_m
            If (allocated(p_lm_even(i)%data)) Then
                DeAllocate(p_lm_even(i)%data)
            Endif
        Enddo
        DeAllocate(p_lm_even)
    Endif
    If (allocated(n_l_even)) DeAllocate(n_l_even)
    If (allocated(n_l_odd)) DeAllocate(n_l_odd)
    If (allocated(lvals)) Then
        Do m = 1, n_m
            If (allocated(lvals(m)%even)) DeAllocate(lvals(m)%even)
            If (allocated(lvals(m)%odd)) DeAllocate(lvals(m)%odd)
        Enddo
        DeAllocate(lvals)
    Endif
    If (allocated(lvalsi)) Then
        Do m = 1, n_m
            If (allocated(lvalsi(m)%even)) DeAllocate(lvalsi(m)%even)
            If (allocated(lvalsi(m)%odd)) DeAllocate(lvalsi(m)%odd)
        Enddo
        DeAllocate(lvalsi)
    Endif

End Subroutine DeAllocate_Parity_Plms

Subroutine Initialize_Legendre(nt,lmax,mval,parity_in,spin_tables_in)
    Implicit None
    Real(kind=qp) :: coloc_min,coloc_max
    Logical, Intent(In) :: parity_in
    Logical, Intent(In), Optional :: spin_tables_in
    Integer, Intent(in) :: nt,lmax, mval(:)
    If (present(spin_tables_in)) build_spin_tables = spin_tables_in


    parity = parity_in
    ! Set up the grid
    coloc_min = -1
    coloc_max = 1
    n_theta = nt
    l_max = lmax
    Allocate(coloc(1:n_theta))
    Allocate(gl_weights(1:n_theta))
    n_m = size(mval)
    Allocate(m_values(1:n_m))
    m_values(:) = mval(:)
    Call Find_Colocation(coloc_min, coloc_max,coloc,gl_weights,n_theta)
    Call Compute_Plms()
    If (build_spin_tables) Call Compute_Spin_Wlms()
End Subroutine Initialize_Legendre


Subroutine Compute_Plms()
    ! Subroutine Compute_Plms(m_values,n_theta, l_max)
    ! We feed this a list of m_values (presumably distributed across processors)
    ! And also l_max.  This is sufficient to initialize the legendre polynomials
    Implicit None
    Real(kind=qp) ::  x,tmp,factorial_ratio,amp, renorm
    Integer :: i, m, l, mv, ntmax

    n_m = size(m_values)

    Allocate(n_l(1:n_m))
    ! Now calculate the p_lms
    ! y_lm(theta) is stored as ylm(m)%(theta,el)
    ! Forget about de-aliasing right now

    Allocate(p_lm(1:n_m))
    Allocate(p_lmq(1:n_m))
    Allocate(ip_lm(1:n_m))

    ntmax = n_theta
    If (parity) Then
        Allocate(p_lm_odd(1:n_m))
        Allocate(p_lm_even(1:n_m))
        Allocate(ip_lm_odd(1:n_m))
        Allocate(ip_lm_even(1:n_m))
        Allocate(n_l_even(1:n_m))
        Allocate(n_l_odd(1:n_m))
        Allocate(lvals(1:n_m))
        Allocate(lvalsi(1:n_m))
        ntmax = n_theta/2
    Endif


    ! Compute P_lm(theta) for all l's at each m
    ! Calculation is done in quad precision.  Storage is done in double.
    ! One m at a time to save memory

    Do m = 1, n_m
        n_l(m) = l_max-m_values(m)+1
        Allocate(p_lmq(m)%data(1:ntmax,m_values(m):l_max))




    ! First, fill in the l = m pieces (closed form expression)
    ! and the l = m+1 pieces
    factorial_ratio = 1.0_qp

        mv = m_values(m)
        Call compute_factorial_ratio(mv,factorial_ratio)
        amp = ((mv+0.5_qp)/(2.0_qp*PiQuad))**0.5_qp
        amp = amp*factorial_ratio
        Do i = 1, ntmax
            x = coloc(i)
            tmp = 1.0_qp-x*x

            If (mod(mv,2) .eq. 1) Then
                !odd m
                p_lmq(m)%data(i,mv) = -amp*tmp**(mv/2+0.5_qp)
            Else
                !even m
                p_lmq(m)%data(i,mv) = amp*tmp**(mv/2)
            Endif
            If (mv .lt. l_max) then
                p_lmq(m)%data(i,mv+1) = p_lmq(m)%data(i,mv)*x*(2.0_qp*mv+3)**0.5_qp
            Endif
        Enddo


        !General recursion for l > m+1
        mv = m_values(m)
        Do l = mv+2, l_max
            Do i = 1, ntmax
                x = coloc(i)
                amp = (l-1)**2-mv*mv
                amp = amp/ (4.0_qp*(l-1)**2-1.0_qp)
                amp = amp**0.5_qp
                tmp = p_lmq(m)%data(i,l-1)*x-amp*p_lmq(m)%data(i,l-2)
                amp = (4.0_qp*l*l-1.0_qp)/(l*l-mv*mv)
                p_lmq(m)%data(i,l) = tmp*amp**0.5_qp
            Enddo
        Enddo

        If (parity) Then
            Call parity_resort(m)
        Else
            !Fill in double precision arrays
            ! Add normalization for integration
            Allocate(ip_lm(m)%data(1:ntmax,m_values(m):l_max))
            Allocate(p_lm(m)%data(m_values(m):l_max,1:ntmax))
            mv = m_values(m)
            Do l = mv, l_max
                Do i = 1, ntmax

                    p_lm(m)%data(l,i)  = p_lmq(m)%data(i,l)
                    renorm = 2.0_qp*PiQuad*gl_weights(i)
                    tmp = p_lmq(m)%data(i,l)*renorm
                    ip_lm(m)%data(i,l) = tmp
                Enddo
            Enddo
        Endif
        DeAllocate(p_lmq(m)%data)

    Enddo
    ! DeAllocate data structures that will no longer be used.
    ! "data" attributes have already been deallocated.
    DeAllocate(p_lmq)
    If (parity) Then
        DeAllocate(p_lm)
        DeAllocate(ip_lm)
    Endif

End Subroutine Compute_Plms

Subroutine Parity_Resort(m)
    Implicit None
    Integer, Intent(In) :: m
    Integer :: l, indeven, indodd,partest, i
    Real(kind=qp) :: renorm, tmp
    Real(kind=qp) :: PTS_normalization, STP_normalization
    ! Resort the p_lms into even and odd arrays

        ! We wrap a normalization factor, related to the FFT
        ! into the Legendre weights.
        If (m_values(m) .eq. 0) Then
            PTS_normalization = 1.0d0/(2.0d0*n_theta)
            STP_normalization = 1.0d0
        Else
            PTS_normalization = 1.0d0/(n_theta)
            STP_normalization = 0.5d0
        Endif
        n_l_even(m) = 0
        n_l_odd(m) = 0
        Do l = m_values(m), l_max
            partest = l-m_values(m)
            If (Mod(partest,2) .eq. 1) Then
                n_l_odd(m) = n_l_odd(m)+1
            Else
                n_l_even(m) = n_l_even(m)+1
            Endif
        Enddo


        If (n_l_even(m) .gt. 0) Then
            Allocate(ip_lm_even(m)%data(1:n_theta/2,1:n_l_even(m)))
            Allocate(p_lm_even(m)%data(1:n_l_even(m),1:n_theta/2))
            Allocate(lvals(m)%even(1:n_l_even(m)))
            Allocate(lvalsi(m)%even(1:n_l_even(m)))
        Endif
        If (n_l_odd(m) .gt. 0) Then
            Allocate(ip_lm_odd(m)%data(1:n_theta/2,1:n_l_odd(m)))
            Allocate(p_lm_odd(m)%data(1:n_l_odd(m),1:n_theta/2))
            Allocate(lvals(m)%odd(1:n_l_odd(m)))
            Allocate(lvalsi(m)%odd(1:n_l_odd(m)))
        Endif
        indeven = 1
        indodd = 1
        Do l = m_values(m), l_max
            partest = l-m_values(m)
            If (Mod(partest,2) .eq. 1) Then
                 lvals(m)%odd(indodd) = l
                lvalsi(m)%odd(indodd) = l
                Do i = 1, n_theta/2
                    renorm = 2.0_qp*PiQuad*gl_weights(i)
                    tmp = p_lmq(m)%data(i,l)*renorm
                    ip_lm_odd(m)%data(i,indodd) = tmp*PTS_normalization
                    p_lm_odd(m)%data(indodd,i) = p_lmq(m)%data(i,l)*STP_normalization
                Enddo
                indodd = indodd +1

            Else
                 lvals(m)%even(indeven) = l
                lvalsi(m)%even(indeven) = l
                Do i = 1, n_theta/2
                    renorm = 2.0_qp*PiQuad*gl_weights(i)
                    tmp = p_lmq(m)%data(i,l)*renorm
                    ip_lm_even(m)%data(i,indeven) = tmp*PTS_normalization
                    p_lm_even(m)%data(indeven,i) =   p_lmq(m)%data(i,l)*STP_normalization
                Enddo
                indeven = indeven+1
            Endif
        Enddo



End Subroutine Parity_Resort
Subroutine Find_Colocation(x1,x2,abscissas, weights, order_n)
    Implicit None
    ! Compute the Gauss-Legendre colocation points and discretized weights
    ! appropriate for the interval x1 < x < x2
    ! Based heavily on Numerical Recipes Volume 2
    ! Variables have been renamed for clarity
    ! Legendre polynomial calculation is accomplished in a separate subroutine
    !    to help with readability.
    Real(kind=qp), Intent(Out) :: abscissas(1:), weights(1:)
    Real(kind=qp), Intent(In) :: x1,x2
    Integer, Intent(In) :: order_n
    Real(kind=qp) :: pn, ith_root,deriv_pn, new_guess
    Real(kind=qp) :: midpoint, scaling
    Integer :: i, n_roots
    Logical :: converged
    Real(kind=qp)  :: eps, del
    midpoint = 0.5_qp*(x1+x2)
    scaling  = 0.5_qp*(x2-x1)
    n_roots  = (order_n+1)/2    ! Roots are symmetric - exploit symmetry

    eps = 3.0e-14_qp


    Do i = 1, n_roots
        ith_root = cos(PiQuad*(i-0.25_qp)/(order_n+0.5_qp))
        converged = .false.
        Do While (.not. converged)
            Call nth_legendre(ith_root,order_n,pn,deriv_pn)

            new_guess = ith_root - pn/deriv_pn
            del = abs(ith_root-new_guess)
            ith_root = new_guess
            if (del .le. eps) Then
                converged = .true.
            Endif
            abscissas(i) = midpoint-scaling*ith_root
            abscissas(order_n+1-i) = midpoint+scaling*ith_root
            weights(i) = 2.0_qp*scaling/((1.0_qp-ith_root*ith_root)*deriv_pn*deriv_pn)
            weights(order_n+1-i) = weights(i)
        Enddo
    Enddo
End Subroutine Find_Colocation

Subroutine nth_legendre(x,n,pn,deriv_pn)
    ! Calculates the value of the nth legendre polynomial at location x
    ! Returns x, p_n(x) and d_by_dx (p_n(x))
    Implicit None
    Real(kind=qp), Intent(Out) :: pn, deriv_pn
    Real(kind=qp), Intent(In) :: x
    Integer, Intent(In) :: n
    Integer :: j
    Real(kind=qp) :: pn_minus1, pn_minus2
    pn = 1.0_qp    !p0
    pn_minus1 = 0.0_qp
    ! Use recursion relation to build p_order_n(x)
    Do j = 1, n
        pn_minus2 = pn_minus1
        pn_minus1 = pn
        pn = ( (2.0_qp*j-1.0_qp)*x*pn_minus1 - (j-1.0_qp)*pn_minus2 )/j
    Enddo
    deriv_pn = n*(x*pn-pn_minus1)/(x*x-1.0_qp)
End Subroutine nth_legendre


Subroutine factorial(n,f)
    Integer :: n
    Real(qp), Intent(Out) :: f
    Integer :: i
    f = 1.0_qp
    Do i = 1, n
        f = f*i
    Enddo
End Subroutine factorial

Subroutine compute_factorial_ratio(m,ratio)
        ! build sqrt( (2m-1)!!(2m-1)!!/(2m)!! ) stably
        Integer, Intent(In) :: m
        Real(kind=qp), Intent(Out) :: ratio
        Integer :: i
        ratio = 1.0_qp
        Do    i = 1, m
            ratio = ratio*((i-0.5_qp)/i)**0.5_qp  !ratio = ratio*(2m-1)/(2m)
        Enddo
End Subroutine compute_factorial_ratio
End Module Legendre_Polynomials

