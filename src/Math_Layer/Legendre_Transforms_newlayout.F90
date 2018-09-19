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

! This module contains routines for performing the Legendre transforms.
! These are really just matrix-matrix multiplies that we offload to dgemm.
! The polynomials themselves are computed in Legendre_Polynomials.F90 at startup.
! Parity is accounted for, as well as "Legendre" data structures.
! Notes:
!            (Dec 5, 2013 - 9:15 p.m.) : Verified that this module sees the correct lower bounds of
!                        data_in/out(i)%data, which is allocated data(m_value(i):l_max) in external modules.
!                        Was concerned that it might index array data starting at 1, not m_value(i) [Fortran standard?]

Module Legendre_Transforms
    Use Legendre_Polynomials
    Use Timing
    Use Structures
    !Type, Public :: rmcontainer
    !    Real*8, Allocatable :: data(:,:)
    !End Type rmcontainer
    Interface Legendre_Transform
        Module Procedure PtS_2d_dgpv2, StP_2d_dgp
        Module Procedure PtS_4d_dgpv2 !, StP_4d_dgp
    End Interface
Contains

Subroutine Test_Legendre
    Implicit None
    Real*8, Allocatable :: theta(:),tmp(:)
    Real*8, Allocatable :: ylm(:,:)
    Integer :: i,l,m, nrhs, nl
    Integer, Allocatable :: l_test(:), m_test(:)
    Type(p_lm_array), Allocatable :: ans(:), ans2(:)
    Real*8 :: st, ct, tmp2,chk, alpha, beta
    Logical :: use_dgemm = .true.
    ! Test the physical to spectral legendre transform
    !  using tabulated associated legendre transforms
    nrhs = 6
    Allocate(theta(1:n_theta))
    Allocate(tmp(1:n_theta))
    Allocate(ylm(1:n_theta,1:nrhs)) ! Test the first six spherical harmonics (modulu the exp(i m phi) piece)
    Allocate(ans(1:n_m))
    Allocate(ans2(1:n_m))
    Allocate(l_test(1:nrhs))
    Allocate(m_test(1:nrhs))

    Do m = 1,n_m
        Allocate( ans(m)%data(m_values(m):l_max,nrhs))
        Allocate(ans2(m)%data(m_values(m):l_max,nrhs))
         ans(m)%data(:,:) = 0.0d0
        ans2(m)%data(:,:) = 53.0d0
    Enddo

    Do i = 1, n_theta
        theta(i) = acos(coloc(i))
        ct = cos(theta(i))
        st = sin(theta(i))
        ylm(i,1) = (0.25d0/pi)**0.5d0    !00
        ylm(i,2) = -st*(3.0d0/8.0d0/pi)**0.5d0 !11
        ylm(i,3) = ct*(3.0d0/4.0/pi)**0.5d0 ! 10
        ylm(i,4) = st*st*0.25d0*(15.0d0/2.0d0/pi)**0.5d0 !22
        ylm(i,5) = -st*ct*(15.0d0/8.0d0/pi)**0.5d0 !21
        ylm(i,6) = (1.5d0*ct*ct-0.5d0)*(5.0d0/4.0/pi)**0.5d0    !20
        tmp(i) = 1.0d0/(1-coloc(i)*coloc(i))
    Enddo
    l_test(1) = 0
    l_test(2:3) = 1
    l_test(4:6) = 2
    m_test(1) = 0
    m_test(2) = 1
    m_test(3) = 0
    m_test(4) = 2
    m_test(5) = 1
    m_test(6) = 0

    !Next perform the transform on our test harmonics

    If (parity) Then
        write(6,*)'testing parity transform'
        Call Pts_2D_drp(ylm,ans,nrhs)
        Call PtS_2D_dgp(ylm, ans2,nrhs)
    Else
        Call Pts_2D_dr(ylm, ans,nrhs)
        Call PtS_2D_dg(ylm, ans2,nrhs)
    Endif


    !Verify the orthogonality of m1 = m2, but l1 ne l2
    Do i = 1, nrhs
        Write(6,*)'============================='
        Write(6,*)'     Test mode  '
        Write(6,*)'l = ', l_test(i), 'm = ', m_test(i)
        Do m = 1, n_m
            if (m_test(i) .eq. m_values(m)) then
                Do l = m_values(m), l_max
                    Write(6,*)l, m_values(m), ans(m)%data(l,i), ans2(m)%data(l,i)
                Enddo
            Endif
        Enddo
    Enddo

    ! Next verify the orthogonality of m1 ne m2, but l1 eq l2
    !Next perform the transform on our test harmonics
    Do i = 1, nrhs
        ylm(:,i) = ylm(:,i)*tmp(:)
    Enddo


    If (parity) Then
        Call Pts_2D_drp(ylm, ans,nrhs)
        Call PtS_2D_dgp(ylm, ans2,nrhs)
    Else
        Call Pts_2D_dr(ylm, ans,nrhs)
        Call PtS_2D_dg(ylm, ans2,nrhs)
    Endif


    Do i = 1, nrhs
        Write(6,*)'============================='
        Write(6,*)'     Test 2   '
        Write(6,*)'l = ', l_test(i), 'm = ', m_test(i)
        Do m = 1, n_m
            Do l = m_values(m), l_max
                if (l_test(i) .eq. l) then
                    Write(6,*)l, m_values(m), ans(m)%data(l,i), ans2(m)%data(l,i)
                    if(m_test(i) .eq. m_values(m)) then
                        chk = (2.0d0*l+1.0d0)/2.0d0 !4.0d0/pi
                        chk = chk/m_values(m)
                        write(6,*)'check: ', chk
                    endif
                Endif
            Enddo
        Enddo
    Enddo



103 format(d7.2)

    DeAllocate(theta)
    DeAllocate(ylm)
    DeAllocate(tmp)
    Do m = 1,n_m
        DeAllocate(ans(m)%data)
        DeAllocate(ans2(m)%data)
    Enddo
    DeAllocate(ans,ans2)
    DeAllocate(l_test,m_test)
End Subroutine Test_Legendre

Subroutine Test_Simple_Dgemm
    Real*8, Allocatable :: a(:,:), b(:,:), c(:,:)
    Integer :: m, n, k, i, ntimes
    Real*8 :: alpha

    ! Just to make sure I've done dgemm correctly
    ! a = | 1 2 |  b = |5|  a#b = c = | 17 |
    !     = | 3 4 |      |6|              = | 39 |

    m = 2
    n = 2
    k = 2
    Allocate(a(1:2,1:2))
    Allocate(b(1:2,1:2))
    Allocate(c(1:2,1:2))

    a(1,1) = 1.0d0
    a(2,1) = 2.0d0
    a(1,2) = 3.0d0
    a(2,2) = 4.0d0
    b(1,1) = 5.0d0
    b(2,1) = 6.0d0
    b(1,2) = 7.0d0
    b(2,2) = 8.0d0
    c(:,:) = 0.0d0
    write(6,*)'calling dgemm'
    alpha = 1.0d0
    beta = 0.0d0
    ntimes = 100
    do i = 1, ntimes
        CALL DGEMM('T','N',m,n,k, alpha, a, m, b, k, beta,c,m)
    enddo
    Write(6,*)'c is : ', c

    a(1,1) = 1.0d0
    a(1,2) = 2.0d0
    a(2,1) = 3.0d0
    a(2,2) = 4.0d0
    CALL DGEMM('N','N',m,n,k, alpha, a, m, b, k, beta,c,m)
    Write(6,*)'c is : ', c

    b(1,1) = 5.0d0
    b(1,2) = 6.0d0
    b(2,1) = 7.0d0
    b(2,2) = 8.0d0
    CALL DGEMM('N','T',m,n,k, alpha, a, m, b, k, beta,c,m)
    Write(6,*)'c is : ', c

    a(1,1) = 1.0d0
    a(2,1) = 2.0d0
    a(1,2) = 3.0d0
    a(2,2) = 4.0d0
    CALL DGEMM('T','T',m,n,k, alpha, a, m, b, k, beta,c,m)
    Write(6,*)'c is : ', c

    DeAllocate(a,b,c)



End Subroutine Test_simple_dgemm


Subroutine Test_Simple_Dgemm2
    Real*8, Allocatable :: a(:,:), b(:,:), c(:,:)
    Integer :: m, n, k, i, ntimes
    Real*8 :: alpha

    !//////////////////////////////////////////////
    ! Just to make sure I've done dgemm correctly
    !
    ! a = | 1 2 |  b = |7  9|  a#b = c = | 23  29 |
    !       | 3 4 |      |8 10|                | 53  67 |
    !        | 5 6 |                        | 83 105 |


    !//////////////////////////////////
    !  Test 1.  Columns run fastest.
    m = 3
    n = 2
    k = 2
    Allocate(a(1:m,1:k))

    Allocate(b(1:k,1:n))

    Allocate(c(1:m,1:n))



    a(1,1) = 1.0d0
    a(1,2) = 2.0d0
    a(2,1) = 3.0d0
    a(2,2) = 4.0d0
    a(3,1) = 5.0d0
    a(3,2) = 6.0d0




    b(1,1) = 7.0d0
    b(1,2) = 9.0d0
    b(2,1) = 8.0d0
    b(2,2) = 10.0d0

    c(:,:) = 0.0d0
    write(6,*)'calling dgemm'
    alpha = 1.0d0
    beta = 0.0d0
    ntimes = 1
    do i = 1, ntimes
        CALL DGEMM('N','N',m,n,k, alpha, a, m, b, k, beta,c,m)
    enddo
    Write(6,*)'c is : ', c



    DeAllocate(a,b,c)

    !//////////////////////////////////
    !  Test 2.  Rows of A run fastest.
    !  B and C remain the same
    m = 3
    n = 2
    k = 2
    Allocate(a(1:k,1:m))

    Allocate(b(1:k,1:n))

    Allocate(c(1:m,1:n))



    a(1,1) = 1.0d0
    a(2,1) = 2.0d0
    a(1,2) = 3.0d0
    a(2,2) = 4.0d0
    a(1,3) = 5.0d0
    a(2,3) = 6.0d0




    b(1,1) = 7.0d0
    b(1,2) = 9.0d0
    b(2,1) = 8.0d0
    b(2,2) = 10.0d0

    c(:,:) = 0.0d0
    write(6,*)'calling dgemm'
    alpha = 1.0d0
    beta = 0.0d0
    ntimes = 1
    do i = 1, ntimes
        CALL DGEMM('T','N',m,n,k, alpha, a, k, b, k, beta,c,m)
    enddo
    Write(6,*)'c is : ', c



    DeAllocate(a,b,c)



End Subroutine Test_simple_dgemm2



Subroutine Time_MatMult(sz,ntimes)
    Implicit None
    Real*8, Allocatable :: a(:,:), b(:,:), c(:,:)
    Integer, Intent(in) :: sz, ntimes
    Integer :: m, n, k, i,j
    Real*8 :: alpha, beta, t1, delta


    m = sz
    n = sz
    k = sz
    Allocate(a(1:sz,1:sz))
    Allocate(b(1:sz,1:sz))
    Allocate(c(1:sz,1:sz))

    Do i = 1, sz
        Do j = 1, sz
            a(i,j) = i+j
            b(i,j) = i-j
        Enddo
    Enddo
    c(:,:) = 0.0d0

    alpha = 1.0d0
    beta = 0.0d0
    Write(6,*)'Ntimes = ', ntimes
    t1=secnds(0.0)
    do i = 1, ntimes
        CALL DGEMM('T','N',m,n,k, alpha, a, m, b, k, beta,c,m)
    enddo
    delta = secnds(t1)
    Write(6,*)'T N time: ', delta

    t1=secnds(0.0)
    do i = 1, ntimes
        CALL DGEMM('T','T',m,n,k, alpha, a, m, b, k, beta,c,m)
    enddo
    delta = secnds(t1)
    Write(6,*)'T T time: ', delta

    t1=secnds(0.0)
    do i = 1, ntimes
        CALL DGEMM('N','T',m,n,k, alpha, a, m, b, k, beta,c,m)
    enddo
    delta = secnds(t1)
    Write(6,*)'N T time: ', delta

    t1=secnds(0.0)
    do i = 1, ntimes
        CALL DGEMM('N','N',m,n,k, alpha, a, m, b, k, beta,c,m)
    enddo
    delta = secnds(t1)
    Write(6,*)'N N time: ', delta


    DeAllocate(a,b,c)



End Subroutine Time_MatMult




Subroutine PtS_2d_dg(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...DGEMM
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Real*8 :: alpha, beta
    Integer :: m,nl

    alpha = 1.0d0
    beta = 0.0d0
    !Do m = 0, l_max,m_mod
    Do m = 1, n_m
            nl = l_max-m_values(m)+1
            CALL DGEMM('T','N',nl,nrhs,n_theta, alpha, p_lm(m)%data, n_theta,data_in , n_theta, beta,data_out(m)%data,nl)
            !CALL DGEMM('T','N',nl(m),nrhs,n_theta, alpha, p_lm(m)%data, n_theta,data_in , n_theta, beta,data_out(m)%data,nl(m))
    Enddo

End Subroutine PtS_2d_dg

Subroutine PtS_2d_dr(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...Direct
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Integer :: m,i,l
    ! Need to try a version 2 of this where the sum function is not invoked
    !Do m = 0, l_max,m_mod
    Do m = 1, n_m
        Do i = 1, nrhs
            Do l = m_values(m), l_max
                data_out(m)%data(l,i) = Sum(data_in(:,i)*p_lm(m)%data(:,l))
            Enddo
        Enddo
    Enddo

End Subroutine PtS_2d_dr

Subroutine PtS_2d_drp(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...Direct...Parity
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Real*8, Allocatable :: feven(:,:), fodd(:,:)
    Integer :: m,i,j,l,nt1,offset

    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs))
    Allocate(fodd(1:n_theta/2,1:nrhs))
    feven(1:n_theta/2,:) = data_in(1:n_theta/2,:)
    fodd(1:n_theta/2,:) = data_in(1:n_theta/2,:)
    nt1 = n_theta+1
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i) = feven(j,i)+data_in(nt1-j,i)
             fodd(j,i) =  fodd(j,i)-data_in(nt1-j,i)
        Enddo
    enddo

    ! One possible way - involves a lot of indexing.
    Do m = 1, n_m
        Do i = 1, nrhs
            Do j = 1, n_l_odd(m)
                l = lvals(m)%odd(j)
                data_out(m)%data(l,i) = Sum(fodd(:,i)*p_lm_odd(m)%data(:,j))
            Enddo
        Enddo
    Enddo
    Do m = 1, n_m
        Do i = 1, nrhs
            Do j = 1, n_l_even(m)
                l = lvals(m)%even(j)
                data_out(m)%data(l,i) = Sum(feven(:,i)*p_lm_even(m)%data(:,j))
            Enddo
        Enddo
    Enddo

    ! A possible alternative is to simply store the odds up front and evens later
    ! This would require re-interleaving in the end
    !Do m = 1, n_m
    !    Do i = 1, nrhs
    !        Do j = 1, n_ell_odd(m)
    !            !l = lvals(m)%odd(j)
    !            data_out(m)%data(j,i) = Sum(fodd(:,i)*p_lm_odd(m)%data(:,j))
    !        Enddo
    !        offset = n_ell_odd(m)
    !        Do j = 1, n_ell_even(m)
    !            data_out(m)%data(j+offset,i) = Sum(feven(:,i)*p_lm_even(m)%data(:,j))
    !        Enddo
    !    Enddo
    !Enddo

    DeAllocate(feven,fodd)

End Subroutine PtS_2d_drp

Subroutine PtS_2d_dgp(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...DGEMM.... Parity
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),fodd(:,:), feven(:,:)
    Integer :: m,nl,nt1,i,j,l, nt2



    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs))
    Allocate(fodd(1:n_theta/2,1:nrhs))
    feven(1:n_theta/2,:) = data_in(1:n_theta/2,:)
    fodd(1:n_theta/2,:) = data_in(1:n_theta/2,:)
    nt1 = n_theta+1
    nt2 = n_theta/2
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i) = feven(j,i)+data_in(nt1-j,i)
             fodd(j,i) =  fodd(j,i)-data_in(nt1-j,i)
        Enddo
    enddo



    alpha = 1.0d0
    beta = 0.0d0
    Do m = 1, n_m

            If (n_l_even(m) .gt. 0) then
                Allocate(temp(1:n_l_even(m),1:nrhs))
                CALL DGEMM('T','N',n_l_even(m),nrhs,nt2, alpha, p_lm_even(m)%data, nt2,feven , nt2, beta,temp,n_l_even(m))
                Do i =1, nrhs
                Do j = 1, n_l_even(m)
                    l = lvals(m)%even(j)
                    data_out(m)%data(l,i) = temp(j,i)
                Enddo
                Enddo
                DeAllocate(temp)
            Endif

            If (n_l_odd(m) .gt. 0) then
                Allocate(temp(1:n_l_odd(m),1:nrhs))
                CALL DGEMM('T','N',n_l_odd(m),nrhs,nt2, alpha, p_lm_odd(m)%data, nt2,fodd , nt2, beta,temp,n_l_odd(m))
                Do i = 1, nrhs
                Do j = 1, n_l_odd(m)
                    l = lvals(m)%odd(j)
                    data_out(m)%data(l,i) = temp(j,i)
                Enddo
                Enddo
                DeAllocate(temp)
            Endif


    Enddo

    DeAllocate(feven)
    DeAllocate(fodd)

End Subroutine PtS_2d_dgp

Subroutine PtS_2d_dgpv2(data_in, data_out)
    ! Physical-to_Spectral ..2D...DGEMM.... Parity
    ! Data_in is a simple 2-d array (theta is first index)
    ! version 2 --- data_in is a 3D array (ntheta,nrhs,n_m)
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(rmcontainer), Intent(InOut) :: data_out(1:)
    Real*8, Intent(In) :: data_in(:,:,:)
    Integer  :: nrhs
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),fodd(:,:,:), feven(:,:,:)
    Integer :: m,nl,nt1,i,j,l, nt2,ddims(3),k

    ddims = shape(data_in)
    n_m = ddims(3)
    nrhs = ddims(2)
    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs,1:n_m))
    Allocate(fodd(1:n_theta/2,1:nrhs,1:n_m))
    feven(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    fodd(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    nt1 = n_theta+1
    nt2 = n_theta/2
    Do k = 1, n_m
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i,k) = feven(j,i,k)+data_in(nt1-j,i,k)
             fodd(j,i,k) =  fodd(j,i,k)-data_in(nt1-j,i,k)
        Enddo
    enddo
    Enddo


    alpha = 1.0d0
    beta = 0.0d0
    Do m = 1, n_m

            If (n_l_even(m) .gt. 0) then
                Allocate(temp(1:n_l_even(m),1:nrhs))
                CALL DGEMM('T','N',n_l_even(m),nrhs,nt2, alpha, ip_lm_even(m)%data, nt2,feven(:,:,m) , nt2, beta,temp,n_l_even(m))
                Do i =1, nrhs
                Do j = 1, n_l_even(m)
                    l = lvals(m)%even(j)
                    data_out(m)%data(l,i) = temp(j,i)
                Enddo
                Enddo
                DeAllocate(temp)
            Endif

            If (n_l_odd(m) .gt. 0) then
                Allocate(temp(1:n_l_odd(m),1:nrhs))
                CALL DGEMM('T','N',n_l_odd(m),nrhs,nt2, alpha, ip_lm_odd(m)%data, nt2,fodd(:,:,m) , nt2, beta,temp,n_l_odd(m))
                Do i = 1, nrhs
                Do j = 1, n_l_odd(m)
                    l = lvals(m)%odd(j)
                    data_out(m)%data(l,i) = temp(j,i)
                Enddo
                Enddo
                DeAllocate(temp)
            Endif


    Enddo

    DeAllocate(feven)
    DeAllocate(fodd)

End Subroutine PtS_2d_dgpv2


Subroutine StP_2d_dgp(data_in, data_out)
    ! Physical-to_Spectral ..2D...DGEMM.... Parity
    ! Data_in is a spectral structure data_in(m)%data(l,i) ! i is radius or what have you
    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(rmcontainer), Intent(In) :: data_in(:)
    Real*8, Intent(InOut) :: data_out(:,:,:)
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),temp2(:,:)
    Integer :: m,nl,nt1,i,j,l, nt2, ddims(3), nrhs
    ddims = shape(data_out)
    n_m = ddims(3)
    nrhs = ddims(2)
    nt1 = n_theta+1
    nt2 = n_theta/2
    data_out(:,:,:) = 0.0d0


    !////////////////////////////////////
    ! In progress
    alpha = 1.0d0
    beta = 0.0d0
    Allocate(temp(1:nt2,1:nrhs))
    ! Solve for odd and even functions
    Do m = 1, n_m

        If (n_l_even(m) .gt. 0) then
            ! This feels unnecessarily clunky.  Might want to consider storing spectral data as even/odd modes.
            ! Just get it running for now
            Allocate(temp2(1:n_l_even(m),1:nrhs))
            Do i =1, nrhs
            Do j = 1, n_l_even(m)
                l = lvals(m)%even(j)
                temp2(j,i) = data_in(m)%data(l,i)
            Enddo
            Enddo
            CALL DGEMM('T','N',nt2,nrhs,n_l_even(m), alpha, p_lm_even(m)%data, n_l_even(m),temp2 , n_l_even(m), beta,temp,nt2)
            data_out(1:nt2,:,m) = temp    ! store symmetric part in data_out
            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(nt1-j,i,m) = temp(j,i)    ! reflect even modes about equator
                Enddo
            Enddo
            DeAllocate(temp2)
        Endif

        If (n_l_odd(m) .gt. 0) then
            Allocate(temp2(1:n_l_odd(m),1:nrhs))
            Do i =1, nrhs
            Do j = 1, n_l_odd(m)
                l = lvals(m)%odd(j)
                temp2(j,i) = data_in(m)%data(l,i)
            Enddo
            Enddo
            CALL DGEMM('T','N',nt2,nrhs,n_l_odd(m), alpha, p_lm_odd(m)%data, n_l_odd(m),temp2 , n_l_odd(m), beta,temp,nt2)
            Do i = 1, nrhs
                Do j = 1, nt2
                    data_out(j,i,m) = data_out(j,i,m)+temp(j,i)
                    data_out(nt1-j,i,m) = data_out(nt1-j,i,m)-temp(j,i)    ! antisymmetric about equator
                Enddo
            Enddo
            DeAllocate(temp2)
        Endif
    Enddo


    ! Note - not sure if it's faster to make a variable named nt2j1 = nt2-j+1 or just let it compute on the fly
    DeAllocate(temp)


End Subroutine StP_2d_dgp



Subroutine PtS_2d_dr2(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...Direct
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Integer :: m,i,l,th
    ! Does not use the sum function
    !Do m = 0, l_max
    Do m = 1,n_m
        Do i = 1, nrhs
            Do l = m_values(m), l_max
                data_out(m)%data(l,i) =   data_in(1,i)*p_lm(m)%data(1,l)
                Do th = 2, n_theta
                    data_out(m)%data(l,i) = data_out(m)%data(l,i)+data_in(th,i)*p_lm(m)%data(th,l)
                Enddo
            Enddo
        Enddo
    Enddo

End Subroutine PtS_2d_dr2

Subroutine PtS_2d_dr3(data_in, data_out, nrhs)
    ! Physical-to_Spectral ..2D...Direct
    ! Data_in is a simple 2-d array (theta is first index)
    Implicit None
    Type(p_lm_array), Intent(InOut) :: data_out(:)
    Real*8, Intent(In) :: data_in(:,:)
    Integer, Intent(In) :: nrhs
    Integer :: m,i,l,th
    ! Does not use the sum function
    ! Exhanged l and i loops
    !Do m = 0, l_max
    Do m = 1,n_m
        Do l = m_values(m), l_max
        Do i = 1, nrhs
            !Do l = m_values(m), l_max
                data_out(m)%data(l,i) =   data_in(1,i)*p_lm(m)%data(1,l)
                Do th = 2, n_theta
                    data_out(m)%data(l,i) = data_out(m)%data(l,i)+data_in(th,i)*p_lm(m)%data(th,l)
                Enddo
            Enddo
        Enddo
    Enddo

End Subroutine PtS_2d_dr3

!###############################
!  4-D modifications
!  data in spectral space will look like data(m)%(ell,real/imag,radius, field id)
!    data in physical space will look the same (since it is temporary config anyway)
Subroutine PtS_4d_dgpv2(data_in, data_out)
    ! Physical-to_Spectral ..2D...DGEMM.... Parity
    ! Data_in is a simple 2-d array (theta is first index)
    ! version 2 --- data_in is a 3D array (ntheta,nrhs,n_m)

    Implicit None
    ! data in is dimensioned (theta,nrhs)
    !Type(hybrid_m), Intent(InOut) :: data_in(:)
    !Type(spectral_m), Intent(InOut) :: data_out(:)
    Type(rmcontainer4D), Intent(InOut) :: data_out(1:)
    Real*8, Intent(In) :: data_in(:,:,:)
    Integer  :: nrhs
    Real*8 :: alpha, beta
    Real*8, Allocatable :: temp(:,:),fodd(:,:,:), feven(:,:,:)
    Integer :: m,nl,nt1,i,j,l, nt2,ddims(3),k , oddims(4), nfield
    Integer :: rmn, rmx, nr

    oddims = shape(data_out(1)%data)
    nfield = oddims(4)
    rmn = LBOUND(data_out(1)%data,2)
    rmx = UBOUND(data_out(1)%data,2)
    nr = rmx-rmn+1

    ddims = shape(data_in)
    n_m = ddims(3)
    nrhs = ddims(2)
    !To exploit parity, we first need
    !to build the even and odd functions
    Allocate(feven(1:n_theta/2,1:nrhs,1:n_m))
    Allocate(fodd(1:n_theta/2,1:nrhs,1:n_m))
    feven(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    fodd(1:n_theta/2,:,:) = data_in(1:n_theta/2,:,:)
    nt1 = n_theta+1
    nt2 = n_theta/2
    Do k = 1, n_m
    Do i = 1, nrhs
        Do j = 1, n_theta/2
            feven(j,i,k) = feven(j,i,k)+data_in(nt1-j,i,k)
             fodd(j,i,k) =  fodd(j,i,k)-data_in(nt1-j,i,k)
        Enddo
    enddo
    Enddo


    alpha = 1.0d0
    beta = 0.0d0
    Do m = 1, n_m

            If (n_l_even(m) .gt. 0) then
                Allocate(temp(1:n_l_even(m),1:nrhs))
                CALL DGEMM('T','N',n_l_even(m),nrhs,nt2, alpha, ip_lm_even(m)%data, nt2,feven(:,:,m) , nt2, beta,temp,n_l_even(m))

                ! OLD DATA FORMAT FOR REFERENCE
                !Do i =1, nrhs
                !Do j = 1, n_l_even(m)
                !    l = lvals(m)%even(j)
                !    data_out(m)%data(l,i) = temp(j,i)
                !Enddo
                !Enddo

                istart = 1
                iend = nr
                Do f = 1, nfield
                Do im =1, 2
                Do j = 1, n_l_even(m)
                    l = lvals(m)%even(j)
                    data_out(m)%data(l,rmn:rmx,im,f) = temp(j,istart:iend)
                Enddo
                istart = istart+nr
                iend = iend+nr
                Enddo
                Enddo
                DeAllocate(temp)
            Endif

            If (n_l_odd(m) .gt. 0) then
                Allocate(temp(1:n_l_odd(m),1:nrhs))
                CALL DGEMM('T','N',n_l_odd(m),nrhs,nt2, alpha, ip_lm_odd(m)%data, nt2,fodd(:,:,m) , nt2, beta,temp,n_l_odd(m))
                ! AGAIN, for reference
                !Do i = 1, nrhs
                !Do j = 1, n_l_odd(m)
                !    l = lvals(m)%odd(j)
                !    data_out(m)%data(l,i) = temp(j,i)
                !Enddo
                !Enddo
                istart = 1
                iend = nr
                Do f = 1, nfield
                Do im =1, 2
                Do j = 1, n_l_odd(m)
                    l = lvals(m)%odd(j)
                    data_out(m)%data(l,1:nr,im,f) = temp(j,istart:iend)
                Enddo
                istart = istart+nr
                iend = iend+nr
                Enddo
                Enddo

                DeAllocate(temp)
            Endif


    Enddo

    DeAllocate(feven)
    DeAllocate(fodd)

End Subroutine PtS_4d_dgpv2


End Module Legendre_Transforms
