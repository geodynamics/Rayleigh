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

Module Spectral_Derivatives
    Use Structures
    Implicit None
    ! This module contains d_by_dtheta and d_by_dphi routines
    Type(rmcontainer), Allocatable :: deriv_coefs(:), sderiv_coefs(:)
    Integer, Allocatable :: mlocal(:)
    Integer :: nm_local
    Integer, Private :: lmax, tnrl,nrl

    !These routines compute sin(theta)dA_by_dtheta
    Interface d_by_dtheta
        Module Procedure d_dtheta_single,d_dtheta_buffer, d_dtheta_buff2arr
        Module Procedure d_dtheta_single3D,d_dtheta_buffer4d, d_dtheta_4dbuff2arr
    End Interface

    !These routines compute 1/sin(theta) d(sin^2(theta) A)_by_dtheta
    Interface d_by_sdtheta
        Module Procedure d_sdtheta_single,d_sdtheta_buffer, d_sdtheta_buff2arr
        Module Procedure d_sdtheta_single3D,d_sdtheta_buffer4d, d_sdtheta_4dbuff2arr
    End Interface

    Interface d_by_dphi
        Module Procedure d_by_dphi3D, d_by_dphi_buff2arr, d_by_dphi_rlmbuff, d_by_dphi3Dbuff
        Module Procedure d_by_dphi_4dbuff2arr, d_by_dphi_rlmbuff4d, d_by_dphi_arr2arr
    End Interface

Contains

Subroutine Initialize_Angular_Derivatives(mvals,maxl,nrl0)
    Implicit None
    Integer, Intent(in) :: mvals(1:)
    Integer, Intent(in) :: maxl,nrl0
    Integer :: i, m, l
    Real*8 :: num, denom
    nm_local = size(mvals)
    Allocate(mlocal(1:nm_local))
    mlocal(:) = mvals(:)
    lmax = maxl
    nrl = nrl0
    tnrl = 2*nrl    ! 2 *nr_local
    Allocate(deriv_coefs(1:nm_local))    ! sin(theta)dA_by_dtheta coefs
    Allocate(sderiv_coefs(1:nm_local))  ! 1/sin(theta) d(sin^2(theta) A)_by_dtheta coefs

    Do i =1, nm_local
        m = mlocal(i)
        Allocate(deriv_coefs(i)%data(1:2,m:lmax))
        Allocate(sderiv_coefs(i)%data(1:2,m:lmax))
        deriv_coefs(i)%data(:,:) = 0.0d0
        sderiv_coefs(i)%data(:,:) = 0.0d0
        Do l = m, lmax
            num = (l+m)*(l-m)*1.0d0
            denom = (2*l+1)*(2*l-1)*1.0d0
            deriv_coefs(i)%data(1,l) = (l-1)*sqrt(num/denom)    ! l-1 coefficient
            sderiv_coefs(i)%data(1,l) = (l+1)*sqrt(num/denom)

            num = (l+1+m)*(l+1-m)*1.0d0
            denom = (2*l+1)*(2*l+3)*1.0d0
            deriv_coefs(i)%data(2,l) = -(l+2)*sqrt(num/denom)    ! l+1 coefficient
            sderiv_coefs(i)%data(2,l) = -l*sqrt(num/denom)
        Enddo
    Enddo

End Subroutine Initialize_Angular_Derivatives



!////////////////////////////
! Theta derivatives



! Note that in ASH, derivative arrays are calculated out to lmax+1, and I would have extra
! lines below that look like:
! B(i)%data(lmax+1,k) = A(i)%data(lmax,k)*deriv_coefs(i)%data(1,lmax+1) - not doing this yet

Subroutine d_dtheta_single(A,B)
    Implicit None
    Type(rmcontainer), Intent(InOut) :: A(1:),B(1:)
    Integer :: i, m, l, k
    ! Computes B = sin(theta)dA_by_d_theta
    !$OMP PARALLEL DO PRIVATE(i,m,k,l)
    Do i = 1, nm_local
        m = mlocal(i)
        If (m .ne. lmax) Then
            Do k = 1, tnrl
                Do l = m+1, lmax-1
                    B(i)%data(l,k) = A(i)%data(l-1,k)*deriv_coefs(i)%data(1,l)   +   A(i)%data(l+1,k)*deriv_coefs(i)%data(2,l)
                Enddo
                B(i)%data(m,k)    = A(i)%data(m+1   ,k)*deriv_coefs(i)%data(2,m)
                B(i)%data(lmax,k) = A(i)%data(lmax-1,k)*deriv_coefs(i)%data(1,lmax)

            Enddo
        Else
            Do k = 1, tnrl
                B(i)%data(lmax,k) = 0.0d0
            Enddo
        Endif
    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_dtheta_single

!//////////////////////////////
!    For New Layout
Subroutine d_dtheta_single3D(A,B)
    Implicit None
    Type(rmcontainer3D), Intent(InOut) :: A(1:),B(1:)
    Integer :: i, m, l
    ! Computes B = sin(theta)dA_by_d_theta
    !$OMP PARALLEL DO PRIVATE(i,m,l)
       Do i = 1, nm_local
        m = mlocal(i)
        If (m .ne. lmax) Then
            Do l = m+1, lmax-1
                B(i)%data(l,:,:) = A(i)%data(l-1,:,:)*deriv_coefs(i)%data(1,l)   +   A(i)%data(l+1,:,:)*deriv_coefs(i)%data(2,l)
            Enddo
            B(i)%data(m,:,:)    = A(i)%data(m+1   ,:,:)*deriv_coefs(i)%data(2,m)
            B(i)%data(lmax,:,:) = A(i)%data(lmax-1,:,:)*deriv_coefs(i)%data(1,lmax)

        Else
            B(i)%data(lmax,:,:) = 0.0d0
        Endif
    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_dtheta_single3D

Subroutine d_sdtheta_single(A,B)
    Implicit None
    Type(rmcontainer), Intent(InOut) :: A(1:),B(1:)
    Integer :: i, m, l, k
    ! Computes B = 1/sin(theta)d(sin^2 A)_by_d_theta
    !$OMP PARALLEL DO PRIVATE(i,m,k,l)
    Do i = 1, nm_local
        m = mlocal(i)
        If (m .ne. lmax) Then
            Do k = 1, tnrl
                Do l = m+1, lmax-1
                    B(i)%data(l,k) = A(i)%data(l-1,k)*sderiv_coefs(i)%data(1,l)   +   A(i)%data(l+1,k)*sderiv_coefs(i)%data(2,l)
                Enddo
                B(i)%data(m,k)    = A(i)%data(m+1   ,k)*sderiv_coefs(i)%data(2,m)
                B(i)%data(lmax,k) = A(i)%data(lmax-1,k)*sderiv_coefs(i)%data(1,lmax)

            Enddo
        Else
            Do k = 1, tnrl
                B(i)%data(lmax,k) = 0.0d0
            Enddo
        Endif
    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_sdtheta_single

!/////// For new layout
subroutine d_sdtheta_single3D(A,B)
    Implicit None
    Type(rmcontainer3D), Intent(InOut) :: A(1:),B(1:)
    Integer :: i, m, l
    ! Computes B = 1/sin(theta)d(sin^2 A)_by_d_theta
    !$OMP PARALLEL DO PRIVATE(i,m,l)
    Do i = 1, nm_local
        m = mlocal(i)
        If (m .ne. lmax) Then
            Do l = m+1, lmax-1
                B(i)%data(l,:,:) = A(i)%data(l-1,:,:)*sderiv_coefs(i)%data(1,l)   +   A(i)%data(l+1,:,:)*sderiv_coefs(i)%data(2,l)
            Enddo
            B(i)%data(m,:,:)    = A(i)%data(m+1   ,:,:)*sderiv_coefs(i)%data(2,m)
            B(i)%data(lmax,:,:) = A(i)%data(lmax-1,:,:)*sderiv_coefs(i)%data(1,lmax)

        Else
            B(i)%data(lmax,:,:) = 0.0d0
        Endif
    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_sdtheta_single3D


Subroutine d_dtheta_buffer(A,fin,fout)
    Type(rmcontainer), Intent(InOut) :: A(1:)
    Integer, Intent(In) :: fin, fout
    Integer :: i, m, k, l
    Integer :: ind1, ind2
    ind1 = (fin-1)*tnrl
    ind2 = (fout-1)*tnrl
    !$OMP PARALLEL DO PRIVATE(i,m,k,l)
    Do i = 1, nm_local
        m = mlocal(i)
        If ( m .ne. lmax) Then
        Do k = 1,tnrl
        Do l = m+1, lmax-1
            A(i)%data(l,k+ind2) = A(i)%data(l-1,k+ind1)*deriv_coefs(i)%data(1,l) &
                                  + A(i)%data(l+1,k+ind1)*deriv_coefs(i)%data(2,l)
        Enddo
        A(i)%data(m,k+ind2)    = A(i)%data(m+1   ,k+ind1)*deriv_coefs(i)%data(2,m)
        A(i)%data(lmax,k+ind2) = A(i)%data(lmax-1,k+ind1)*deriv_coefs(i)%data(1,lmax)
        Enddo
        Else
        A(i)%data(m,ind2+1:ind2+tnrl) = 0.0d0
        Endif

    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_dtheta_buffer

!///// For New layout
Subroutine d_dtheta_buffer4D(A,fin,fout)
    Type(rmcontainer4d), Intent(InOut) :: A(1:)
    Integer, Intent(In) :: fin, fout
    Integer :: i, m, l
    Integer :: ind1, ind2
    ind1 = (fin-1)*tnrl
    ind2 = (fout-1)*tnrl
    !$OMP PARALLEL DO PRIVATE(i,m,l)
    Do i = 1, nm_local
        m = mlocal(i)
        If ( m .ne. lmax) Then

        Do l = m+1, lmax-1
            A(i)%data(l,:,:,fout) = A(i)%data(l-1,:,:,fin)*deriv_coefs(i)%data(1,l) &
            + A(i)%data(l+1,:,:,fin)*deriv_coefs(i)%data(2,l)
        Enddo
        A(i)%data(m,:,:,fout)    = A(i)%data(m+1   ,:,:,fin)*deriv_coefs(i)%data(2,m)
        A(i)%data(lmax,:,:,fout) = A(i)%data(lmax-1,:,:,fin)*deriv_coefs(i)%data(1,lmax)

        Else
            A(i)%data(m,:,:,fout) = 0.0d0
        Endif

    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_dtheta_buffer4D

Subroutine d_sdtheta_buffer(A,fin,fout)
    Type(rmcontainer), Intent(InOut) :: A(1:)
    Integer, Intent(In) :: fin, fout
    Integer :: i, m, k, l
    Integer :: ind1, ind2
    ind1 = (fin-1)*tnrl
    ind2 = (fout-1)*tnrl
    !$OMP PARALLEL DO PRIVATE(i,m,k,l)
    Do i = 1, nm_local
        m = mlocal(i)
        If (m .ne. lmax) Then
        Do k = 1,tnrl
        Do l = m+1, lmax-1
            A(i)%data(l,k+ind2) = A(i)%data(l-1,k+ind1)*sderiv_coefs(i)%data(1,l) &
                                  + A(i)%data(l+1,k+ind1)*sderiv_coefs(i)%data(2,l)
        Enddo
        A(i)%data(m   ,k+ind2) = A(i)%data(m+1   ,k+ind1)*sderiv_coefs(i)%data(2,m)
        A(i)%data(lmax,k+ind2) = A(i)%data(lmax-1,k+ind1)*sderiv_coefs(i)%data(1,lmax)
        Enddo
        Else
            A(i)%data(m,ind2+1:ind2+tnrl) = 0.0d0
        Endif
    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_sdtheta_buffer

!//////////////////////
!  For new layout
subroutine d_sdtheta_buffer4d(A,fin,fout)
    Type(rmcontainer4d), Intent(InOut) :: A(1:)
    Integer, Intent(In) :: fin, fout
    Integer :: i, m, l
    !$OMP PARALLEL DO PRIVATE(i,m,l)
    Do i = 1, nm_local
        m = mlocal(i)
        If (m .ne. lmax) Then
        Do l = m+1, lmax-1
            A(i)%data(l,:,:,fout) = A(i)%data(l-1,:,:,fin)*sderiv_coefs(i)%data(1,l) &
                                    + A(i)%data(l+1,:,:,fin)*sderiv_coefs(i)%data(2,l)
        Enddo
        A(i)%data(m   ,:,:,fout) = A(i)%data(m+1   ,:,:,fin)*sderiv_coefs(i)%data(2,m)
        A(i)%data(lmax,:,:,fout) = A(i)%data(lmax-1,:,:,fin)*sderiv_coefs(i)%data(1,lmax)
        Else
            A(i)%data(m,:,:,fout) = 0.0d0
        Endif
    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_sdtheta_buffer4d


Subroutine d_dtheta_buff2arr(A,fin,arr)
    Type(rmcontainer), Intent(InOut) :: A(1:), arr(1:)
    Integer, Intent(In) :: fin
    Integer :: i, m, k, l
    Integer :: ind1
    ind1 = (fin-1)*tnrl
    !$OMP PARALLEL DO PRIVATE(i,m,k,l)
    Do i = 1, nm_local
        m = mlocal(i)
        if (m .ne. lmax) then
            Do k = 1,tnrl
            Do l = m+1, lmax-1
                arr(i)%data(l,k) = A(i)%data(l-1,k+ind1)*deriv_coefs(i)%data(1,l) &
                                   + A(i)%data(l+1,k+ind1)*deriv_coefs(i)%data(2,l)
            Enddo
            arr(i)%data(m,k)    = A(i)%data(m+1   ,k+ind1)*deriv_coefs(i)%data(2,m)
            arr(i)%data(lmax,k) = A(i)%data(lmax-1,k+ind1)*deriv_coefs(i)%data(1,lmax)
            Enddo

        else
            arr(i)%data(m,:) = 0.0d0
        endif

    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_dtheta_buff2arr

!///// For New layout
Subroutine d_dtheta_4dbuff2arr(A,fin,arr)
    Type(rmcontainer4d), Intent(InOut) :: A(1:)
    Type(rmcontainer3d), Intent(InOut) :: arr(1:)
    Integer, Intent(In) :: fin
    Integer :: i, m, l
    !$OMP PARALLEL DO PRIVATE(i,m,l)
    Do i = 1, nm_local
        m = mlocal(i)
        if (m .ne. lmax) then
            Do l = m+1, lmax-1
                arr(i)%data(l,:,:) = A(i)%data(l-1,:,:,fin)*deriv_coefs(i)%data(1,l)& 
                                     + A(i)%data(l+1,:,:,fin)*deriv_coefs(i)%data(2,l)
            Enddo
            arr(i)%data(m,:,:)    = A(i)%data(m+1   ,:,:,fin)*deriv_coefs(i)%data(2,m)
            arr(i)%data(lmax,:,:) = A(i)%data(lmax-1,:,:,fin)*deriv_coefs(i)%data(1,lmax)
        else
            arr(i)%data(m,:,:) = 0.0d0
        endif

    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_dtheta_4dbuff2arr




Subroutine d_sdtheta_buff2arr(A,fin,arr)
    Type(rmcontainer), Intent(InOut) :: A(1:), arr(1:)
    Integer, Intent(In) :: fin
    Integer :: i, m, k, l
    Integer :: ind1
    ind1 = (fin-1)*tnrl
    !$OMP PARALLEL DO PRIVATE(i,m,k,l)
    Do i = 1, nm_local
        m = mlocal(i)
        if (m .ne. lmax) then
        Do k = 1,tnrl
        Do l = m+1, lmax-1
            arr(i)%data(l,k) = A(i)%data(l-1,k+ind1)*sderiv_coefs(i)%data(1,l)   +   A(i)%data(l+1,k+ind1)*sderiv_coefs(i)%data(2,l)
        Enddo
        arr(i)%data(m,k)    = A(i)%data(m+1   ,k+ind1)*sderiv_coefs(i)%data(2,m)
        arr(i)%data(lmax,k) = A(i)%data(lmax-1,k+ind1)*sderiv_coefs(i)%data(1,lmax)
        Enddo
        Else
            arr(i)%data(m,:) = 0.0d0
        Endif
    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_sdtheta_buff2arr

!///////////////////////////////////
! For new layout
Subroutine d_sdtheta_4dbuff2arr(A,fin,arr)
    Type(rmcontainer4d), Intent(InOut) :: A(1:)
    Type(rmcontainer3d), Intent(InOut) :: arr(1:)
    Integer, Intent(In) :: fin
    Integer :: i, m, l
    !$OMP PARALLEL DO PRIVATE(i,m,l)
    Do i = 1, nm_local
        m = mlocal(i)
        if (m .ne. lmax) then
        Do l = m+1, lmax-1
            arr(i)%data(l,:,:) = A(i)%data(l-1,:,:,fin)*sderiv_coefs(i)%data(1,l) &
                                 + A(i)%data(l+1,:,:,fin)*sderiv_coefs(i)%data(2,l)
        Enddo
        arr(i)%data(m,:,:)    = A(i)%data(m+1   ,:,:,fin)*sderiv_coefs(i)%data(2,m)
        arr(i)%data(lmax,:,:) = A(i)%data(lmax-1,:,:,fin)*sderiv_coefs(i)%data(1,lmax)
        Else
            arr(i)%data(m,:,:) = 0.0d0
        Endif
    Enddo
    !$OMP END PARALLEL DO
End Subroutine d_sdtheta_4dbuff2arr
!///////////////////////////
! Phi derivatives
    !/////////////////////////////////////////////
    !  d_by_dm3D
    !  Computes derivative of arrin along first dimension
    !  Arrin is assumed to be in spectral space in dimension 1
    !  Storage of m values is assumed to be in FFTW's r2c in place format
    !    Note that this has no spatial scale factor (d_by_dphi vs d_by_dx)
    Subroutine d_by_dphi3D(arrin,arrout)
        Implicit None
        Integer :: i, j, k, m, ashape(1:3),ni,nj,nk
        Real*8, Intent(InOut) :: arrin(:,1:,1:)
        Real*8, Intent(InOut) :: arrout(:,1:,1:)
        ashape = shape(arrin)
        ni = ashape(2)
        nj = ashape(3)
        nk = ashape(1)
        !$OMP PARALLEL DO PRIVATE(i,j,k,m)
        Do j = 1, nj
            Do i = 1, ni
                Do k = 1, nk,2
                    m = (k-1)/2
                    arrout(k,i,j) = -m*arrin(k+1,i,j)
                    arrout(k+1,i,j) = m*arrin(k,i,j)
                Enddo
            Enddo
        Enddo
        !$OMP END PARALLEL DO
    End Subroutine d_by_dphi3D

    Subroutine d_by_dphi3Dbuff(arrin,from_ind,to_ind)
        Implicit None
        Integer :: i, j, k, m, ashape(1:4),ni,nj,nk
        Real*8, Intent(InOut) :: arrin(:,1:,1:,1:)
        Integer, Intent(In) :: from_ind, to_ind
        Real*8 :: tmp
        ashape = shape(arrin)
        ni = ashape(2)
        nj = ashape(3)
        nk = ashape(1)
        If (from_ind .ne. to_ind) Then
            !$OMP PARALLEL DO PRIVATE(i,j,k,m)
            Do j = 1, nj
                Do i = 1, ni
                    Do k = 1, nk,2
                        m = (k-1)/2
                        arrin(k,i,j,to_ind) = -m*arrin(k+1,i,j,from_ind)
                        arrin(k+1,i,j,to_ind) = m*arrin(k,i,j,from_ind)
                    Enddo
                Enddo
            Enddo
            !$OMP END PARALLEL DO
        Else
            !$OMP PARALLEL DO PRIVATE(i,j,k,m,tmp)
            Do j = 1, nj
                Do i = 1, ni
                    Do k = 1, nk,2
                        m = (k-1)/2
                        tmp = arrin(k,i,j,to_ind)
                        arrin(k,i,j,to_ind) = -m*arrin(k+1,i,j,from_ind)
                        arrin(k+1,i,j,to_ind) = m*tmp
                    Enddo
                Enddo
            Enddo
            !$OMP END PARALLEL DO
        Endif
    End Subroutine d_by_dphi3Dbuff


    !///////////////////////////////////////////////////////
    !  Computes phi derivative of variable from_ind stored in buff
    !  Stores result in arr.  Arr and buff are assumed to be in rlm config
    !  For rlm config, 2nd index runs 1:nr (real part) and nr+1:2nr (imaginary_part)
    Subroutine d_by_dphi_buff2arr(buff,from_ind, arr)
        Implicit None
        Type(rmcontainer), Intent(InOut) :: buff(1:), arr(1:)
        Integer, Intent(In) :: from_ind
        Integer :: rrstart, rmid, irend, i, m
        rrstart = (from_ind-1)*tnrl+1
        rmid = rrstart+nrl
        irend = rrstart+tnrl-1
        !$OMP PARALLEL DO PRIVATE(i,m)
        Do i = 1, nm_local
                m = mlocal(i)
                arr(i)%data(:,1:nrl)      = -m*buff(i)%data(:,rmid:irend)
                arr(i)%data(:,nrl+1:tnrl) =  m*buff(i)%data(:,rrstart:rmid-1)
        Enddo
        !$OMP END PARALLEL DO

    End Subroutine d_by_dphi_buff2arr

    !/////////////////
    ! For new layout
    Subroutine d_by_dphi_4dbuff2arr(buff,from_ind, arr)
        Implicit None
        Type(rmcontainer4d), Intent(InOut) :: buff(1:)
        Type(rmcontainer3d), Intent(InOut) :: arr(1:)
        Integer, Intent(In) :: from_ind
        Integer :: i, m
        !$OMP PARALLEL DO PRIVATE(i,m)
        Do i = 1, nm_local
                m = mlocal(i)
                arr(i)%data(:,:,1) = -m*buff(i)%data(:,:,2,from_ind)
                arr(i)%data(:,:,2) =  m*buff(i)%data(:,:,1,from_ind)
        Enddo
        !$OMP END PARALLEL DO

    End Subroutine d_by_dphi_4dbuff2arr

    Subroutine d_by_dphi_arr2arr(arr1, arr2)
        Implicit None

        Type(rmcontainer3d), Intent(InOut) :: arr1(1:), arr2(1:)

        Integer :: i, m
        !$OMP PARALLEL DO PRIVATE(i,m)
        Do i = 1, nm_local
                m = mlocal(i)
                arr2(i)%data(:,:,1) = -m*arr1(i)%data(:,:,2)
                arr2(i)%data(:,:,2) =  m*arr1(i)%data(:,:,1)
        Enddo
        !$OMP END PARALLEL DO

    End Subroutine d_by_dphi_arr2arr

    !///////////////////////////////////////////////////////
    !  Computes phi derivative of variable from_ind stored in buff
    !  Stores result in variable to_ind of buff arr.
    !  Buff is assumed to be in rlm config
    !  from_ind and to_ind can be the same (in-place)
    Subroutine d_by_dphi_rlmbuff(buff,from_ind, to_ind)
        Implicit None
        Type(rmcontainer), Intent(InOut) :: buff(1:)
        Integer, Intent(In) :: from_ind, to_ind
        Integer :: rrstart, rmid, irend, rrstart2, rmid2, irend2, i, m
        Real*8, Allocatable :: temp(:,:)

        rrstart = (from_ind-1)*tnrl+1
        rmid = rrstart+nrl
        irend = rrstart+tnrl-1

        If (from_ind .ne. to_ind) Then
            rrstart2 = (to_ind-1)*tnrl+1
            rmid2 = rrstart2+nrl
            irend2 = rrstart2+tnrl-1
            !$OMP PARALLEL DO PRIVATE(i,m)
            Do i = 1, nm_local
                m = mlocal(i)
                buff(i)%data(:,rrstart2:rmid2-1) = -m*buff(i)%data(:,rmid:irend)
                buff(i)%data(:,rmid2:irend2) = m*buff(i)%data(:,rrstart:rmid-1)
            Enddo
            !$OMP END PARALLEL DO
        Else
            ! in-place
            Allocate(temp(0:lmax,1:nrl))
            !$OMP PARALLEL DO PRIVATE(i,m,temp)
            Do i = 1, nm_local
                m = mlocal(i)
                temp(m:lmax,1:nrl) = buff(i)%data(m:lmax,rrstart:rmid-1)    ! save the real piece
                buff(i)%data(:,rrstart:rmid-1) = -m*buff(i)%data(:,rmid:irend) !overwrite it with new real piece
                buff(i)%data(m:lmax,rmid:irend) = m*temp(m:lmax,1:nrl)    ! build new imaginary piece

            Enddo
            !OMP END PARALLEL DO
            DeAllocate(temp)
        Endif

    End Subroutine d_by_dphi_rlmbuff

    !///////////////////////////
    ! For new layout
    Subroutine d_by_dphi_rlmbuff4d(buff,from_ind, to_ind)
        Implicit None
        Type(rmcontainer4d), Intent(InOut) :: buff(1:)
        Integer, Intent(In) :: from_ind, to_ind
        Integer ::  i, m
        Real*8, Allocatable :: temp(:,:)


        If (from_ind .ne. to_ind) Then
            !$OMP PARALLEL DO PRIVATE(i,m)
            Do i = 1, nm_local
                m = mlocal(i)
                buff(i)%data(:,:,1,to_ind) = -m*buff(i)%data(:,:,2,from_ind)
                buff(i)%data(:,:,2,to_ind) =  m*buff(i)%data(:,:,1,from_ind)
            Enddo
            !$OMP END PARALLEL DO
        Else
            ! in-place
            Allocate(temp(0:lmax,1:nrl))
            !$OMP PARALLEL DO PRIVATE(i,m,temp)
            Do i = 1, nm_local
                m = mlocal(i)
                temp(m:lmax,:) = buff(i)%data(m:lmax,:,1,from_ind)
                buff(i)%data(:,:,1,from_ind) = -m*buff(i)%data(:,:,2,from_ind)
                buff(i)%data(m:lmax,:,2,from_ind) = m*temp(m:lmax,:)
            Enddo
            !$OMP END PARALLEL DO
            DeAllocate(temp)
        Endif

    End Subroutine d_by_dphi_rlmbuff4d

End Module Spectral_Derivatives
