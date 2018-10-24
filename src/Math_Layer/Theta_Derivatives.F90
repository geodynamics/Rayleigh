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

Module Theta_Derivatives
    Use Structures
    Implicit None
    Type(rmcontainer), Allocatable :: deriv_coefs(:)
    Integer, Allocatable :: mlocal(:)
    Integer :: nm_local
    Integer, Private :: lmax, tnrl
    Interface d_by_d_theta
        Module Procedure d_dtheta_single,d_dtheta_buffer
    End Interface
Contains

Subroutine Initialize_Theta_Derivatives(mvals,maxl,nrl)
    Implicit None
    Integer, Intent(in) :: mvals(1:)
    Integer, Intent(in) :: maxl,nrl
    Integer :: i, m, l
    Real*8 :: num, denom
    nm_local = size(mvals)
    Allocate(mlocal(1:nm_local))
    mlocal(:) = mvals(:)
    lmax = maxl
    tnrl = 2*nrl    ! 2 *nr_local
    Allocate(deriv_coefs(1:nm_local))

    Do i =1, nm_local
        m = mlocal(i)
        Allocate(deriv_coefs(i)%data(1:2,m:lmax))
        deriv_coefs(i)%data(:,:) = 0.0d0
        Do l = m, lmax
            num = (l+m)*(l-m)*1.0d0
            denom = (2*l+1)*(2*l-1)*1.0d0
            deriv_coefs(i)%data(1,l) = (l-1)*sqrt(num/denom)    ! l-1 coefficient

            num = (l+1+m)*(l+1-m)*1.0d0
            denom = (2*l+1)*(2*l+3)*1.0d0
            deriv_coefs(i)%data(2,l) = -(l+2)*sqrt(num/denom)    ! l+1 coefficient
        Enddo
    Enddo

End Subroutine Initialize_Theta_Derivatives

! Note that in ASH, derivative arrays are calculated out to lmax+1, and I would have extra
! lines below that look like:
! B(i)%data(lmax+1,k) = A(i)%data(lmax,k)*deriv_coefs(i)%data(1,lmax+1) - not doing this yet

Subroutine d_dtheta_single(A,B)
    Implicit None
    Type(rmcontainer), Intent(InOut) :: A(1:),B(1:)
    Integer :: i, m, l, k
    ! Computes B = sin(theta)dA_by_d_theta
    Do i = 1, nm_local
        m = mlocal(i)
        If (m .ne. lmax) Then
            Do k = 1, tnrl
                Do l = m+1, lmax-1
                    B(i)%data(l,k) = A(i)%data(l-1,k)*deriv_coefs(i)%data(1,l)   +    &
                                     A(i)%data(l+1,k)*deriv_coefs(i)%data(2,l)
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
End Subroutine d_dtheta_single


Subroutine d_dtheta_buffer(A,fin,fout)
    Type(rmcontainer), Intent(InOut) :: A(1:)
    Integer, Intent(In) :: fin, fout
    Integer :: i, m, k, l
    Integer :: ind1, ind2
    ind1 = (fin-1)*tnrl
    ind2 = (fout-1)*tnrl
    Do i = 1, nm_local
        m = mlocal(i)
        Do k = 1,tnrl
        Do l = m+1, lmax
           A(i)%data(l,k+ind2) = A(i)%data(l-1,k+ind1)*deriv_coefs(i)%data(1,l)   +   &
                                 A(i)%data(l+1,k+ind1)*deriv_coefs(i)%data(2,l)
        Enddo
        A(i)%data(m,k+ind2)    = A(i)%data(m+1   ,k+ind1)*deriv_coefs(i)%data(2,m)
        A(i)%data(lmax,k+ind2) = A(i)%data(lmax-1,k+ind1)*deriv_coefs(i)%data(1,lmax)
        Enddo
    Enddo
End Subroutine d_dtheta_buffer

End Module Theta_Derivatives
