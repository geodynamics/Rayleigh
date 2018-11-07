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

Module Math_Utility
    ! This module contains potentially useful stand alone routines
Contains

    Subroutine Indefinite_Integral(f,intf,xgrid)
        Implicit None
        Real*8, Intent(In) :: f(1:)
        Real*8, Intent(InOut) :: intf(1:), xgrid(1:)
        Real*8 :: trap_area
        Integer :: xsize(1), i, nx
        !Computes the indefinite integral of f using trapezoidal rule
        !Integrates from bottom of grid (nx up to 1 --i.e., assumes reverse grid)
        xsize = size(xgrid)
        nx = xsize(1)
        intf(nx) = 0.0d0
        Do i = nx-1, 1, -1
            trap_area = (xgrid(i)-xgrid(i+1))*(f(i)+f(i+1))*0.5d0
            intf(i) = intf(i+1)+trap_area
        Enddo
    End Subroutine Indefinite_Integral

    Subroutine tanh_profile(x,y,flip)
        Real*8, Intent(In) :: x(1:)
        Real*8, Intent(InOut) :: y(1:)
        Logical, Intent(In), Optional :: flip
        Real*8 :: flip_factor = 1.0d0
        Integer :: xsize(1), nx, i
        xsize = size(x)
        nx = xsize(1)
        if (present(flip)) Then
            if(flip) flip_factor = -1.0d0
        endif

        Do i = 1, nx
            y(i) = 0.5d0*(1.0d0+flip_factor*tanh(x(i)))
        Enddo

    End Subroutine tanh_profile

    Subroutine Spline_Interpolate(profile_in, grid_in, profile_out, grid_out)
        Implicit None
        Real*8, Intent(In) :: profile_in(1:), grid_in(1:), grid_out(1:)
        Real*8, Intent(InOut) :: profile_out(1:)
        Integer :: ngrid_in, ngrid_out, r
        Real*8, Allocatable :: spy2(:)
        Real*8 :: splx, sply, max_grid_in, min_grid_in, plower, pupper



        ngrid_in = size(grid_in)
        ngrid_out = size(grid_out)

        max_grid_in = grid_in(1)
        min_grid_in = grid_in(ngrid_in)

        pupper = profile_in(1)
        plower = profile_in(ngrid_in)

        Allocate(spy2(1:ngrid_in))
        spy2(1:ngrid_in) = 0.0d0
        profile_out(1:ngrid_out) = 0.0d0
        Call Spline(grid_in, profile_in, ngrid_in, 2.0D30, 2.0D30, spy2)
        Do r = 1, ngrid_out
         If ( (grid_out(r) .le. max_grid_in) .and. (grid_out(r) .ge. min_grid_in) ) Then
            splx = grid_out(r)
            Call Splint(grid_in, profile_in,spy2,ngrid_in, splx, sply)
            profile_out(r) = sply
         Endif
        Enddo

        !Handle values that lie outside the bounds of the in_grid
        Do r = 1, ngrid_out
            If (grid_out(r) .ge. max_grid_in) Then
                profile_out(r) = pupper
            Endif
            If (grid_out(r) .le. min_grid_in) Then
                profile_out(r) = plower
            Endif
        Enddo


        DeAllocate(spy2)
    End Subroutine Spline_Interpolate

    ! Numerical Recipes Routines for Spline Interpolation
    Subroutine Spline(x,y,n,yp1,ypn,y2)
        ! From Numerical Recipes in Fortran
        Integer:: n, NMAX
        Real(8) :: yp1, ypn, x(n), y(n), y2(n)
        PARAMETER (NMAX = 10000)
        Integer :: i, k
        Real(8) :: p, qn, sig, un, u(NMAX)

        If (yp1 .gt. 0.99D30) Then
            y2(1) = 0.0D0
            u(1) = 0.0D0
        Else
            y2(1) = -0.5D0
            u(1) = ( 3.0D0 / ( x(2)-x(1) ) ) * ( (y(2)-y(1))/(x(2)-x(1)) -yp1 )
        Endif

        Do i = 2, n-1
            sig = (x(i)-x(i-1)) / (x(i+1)-x(i-1))
            p = sig*y2(i-1)+2.0D0
            y2(i) = (sig-1.0D0)/p
            u(i) = (6.0D0*( (y(i+1)-y(i)) / (x(i+1)-x(i)) - (y(i)-y(i-1)) &
                & /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
        Enddo

        If (ypn .gt. 0.99D30) Then
            qn = 0.0D0
            un = 0.0D0
        Else
            qn = 0.5D0
            un = (3.0D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
        Endif

        y2(n) = (un-qn*u(n-1))/(qn*y2(n-1)+1.0D0)

        Do k = n-1, 1, -1
            y2(k) = y2(k)*y2(k+1)+u(k)
        Enddo
        Return
    End Subroutine Spline

    Subroutine Splint(xa,ya,y2a,n,x,y)
        ! From Numerical Recipes in Fortan
        Integer :: n
        Real(8) x, y, xa(n), y2a(n), ya(n)
        Integer :: k, khi, klo
        Real(8) a,b,h

        klo = 1
        khi = n
1   If ( (khi-klo) .gt. 1) Then
            k = (khi+klo)/2
            If (xa(k) .lt. x) Then        ! if xa is in ascending order, change lt to gt
                khi = k
            else
             klo = k
            endif
            Goto 1
        Endif

        h = xa(khi)-xa(klo)
        If (h .eq. 0.0D0) Then
            Write(6,*) 'bad xa input in splint'
            STOP
        Endif
        a = (xa(khi)-x)/h
        b = (x-xa(klo))/h

        y = a*ya(klo)+ b*ya(khi)+ &
            & ( (a**3-a)*y2a(klo)+(b**3-b)*y2a(khi) )*(h**2)/6.0D0

        Return
    End Subroutine Splint

End Module Math_Utility
