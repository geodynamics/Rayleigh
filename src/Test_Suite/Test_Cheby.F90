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

Module Test_Cheby
    Use ProblemSize

    Implicit None
    Real*8, Allocatable :: colocx(:)
Contains

    Subroutine Test_Chebyshev_Transforms()
        Implicit None
        Allocate(colocx(1:N_R))
        colocx(:) = radius(:)
        !Call Initialize_Chebyshev(colocx,1.0d0,1.5d0)
        !scaling = colocx(1)-colocx(N_R))
        !colocx = 0.5d0*(colocx)
        !write(6,*)'colocx: ', colocx
        Call Test_Transform_1d()
        Call Test_Transform_2d()
        Call Test_Transform_3d()
        Call Test_Transform_4d()
        Call Test_Derivatives()
        DeAllocate(colocx)
    End Subroutine Test_Chebyshev_Transforms

    Subroutine Test_Transform_1d
        Real*8 :: rdiff
        Real*8, Allocatable :: f(:), c(:),fcheck(:)


        Allocate(f(1:N_R))
        Allocate(fcheck(1:N_R))
        Allocate(c(1:N_R))
        f = exp(-4.0d0*colocx**2)    ! Even function
        f = f+exp(-4.0d0*(colocx+0.25)**2)    ! odd function
        f = f-exp(-4.0d0*(colocx-0.25)**2)


        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev transform of 1-D array.'
            Write(6,*)'Parity is off'
            !Call gridcp%to_spectral(f,c)
            !Call gridcp%from_Spectral(c,fcheck)
            Call reldiff1d(f,fcheck,rdiff)
            Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff
        Endif


        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev transform of 1-D array.'
            Write(6,*)'Parity is on'
            c(:) = 0.0d0
            fcheck(:) = 0.0d0
            !Call gridcp%to_spectral(f,c)
            !Call gridcp%from_spectral(c,fcheck)
            Call reldiff1d(f,fcheck,rdiff)
            Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff
        Endif

        DeAllocate(f,c,fcheck)

    End Subroutine Test_Transform_1d

    Subroutine Reldiff1d(f_good,f_check,rdiff)
        Implicit None
        Real*8, Intent(In) :: f_good(:), f_check(:)
        Real*8, Intent(InOut) :: rdiff
        Real*8 :: test, eps
        Integer :: i, n
        rdiff = 0.0d0
        n = size(f_good)
        eps = 1.0d-15
        rdiff = 0.0d0
        do i = 1, n
            test = abs( (f_check(i)-f_good(i))/f_good(i) )
            if (test .gt. rdiff) rdiff = test
        Enddo
    End Subroutine Reldiff1d


    Subroutine Reldiff2d(f_good,f_check,rdiff)
        Implicit None
        Real*8, Intent(In) :: f_good(:,:), f_check(:,:)
        Real*8, Intent(InOut) :: rdiff
        Real*8 :: test, eps
        Integer :: i, j, n, nn, dims(2)
        rdiff = 0.0d0
        dims = shape(f_good)
        n = dims(1)
        nn = dims(2)
        eps = 1.0d-15
        rdiff = 0.0d0
        Do j = 1, nn
        Do i = 1, n
            test = abs( (f_check(i,j)-f_good(i,j))/f_good(i,j) )
            if (test .gt. rdiff) rdiff = test
        Enddo
        Enddo
    End Subroutine Reldiff2d


    Subroutine Test_Transform_2d
        Integer :: i,  n2
        Real*8 :: rdiff
        Real*8, Allocatable :: f(:,:), c(:,:),fcheck(:,:)


        n2 = 3
        Allocate(f(1:N_R,1:n2))
        Allocate(fcheck(1:N_R,1:n2))
        Allocate(c(1:N_R,1:n2))

        Do i = 1, n2
            f(:,i) = exp(-4.0d0*(colocx)**2)
            f(:,i) = f(:,i)+exp(-i*(colocx-0.25)**2) -exp(-i*(colocx-0.25)**2)
        Enddo


        !Call gridcp%to_spectral(f,c)
        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev transform of 2-D array.'
            Write(6,*)'Parity is off'
            !Call gridcp%to_spectral(f,c)
            !Call gridcp%from_spectral(c,fcheck)
            Call reldiff2d(f,fcheck,rdiff)
            Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff
        Endif

        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev transform of 2-D array.'
            Write(6,*)'Parity is on'
            fcheck(:,:) = 0.0d0
            c(:,:) = 0.0d0
            !Call gridcp%to_spectral(f,c)
            !Call gridcp%from_spectral(c,fcheck)
            Call reldiff2d(f,fcheck,rdiff)
            Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff
        Endif

        DeAllocate(f,c,fcheck)

    End Subroutine Test_Transform_2d

    Subroutine Test_Transform_3d
        Integer :: i, j, n2, n3
        Real*8 :: rdiff
        Real*8, Allocatable :: f(:,:,:), c(:,:,:),fcheck(:,:,:)

        n2 = 3
        n3 = 3
        Allocate(f(1:N_R,1:n2,1:n3))
        Allocate(fcheck(1:N_R,1:n2,1:n3))
        Allocate(c(1:N_R,1:n2,1:n3))

        Do j = 1, n3
        Do i = 1, n2
            f(:,i,j) = exp(-4.0d0*(colocx)**2)
            f(:,i,j) = f(:,i,j)+exp(-i*(colocx-0.25)**2) -exp(-i*(colocx-0.25)**2)
            f(:,i,j) = f(:,i,j)+exp(-j*colocx**2)
        Enddo
        Enddo

        !Call gridcp%to_spectral(f,c)
        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev transform of 3-D array.'
            Write(6,*)'Parity is off'
            !Call gridcp%to_spectral(f,c)
            !Call gridcp%from_spectral(c,fcheck)
            Call reldiff3d(f,fcheck,rdiff)
            Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff
        Endif


        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev transform of 3-D array.'
            Write(6,*)'Parity is on'
            !Call gridcp%to_spectral(f,c)
            !Call gridcp%from_spectral(c,fcheck)
            Call reldiff3d(f,fcheck,rdiff)
            Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff
        Endif
        DeAllocate(f,c,fcheck)
    End Subroutine Test_Transform_3d

    Subroutine Reldiff3d(f_good,f_check,rdiff)
        Implicit None
        Real*8, Intent(In) :: f_good(:,:,:), f_check(:,:,:)
        Real*8, Intent(InOut) :: rdiff
        Real*8 :: test, eps
        Integer :: i, j, k, n, nn, nnn, dims(3)
        rdiff = 0.0d0
        dims = shape(f_good)
        n = dims(1)
        nn = dims(2)
        nnn = dims(3)
        eps = 1.0d-8
        rdiff = 0.0d0
        Do k = 1, nnn
        Do j = 1, nn
        Do i = 1, n
            test = abs( (f_check(i,j,k)-f_good(i,j,k))/(f_good(i,j,k)+eps) )
            if (test .gt. rdiff) rdiff = test
        Enddo
        Enddo
        Enddo
    End Subroutine Reldiff3d

    Subroutine Test_Transform_4d
        Integer :: i, j, n2, n3,n4
        Real*8 :: rdiff
        Real*8, Allocatable :: f(:,:,:,:), c(:,:,:,:),fcheck(:,:,:,:)

        n2 = 3
        n3 = 3
        n4 = 2
        Allocate(f(1:N_R,1:n2,1:n3,1:n4))
        Allocate(fcheck(1:N_R,1:n2,1:n3,1:n4))
        Allocate(c(1:N_R,1:n2,1:n3,1:n4))

        Do j = 1, n3
        Do i = 1, n2
            f(:,i,j,1) = exp(-4.0d0*(colocx)**2)
            f(:,i,j,1) = f(:,i,j,1)+exp(-i*(colocx-0.25)**2) -exp(-i*(colocx-0.25)**2)
            f(:,i,j,1) = f(:,i,j,1)+exp(-j*colocx**2)
        Enddo
        Enddo
        f(:,:,:,2) = 2.0d0*f(:,:,:,1)


        Call gridcp%to_spectral(f,c)
        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev transform of 4-D array.'
            Write(6,*)'Parity is off'
            Call gridcp%to_spectral(f,c)
            Call gridcp%from_spectral(c,fcheck)
            Call reldiff4d(f,fcheck,rdiff)
            Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff
        Endif


        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev transform of 4-D array.'
            Write(6,*)'Parity is on'
            fcheck = 0.0d0
            c = 0.0d0
            Call gridcp%to_spectral(f,c)
            Call gridcp%from_spectral(c,fcheck)
            Call reldiff4d(f,fcheck,rdiff)
            Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff
        Endif
        DeAllocate(f,c,fcheck)
    End Subroutine Test_Transform_4d

    Subroutine Reldiff4d(f_good,f_check,rdiff)
        Implicit None
        Real*8, Intent(In) :: f_good(:,:,:,:), f_check(:,:,:,:)
        Real*8, Intent(InOut) :: rdiff
        Real*8 :: test, eps
        Integer :: i, j, k, kk, n, nn, nnn, nnnn, dims(4)
        rdiff = 0.0d0
        dims = shape(f_good)
        n = dims(1)
        nn = dims(2)
        nnn = dims(3)
        nnnn = dims(4)
        eps = 1.0d-12
        rdiff = 0.0d0
        Do kk = 1, nnnn
        Do k = 1, nnn
        Do j = 1, nn
        Do i = 1, n
            test = abs( (f_check(i,j,k,kk)-f_good(i,j,k,kk))/(f_good(i,j,k,kk)+eps) )
            if (test .gt. rdiff) rdiff = test
        Enddo
        Enddo
        Enddo
        Enddo
    End Subroutine Reldiff4d

    Subroutine Test_Derivatives()
        ! Radial derivatives are taken on the first index of 4-D arrays in the main code
        ! We test that functionality here
        Implicit None
        Integer :: i, j, k, kk, n2, n3,n4,nn, dorder,ind,d
        Real*8 :: rdiff
        Real*8, Allocatable :: f(:,:,:,:), c(:,:,:,:),ans(:,:,:,:)
        Real*8 :: worst, worsta
        n2 = 2
        n3 = 2
        n4 = 4
        Allocate(f(1:N_R,1:n2,1:n3,1:n4))
        Allocate(ans(1:N_R,1:n2,1:n3,1:n4))
        Allocate(c(1:N_R,1:n2,1:n3,1:n4))
        f = 0.0d0
        Do j = 1, n3
        Do i = 1, n2
            !f(:,i,j,1) = i*j*exp(-4.0d0*(colocx)**2)
            f(:,i,j,1) = exp(-1.0d0*colocx)
        Enddo
        Enddo
        ans = 0.0d0
        ans(:,:,:,1) = f(:,:,:,1)
        Do j = 1, n3
        Do i = 1, n2
            ans(:,i,j,2) = -8.0d0*colocx*f(:,i,j,1)    ! 1st derivative
            ans(:,i,j,3) = f(:,i,j,1)*(64.0d0*colocx**2-8.0d0)            ! 2nd derivative
            ans(:,i,j,4) = f(:,i,j,1)*(192.0d0*colocx-512.0d0*colocx**3)    ! 3rd derivative
            ans(:,i,j,2) = -f(:,i,j,1)
            ans(:,i,j,3) = f(:,i,j,1)
            ans(:,i,j,4) = -f(:,i,j,1)
        Enddo
        Enddo

        dorder = 1
        Call gridcp%to_spectral(f,c)
        !Out of place derivatives
        Call gridcp%d_by_dr_cp(1,2,c,1)
        Call gridcp%d_by_dr_cp(1,3,c,2)
        Call gridcp%d_by_dr_cp(1,4,c,3)

        Call gridcp%from_spectral(c,f)

        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev Derivative of 4-D array.'
            Do i = 1, 3
                Call reldiff3d2(f(:,:,:,i+1),ans(:,:,:,i+1),rdiff,ind,worst,worsta)
                Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff, maxval(f), maxval(ans)
                !Write(6,*)'Worst results at ind: ', ind, colocx(ind), worst, worsta
            Enddo
            ! Can uncomment this to quickly check visually in IDL.
             Open(unit=15,file='cheby_derivs0',form='unformatted', status='replace')
         Write(15)n_r,n2,n3,n4
            Write(15)((((f(i,j,k,kk),i=1,n_r),j=1,n2),k=1,n3),kk=1,n4)
            Write(15)((((ans(i,j,k,kk),i=1,n_r),j=1,n2),k=1,n3),kk=1,n4)
            Close(15)

        Endif
        !////////////////////////////////////////////////////////////////////////////////////////
        ! Next we peform a similiar test using the dcheby arrays (used for the implicit solve)
        f(:,:,:,2:4) = 0.0d0
        Do d = 1,3
        Do k = 1, n3
            Do j = 1, n2
                Do i = 1, n_r
                    Do nn = 1, n_r
                    f(i,j,k,d+1) = 0 !f(i,j,k,d+1)+dcheby(i,nn,d)*c(nn,j,k,1)
                    Enddo
                Enddo
            Enddo
        Enddo
        Enddo


        If (my_rank .eq. 0) Then
            Write(6,*)'Checking Chebyshev Derivative Array.'
            Do i = 1, 3
                Call reldiff3d2(f(:,:,:,i+1),ans(:,:,:,i+1),rdiff,ind,worst,worsta)
                Write(6,*)'Maximum relative difference (fcheck-f)/f is : ', rdiff
                !Write(6,*)'Worst results at ind: ', ind, colocx(ind), worst, worsta
            Enddo
            ! Can uncomment this to quickly check visually in IDL.
             Open(unit=15,file='cheby_derivs',form='unformatted', status='replace')
         Write(15)n_r,n2,n3,n4
            Write(15)((((f(i,j,k,kk),i=1,n_r),j=1,n2),k=1,n3),kk=1,n4)
            Write(15)((((ans(i,j,k,kk),i=1,n_r),j=1,n2),k=1,n3),kk=1,n4)
            Close(15)
        Endif
        DeAllocate(f,c,ans)
    End Subroutine Test_Derivatives
    Subroutine Reldiff3d2(f_good,f_check,rdiff,ii,worst, worsta)
        Implicit None
        Real*8, Intent(In) :: f_good(:,:,:), f_check(:,:,:)
        Real*8, Intent(InOut) :: rdiff
        Real*8 :: test, eps
        Integer :: i, j, k, n, nn, nnn, dims(3)
        Integer, Intent(Out) :: ii
        Real*8, Intent(Out) :: worst, worsta
        rdiff = 0.0d0
        dims = shape(f_good)
        n = dims(1)
        nn = dims(2)
        nnn = dims(3)
        eps = 1.0d-8
        rdiff = 0.0d0
        Do k = 1, nnn
        Do j = 1, nn
        Do i = 1, n
            test = abs( (f_check(i,j,k)-f_good(i,j,k))/( abs(f_good(i,j,k))+eps) )
            if (test .gt. rdiff) then
                rdiff = test
                ii = i
                worst = f_check(i,j,k)
                worsta = f_good(i,j,k)
            Endif
        Enddo
        Enddo
        Enddo
    End Subroutine Reldiff3d2

End Module Test_Cheby
