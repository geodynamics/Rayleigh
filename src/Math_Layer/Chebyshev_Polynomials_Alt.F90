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

Module Chebyshev_Polynomials_Alt
    ! Module for computing Chebyshev Polynomial Arrays
    ! Uses objects approach to facilitate interpolation between radial grids of different dimensions
    Implicit None
    Type, Public :: Cheby_Transform_Interface
        Integer :: N_max, N_even, N_odd, n_x

        Real*8, Allocatable :: x(:)  ! The colocation points
        Real*8, Allocatable :: cheby(:,:)        ! cheby(r,k) is chebyshev polynomial of degree k-1 at radius r
        Real*8, Allocatable :: cheby_even(:,:), cheby_odd(:,:) ! even and odd chebyshev arrays
        Real*8, Allocatable :: dcheby(:,:,:)    ! The Chebyshev Derivative Arrays
        Real*8, Allocatable :: integration_weights(:)
        Real*8  :: pi_over_N
        Logical :: Parity = .true.
        Logical :: initialized = .false.
        Real*8  :: scaling ! x runs from -0.5 to 0.5 by default

        Contains

        Procedure :: Init => Init_Cheby_Interface
        Procedure :: gen_colocation_points
        Procedure :: gen_Tn
        Procedure :: gen_Tn_Deriv_arrays
        Procedure :: Destroy
        Procedure :: tospec4d => To_Spectral_4D
    End Type Cheby_Transform_Interface
    Real*8, Private ::    Pi  = 3.1415926535897932384626433832795028841972d+0


    !Interface Cheby_To_Spectral
    !    Module Procedure To_Spectral_1D, To_Spectral_2D, To_Spectral_3D, To_Spectral_4D
    !End Interface


    !Interface Cheby_From_Spectral
    !    Module Procedure From_Spectral_1D, From_Spectral_2D, From_Spectral_3D, From_Spectral_4D
    !End Interface

    !Interface d_by_dr_cpalt
    !    Module Procedure Cheby_Deriv_Buffer_4Dalt
    !End Interface

Contains
    Subroutine Destroy(self)
        !Cleanup routine that frees all memory used by the Chebyshev Interface object
        Implicit None
        Class(Cheby_Transform_Interface) :: self
        DeAllocate(self%x,self%cheby)
        If (allocated(self%cheby_even)) DeAllocate(self%cheby_even)
        If (allocated(self%cheby_odd)) DeAllocate(self%cheby_odd)
        If (allocated(self%dcheby)) DeAllocate(self%dcheby)
        If (allocated(self%integration_weights)) DeAllocate(self%integration_weights)
        self%initialized = .false.
    End Subroutine Destroy
    Subroutine Init_Cheby_Interface(self,grid, xmin,xmax, integral, derivatives)
        Implicit None
        Class(Cheby_Transform_Interface) :: self
        Real*8, Intent(InOut) :: grid(:)
        Real*8, Intent(In), Optional :: xmin, xmax
        Real*8 ::delta, gmin, tmp, xx
        Integer :: r, n_max
        Logical, Intent(In), Optional :: integral, derivatives
        n_max = size(grid)
        If (.not. self%initialized) Then
            self%N_max = size(grid)
            Call self%gen_colocation_points()
            grid(:) = self%x(:)
            Call self%Gen_Tn()
            If (present(derivatives)) Then
                If (derivatives) Then
                    Call self%Gen_Tn_Deriv_Arrays(3)
                Endif
            Endif
            If (present(xmin)) Then
                delta = xmax-xmin
                self%scaling = (self%x(1)-self%x(N_max))/delta

                If (present(derivatives)) Then
                    If (derivatives) Then
                        self%dcheby(:,:,1) = self%dcheby(:,:,1)*self%scaling
                        self%dcheby(:,:,2) = self%dcheby(:,:,2)*(self%scaling**2)
                        self%dcheby(:,:,3) = self%dcheby(:,:,3)*(self%scaling**3)
                        grid(:) = grid(:)/self%scaling
                        gmin = grid(N_max)
                        grid(:) = grid(:)-gmin+xmin
                    Endif
                Endif
            Endif
            self%initialized = .true.
        Else
            grid(:) = self%x(:)
        Endif


        If (present(integral)) Then
            If (integral) Then
                Allocate(self%integration_weights(1:self%n_max))
                self%integration_weights = 0.0d0
                tmp = 1.5d0*Pi * (grid(1)-grid(N_max)) / &
                    & ( (grid(1)**3 - grid(N_max)**3) * N_max )
                Do r=1,N_max
                    xx = (2.0d0*grid(r)-grid(N_max)-grid(1))/(grid(1)-grid(N_max))
                    self%integration_weights(r) = grid(r)**2 * tmp * sqrt(1.0d0-xx*xx)
                Enddo
            Endif
        Endif

    End Subroutine Init_Cheby_Interface

    Subroutine Rescale_Grid_CP(self,length_scale)
        Implicit None
        Class(Cheby_Transform_Interface) :: self
        Real*8, Intent(In) :: length_scale
        ! Following initialization, we can rescale the chebyshev arrays if we choose
        ! This is useful when nondimensionalizing after the reference state has been set up
        ! (which typically requires a radial grid to have been established)

        self%dcheby(:,:,1) = self%dcheby(:,:,1)*length_scale
        self%dcheby(:,:,2) = self%dcheby(:,:,2)*(length_scale**2)
        self%dcheby(:,:,3) = self%dcheby(:,:,3)*(length_scale**3)

    End Subroutine Rescale_Grid_CP

    Subroutine Gen_Colocation_Points(self)
        Implicit None
        Integer :: i
        Real*8 :: arg
        Class(Cheby_Transform_Interface) :: self
        Allocate(self%x(1:self%N_max))
        self%pi_over_N = pi/(self%N_max*1.0d0)
        arg = (0.5d0)*self%pi_over_N
        Do i = 1, self%N_Max
            self%x(i) = cos(arg)
            arg = arg+self%pi_over_N
        Enddo
    End Subroutine Gen_Colocation_Points

    Subroutine Gen_Tn(self)
        Implicit None
        Class(Cheby_Transform_Interface) :: self
        Integer :: i, k, r, n_max, n_odd, n_even, n_x
        Real*8 :: acx, arg
        n_max = self%n_max
        Allocate(self%cheby(1:N_max,1:N_max))
        Do r = 1, N_max
            acx = self%pi_over_n*(r-1+0.5d0)
            Do i = 1, N_max
                k = i -1
                arg = k*acx
                self%cheby(r,i) = cos(arg)
            Enddo
        Enddo
        If (self%parity) Then
            n_odd = N_max/2
            n_even = n_odd+mod(N_max,2)
            self%n_even = n_even
            self%n_odd = n_odd
            n_x = n_even
            self%n_x = n_x
               Allocate(self%cheby_even(1:N_x,1:N_even))
            Allocate(self%cheby_odd(1:N_x,1:N_odd))
            Do i = 1, n_even
                self%cheby_even(1:N_x,i) = self%cheby(1:N_x,2*i-1)
            Enddo
            Do i = 1,n_odd
                self%cheby_odd(1:N_x,i) = self%cheby(1:N_x,2*i)
            Enddo
            If (n_x .ne. n_odd) Then
                ! We actually have an x = 0 point
                ! We must be careful not to double count the power here
                ! when exploiting parity to speed up the transforms.
                ! No such adjustment need be made for the regular cheby array
                self%cheby_even(n_x,:) = 0.5d0*self%cheby_even(n_x,:)
                self%cheby_odd(n_x,:) = 0.5d0*self%cheby_odd(n_x,:)
                ! There is really no need to modify cheby_odd, but this way it is consistently stored.
            Endif
        Endif
    End Subroutine Gen_Tn

    Subroutine Gen_Tn_Deriv_Arrays(self,dmax)
        Implicit None
        Class(Cheby_Transform_Interface) :: self
        Integer, Intent(In) :: dmax
        Integer :: i, k,n,d, n_max
        Real*8, Allocatable :: alpha(:,:)

        ! sum_n (alpha_kn  c_n) = c'_k
        n_max = self%n_max
        Allocate(alpha(1:N_max,1:N_max))
        alpha(:,:) = 0.0d0
        alpha(N_max,:) = 0.0d0
        alpha(N_max-1,N_max) = 2.0d0*(N_max-1)
        Do k = N_max-2, 1, -1
            alpha(k,k+1) = 2.0d0*k
            alpha(k,:) = alpha(k,:)+alpha(k+2,:)
        Enddo
        Allocate(self%dcheby(1:N_max,1:N_max,0:dmax))
        self%dcheby(:,:,:) = 0.0d0
        self%dcheby(:,:,0) = self%cheby(:,:)
        self%dcheby(:,1,0) = self%dcheby(:,1,0) - 0.5d0    ! This accounts for the -1/2c_0 in going from_spectral
        If (dmax .ge. 1) Then
            Do d = 1, dmax
            Do n = 1, N_Max
            Do k = 1, N_max
                Do i = 1, N_max
                    self%dcheby(k,n,d) = self%dcheby(k,n,d) + self%dcheby(k,i,d-1)*alpha(i,n)
                Enddo
            Enddo
            Enddo
            Enddo
        Endif
        DeAllocate(alpha)
    End Subroutine Gen_Tn_Deriv_Arrays

    Subroutine To_Spectral_1D(self,f_in,c_out)
        Implicit None
        Class(Cheby_Transform_Interface) :: self
        Real*8, Intent(In) :: f_in(:)
        Real*8, Intent(InOut) :: c_out(:)
        Real*8 :: alpha, beta
        Real*8, Allocatable :: f_even(:), f_odd(:), c_temp(:)
        Integer :: i, N_x, n_even, n_odd, n_max
        alpha = 2.0d0/n_max
        beta = 0.0d0

        n_x = self%n_x
        n_even = self%n_even
        n_odd = self%n_odd
        n_max = self%n_max

        If (self%parity) Then
            Allocate(c_temp(1:n_even ))
            Allocate(f_even(1:n_x ))
            Allocate( f_odd(1:n_x  ))
            Do i = 1, N_x
                f_even(i) = f_in(i)+f_in(N_max-i+1)
                f_odd(i)  = f_in(i)-f_in(N_max-i+1)
            Enddo
            CALL DGEMM('T','N',N_even,1,N_x, alpha, self%cheby_even, N_x,f_even , N_x, beta,c_temp,N_even)
            Do i = 1, n_even
                c_out(2*i-1) = c_temp(i)
            Enddo
            If (n_even .ne. n_odd) Then
                DeAllocate(c_temp)
                Allocate(c_temp(1:n_odd))
            Endif
            Call DGEMM('T','N',N_odd,1,N_x, alpha, self%cheby_odd, N_x,f_odd , N_x, beta,c_temp,N_odd)
            Do i = 1, N_odd
                c_out(2*i) = c_temp(i)
            Enddo
            DeAllocate(f_even,f_odd,c_temp)

        Else
            CALL DGEMM('T','N',N_max,1,N_Max, alpha, self%cheby, N_max,f_in , N_Max, beta,c_out,N_max)
        Endif
    End Subroutine To_Spectral_1D

    Subroutine To_Spectral_2D(self, f_in,c_out)
        Implicit None
        Real*8, Intent(In) :: f_in(:,:)
        Real*8, Intent(InOut) :: c_out(:,:)
        Real*8 :: alpha, beta
        Real*8, Allocatable :: f_even(:,:), f_odd(:,:), c_temp(:,:)
        Integer :: i, j, n2, dims(2), n_x, n_even, n_odd, n_max
        Class(Cheby_Transform_Interface) :: self
        n_x = self%n_x
        n_even = self%n_even
        n_odd = self%n_odd
        n_max = self%n_max

        alpha = 2.0d0/n_max
        beta = 0.0d0
        dims = shape(f_in)
        n2 = dims(2)

        If (self%parity) Then
            Allocate(c_temp(1:n_even, 1:n2 ))
            Allocate(f_even(1:n_x   , 1:n2 ))
            Allocate( f_odd(1:n_x   , 1:n2 ))
            Do j = 1, n2
            Do i = 1, N_x
                f_even(i,j) = f_in(i,j)+f_in(N_max-i+1,j)
                f_odd(i,j)  = f_in(i,j)-f_in(N_max-i+1,j)
            Enddo
            Enddo
            CALL DGEMM('T','N',N_even,n2,N_x, alpha, self%cheby_even, N_x,f_even , N_x, beta,c_temp,N_even)
            Do j = 1, n2
            Do i = 1, n_even
                c_out(2*i-1,j) = c_temp(i,j)
            Enddo
            Enddo
            If (n_even .ne. n_odd) Then
                DeAllocate(c_temp)
                Allocate(c_temp(1:n_odd,1:n2))
            Endif
            Call DGEMM('T','N',N_odd,n2,N_x, alpha, self%cheby_odd, N_x,f_odd , N_x, beta,c_temp,N_odd)
            Do j = 1, n2
            Do i = 1, N_odd
                c_out(2*i,j) = c_temp(i,j)
            Enddo
            Enddo
            DeAllocate(f_even,f_odd,c_temp)

        Else
            CALL DGEMM('T','N',N_max,n2,N_Max, alpha, self%cheby, N_max,f_in , N_Max, beta,c_out,N_max)
        Endif
    End Subroutine To_Spectral_2D

    Subroutine To_Spectral_3D(self, f_in,c_out)
        Implicit None
        ! Play a sneaky trick on dgemm and see if it sticks without complaining
        Real*8, Intent(In) :: f_in(:,:,:)
        Real*8, Intent(InOut) :: c_out(:,:,:)
        Real*8 :: alpha, beta
        Real*8, Allocatable :: f_even(:,:,:), f_odd(:,:,:), c_temp(:,:,:)
        Integer :: i, j, k, n2, n3, dims(3)
        Integer :: n_x, n_even, n_odd, n_max
        Class(Cheby_Transform_Interface) :: self
        n_x = self%n_x
        n_even = self%n_even
        n_odd = self%n_odd
        n_max = self%n_max

        alpha = 2.0d0/n_max
        beta = 0.0d0
        dims = shape(f_in)
        n2 = dims(2)
        n3 = dims(3)
        If (self%parity) Then
            Allocate(c_temp(1:n_even, 1:n2, 1:n3 ))
            Allocate(f_even(1:n_x   , 1:n2, 1:n3 ))
            Allocate( f_odd(1:n_x   , 1:n2, 1:n3 ))
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, N_x
                f_even(i,j,k) = f_in(i,j,k)+f_in(N_max-i+1,j,k)
                f_odd(i,j,k)  = f_in(i,j,k)-f_in(N_max-i+1,j,k)
            Enddo
            Enddo
            Enddo
            CALL DGEMM('T','N',N_even,n2*n3,N_x, alpha, self%cheby_even, N_x,f_even , N_x, beta,c_temp,N_even)
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_even
                c_out(2*i-1,j,k) = c_temp(i,j,k)
            Enddo
            Enddo
            Enddo
            If (n_even .ne. n_odd) Then
                DeAllocate(c_temp)
                Allocate(c_temp(1:n_odd,1:n2,1:n3))
            Endif
            Call DGEMM('T','N',N_odd,n2*n3,N_x, alpha, self%cheby_odd, N_x,f_odd , N_x, beta,c_temp,N_odd)
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, N_odd
                c_out(2*i,j,k) = c_temp(i,j,k)
            Enddo
            Enddo
            Enddo
            DeAllocate(f_even,f_odd,c_temp)

        Else
            CALL DGEMM('T','N',N_max,n2*n3,N_Max, alpha, self%cheby, N_max,f_in , N_Max, beta,c_out,N_max)
        Endif
    End Subroutine To_Spectral_3D

    Subroutine To_Spectral_4D(self,f_in,c_out)
        Implicit None
        ! Play a sneaky trick on dgemm and see if it sticks without complaining
        Real*8, Intent(In) :: f_in(:,:,:,:)
        Real*8, Intent(InOut) :: c_out(:,:,:,:)
        Real*8 :: alpha, beta
        Real*8, Allocatable :: f_even(:,:,:,:), f_odd(:,:,:,:), c_temp(:,:,:,:)
        Integer :: i, j, k, kk, n2, n3, n4, dims(4)
        Integer :: n_x, n_even, n_odd, n_max
        Class(Cheby_Transform_Interface) :: self
        n_x = self%n_x
        n_even = self%n_even
        n_odd = self%n_odd
        n_max = self%n_max

        alpha = 2.0d0/n_max
        beta = 0.0d0
        dims = shape(f_in)
        n2 = dims(2)
        n3 = dims(3)
        n4 = dims(4)
        If (self%parity) Then
            Allocate(c_temp(1:n_even, 1:n2, 1:n3, 1:n4 ))
            Allocate(f_even(1:n_x   , 1:n2, 1:n3, 1:n4 ))
            Allocate( f_odd(1:n_x   , 1:n2, 1:n3, 1:n4 ))
            Do kk= 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, N_x
                f_even(i,j,k, kk) = f_in(i,j,k,kk)+f_in(N_max-i+1,j,k,kk)
                f_odd(i,j,k, kk)  = f_in(i,j,k,kk)-f_in(N_max-i+1,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo
            CALL DGEMM('T','N',N_even,n2*n3*n4,N_x, alpha, self%cheby_even, N_x,f_even , N_x, beta,c_temp,N_even)
            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_even
                c_out(2*i-1,j,k,kk) = c_temp(i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo
            If (n_even .ne. n_odd) Then
                DeAllocate(c_temp)
                Allocate(c_temp(1:n_odd,1:n2,1:n3,1:n4))
            Endif
            Call DGEMM('T','N',N_odd,n2*n3*n4,N_x, alpha, self%cheby_odd, N_x,f_odd , N_x, beta,c_temp,N_odd)
            Do kk =1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, N_odd
                c_out(2*i,j,k,kk) = c_temp(i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo
            DeAllocate(f_even,f_odd,c_temp)

        Else
            CALL DGEMM('T','N',N_max,n2*n3*n4,N_Max, alpha, self%cheby, N_max,f_in , N_Max, beta,c_out,N_max)
        Endif
    End Subroutine To_Spectral_4D

    Subroutine From_Spectral_1D(self, c_in,f_out)
        Implicit None
        Real*8, Intent(In) :: c_in(:)
        Real*8, Intent(InOut) :: f_out(:)
        Real*8, Allocatable :: c_temp(:), f_temp(:)
        Real*8 :: alpha, beta
        Integer :: i
        Integer :: n_x, n_even, n_odd, n_max
        Class(Cheby_Transform_Interface) :: self
        n_x = self%n_x
        n_even = self%n_even
        n_odd = self%n_odd
        n_max = self%n_max

        alpha = 1.0d0
        beta = 0.0d0



        If (self%parity) Then
            Allocate(c_temp(1:n_even))
            Allocate(f_temp(1:n_x))
            Do i = 1, n_even
                c_temp(i) = c_in(2*i-1)
            Enddo
            CALL DGEMM('N','N',N_x,1,N_even, alpha, self%cheby_even, N_x,c_temp , N_even, beta,f_temp,N_x)
            Do i = 1, n_even
                f_out(i) = f_temp(i)
                f_out(N_max-i+1) = f_temp(i)
            Enddo
            If (n_even .ne. n_odd) Then
                f_out(n_x) = f_out(n_x)*2.0d0
                DeAllocate(c_temp)
                Allocate(c_temp(1:n_odd))
            Endif
            Do i = 1, n_odd
                c_temp(i) = c_in(2*i)
            Enddo
            CALL DGEMM('N','N',N_x,1,N_odd, alpha, self%cheby_odd, N_x,c_temp , N_odd, beta,f_temp,N_x)
            Do i = 1, n_odd
                f_out(i) = f_out(i) + f_temp(i)
                f_out(N_max-i+1) = f_out(N_max-i+1)-f_temp(i)
            Enddo
            If (n_even .ne. n_odd) Then
                f_out(n_x) = f_out(n_x) + f_temp(n_x)*2.0d0
            Endif
            DeAllocate(c_temp, f_temp)
        Else
            CALL DGEMM('N','N',N_max,1,N_Max, alpha, self%cheby, N_max,c_in , N_Max, beta,f_out,N_max)

        Endif
        f_out(:) = f_out(:) -c_in(1)/2.0d0
    End Subroutine From_Spectral_1D

    Subroutine From_Spectral_2D(self, c_in,f_out)
        Implicit None
        Real*8, Intent(In) :: c_in(:,:)
        Real*8, Intent(InOut) :: f_out(:,:)
        Real*8, Allocatable :: c_temp(:,:), f_temp(:,:)
        Real*8 :: alpha, beta
        Integer :: i, j, n2, dims(2)
        Integer :: n_x, n_even, n_odd, n_max
        Class(Cheby_Transform_Interface) :: self
        n_x = self%n_x
        n_even = self%n_even
        n_odd = self%n_odd
        n_max = self%n_max

        alpha = 1.0d0
        beta = 0.0d0

        dims = shape(c_in)
        n2 = dims(2)

        If (self%parity) Then
            Allocate(c_temp(1:n_even,1:n2))
            Allocate(f_temp(1:n_x,1:n2))
            Do j = 1, n2
            Do i = 1, n_even
                c_temp(i,j) = c_in(2*i-1,j)
            Enddo
            enddo

            CALL DGEMM('N','N',N_x,n2,N_even, alpha, self%cheby_even, N_x,c_temp , N_even, beta,f_temp,N_x)
            Do j = 1, n2
            Do i = 1, n_even
                f_out(i,j) = f_temp(i,j)
                f_out(N_max-i+1,j) = f_temp(i,j)
            Enddo
            Enddo
            If (n_even .ne. n_odd) Then
                Do j = 1, n2
                    f_out(n_x,j) = f_out(n_x,j)*2.0d0
                Enddo
                DeAllocate(c_temp)
                Allocate(c_temp(1:n_odd,1:n2))
            Endif
            Do j = 1, n2
            Do i = 1, n_odd
                c_temp(i,j) = c_in(2*i,j)
            Enddo
            Enddo
            CALL DGEMM('N','N',N_x,n2,N_odd, alpha, self%cheby_odd, N_x,c_temp , N_odd, beta,f_temp,N_x)

            Do j = 1, n2
            Do i = 1, n_odd
                f_out(i,j) = f_out(i,j) + f_temp(i,j)
                f_out(N_max-i+1,j) = f_out(N_max-i+1,j)-f_temp(i,j)
            Enddo
            Enddo
            If (n_even .ne. n_odd) Then
                f_out(n_x,:) = f_out(n_x,:) + f_temp(n_x,:)*2.0d0
            Endif
            DeAllocate(c_temp, f_temp)
        Else
            CALL DGEMM('N','N',N_max,n2,N_Max, alpha, self%cheby, N_max,c_in , N_Max, beta,f_out,N_max)

        Endif
        Do j = 1, n2
            f_out(:,j) = f_out(:,j) -c_in(1,j)/2.0d0
        Enddo
    End Subroutine From_Spectral_2D

    Subroutine From_Spectral_3D(self, c_in,f_out)
        Implicit None
        Real*8, Intent(In) :: c_in(:,:,:)
        Real*8, Intent(InOut) :: f_out(:,:,:)
        Real*8, Allocatable :: c_temp(:,:,:), f_temp(:,:,:)
        Real*8 :: alpha, beta
        Integer :: i, j, k, n2, n3, dims(3)
        Integer :: n_x, n_even, n_odd, n_max
        Class(Cheby_Transform_Interface) :: self
        n_x = self%n_x
        n_even = self%n_even
        n_odd = self%n_odd
        n_max = self%n_max

        alpha = 1.0d0
        beta = 0.0d0

        dims = shape(c_in)
        n2 = dims(2)
        n3 = dims(3)
        If (self%parity) Then
            Allocate(c_temp(1:n_even,1:n2, 1:n3))
            Allocate(f_temp(1:n_x   ,1:n2, 1:n3))
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_even
                c_temp(i,j,k) = c_in(2*i-1,j,k)
            Enddo
            Enddo
            Enddo

            CALL DGEMM('N','N',N_x,n2*n3,N_even, alpha, self%cheby_even, N_x,c_temp , N_even, beta,f_temp,N_x)
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_even
                f_out(i,j,k) = f_temp(i,j,k)
                f_out(N_max-i+1,j,k) = f_temp(i,j,k)
            Enddo
            Enddo
            Enddo

            If (n_even .ne. n_odd) Then
                Do k = 1, n3
                Do j = 1, n2
                    f_out(n_x,j,k) = f_out(n_x,j,k)*2.0d0
                Enddo
                Enddo
                DeAllocate(c_temp)
                Allocate(c_temp(1:n_odd,1:n2,1:n3))
            Endif
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_odd
                c_temp(i,j,k) = c_in(2*i,j,k)
            Enddo
            Enddo
            Enddo
            CALL DGEMM('N','N',N_x,n2*n3,N_odd, alpha, self%cheby_odd, N_x,c_temp , N_odd, beta,f_temp,N_x)

            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_odd
                f_out(i,j,k) = f_out(i,j,k) + f_temp(i,j,k)
                f_out(N_max-i+1,j,k) = f_out(N_max-i+1,j,k)-f_temp(i,j,k)
            Enddo
            Enddo
            Enddo
            If (n_even .ne. n_odd) Then
                f_out(n_x,:,:) = f_out(n_x,:,:) + f_temp(n_x,:,:)*2.0d0
            Endif
            DeAllocate(c_temp, f_temp)
        Else
            CALL DGEMM('N','N',N_max,n2*n3,N_Max, alpha, self%cheby, N_max,c_in , N_Max, beta,f_out,N_max)

        Endif
        Do k = 1, n3
        Do j = 1, n2
            f_out(:,j,k) = f_out(:,j,k) -c_in(1,j,k)/2.0d0
        Enddo
        Enddo
    End Subroutine From_Spectral_3D

    Subroutine From_Spectral_4D(self,c_in,f_out)
        Implicit None
        Real*8, Intent(In) :: c_in(:,:,:,:)
        Real*8, Intent(InOut) :: f_out(:,:,:,:)
        Real*8, Allocatable :: c_temp(:,:,:,:), f_temp(:,:,:,:)
        Real*8 :: alpha, beta
        Integer :: i, j, k, kk, n2, n3, n4, dims(4)
        Integer :: n_x, n_even, n_odd, n_max
        Class(Cheby_Transform_Interface) :: self
        n_x = self%n_x
        n_even = self%n_even
        n_odd = self%n_odd
        n_max = self%n_max

        alpha = 1.0d0
        beta = 0.0d0

        dims = shape(c_in)
        n2 = dims(2)
        n3 = dims(3)
        n4 = dims(4)
        If (self%parity) Then
            Allocate(c_temp(1:n_even,1:n2, 1:n3, 1:n4))
            Allocate(f_temp(1:n_x   ,1:n2, 1:n3, 1:n4))
            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_even
                c_temp(i,j,k,kk) = c_in(2*i-1,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo

            CALL DGEMM('N','N',N_x,n2*n3*n4,N_even, alpha, self%cheby_even, N_x,c_temp , N_even, beta,f_temp,N_x)

            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_even
                f_out(i,j,k,kk) = f_temp(i,j,k,kk)
                f_out(N_max-i+1,j,k,kk) = f_temp(i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo

            If (n_even .ne. n_odd) Then
                Do kk = 1, n4
                Do k = 1, n3
                Do j = 1, n2
                    f_out(n_x,j,k,kk) = f_out(n_x,j,k,kk)*2.0d0
                Enddo
                Enddo
                Enddo
                DeAllocate(c_temp)
                Allocate(c_temp(1:n_odd,1:n2,1:n3,1:n4))
            Endif
            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_odd
                c_temp(i,j,k,kk) = c_in(2*i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo
            CALL DGEMM('N','N',N_x,n2*n3*n4,N_odd, alpha, self%cheby_odd, N_x,c_temp , N_odd, beta,f_temp,N_x)

            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_odd
                f_out(i,j,k,kk) = f_out(i,j,k,kk) + f_temp(i,j,k,kk)
                f_out(N_max-i+1,j,k,kk) = f_out(N_max-i+1,j,k,kk)-f_temp(i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo
            If (n_even .ne. n_odd) Then
                f_out(n_x,:,:,:) = f_out(n_x,:,:,:) + f_temp(n_x,:,:,:)*2.0d0
            Endif
            DeAllocate(c_temp, f_temp)
        Else
            CALL DGEMM('N','N',N_max,n2*n3*n4,N_Max, alpha, self%cheby, N_max,c_in , N_Max, beta,f_out,N_max)

        Endif
        Do kk = 1, n4
        Do k = 1, n3
        Do j = 1, n2
            f_out(:,j,k,kk) = f_out(:,j,k,kk) -c_in(1,j,k,kk)/2.0d0
        Enddo
        Enddo
        Enddo
    End Subroutine From_Spectral_4D

    Subroutine Cheby_Deriv_Buffer_4Dalt(self,ind,dind,buffer,dorder)
        ! This is exactly the same procedure as in the normal chebyshev polynomials module
        ! I just changed the name for now to avoid any conflicts
        Implicit None
        Real*8,  Intent(InOut) :: buffer(0:,1:,1:,1:)    ! Makes it easier to reconcile with my IDL code
        Integer, Intent(In)    :: ind, dind, dorder
        Real*8, Allocatable :: dbuffer(:,:)
        Integer :: dims(4), n,n2,n3, i,j,k, order
        Real*8 :: scaling
        Class(Cheby_Transform_Interface) :: self
        scaling = self%scaling
        dims = shape(buffer)
        n = dims(1)
        n2 = dims(2)
        n3 = dims(3)
        If (ind .ne. dind) Then
            Do k = 1, n3
                Do j = 1, n2
                    buffer(n-1,j,k,dind) = 0.0d0
                    buffer(n-2,j,k,dind) = 2.0d0*(n-1)*buffer(n-1,j,k,ind)*scaling
                    Do i = n-3,0, -1
                        buffer(i,j,k,dind) = buffer(i+2,j,k,dind)+2.0d0*(i+1)*buffer(i+1,j,k,ind)*scaling
                    Enddo
                Enddo
            Enddo
            If (dorder .gt. 1) Then
                Allocate(dbuffer(0:n-1,1:dorder))
                Do k = 1, n3
                    Do j = 1, n2
                        dbuffer(:,1) = buffer(:,j,k,dind)
                        Do order = 2, dorder
                            dbuffer(n-1,order) = 0.0d0
                            dbuffer(n-2,order) = 2.0d0*(n-1)*dbuffer(n-1,order-1)*scaling
                            Do i = n -3, 0, -1
                                dbuffer(i,order) = dbuffer(i+2,order)+2.0d0*(i+1)*dbuffer(i+1,order-1)*scaling
                            Enddo
                        Enddo
                        buffer(:,j,k,dind) = dbuffer(:,dorder)
                    Enddo
                Enddo

                DeAllocate(dbuffer)
            Endif
        Else
            ! In-place
        Endif
        buffer(:,:,:,dind) = buffer(:,:,:,dind)
    End Subroutine Cheby_Deriv_Buffer_4Dalt





End Module Chebyshev_Polynomials_Alt
