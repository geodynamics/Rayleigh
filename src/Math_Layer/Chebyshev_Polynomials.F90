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

Module Chebyshev_Polynomials
    Use Structures
    ! Module for computing Chebyshev Polynomial Arrays and the associated Derivative Arrays
    Implicit None
    Integer :: cp_nthreads
    Real*8, private ::    Pi  = 3.1415926535897932384626433832795028841972d+0
    Real*8, private :: pi_over_N
    Logical :: DeAlias = .true.
    Logical :: Parity = .true.
    Logical :: initialized = .false.
    Logical, Private :: use_extrema = .false. ! Set to true to use extrema points instead of zeros for T_Nmax

    Type, Public :: Cheby_Grid
        Integer :: domain_count = 0  ! Number of subdomains for this grid object
        Integer ::    max_npoly = 0  ! maximum number of (aliased) polynomials within any domain
        Integer ::       ntotal = 0  ! Number of colocation points for GLOBAL grid

        Integer, Allocatable :: npoly(:)  ! Number of polynomials within each subdomain
        Integer, Allocatable :: n_x(:)    ! Number of colocation points within each subdomain
        Integer, Allocatable :: rda(:)    ! Deliasing index for each subdomain
        Real*8, Allocatable  :: domain_bounds(:) ! The bounds that divide each subdomain
        Real*8, Allocatable  :: n_even(:), n_odd(:) ! # of even and odd polynomials within each domain
        Real*8, Allocatable  :: x(:,:), theta(:,:) ! colocation points and angle for each subdomain
        Real*8, Allocatable  :: deriv_scaling(:)

        Type(rmcontainer)  , Allocatable :: cheby_even(:), cheby_odd(:)
        Type(rmcontainer3d), Allocatable :: dcheby(:)

        Contains

            Procedure :: Init => Initialize_Cheby_Grid
            Procedure :: gen_Tn
            Procedure :: gen_colocation_points
            Procedure :: gen_Tn_Deriv_arrays
            Procedure :: dealias_buffer => dealias_buffer4d
            Procedure :: to_spectral => to_spectral4d
            Procedure :: from_spectral => from_spectral4d
            Procedure :: d_by_dr_cp => Cheby_Deriv_Buffer_4D
    End Type Cheby_Grid

    Type(Cheby_Grid), Public :: main_grid  ! Publically accessible grid object (alleviates need for pointers)

Contains

    Subroutine Initialize_Cheby_Grid(self,grid,integration_weights,ndomains,npoly, &
                & bounds,dmax,nthread,dealias_by, verbose)
        Implicit None
        Class(Cheby_Grid) :: self
        Real*8, Intent(InOut) :: grid(1:), integration_weights(1:)
        Real*8, Intent(In) :: bounds(1:)
        Integer, Intent(In) :: ndomains, npoly(1:)
        Integer, Intent(In), Optional :: nthread, dmax
        Integer, Intent(In), Optional :: dealias_by(1:)
        Logical, Intent(In), Optional :: verbose

        Real*8 :: domain_delta, scaling, upperb, lowerb, xmin,xmax
        Real*8 :: gmax, gmin, xx,int_scale
        Integer :: r, i,n, domain_count
        Integer :: ind, ind2, n_max, dmx, gindex, db,scheck
        Logical :: custom_dealiasing = .false.
        Logical :: report
        report = .false.

        If (present(verbose)) Then
            If (verbose) report = .true.
        Endif

        If (report) Then
            write(6,*)'Ndomains: ', ndomains
            write(6,*)'npoly   : ', npoly(1:ndomains)
            write(6,*)'bounds  : ', bounds(1:ndomains+1)
        Endif
        If (present(dealias_by)) Then
            scheck = size(dealias_by)
            If ( (scheck .ge. ndomains) .and. report ) Then
                Write(6,*)'dealias_by: ', dealias_by(1:ndomains)
            Endif
        Endif

        !Note that bounds and npoly are assumed to be provided in ascending order
        !We reverse them so that the subdomains agree with a globally reversed grid
        dmx = 3
        If (present(dmax)) Then
            If (dmax .ge. 1) dmx = dmax
        Endif


        If (present(dealias_by)) Then
            scheck = size(dealias_by)
            if (scheck .ge. ndomains) custom_dealiasing = .true.
        Endif

        Allocate(self%cheby_even(1:ndomains))
        Allocate(self%cheby_odd(1:ndomains))
        Allocate(self%dcheby(1:ndomains))
        Allocate(self%n_odd(1:ndomains))
        Allocate(self%n_even(1:ndomains))


        domain_count = ndomains
        self%domain_count = ndomains
        Allocate(self%domain_bounds(1:domain_count+1))
        Allocate(self%npoly(1:domain_count))
        Allocate(self%rda(1:domain_count))
        Do i = 1, domain_count
            n = npoly(domain_count-i+1)
            self%npoly(i)  =  n
            db = (2*n)/3+1
            self%rda(i) = db
            If (custom_dealiasing) Then
                db = dealias_by(domain_count+1-i)

                If ((db .ge. 1) .and. (db .lt. n) ) Then
                    self%rda(i) = n-db+1
                Endif

            Endif
        Enddo

        Do i = 1, domain_count+1
            self%domain_bounds(i) = bounds(domain_count-i+2)
        Enddo

        gmax = self%domain_bounds(1)
        gmin = self%domain_bounds(domain_count+1)

        self%max_npoly = maxval(self%npoly)
        self%ntotal = sum(self%npoly)
        Allocate(self%x(1:self%max_npoly, 1:domain_count))
        Allocate(self%theta(1:self%max_npoly,1:domain_count))
        Allocate(self%deriv_scaling(1:self%ntotal))

        Call self%gen_colocation_points()
        Call self%gen_Tn()
        Call self%gen_Tn_Deriv_Arrays(dmx)

        !Compute the global grid and rescale the chebyshev derivative arrays
        ind = 1
        integration_weights(:) = 0.0d0
        grid(:) = 0.0d0
        Do n = 1, self%domain_count
            n_max = self%npoly(n)
            ind2 = ind+self%npoly(n)-1

            upperb = self%domain_bounds(n)
            lowerb = self%domain_bounds(n+1)
            xmin = self%x(n_max,n)
            xmax = self%x(1,n)

            domain_delta = upperb-lowerb

            scaling = domain_delta/(xmax-xmin)
            grid(ind:ind2) = (self%x(1:n_max,n) -xmin)*scaling +lowerb


            int_scale = (3*Pi* scaling)/( (gmax**3 - gmin**3) * n_max )

            scaling = 1.0d0/scaling


            Do i = 1, dmx
                self%dcheby(n)%data(:,:,i) = &
                    & self%dcheby(n)%data(:,:,i)*(scaling**i)
            Enddo

            Do i=1,N_max
                gindex = ind+i-1
                 xx = self%x(i,n)
                integration_weights(gindex) = &
                    & int_scale*grid(gindex)**2  * sqrt(1.0d0-xx*xx)
                self%deriv_scaling(gindex) = scaling
            Enddo
            integration_weights(ind) = integration_weights(ind)*0.5d0      !Boundaries x 1/2 (since on zero points)
            integration_weights(ind2) = integration_weights(ind2)*0.5d0



            ind = ind2+1
        Enddo

        cp_nthreads = 1
        If (present(nthread)) Then
            If (report) write(6,*)'nthreads: ', nthread
            If (nthread .gt. 1) cp_nthreads = nthread
        Endif

    End Subroutine Initialize_Cheby_Grid


    Subroutine Gen_Colocation_Points(self)
        Implicit None
        Class(Cheby_Grid) :: self
        Integer           :: i,n,n_max
        Real*8            :: arg, dtheta, theta0
        ! Calculate the colocation points X { -1 , 1}
        ! Also calculate the theta grid theta { 0 , pi }

        Do n = 1,self%domain_count
            n_max = self%npoly(n)

            If (use_extrema) Then  !Use the extrema of T_{N_max-1}

                dtheta = pi/(N_max-1)
                theta0 = 0.0d0
                arg = theta0

                Do i = 1, N_max
                    self%theta(i,n) = arg
                    self%x(i,n) = cos(arg)
                    arg = arg+dtheta
                Enddo

            Else !Use the zeroes of T_{N_max}

                dtheta = pi/N_max
                theta0 = dtheta*0.5d0
                arg = theta0

                Do i = 1, N_Max
                    self%theta(i,n) = arg
                    self%x(i,n) = cos(arg)
                    arg = arg+dtheta
                Enddo

            Endif

        Enddo
    End Subroutine Gen_Colocation_Points


    Subroutine Gen_Tn(self)
        Implicit None
        Class (Cheby_Grid)  :: self
        Integer             :: i, k,n,n_even, n_max, n_odd, n_x, r
        Real*8              :: arg
        Real*8, Allocatable :: cheby(:,:) ! workspace
        Allocate(cheby(1:self%max_npoly,1:self%max_npoly))
        Do n = 1,self%domain_count
            n_max = self%npoly(n)
            Do r = 1, N_max
                Do i = 1, N_max
                    k = i -1
                    arg = k*self%theta(r,n)
                    cheby(r,i) = cos(arg)
                Enddo
            Enddo

            n_odd = N_max/2
            n_even = n_odd+mod(N_max,2)
            n_x = n_even
            self%n_odd(n) = n_odd
            self%n_even(n) = n_even

            Allocate(self%cheby_even(n)%data(1:N_x,1:N_even))
            Allocate(self%cheby_odd( n)%data(1:N_x,1:N_odd))
            Do i = 1, n_even
                self%cheby_even(n)%data(1:N_x,i) = &
                    & cheby(1:N_x,2*i-1)
            Enddo
            Do i = 1,n_odd
                self%cheby_odd(n)%data(1:N_x,i) = &
                    & cheby(1:N_x,2*i)
            Enddo
            If (n_x .ne. n_odd) Then
                ! We actually have an x = 0 point
                ! We must be careful not to double count the power here
                ! when exploiting parity to speed up the transforms.
                ! No such adjustment need be made for the regular cheby array
                self%cheby_even(n)%data(n_x,:) = &
                    & 0.5d0*self%cheby_even(n)%data(n_x,:)
                self%cheby_odd(n )%data(n_x,:) = &
                    & 0.5d0*self%cheby_odd( n)%data(n_x,:)
                ! There is really no need to modify cheby_odd, but this way it is consistently stored.
            Endif
        Enddo
        DeAllocate(cheby)
    End Subroutine Gen_Tn

    Subroutine Gen_Tn_Deriv_Arrays(self,dmax)
        Implicit None
        Class (Cheby_Grid)  :: self
        Integer, Intent(In) :: dmax
        Integer :: d,i, k, m, n, n_max ,r
        Real*8 :: arg
        Real*8, Allocatable :: alpha(:,:)

        ! sum_n (alpha_kn  c_n) = c'_k
        Allocate(alpha(1:self%max_npoly,self%max_npoly))

        Do m = 1,self%domain_count

            N_max = self%npoly(m)

            alpha(:,:) = 0.0d0
            alpha(N_max,:) = 0.0d0
            alpha(N_max-1,N_max) = 2.0d0*(N_max-1)

            Do k = N_max-2, 1, -1
                alpha(k,k+1) = 2.0d0*k
                alpha(k,:) = alpha(k,:)+alpha(k+2,:)
            Enddo

            Allocate(self%dcheby(m)%data(1:N_max,1:N_max,0:dmax))
            self%dcheby(m)%data(:,:,:) = 0.0d0

            Do r = 1, N_max
                Do i = 1, N_max
                    k = i -1
                    arg = k*self%theta(r,m)
                    self%dcheby(m)%data(r,i,0) = cos(arg)
                Enddo
            Enddo


            self%dcheby(m)%data(:,1,0) = &
                self%dcheby(m)%data(:,1,0) - 0.5d0    ! This accounts for the -1/2c_0 in going from_spectral
            If (dmax .ge. 1) Then
                Do d = 1, dmax
                Do n = 1, N_max
                Do k = 1, N_max
                Do i = 1, N_max
                    self%dcheby(m)%data(k,n,d) = self%dcheby(m)%data(k,n,d) &
                        & + self%dcheby(m)%data(k,i,d-1)*alpha(i,n)
                Enddo
                Enddo
                Enddo
                Enddo
            Endif
        Enddo
        DeAllocate(alpha)
    End Subroutine Gen_Tn_Deriv_Arrays

    Subroutine dealias_buffer4d(self, buffer)
        Implicit None
        Class(Cheby_Grid) :: self
        Real*8, Intent(InOut) :: buffer(1:,1:,1:,1:)
        Integer :: bsize(4), i,j,k,n,offset, ind1,ind2

        bsize = shape(buffer)
        !do n = 1, self%domain_count
        !    WRite(6,*)'da check: ', self%rda(n), self%npoly(n)
        !Enddo
        Do k = 1, bsize(4)
        Do j = 1, bsize(3)
        Do i = 1, bsize(2)
        offset = 0
        Do n = 1, self%domain_count
            ind1 = self%rda(n)+offset
            ind2 = self%npoly(n)+offset
            buffer(ind1:ind2,i,j,k) = 0.0d0
            offset = offset+self%npoly(n)
        Enddo
        Enddo
        Enddo
        Enddo
    End Subroutine dealias_buffer4d
    Subroutine To_Spectral4D(self,f_in,c_out)
        Implicit None
        Class(Cheby_Grid) :: self
        Real*8, Intent(In) :: f_in(:,:,:,:)
        Real*8, Intent(InOut) :: c_out(:,:,:,:)
        Real*8 :: alpha, beta
        Real*8, Allocatable :: f_even(:,:,:,:), f_odd(:,:,:,:), c_temp(:,:,:,:)
        Integer :: i, j, k, kk, n2, n3, n4, dims(4),nsub,nglobal, hoff, hh
        Integer :: istart, iend, n_even, n_odd, n_max, n_x


        beta = 0.0d0
        dims = shape(f_in)
        nglobal = dims(1)
        n2 = dims(2)
        n3 = dims(3)
        n4 = dims(4)
        nsub = self%domain_count
        hoff = 0
        DO hh = 1, nsub
            n_even = self%n_even(hh)
            n_odd = self%n_odd(hh)
            n_max = self%npoly(hh)
            n_x = self%n_even(hh)
            alpha = 2.0d0/n_max
            Allocate(c_temp(1:n_even, 1:n2, 1:n3, 1:n4 ))
            Allocate(f_even(1:n_x   , 1:n2, 1:n3, 1:n4 ))
            Allocate( f_odd(1:n_x   , 1:n2, 1:n3, 1:n4 ))

            Do kk= 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, N_x
                f_even(i,j,k, kk) = f_in(hoff+i,j,k,kk)+f_in(hoff+N_max-i+1,j,k,kk)
                f_odd(i,j,k, kk)  = f_in(hoff+i,j,k,kk)-f_in(hoff+N_max-i+1,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo

            CALL DGEMM('T','N',N_even,n2*n3*n4,N_x, alpha, &
                & self%cheby_even(hh)%data, N_x,f_even , N_x, beta,c_temp,N_even)

            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_even
                c_out(hoff+2*i-1,j,k,kk) = c_temp(i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo

            If (n_even .ne. n_odd) Then
                DeAllocate(c_temp)
                Allocate(c_temp(1:n_odd,1:n2,1:n3,1:n4))
            Endif
            Call DGEMM('T','N',N_odd,n2*n3*n4,N_x, alpha, &
                & self%cheby_odd(hh)%data, N_x,f_odd , N_x, beta,c_temp,N_odd)

            Do kk =1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, N_odd
                c_out(hoff+2*i,j,k,kk) = c_temp(i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo
            DeAllocate(f_even,f_odd,c_temp)
            hoff = hoff+self%npoly(hh)

        ENDDO
    End Subroutine To_Spectral4D

    Subroutine From_Spectral4D(self,c_in,f_out)
        Implicit None
        ! Parity is assumed to be true for NOW
        ! This may be inefficient, but I just want it to work first.  Will optimize later
        Class(Cheby_Grid) :: self
        Real*8, Intent(In) :: c_in(:,:,:,:)
        Real*8, Intent(InOut) :: f_out(:,:,:,:)
        Real*8, Allocatable :: c_temp(:,:,:,:), f_temp(:,:,:,:)
        Real*8 :: alpha, beta
        Integer :: i, j, k, kk, n2, n3, n4, dims(4),nsub,nglobal, hoff, hh
        Integer :: istart, iend, n_even, n_odd, n_max, n_x
        alpha = 1.0d0
        beta = 0.0d0

        dims = shape(c_in)
        nglobal = dims(1)
        n2 = dims(2)
        n3 = dims(3)
        n4 = dims(4)
        nsub = self%domain_count

        hoff = 0
        DO hh = 1, nsub
            n_even = self%n_even(hh)
            n_odd = self%n_odd(hh)
            n_max = self%npoly(hh)
            n_x = self%n_even(hh)

            Allocate(c_temp(1:n_even,1:n2, 1:n3, 1:n4))
            Allocate(f_temp(1:n_x   ,1:n2, 1:n3, 1:n4))
            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_even
                c_temp(i,j,k,kk) = c_in(hoff+2*i-1,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo

            CALL DGEMM('N','N',N_x,n2*n3*n4,N_even, alpha, &
                & self%cheby_even(hh)%data, N_x,c_temp , N_even, beta,f_temp,N_x)

            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_even
                f_out(i+hoff,j,k,kk) = f_temp(i,j,k,kk)
                f_out(hoff+N_max-i+1,j,k,kk) = f_temp(i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo

            If (n_even .ne. n_odd) Then
                Do kk = 1, n4
                Do k = 1, n3
                Do j = 1, n2
                    f_out(hoff+n_x,j,k,kk) = f_out(hoff+n_x,j,k,kk)*2.0d0
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
                c_temp(i,j,k,kk) = c_in(hoff+2*i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo
            CALL DGEMM('N','N',N_x,n2*n3*n4,N_odd, alpha, &
                & self%cheby_odd(hh)%data, N_x,c_temp , N_odd, beta,f_temp,N_x)

            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
            Do i = 1, n_odd
                f_out(hoff+i,j,k,kk) = f_out(hoff+i,j,k,kk) + f_temp(i,j,k,kk)
                f_out(hoff+N_max-i+1,j,k,kk) = f_out(hoff+N_max-i+1,j,k,kk)-f_temp(i,j,k,kk)
            Enddo
            Enddo
            Enddo
            Enddo
            If (n_even .ne. n_odd) Then
                f_out(hoff+n_x,:,:,:) = f_out(hoff+n_x,:,:,:) + f_temp(hoff+n_x,:,:,:)*2.0d0
            Endif
            DeAllocate(c_temp, f_temp)
            hoff = hoff+self%npoly(hh)
        ENDDO ! HH

        hoff = 0
        Do hh = 1, nsub
            n_max = self%npoly(hh)
            istart = 1+hoff
            iend = istart+n_max-1
            Do kk = 1, n4
            Do k = 1, n3
            Do j = 1, n2
                f_out(istart:iend,j,k,kk) = f_out(istart:iend,j,k,kk) -c_in(hoff+1,j,k,kk)/2.0d0
            Enddo
            Enddo
            Enddo
            hoff = hoff+self%npoly(hh)
        Enddo
    End Subroutine From_Spectral4D

    Subroutine Cheby_Deriv_Buffer_4D(self,ind,dind,buffer,dorder)
#ifdef useomp
        Use Omp_lib
#endif
        Implicit None
        Class (Cheby_Grid) :: self
        Real*8,  Intent(InOut) :: buffer(0:,1:,1:,1:)    ! Makes it easier to reconcile with my IDL code
        Integer, Intent(In)    :: ind, dind, dorder
        Real*8, Allocatable :: dbuffer(:,:,:)
        Integer :: dims(4), n,n2,n3, i,j,k, order
        Integer :: kstart, kend, nthr,trank
        Integer :: nglobal, nsub, hoff, hh
        dims = shape(buffer)
        nglobal = dims(1)
        nsub = self%domain_count

        n2 = dims(2)
        n3 = dims(3)
        If (ind .ne. dind) Then
        !$OMP PARALLEL DO PRIVATE(i,j,k,hh,hoff,n)
            Do k = 1, n3
                Do j = 1, n2
                hoff = 0
                DO hh = 1, nsub
                    !scaling = self%scaling(hh)
                    n = self%npoly(hh)

                    buffer(hoff+n-1,j,k,dind) = 0.0d0
                    buffer(hoff+n-2,j,k,dind) = 2.0d0*(n-1)*buffer(hoff+n-1,j,k,ind) !*scaling
                    Do i = n-3,0, -1
                        buffer(hoff+i,j,k,dind) = buffer(hoff+i+2,j,k,dind) &
                            & +2.0d0*(i+1)*buffer(hoff+i+1,j,k,ind) !*scaling
                    Enddo
                    hoff = hoff+self%npoly(hh)
                ENDDO !hh
                Enddo
            Enddo
         !$OMP END PARALLEL DO
            If (dorder .gt. 1) Then
                Allocate(dbuffer(0:nglobal-1,1:dorder,0:cp_nthreads-1))
            !$OMP PARALLEL PRIVATE(i,j,k,trank,order,kstart,kend,nthr,n,hh,hoff)
#ifdef useomp
                trank = omp_get_thread_num()
                  nthr  = omp_get_num_threads()
                  kstart = (trank*n3)/nthr+1
                  kend = ((trank+1)*n3)/nthr
#else
                trank = 0
                kstart = 1
                kend = n3
#endif
                Do k = kstart,kend

                    Do j = 1, n2
                        dbuffer(:,1,trank) = buffer(:,j,k,dind)
                        Do order = 2, dorder
                            hoff = 0
                            DO hh = 1, nsub
                                n = self%npoly(hh)

                            dbuffer(hoff+n-1,order,trank) = 0.0d0
                            dbuffer(hoff+n-2,order,trank) = 2.0d0*(n-1)*dbuffer(hoff+n-1,order-1,trank)
                            Do i = n -3, 0, -1
                                dbuffer(hoff+i,order,trank) = dbuffer(hoff+i+2,order,trank)+ &
                                    & 2.0d0*(i+1)*dbuffer(hoff+i+1,order-1,trank)
                            Enddo
                            hoff = hoff+self%npoly(hh)
                            ENDDO
                        Enddo
                        buffer(:,j,k,dind) = dbuffer(:,dorder,trank)
                    Enddo
                Enddo
            !$OMP END PARALLEL
                DeAllocate(dbuffer)
            Endif
        Else
            ! In-place -- Needs developing
        Endif
        Do k = 1, n3
        Do j = 1, n2
        buffer(:,j,k,dind) = buffer(:,j,k,dind)*(self%deriv_scaling(:)**dorder)
        Enddo
        Enddo
    End Subroutine Cheby_Deriv_Buffer_4D

End Module Chebyshev_Polynomials
