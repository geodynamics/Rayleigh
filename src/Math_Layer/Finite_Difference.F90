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

!///////////////////////////////////////////////////
!
!        Finite-Difference Derivative Module
!            (This is not compact finite-difference)
!            (Based on Taylor-exansion)
!
!//////////////////////////////////////////////////
Module Finite_Difference
    Use Math_Constants, Only : zero, one
    Implicit None
    Integer, Private :: n_x_fd
    Real*8, Private, Allocatable :: x_fd(:)
    Integer, Private, save :: N_s = 5  ! number of stencil points
    Integer, Private, save :: kd ! location of the diagonal in the stencil
    !Real*8 :: one=1.0d0
    !Real*8 :: zero = 0.0d0

    !***********************************************************************************
    !            Variables used for the new radial derivative coefficient generation scheme
    !                (there may be good reasons for these to be private... but public for now...)
    real*8, Public, Allocatable, Target :: d_coefs(:,:), dd_coefs(:,:), ddd_coefs(:,:)! Derivative coefficients
    Integer, Public :: n_stencil_1 = 5        ! stencil size for first order radial derivative - 4th order accurate by default
    Integer, Public :: n_stencil_2 = 5        ! stencil size for 2nd order radial derivative - 4th order accurate by default
    Integer, Public :: n_stencil_3 = 5        ! stencil size for 3rd order radial derivative - 2th order accurate by default
    Integer, Public :: n_stencil_2_nl = 9         ! stencil size for special 2nd order radial derivative
    Integer, Public :: boundary_accuracy_1 = 2        ! boundary accuracy for 1st order radial derivative
    Integer, Public :: boundary_accuracy_2 = 2        ! boundary accuracy for 2nd order radial derivative
    Integer, Public :: boundary_accuracy_3 = 2        ! boundary accuracy for 3rd order radial derivative
    Integer, Public :: boundary_accuracy_2_nl = 7        ! boundary accuracy for special 2nd order radial derivative - set so that full stencil is used (even if it has zeroes)
    Integer, Public :: derivative_transition_type = 1    ! Use asymmetric mesh for near-boundary points by default (set to 2 to use forward/backward differences)
    Integer, Public :: stencil_sizes(3), max_stencil_size
    !Namelist /Derivatives_Namelist/ n_stencil_1, n_stencil_2, n_stencil_3, boundary_accuracy_1, boundary_accuracy_2, boundary_accuracy_3
    !********************************************************

    Interface d_by_dx
        Module Procedure d_by_dx1d
        Module Procedure d_by_dx3d
        Module Procedure d_by_dx3d2
        Module Procedure d_by_dx3d3
    End Interface

Contains

    !=========================================================================
    !                         Section (I)  :   General Intialization Routines



    Subroutine Initialize_Derivatives(xin, integration_weights) ! also includes other geometrical quantities
        Real*8, Intent(In) :: xin(1:)
        Real*8, Intent(InOut) :: integration_weights(1:)
        Integer :: i
        Real*8 :: delr, int_sum
        n_x_fd = size(xin)
        Allocate(x_fd(1:n_x_fd))
        x_fd(1:n_x_fd) = xin(1:n_x_fd)

        Call Initialize_dcoefs()

        integration_weights(:) = 0.0d0
        Do i = 2, n_x_fd-1
            delr = (xin(i-1)-xin(i+1))/2.0d0
            integration_weights(i) = delr*xin(i)**2
        Enddo
        delr = ( xin(1)-xin(2) )/ 2.0d0
        integration_weights(1) = delr*xin(1)**2

        delr = (xin(n_x_fd-1)-xin(n_x_fd))/2.0d0
        integration_weights(n_x_fd) = delr*xin(n_x_fd)**2

        int_sum = sum(integration_weights)
        integration_weights = integration_weights/int_sum
    End Subroutine Initialize_Derivatives


    Subroutine Rescale_Grid_FD(length_scale)
        Implicit None
        Real*8, Intent(In) :: length_scale
        ! Following initialization, we can rescale the FD arrays if we choose
        ! This is useful when nondimensionalizing after the reference state has been set up
        ! (which typically requires a radial grid to have been established)
        x_fd(:) = x_fd(:)/length_scale
        d_coefs(:,:) = d_coefs(:,:)*length_scale
        dd_coefs(:,:) = dd_coefs(:,:)*(length_scale**2)
        ddd_coefs(:,:) = ddd_coefs(:,:)*(length_scale**3)
    End Subroutine Rescale_Grid_FD

    !=========================================================================
    !
    !      Section (III)  :   Finite Difference Derivative Routines
    !
    !=========================================================================


    !************************************************************************
    Subroutine d_by_dx1d(field, d_field, dorder)
        !  Derivative of 1-D Real Array (x in-processor)
        real*8 :: field(:)
        real*8 :: d_field(:)
        real*8, Allocatable :: temp_field(:)
        real*8, Pointer :: coefs(:,:)
        Integer :: dorder,r,stencil_radius, stencil_size, nbpts,istart,iend, stencil_center, offset
        Integer :: i

        Allocate(temp_field(1:n_x_fd))
        temp_field(:) = zero
        Nullify(coefs)
        Select Case (dorder)

            Case (1)
                coefs => d_coefs
                stencil_size = stencil_sizes(1)
                nbpts = dorder+boundary_accuracy_1
            Case (2)
                coefs => dd_coefs
                stencil_size = stencil_sizes(2)
                nbpts = dorder+boundary_accuracy_2
            Case (3)
                coefs => ddd_coefs
                stencil_size = stencil_sizes(3)
                nbpts = dorder+boundary_accuracy_3
            Case Default
                coefs => d_coefs
                stencil_size = stencil_sizes(1)
                nbpts = 1+boundary_accuracy_1
        End Select
        stencil_radius = stencil_size/2
        stencil_center = stencil_radius+1
        offset = (N_S-stencil_size)/2

        !boundaries first
        Do i = 1, nbpts
            temp_field(1) = temp_field(1)+field(i)*coefs(i,1)
            temp_field(n_x_fd) = temp_field(n_x_fd)+field(n_x_fd-i+1)*coefs(N_S-i+1,n_x_fd)
        Enddo
        ! Next do low r values that don't get the full stencil
        Do r = 2, stencil_radius
            istart = stencil_center-r+1
            Do i = istart, stencil_size
                temp_field(r) = temp_field(r)+field(r+i-stencil_center)*coefs(i+offset,r)
            Enddo
        Enddo
        ! Now the full stencil points
        Do r = stencil_center,n_x_fd-stencil_radius
            Do i = 1, stencil_size
                temp_field(r) = temp_field(r)+field(r+i-stencil_center)*coefs(i+offset,r)
            Enddo
        Enddo
        ! And finally the upper r values that don't get the full stencil
        Do r = n_x_fd-stencil_radius+1,n_x_fd-1
            iend = stencil_center+(n_x_fd-r)
            Do i =1, iend
                temp_field(r) = temp_field(r)+field(r+i-stencil_center)*coefs(i+offset,r)
            Enddo
        Enddo


        d_field = temp_field

        DeAllocate(temp_field)
        Nullify(coefs)
    End Subroutine d_by_dx1d




!************************************************************************
    Subroutine d_by_dx3d(field, dfield)
!************************************************************************
        Implicit None
        real*8, Intent(In) :: field(:,:,:)
        real*8, Intent(Out) :: dfield(:,:,:,:)
        real*8, Pointer :: coefs(:,:)
        Integer :: dorder,r,stencil_radius, stencil_size, nbpts,istart,iend, stencil_center, offset
        Integer :: i, ii, jj, iimin, iimax, jjmin, jjmax, dmax,dmin,d

        !************************************************************************
        ! Compute the radial derivatives along dimension 1 of a 3-D real array
        ! R is assumed to be in processor, so dimension 1 is assumed to run from 1 to n_x_fd
        ! Dimensions 2 and 3 can be of arbitrary size and have arbitrary starting indices.
        ! Number of derivatives computed is determined by the size of dfield's fourth dimension

        iimin = Lbound(field,2)
        iimax = Ubound(field,2)
        jjmin = Lbound(field,3)
        jjmax = Ubound(field,3)
         dmin = Lbound(dfield,4)
         dmax = Ubound(dfield,4)
        dfield(:,:,:,:) = zero
        Nullify(coefs)

        Do d = dmin, dmax
            dorder = d-dmin+1
            Select Case (dorder)

            Case (1)
                coefs => d_coefs
                stencil_size = stencil_sizes(1)
                nbpts = dorder+boundary_accuracy_1
            Case (2)
                coefs => dd_coefs
                stencil_size = stencil_sizes(2)
                nbpts = dorder+boundary_accuracy_2
            Case (3)
                coefs => ddd_coefs
                stencil_size = stencil_sizes(3)
                nbpts = dorder+boundary_accuracy_3
            Case Default
                coefs => d_coefs
                stencil_size = stencil_sizes(1)
                nbpts = 1+boundary_accuracy_1
            End Select
            stencil_radius = stencil_size/2
            stencil_center = stencil_radius+1
            offset = (N_S-stencil_size)/2
            Write(6,*)'I am in and dorder is : ', d
            Do jj = jjmin, jjmax
                Do ii = iimin, iimax

                    !boundaries first
                    Do i = 1, nbpts
                        dfield(  1,ii,jj,dorder) = dfield(  1,ii,jj,dorder)+field(      i,ii,jj)*coefs(i,1)
                        dfield(n_x_fd,ii,jj,dorder) = dfield(n_x_fd,ii,jj,dorder)+field(n_x_fd-i+1,ii,jj)*coefs(N_S-i+1,n_x_fd)
                    Enddo

                    ! Next do low r values that don't get the full stencil
                    Do r = 2, stencil_radius
                        istart = stencil_center-r+1
                        Do i = istart, stencil_size
                            dfield(r,ii,jj,dorder) = dfield(r,ii,jj,dorder)+field(r+i-stencil_center,ii,jj)*coefs(i+offset,r)
                        Enddo
                    Enddo

                    ! Now the full stencil points
                    Do r = stencil_center,n_x_fd-stencil_radius
                        Do i = 1, stencil_size
                            dfield(r,ii,jj,dorder) = dfield(r,ii,jj,dorder)+field(r+i-stencil_center,ii,jj)*coefs(i+offset,r)
                        Enddo
                    Enddo

                    ! And finally the upper r values that don't get the full stencil
                    Do r = n_x_fd-stencil_radius+1,n_x_fd-1
                        iend = stencil_center+(n_x_fd-r)
                        Do i =1, iend
                            dfield(r,ii,jj,dorder) = dfield(r,ii,jj,dorder)+field(r+i-stencil_center,ii,jj)*coefs(i+offset,r)
                        Enddo
                    Enddo

                Enddo
            Enddo

            Nullify(coefs)
        Enddo

    End Subroutine d_by_dx3d





!************************************************************************
    Subroutine d_by_dx3d2(dfield,dmax)
!************************************************************************
        Implicit None
        Real*8, Intent(InOut) :: dfield(:,:,:,0:)
        Integer, Intent(In) :: dmax
        Real*8, Pointer :: coefs(:,:)
        Integer :: dorder,r,stencil_radius, stencil_size, nbpts,istart,iend, stencil_center, offset
        Integer :: i, ii, jj, iimin, iimax, jjmin, jjmax, dmin,d

        !************************************************************************
        ! Compute the radial derivatives along dimension 1 of a 3-D real array
        ! R is assumed to be in processor, so dimension 1 is assumed to run from 1 to n_x_fd
        ! Dimensions 2 and 3 can be of arbitrary size and have arbitrary starting indices.
        ! Number of derivatives computed is determined by the size of dfield's fourth dimension

        iimin = Lbound(dfield,2)
        iimax = Ubound(dfield,2)
        jjmin = Lbound(dfield,3)
        jjmax = Ubound(dfield,3)
         dmin = 1


        dfield(:,:,:,1:dmax) = zero
        Nullify(coefs)

        Do d = dmin, dmax
            dorder = d-dmin+1
            Select Case (dorder)

            Case (1)
                coefs => d_coefs
                stencil_size = stencil_sizes(1)
                nbpts = dorder+boundary_accuracy_1
            Case (2)
                coefs => dd_coefs
                stencil_size = stencil_sizes(2)
                nbpts = dorder+boundary_accuracy_2
            Case (3)
                coefs => ddd_coefs
                stencil_size = stencil_sizes(3)
                nbpts = dorder+boundary_accuracy_3
            Case Default
                coefs => d_coefs
                stencil_size = stencil_sizes(1)
                nbpts = 1+boundary_accuracy_1
            End Select
            stencil_radius = stencil_size/2
            stencil_center = stencil_radius+1
            offset = (N_S-stencil_size)/2

            Do jj = jjmin, jjmax
                Do ii = iimin, iimax

                    !boundaries first
                    Do i = 1, nbpts
                        dfield(  1,ii,jj,dorder) = dfield(  1,ii,jj,dorder)+dfield(      i,ii,jj,0)*coefs(i,1)
                        dfield(n_x_fd,ii,jj,dorder) = dfield(n_x_fd,ii,jj,dorder)+dfield(n_x_fd-i+1,ii,jj,0)*coefs(N_S-i+1,n_x_fd)
                    Enddo

                    ! Next do low r values that don't get the full stencil
                    Do r = 2, stencil_radius
                        istart = stencil_center-r+1
                        Do i = istart, stencil_size
                            dfield(r,ii,jj,dorder) = dfield(r,ii,jj,dorder)+dfield(r+i-stencil_center,ii,jj,0)*coefs(i+offset,r)
                        Enddo
                    Enddo

                    ! Now the full stencil points
                    Do r = stencil_center,n_x_fd-stencil_radius
                        Do i = 1, stencil_size
                            dfield(r,ii,jj,dorder) = dfield(r,ii,jj,dorder)+dfield(r+i-stencil_center,ii,jj,0)*coefs(i+offset,r)
                        Enddo
                    Enddo

                    ! And finally the upper r values that don't get the full stencil
                    Do r = n_x_fd-stencil_radius+1,n_x_fd-1
                        iend = stencil_center+(n_x_fd-r)
                        Do i =1, iend
                            dfield(r,ii,jj,dorder) = dfield(r,ii,jj,dorder)+dfield(r+i-stencil_center,ii,jj,0)*coefs(i+offset,r)
                        Enddo
                    Enddo

                Enddo
            Enddo

            Nullify(coefs)
        Enddo

    End Subroutine d_by_dx3d2


!************************************************************************
    Subroutine d_by_dx3d3(ind,dind, buffer, dorder)
!************************************************************************
        Implicit None
        real*8, Intent(InOut) :: buffer(:,:,:,:)
        real*8, Pointer :: coefs(:,:)
        Integer, Intent(In) :: ind, dind
        Integer :: dorder,r,stencil_radius, stencil_size, nbpts,istart,iend, stencil_center, offset
        Integer :: i, ii, jj, iimin, iimax, jjmin, jjmax
        !************************************************************************
        ! Compute the radial derivatives along dimension 1 of a 3-D real array
        ! R is assumed to be in processor, so dimension 1 is assumed to run from 1 to n_x_fd
        ! Dimensions 2 and 3 can be of arbitrary size and have arbitrary starting indices.
        ! Derivative order specified by dorder
        ! derivative is taken on buffer(:,:,:,ind) and stored in buffer(:,:,:,dind)

        iimin = Lbound(buffer,2)
        iimax = Ubound(buffer,2)
        jjmin = Lbound(buffer,3)
        jjmax = Ubound(buffer,3)

        buffer(:,:,:,dind) = zero
        Nullify(coefs)


        Select Case (dorder)

            Case (1)
                coefs => d_coefs
                stencil_size = stencil_sizes(1)
                nbpts = dorder+boundary_accuracy_1
            Case (2)
                coefs => dd_coefs
                stencil_size = stencil_sizes(2)
                nbpts = dorder+boundary_accuracy_2
            Case (3)
                coefs => ddd_coefs
                stencil_size = stencil_sizes(3)
                nbpts = dorder+boundary_accuracy_3
            Case Default
                coefs => d_coefs
                stencil_size = stencil_sizes(1)
                nbpts = 1+boundary_accuracy_1
        End Select
        stencil_radius = stencil_size/2
        stencil_center = stencil_radius+1
        offset = (N_S-stencil_size)/2


        Do jj = jjmin, jjmax
            Do ii = iimin, iimax

                !boundaries first
                Do i = 1, nbpts
                    buffer(     1,ii,jj,dind) = buffer(     1,ii,jj,dind)+buffer(         i,ii,jj,ind)*coefs(i,1)
                    buffer(n_x_fd,ii,jj,dind) = buffer(n_x_fd,ii,jj,dind)+buffer(n_x_fd-i+1,ii,jj,ind)*coefs(N_S-i+1,n_x_fd)
                Enddo

                ! Next do low r values that don't get the full stencil
                Do r = 2, stencil_radius
                    istart = stencil_center-r+1
                    Do i = istart, stencil_size
                        buffer(r,ii,jj,dind) = buffer(r,ii,jj,dind)+buffer(r+i-stencil_center,ii,jj,ind)*coefs(i+offset,r)
                    Enddo
                Enddo

                ! Now the full stencil points
                Do r = stencil_center,n_x_fd-stencil_radius
                    Do i = 1, stencil_size
                        buffer(r,ii,jj,dind) = buffer(r,ii,jj,dind)+buffer(r+i-stencil_center,ii,jj,ind)*coefs(i+offset,r)
                    Enddo
                Enddo

                ! And finally the upper r values that don't get the full stencil
                Do r = n_x_fd-stencil_radius+1,n_x_fd-1
                    iend = stencil_center+(n_x_fd-r)
                    Do i =1, iend
                        buffer(r,ii,jj,dind) = buffer(r,ii,jj,dind)+buffer(r+i-stencil_center,ii,jj,ind)*coefs(i+offset,r)
                    Enddo
                Enddo

            Enddo
        Enddo

        Nullify(coefs)


    End Subroutine d_by_dx3d3



    Subroutine Initialize_dcoefs()
        Integer :: offset
        real*8, Allocatable :: temp_coefs(:,:)
        n_stencil_2_nl = 2*n_stencil_1-1
        boundary_accuracy_2_nl = n_stencil_2_nl-2
        stencil_sizes(1) = n_stencil_1
        stencil_sizes(2) = n_stencil_2
        stencil_sizes(3) = n_stencil_3
        max_stencil_size = MAXVAL(stencil_sizes)
        N_s = max_stencil_size
        kd = (N_s+1)/2
        Allocate(d_coefs(N_s,n_x_fd))
        Allocate(dd_coefs(N_s,n_x_fd))
        Allocate(ddd_coefs(N_s,n_x_fd))

        d_coefs(:,:) = zero
        dd_coefs(:,:) = zero
        ddd_coefs(:,:) = zero

        ! Generate 1st derivative coefficients
        Allocate(temp_coefs(n_stencil_1,n_x_fd))
        offset = (N_s-n_stencil_1)/2

        Call Generate_Coefs(x_fd,n_stencil_1,temp_coefs,1, boundary_accuracy_1)

        d_coefs(1+offset:N_s-offset,2:n_x_fd-1) = temp_coefs(1:n_stencil_1,2:n_x_fd-1)
        d_coefs(1:n_stencil_1,1) = temp_coefs(1:n_stencil_1,1)    ! offset the boundaries
        d_coefs(N_s-n_stencil_1+1:N_s,n_x_fd) = temp_coefs(1:n_stencil_1,n_x_fd)
        DeAllocate(temp_coefs)

        ! 2nd Derivative
        ! Note:  Much of this work will be undone when build_integration_matrix is called.
        !             Due to the strange behavior of the second derivative, we generate a set of
        !            second derivative coefficients for the interior points from a double
        !            application of the 1st derivative.  This is not ideal, but it is stable.
        Allocate(temp_coefs(n_stencil_2,n_x_fd))
        offset = (N_s-n_stencil_2)/2
        Call Generate_Coefs(x_fd,n_stencil_2,temp_coefs,2, boundary_accuracy_2)
        dd_coefs(1+offset:N_s-offset,2:n_x_fd-1) = temp_coefs(1:n_stencil_2,2:n_x_fd-1)
        dd_coefs(1:n_stencil_2,1) = temp_coefs(1:n_stencil_2,1)    ! offset the boundaries
        dd_coefs(N_s-n_stencil_2+1:N_s,n_x_fd) = temp_coefs(1:n_stencil_2,n_x_fd)

        DeAllocate(temp_coefs)

        !3rd Derivative
        Allocate(temp_coefs(n_stencil_3,n_x_fd))
        offset = (N_s-n_stencil_3)/2

        Call Generate_Coefs(x_fd,n_stencil_3,temp_coefs,3, boundary_accuracy_3)

        ddd_coefs(1+offset:N_s-offset,2:n_x_fd-1) = temp_coefs(1:n_stencil_3,2:n_x_fd-1)
        ddd_coefs(1:n_stencil_3,1) = temp_coefs(1:n_stencil_3,1)    ! offset the boundaries
        ddd_coefs(N_s-n_stencil_3+1:N_s,n_x_fd) = temp_coefs(1:n_stencil_3,n_x_fd)

        DeAllocate(temp_coefs)


    End Subroutine Initialize_dcoefs



    Subroutine Generate_Coefs(grid,n_stencil,coefs,deriv_order, boundary_accuracy)
        Implicit None
        real*8 :: grid(:), coefs(:,:)
        real*8 :: dmean, jfactorial, rescaling, dsum
        real*8, Allocatable :: temp(:,:), dees(:), work(:)
        Integer :: i,j,k, istart, iend, jstart, jend, itemp(2), nxy, ngrid
        Integer :: n_stencil, deriv_order, stencil_center, stencil_radius
        Integer :: boundary_accuracy, transition_type
        Integer :: n_boundary_pts, diff_dir, this_pt, info
        Integer, Allocatable :: ipiv(:), boundary_pts(:), directions(:)
      ! The stencil for the (fully) interior points is assumed to be centered (i.e. odd number of points)
        ! Note that this can be done for any grid - not just radial since variables like "r" and "radius" are never used
        !*******************************************************************************************************
        ! Transition Type:
        !      The transition_type determines how to handle points which are not boundary points but which are
        !     too close to the boundary to use the full stencil.
        ! 1:   Non-symmetric stencils are applied as the boundaries are approached
        ! 2:   Forward/backward methods are used as the left/right boundaries are approached (with same accuracy as the boundary accuracy)
        !
        !    Currently only type 1 is supported by the derivative routines and Implicit.F
        !*****************************************************************************************************
        !Write(6,*)deriv_order
         ngrid = Size(grid)
        stencil_radius = n_stencil/2
        stencil_center = stencil_radius+1
        coefs(:,:) = zero
        transition_type = derivative_transition_type

        If (transition_type .eq. 1) Then
            istart = 2
            iend = ngrid-1
        Endif
        If (transition_type .eq. 2) Then
            istart = stencil_center
            iend = ngrid-stencil_radius
        Endif
        n_boundary_pts = istart-1+ngrid-iend
        Allocate(boundary_pts(1:n_boundary_pts))
        Allocate(directions(1:n_boundary_pts))    ! 1 => left side, forward difference.  -1 => right side, backward difference

        Do i = 1, istart -1
            boundary_pts(i) = i
            directions(i) = 1
        Enddo
        Do i = istart, n_boundary_pts
            boundary_pts(i) = ngrid+(i-n_boundary_pts)
            directions(i) = -1
        Enddo

        Allocate(work(1:n_stencil))
        Allocate(ipiv(1:n_stencil))
        Allocate(dees(1:n_stencil))    ! The distance between the point we are taking the derivative at
                                                 ! and all points within the stencil

        ! Handle interior points first (symmetric and non-symmetric (if desired) stencils
        Do i = istart, iend

            dees(:) = zero

            ! Correct for the full stencil extending across the boundaries if necessary
            itemp(1) = 1
            itemp(2) = stencil_center-i+1
            jstart = MAXVAL(itemp)

            itemp(1) = n_stencil
            itemp(2) = (ngrid-i)+stencil_center
            jend = MINVAL(itemp)

            nxy = jend-jstart+1
            dsum = zero

            Do j = jstart, jend
                !Write(6,*)i,jstart,jend,stencil_center
                dees(j) = grid(i+j-stencil_center)-grid(i)
                dsum = dsum+ABS(dees(j))
            Enddo

            dees(stencil_center) = zero ! always true

            dmean = dsum/DBLE(nxy-1)
            dees(:) = dees(:)/dmean         ! Here we effectively do some matrix preconditioning so that we don't have one column
                                                 ! with values of order delta_r^4 and another of delta_r.  Will rescale the answer later.

            Allocate(temp(1:nxy,1:nxy))
            temp(:,:) = zero
            temp(:,1) = one
            Do j = 1, nxy
                Do k = 2, nxy
                    temp(j,k) = dees(j+jstart-1)**(k-1)
                Enddo
            Enddo
            ! Invert the matrix

            Call dgetrf(nxy,nxy,temp,nxy,ipiv,info)
            Call dgetri(nxy,temp,nxy,ipiv,work,nxy,info)    ! temp is now the inverse matrix

            jfactorial = 1.0D0
            Do j = 2, nxy
                jfactorial = jfactorial*DBLE(j-1)
                rescaling = (1.0D0/dmean)**(j-1)
                temp(j,:) = temp(j,:)*jfactorial*rescaling
            Enddo

            ! The coefficients for derivative of order n are contained in row n+1
            coefs(jstart:jend,i) = temp(deriv_order+1,1:nxy)

            DeAllocate(temp)
        Enddo
        DeAllocate(dees)
        !Now handle the boundary points
        nxy = deriv_order+boundary_accuracy
        Allocate(dees(1:nxy))
        Allocate(temp(1:nxy,1:nxy))

        Do i = 1,n_boundary_pts
            diff_dir = directions(i)
            this_pt = boundary_pts(i)
            dsum = 0.0D0
            dees(:) = zero
            If (diff_dir .eq. 1) Then
                Do j = 2, nxy
                    dees(j) = grid(this_pt+j-1)-grid(this_pt)
                    dsum = dsum+dees(j)
                Enddo
            Else    ! diff_dir = -1
                Do j = 1, nxy-1
                    dees(j) = grid(this_pt+(j-nxy))-grid(this_pt)
                    dsum = dsum+dees(j)
                Enddo
            Endif

            dmean = dsum/DBLE(nxy-1)
            dees(:) = dees(:)/dmean
            temp(:,:) = 0.0D0
            temp(:,1) = one

            Do j = 1, nxy
                Do k = 2, nxy
                    temp(j,k) = dees(j)**(k-1)
                Enddo
            Enddo

            Call dgetrf(nxy,nxy,temp,nxy,ipiv,info)
            Call dgetri(nxy,temp,nxy,ipiv,work,nxy,info)    ! temp is now the inverse matrix

            jfactorial = 1.0D0
            Do j = 2, nxy
                jfactorial = jfactorial*DBLE(j-1)
                rescaling = (1.0D0/dmean)**(j-1)
                temp(j,:) = temp(j,:)*jfactorial*rescaling
            Enddo

            If (diff_dir .eq. 1) Then
                coefs(1:nxy,this_pt) = temp(deriv_order+1,1:nxy)
            Else    !diff_dir = -1
                coefs(n_stencil-nxy+1:n_stencil, this_pt) = temp(deriv_order+1,1:nxy)
            Endif

        Enddo

        DeAllocate(temp)
        DeAllocate(dees)
        DeAllocate(work)
        DeAllocate(ipiv)

    End Subroutine Generate_Coefs


    !=============================
    ! Implicit.F row-loading routines  -- this should maybe be moved out later, but OK for now
    Subroutine Load_Interior_Rows(row,col,amp,dorder,mpointer)
        Integer, Intent(In) :: row, col, dorder ! ,rr,rp,n,np
        real*8, Intent(In) :: amp(:)
        real*8, Pointer, Dimension(:,:), Intent(In) :: mpointer

                Call Load_Interior_RowsFD(row,col,amp,dorder,mpointer)



    End Subroutine Load_Interior_Rows
    Subroutine Load_Single_Row(r,row,col,amp,dorder,mpointer, clear_row, amp_arr, boundary)
        Integer, Intent(In) :: r,row, col, dorder ! ,rr,rp,n,np
        real*8, Intent(In) :: amp
        real*8, Intent(In), Optional :: amp_arr(:)    ! lets us load every column of the row
        real*8, Pointer, Dimension(:,:), Intent(InOut) :: mpointer
        Logical, Intent(In), Optional :: clear_row, boundary


                Call Load_Single_RowFD(r,row,col,amp,dorder,mpointer, clear_row, amp_arr, boundary)



    End Subroutine Load_Single_Row
    Subroutine Load_Interior_RowsFD(row,col,amp,dorder,mpointer)
        Integer, Intent(In) :: row, col, dorder ! ,rr,rp,n,np
        Integer :: itemp(1:2), k, kk, kstart, kend, nbpts
        Integer :: r, stencil_size
        real*8, Intent(In) :: amp(:)
        real*8, Pointer, Dimension(:,:), Intent(In) :: mpointer
        real*8, Pointer :: coefs(:,:)
        Logical :: fd_ash = .true.

        If (fd_ash) Then
            If (dorder .gt. 0) Then
                nullify(coefs)
                Select Case (dorder)
                    Case (1)
                        coefs => d_coefs
                        stencil_size = stencil_sizes(1)
                        nbpts = dorder+boundary_accuracy_1
                    Case (2)
                        coefs => dd_coefs
                        stencil_size = stencil_sizes(2)
                        nbpts = dorder+boundary_accuracy_2
                    Case (3)
                        coefs => ddd_coefs
                        stencil_size = stencil_sizes(3)
                        nbpts = dorder+boundary_accuracy_3
                    Case Default
                        coefs => d_coefs
                        stencil_size = stencil_sizes(1)
                        nbpts = 1+boundary_accuracy_1
                End Select
                ! Boundaries first
                Do k = 1,nbpts
                    mpointer(1+row,k+col) = mpointer(1+row,k+col)+amp(1)*coefs(k,1)
                Enddo
                Do k = max_stencil_size-nbpts+1, max_stencil_size
                    mpointer(n_x_fd+row,col+n_x_fd+(k-max_stencil_size)) = &
                        mpointer(n_x_fd+row,col+n_x_fd+(k-max_stencil_size))+amp(n_x_fd)*coefs(k,n_x_fd)
                Enddo
                Do r = 2, n_x_fd-1
                    itemp(1) = 1
                    itemp(2) = kd-r+1
                    kstart = MAXVAL(itemp)

                    itemp(1) = stencil_size
                    itemp(2) = (n_x_fd-r)+kd
                    kend = MINVAL(itemp)
                    Do kk = kstart,kend        ! kk iterates over the stencil
                        k = r + kk - kd        ! k is how far from the diagonal of the matrix we are
                                                    ! Note that k is zero when we are at the stencil center (kk = kd)
                        mpointer(row+r,k+col) = mpointer(row+r,k+col) + amp(r)*coefs(kk,r)
                    Enddo
                Enddo
            Else    ! add amp to the diagonal
                !Do r = 2, n_x_fd-1
                Do r = 1, n_x_fd
                    k = r
                    mpointer(row+r,k+col) = mpointer(row+r,k+col) + amp(r)
                Enddo
            Endif
        Endif
    End Subroutine Load_Interior_RowsFD




    Subroutine Load_Single_RowFD(r,row,col,amp,dorder,mpointer, clear_row, amp_arr, boundary)
        Integer, Intent(In) :: r,row, col, dorder ! ,rr,rp,n,np
        Integer :: itemp(1:2), k,  nbpts, kstart,kend,kk
        Integer :: stencil_size
        real*8, Intent(In) :: amp
        real*8, Intent(In), Optional :: amp_arr(:)    ! lets us load every column of the row
        real*8, Pointer, Dimension(:,:), Intent(InOut) :: mpointer
        real*8, Pointer :: coefs(:,:)
        Logical, Intent(In), Optional :: clear_row, boundary
        Logical :: fd_ash = .true.
        If (present(clear_row)) Then
            ! clear everything in this row
            mpointer(r+row,:) = 0.0d0
        Endif
        If (fd_ash) Then
            If (dorder .gt. 0) Then
                nullify(coefs)
                Select Case (dorder)
                    Case (1)
                        coefs => d_coefs
                        stencil_size = stencil_sizes(1)
                        nbpts = dorder+boundary_accuracy_1
                    Case (2)
                        coefs => dd_coefs
                        stencil_size = stencil_sizes(2)
                        nbpts = dorder+boundary_accuracy_2
                    Case (3)
                        coefs => ddd_coefs
                        stencil_size = stencil_sizes(3)
                        nbpts = dorder+boundary_accuracy_3
                    Case Default
                        coefs => d_coefs
                        stencil_size = stencil_sizes(1)
                        nbpts = 1+boundary_accuracy_1
                End Select

                If ((r .ne. n_x_fd) .and. (r .ne. 1) ) Then
                    If (.not. present(boundary)) Then
                        itemp(1) = 1
                        itemp(2) = kd-r+1
                        kstart = MAXVAL(itemp)

                        itemp(1) = stencil_size
                        itemp(2) = (n_x_fd-r)+kd
                        kend = MINVAL(itemp)
                        Do kk = kstart,kend        ! kk iterates over the stencil
                            k = r + kk - kd        ! k is how far from the diagonal of the matrix we are
                                                    ! Note that k is zero when we are at the stencil center (kk = kd)
                            mpointer(row+r,k+col) = mpointer(row+r,k+col) + amp*coefs(kk,r)
                        Enddo
                    Endif
                Endif
                If (r .eq. 1) Then
                    Do k = 1,nbpts
                        mpointer(r+row,k+col) = mpointer(r+row,k+col)+amp*coefs(k,1)
                    Enddo
                Endif
                If (present(boundary)) Then
                If (r .eq. 2) Then            ! Some boundary conditions use level 2 with a derivative evaluated at r =1
                    Do k = 1,nbpts
                        mpointer(r+row,k+col) = mpointer(r+row,k+col)+amp*coefs(k,1)
                    Enddo
                Endif
                If (r .eq. (n_x_fd-1)) Then    ! Some boundary conditions use level n_x_fd-1 with a derivative evaluated at level n_x_fd
                    Do k = max_stencil_size-nbpts+1, max_stencil_size
                        mpointer(r+row,col+n_x_fd+(k-max_stencil_size)) = mpointer(r+row,col+n_x_fd+(k-max_stencil_size)) &
                                                                          + amp*coefs(k,n_x_fd)
                    Enddo
                Endif
                Endif
                If (r .eq. n_x_fd) Then
                    Do k = max_stencil_size-nbpts+1, max_stencil_size
                        mpointer(r+row,col+r+(k-max_stencil_size)) = mpointer(r+row,col+r+(k-max_stencil_size))+amp*coefs(k,n_x_fd)
                    Enddo
                Endif

            Else    ! add amp to the diagonal
                If (.not. present(amp_arr)) Then
                    mpointer(row+r,r+col) = mpointer(row+r,r+col) + amp
                Else
                    mpointer(row+r,col+1:col+n_x_fd) = mpointer(row+r,col+1:col+n_x_fd) + amp_arr(1:n_x_fd)
                Endif
            Endif
        Endif
    End Subroutine Load_Single_RowFD




End Module Finite_Difference
