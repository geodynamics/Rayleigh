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

Module Linear_Solve
    Use Finite_Difference
    Use Chebyshev_Polynomials, Only : Cheby_Grid
    !==========================================================================
    ! Generalized Implicit Time-stepping
    ! Currently assumes that implicit time stepping can be done with 1 dimension (only)
    ! in processor.  More general, 2 dimensional version could be added to
    ! (for example) implicitly advance Coriolis terms.
    Implicit None
    Integer, Save, Private :: n_equations, n_vars, n_modes    ! Number of variables and equations to be solved
    Integer, Save, Private :: n_modes_total
    Integer, Save, Private, Allocatable :: nsub_modes(:)

    Integer, Save, Private :: ndim1, ndim2, n_links, maximum_deriv_order        ! Variable used to keep track of linked equations (i.e. WPS)

    real*8, Save, Private :: LHS_time_factor, RHS_time_factor    ! Forward and Backward time-weighting of the implicit scheme.
    real*8, Allocatable :: dfield(:,:,:,:)
    Logical :: band_solve = .false.
    Logical :: sparse_solve = .false.
    Logical, Private :: chebyshev = .false.
    Real*8, Allocatable, Private :: temp_rhs(:,:,:)
    Type(Cheby_Grid), Pointer, Private :: cpgrid
    Type Data_arrays        ! support structure for the equation structure
        real*8, Allocatable :: data(:,:,:)    ! dimensioned (r, mode, derivative_order)
    End Type Data_arrays
    Type Coefficient_array        ! support structure for the equation structure
        real*8, Allocatable :: data(:,:)    ! dimensioned (r, mode, derivative_order)
    End Type Coefficient_Array
    !===============================================================
    !   Equation Structure
    Type Equation
        Type(Coefficient_Array), Allocatable :: coefs(:)
        real*8, Allocatable :: lhs(:,:)        ! The matrix to be inverted for this equation
                                                                    ! If an equation set is linked, this is only allocated for the primary
                                                                    ! equation in the link.
        real*8, Pointer, dimension(:,:) :: mpointer    ! Points to LHS of the primary equation in a linked set (or LHS if not linked)

        real*8, Allocatable:: rhs(:,:,:)    ! RHS array for all modes of a given equation.  Only allocated for mode 1
        real*8, Pointer, dimension(:,:,:) :: rhs_pointer    ! Points to the appropriate portion of the joint RHS array
                                                                                        ! for this equation dimensioned r,real/imag, sub-mode
        Integer, Allocatable :: pivot(:)        ! Pivot array corresponding to the LHS matrix

        Integer, Allocatable :: links(:)                    ! Array containing indices of other equations that this equation is linked to.
        Integer :: nlinks = 1                                ! Number of equations this equation is linked to (including itself).

        Logical :: primary = .true.                        ! True if this equation is the primary equation of a linked set
        Logical :: linked = .false.                        ! True if this equation is a member of a linked set
        Logical :: solvefor = .false.                        ! Set to true if this mode is solved for in this equation

        Integer, Allocatable :: colblock(:)                ! The column block for each variable reference by this equation
        Integer :: rowblock = 0                                ! The rowblock of this equation in the LHS matrix


        !//////////////////////////////////////////////////////////////////\
        !Additional attributes used if sparse solve is active
        Integer, Allocatable :: sparse_ia(:), sparse_ja(:)
        Real*8, Allocatable :: sparse_mat(:)
        Integer*8 :: pt(64)
        Integer :: mxfct, mnum, mtype, phase, n, nrhs, error, msglvl
        Integer :: solver, idum, nrows
        Integer :: iparm(64)
        Real*8  :: dparm(64), ddum
        logical :: sparse_initialized = .false.

        Contains
        Procedure :: LU_Solve_Sparse
        Procedure :: LU_Decompose_Sparse

    End Type Equation




    !========================================================
    !    Variable Structure
    Type cdarrays
        Complex*16, allocatable :: data(:,:)        ! Dimensioned 1:ndim1, 1:n_modes_total
    End Type cdarrays

    Type ddarrays
        real*8, Allocatable :: data(:,:,:)        ! Dimensioned 1:ndim1, 1:2 (real/complex), 1:n_modes_total
    End Type ddarrays


    Type Variable
        Integer :: var_start        ! The first dimension starting index of this variable in the RHS array of equation_set(1,equ_ind)
        Integer :: var_end        ! The first dimension ending index
        Integer :: equ_ind        ! The equation responsible for the RHS array this variable is stored in following the solve.
        real*8, pointer, dimension(:,:,:) :: data    ! The portion of the RHS array of equation_set(1,equ_ind) that holds this variable
        Logical, Allocatable :: in_equation(:,:,:)
        Integer :: max_dorder ! The maximum derivative order calculated for this variable (for linear or nonlinear terms)
        Integer :: dmax    ! The maximum derivative order that needs to be saved for this variable
                                                ! when leaving the implicit configuration.
        Integer :: dsave_type = 2
        Type(cdarrays), Allocatable :: cderivs(:)
        Type(ddarrays), Allocatable :: derivs(:)
    End Type Variable
    !=================================================

    Type(equation), Allocatable, Target :: equation_set(:,:)
    Type(data_arrays), Allocatable :: Implicit_RHS(:)
    Type(variable), Allocatable, Target :: var_set(:)

    Contains
    Subroutine Use_Chebyshev(ingrid)
        Type(Cheby_Grid), Pointer, Intent(InOut) :: ingrid
        chebyshev = .true.
        nullify(cpgrid)
        cpgrid => ingrid
    End Subroutine Use_Chebyshev
    Subroutine Use_BandSolve()
        band_solve = .true.
    End Subroutine Use_BandSolve
    Subroutine Use_SparseSolve()
        sparse_solve = .true.
    End Subroutine Use_SparseSolve

    Subroutine DeAllocate_Derivatives()
        Implicit None
        Integer :: i, j, dmax, dtype
        ! DeAllocate Every Variable's derivative array
        Do i = 1, n_vars
            dmax = var_set(i)%dmax
            dtype = var_set(i)%dsave_type
            If (dmax .ge. 0) Then
                If (dtype .eq. 1) Then
                    Do j = 0, dmax
                        If (allocated(var_set(i)%cderivs(j)%data)) DeAllocate(var_set(i)%cderivs(j)%data)
                    Enddo
                Else
                    Do j = 0, dmax
                        If (allocated(var_set(i)%derivs(j)%data)) DeAllocate(var_set(i)%derivs(j)%data)
                    Enddo
                Endif
            Endif
        Enddo
    End Subroutine DeAllocate_Derivatives

    Subroutine Set_Deriv_Save(varind,maxd)
        Implicit None
        Integer, Intent(In) :: maxd, varind
        Integer :: dtype
        var_set(varind)%dmax = maxd
        dtype = var_set(varind)%dsave_type
        If (dtype .eq. 1) Then
            Allocate(var_set(varind)%cderivs(0:maxd))
        Else
            Allocate(var_set(varind)%derivs(0:maxd))
        Endif
    End Subroutine Set_Deriv_Save

    Subroutine Set_Time_Factors(lhs_factor,rhs_factor)
        Implicit None
        real*8, Intent(In) :: lhs_factor, rhs_factor
        lhs_time_factor = lhs_factor
        rhs_time_factor = rhs_factor
    End Subroutine Set_Time_Factors

    !=======================================================
    !                Initialization Routines
    Subroutine Initialize_Equation_Set(neq,nvar,ndim,nmode, nsub,nphase)
        Implicit None
        Integer, Intent(In) :: neq, nvar, nmode, ndim, nsub(:)
        Integer :: i,j, nphase
        maximum_deriv_order = 0
        n_modes = nmode
        n_modes_total = Sum(nsub)

        Allocate(nsub_modes(1:n_modes))
        nsub_modes = nsub
        n_equations = neq
        n_vars = nvar
        ndim1 = ndim
        ndim2 = nphase
        Allocate(equation_set(1:n_modes,1:n_equations))
        Allocate(Implicit_RHS(1:n_equations))
        Do j = 1, n_equations
            Do i = 1, n_modes
                Allocate(equation_set(i,j)%coefs(1:n_vars))
                Allocate(equation_set(i,j)%colblock(1:n_vars))
                equation_set(i,j)%colblock(:) =0
            Enddo
        Enddo
        Allocate(var_set(1:nvar))
        Do i = 1, nvar
            var_set(i)%var_start = 1
            var_set(i)%var_end = ndim1
            var_set(i)%equ_ind = i
            var_set(i)%max_dorder = 0
        Enddo
    End Subroutine Initialize_Equation_Set

    Subroutine Link_Equations(eq_links, var_links,nlinks,mode)
        Implicit None
        Integer, Intent(In) :: nlinks, mode, var_links(:), eq_links(:)
        Integer :: i,j
        Do i = 1, nlinks
            ! Record which equations are linked with the equation structure
            Allocate(equation_set(mode,eq_links(i))%links(1:nlinks))
            equation_set(mode,eq_links(i))%links(:) = eq_links
            equation_set(mode,eq_links(i))%linked = .true.
            equation_set(mode,eq_links(i))%nlinks = nlinks
            ! Establish the row block for each equation within the equation matrix
            equation_set(mode,eq_links(i))%rowblock = (i-1)*ndim1
            If (i .gt. 1) Then
                equation_set(mode,eq_links(i))%primary = .false.
            Endif

            Do j = 1, nlinks
                ! Establish the colums within the matrix pertaining to each variable
                equation_set(mode,eq_links(i))%colblock(var_links(j)) = (j-1)*ndim1

            Enddo
            var_set(var_links(i))%var_start = (i-1)*ndim1+1        ! mode independent for now
            var_set(var_links(i))%var_end = i*ndim1
            var_set(i)%equ_ind = eq_links(1)
        Enddo
    End Subroutine Link_Equations

    Subroutine Initialize_Equation_Coefficients(eq, var,dorder, mode)
        Implicit None
        Integer, Intent(In) :: mode,eq, var, dorder

        Allocate( equation_set(mode,eq)%coefs(var)%data(1:ndim1,0:dorder)  )
        equation_set(mode,eq)%coefs(var)%data(:,:) = 0.0d0
        equation_set(mode,eq)%solvefor = .true.
        If (dorder .gt. maximum_deriv_order) maximum_deriv_order = dorder
        If (dorder .gt. var_set(var)%max_dorder) var_set(var)%max_dorder = dorder
    End Subroutine Initialize_Equation_Coefficients

    Subroutine Reset_Equation_Coefficients()
        Implicit None
        Integer :: i, j, k
        Do k = 1, n_equations
            Do j = 1, n_modes
                If (allocated(equation_set(j,k)%lhs)) equation_set(j,k)%lhs(:,:) = 0.0d0
                if (band_solve .or. sparse_solve) then
                    If (allocated(equation_set(j,k)%lhs) .and. equation_set(1,k)%primary) Then
                        DeAllocate(equation_set(j,k)%lhs)
                    Endif
                endif

                Do i = 1, n_vars
                    If(allocated(equation_set(j,k)%coefs(i)%data)) Then
                        equation_set(j,k)%coefs(i)%data(:,:) = 0.0d0
                    Endif
                Enddo
            Enddo
        Enddo
    End Subroutine Reset_Equation_Coefficients

    Subroutine Finalize_Equations()
        Implicit None
        Integer :: k,j

        Do k = 1, n_vars
            j = var_set(k)%max_dorder
            Allocate(var_set(k)%in_equation(1:n_modes,1:n_equations,0:j))
            var_set(k)%in_equation(:,:,:) = .false.
        Enddo
    End Subroutine Finalize_Equations

    Subroutine DeAllocate_LHS(mode_ind)
        Implicit None
        Integer, Intent(In) :: mode_ind
        Integer :: k,j

        j = mode_ind
        Do k = 1, n_equations
            If (equation_set(j,k)%solvefor) Then
                ! Only DeAllocate matrix information for modes we actually solve for
                If (equation_set(j,k)%primary) Then
                    If (allocated(equation_set(j,k)%lhs)) Then
                        DeAllocate(equation_set(j,k)%lhs)
                    Endif
                    If (allocated(equation_set(j,k)%pivot)) Then
                        DeAllocate(equation_set(j,k)%pivot)
                    Endif
                    nullify(equation_set(j,k)%mpointer)
                Else
                    nullify(equation_set(j,k)%mpointer)
                Endif
            Endif
        Enddo

    End Subroutine DeAllocate_LHS

    Subroutine Allocate_LHS(mode_ind)
        Implicit None
        Integer, Intent(In) :: mode_ind
        Integer :: k,j, ndim, ind

        ! Check to see if each equation's matrix has already been allocated.
        ! If not, allocate it.  Account for linked equations.
        j = mode_ind
        Do k = 1, n_equations
            If (equation_set(j,k)%solvefor) Then        ! Only Allocate matrix information for modes we actually solve for
                If (equation_set(j,k)%primary) Then
                    ndim = ndim1*equation_set(j,k)%nlinks
                    If (.not. allocated(equation_set(j,k)%lhs)) Then
                        Allocate(equation_set(j,k)%lhs(1:ndim,1:ndim))
                    Endif
                    equation_set(j,k)%lhs(1:ndim,1:ndim) = 0.0d0
                    If (.not. allocated(equation_set(j,k)%pivot)) Then
                        Allocate(equation_set(j,k)%pivot(1:ndim))
                    Endif
                    equation_set(j,k)%mpointer => equation_set(j,k)%lhs
                Else
                    ind = equation_set(j,k)%links(1)
                    equation_set(j,k)%mpointer => equation_set(j,ind)%lhs
                Endif
            Endif
        Enddo

    End Subroutine Allocate_LHS

    Subroutine Allocate_RHS(zero_rhs)
        Implicit None
        Integer :: k,j, ndim, ind, indx, nsub
        Logical, Intent(In), Optional :: zero_rhs
        ! We use one large RHS for each set of equations, with the first mode's
        ! equation object holding that RHS. Each mode's equation object points to
        ! the appropriate parts of that RHS space.  The idea is that the RHS can be
        ! accessed more efficiently when adding nonlinear and CN terms.

        Do k = 1, n_equations

            !/// Allocation
            If (equation_set(1,k)%primary) Then
                ndim = ndim1*equation_set(1,k)%nlinks
                Allocate(equation_set(1,k)%rhs(1:ndim,1:ndim2,1:n_modes_total))
                If(present(zero_rhs)) Then
                    If (zero_rhs .eqv. .true.) Then
                        equation_set(1,k)%rhs(:,:,:) = 0.0d0
                    Endif
                Endif
            Endif

            !/// Pointing
            indx = 1
            Do j = 1, n_modes
                nsub = nsub_modes(j)

                If (equation_set(j,k)%primary) Then
                    equation_set(j,k)%rhs_pointer => equation_set(1,k)%rhs(:,:,indx:indx+nsub-1)
                Else
                    ind = equation_set(j,k)%links(1)
                    equation_set(j,k)%rhs_pointer => equation_set(1,ind)%rhs(:,:,indx:indx+nsub-1)

                Endif
                indx = indx+nsub

            Enddo

        Enddo

    End Subroutine Allocate_RHS

    Subroutine DeAllocate_RHS()
        Implicit None
        Integer :: k, j
        ! DeAllocate the RHS for each equation
        Do k = 1, n_equations

            Do j = 1, n_modes
                If (Allocated(equation_set(j,k)%rhs)) DeAllocate(equation_set(j,k)%rhs)
                Nullify(equation_set(j,k)%rhs_pointer)
            Enddo
        Enddo
    End Subroutine DeAllocate_RHS



    !===========================
    !  Matrix Solve Routine
    Subroutine Implicit_Solve
        Implicit None
        Integer :: j, k,maxlink
        If (band_solve) Then
            maxlink = 1
            Do k = 1, n_equations
                maxlink = max(maxlink,equation_set(1,k)%nlinks)
            Enddo
            if (maxlink .gt. 1) Allocate(temp_rhs(1:ndim1*maxlink,1:ndim2,1:n_modes_total))
        Endif
        Do k = 1, n_equations
            If (band_solve .and. (equation_set(1,k)%nlinks .gt. 1) .and. equation_set(1,k)%primary) Then
                Call Band_Arrange_RHS(k)
            Endif
            Do j = 1, n_modes
                If (equation_set(j,k)%primary .and. equation_set(j,k)%solvefor) Then
                    If (band_solve) Then
                        Call lu_solve_band(Equation_set(j,k)%LHS , Equation_Set(j,k)%rhs_pointer , Equation_Set(j,k)%Pivot)
                    Else If (sparse_solve) Then
                        !Call LU_Solve_Sparse(j,k)
                        Call Equation_set(j,k)%LU_Solve_Sparse()
                    Else
                        Call lu_solve_full(Equation_set(j,k)%LHS , Equation_Set(j,k)%rhs_pointer , Equation_Set(j,k)%Pivot)
                    Endif
                Endif
            Enddo
            If (band_solve .and. (equation_set(1,k)%nlinks .gt. 1) .and. equation_set(1,k)%primary) Then
                Call Band_ReArrange_RHS(k)
            Endif
        Enddo
        If (band_solve) Then
            if (maxlink .gt. 1) DeAllocate(temp_rhs)
        Endif

    End Subroutine Implicit_Solve

    Subroutine LU_Decompose_Matrices()
        Implicit None
        Integer :: j, k
        Do k = 1, n_equations
            Do j = 1, n_modes
                If (equation_set(j,k)%primary .and. equation_set(j,k)%solvefor) Then
                    if (band_solve) Then

                        Call lu_decompose_band(Equation_set(j,k)%LHS , Equation_Set(j,k)%Pivot)

                    else if (sparse_solve) Then
                        !Call lu_decompose_sparse(j,k)
                        Call Equation_set(j,k)%lu_decompose_sparse()
                    else
                        Call lu_decompose_full(Equation_set(j,k)%LHS , Equation_Set(j,k)%Pivot)
                    Endif
                Endif
            Enddo
        Enddo


    End Subroutine LU_Decompose_Matrices

    Subroutine Point_Variables()
        Implicit None
        Integer :: vstart, vend, i, equ_ind
        Do i = 1, n_vars
            vstart = var_set(i)%var_start
            vend   = var_set(i)%var_end
            equ_ind = var_set(i)%equ_ind
            var_set(i)%data => equation_set(1,equ_ind)%rhs(vstart:vend,:,:)
        Enddo
    End Subroutine Point_Variables


    Subroutine Add_Implicit_Term(eqind,varind,dorder,amp,mode, static)
        ! This subroutine adds derivative coefficients into the implicit matrix
        ! and also saves them into the correct part of the equation coefficient arrays
        ! If static is set, it means that this term is a piece of:
        ! (W^n+1-W^n)/dt  and gets added to the RHS and LHS with no dt factor
        ! and with the same amplitude.
        Implicit None
        Integer, Intent(In) :: eqind, varind, dorder, mode
        real*8, Intent(InOut) :: amp(:)
        real*8 :: time_amp
        Integer :: rowblock, colblock
        Logical, Intent(In), optional :: static
        real*8, Pointer, Dimension(:,:) :: mpointer
        rowblock = equation_set(mode,eqind)%rowblock
        colblock = equation_set(mode,eqind)%colblock(varind)
        mpointer => equation_set(mode,eqind)%mpointer

        If(present(static)) Then
            equation_set(mode,eqind)%coefs(varind)%data(:,dorder) = amp + &
                & equation_set(mode,eqind)%coefs(varind)%data(:,dorder)
        Else
            equation_set(mode,eqind)%coefs(varind)%data(:,dorder) = amp*RHS_Time_Factor+ &
                & equation_set(mode,eqind)%coefs(varind)%data(:,dorder)
        Endif

        If (present(static)) Then
            If (chebyshev) Then
                Call Load_Interior_Rows_Cheby(rowblock, colblock,amp,dorder,mpointer)
            Else
                Call Load_Interior_Rows(rowblock, colblock,amp,dorder,mpointer) ! Uses existing version of load rows
            Endif
        Else
            !//// This was here before I added the if block above
            !equation_set(mode,eqind)%coefs(varind)%data(:,dorder) = amp
            time_amp = LHS_time_factor
            amp = amp*time_amp
            If (chebyshev) Then
                Call Load_Interior_Rows_Cheby(rowblock, colblock,amp,dorder,mpointer)
            Else
                Call Load_Interior_Rows(rowblock, colblock,amp,dorder,mpointer) ! Uses existing version of load rows
            Endif
            amp = amp/time_amp
            var_set(varind)%in_equation(mode,eqind,dorder) = .true.
        Endif

    End Subroutine Add_Implicit_Term

    Subroutine Compute_Implicit_RHS()
        Implicit None
        Integer :: i,j,k,d,nsub, indx, djmax, ii, jj
        real*8, Pointer, Dimension(:) :: coefs

        nullify(coefs)
        Call point_variables()    ! Point the variable data pointers to the appropriate part of the RHS arrays
        ! Zero out the RHS
        Do j = 1, n_equations
            ! need to allocate somewhere
            If (.not. Allocated(Implicit_RHS(j)%data)) Then
                Allocate(Implicit_RHS(j)%data(1:ndim1,1:ndim2,1:n_modes_total))
            Endif
            Implicit_RHS(j)%data = 0.0d0
        Enddo

        ! Allocate the array that will hold each variable and its derivatives
        Allocate(dfield(1:ndim1,1:ndim2,1:n_modes_total,0:maximum_deriv_order))


        Do k = 1, n_vars    ! Iterate over each variable

            dfield(:,:,:,0) = var_set(k)%data    ! Copy variable from the RHS array
            !dpointer => dfield(:,:,:,0)
            !Call Copy_Variable(k,dpointer)
            djmax = var_set(k)%max_dorder
            Call d_by_dx(dfield,djmax)    ! Compute all radial derivatives needed for this variable.
            Call Save_Derivatives(k) ! save this derivative for transposition later

            Do d = 0, djmax

                Do j = 1, n_equations

                        indx = 1
                        Do i = 1, n_modes
                        If (var_set(k)%in_equation(i,j,d))    Then    ! Add this derivative to this equation
                            coefs => equation_set(i,j)%coefs(k)%data(:,d)    ! might want to restructure how this is stored

                            nsub = nsub_modes(i)-1
                            Do jj = indx, indx+nsub
                                Do ii = 1, ndim2
                                    Implicit_RHS(j)%data(:,ii,jj) = Implicit_RHS(j)%data(:,ii,jj)+dfield(:,ii,jj,d)*coefs
                                Enddo
                            Enddo


                        Endif
                        indx = indx+nsub_modes(i)
                        Enddo

                Enddo
            Enddo
        Enddo

        DeAllocate(dfield)
    End Subroutine Compute_Implicit_RHS

    Subroutine Add_Implicit_RHS()
        Implicit None
        Integer :: i,j,k,istart, ii
        Do k = 1, n_equations
            If (equation_set(1,k)%primary) Then
                If (equation_set(1,k)%nlinks .gt. 1) Then
                    Do j = 1, equation_set(1,k)%nlinks
                        i = equation_set(1,k)%links(j)

                        istart = (j-1)*ndim1
                        !iend   = istart+ndim1
                        Do ii = 1, ndim1
                        equation_set(1,k)%rhs(istart+ii,:,:) = equation_set(1,k)%rhs(istart+ii,:,:) +&
                            & RHS_Time_Factor*Implicit_RHS(i)%data(ii,:,:)
                        Enddo
                    Enddo
                Else
                    equation_set(1,k)%rhs = equation_set(1,k)%rhs+RHS_Time_Factor*Implicit_RHS(k)%data
                Endif
            Endif
        Enddo
    End Subroutine Add_Implicit_RHS

    Subroutine Set_RHS(eqid,set_to)
        ! Set the RHS of equation eqid to the value of set_to
        Implicit None
        Real*8, Intent(InOut) :: set_to(:,:,:)
        Integer, Intent(In) :: eqid
        Integer :: istart,iend, ind


        If (n_modes .gt. 0) Then
            !Primary equation object always has the full, allocated rhs.
            ind = eqid
            if (.not. equation_set(1,eqid)%primary) then
                ind = equation_set(1,eqid)%links(1)
            endif
            ! Individual RHS's inhabit row ranges defined by rowblock and ndim1
            istart = equation_set(1,eqid)%rowblock+1
            iend = istart+ndim1-1


            equation_set(1,ind)%rhs(istart:iend,:,:) = set_to(1:ndim1,:,:)
        Endif
    End Subroutine Set_RHS

    Subroutine Zero_RHS(eqid)
        ! Set the RHS of equation eqid to zero
        ! This is useful for when we want to remove explicity time-dependence
        ! from the equation (as with mantle convection problems).
        Implicit None
        Integer, Intent(In) :: eqid
        Integer :: istart,iend, ind


        If (n_modes .gt. 0) Then
            !Primary equation object always has the full, allocated rhs.
            ind = eqid
            if (.not. equation_set(1,eqid)%primary) then
                ind = equation_set(1,eqid)%links(1)
            endif
            ! Individual RHS's inhabit row ranges defined by rowblock and ndim1
            istart = equation_set(1,eqid)%rowblock+1
            iend = istart+ndim1-1


            equation_set(1,ind)%rhs(istart:iend,:,:) = 0.0d0
        Endif
    End Subroutine Zero_RHS

    Subroutine Get_All_RHS(buffer)
        Implicit None
        ! Copy equation structure RHSs to the buffer (dlink RHSs)
        Real*8, Intent(InOut) :: buffer(:,:,:,1:)
        Integer :: i, ind,istart,iend
        If (n_modes .gt. 0) Then
            Do i = 1, n_equations
                ! Individual RHS's inhabit row ranges defined by rowblock and ndim1
                istart = equation_set(1,i)%rowblock+1
                iend = istart+ndim1-1
                !Primary equation object always has the full, allocated rhs.
                ind = i
                if (.not. equation_set(1,i)%primary) then
                    ind = equation_set(1,i)%links(1)
                endif

                buffer(1:ndim1,:,:,i) = equation_set(1,ind)%rhs(istart:iend,:,:)
            Enddo
        Endif
    End Subroutine Get_All_RHS

    Subroutine Set_All_RHS(buffer)
        ! Copy RHS's from the buffer into the equation structure
        ! Buffer RHS's are assumed to be unlinked.
        Implicit None
        Real*8, Intent(InOut) :: buffer(:,:,:,1:)
        Integer :: i, ind,istart,iend

            Do i = 1, n_equations
                ! Individual RHS's inhabit row ranges defined by rowblock and ndim1
                istart = equation_set(1,i)%rowblock+1
                iend = istart+ndim1-1
                !Primary equation object always has the full, allocated rhs.
                ind = i
                if (.not. equation_set(1,i)%primary) then
                    ind = equation_set(1,i)%links(1)
                endif

                equation_set(1,ind)%rhs(istart:iend,:,:) = buffer(1:ndim1,:,:,i)
            Enddo
    End Subroutine Set_All_RHS

    Subroutine Add_To_All_RHS(buffer,mfactor)
        ! Add mfactor*buffer to the equation set RHS's
        ! buffer is assumed to be in unlinked format
        ! Right now this is just an easy way of adding an AB term
        Implicit None
        Real*8, Intent(InOut) :: buffer(:,:,:,1:)
        Real*8, Intent(In) :: mfactor
        Integer :: i, ind,istart,iend

            Do i = 1, n_equations
                ! Individual RHS's inhabit row ranges defined by rowblock and ndim1
                istart = equation_set(1,i)%rowblock+1
                iend = istart+ndim1-1
                !Primary equation object always has the full, allocated rhs.
                ind = i
                if (.not. equation_set(1,i)%primary) then
                    ind = equation_set(1,i)%links(1)
                endif

                equation_set(1,ind)%rhs(istart:iend,:,:) = equation_set(1,ind)%rhs(istart:iend,:,:)+mfactor*buffer(1:ndim1,:,:,i)
            Enddo
    End Subroutine Add_to_all_RHS


    !///////// TODO
    Subroutine Add_Derivative(eqid,varid, dorder,addto,dfield,dind)
        ! Adds dfield multplied by the coefficients appropriate for
        ! equation eqid to the array addto.  dfield is assumed to be the derivative of
        ! variable varid
        Implicit None
        Real*8, Intent(InOut) :: addto(:,:,:,:), dfield(:,:,:,:)
        Integer, Intent(In) :: eqid, varid, dorder,dind
        Integer :: i,nsub, indx, ii, jj
        real*8, Pointer, Dimension(:) :: coefs
        Integer :: itmp
        nullify(coefs)
        indx = 1
        Do i = 1, n_modes
            nsub = nsub_modes(i)-1 !like n_m locally stored for ell(i)

            ! TODO:   N.F. June 27, 2017
            ! This logic should be streamlined.  Possibly best to create a support array
            ! that indicates the maximum dorder for a particular mode....
            If (allocated(equation_set(i,eqid)%coefs(varid)%data)) then
            itmp =UBOUND( equation_set(i,eqid)%coefs(varid)%data,2 )
            If (dorder .le. itmp) Then
            coefs => equation_set(i,eqid)%coefs(varid)%data(:,dorder)
            Do jj = indx, indx+nsub
                Do ii = 1,ndim2
                    addto(:,ii,jj,eqid) = addto(:,ii,jj,eqid)+dfield(:,ii,jj,dind)*coefs
                    ! if we didn't want to use the pointer (worth testing), we could just write this...
                    !addto(:,ii,jj) = addto(:,ii,jj)+dfield(:,ii,jj)*equation_set(i,eqid)%coefs(varid)%data(:,dorder)
                Enddo
            Enddo
            Endif
            Endif
            indx = indx+nsub_modes(i)
        Enddo
        nullify(coefs)
        ! It might also be worth making a coef array like coefs(:,dorder,i,eqid), but it would be memory wasteful
        !  since no variable but W has a third derivative, and since we only need the first derivative for P
        ! Also, all variables do not appear in equation equation
        ! More notes:  The number of elements of that coef array would then be N_r*nel_local*3*n_eq*n_var = N_r*nel_local*48 for hydro
    End Subroutine Add_Derivative




    Subroutine Save_Derivatives(varind)
        Implicit None
        Integer, Intent(in) :: varind
        Integer :: i, dmax, dtype
        dtype = var_set(varind)%dsave_type
        dmax = var_set(varind)%dmax
        If (dtype .eq. 1) Then
            If (dmax .gt. 0) Then        ! changed ge to gt here and just below.  NF Oct 18, 2013
                Do i = 0, dmax
                    If (.not. Allocated(var_set(varind)%cderivs(i)%data)) Then
                        Allocate(var_set(varind)%cderivs(i)%data(1:ndim1,1:n_modes_total))
                    Endif
                        !var_set(varind)%cderivs(i)%data = Cmplx(dfield(:,1,:,i), dfield(:,2,:,i),Precision)

                Enddo
            Endif
        Else
            If (dmax .gt. 0) Then
                Do i = 0, dmax
                    If (.not. Allocated(var_set(varind)%derivs(i)%data)) Then
                        Allocate(var_set(varind)%derivs(i)%data(1:ndim1,1:ndim2,1:n_modes_total))
                    Endif
                        var_set(varind)%derivs(i)%data = dfield(:,:,:,i)

                Enddo
            Endif

        Endif
    End Subroutine Save_Derivatives

    Subroutine write_matrix(mode,eqind)
        Integer, Intent(In) :: mode, eqind
        Integer :: i, j, nlinks
        Integer :: sig = 314
        real*8, Pointer, Dimension(:,:) :: mpointer
        Character*120 :: filename
        filename = 'matrix_dump'
        mpointer => equation_set(mode,eqind)%mpointer
        nlinks = equation_set(mode,eqind)%nlinks
        Open(unit=15,file=filename,form='unformatted', status='replace',access='stream')
        Write(15)sig
        Write(15)ndim1*nlinks
        Write(15)ndim1*nlinks
        Write(15)((mpointer(i,j), i = 1, ndim1*3), j = 1, ndim1*3)
        Close(15)

    End Subroutine write_matrix

    Subroutine print_row(mode,row,eqind)
        Integer, Intent(In) :: mode, row, eqind
        Integer :: rowblock,i
        real*8, Pointer, Dimension(:,:) :: mpointer
        mpointer => equation_set(mode,eqind)%mpointer
        rowblock = equation_set(mode,eqind)%rowblock
        Do i = 1, ndim1*3
            Write(6,*)i, mpointer(rowblock+row,i)
        Enddo
    End Subroutine print_row

    Subroutine print_column(mode,col,eqind)
        Integer, Intent(In) :: mode, col, eqind
        Integer :: rowblock,i
        real*8, Pointer, Dimension(:,:) :: mpointer
        mpointer => equation_set(mode,eqind)%mpointer
        Do i = 1, ndim1*3
            Write(6,*)i, mpointer(i,col)
        Enddo
    End Subroutine print_column


    Subroutine Load_BC(mode,row,eqind,varind,amp,dorder,integral)
        Implicit None
        Integer, Intent(In) :: mode, row, eqind,varind,dorder
        Integer :: colblock, rowblock
        real*8, Intent(In) :: amp
        real*8, Intent(In), Optional :: integral(:)
        real*8, Pointer, Dimension(:,:) :: mpointer
        mpointer => equation_set(mode,eqind)%mpointer
        colblock = equation_set(mode,eqind)%colblock(varind)
        rowblock = equation_set(mode,eqind)%rowblock

        If (present(integral)) Then
            mpointer(rowblock+row,colblock+1:colblock+ndim1) = integral(1:ndim1)

        Else
            If (chebyshev) Then
                Call Load_Single_Row_Cheby(row,rowblock,colblock,amp,dorder,mpointer, boundary = .true.)
            Else
                Call Load_Single_Row(row,rowblock,colblock,amp,dorder,mpointer, boundary = .true.)
            Endif
        Endif
    End Subroutine Load_BC


    Subroutine FEContinuity(eqind, mode, varind,row,dorder)
        Implicit None
        Integer, Intent(In) :: eqind, varind, mode, row, dorder
        Integer :: colblock, rowblock
        real*8, Pointer, Dimension(:,:) :: mpointer
        mpointer => equation_set(mode,eqind)%mpointer
        colblock = equation_set(mode,eqind)%colblock(varind)
        rowblock = equation_set(mode,eqind)%rowblock

        Call Cheby_Continuity(row,rowblock,colblock,dorder,mpointer)

    End Subroutine FEContinuity
    Subroutine Clear_Row(eqind, mode,row)
        Implicit None
        Integer, Intent(In) :: eqind, mode, row
        Integer :: rowblock
        rowblock = equation_set(mode,eqind)%rowblock
        equation_set(mode,eqind)%mpointer(rowblock+row,:) = 0.0d0

    End Subroutine Clear_Row



    Subroutine LU_Solve_Sparse(self)
        Implicit None
        Class(Equation) :: self

        Integer :: error, mtype, phase, mxfct, mnum, nrows
        Integer :: msglvl, nrhs
        Real*8, Allocatable :: faux_rhs(:,:), faux_x(:,:)
        self%mtype = 11
        phase  = 33      ! only solve
        msglvl = 0

        self%nrows = self%nlinks*ndim1    ! and these
        self%nrhs = size(self%rhs_pointer)/(self%nrows)  ! Add these into equation structure
        !Write(6,*)'scheck: ', shape(self%rhs), self%nrhs, self%nrows
        !Write(6,*)'scheck2: ', size(self%sparse_mat), size(self%sparse_ia), size(self%sparse_ja)
        !Write(6,*)'iparm: ', self%iparm(5)
        self%mxfct = 1
        self%mnum = 1
        self%iparm(8) = 1

        Allocate(faux_rhs(1:self%nrows,1:self%nrhs))
        Allocate(faux_x(1:self%nrows,1:self%nrhs))
        faux_rhs(:,:) = 0.0d0
        faux_x(:,:) = 0.0d0

        !Write(6,*)'scheck: ', self%nrows*self%nrhs, size(faux_rhs), self%nrows, self%nrhs
        !Write(6,*)'scheck2: ', size(self%rhs_pointer), size(self%rhs_pointer)/self%nrows
        !Write(6,*)'scheck3: ', shape(self%rhs_pointer)
        !CALL pardiso (self%pt, self%mxfct, self%mnum, self%mtype, phase, self%nrows, &
        !            & self%sparse_mat, self%sparse_ia, self%sparse_ja, &
         !           & self%idum, self%nrhs, self%iparm, msglvl, self%rhs_pointer, &
          !          & faux_x, error, self%dparm)


        IF (error .NE. 0) THEN
            WRITE(*,*) '33:  The following ERROR was detected: ', error
            STOP
        END IF
        DeAllocate(faux_rhs)
        DeAllocate(faux_x)
    End Subroutine LU_Solve_Sparse



    Subroutine LU_Decompose_Sparse(self)
        Implicit None
        Class(Equation) :: self

        Integer :: error, mtype, solver,phase, mxfct, mnum, nrows
        Integer :: msglvl, nrhs
        self%nrows = self%nlinks*ndim1
        self%nrhs = size(self%rhs_pointer)/self%nrows

        self%mtype = 11  !unsymmetric matrix
        solver = 0  !sparse direct method
        msglvl = 0      ! with statistical information
        self%mxfct  = 1      ! msglvl = 1 will generate output to the screen
        self%mnum   = 1

        If (.not. self%sparse_initialized) Then
            self%sparse_initialized = .true.
            self%iparm(:) = 0
            !Call pardisoinit(self%pt, self%mtype, solver, self%iparm, &
            !    & self%dparm, error)
            Write(6,*)'iparm(2) = ', self%iparm(2), self%iparm(28)
            IF (error .NE. 0) THEN
                IF (error.EQ.-10 ) WRITE(*,*) 'No license file found'
                IF (error.EQ.-11 ) WRITE(*,*) 'License is expired'
                IF (error.EQ.-12 ) WRITE(*,*) 'Wrong username or hostname'
                STOP
            ENDIF



            phase  = 11      ! only reordering and symbolic factorization

            !Write(6,*)'Checkcheck: ', self%iparm(3), self%iparm(6)
            self%iparm(3) = 1
            self%iparm(6) = 1
            !CALL pardiso (self%pt, self%mxfct, self%mnum, self%mtype, phase, self%nrows, &
            !            & self%sparse_mat, self%sparse_ia, self%sparse_ja, &
            !            & self%idum, self%nrhs, self%iparm, msglvl, self%ddum, &
            !            & self%ddum, error, self%dparm)
            !WRITE(*,*) 'Reordering completed ... '

            IF (error .NE. 0) THEN
                WRITE(*,*) '11: The following ERROR was detected: ', error
                STOP
            END IF
        Endif

        phase  = 22      ! only factorization

        !CALL pardiso (self%pt, self%mxfct, self%mnum, self%mtype, phase, self%nrows, &
        !            & self%sparse_mat, self%sparse_ia, self%sparse_ja, &
        !            & self%idum, self%nrhs, self%iparm, msglvl, self%ddum, &
        !            & self%ddum, error, self%dparm)
        IF (error .NE. 0) THEN
            WRITE(*,*) '22:  The following ERROR was detected: ', error
            STOP
        END IF

    End Subroutine LU_Decompose_Sparse



    ! These last two routines are just wrappers for lapack routines  - possibly unnecessary
    Subroutine LU_Decompose_full(mat, pvt)
        Real*8, Intent(InOut) :: mat(:,:)
        Integer, Intent(Inout) :: pvt(:)
        Integer :: n,info

        n = Size(mat,1)
        Call Dgetrf(n, n, mat, n, pvt, info)
    !  Write(6,*)'info : ',info,n
    End Subroutine LU_Decompose_full




    Subroutine LU_Decompose_band(mat, pvt)
      Real*8, Intent(InOut) :: mat(:,:)
      Integer, Intent(out) :: pvt(:)
      Integer :: n, info,ku,kl,lda

      n = Size(mat,2)

      lda = Size(mat,1)

      kl = (lda - 1)/3
      ku = kl

      Call Dgbtrf(n, n, kl, ku, mat, lda, pvt, info)

    End Subroutine LU_Decompose_band
    Subroutine LU_Solve_Band(mat, rhs, pvt, na, nb)
          Real*8,intent(in) :: mat(:,:)
          Real*8,Intent(inout) :: rhs(:,:,:)
          Integer, Intent(in) :: pvt(:)
          Integer, Optional :: na, nb
          Integer :: ma, mb, lda, ku, kl,info

          If (Present(na)) Then
             ma = na
          Else
             ma = Size(mat,2)
          End If

          If (Present(nb)) Then
             mb = nb
          Else
             mb = Size(rhs)/Size(rhs,1)
        !     mb = Size(rhs)/Size(rhs,3)

          End If

          lda = Size(mat,1)

          kl = (lda - 1)/3
          ku = kl
        !  Call dgbtrs('N', ma, kl, ku, mb, mat, lda, pvt, rhs, Size(rhs,3), info)

          Call dgbtrs('N', ma, kl, ku, mb, mat, lda, pvt, rhs, Size(rhs,1), info)


    End Subroutine LU_Solve_Band

    Subroutine LU_Solve_full(mat, rhs, pvt, na, nb)
        Real*8,intent(in) :: mat(:,:)
        Real*8,Intent(inout) :: rhs(:,:,:)
        Integer, Intent(in) :: pvt(:)
        Integer, Optional :: na, nb
        Integer :: ma, mb, info

        If (Present(na)) Then
            Write(6,*)'na specified: ', na, Size(rhs,1)
            ma = na
        Else
            ma = Size(mat,1)
        End If

        If (Present(nb)) Then
            mb = nb
            Write(6,*)'mb specified: ', mb, Size(rhs)/Size(rhs,1)
        Else
            mb = Size(rhs)/Size(rhs,1)
        End If

        Call dgetrs('N', ma, mb, mat, Size(mat,1), pvt, rhs, Size(rhs,1), info)


        !If(Present(nb)) Then
        If (info .ne. 0) Then
            Write(6,*)'Problem is solve!  info is ', info
        Endif

    End Subroutine LU_Solve_full

    Subroutine Band_Arrange(equ,mode)
            Integer, Intent(In) :: equ, mode
            Integer :: link,i,nupper
            If (equation_set(mode,equ)%solvefor) Then
                If (equ .eq. 1) Then
                    nupper = 3*cpgrid%max_npoly
                    !Call Band_Load_Single(mode,equ,11)
                     Call Band_Load_Single(mode,equ,nupper)

                Else
                    nupper = cpgrid%max_npoly
                    Call Band_Load_Single(mode,equ,nupper) ! 3
                Endif


            equation_set(mode,equ)%mpointer => equation_set(mode,equ)%lhs
            if (equation_set(mode,equ)%nlinks .gt. 1) Then
                Do i = 1, equation_set(mode,equ)%nlinks
                    link = equation_set(mode,equ)%links(i)
                    equation_set(mode,link)%mpointer => equation_set(mode,equ)%lhs
                Enddo
            Endif

            Endif
    End Subroutine Band_Arrange


    Subroutine Band_Arrange_RHS(k)
        Integer, Intent(In) :: k
        Integer :: nlinks, i, r,j, nrow,rind

        ! Take a normal linked RHS and arrange it so the variables are interleaved
        nlinks = equation_set(1,k)%nlinks
        nrow = ndim1*nlinks

        Do j = 1, n_modes_total
            Do i = 1, ndim2
                Do r = 1, nrow
                    temp_rhs(r,i,j) = equation_set(1,k)%rhs(r,i,j)
                Enddo
            Enddo
        Enddo
        Do r = 1, ndim1
            rind = nlinks*(r-1)+1
            Do i = 0, nlinks-1
                equation_set(1,k)%rhs(rind+i,:,:) = temp_rhs(r+i*ndim1,:,:)
            Enddo
        Enddo


    End Subroutine Band_Arrange_RHS

    Subroutine Band_ReArrange_RHS(k)
        Integer, Intent(In) :: k
        Integer :: nlinks, i, r, j, nrow,rind

        ! Take a normal linked RHS and arrange it so the variables are interleaved
        nlinks = equation_set(1,k)%nlinks
        nrow = ndim1*nlinks

        Do j = 1, n_modes_total
            Do i = 1, ndim2
                Do r = 1, nrow
                    temp_rhs(r,i,j) = equation_set(1,k)%rhs(r,i,j)
                Enddo
            Enddo
        Enddo
        Do r = 1, ndim1
            rind = nlinks*(r-1)+1
            Do i = 0, nlinks-1
                equation_set(1,k)%rhs(r+i*ndim1,:,:) = temp_rhs(rind+i,:,:)
            Enddo
        Enddo


    End Subroutine Band_ReArrange_RHS

    Subroutine Sparse_Load(eind,mind)
        Implicit None
        !Equation_set(j,k)%LHS  j is equation, k is mode
        ! Oct 5, 2016:  This version makes some assumptions concerning the structure of the matrix
        Integer, Intent(In) :: eind, mind  ! Equation and mode indices
        Integer :: i,ii,j, k, n, npoly, npoly2, sind1, sind2,mcind1,mcind2
        Integer :: c_offset, s_offset, rindex !column and sparse matrix offset & row index
        Integer :: nlinks, nsub, element_count, n_rows, v_offset

        Real*8, Allocatable :: sparse_mat(:)
        Integer, Allocatable :: ia(:), ja(:)

        If (equation_set(mind,eind)%solvefor) Then

        !Write(6,*)'Loading sparse matrix'

        nlinks = equation_set(mind,eind)%nlinks
        nsub = cpgrid%domain_count
        N_rows = ndim1*nlinks
        element_count = 0

        Do i = 1, nsub
            npoly = cpgrid%npoly(i)
            element_count = element_count + cpgrid%npoly(i)*cpgrid%npoly(i)
            if (i .ne. 1) element_count = element_count+cpgrid%npoly(i-1)
            if (i .ne. nsub) element_count = element_count+cpgrid%npoly(i+1)
        Enddo

        element_count = element_count*nlinks*nlinks
        ! A little clunky, but works for now (copy into eq_set arrays at end)
        If (.not. allocated(equation_set(mind,eind)%sparse_mat)) Then
            Allocate(equation_set(mind,eind)%sparse_mat(1:element_count))
        Endif
        If (.not. allocated(equation_set(mind,eind)%sparse_ja)) Then
            Allocate(equation_set(mind,eind)%sparse_ja(1:element_count))
        Endif
        If (.not. allocated(equation_set(mind,eind)%sparse_ia)) Then
            Allocate(equation_set(mind,eind)%sparse_ia(1:n_rows+1))
        Endif

        Allocate(sparse_mat(1:element_count))
        Allocate(ja(1:element_count))
        Allocate(ia(1:n_rows+1))

        sparse_mat(:) = 0.D0
        ja(:) = 0
        ia(:) = 0

        rindex = 1
        s_offset = 0
        !Write(6,*)'INSIDE: ', eind, mind, nlinks
        Do n = 1, nlinks
            c_offset = 0
            Do i = 1, nsub
                npoly = cpgrid%npoly(i)


                ia(rindex) = s_offset+1
                If (i .eq. 1) Then  ! Treat lower domain boundaries as appropriate

                    Do k = 1, nlinks  ! Iterate across columns associated with ALL variables
                        v_offset = (k-1)*ndim1
                        sind1 = s_offset+1
                        sind2 = s_offset+npoly
                        mcind1 = v_offset+1
                        mcind2 = v_offset+npoly
                        sparse_mat(sind1:sind2) = &
                            & Equation_set(mind,eind)%LHS(rindex,mcind1:mcind2)
                        Do j = 1,npoly
                            ja(s_offset+j) = j+v_offset
                        Enddo
                        s_offset = s_offset+npoly
                    Enddo
                Else
                    npoly2 = cpgrid%npoly(i-1)
                    Do k = 1, nlinks
                        v_offset = (k-1)*ndim1
                        sind1 = s_offset+1
                        sind2 = s_offset+npoly+npoly2

                        mcind1 = v_offset+c_offset-npoly2+1
                        mcind2 = v_offset+c_offset+npoly
                        sparse_mat(sind1:sind2) = &
                            & Equation_set(mind,eind)%LHS(rindex,mcind1:mcind2)

                        Do j = 1,npoly+npoly2
                            ja(s_offset+j) = mcind1-1+j
                        Enddo

                        s_offset = s_offset+npoly+npoly2
                    Enddo
                Endif
                rindex = rindex+1
                ia(rindex) = s_offset+1
                ! Next, load the interior rows
                Do k = 2, npoly-1
                    ia(rindex) = s_offset+1

                    Do ii =1, nlinks
                        v_offset = (ii-1)*ndim1
                        mcind1 = v_offset+c_offset+1
                        mcind2 = v_offset+c_offset+npoly
                        sparse_mat(s_offset+1:s_offset+npoly) = &
                            & Equation_set(mind,eind)%LHS(rindex,mcind1:mcind2)
                        Do j = 1,npoly
                            ja(s_offset+j) = mcind1-1+j
                        Enddo
                        s_offset = s_offset+npoly
                    Enddo
                    rindex = rindex+1
                Enddo

                ia(rindex) = s_offset+1
                If ( i .eq. nsub) Then  ! Treat upper domain boundaries as appropriate

                    Do k = 1, nlinks
                        v_offset = (k-1)*ndim1
                        sind1 = s_offset+1
                        sind2 = s_offset+npoly
                        mcind1 = v_offset+c_offset+1
                        mcind2 = v_offset+c_offset+npoly
                        sparse_mat(sind1:sind2) = &
                            & Equation_set(mind,eind)%LHS(rindex,mcind1:mcind2)
                        Do j = 1,npoly
                            ja(s_offset+j) = mcind1-1+j
                        Enddo

                        s_offset = s_offset+npoly
                    Enddo
                Else
                    npoly2 = cpgrid%npoly(i+1)
                    Do k = 1, nlinks
                        v_offset = (k-1)*ndim1
                        sind1 = s_offset+1
                        sind2 = s_offset+npoly+npoly2

                        mcind1 = v_offset+c_offset+1
                        mcind2 = v_offset+c_offset+npoly+npoly2

                        sparse_mat(sind1:sind2) = &
                            & Equation_set(mind,eind)%LHS(rindex,mcind1:mcind2)
                        Do j = 1,npoly+npoly2
                            ja(s_offset+j) = mcind1-1+j
                        Enddo

                        s_offset = s_offset+npoly+npoly2
                        Enddo
                Endif

                rindex = rindex+1
                c_offset = c_offset+npoly
            Enddo
        Enddo
        ia(rindex) = s_offset+1
        !Write(6,*)'ROW CHECK: ', rindex-1, n_rows
        if (s_offset .ne. element_count) Write(6,*)'ELEMENT COUNT INCONSISTENT', s_offset, element_count

        Equation_set(mind,eind)%sparse_mat(:) = sparse_mat(:)
        Equation_set(mind,eind)%sparse_ia(:) = ia(:)
        Equation_set(mind,eind)%sparse_ja(:) = ja(:)
        DeAllocate(sparse_mat,ja,ia)
        DeAllocate(Equation_Set(mind,eind)%LHS)
        !DeAllocate(equation_set(mind,eind)%sparse_mat)
        !DeAllocate(equation_set(mind,eind)%sparse_ja)
        !DeAllocate(equation_set(mind,eind)%sparse_ia)
        !Write(6,*)'Outside: ', eind, mind, nlinks
        Endif
    End Subroutine Sparse_Load



    Subroutine Band_Load_Single(j,k,n_upper)
          !Equation_set(j,k)%LHS  j is equation, k is mode
          Real*8, Allocatable :: band_matrix(:,:), temp_rows(:,:)
          Integer :: r, i, row_diag, rput, N_Rows,j,k
          Integer :: istart, iend, n_upper, n_lower,nlinks,rind
          nlinks = equation_set(j,k)%nlinks

            N_rows = ndim1*nlinks
          n_lower = n_upper
          row_diag = n_lower+n_upper+1
          Allocate(band_matrix(2*n_lower+n_upper+1, N_Rows))
          band_matrix(:,:) = 0.D0

            if (nlinks .gt. 1) Then
                ! If we want to band solve a linked equation, we have to
                ! intereave the rows and columns to make it banded
                Allocate(temp_rows(1:N_Rows,1:N_rows))

                Do r = 1, ndim1    ! rows
                    rind = nlinks*(r-1)+1
                    Do i = 0, nlinks-1
                        temp_rows(rind+i,:) = equation_set(j,k)%LHS(r+i*ndim1,:)
                    Enddo
                Enddo
                Equation_set(j,k)%LHS(:,:) = temp_rows(:,:)

                Do r = 1, ndim1    ! columns
                    rind = nlinks*(r-1)+1
                    Do i = 0, nlinks-1
                        temp_rows(:,rind+i) = equation_set(j,k)%LHS(:,r+i*ndim1)
                    Enddo
                Enddo
                Equation_set(j,k)%LHS(:,:) = temp_rows(:,:)

                DeAllocate(temp_rows)
            Endif

          Do r = 1, N_Rows
             istart = r-n_lower
             iend = r+n_upper
             If (istart .lt. 1) istart = 1
             If (iend .gt. N_Rows) iend = N_Rows
             Do i = istart, iend
                rput = row_diag+(r-i)
                band_matrix(rput ,i ) = Equation_Set(j,k)%LHS(r,i)
             Enddo
          Enddo
          DeAllocate(Equation_Set(j,k)%LHS)
          Allocate(Equation_Set(j,k)%LHS(2*n_lower+n_upper+1, N_Rows))
          Equation_Set(j,k)%LHS(:,:) = band_matrix(:,:)
          DeAllocate(band_matrix)
    End Subroutine Band_Load_Single


!Matrix row-loading routines
    Subroutine Cheby_Continuity(rind,row,col,dorder,mpointer) !, clear_row, boundary)
        Implicit None
        Integer, Intent(In) :: rind,row, col, dorder
        Integer :: n, offleft, offright, r
        Integer :: hh, ind, nsub
        real*8, Pointer, Dimension(:,:), Intent(InOut) :: mpointer

        nsub = cpgrid%domain_count


        r = cpgrid%npoly(1)         !Decide whether we are replacing rows at
        if (rind .eq. 1) r = r+1    !either the top or bottom of each domain

        offleft = 0                 !Column offset for left-side domain
        offright = cpgrid%npoly(1)  !Column offset for right-side domain
        !write(6,*)'rind check: ', rind, r
        Do hh = 1, nsub -1
            !write(6,*)'rcheck 2: ', r
            ! Clear this row completely
            mpointer(r+row,:) = 0.0d0

            ! Load the left-side domain first (using its "upper" boundary values)
            ind = cpgrid%npoly(hh)
            Do n = 1, cpgrid%rda(hh) -1     ! De-Alias at boundaries (single rows are really just for boundaries)
                mpointer(row+r,col+n+offleft) = mpointer(row+r,col+n+offleft) &
                    &  + cpgrid%dcheby(hh)%data(ind,n,dorder)
            Enddo

            !Load the right-side domain (using its "lower" boundary values)
            Do n = 1, cpgrid%rda(hh+1)-1
                mpointer(row+r,col+n+offright) = mpointer(row+r,col+n+offright) &
                    & - cpgrid%dcheby(hh+1)%data(1,n,dorder)
            Enddo

            !Advance row and offset indices to next subdomain


            r = r+cpgrid%npoly(hh+1)
            offleft = offleft+cpgrid%npoly(hh)
            offright = offright+cpgrid%npoly(hh+1)
        Enddo


    End Subroutine Cheby_Continuity


    Subroutine Load_Interior_Rows_Cheby(row,col,amp,dorder,mpointer)
        Implicit None
        Integer, Intent(In) :: row, col, dorder
        Integer :: r, n, off1,npoly, hh, i, nsub
        real*8, Intent(In) :: amp(:)
        real*8, Pointer, Dimension(:,:), Intent(In) :: mpointer

        r    = 1
        off1 = 0
        nsub = cpgrid%domain_count
        Do hh = 1, nsub
            npoly = cpgrid%npoly(hh)
            Do i = 1, npoly
                Do n = 1, npoly
                    mpointer(row+r,col+n+off1) = mpointer(row+r,col+n+off1) &
                        & +amp(r)*cpgrid%dcheby(hh)%data(i,n,dorder)
                Enddo
                r = r+1
            Enddo
            off1 = off1+npoly
        Enddo
        !write(6,*)'Rout: ', r-1

    End Subroutine Load_Interior_Rows_Cheby

    Subroutine Load_Single_Row_Cheby(r,row,col,amp,dorder,mpointer, clear_row, boundary)
        Implicit None
        Integer, Intent(In) :: r,row, col, dorder
        Integer :: n, off1, local_index, domain, nsub,hh,rupper
        real*8, Intent(In) :: amp
        real*8, Pointer, Dimension(:,:), Intent(InOut) :: mpointer
        Logical, Intent(In), Optional :: clear_row, boundary
        Logical :: bjunk
        If (present(clear_row)) Then
            ! clear everything in this row
            mpointer(r+row,:) = 0.0d0
        Endif
        If (present(boundary)) Then
            ! Do nothing at the moment
            bjunk = boundary    ! Placeholder to avoid Intel compiler warnings (unused vars)
        Endif

        !Before loading the row, establish which domain r corresponds to.

        nsub = cpgrid%domain_count
        domain = 1                  ! initial guess at the r's domain
        off1 = 0                    ! column offset associated with that domain
        rupper = cpgrid%npoly(1)    ! uppermost radial index that belongs to domain being considered
        local_index = r             ! local index within r's domain
        !Write(6,*)'rcheck : ', local_index
        Do hh = 1, nsub-1
            If (local_index .gt. rupper) Then
                off1 = off1+cpgrid%npoly(hh)
                domain = domain+1
                local_index = local_index-rupper
            Endif
            rupper = cpgrid%npoly(hh+1)
        Enddo
        !Write(6,*)'local_index, r, domain : ', local_index, r, domain
        Do n = 1, cpgrid%rda(domain)-1     ! De-Alias at boundaries (single rows are really just for boundaries)
            mpointer(row+r,col+n+off1) = mpointer(row+r,col+n+off1) &
                & + amp*cpgrid%dcheby(domain)%data(local_index,n,dorder)
        Enddo

    End Subroutine Load_Single_Row_Cheby



End Module Linear_Solve

