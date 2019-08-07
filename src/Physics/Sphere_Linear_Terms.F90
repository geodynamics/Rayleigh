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

Module Sphere_Linear_Terms
    Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values, my_lm_min, my_nl_lm, my_nm_lm, my_lm_lval, my_lm_max
    Use Controls
    Use ProblemSize
    Use Linear_Solve
    Use Fields
    Use BoundaryConditions
    Use Timers
    Use ClockInfo
    Use PDE_Coefficients
    Use Math_Constants
    Implicit None
    Real*8, Allocatable :: Lconservation_weights(:)

Contains
    Subroutine Linear_Init()
        Implicit None
        Real*8 :: amp, T,arg
        Integer :: n, r, m, nm
        !Depending on process layout, some ranks may not participate in the solve
        If (my_num_lm .gt. 0) Then 
 
            Call Initialize_Linear_System()

            If (strict_L_conservation) Then

                Allocate(Lconservation_weights(1:N_R))
                Lconservation_weights(1:N_R) = 0.0d0
                nm = 0
                Do m = 1, gridcp%domain_count
                    Do n = 1, gridcp%npoly(m)
                        Do r = 1, gridcp%npoly(m)
                            T = gridcp%dcheby(m)%data(r,n,0)
                            Lconservation_weights(n+nm) = Lconservation_weights(n+nm) + radial_integral_weights(r+nm) * T
                        Enddo
                    Enddo
                    Lconservation_weights( nm+(2*gridcp%npoly(m))/3+1:nm+gridcp%npoly(m) ) = 0.0d0  ! De-Alias
                    nm = nm + gridcp%npoly(m)
                Enddo
            Endif
        Endif
    End Subroutine Linear_Init

    Subroutine Reset_Linear_Equations()
        Implicit None
        Real*8 :: rhs_factor, lhs_factor
        Call Reset_Equation_Coefficients()    ! Zero out all coefficients
        lhs_factor = -deltat*alpha_implicit    ! Crank Nicolson scheme - alpha = 0.5
        rhs_factor = deltat*(1.0d0-alpha_implicit)
        Call Set_Time_Factors(lhs_factor,rhs_factor)
        Call Load_Linear_Coefficients()
        Call LU_Decompose_Matrices()
    End Subroutine Reset_Linear_Equations


    Subroutine Initialize_Linear_System
        Implicit None
        Integer :: lp, l, nlinks
        Integer, Allocatable :: eq_links(:), var_links(:)
        Type(Cheby_Grid), Pointer :: gridpointer
        nullify(gridpointer)
        gridpointer => gridcp
        If (chebyshev) Call Use_Chebyshev(gridpointer)    ! Turns chebyshev mode to "on" for the linear solve

        Call Initialize_Equation_Set(n_equations,n_variables,N_R,my_nl_lm, my_nm_lm,2)

        Do lp = 1, my_nl_lm
            l = my_lm_lval(lp)
            If (l .eq.0) Then

                nlinks = 3
                Allocate(eq_links(1:3))
                Allocate(var_links(1:3))

                eq_links(1) = weq
                eq_links(2) = peq
                eq_links(3) = teq

                var_links(1) = wvar
                var_links(2) = pvar
                var_links(3) = tvar

                ! Not optimal (right now if variables are linked for one mode, they are linked for all).
                Call link_equations(eq_links, var_links,nlinks,lp)

                ! ell = 0 pressure equation is hydrostatic balance
                Call Initialize_Equation_Coefficients(weq, wvar, 0,lp)        ! identity matrix for W
                Call Initialize_Equation_Coefficients(peq, pvar, 1,lp)
                Call Initialize_Equation_Coefficients(peq, tvar, 0,lp)
                Call Initialize_Equation_Coefficients(teq, tvar, 2,lp)

                DeAllocate(eq_links)
                DeAllocate(var_links)
            Else
                ! W equation
                Call Initialize_Equation_Coefficients(weq,wvar,2,lp)
                Call Initialize_Equation_Coefficients(weq,pvar,1,lp)
                Call Initialize_Equation_Coefficients(weq,tvar,0,lp)


                ! P equation
                Call Initialize_Equation_Coefficients(peq,wvar, 3,lp)
                Call Initialize_Equation_Coefficients(peq,pvar, 0,lp)

                ! T equation
                Call Initialize_Equation_Coefficients(teq,tvar, 2,lp)

                If (advect_reference_state) Then
                    Call Initialize_Equation_Coefficients(teq,wvar, 0,lp)
                Endif

                ! Z equation
                Call Initialize_Equation_Coefficients(zeq,zvar, 2,lp)

                If (magnetism) Then
                    Call Initialize_Equation_Coefficients(ceq,cvar, 2,lp)
                    Call Initialize_Equation_Coefficients(aeq,avar, 2,lp)
                Endif

                nlinks = 3
                Allocate(eq_links(1:3))
                Allocate(var_links(1:3))

                eq_links(1) = weq
                eq_links(2) = peq
                eq_links(3) = teq

                var_links(1) = wvar
                var_links(2) = pvar
                var_links(3) = tvar
                Call link_equations(eq_links, var_links,nlinks,lp)



                DeAllocate(eq_links)
                DeAllocate(var_links)
            Endif
        Enddo
        Call Finalize_Equations()
        If (bandsolve) Call Use_BandSolve()
        If (sparsesolve) Call Use_SparseSolve()

    End Subroutine Initialize_Linear_System

    Subroutine Load_Linear_Coefficients()
        Implicit None

        Real*8, Allocatable :: H_Laplacian(:), amp(:)
        Integer :: l, lp
        Real*8 :: diff_factor,ell_term
        !rmin_norm
        diff_factor = 1.0d0 ! hyperdiffusion factor (if desired, 1.0d0 is equivalent to no hyperdiffusion)
        Allocate(amp(1:N_R))
        Allocate(H_Laplacian(1:N_R))
        Do lp = 1, my_nl_lm
            If (bandsolve) Call DeAllocate_LHS(lp)
            Call Allocate_LHS(lp)
            l = my_lm_lval(lp)

            If (hyperdiffusion) Then
                ell_term = ((l-1.0d0)/(l_max-1.0d0))**hyperdiffusion_beta
                diff_factor = 1.0d0+hyperdiffusion_alpha*ell_term
            Endif
            H_Laplacian = - l_l_plus1(l) * OneOverRSquared
            If (l .eq. 0) Then
                !====================================================
                !            Temperature Equation
                ! T only
                amp = 1.0d0
                Call add_implicit_term(teq,tvar, 0, amp,lp, static = .true.)    ! Time independent part

                !amp = 2.0d0/radius/Pr
                amp = 2.0d0/radius*kappa
                Call add_implicit_term(teq,tvar, 1, amp,lp)
                !amp = 1.0d0/Pr
                amp = 1.0d0*kappa
                Call add_implicit_term(teq,tvar, 2, amp,lp)

                ! Kappa,rho, T variation in radius
                amp = S_Diffusion_Coefs_1
                Call add_implicit_term(teq,tvar,1,amp,lp)

                !=======================================
                !   Hydrostatic balance
                !      This part of the equation is static (i.e. not time-evolving)

                amp = -ref%Buoyancy_Coeff
                Call add_implicit_term(peq, tvar, 0, amp,lp, static = .true.)            ! Gravity    --- Need LHS_Only Flag

                amp = ref%dpdr_W_term
                Call add_implicit_term(peq, pvar, 1, amp,lp, static = .true.)            ! dPdr     --- Here too


                ! Pressure equation is replaced with a statement that the ell=0 W is zero
                amp = 1.0d0
                Call add_implicit_term(weq, wvar, 0, amp,lp, static = .true.)            ! --- any maybe here

                !  If band solve, do redefinition here
            Else

                !==================================================
                !                Radial Momentum Equation



                ! Temperature
                amp = -ref%Buoyancy_Coeff/H_Laplacian
                Call add_implicit_term(weq, tvar, 0, amp,lp)            ! Gravity


                ! Pressure
                !amp = 1.0d0/(Ek*H_Laplacian)*ref%density        ! dPdr
                amp = ref%dpdr_W_term/H_Laplacian
                Call add_implicit_term(weq,pvar, 1, amp,lp)

                ! W

                If (inertia) Then
                    ! (from u_{t+1} in CN method -- no dt factor)
                    amp = 1.0d0
                    Call add_implicit_term(weq,wvar, 0, amp,lp,static = .true.)
                Endif

                !amp = H_Laplacian        ! Diffusion
                amp = H_Laplacian*nu*diff_factor
                Call add_implicit_term(weq,wvar, 0, amp,lp)
                !amp = 1.0d0
                amp = nu*diff_factor
                Call add_implicit_term(weq,wvar, 2, amp,lp)

                ! These two diffusion bits are different
                ! depending on variation of rho and nu
                amp = W_Diffusion_Coefs_0*diff_factor
                Call add_implicit_term(weq,wvar, 0, amp,lp)
                amp = W_Diffusion_Coefs_1*diff_factor
                Call add_implicit_term(weq,wvar, 1, amp,lp)

                !==================================================
                !                Pressure (dWdr) Equation

                ! Pressure
                !amp = -(1.0d0)/Ek*ref%density
                amp = ref%pressure_dwdr_term
                Call add_implicit_term(peq,pvar, 0, amp,lp)

                ! W
                If (inertia) Then
                    ! (from u_{t+1} in CN method -- no dt factor)
                    amp = 1.0d0
                    Call add_implicit_term(peq,wvar, 1, amp,lp, static = .true.)
                Endif

                !amp =-H_Laplacian*2.0d0/radius
                amp =-nu*H_Laplacian*2.0d0/radius*diff_factor
                Call add_implicit_term(peq,wvar, 0, amp,lp)
                !amp = H_Laplacian
                amp = H_Laplacian*nu*diff_factor
                Call add_implicit_term(peq,wvar, 1, amp,lp)
                !amp = 1.0d0
                amp = nu*diff_factor
                Call add_implicit_term(peq,wvar, 3, amp,lp)


                ! Again, these two bits depend on radial variation of rho and nu
                amp = dW_Diffusion_Coefs_0*H_Laplacian*diff_factor
                Call add_implicit_term(peq,wvar, 0, amp,lp)

                amp = dW_Diffusion_Coefs_1*diff_factor
                Call add_implicit_term(peq,wvar, 1, amp,lp)

                amp = dW_Diffusion_Coefs_2*diff_factor
                Call add_implicit_term(peq,wvar, 2, amp,lp)
                !====================================================
                !            Temperature Equation

                ! T
                amp = 1.0d0
                Call add_implicit_term(teq,tvar, 0, amp,lp, static = .true.)        ! Time independent term

                !amp = H_Laplacian/Pr        ! Diffusion
                amp = H_Laplacian*kappa*diff_factor
                Call add_implicit_term(teq,tvar, 0, amp,lp)
                !amp = 2.0d0/radius/Pr
                amp = 2.0d0/radius*kappa*diff_factor
                Call add_implicit_term(teq,tvar, 1, amp,lp)

                ! amp = 1.0d0/Pr
                amp = kappa*diff_factor
                Call add_implicit_term(teq,tvar, 2, amp,lp)

                ! Kappa,rho, T variation in radius
                amp = S_Diffusion_Coefs_1*diff_factor
                Call add_implicit_term(teq,tvar,1,amp,lp)

                !Reference State Advection (only do this if reference state is non-adiabatic)
                If (advect_reference_state) Then
                    amp = (H_Laplacian/ref%density)*ref%dsdr
                    Call add_implicit_term(teq,wvar,0,amp,lp)
                Endif

                !=====================================================
                !    Z Equation

                If (inertia) Then
                    ! (from u_{t+1} in CN method -- no dt factor)
                    amp = 1.0d0
                    Call add_implicit_term(zeq,zvar, 0, amp,lp, static = .true.)    ! Time-independent piece
                Endif

                !amp = H_Laplacian
                amp = H_Laplacian*nu*diff_factor
                Call add_implicit_term(zeq,zvar, 0, amp,lp)
                !amp = 1.0d0
                amp = nu
                Call add_implicit_term(zeq,zvar, 2, amp,lp)

                ! Variation of rho and nu
                amp = Z_Diffusion_Coefs_0*diff_factor
                Call add_implicit_term(zeq,zvar, 0, amp,lp)
                amp = Z_Diffusion_Coefs_1*diff_factor
                Call add_implicit_term(zeq,zvar, 1, amp,lp)
                If (magnetism) Then
                    !=========================================
                    !  Btor Equation
                    amp = 1.0d0
                    Call add_implicit_term(aeq,avar, 0, amp,lp, static = .true.)    ! Time-independent piece

                    amp = H_Laplacian*eta*diff_factor
                    Call add_implicit_term(aeq,avar, 0, amp,lp)

                    amp = 1.0d0*eta*diff_factor
                    Call add_implicit_term(aeq,avar, 2, amp,lp)

                    ! Eta variation in radius
                    amp = A_Diffusion_Coefs_1*diff_factor
                    Call add_implicit_term(aeq,avar,1,amp,lp)

                    !=========================================
                    !  Bpol Equation
                    amp = 1.0d0
                    Call add_implicit_term(ceq,cvar, 0, amp,lp, static = .true.)    ! Time-independent piece

                    amp = H_Laplacian*eta*diff_factor
                    Call add_implicit_term(ceq,cvar, 0, amp,lp)

                    amp = 1.0d0*eta*diff_factor
                    Call add_implicit_term(ceq,cvar, 2, amp,lp)
                Endif

                ! If band solve, do the redefinition of the matrix here

            Endif
            Call Set_Boundary_Conditions(lp)
            If (sparsesolve) Then
                !Write(6,*)'matrix: ', weq,lp, my_rank, l
                Call Sparse_Load(weq,lp)
                !Write(6,*)'matrix: ', zeq,lp,my_rank, l
                Call Sparse_Load(zeq,lp)
                If (magnetism) Then
                    Call Sparse_Load(aeq,lp)
                    Call Sparse_Load(ceq,lp)
                Endif
            Endif


            If (bandsolve) Then
                Call Band_Arrange(weq,lp)
                Call Band_Arrange(zeq,lp)
                If (magnetism) Then
                    Call Band_Arrange(aeq,lp)
                    Call Band_Arrange(ceq,lp)
                Endif
            Endif

        Enddo
        DeAllocate(amp)
        DeAllocate(H_Laplacian)
    End Subroutine Load_Linear_Coefficients

    Subroutine Set_Boundary_Conditions(mode_ind)
        ! Modified version of set_boundary_conditions
        ! Designed to work with more memory friendly logic
        ! Sets boundary condition of indicated l-value (index lp)
        ! only.  Does not loop over lp.
        Implicit None
        Real*8 :: samp,one
        Integer, Intent(In) :: mode_ind
        Integer :: l, r,lp
        one = 1.0d0
        lp = mode_ind

        l = my_lm_lval(lp)

        If (l .eq. 0) Then
            Call Clear_Row(peq,lp,1)            ! Pressure only has one boundary condition
            Call Clear_Row(teq,lp,1)
            Call Clear_Row(teq,lp,N_R)

            ! "1" denotes linking at index 1, starting in domain 2
            ! "2" denotes linking at index npoly, starting in domain 1
            Call FEContinuity(peq,lp,pvar,1,0) ! Pressure is continuous
            Call FEContinuity(teq,lp,tvar,1,1) ! T' / S' is continuous (for ell = 0)
            Call FEContinuity(teq,lp,tvar,2,0) ! T/S is continuous

            !*******************************************************
            ! Entropy Boundary Conditions

            ! Temperature Boundary Conditions (T fixed bottom and top)
            r = 1
            If (fix_tvar_top) Then
                Call Load_BC(lp,r,teq,tvar,one,0)    !upper boundary
            Endif
            If (fix_dtdr_top) Then
                Call Load_BC(lp,r,teq,tvar,one,1)
            Endif

            r = N_R
            If (fix_tvar_bottom) Then
                Call Load_BC(lp,r,teq,tvar,one,0)    ! lower boundary
            Endif
            If (fix_dtdr_bottom) Then
                Call Load_BC(lp,r,teq,tvar,one,1)
            Endif


            ! The ell=0 pressure is really a diagnostic of the system.
            ! It doesn't drive anything.  The simplist boundary condition
            ! is to enforce a pressure node at the top.
            r = 1
            Call Load_BC(lp,r,peq,pvar,one,0)


        Else

            !*******************************************************
            !        Clear the boundary rows
            Call Clear_Row(weq,lp,1)
            Call Clear_Row(weq,lp,N_R)
            Call Clear_Row(peq,lp,1)
            Call Clear_Row(peq,lp,N_R)
            Call Clear_Row(teq,lp,1)
            Call Clear_Row(teq,lp,N_R)
            Call Clear_Row(zeq,lp,1)
            Call Clear_Row(zeq,lp,N_R)

            ! "1" denotes linking at index 1, starting in domain 2
            ! "2" denotes linking at index npoly, starting in domain 1
            Call FEContinuity(peq,lp,pvar,2,0)   ! Pressure is continuous
            Call FEContinuity(peq,lp,wvar,1,2)          ! W'' is continuous

            Call FEContinuity(weq,lp,wvar,2,0)   ! W is continuous
            Call FEContinuity(weq,lp,wvar,1,1)          ! W' is continuous

            Call FEContinuity(teq,lp,tvar,2,0)   ! T/S is continuous
            Call FEContinuity(teq,lp,tvar,1,1)          ! T' / S' is continuous (for ell = 0)


            Call FEContinuity(zeq,lp,zvar,2,0)   ! Z is continuous
            Call FEContinuity(zeq,lp,zvar,1,1)          ! Z' is continuous


            !*******************************************************
            ! Entropy Boundary Conditions

            ! Temperature Boundary Conditions (T fixed bottom and top)
            r = 1
            If (fix_tvar_top) Then
                If (.not. fix_divrfc_top) Then
                    Call Load_BC(lp,r,teq,tvar,one,0)    !upper boundary
                Endif
            Endif
            If (fix_dtdr_top) Then
                Call Load_BC(lp,r,teq,tvar,one,1)
            Endif

            r = N_R
            If (fix_tvar_bottom) Then
                Call Load_BC(lp,r,teq,tvar,one,0)    ! lower boundary
            Endif
            If (fix_dtdr_bottom) Then
                Call Load_BC(lp,r,teq,tvar,one,1)
            Endif

            !///////////////////////////////////////////////////////////////////
            ! These four different boundary conditions are similar, though
            ! slightly different, in nature.  The idea is to try some boundary conditions
            ! that allow entropy and it's derivatives to vary on the boundary
            ! Either Del dot Grad S, Del dot F_conductive, or Del_r dot Grad S or F_conductive is zero

            If (fix_divrt_top) Then
                r = 1
                !d2sdr2 + 2/r dsdr = 0 at r = r_top
                Call Load_BC(lp,r,teq,tvar,one,2)
                samp = 2.0d0/radius(r)
                Call Load_BC(lp,r,teq,tvar,samp,1)
            Endif

            If (fix_divt_top) Then
                r = 1
                !d2sdr2 + 2/r dsdr -l(l+1)/r^2= 0 at r = r_top
                Call Load_BC(lp,r,teq,tvar,one,2)
                samp = 2.0d0/radius(r)
                Call Load_BC(lp,r,teq,tvar,samp,1)
                samp = -l*1.0d0*(l*1.0d0+1)/radius(r)**2
                Call Load_BC(lp,r,teq,tvar,samp,0)
            Endif


            If (fix_divrFc_top) Then
                r = 1
                !d2sdr2 + 2/r dsdr = 0 at r = r_top
                Call Load_BC(lp,r,teq,tvar,one,2)
                samp = 2.0d0/radius(r)+(dlnkappa(r)+ref%dlnrho(r)+ref%dlnT(r))
                Call Load_BC(lp,r,teq,tvar,samp,1)
            Endif
            If (fix_divFc_top) Then
                r = 1
                if ( l .gt. 1) Then
                    !d2sdr2 + 2/r dsdr -l(l+1)/r^2= 0 at r = r_top
                    Call Load_BC(lp,r,teq,tvar,one,2)
                    samp = 2.0d0/radius(r)+(dlnkappa(r)+ref%dlnrho(r)+ref%dlnT(r))
                    Call Load_BC(lp,r,teq,tvar,samp,1)
                    samp = -l*1.0d0*(l*1.0d0+1)/radius(r)**2
                    Call Load_BC(lp,r,teq,tvar,samp,0)
                Else
                    ! This boundary condition generally works well, but it can develop an ell=1
                    ! convection mode that causes undesired (and probably unphysical) behavior
                    samp = 1.0d0
                    Call Load_BC(lp,r,teq,tvar,samp,0)
                Endif
            Endif


            !************************************************************
            ! Velocity Boundary Conditions

            ! Impenetrable top and bottom
            ! W vanishes at the boundaries
            r = 1
            Call Load_BC(lp,r,weq,wvar,one,0)
            r = N_R
            Call Load_BC(lp,r,weq,wvar,one,0)


            If (no_slip_top) Then
                r = 1

                Call Load_BC(lp,r,zeq,zvar,one,0)
                Call Load_BC(lp,r,peq,wvar,one,1)
            Else
                ! Else stress-free
                r = 1
                samp = -(2.0d0/radius(r)+ref%dlnrho(r))
                Call Load_BC(lp,r,peq,wvar,one,2)
                Call Load_BC(lp,r,peq,wvar,samp,1)


                Call Load_BC(lp,r,zeq,zvar,one,1)
                Call Load_BC(lp,r,zeq,zvar,samp,0)
            Endif

            If (no_slip_bottom) Then
                r = N_R
                Call Load_BC(lp,r,peq,wvar,one,1)
                Call Load_BC(lp,r,zeq,zvar,one,0)

            Else
                !stress_free_bottom
                r = N_R
                samp = -(2.0d0/radius(r)+ref%dlnrho(r))
                Call Load_BC(lp,r,peq,wvar,one,2)
                Call Load_BC(lp,r,peq,wvar,samp,1)
                Call Load_BC(lp,r,zeq,zvar,one,1)
                Call Load_BC(lp,r,zeq,zvar,samp,0)
            Endif

            If ((l .eq. 1) .and. (strict_L_Conservation) ) then
                Call Clear_Row(zeq,lp,1)
                Call Load_BC(lp,1,zeq,zvar,one,0,integral = Lconservation_weights)
            Endif


            !*******************************************************
            !        Magnetic Boundary Conditions

            If (Magnetism) Then
                !  Clear the boundary rows
                Call Clear_Row(ceq,lp,1)
                Call Clear_Row(ceq,lp,N_R)
                Call Clear_Row(aeq,lp,1)
                Call Clear_Row(aeq,lp,N_R)


                Call FEContinuity(aeq,lp,avar,2,0)   ! A is continuous
                Call FEContinuity(aeq,lp,avar,1,1)          ! A' is continuous

                Call FEContinuity(ceq,lp,cvar,2,0)   ! C is continuous
                Call FEContinuity(ceq,lp,cvar,1,1)          ! C' is continuous


                ! Match to a potential field at top and bottom
                ! Btor = 0 at top and bottom
                r = 1
                Call Load_BC(lp,r,aeq,avar,one,0)
                r = N_R
                Call Load_BC(lp,r,aeq,avar,one,0)

                ! dBpol/dr+ell*Bpol/r = 0 at outer boundary
                r = 1
                Call Load_BC(lp,r,ceq,cvar,one,1)
                samp = my_lm_lval(lp)*one_over_r(r)
                Call Load_BC(lp,r,ceq,cvar,samp,0)

                ! dBpol/dr-(ell+1)*Bpol/r = 0 at inner boundary
                r = N_R
                Call Load_BC(lp,r,ceq,cvar,one,1)
                samp = - (l+1)*One_Over_R(r)
                Call Load_BC(lp,r,ceq,cvar,samp,0)


                If (fix_poloidalfield_top) Then
                    Call Clear_Row(ceq,lp,1)
                    Call Clear_Row(aeq,lp,1)
                    r = 1
                    Call Load_BC(lp,r,ceq,cvar,one,0)
                    Call Load_BC(lp,r,aeq,avar,one,0)
                Endif
                If (fix_poloidalfield_bottom) Then
                    Call Clear_Row(aeq,lp,N_R)
                    Call Clear_Row(ceq,lp,N_R)
                    r = N_R
                    Call Load_BC(lp,r,ceq,cvar,one,0)
                    Call Load_BC(lp,r,aeq,avar,one,0)
                Endif

            Endif    ! Magnetism


        Endif ! l = 0 or not


    End Subroutine Set_Boundary_Conditions

    Subroutine Enforce_Boundary_Conditions
        Implicit None
        Real*8 :: bc_val
        Integer :: uind, lind
        Integer :: real_ind, imag_ind

        Call Apply_Boundary_Mask(bc_values)
        Call Domain_Continuity()

    End Subroutine Enforce_Boundary_Conditions

    Subroutine Domain_Continuity()
        Implicit None
        Integer :: l, indx, ii,lp, j, k,n
        ! start applying the boundary and continuity conditions by setting
        ! the appropriate right hand sides.

        ! Will wrap this more into Linear_Solve.F90 soon
        ii = 2*N_R


        indx = 1

        Do lp = 1, my_nl_lm
            n_m = my_nm_lm(lp)-1    ! really n_m, but the indexing below is from the old implicit solv
            l = my_lm_lval(lp)

            If (chebyshev) Then
                ! Apply continuity conditions across the subdomains

                j = 0
                do n = 1, ndomains-1
                    j = j+gridcp%npoly(n)
                    equation_set(1,weq)%RHS(      j, : , indx:indx+n_m) = 0.0d0
                    equation_set(1,weq)%RHS(2*N_R+j, : , indx:indx+n_m) = 0.0d0
                    if (l .ne. 0) then
                        equation_set(1,zeq)%RHS(      j, : , indx:indx+n_m) = 0.0d0
                        equation_set(1,weq)%RHS(N_R  +j, : , indx:indx+n_m) = 0.0d0 ! dpdr continuity - only for ell =/ 0
                    endif
                Enddo


                j = 1
                Do n = 1, ndomains-1
                    j = j+gridcp%npoly(n)
                    equation_set(1,weq)%RHS(      j  ,:, indx:indx+n_m) = 0.0d0
                    equation_set(1,weq)%RHS(  N_R+j  ,:, indx:indx+n_m) = 0.0d0
                    equation_set(1,weq)%RHS(2*N_R+j  ,:, indx:indx+n_m) = 0.0d0
                    if (l .ne. 0) equation_set(1,zeq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                Enddo

                If (Magnetism) Then
                    if (l .ne. 0) then
                        j = 0
                        Do n = 1, ndomains-1
                            j = j+gridcp%npoly(n)
                            equation_set(1,aeq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                            equation_set(1,ceq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                        Enddo

                        j = 1
                        Do n = 1, ndomains-1
                            j = j+gridcp%npoly(n)
                            equation_set(1,aeq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                            equation_set(1,ceq)%RHS(      j,: , indx:indx+n_m) = 0.0d0
                        Enddo
                    endif
                Endif


            Endif

            indx = indx + n_m + 1
        Enddo
    End Subroutine Domain_Continuity

End Module Sphere_Linear_Terms
