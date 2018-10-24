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

Module Sphere_Spectral_Space
    Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values, my_lm_min, my_nl_lm, my_nm_lm, my_lm_lval, my_lm_max
    Use Parallel_Framework
    Use Controls
    Use ProblemSize
    Use Finite_Difference, Only : d_by_dx
    Use Linear_Solve
    Use Fields
    Use BoundaryConditions
    Use ClockInfo
    Use Timers
    Use Sphere_Linear_Terms
    Implicit None
    Type(SphericalBuffer) :: ctemp ! workspace
Contains



    Subroutine Post_Solve()
        Implicit None
        Integer :: m, i
        Character*12 :: tstring, otstring

        ! wsp%p1b is assumed to be allocated
        Call StopWatch(psolve_time)%startclock()
        Call wsp%construct('p1a')
        wsp%config = 'p1a'

        old_deltat = deltat
        If (new_timestep) Then
            deltat = new_deltat
            new_timestep = .false.
            If (my_rank .eq. 0) Then
                Write(otstring,t_ofmt)old_deltat
                Write(tstring,t_ofmt)deltat
                Call stdout%print(' Timestep has changed from '//Trim(otstring)//' to '//Trim(tstring)//'.')
                Call stdout%partial_flush()  ! Make SURE that a changing timestep is recorded ...
                                             ! ... even at the expense of additional file I/O for redirected stdout

            Endif

            Call Reset_Linear_Equations()
        Endif

        if (iteration .eq. 1) then
                !Euler Step
                new_ab_factor = deltat
                old_ab_factor = 0.0d0
        else
                new_ab_factor = 0.5d0*deltat*(2 + deltat/old_deltat)
                old_ab_factor = -0.5d0*deltat**2/old_deltat
        endif

        wsp%p1b = wsp%p1b*old_ab_factor

        !Copy each variable out of the RHS into the top part of the buffer
        ! These variables are in spectral space radially
        !!!DDDD Write(6,*)'I am getting the new rhs: ', my_rank
        Call Get_All_RHS(wsp%p1a)



        Call gridcp%dealias_buffer(wsp%p1a)    ! de-alias


        ! This is terribly inefficient, but I just want to test the stability of Chebyshev vs. FD for not..
        ! We'll create a new buffer.  ctemp
        ! Store all the permanent derivatives there - in c space
        ctemp%nf1a = 4
        ctemp%nf1b = 4
        If (magnetism) then
            ctemp%nf1a = 5
            ctemp%nf1b = 5
        Endif
        Call ctemp%construct('p1a')
        ! W..
        Call gridcp%d_by_dr_cp(wvar,d3wdr3,wsp%p1a,3)
        ctemp%p1a(:,:,:,1) = wsp%p1a(:,:,:,d3wdr3)

        Call gridcp%d_by_dr_cp(wvar,dwdr   ,wsp%p1a,1)
        Call gridcp%d_by_dr_cp(wvar,d2wdr2 ,wsp%p1a,2)
        ! P....n
        Call gridcp%d_by_dr_cp(pvar,dpdr1,wsp%p1a,1)
        ctemp%p1a(:,:,:,2) = wsp%p1a(:,:,:,dpdr1)
        ! T
        Call gridcp%d_by_dr_cp(tvar,d2tdr2,wsp%p1a,2)
        ctemp%p1a(:,:,:,3) = wsp%p1a(:,:,:,d2tdr2)
        Call gridcp%d_by_dr_cp(tvar,dtdr,wsp%p1a,1)
        ! Z..
        Call gridcp%d_by_dr_cp(zvar,d2zdr2,wsp%p1a,2)
        ctemp%p1a(:,:,:,4) = wsp%p1a(:,:,:,d2zdr2)
        Call gridcp%d_by_dr_cp(zvar,dzdr,wsp%p1a,1)

        ! Magnetism
        If (magnetism) Then
            Call gridcp%d_by_dr_cp(avar,d2adr2,wsp%p1a,2)
            ctemp%p1a(:,:,:,5) = wsp%p1a(:,:,:,d2adr2)
            Call gridcp%d_by_dr_cp(avar,dadr  ,wsp%p1a,1)
            Call gridcp%d_by_dr_cp(cvar,dcdr  ,wsp%p1a,1)
            Call gridcp%d_by_dr_cp(cvar,d2cdr2,wsp%p1a,2)
        Endif

        !//////////////////////////////////////////////////////////////////////////
        ! Now everything we need is in the wsp or ctemp buffer
        ! The ctemp terms are those terms that do not leave this configuration
        ! transform them now & add them to appropriate equations
        Call ctemp%construct('p1b')
        Call gridcp%dealias_buffer(ctemp%p1a)    ! de-alias

        Call gridcp%From_Spectral(ctemp%p1a,ctemp%p1b)
        If (output_iteration) Then
            ! Grab dpdr
            Call cobuffer%construct('p1a')
            cobuffer%p1a(:,:,:,dpdr_cb) = ctemp%p1b(:,:,:,2)
        Endif


        Call Add_Derivative(peq,wvar,3,wsp%p1b,ctemp%p1b,1)
        Call Add_Derivative(weq,pvar,1,wsp%p1b,ctemp%p1b,2)
        Call Add_Derivative(teq,tvar,2,wsp%p1b,ctemp%p1b,3)
        Call Add_Derivative(zeq,zvar,2,wsp%p1b,ctemp%p1b,4)
        If (magnetism) Then
            Call Add_Derivative(aeq,avar,2,wsp%p1b,ctemp%p1b,5)
        Endif
        Call ctemp%deconstruct('p1a')
        Call ctemp%deconstruct('p1b')

        !//////////////////////////////////////////////
        !  Next, we reconstruct ctemp%p1a and copy wsp%p1a into it
        ctemp%nf1a = wsp%nf1a
        Call ctemp%construct('p1a')
        ctemp%p1a(:,:,:,:) = wsp%p1a(:,:,:,:)
        Call gridcp%dealias_buffer(ctemp%p1a) !De-Alias
        wsp%p1a(:,:,:,:) = 0.0d0    ! Shouldn't need to do this, but just to be sure
        Call gridcp%From_Spectral(ctemp%p1a,wsp%p1a)
        Call ctemp%deconstruct('p1a')

        !/////////////////////////////////////////////////////////////////
        !  The rest of the code can remain unchanged
        !/////////////////////////////////////////////////////////////////
        !Load the W derivatives into the appropriate RHS's

        Call Add_Derivative(peq,wvar,0,wsp%p1b,wsp%p1a,wvar)
        Call Add_Derivative(peq,wvar,1,wsp%p1b,wsp%p1a,dwdr)
        Call Add_Derivative(peq,wvar,2,wsp%p1b,wsp%p1a,d2wdr2)

        Call Add_Derivative(weq,wvar,0,wsp%p1b,wsp%p1a,wvar)
        Call Add_Derivative(weq,wvar,1,wsp%p1b,wsp%p1a,dwdr)
        Call Add_Derivative(weq,wvar,2,wsp%p1b,wsp%p1a,d2wdr2)



        !//////////////////////////////
        !  P Terms

        Call Add_Derivative(peq,pvar,0,wsp%p1b,wsp%p1a,pvar)


        !///////////////////////////////
        ! T Terms

        Call Add_Derivative(teq,tvar,1,wsp%p1b,wsp%p1a,dtdr)

        Call Add_Derivative(teq,tvar,0, wsp%p1b,wsp%p1a,tvar)
        Call Add_Derivative(weq,tvar,0, wsp%p1b,wsp%p1a,tvar)    ! gravity

        ! Convert temperature to temperature/r (will take derivatives of this for advection)
        ! As of April 30, 2017, we multiply by 1/r in physical space instead
        !Do m = 1, my_num_lm
        !    Do i = 1, 2
        !        wsp%p1a(:,i,m,tvar) = wsp%p1a(:,i,m,tvar)/radius(:)
        !    Enddo
        !Enddo


        !///////////////////////////////
        !  Z Terms


        Call Add_Derivative(zeq,zvar,0,wsp%p1b,wsp%p1a,zvar)
        Call Add_Derivative(zeq,zvar,1,wsp%p1b,wsp%p1a,dzdr)

        !///////////////////////////////////////
        !  Magnetic Terms
        If (magnetism) Then
            !//////////////
            ! A-terms (Toroidal magnetic field)

            Call Add_Derivative(aeq,avar,0,wsp%p1b,wsp%p1a,avar)

            !///////////////////
            ! C-terms (Poloidal magnetic field)

            Call Add_Derivative(ceq,cvar,2,wsp%p1b,wsp%p1a,d2cdr2)

            Call Add_Derivative(ceq,cvar,0,wsp%p1b,wsp%p1a,cvar)

        Endif


        !Load the old ab array into the RHS
        Call Set_All_RHS(wsp%p1b)    ! RHS now holds old_AB+CN factors


        Call wsp%deconstruct('p1b')
        Call StopWatch(psolve_time)%increment()

        Call StopWatch(ctranspose_time)%startclock()




        If (output_iteration) Then
            !Convert p/rho to p
            ! We already took d/dr(p/rho), so we'll fix that later
            Do m = 1, my_num_lm
                Do i = 1, 2
                    wsp%p1a(:,i,m,pvar) = wsp%p1a(:,i,m,pvar)*ref%density(:)
                Enddo
            Enddo
            Call cobuffer%reform()
        Endif
        Call wsp%reform()    ! move from p1a to s2a
        Call StopWatch(ctranspose_time)%increment()

    End Subroutine Post_Solve

    Subroutine AdvanceTime
        Implicit None
        ! wsp will be in 'p1b' config
        ! p1b contains the new adams bashforth term
        !Call print_max_spec2(pvar)

        if (.not. nonlinear) then
            wsp%p1b(:,:,:,:) = 0.0d0
        endif
        if (magnetism) then
            Call Finalize_EMF()
        endif
        Call Add_to_All_RHS(wsp%p1b,new_ab_factor)
        Call Enforce_Boundary_Conditions()
        Call StopWatch(solve_time)%startclock()
        Call Implicit_Solve()
        Call StopWatch(solve_time)%increment()
        simulation_time = simulation_time+deltat
        ! The righthand side of the equation set structure
        ! Now contains the updated fields.
    End Subroutine AdvanceTime
    Subroutine Finalize_EMF()
        Implicit None
        Integer m, i,j,jstart,jend
        ! we need to take one last radial derivative and combine terms

        If (chebyshev) Then
            ! Again, terribly inefficient, but we are looking to check the MHD right now.
            ! Will optimize this later.
            ctemp%nf1a = 2
            ctemp%nf1b = 2
            Call ctemp%construct('p1a')
            Call ctemp%construct('p1b')
            ctemp%p1a(:,:,:,:) = 0.0d0
            Do m = 1, my_num_lm
                Do i = 1, 2
                    ctemp%p1a(:,i,m,1) = wsp%p1b(:,i,m,emfphi)
                Enddo
            Enddo

            Call gridcp%to_spectral(ctemp%p1a,ctemp%p1b)
            Call gridcp%d_by_dr_cp(1,2,ctemp%p1b,1)
            Call gridcp%dealias_buffer(ctemp%p1b)
            Call gridcp%from_spectral(ctemp%p1b,ctemp%p1a)



            Do m = 1, my_num_lm
                Do i = 1, 2
                    wsp%p1b(:,i,m,avar) = wsp%p1b(:,i,m,avar) + ctemp%p1a(:,i,m,2)
                Enddo
            Enddo
            Call ctemp%deconstruct('p1a')
            Call ctemp%deconstruct('p1b')
        Else
            ctemp%nf1a = 1
            Call ctemp%construct('p1a')
            ctemp%p1a(:,:,:,:) = 0.0d0
            Do m = 1, my_num_lm
                Do i = 1, 2
                    ctemp%p1a(:,i,m,1) = wsp%p1b(:,i,m,avar)
                Enddo
            Enddo
            Call d_by_dx(emfphi,avar,wsp%p1b,1)

            Do m = 1, my_num_lm
                Do i = 1, 2
                    wsp%p1b(:,i,m,avar) = ctemp%p1a(:,i,m,1) + wsp%p1b(:,i,m,avar)
                Enddo
            Enddo
            Call ctemp%deconstruct('p1a')
        Endif

    End Subroutine Finalize_EMF

End Module Sphere_Spectral_Space
