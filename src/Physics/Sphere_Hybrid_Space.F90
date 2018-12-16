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

#include "indices.F"
Module Sphere_Hybrid_Space

    ! NOTE: WE NEED a 1/density variable
    Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values, my_lm_min, my_nl_lm, my_nm_lm, my_lm_lval, my_lm_max
    Use Parallel_Framework
    Use Controls
    Use ProblemSize
    Use Legendre_Polynomials, Only : p_lm_array
    Use Legendre_Transforms, Only : Legendre_Transform
    Use Spectral_Derivatives
    Use Fields
    Use Timers
    Use ClockInfo
    Use ReferenceState

    Implicit None
    Real*8, Allocatable :: over_rhor(:), over_rhorsq(:), drho_term(:)

    Type(rmcontainer3D), Allocatable :: ftemp1(:), ftemp2(:),ftemp3(:), ftemp4(:)
Contains


    Subroutine Hybrid_Init()
        Integer :: r1, r2
        ! Allocate a few useful arrays that prevent extra mult/adds
        Allocate(over_rhor(my_r%min:my_r%max))
        Allocate(over_rhorsq(my_r%min:my_r%max))
        Allocate(drho_term(my_r%min:my_r%max))
        r1 = my_r%min
        r2 = my_r%max
        over_rhor(r1:r2) = one_over_r(r1:r2)/ref%density(r1:r2)
        over_rhorsq(r1:r2) = OneOverRSquared(r1:r2)/ref%density(r1:r2)
        drho_term(r1:r2) = ref%dlnrho(r1:r2)+one_over_r(r1:r2)

    End Subroutine Hybrid_Init

    Subroutine rlm_spacea()
        Implicit None
        Integer :: mp,r,imi,m,l, ind_top
        Call StopWatch(rlma_time)%startclock()

        ! Zero out l_max mode
        Do mp = my_mp%min, my_mp%max
            SBUFFA(l_max,:,:,:) = 0.0d0
        Enddo


        ! Allocate two work arrays
        Call Allocate_rlm_Field(ftemp1)
        Call Allocate_rlm_Field(ftemp2)

        If (output_iteration) Call Hybrid_Output_Initial()

        Call Velocity_Components()
        Call Velocity_Derivatives()
        Call d_by_dtheta(wsp%s2a,tvar,dtdt)


        If (magnetism) Call compute_BandCurlB()

        If (output_iteration) Call Hybrid_Output_Final()



        Call DeAllocate_rlm_Field(ftemp1)
        Call DeAllocate_rlm_Field(ftemp2)

        ! Zero out l_max mode
        Do mp = my_mp%min, my_mp%max
            wsp%s2a(mp)%data(l_max,:,:,:) = 0.0d0
        Enddo

        Call StopWatch(rlma_time)%increment()

        ! Legendre Transform and transpose the buffer
        Call wsp%construct('p2a')
        Call StopWatch(legendre_time)%startclock()
        Call Legendre_Transform(wsp%s2a,wsp%p2a)
        Call StopWatch(legendre_time)%increment()
        Call wsp%deconstruct('s2a')
        wsp%config = 'p2a'

        Call StopWatch(rtranspose_time)%startclock()

        If (output_iteration) Then
            Call wsp%reform(nextra_recv = output_nextra)
        Else
            Call wsp%reform()    ! We are now in p3a
        Endif

        Call StopWatch(rtranspose_time)%increment()
    End Subroutine rlm_spacea

    Subroutine rlm_spaceb()
        Implicit None
        Integer :: m, mp, r, imi
        ! Upon entry into this routine, we have the following quantities
        ! Tvar : RHS for the T equation
        ! Wvar : l(l+1)*RHS for the W equation
        ! Pvar : r/sintheta * [u dot grad u]_theta
        ! Zvar : r/sintheta * [u dot grad u]_phi

        ! The RHS for T is ready to go
        ! The W, Z and dWdr RHS's need a little work

        ! Transform
        Call wsp%construct('s2b')

        Call StopWatch(legendre_time)%startclock()
        Call Legendre_Transform(wsp%p2b,wsp%s2b)
        Call StopWatch(legendre_time)%increment()

        Call wsp%deconstruct('p2b')
        wsp%config = 's2b'
        Call StopWatch(rlmb_time)%startclock()



        ! The NL RHS for W is r^2/(l(l+1)) * the NL RHS for Ur
        ! We already have the r^2 taken care of.  Now for the l(l+1)

        DO_IDX2
            SBUFFB(IDX2,wvar) = SBUFFB(IDX2,wvar)*over_l_l_plus1(m:l_max)
        END_DO

        ! Now for the Z RHS, formed from the radial component of the curl of u dot grad u

        Call Allocate_rlm_Field(ftemp1)
        Call Allocate_rlm_Field(ftemp2)

        Call d_by_sdtheta(wsp%s2b, zvar,ftemp1)    ! need to be sure we have this indexing correct
        Call d_by_dphi(wsp%s2b,pvar,ftemp2)

        DO_IDX2
            ftemp1(mp)%data(IDX2) = ( ftemp2(mp)%data(IDX2)- &
                & ftemp1(mp)%data(IDX2) )*over_l_l_plus1(m:l_max)
        END_DO



        Call d_by_dphi(wsp%s2b,zvar,ftemp2)
        Call d_by_sdtheta(wsp%s2b,pvar,zvar)

        DO_IDX2
            SBUFFB(IDX2,pvar) = ( SBUFFB(IDX2,zvar)+ &
                & ftemp2(mp)%data(IDX2) )*over_l_l_plus1(m:l_max)
        END_DO
        ! dwdr RHS (p equation) is now loaded


        DO_IDX2
            SBUFFB(IDX2,zvar) = ftemp1(mp)%data(IDX2)
        END_DO
        ! Z RHS is now loaded



        ! The ell =0 w and p and z equations have zero RHS
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            if (m .eq. 0) then
                SBUFFB(0,my_r%min:my_r%max,1:2, pvar) = 0.0d0
                SBUFFB(0,my_r%min:my_r%max,1:2, wvar) = 0.0d0
                SBUFFB(0,my_r%min:my_r%max,1:2, zvar) = 0.0d0
            endif
        Enddo


        If (magnetism) Call adjust_emf()

        Call DeAllocate_rlm_Field(ftemp1)
        Call DeAllocate_rlm_Field(ftemp2)


        ! Zero out l_max mode
        Do mp = my_mp%min, my_mp%max
            SBUFFB(l_max,:,:,:) = 0.0d0
        Enddo

        Call StopWatch(rlmb_time)%increment()

        Call StopWatch(ctranspose_time)%startclock()
        Call wsp%reform() ! move to the solve space
        Call StopWatch(ctranspose_time)%increment()

        Call Adjust_TimeStep()


    End Subroutine rlm_spaceb

    Subroutine Hydro_Output_Derivatives()
        Implicit None
        Integer :: r, l, m, mp, imi
        ! Compute sin(theta) dP/dtheta and
        ! place it in the cobuffer
        Call d_by_dtheta(wsp%s2a,pvar,ftemp1)
        DO_IDX2
            ASBUFFA(IDX2,dpdt_cb) = ftemp1(mp)%data(IDX2)
        END_DO
    End Subroutine Hydro_Output_Derivatives

    Subroutine Velocity_Components()
        Implicit None
        Integer ::  m, mp,  r, imi


        ! Compute the velocity vield

        ! vr    overwrites w
        DO_IDX2
            SBUFFA(IDX2,vr) = l_l_plus1(m:l_max)*SBUFFA(IDX2,vr)*Over_RhoRSQ(r)
        END_DO

        ! We compute sintheta v_theta
        Call d_by_dtheta(wsp%s2a,dwdr,ftemp1)
        Call d_by_dphi(wsp%s2a,zvar,    ftemp2)

        DO_IDX2
            ftemp1(mp)%data(IDX2) = ftemp1(mp)%data(IDX2)+ftemp2(mp)%data(IDX2)
        END_DO

        DO_IDX2
                SBUFFA(IDX2,vtheta) = ftemp1(mp)%data(IDX2)*Over_RhoR(r)
        END_DO

        ! Now sintheta v_phi
        Call   d_by_dphi(wsp%s2a,dwdr,    ftemp1)
        Call d_by_dtheta(wsp%s2a,zvar,ftemp2)
        DO_IDX2
            ftemp1(mp)%data(IDX2) = ftemp1(mp)%data(IDX2)-ftemp2(mp)%data(IDX2)
        END_DO

        DO_IDX2
            SBUFFA(IDX2,vphi) = ftemp1(mp)%data(IDX2)*Over_RhoR(r)
        END_DO

    End Subroutine Velocity_Components


    Subroutine Velocity_Derivatives()
        Implicit None
        Integer :: r, l, m, mp, imi
        !/////////////////////////////////
        ! sintheta dv theta dr
        Call d_by_dtheta(wsp%s2a,d2wdr2,ftemp1)    ! Store sintheta dwdtheta there for now.  We're going to use it a bit anyway.
        Call d_by_dphi(wsp%s2a,dzdr,    ftemp2)       ! Will overwrite this with dTdtheta shortly


        DO_IDX2
            ftemp1(mp)%data(IDX2) = ftemp1(mp)%data(IDX2)+ftemp2(mp)%data(IDX2)
        END_DO

        DO_IDX2
            SBUFFA(IDX2,dvtdr) = ftemp1(mp)%data(IDX2)*Over_RhoR(r)
        END_DO

        ! .... Small correction for density variation  :  - u_theta*dlnrhodr (added -u_theta/r as well here)
        ! Notice that there is a -u_theta/r term above.  These should be combined
        ! for efficiency later

        DO_IDX2
            SBUFFA(IDX2,dvtdr) = SBUFFA(IDX2,dvtdr)- &
                & SBUFFA(IDX2,vtheta)*drho_term(r)
        END_DO

        !/////////////////////////////////
        ! sinphi dv phi dr
        Call d_by_dphi(wsp%s2a,d2wdr2,ftemp1)    ! Store sintheta dwdtheta there for now.  We're going to use it a bit anyway.
        Call d_by_dtheta(wsp%s2a,dzdr,    ftemp2)       ! Will overwrite this with dTdtheta shortly

        DO_IDX2
            ftemp1(mp)%data(IDX2) = ftemp1(mp)%data(IDX2)-ftemp2(mp)%data(IDX2)
        END_DO

        DO_IDX2
            SBUFFA(IDX2,dvpdr) = ftemp1(mp)%data(IDX2)*Over_RhoR(r)
        END_DO

        ! .... Small correction for density variation  :  - u_phi*dlnrhodr
        ! .... moved -u_phi/r here as well
        DO_IDX2
            SBUFFA(IDX2,dvpdr) = SBUFFA(IDX2,dvpdr)- &
                &  SBUFFA(IDX2,vphi)*drho_term(r)
        END_DO
        !/////////////////////////////////////////
        ! dvrdr    overwrites dwdr

        DO_IDX2
            SBUFFA(IDX2,dvrdr) = l_l_plus1(m:l_max)* &
                & SBUFFA(IDX2,dvrdr)*Over_RhoRSQ(r)
        END_DO


        DO_IDX2
            SBUFFA(IDX2,dvrdr) = SBUFFA(IDX2,dvrdr)- &
                & SBUFFA(IDX2,vr)*Two_Over_R(r)
        END_DO

        ! .... Small correction for density variation  :  - u_r*dlnrhodr
        DO_IDX2
            SBUFFA(IDX2,dvrdr) = SBUFFA(IDX2,dvrdr)- &
                & SBUFFA(IDX2,vr)*ref%dlnrho(r)
        END_DO


        Call d_by_dtheta(wsp%s2a,vr,dvrdt)


        ! Convert Z to ell(ell+1) Z/r^2  (i.e. omega_r)
        DO_IDX2
            SBUFFA(IDX2,zvar) = l_l_plus1(m:l_max)*SBUFFA(IDX2,zvar)*Over_RhoRSQ(r)
        END_DO
    End Subroutine Velocity_Derivatives

    Subroutine Compute_BandCurlB()
        Implicit None
        Integer :: imi, m, mp, r

        ! This routine computes B and Del X B


        !/////////////// BR /////////////////////
        ! First convert C to Br  ! Br overwrites C
        DO_IDX2
            SBUFFA(IDX2,Br) = l_l_plus1(m:l_max)*SBUFFA(IDX2,Br)*OneOverRSquared(r)
        END_DO

        !////////////////// [Del x B]_r ///////////////////////////
        ! (does not overwrite any existing fields)
        DO_IDX2
            SBUFFA(IDX2,curlbr) = l_l_plus1(m:l_max) &
               *SBUFFA(IDX2,Avar)*OneOverRSquared(r)
        END_DO

        ! Convert d2cdr2 to d2cdr2-Br (br = cl(l+1)/r^2
        DO_IDX2
            SBUFFA(IDX2,d2cdr2) = SBUFFA(IDX2,d2cdr2)-SBUFFA(IDX2,Br)
        END_DO

        ! Free up the dAdr space -- get its two angular derivatives
        Call d_by_dtheta(wsp%s2a,dadr,ftemp1)
        Call d_by_dphi(  wsp%s2a,dadr,ftemp2)

        !////////// [Del x B]_phi //////////////////////////
        ! overwrite d_a_dr with d_d_phi(d_a_dr)
        DO_IDX2
            SBUFFA(IDX2,dadr) = ftemp2(mp)%data(IDX2)
        END_DO
        ! overwrite ftemp2 with d_d_theta (d2cdr2-br)
        Call d_by_dtheta(  wsp%s2a,d2cdr2,ftemp2)

        ! Add this term to d_d_phi(d_a_dr) to build rsintheta [del x b]_phi (overwrite dadr)
        DO_IDX2
            SBUFFA(IDX2,curlbphi) = SBUFFA(IDX2,curlbphi)+ftemp2(mp)%data(IDX2)
            SBUFFA(IDX2,curlbphi) = SBUFFA(IDX2,curlbphi)
        END_DO

        !/////////////[Del x B]_theta ///////////////////////
        Call d_by_dphi(  wsp%s2a,d2cdr2,ftemp2)       ! get phi derivative of d2cdr2-Br

        ! Combine with ftemp1 to build rsintheta [del x B]_theta (overwrites d2cdr2)
        DO_IDX2
            SBUFFA(IDX2,curlbtheta) = (ftemp1(mp)%data(IDX2)-ftemp2(mp)%data(IDX2))
        END_DO


        !////////////B Theta
        ! Free up the A space -- get its two angular derivatives
        Call d_by_dtheta(wsp%s2a,avar,ftemp1)
        Call d_by_dphi(  wsp%s2a,avar,ftemp2)


        ! overwrite A with dA_d_phi
        DO_IDX2
            SBUFFA(IDX2,Avar) = ftemp2(mp)%data(IDX2)
        END_DO

        ! overwrite ftemp2 with d_d_theta (dcdr)
        Call d_by_dtheta(  wsp%s2a,dcdr,ftemp2)

        ! Add this term to dA_d_phi to build rsintheta B_theta
        DO_IDX2
            SBUFFA(IDX2,Avar) = SBUFFA(IDX2,Avar)+ftemp2(mp)%data(IDX2)
        END_DO

        !///////////// Bphi
        Call d_by_dphi(  wsp%s2a,dcdr,ftemp2)       ! get phi derivative of dcdr

        ! Combine with ftemp1 to build rsintheta B_phi
        DO_IDX2
            SBUFFA(IDX2,dcdr) = ftemp2(mp)%data(IDX2)-ftemp1(mp)%data(IDX2)
        END_DO
    End Subroutine Compute_BandCurlB

    Subroutine Bfield_Derivatives()
        Implicit None
        Integer :: r, l, m, mp, imi

        ! These terms are only needed if we want to output
        ! inductions terms in the diagnostics


        !/////////////////////////////////
        ! sintheta dB theta dr
        Call d_by_dtheta(ftemp3,ftemp1)
        Call d_by_dphi(ftemp4,    ftemp2)


        DO_IDX2
            ftemp1(mp)%data(IDX2) = ftemp1(mp)%data(IDX2)+ftemp2(mp)%data(IDX2)
        END_DO

        DO_IDX2
            ASBUFFA(IDX2,dbtdr_cb) = ftemp1(mp)%data(IDX2)*one_over_r(r)
        END_DO

        DO_IDX2
            ASBUFFA(IDX2,dbtdr_cb) = ASBUFFA(IDX2,dbtdr_cb)- &
                & SBUFFA(IDX2,btheta)*OneOverRSquared(r)  ! (take care) btheta is really rsintheta btheta
        END_DO                                              ! hence 1/r^2 instead of 1/r

        !/////////////////////////////////
        ! sintheta dB phi dr
        Call d_by_dphi(ftemp3,ftemp1)
        Call d_by_dtheta(ftemp4,    ftemp2)

        DO_IDX2
            ftemp1(mp)%data(IDX2) = ftemp1(mp)%data(IDX2)-ftemp2(mp)%data(IDX2)
        END_DO

        DO_IDX2
            ASBUFFA(IDX2,dbpdr_cb) = ftemp1(mp)%data(IDX2)*one_over_r(r)
        END_DO

        DO_IDX2
            ASBUFFA(IDX2,dbpdr_cb) = ASBUFFA(IDX2,dbpdr_cb)- &
                &  SBUFFA(IDX2,bphi)*OneOverRSquared(r) ! (take care) bphi is really rsinthetabphi
        END_DO

        !/////////////////////////////////////////
        ! dB r dr  (dbrdr_cb holds dcdr up until this point)

        DO_IDX2
            ASBUFFA(IDX2,dbrdr_cb) = l_l_plus1(m:l_max)* &
                & ASBUFFA(IDX2,dbrdr_cb)*OneOverRSquared(r)
        END_DO


        DO_IDX2
            ASBUFFA(IDX2,dbrdr_cb) = ASBUFFA(IDX2,dbrdr_cb)- &
                & SBUFFA(IDX2,br)*Two_Over_R(r)
        END_DO

        ! sintheta dbrdt
        Call d_by_dtheta(wsp%s2a,br,ftemp1)
        DO_IDX2
            ASBUFFA(IDX2,dbrdt_cb) = ftemp1(mp)%data(IDX2)
        END_DO


    End Subroutine BField_Derivatives

    Subroutine Hybrid_Output_Initial()
        Implicit None
        Integer :: r, l, m, mp, imi
        If (magnetism) Then
            Call Allocate_rlm_Field(ftemp3)
            Call Allocate_rlm_Field(ftemp4)
            ! First we grab a copy of several variables whose
            ! values will be overwritten in B and J are computed
        
            ! Convert A to ell(ell+1) A/r^2  (i.e. [curl B]_r)
            DO_IDX2
                ASBUFFA(IDX2,avar_cb) = l_l_plus1(m:l_max)* &
                                        SBUFFA(IDX2,avar)*OneOverRSquared(r)
            END_DO

            DO_IDX2
                ASBUFFA(IDX2,dbrdr_cb) = SBUFFA(IDX2,dcdr)
            END_DO

            DO_IDX2
                ftemp3(mp)%data(IDX2)  = SBUFFA(IDX2,d2cdr2)
            END_DO

            DO_IDX2
                ftemp4(mp)%data(IDX2) = SBUFFA(IDX2,dadr)
            END_DO

        Endif
    End Subroutine Hybrid_Output_Initial

    Subroutine Hybrid_Output_Final()
        Implicit None
        Integer :: r, l, m, mp, imi
        Do mp = my_mp%min, my_mp%max
            ASBUFFA(l_max,:,:,:) = 0.0d0
        Enddo
        Call Hydro_Output_Derivatives()
        If (magnetism) Then
            ! We compute some derivatives of B as well
            Call BField_Derivatives()
            Call Deallocate_rlm_Field(ftemp3)
            Call Deallocate_rlm_Field(ftemp4)
        Endif
        Call cobuffer%construct('p2a')
        cobuffer%config = 'p2a'
        Call Legendre_Transform(cobuffer%s2a,cobuffer%p2a)
        Call cobuffer%deconstruct('s2a')

        Call cobuffer%reform()

    End Subroutine Hybrid_Output_Final

    Subroutine Adjust_Emf()
        Implicit None
        Integer :: m, mp, r,imi

        Call d_by_sdtheta(wsp%s2b, emfphi,ftemp1)
        Call d_by_dphi(wsp%s2b,emftheta,ftemp2)

        Call Allocate_rlm_Field(ftemp3)
        ! Copy out emf_theta before we overwrite it
        DO_IDX2
            ftemp3(mp)%data(IDX2) = SBUFFB(IDX2,emftheta)
        END_DO

        ! Now for the C RHS, formed from the radial component of the curl of the emf
        ! cvar overwrites emftheta
        DO_IDX2
            SBUFFB(IDX2,Cvar) = ( ftemp1(mp)%data(IDX2)- &
                & ftemp2(mp)%data(IDX2) )*over_l_l_plus1(m:l_max)
        END_DO


        Call d_by_dphi(wsp%s2b,emfphi,ftemp2)
        ! Move ftemp3 (emftheta) into emfphi's old spot
        DO_IDX2
            SBUFFB(IDX2,emfphi)=ftemp3(mp)%data(IDX2)
        END_DO
        Call d_by_sdtheta(wsp%s2b, emfphi,ftemp1)

        DO_IDX2
            SBUFFB(IDX2,emfphi) = ( ftemp2(mp)%data(IDX2)+ &
                & ftemp1(mp)%data(IDX2) )*over_l_l_plus1(m:l_max)
        END_DO
        Call DeAllocate_rlm_Field(ftemp3)
        ! Ensure there is no ell=0 emf  -- should I do this?
        !rmn1 = (emfr-1)    *tnr+1
        !rmn2 = (emftheta-1)*tnr+1
        !rmn3 = (emfphi-1)  *tnr+1
        !Do mp = my_mp%min, my_mp%max
        !    m = m_values(mp)
        !    if (m .eq. 0) then
        !        wsp%s2b(mp)%data(0,rmn1:rmn1+tnr-1) = 0.0d0
        !        wsp%s2b(mp)%data(0,rmn2:rmn2+tnr-1) = 0.0d0
        !        wsp%s2b(mp)%data(0,rmn3:rmn3+tnr-1) = 0.0d0
        !    endif
        !Enddo
    End Subroutine Adjust_EMF

    Subroutine Adjust_TimeStep()
        Implicit None
        Real*8 :: maxt2, maxt
        Character*8 :: dtfmt ='(ES10.4)'
        Character*14 :: tmstr, tmstr2

        Call wsp%unload_cargo(global_msgs)


        maxt2 = global_msgs(1)
        if (maxt2 .gt. 0.0d0) Then
            maxt = 1.0d0/sqrt(maxt2)

            if (deltat .lt. maxt*cflmin) then
                ! we can increase our timestep
                new_deltat = Min(cflmax*maxt,max_time_step)

            elseif (deltat .gt. (maxt*cflmax)) then
                new_deltat = cflmax*maxt
                if (new_deltat .gt. deltat*(1.0d0-min_dt_change)) then
                    ! As much as possible, we would like to avoid
                    ! changing the timestep (slow process).  When we do change it,
                    ! make sure we give it a good bump.
                    new_deltat = deltat*(1.0d0-min_dt_change)
                endif
            endif
        Endif
        if (new_deltat .gt. (max_time_step*1.000001d0)) Then
            new_deltat = max_time_step
        Endif
        If (new_deltat .ne. deltat) Then
            new_timestep = .true.
        Endif
        If (new_deltat .lt. min_time_step) Then
            If (my_rank .eq. 0) Then
                Call stdout%print('Time step became too small.')
                Write(tmstr,dtfmt)new_deltat
                Write(tmstr2,dtfmt)min_time_step
                Call stdout%print(' DeltaT became : '//tmstr//'  Min DeltaT Allowed:   '//tmstr2)
                Call stdout%partial_flush()
            Endif
            Call pfi%exit()
            Stop
        Endif

    End Subroutine Adjust_TimeStep

    Subroutine Allocate_rlm_Field(arr)
        Implicit None
        Type(rmcontainer3D), Intent(InOut), Allocatable :: arr(:)
        Integer :: mp,m


        Allocate(arr(my_mp%min:my_mp%max))
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            Allocate(arr(mp)%data(m:l_max,my_r%min:my_r%max,1:2))
            arr(mp)%data(:,:,:) = 0.0d0
        Enddo
    End Subroutine Allocate_rlm_Field

    Subroutine DeAllocate_rlm_Field(arr)
        Implicit None
        Type(rmcontainer3D), Intent(InOut), Allocatable :: arr(:)
        Integer :: mp
        Do mp = my_mp%min, my_mp%max
            DeAllocate(arr(mp)%data)
        Enddo
        DeAllocate(arr)
    End Subroutine DeAllocate_rlm_Field
End Module Sphere_Hybrid_Space
