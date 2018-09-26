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
!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_SECOND_DERIVATIVES
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Second_Derivatives
    Use Diagnostics_Base
    Use Structures
    Use Spectral_Derivatives
    Implicit None


    Integer, Allocatable :: ddindmap(:,:)
    Integer :: nddfields
    Logical :: compute_vr_dd   = .false.
    Logical :: compute_vt_dd   = .false.
    Logical :: compute_vp_dd   = .false.

    Logical :: compute_pvar_dd = .false.
    Logical :: compute_tvar_dd = .false.

    Logical :: compute_br_dd   = .false.
    Logical :: compute_bt_dd   = .false.
    Logical :: compute_bp_dd   = .false.
Contains

    Subroutine Init_Derivative_Logic()
        IMPLICIT NONE
        Integer :: i
        !///////////////////////////////////////////////////////
        ! Check to see if the user has specified any of the second
        ! derivatives individually
        do i = dv_r_d2r, dvm_r_d2tp,3
            if (sometimes_compute(i)) compute_vr_dd = .true.
        enddo

        do i = dv_theta_d2r, dvm_theta_d2tp,3
            if (sometimes_compute(i)) compute_vt_dd = .true.
        enddo
        do i = dv_phi_d2r, dvm_phi_d2tp,3
            if (sometimes_compute(i)) compute_vp_dd = .true.
        enddo

        do i = db_r_d2r, dbm_r_d2tp,3
            if (sometimes_compute(i)) compute_br_dd = .true.
        enddo

        do i = db_theta_d2r, dbm_theta_d2tp,3
            if (sometimes_compute(i)) compute_bt_dd = .true.
        enddo

        do i = db_phi_d2r, dbm_phi_d2tp,3
            if (sometimes_compute(i)) compute_bp_dd = .true.
        enddo

        do i = entropy_d2r, entropy_m_d2tp,2
            if (sometimes_compute(i)) compute_tvar_dd = .true.
        enddo

        do i = pressure_d2r, pressure_m_d2tp,2
            if (sometimes_compute(i)) compute_pvar_dd = .true.
        enddo


        !//////////////////////////////////////////////////////////////////
        ! Terms related to viscosity
        If (sometimes_compute(viscous_force_r))  compute_vr_dd = .true.
        If (sometimes_compute(viscous_pforce_r)) compute_vr_dd = .true.
        If (sometimes_compute(viscous_mforce_r)) compute_vr_dd = .true.

        If (sometimes_compute(viscous_force_theta))  compute_vt_dd = .true.
        If (sometimes_compute(viscous_pforce_theta)) compute_vt_dd = .true.
        If (sometimes_compute(viscous_mforce_theta)) compute_vt_dd = .true.

        If (sometimes_compute(viscous_force_phi))  compute_vp_dd = .true.
        If (sometimes_compute(viscous_pforce_phi)) compute_vp_dd = .true.
        If (sometimes_compute(viscous_mforce_phi)) compute_vp_dd = .true.

        If (sometimes_compute(visc_work) .or. sometimes_compute(visc_work_pp) &
            .or. sometimes_compute(visc_work_mm) ) Then
            compute_vr_dd = .true.
            compute_vt_dd = .true.
            compute_vp_dd = .true.
        Endif

        Do i = visc_flux_r, visc_fluxmm_r, 3
            if (sometimes_compute(i)) compute_vr_dd = .true.
        Enddo
        Do i = visc_flux_theta, visc_fluxmm_theta, 3
            if (sometimes_compute(i)) compute_vt_dd = .true.
        Enddo
        Do i = visc_flux_phi, visc_fluxmm_phi, 3
            if (sometimes_compute(i)) compute_vp_dd = .true.
        Enddo


        !/////////////////////////////////////////////////////////////////
        ! Check to see if we are computing thermal diffusion terms
        If (sometimes_compute(s_diff) .or. sometimes_compute(sp_diff) &
            .or. sometimes_compute(sm_diff) ) Then
            compute_tvar_dd = .true.
            compute_pvar_dd = .true.
        Endif

        do i = s_diff_r, sm_diff_phi
            if (sometimes_compute(i)) compute_tvar_dd = .true.
        enddo

        !//////////////////////////////////////////////////////
        ! Are we computing magnetic diffusion terms?
        If (sometimes_compute(induct_diff_r) .or. sometimes_compute(induct_diff_bm_r) &
            .or. sometimes_compute(induct_diff_bp_r) ) Then
            compute_br_dd=.true.
        Endif

        If (sometimes_compute(induct_diff_theta) .or. sometimes_compute(induct_diff_bm_theta) &
            .or. sometimes_compute(induct_diff_bp_theta) ) Then
            compute_bt_dd=.true.
        Endif

        If (sometimes_compute(induct_diff_phi) .or. sometimes_compute(induct_diff_bm_phi) &
            .or. sometimes_compute(induct_diff_bp_phi) ) Then
            compute_bp_dd=.true.
        Endif


        If (sometimes_compute(idiff_work) .or. sometimes_compute(idiff_work_pp) &
            .or. sometimes_compute(idiff_work_mm) ) Then
            compute_br_dd = .true.
            compute_bt_dd = .true.
            compute_bp_dd = .true.
        Endif

        ! Execute a lot of compute_q logic here to see if the
        ! different compute_xx_dd variables should be set to true.

        If (compute_vr_dd) need_second_derivatives = .true.
        If (compute_vt_dd) need_second_derivatives = .true.
        If (compute_vp_dd) need_second_derivatives = .true.

        If (compute_tvar_dd) need_second_derivatives = .true.
        If (compute_pvar_dd) need_second_derivatives = .true.

        If (compute_br_dd) need_second_derivatives = .true.
        If (compute_bt_dd) need_second_derivatives = .true.
        If (compute_bp_dd) need_second_derivatives = .true.

        ! Turbulent KE generation
        If (sometimes_compute(production_buoyant_pKE)) need_second_derivatives = .true.
        If (sometimes_compute(production_shear_pKE)) need_second_derivatives = .true.

        If (sometimes_compute(dissipation_viscous_pKE)) need_second_derivatives = .true.

        If (sometimes_compute(transport_pressure_pKE)) need_second_derivatives = .true.
        If (sometimes_compute(transport_viscous_pKE)) Then
            need_second_derivatives = .true.
            compute_vr_dd = .true.
            compute_vt_dd = .true.
            compute_vp_dd = .true.
        Endif
        If (sometimes_compute(transport_turbadvect_pKE)) need_second_derivatives = .true.
        If (sometimes_compute(transport_meanadvect_pKE)) need_second_derivatives = .true.

        If (sometimes_compute(rflux_pressure_pKE)) need_second_derivatives = .true.
        If (sometimes_compute(rflux_viscous_pKE)) need_second_derivatives = .true.
        If (sometimes_compute(rflux_turbadvect_pKE)) need_second_derivatives = .true.
        If (sometimes_compute(rflux_meanadvect_pKE)) need_second_derivatives = .true.

        If (sometimes_compute(thetaflux_pressure_pKE)) need_second_derivatives = .true.
        If (sometimes_compute(thetaflux_viscous_pKE)) need_second_derivatives = .true.
        If (sometimes_compute(thetaflux_turbadvect_pKE)) need_second_derivatives = .true.
        If (sometimes_compute(thetaflux_meanadvect_pKE)) need_second_derivatives = .true.

    End Subroutine Init_Derivative_Logic

    Subroutine Initialize_Second_Derivatives()
        ! Initializes all the indexing related to computing and
        ! accessing second derivatives at output time.
        ! Most of the actual indexing is handled by Set_DD_Indices
        IMPLICIT NONE
        INTEGER :: ndind
        INTEGER :: ddfcount(3,2)



        Call Init_Derivative_Logic()


        nddfields = 0   ! Number of fields whose second derivatives we want
        ndind = 0       ! internal indexing variable

        IF (compute_vr_dd) nddfields = nddfields +1
        IF (compute_vt_dd) nddfields = nddfields +1
        IF (compute_vp_dd) nddfields = nddfields +1

        IF (compute_tvar_dd) nddfields = nddfields +1
        IF (compute_pvar_dd) nddfields = nddfields +1

        IF (magnetism) THEN
            IF (compute_br_dd) nddfields = nddfields +1
            IF (compute_bt_dd) nddfields = nddfields +1
            IF (compute_bp_dd) nddfields = nddfields +1
        ENDIF




        Allocate(ddindmap(4,nddfields*2))
        ddindmap(:,:) = -1


        If (compute_vr_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dvrdrdr, dvrdrdt, dvrdrdp, &
                                      dvrdtdt, dvrdtdp, dvrdpdp, &
                                      dvrdr  , dvrdt  , dvrdp)
        Endif

        IF (compute_vt_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dvtdrdr, dvtdrdt, dvtdrdp, &
                                      dvtdtdt, dvtdtdp, dvtdpdp, &
                                      dvtdr  , dvtdt  , dvtdp)
        ENDIF

        IF (compute_vp_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dvpdrdr, dvpdrdt, dvpdrdp, &
                                      dvpdtdt, dvpdtdp, dvpdpdp, &
                                      dvpdr  , dvpdt  , dvpdp)
        ENDIF

        IF (compute_tvar_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dtdrdr, dtdrdt, dtdrdp, &
                                      dtdtdt, dtdtdp, dtdpdp, &
                                      dtdr  , dtdt  , dtdp)
        ENDIF

        IF (compute_pvar_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, ddindmap, &
                                      dpdrdr, dpdrdt, dpdrdp, &
                                      dpdtdt, dpdtdp, dpdpdp, &
                                      dpdr  , dpdt  , dpdp)
        ENDIF

        If (magnetism) THEN
            If (compute_br_dd) THEN
                ndind = ndind+1
                Call set_dd_indices(ndind,nddfields, ddindmap, &
                                          dbrdrdr, dbrdrdt, dbrdrdp, &
                                          dbrdtdt, dbrdtdp, dbrdpdp, &
                                          dbrdr  , dbrdt  , dbrdp)
            Endif

            IF (compute_bt_dd) THEN
                ndind = ndind+1
                Call set_dd_indices(ndind,nddfields, ddindmap, &
                                          dbtdrdr, dbtdrdt, dbtdrdp, &
                                          dbtdtdt, dbtdtdp, dbtdpdp, &
                                          dbtdr  , dbtdt  , dbtdp)
            ENDIF

            IF (compute_bp_dd) THEN
                ndind = ndind+1
                Call set_dd_indices(ndind,nddfields, ddindmap,      &
                                          dbpdrdr, dbpdrdt, dbpdrdp, &
                                          dbpdtdt, dbpdtdp, dbpdpdp, &
                                          dbpdr  , dbpdt  , dbpdp)
            ENDIF
        ENDIF



        ddfcount(1,1) = nddfields*4 ! config 1a
        ddfcount(2,1) = nddfields*4 ! 2a
        ddfcount(3,1) = nddfields*7 ! 3a
        ddfcount(3,2) = nddfields*2 ! 3b
        ddfcount(2,2) = nddfields*2 ! 2b
        ddfcount(1,2) = nddfields*4 ! 1b


        Call d2buffer%init(field_count = ddfcount, config = 'p3b')

        Call d2buffer%construct('p3a')
        !WRITE(6,*)'BCHECK: ', shape(d2buffer%p3a)
        Call d2buffer%deconstruct('p3a')
    End Subroutine Initialize_Second_Derivatives

    Subroutine Compute_Second_Derivatives(inbuffer)
        Implicit None
        INTEGER :: i,j, imi, mp, m
        INTEGER :: r,k, t
        Real*8, Intent(InOut) :: inbuffer(1:,my_r%min:,my_theta%min:,1:)
        Type(rmcontainer3D), Allocatable :: ddtemp(:)






        ! Here were compute all second derivatives for N variables
        ! Outline of the process:
        ! 1)  Intialize p3b space of d2buffer
        ! 2)  Load slots [1 : N] of d2buffer with d_by_dr for each variable
        ! 3)  Load slots [N+1 : 2N] with d_by_dtheta for each variable
        ! 4)  Move to p1b/p1a configuration
        ! 5)  Load slots [2N+1 : 3N] with d2_by_dr2 for each variable
        ! 6)  Load slots [3N+1 : 4N] with d2_by_drdt for each variable; move to s2a
        ! 7)  Move to s2a and load sintheta*dxdtdt into slots [N+1 : 2N] (ovewriting dxdt)
        ! 8)  Move to p3a and load d2_by_drdphi into slots [4N+1: 5N] for each variable
        ! 9)  Load d2_by_dphi2 into slots [5N+1: 6N] for each variable
        ! 10) Load d2_by_dtdphi into slots [6N+1 : 7N] for each variable
        ! 11) FFT, correct for sintheta factor, compute means and fluctuations

        ! When this routine is complete, the contents of d2buffer%p3a will be
        ! [   1 : N  ] -- workspace;
        ! [ N+1 : 2N ] -- dxdtdt
        ! [2N+1 : 3N ] -- dxdrdr
        ! [3N+1 : 4N ] -- dxdrdt
        ! [4N+1 : 5N ] -- dxdrdp
        ! [5N+1 : 6N ] -- dxdpdp
        ! [6N+1 : 7N ] -- dxdtdp

        !///////////////////////////////////////////////////////////
        ! Step 1: Initialize p3b portion of d2buffer
        !Write(6,*)'Constructing the buffer...', d2buffer%nf3b

        Call d2buffer%construct('p3b')
        d2buffer%config = 'p3b'

        !Write(6,*)'Complete...'

        !///////////////////////////////////////////////////////////
        ! Steps 2-3:  Load radial and theta derivatives
        Do i = 1, nddfields*2
            d2buffer%p3b(:,:,:,ddindmap(1,i)) = inbuffer(:,:,:,ddindmap(2,i))
        Enddo

        !Write(6,*)'Steps 1-3 complete.'

        !////////////////////////////////////////////////////////////////
        ! Step 4:  Move to p1b/p1a configuration
        Call fft_to_spectral(d2buffer%p3b, rsc = .true.)
        Call d2buffer%reform()
        Call d2buffer%construct('s2b')
        Call Legendre_Transform(d2buffer%p2b,d2buffer%s2b)
        Call d2buffer%deconstruct('p2b')
        d2buffer%config ='s2b'



        !Move to p1b configuration
        Call d2buffer%reform()
        Call d2buffer%construct('p1a')
        Call gridcp%To_Spectral(d2buffer%p1b,d2buffer%p1a)
        Call gridcp%dealias_buffer(d2buffer%p1a)
        d2buffer%p1b = 0.0
        d2buffer%config='p1a'

        !Write(6,*)'Steps 4 complete.'

        !////////////////////////////////////////////////////////////
        ! Steps 5-6:  Load d2_by_dr2 and d2_by_drdt into the buffer
        Do i = 1, nddfields*2
            j = i+nddfields*2
            Call gridcp%d_by_dr_cp(i,j,d2buffer%p1a,1)
        Enddo
        Call gridcp%From_Spectral(d2buffer%p1a,d2buffer%p1b)
        d2buffer%p1a=d2buffer%p1b
        Call d2buffer%deconstruct('p1b')

        ! Ordering of fields in buffer is now dxdr, dxdt, dxdrdr, dxdrdt
        !Write(6,*)'Steps 5-6 complete.'

        !///////////////////////////////////////////////////////////////
        ! Step 7:  Move to s2a & overwrite dxdt with sintheta*{dxdtdt}
        Call d2buffer%reform()

        !Write(6,*)'Step 7 reformation complete'

        Call Allocate_rlm_Field(ddtemp)

        !Write(6,*)'Step 7 allocation complete', nddfields

        Do i = nddfields+1,nddfields*2
            ! We overwrite dxdt with sintheta* {dxdtdt}
            Call d_by_dtheta(d2buffer%s2a,i,ddtemp)
            DO_IDX2
                d2buffer%s2a(mp)%data(IDX2,i) = ddtemp(mp)%data(IDX2)
            END_DO
        Enddo

        !Write(6,*)'Step 7 derivatves complete'

        Call DeAllocate_rlm_Field(ddtemp)

        !Write(6,*)'Step 7 deallocation complete'

        Call d2buffer%construct('p2a')
        Call Legendre_Transform(d2buffer%s2a,d2buffer%p2a)
        Call d2buffer%deconstruct('s2a')
        d2buffer%config = 'p2a'

        !Write(6,*)'Step 7 complete.'

        !/////////////////////////////////////////////////////////////////////
        !  Steps 8-10 : phi derivatives
        Call d2buffer%reform() ! move to p3a

        ! Ordering of fields in buffer is now dxdr, sintheta*{dxdtdt}, dxdrdr, dxdrdt

        !Compute dxdrdp
        Do i = 1, nddfields
            j = i+nddfields*4
            Call d_by_dphi(d2buffer%p3a,i,j)
        Enddo

        ! Grab dxdp and dxdt from inbuffer
        ! Inbuffer is in physical space, so after we copy into p3b, we need to FFT
        Call d2buffer%construct('p3b')
        Do i = 1, nddfields
            d2buffer%p3b(:,:,:,i) = inbuffer(:,:,:,ddindmap(3,i))
            d2buffer%p3b(:,:,:,i+nddfields) = inbuffer(:,:,:,ddindmap(4,i))
        Enddo
        Call fft_to_spectral_rsc(d2buffer%p3b)  ! call this version since dropping in midway through loop

        !Copy dxdp into dxdr space, then calculate dxdpdp
        Do i = 1, nddfields
            j = i+nddfields*5
            d2buffer%p3a(:,:,:,i) = d2buffer%p3b(:,:,:,i)

            Call d_by_dphi(d2buffer%p3a,i,j)
        Enddo

        !Copy dxdt into dxdr space, then calculate dxdtdp
        Do i = 1, nddfields
            j = i+nddfields*6
            d2buffer%p3a(:,:,:,i) = d2buffer%p3b(:,:,:,i+nddfields)
            Call d_by_dphi(d2buffer%p3a,i,j)
        Enddo

        Call d2buffer%deconstruct('p3b')
       ! Write(6,*)'Steps 8-10 complete.'

        !//////////////////////////////////////////
        ! Step 11:   Finalize
        ! FFT
        Call fft_to_physical(d2buffer%p3a,rsc = .true.)

        ! Convert sintheta*{dxdtdt} to dxdtdt
        Do i = nddfields+1,nddfields*2
            DO_PSI
                d2buffer%p3a(PSI,i) = d2buffer%p3a(PSI,i)*csctheta(t)
            END_DO
        Enddo
        !D2buffer is now initialized.  Ordering of fields is:
        ! [dxdt], dxdtdt,dxdrdr,dxdrdt, dxdrdp, dxdpdp, dxdtdp
        ! Note that we have one redundant derivative, dxdt, stored for each field

        ! Now compute the means and fluctuations
        Allocate(d2_ell0(my_r%min:my_r%max,1:nddfields*7))
        Allocate(d2_m0(my_r%min:my_r%max,my_theta%min:my_theta%max,1:nddfields*7))
        Allocate(d2_fbuffer(1:n_phi,my_r%min:my_r%max, &
            my_theta%min:my_theta%max,1:nddfields*7))

        Call ComputeEll0(d2buffer%p3a,d2_ell0)
        Call   ComputeM0(d2buffer%p3a,d2_m0)

        DO j = 1,nddfields*7
            DO_PSI
                d2_fbuffer(PSI,j) = d2buffer%p3a(PSI,j) - d2_m0(PSI2,j)
            END_DO
        ENDDO

        !Write(6,*)'Step 11 complete.'
    End Subroutine Compute_Second_Derivatives

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

    Subroutine Set_DD_Indices(iind, nskip, iindmap, &
                                     dxdrdr, dxdrdt, dxdrdp, &
                                     dxdtdt, dxdtdp, dxdpdp, &
                                       dxdr,   dxdt, dxdp)
        ! Sets indexing within indmap and assigned values to
        ! dxdidj consistent with the logic used in
        ! Compute_Second_Derivatives()
        ! [ N+1 : 2N ] -- dxdtdt
        ! [2N+1 : 3N ] -- dxdrdr
        ! [3N+1 : 4N ] -- dxdrdt
        ! [4N+1 : 5N ] -- dxdrdp
        ! [5N+1 : 6N ] -- dxdpdp
        ! [6N+1 : 7N ] -- dxdtdp
        Implicit None
        INTEGER, Intent(In)    :: iind, nskip
        INTEGER, INTENT(OUT)   :: dxdtdt, dxdrdr, dxdrdt
        INTEGER, INTENT(OUT)   :: dxdrdp, dxdpdp, dxdtdp
        INTEGER, INTENT(IN)    :: dxdr, dxdt, dxdp
        INTEGER, INTENT(INOUT) :: iindmap(:,:)
        dxdtdt = iind+nskip
        dxdrdr = iind+nskip*2
        dxdrdt = iind+nskip*3
        dxdrdp = iind+nskip*4
        dxdpdp = iind+nskip*5
        dxdtdp = iind+nskip*6

        iindmap(1,iind          ) = iind
        iindmap(2,iind          ) = dxdr
        iindmap(3,iind          ) = dxdp
        iindmap(4,iind          ) = dxdt
        !iindmap(1,iind          ) = iind+nskip
        iindmap(1,iind+nskip) = iind+nskip
        iindmap(2,iind+nskip) = dxdt


    End Subroutine Set_DD_Indices
End Module Diagnostics_Second_Derivatives
