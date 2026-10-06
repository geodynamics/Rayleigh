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
    Use Finite_Difference, Only : d_by_dx3d3
    Implicit None


    ! Each field is either a scalar (v_r, T, P, B_r) or one horizontal component of a vector
    ! (v_theta, v_phi, B_theta, B_phi).  See Compute_Second_Derivatives.
    Integer, Parameter :: dd_scalar = 1, dd_theta = 2, dd_phi = 3
    Integer, Allocatable :: dd_type(:)      ! dd_scalar, dd_theta or dd_phi
    Integer, Allocatable :: dd_src(:,:)     ! (1:8,i): buffer indices of x, dxdr, dxdt, dxdp and, for horizontal
                                            !   components, of the other component y, dydr, dydt, dydp
    Integer, Allocatable :: dd_slot(:,:)    ! (1:2,i): d2buffer p3b slots of the fields transformed for field i
    Integer, Allocatable :: dd_scr(:)       ! horizontal components: index of the p3a scratch pair holding
                                            !   d2y/drdphi and d2y/dthetadphi
    Real*8, Allocatable  :: dd_sgn(:)       ! +1 for theta components (Q = horizontal divergence), 
                                            ! -1 for phi components (Q = radial vorticity)
    Integer :: nddfields, ndd_trans, ndd_horiz

Contains

    Subroutine Second_Derivative_Logic(check, l_compute_vr_dd, l_compute_vt_dd, l_compute_vp_dd, &
                                 l_compute_tvar_dd, l_compute_pvar_dd, &
                                 l_compute_br_dd, l_compute_bt_dd, l_compute_bp_dd, need_dd)
        ! Trigger-code logic shared between the once-at-startup buffer-sizing
        ! pass (check => Sometimes_Compute, decides which fields ever need
        ! second derivatives, for dd_src/buffer sizing) and the per-iteration
        ! recheck of whether Compute_Second_Derivatives needs to run this
        ! iteration (check => Compute_Quantity).
        IMPLICIT NONE
        Procedure(Quantity_Check_If) :: check
        Logical, Intent(Out) :: l_compute_vr_dd, l_compute_vt_dd, l_compute_vp_dd
        Logical, Intent(Out) :: l_compute_tvar_dd, l_compute_pvar_dd
        Logical, Intent(Out) :: l_compute_br_dd, l_compute_bt_dd, l_compute_bp_dd
        Logical, Intent(Out) :: need_dd
        Integer :: i

        l_compute_vr_dd   = .false.
        l_compute_vt_dd   = .false.
        l_compute_vp_dd   = .false.
        l_compute_tvar_dd = .false.
        l_compute_pvar_dd = .false.
        l_compute_br_dd   = .false.
        l_compute_bt_dd   = .false.
        l_compute_bp_dd   = .false.
        need_dd      = .false.

        !///////////////////////////////////////////////////////
        ! Check to see if the user has specified any of the second
        ! derivatives individually
        do i = dv_r_d2r, dvm_r_d2tp,3
            if (check(i)) l_compute_vr_dd = .true.
        enddo

        do i = dv_theta_d2r, dvm_theta_d2tp,3
            if (check(i)) l_compute_vt_dd = .true.
        enddo
        do i = dv_phi_d2r, dvm_phi_d2tp,3
            if (check(i)) l_compute_vp_dd = .true.
        enddo

        do i = db_r_d2r, dbm_r_d2tp,3
            if (check(i)) l_compute_br_dd = .true.
        enddo

        do i = db_theta_d2r, dbm_theta_d2tp,3
            if (check(i)) l_compute_bt_dd = .true.
        enddo

        do i = db_phi_d2r, dbm_phi_d2tp,3
            if (check(i)) l_compute_bp_dd = .true.
        enddo

        do i = entropy_d2r, entropy_m_d2tp,2
            if (check(i)) l_compute_tvar_dd = .true.
        enddo

        do i = pressure_d2r, pressure_m_d2tp,2
            if (check(i)) l_compute_pvar_dd = .true.
        enddo


        !//////////////////////////////////////////////////////////////////
        ! Terms related to viscosity
        If (check(visc_work) .or. &
            check(viscous_force_r) .or. &
            check(curl_viscous_force_theta) .or. &
            check(curl_viscous_force_theta_squared) .or. &
            check(curl_viscous_force_phi) .or. &
            check(curl_viscous_force_phi_squared) .or. &
            check(viscous_pforce_r) .or. &
            check(curl_viscous_pforce_theta) .or. &
            check(curl_viscous_pforce_phi) .or. &
            check(viscous_mforce_r) .or. &
            check(curl_viscous_mforce_theta) .or. &
            check(curl_viscous_mforce_phi)) Then
            l_compute_vr_dd = .true.
        Endif


        If (check(visc_work_pp) .or. &
            check(viscous_force_theta) .or. &
            check(curl_viscous_force_r) .or. &
            check(curl_viscous_force_r_squared) .or. &
            check(curl_viscous_force_phi) .or. &
            check(curl_viscous_force_phi_squared) .or. &
            check(viscous_pforce_theta) .or. &
            check(curl_viscous_pforce_r) .or. &
            check(curl_viscous_pforce_phi) .or. &
            check(viscous_mforce_theta) .or. &
            check(curl_viscous_mforce_r) .or. &
            check(curl_viscous_mforce_phi)) Then
            l_compute_vt_dd = .true.
        Endif

        If (check(visc_work_mm) .or. &
            check(viscous_force_phi) .or. &
            check(curl_viscous_force_r) .or. &
            check(curl_viscous_force_r_squared) .or. &
            check(curl_viscous_force_theta) .or. &
            check(curl_viscous_force_theta_squared) .or. &
            check(viscous_pforce_phi) .or. &
            check(curl_viscous_pforce_r) .or. &
            check(curl_viscous_pforce_theta) .or. &
            check(viscous_mforce_phi) .or. &
            check(curl_viscous_mforce_r) .or. &
            check(curl_viscous_mforce_theta)) Then
            l_compute_vp_dd = .true.
        Endif

        Do i = visc_flux_r, visc_fluxmm_r, 3
            if (check(i)) l_compute_vr_dd = .true.
        Enddo
        Do i = visc_flux_theta, visc_fluxmm_theta, 3
            if (check(i)) l_compute_vt_dd = .true.
        Enddo
        Do i = visc_flux_phi, visc_fluxmm_phi, 3
            if (check(i)) l_compute_vp_dd = .true.
        Enddo


        !/////////////////////////////////////////////////////////////////
        ! Check to see if we are computing thermal diffusion terms
        If (check(s_diff) .or. check(sp_diff) &
            .or. check(sm_diff) ) Then
            l_compute_tvar_dd = .true.
            l_compute_pvar_dd = .true.
        Endif

        do i = s_diff_r, sm_diff_phi
            if (check(i)) l_compute_tvar_dd = .true.
        enddo

        !//////////////////////////////////////////////////////
        ! Are we computing magnetic diffusion terms?
        If (check(induct_diff_r) .or. check(induct_diff_bm_r) &
            .or. check(induct_diff_bp_r) ) Then
            l_compute_br_dd=.true.
        Endif

        If (check(induct_diff_theta) .or. check(induct_diff_bm_theta) &
            .or. check(induct_diff_bp_theta) ) Then
            l_compute_bt_dd=.true.
        Endif

        If (check(induct_diff_phi) .or. check(induct_diff_bm_phi) &
            .or. check(induct_diff_bp_phi) ) Then
            l_compute_bp_dd=.true.
        Endif


        If (check(idiff_work) .or. check(idiff_work_pp) &
            .or. check(idiff_work_mm) ) Then
            l_compute_br_dd = .true.
            l_compute_bt_dd = .true.
            l_compute_bp_dd = .true.
        Endif

        !//////////////////////////////////////////////////////
        If (check(curl_v_grad_v_r) .or. check(curl_v_grad_v_r_squared) .or. &
            check(curl_v_grad_v_theta) .or. check(curl_v_grad_v_theta_squared) .or. &
            check(curl_v_grad_v_phi) .or. check(curl_v_grad_v_phi_squared) .or. &
            check(curl_v_grad_v_abs) .or. &
            check(curl_vp_grad_vp_r) .or. check(curl_vp_grad_vp_theta) .or. &
            check(curl_vp_grad_vp_phi) .or. &
            check(curl_vm_grad_vm_r) .or. check(curl_vm_grad_vm_theta) .or. &
            check(curl_vm_grad_vm_phi) .or. &
            check(curl_vp_grad_vm_r) .or. check(curl_vp_grad_vm_theta) .or. &
            check(curl_vp_grad_vm_phi) .or. &
            check(curl_vm_grad_vp_r) .or. check(curl_vm_grad_vp_theta) .or. &
            check(curl_vm_grad_vp_phi) ) Then
            l_compute_vr_dd = .true.
            l_compute_vt_dd = .true.
            l_compute_vp_dd = .true.
        Endif

        If (check(curl_j_cross_b_r) .or. check(curl_j_cross_b_r_squared) .or. &
            check(curl_j_cross_b_theta) .or. check(curl_j_cross_b_theta_squared) .or. &
            check(curl_j_cross_b_phi) .or. check(curl_j_cross_b_phi_squared) .or. &
            check(curl_j_cross_b_abs) .or. &
            check(curl_jp_cross_bp_r) .or. check(curl_jp_cross_bp_theta) .or. &
            check(curl_jp_cross_bp_phi) .or. &
            check(curl_jm_cross_bm_r) .or. check(curl_jm_cross_bm_theta) .or. &
            check(curl_jm_cross_bm_phi) .or. &
            check(curl_jp_cross_bm_r) .or. check(curl_jp_cross_bm_theta) .or. &
            check(curl_jp_cross_bm_phi) .or. &
            check(curl_jm_cross_bp_r) .or. check(curl_jm_cross_bp_theta) .or. &
            check(curl_jm_cross_bp_phi) ) Then
            l_compute_br_dd = .true.
            l_compute_bt_dd = .true.
            l_compute_bp_dd = .true.
        Endif

        ! Execute a lot of compute_q logic here to see if the
        ! different compute_xx variables should be set to true.

        If (l_compute_vr_dd) need_dd = .true.
        If (l_compute_vt_dd) need_dd = .true.
        If (l_compute_vp_dd) need_dd = .true.

        If (l_compute_tvar_dd) need_dd = .true.
        If (l_compute_pvar_dd) need_dd = .true.

        If (l_compute_br_dd) need_dd = .true.
        If (l_compute_bt_dd) need_dd = .true.
        If (l_compute_bp_dd) need_dd = .true.

        ! Turbulent KE generation
        If (check(production_buoyant_pKE)) need_dd = .true.
        If (check(production_shear_pKE)) need_dd = .true.

        If (check(dissipation_viscous_pKE)) need_dd = .true.

        If (check(transport_pressure_pKE)) need_dd = .true.
        If (check(transport_viscous_pKE)) Then
            need_dd = .true.
            l_compute_vr_dd = .true.
            l_compute_vt_dd = .true.
            l_compute_vp_dd = .true.
        Endif
        If (check(transport_turbadvect_pKE)) need_dd = .true.
        If (check(transport_meanadvect_pKE)) need_dd = .true.

        If (check(rflux_pressure_pKE)) need_dd = .true.
        If (check(rflux_viscous_pKE)) need_dd = .true.
        If (check(rflux_turbadvect_pKE)) need_dd = .true.
        If (check(rflux_meanadvect_pKE)) need_dd = .true.

        If (check(thetaflux_pressure_pKE)) need_dd = .true.
        If (check(thetaflux_viscous_pKE)) need_dd = .true.
        If (check(thetaflux_turbadvect_pKE)) need_dd = .true.
        If (check(thetaflux_meanadvect_pKE)) need_dd = .true.

    End Subroutine Second_Derivative_Logic

    Function Second_Derivatives_Needed() result(needed)
        ! Per-iteration determination of whether Compute_Second_Derivatives
        ! needs to run, via compute_quantity (this iteration's menu) rather
        ! than sometimes_compute. The per-field compute_xx flags are
        ! local/transient here -- they do not affect the compute_xx_dd locals
        ! in Initialize_Second_Derivatives, which are fixed once, at startup,
        ! since they size dd_src/d2buffer.
        Implicit None
        Logical :: needed
        Logical :: l_compute_vr_dd, l_compute_vt_dd, l_compute_vp_dd
        Logical :: l_compute_tvar_dd, l_compute_pvar_dd
        Logical :: l_compute_br_dd, l_compute_bt_dd, l_compute_bp_dd

        Call Second_Derivative_Logic(Compute_Quantity, l_compute_vr_dd, l_compute_vt_dd, l_compute_vp_dd, &
                               l_compute_tvar_dd, l_compute_pvar_dd, &
                               l_compute_br_dd, l_compute_bt_dd, l_compute_bp_dd, needed)

    End Function Second_Derivatives_Needed

    Subroutine Initialize_Second_Derivatives()
        ! Initializes all the indexing related to computing and
        ! accessing second derivatives at output time.
        ! Most of the actual indexing is handled by Set_DD_Indices
        IMPLICIT NONE
        INTEGER :: ndind, i
        INTEGER :: ddfcount(3,2)
        Logical :: compute_vr_dd, compute_vt_dd, compute_vp_dd
        Logical :: compute_tvar_dd, compute_pvar_dd
        Logical :: compute_br_dd, compute_bt_dd, compute_bp_dd
        Logical :: need_dd_at_init

        ! Once-at-startup pass: decides, via sometimes_compute (the menu across
        ! the whole run), which fields ever need second derivatives taken, for
        ! sizing dd_src/the d2buffer below. The resulting need_dd_at_init is
        ! discarded; whether Compute_Second_Derivatives runs is refreshed every
        ! iteration by Second_Derivatives_Needed().
        Call Second_Derivative_Logic(Sometimes_Compute, compute_vr_dd, compute_vt_dd, compute_vp_dd, &
                               compute_tvar_dd, compute_pvar_dd, &
                               compute_br_dd, compute_bt_dd, compute_bp_dd, need_dd_at_init)


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




        Allocate(dd_type(nddfields), dd_src(8,nddfields), dd_slot(2,nddfields))
        Allocate(dd_scr(nddfields), dd_sgn(nddfields))
        dd_src(:,:) = -1


        If (compute_vr_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, dd_scalar, &
                                      dvrdrdr, dvrdrdt, dvrdrdp, &
                                      dvrdtdt, dvrdtdp, dvrdpdp, &
                                      vr, dvrdr  , dvrdt  , dvrdp)
        Endif

        IF (compute_vt_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, dd_theta, &
                                      dvtdrdr, dvtdrdt, dvtdrdp, &
                                      dvtdtdt, dvtdtdp, dvtdpdp, &
                                      vtheta, dvtdr  , dvtdt  , dvtdp, &
                                      vphi  , dvpdr  , dvpdt  , dvpdp)
        ENDIF

        IF (compute_vp_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, dd_phi, &
                                      dvpdrdr, dvpdrdt, dvpdrdp, &
                                      dvpdtdt, dvpdtdp, dvpdpdp, &
                                      vphi  , dvpdr  , dvpdt  , dvpdp, &
                                      vtheta, dvtdr  , dvtdt  , dvtdp)
        ENDIF

        IF (compute_tvar_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, dd_scalar, &
                                      dtdrdr, dtdrdt, dtdrdp, &
                                      dtdtdt, dtdtdp, dtdpdp, &
                                      tvar, dtdr  , dtdt  , dtdp)
        ENDIF

        IF (compute_pvar_dd) THEN
            ndind = ndind+1
            Call set_dd_indices(ndind,nddfields, dd_scalar, &
                                      dpdrdr, dpdrdt, dpdrdp, &
                                      dpdtdt, dpdtdp, dpdpdp, &
                                      pvar, dpdr  , dpdt  , dpdp)
        ENDIF

        If (magnetism) THEN
            If (compute_br_dd) THEN
                ndind = ndind+1
                Call set_dd_indices(ndind,nddfields, dd_scalar, &
                                          dbrdrdr, dbrdrdt, dbrdrdp, &
                                          dbrdtdt, dbrdtdp, dbrdpdp, &
                                          br, dbrdr  , dbrdt  , dbrdp)
            Endif

            IF (compute_bt_dd) THEN
                ndind = ndind+1
                Call set_dd_indices(ndind,nddfields, dd_theta, &
                                          dbtdrdr, dbtdrdt, dbtdrdp, &
                                          dbtdtdt, dbtdtdp, dbtdpdp, &
                                          btheta, dbtdr  , dbtdt  , dbtdp, &
                                          bphi  , dbpdr  , dbpdt  , dbpdp)
            ENDIF

            IF (compute_bp_dd) THEN
                ndind = ndind+1
                Call set_dd_indices(ndind,nddfields, dd_phi,      &
                                          dbpdrdr, dbpdrdt, dbpdrdp, &
                                          dbpdtdt, dbpdtdp, dbpdpdp, &
                                          bphi  , dbpdr  , dbpdt  , dbpdp, &
                                          btheta, dbtdr  , dbtdt  , dbtdp)
            ENDIF
        ENDIF

        ! Assign the p3b slots of the fields that are Legendre transformed:
        ! one per scalar (x), two per horizontal component (Q, sin(theta) dx/dr).
        ndd_trans = 0
        ndd_horiz = 0
        dd_slot(:,:) = -1
        dd_scr(:) = -1
        Do i = 1, nddfields
            ndd_trans = ndd_trans+1
            dd_slot(1,i) = ndd_trans
            If (dd_type(i) .ne. dd_scalar) Then
                ndd_trans = ndd_trans+1
                dd_slot(2,i) = ndd_trans
                ndd_horiz = ndd_horiz+1
                dd_scr(i) = ndd_horiz
            Endif
        Enddo

        ! Only ndd_trans fields are transposed on the way in and 3 per variable on the
        ! way back (see Compute_Second_Derivatives).  p3a holds those 3N fields, the 6N
        ! outputs and the 2*ndd_horiz scratch fields.
        ddfcount(1,1) = nddfields*3 ! config 1a
        ddfcount(2,1) = nddfields*3 ! 2a
        ddfcount(3,1) = nddfields*9+ndd_horiz*2 ! 3a
        ddfcount(3,2) = ndd_trans   ! 3b
        ddfcount(2,2) = ndd_trans   ! 2b
        ddfcount(1,2) = ndd_trans   ! 1b


        Call d2buffer%init(field_count = ddfcount, config = 'p3b')

        Call d2buffer%construct('p3a')

        Call d2buffer%deconstruct('p3a')
    End Subroutine Initialize_Second_Derivatives

    Subroutine Compute_Second_Derivatives(inbuffer)
        Implicit None
        INTEGER :: i, j, imi, mp, m
        INTEGER :: r, k, t
        INTEGER :: n, nt, nw, ob, sb, s1, s2, k1, k2, k3, ix, ixr, ixt, iyp
        Real*8 :: sgn
        Real*8, Intent(InOut) :: inbuffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Allocatable :: work(:,:,:,:)
        Type(rmcontainer3D), Allocatable :: ddtemp(:)

        ! Here we compute all second derivatives for N variables.
        !
        ! A field may only be Legendre transformed (forward or inverse) if it is
        ! smooth at the poles, i.e. if its m-th Fourier component behaves like
        ! sin^|m|(theta) * (polynomial in cos(theta)).  Scalars x satisfy this, but
        ! dx/dtheta does not (it has the wrong parity across the pole), and neither
        ! does a horizontal vector component or its radial derivative. 
        ! We only transform smooth quantities and do all division by
        ! sin(theta) pointwise exactly on the grid at the end:
        !
        ! Scalars x in (v_r, T, P, B_r) -- transform C = x only
        !    d2x/dr2       = d2C/dr2
        !    d2x/drdtheta  = [sin(theta) d(dC/dr)/dtheta]/sin(theta), computed spectrally
        !    d2x/dtheta2   = Lap_1(C) - d2x/dphi2/sin^2(theta) - cot(theta) dx/dtheta
        !    where Lap_1 = -l(l+1) is the Laplacian on the unit sphere.
        !
        ! Horizontal components x, with y the other component (v_theta/v_phi, B_theta/B_phi),
        ! sgn = +1 for x = theta component, -1 for x = phi component -- transform
        !    A = sin(theta) dx/dr
        !    Q = dx/dtheta + cot(theta) x + sgn dy/dphi / sin(theta)
        !    Q is r times the horizontal divergence (when x = theta component) or the radial
        !    vorticity (when x = phi component).
        !    Both are smooth scalars for any smooth vector field.
        !    Rearranging,
        !    dx/dtheta = Q - cot(theta) x - sgn dy/dphi / sin(theta), so that
        !    d2x/dr2       = d(A)/dr / sin(theta)
        !    d2x/drdtheta  = dQ/dr - cot(theta) dx/dr - sgn d2y/drdphi / sin(theta)
        !    d2x/dtheta2   = dQ/dtheta + x/sin^2(theta) - cot(theta) dx/dtheta
        !                    + sgn cot(theta) dy/dphi / sin(theta) - sgn d2y/dthetadphi / sin(theta)
        !    with dQ/dtheta = [sin(theta) dQ/dtheta]/sin(theta), computed spectrally.
        !
        ! Truncation at l_max:  multiplying by sin(theta) raises the spectral degree by one,
        ! so A above has an l_max+1 component proportional to the l_max coefficient of x.  The forward
        ! Legendre transform discards it, the result no longer vanishes like sin(theta) at
        ! the poles, and the division by sin(theta) in Step 7 below amplifies the error there.
        ! A is formed in physical space, so it cannot be truncated here; it is exact only
        ! because the inputs come from rlm_spacea (Sphere_Hybrid_Space), where the l_max
        ! mode of every field is zeroed before the transform to physical space.
        ! Any new field added to this routine must satisfy the same condition.
        ! The fields passed to d_by_dtheta (just Q and dC/dr) are truncated at l_max explicitly
        ! (Step 4) for the same reason.
        !
        ! The phi derivatives d2x/drdphi, d2x/dphi2 and d2x/dthetadphi (and d2y/drdphi,
        ! d2y/dthetadphi) are computed from the first derivatives in inbuffer with FFTs
        ! only; no Legendre transform is involved.
        !
        ! Every field in p1a and s2a is transposed so only the three fields needed per variable are kept there;
        ! the radial derivatives are taken in a private work array instead.
        !
        ! Steps:
        ! 1)  Load the p3b slots (see dd_slot) with C, or with Q and A, for each variable
        ! 2)  Move to p1b and copy to the work array (in Chebyshev space if applicable)
        ! 3)  Fill p1a slots [3i-2 : 3i] with, for variable i,
        !        scalar:     C,  dC/dr,  d2C/dr2
        !        horizontal: Q,  dQ/dr,  dA/dr
        ! 4)  Move to s2a; scalar: C -> Lap_1(C), dC/dr -> sin(theta) d(dC/dr)/dtheta;
        !     horizontal: Q -> sin(theta) dQ/dtheta
        ! 5)  Move to p3a
        ! 6)  Load dxdr, dxdp, dxdt into the dxdrdp, dxdpdp, dxdtdp slots (and dydr, dydt into
        !     the scratch slots) and take phi derivatives
        ! 7)  FFT, assemble dxdtdt, dxdrdr and dxdrdt, compute means and fluctuations

        ! When this routine is complete, the contents of d2buffer%p3a will be
        ! [    1 : 3N ] -- transformed fields (workspace)
        ! [3N+1 : 4N ] -- dxdtdt
        ! [4N+1 : 5N ] -- dxdrdr
        ! [5N+1 : 6N ] -- dxdrdt
        ! [6N+1 : 7N ] -- dxdrdp
        ! [7N+1 : 8N ] -- dxdpdp
        ! [8N+1 : 9N ] -- dxdtdp
        ! followed by workspace

        n    = nddfields
        nt   = ndd_trans
        ob   = n*3           ! output base (see Set_DD_Indices)
        sb   = ob+n*6        ! scratch pairs live in [sb+1 : sb+2*ndd_horiz]

        !///////////////////////////////////////////////////////////
        ! Step 1:  Load the fields to be transformed

        Call d2buffer%construct('p3b')
        d2buffer%config = 'p3b'
        d2buffer%p3b(:,:,:,:) = 0.0d0

        Do i = 1, n
            ix  = dd_src(1,i)
            ixr = dd_src(2,i)
            ixt = dd_src(3,i)
            s1  = dd_slot(1,i)
            If (dd_type(i) .eq. dd_scalar) Then
                DO_PSI
                    d2buffer%p3b(PSI,s1) = inbuffer(PSI,ix)
                END_DO
            Else
                s2  = dd_slot(2,i)
                iyp = dd_src(8,i)
                sgn = dd_sgn(i)
                DO_PSI
                    d2buffer%p3b(PSI,s1) = inbuffer(PSI,ixt)+cottheta(t)*inbuffer(PSI,ix) &
                                         + sgn*csctheta(t)*inbuffer(PSI,iyp)
                    d2buffer%p3b(PSI,s2) = inbuffer(PSI,ixr)*sintheta(t)
                END_DO
            Endif
        Enddo


        !////////////////////////////////////////////////////////////////
        ! Step 2:  Move to p1b and copy to the work array
        Call fft_to_spectral(d2buffer%p3b, rsc = .true.)
        Call d2buffer%reform()
        Call d2buffer%construct('s2b')
        Call Legendre_Transform(d2buffer%p2b,d2buffer%s2b)
        Call d2buffer%deconstruct('p2b')
        d2buffer%config ='s2b'

        Call d2buffer%reform() ! move to p1b

        ! The last work slot receives each radial derivative in turn
        nw = nt+1
        Allocate(work(1:size(d2buffer%p1b,1),1:2,1:size(d2buffer%p1b,3),1:nw))
        work(:,:,:,:) = 0.0d0
        If (chebyshev) Then
            Call gridcp%To_Spectral(d2buffer%p1b(:,:,:,1:nt),work(:,:,:,1:nt))
            Call gridcp%dealias_buffer(work(:,:,:,1:nt))
        Else
            work(:,:,:,1:nt) = d2buffer%p1b(:,:,:,1:nt)
        Endif


        !////////////////////////////////////////////////////////////
        ! Step 3:  Fill p1a with the fields needed in s2a and p3a
        Call d2buffer%construct('p1a')
        d2buffer%config='p1a'

        Do i = 1, n
            s1 = dd_slot(1,i)
            k1 = 3*i-2
            k2 = 3*i-1
            k3 = 3*i
            d2buffer%p1a(:,:,:,k1) = work(:,:,:,s1)
            Call radial_derivative(s1, k2, 1)
            If (dd_type(i) .eq. dd_scalar) Then
                Call radial_derivative(s1, k3, 2)
            Else
                Call radial_derivative(dd_slot(2,i), k3, 1)
            Endif
        Enddo

        If (chebyshev) Then
            ! Back to physical space in radius; the work array is no longer needed
            DeAllocate(work)
            Allocate(work(1:size(d2buffer%p1a,1),1:2,1:size(d2buffer%p1a,3),1:size(d2buffer%p1a,4)))
            Call gridcp%From_Spectral(d2buffer%p1a,work)
            d2buffer%p1a = work
        Endif
        DeAllocate(work)
        Call d2buffer%deconstruct('p1b')


        !///////////////////////////////////////////////////////////////
        ! Step 4:  Move to s2a; scalar: C -> Lap_1(C), dC/dr -> sin(theta) d(dC/dr)/dtheta,
        !          horizontal: Q -> sin(theta) dQ/dtheta
        Call d2buffer%reform()

        Call Allocate_rlm_Field(ddtemp)

        Do i = 1, n
            If (dd_type(i) .eq. dd_scalar) Then
                k1 = 3*i-2
                DO_IDX2
                    d2buffer%s2a(mp)%data(IDX2,k1) = -l_l_plus1(m:l_max)*d2buffer%s2a(mp)%data(IDX2,k1)
                END_DO
                Call sintheta_dtheta_inplace(3*i-1)
            Else
                Call sintheta_dtheta_inplace(3*i-2)
            Endif
        Enddo

        Call DeAllocate_rlm_Field(ddtemp)

        Call d2buffer%construct('p2a')
        Call Legendre_Transform(d2buffer%s2a,d2buffer%p2a)
        Call d2buffer%deconstruct('s2a')
        d2buffer%config = 'p2a'


        !/////////////////////////////////////////////////////////////////////
        !  Step 5 : Move to p3a
        Call d2buffer%reform() ! move to p3a

        d2buffer%p3a(:,:,:,ob+1:) = 0.0d0


        !/////////////////////////////////////////////////////////////////////
        !  Step 6 : phi derivatives (FFT only)
        Do i = 1, n
            DO_PSI
                d2buffer%p3a(PSI,ob+i+n*3) = inbuffer(PSI,dd_src(2,i))   ! dxdr -> dxdrdp
                d2buffer%p3a(PSI,ob+i+n*4) = inbuffer(PSI,dd_src(4,i))   ! dxdp -> dxdpdp
                d2buffer%p3a(PSI,ob+i+n*5) = inbuffer(PSI,dd_src(3,i))   ! dxdt -> dxdtdp
            END_DO
            If (dd_type(i) .ne. dd_scalar) Then
                j = sb+dd_scr(i)*2-1
                DO_PSI
                    d2buffer%p3a(PSI,j  ) = inbuffer(PSI,dd_src(6,i)) ! dydr -> dydrdp
                    d2buffer%p3a(PSI,j+1) = inbuffer(PSI,dd_src(7,i)) ! dydt -> dydtdp
                END_DO
            Endif
        Enddo

        ! These are physical; the rest of p3a is already in spectral space (in phi)
        Call fft_to_spectral_rsc(d2buffer%p3a(:,:,:,ob+n*3+1:ob+n*6))
        If (ndd_horiz .gt. 0) Call fft_to_spectral_rsc(d2buffer%p3a(:,:,:,sb+1:sb+ndd_horiz*2))

        Do j = ob+n*3+1, ob+n*6
            Call d_by_dphi(d2buffer%p3a,j,j)
        Enddo
        Do j = sb+1, sb+ndd_horiz*2
            Call d_by_dphi(d2buffer%p3a,j,j)
        Enddo


        !//////////////////////////////////////////
        ! Step 7:   Finalize
        ! FFT
        Call fft_to_physical(d2buffer%p3a,rsc = .true.)

        Do i = 1, n
            ix  = dd_src(1,i)
            ixr = dd_src(2,i)
            ixt = dd_src(3,i)
            k1  = 3*i-2
            k2  = 3*i-1
            k3  = 3*i
            If (dd_type(i) .eq. dd_scalar) Then
                DO_PSI
                    ! d2x/dtheta2   = Lap_1(x) - d2x/dphi2/sin^2(theta) - cot(theta) dx/dtheta
                    d2buffer%p3a(PSI,ob+i    ) = d2buffer%p3a(PSI,k1) &
                                               - csctheta(t)*csctheta(t)*d2buffer%p3a(PSI,ob+i+n*4) &
                                               - cottheta(t)*inbuffer(PSI,ixt)
                    ! d2x/dr2       = d2C/dr2
                    d2buffer%p3a(PSI,ob+i+n  ) = d2buffer%p3a(PSI,k3)
                    ! d2x/drdtheta  = [sin(theta) d(dC/dr)/dtheta]/sin(theta)
                    d2buffer%p3a(PSI,ob+i+n*2) = d2buffer%p3a(PSI,k2)*csctheta(t)
                END_DO
            Else
                iyp = dd_src(8,i)
                sgn = dd_sgn(i)
                j = sb+dd_scr(i)*2-1
                DO_PSI
                    ! d2x/dtheta2 = dQ/dtheta + x/sin^2(theta) - cot(theta) dx/dtheta
                    !         + sgn cot(theta) dy/dphi / sin(theta) - sgn d2y/dthetadphi / sin
                    d2buffer%p3a(PSI,ob+i    ) = d2buffer%p3a(PSI,k1)*csctheta(t) &
                                               + csctheta(t)*csctheta(t)*inbuffer(PSI,ix) &
                                               - cottheta(t)*inbuffer(PSI,ixt) &
                                               + sgn*csctheta(t)*cottheta(t)*inbuffer(PSI,iyp) &
                                               - sgn*csctheta(t)*d2buffer%p3a(PSI,j+1)
                    ! d2x/dr2 = d(sin(theta) dx/dr)/dr / sin(theta)
                    d2buffer%p3a(PSI,ob+i+n  ) = d2buffer%p3a(PSI,k3)*csctheta(t)
                    ! d2x/drdtheta  = dQ/dr - cot(theta) dx/dr - sgn d2y/drdphi / sin(theta)
                    d2buffer%p3a(PSI,ob+i+n*2) = d2buffer%p3a(PSI,k2) &
                                               - cottheta(t)*inbuffer(PSI,ixr) &
                                               - sgn*csctheta(t)*d2buffer%p3a(PSI,j)
                END_DO
            Endif
        Enddo
        !D2buffer is now initialized.  Ordering of output fields is:
        ! dxdtdt, dxdrdr, dxdrdt, dxdrdp, dxdpdp, dxdtdp

        ! Now compute the means and fluctuations
        Allocate(d2_ell0(my_r%min:my_r%max,ob+1:ob+n*6))
        Allocate(d2_m0(my_r%min:my_r%max,my_theta%min:my_theta%max,ob+1:ob+n*6))
        Allocate(d2_fbuffer(1:n_phi,my_r%min:my_r%max, &
                 my_theta%min:my_theta%max,ob+1:ob+n*6))

        Call ComputeEll0(d2buffer%p3a(:,:,:,ob+1:ob+n*6),d2_ell0)
        Call   ComputeM0(d2buffer%p3a(:,:,:,ob+1:ob+n*6),d2_m0)

        DO j = ob+1,ob+n*6
            DO_PSI
                d2_fbuffer(PSI,j) = d2buffer%p3a(PSI,j) - d2_m0(PSI2,j)
            END_DO
        ENDDO

    Contains

        Subroutine radial_derivative(src, dst, dorder)
            ! d^dorder/dr^dorder of work slot src into p1a slot dst
            Integer, Intent(In) :: src, dst, dorder
            If (chebyshev) Then
                Call gridcp%d_by_dr_cp(src, nw, work, dorder)
            Else
                Call d_by_dx3d3(src, nw, work, dorder)
            Endif
            d2buffer%p1a(:,:,:,dst) = work(:,:,:,nw)
        End Subroutine radial_derivative

        Subroutine sintheta_dtheta_inplace(f)
            ! s2a slot f -> sin(theta) d/dtheta of slot f.  That needs modes up to
            ! l_max+1, which s2a cannot hold, so truncate at l_max-1 first so the
            ! result fits exactly and remains divisible by sin(theta) in Step 7
            ! (cf. rlm_spacea in Sphere_Hybrid_Space).
            Integer, Intent(In) :: f
            Do mp = my_mp%min, my_mp%max
                d2buffer%s2a(mp)%data(l_max,:,:,f) = 0.0d0
            Enddo
            Call d_by_dtheta(d2buffer%s2a,f,ddtemp)
            DO_IDX2
                d2buffer%s2a(mp)%data(IDX2,f) = ddtemp(mp)%data(IDX2)
            END_DO
        End Subroutine sintheta_dtheta_inplace

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

    Subroutine Set_DD_Indices(iind, nskip, ftype, &
                                     dxdrdr, dxdrdt, dxdrdp, &
                                     dxdtdt, dxdtdp, dxdpdp, &
                                     x, dxdr, dxdt, dxdp, &
                                     y, dydr, dydt, dydp)
        ! Sets the field type and buffer indices of field iind in dd_type/dd_src and
        ! assigns values to dxdidj consistent with the logic used in
        ! Compute_Second_Derivatives()
        ! [3N+1 : 4N ] -- dxdtdt
        ! [4N+1 : 5N ] -- dxdrdr
        ! [5N+1 : 6N ] -- dxdrdt
        ! [6N+1 : 7N ] -- dxdrdp
        ! [7N+1 : 8N ] -- dxdpdp
        ! [8N+1 : 9N ] -- dxdtdp
        ! (the first 3N slots of d2buffer%p3a hold the transformed fields)
        ! y and its derivatives (the other horizontal component) are required for
        ! horizontal components (ftype = dd_theta or dd_phi) and ignored for scalars.
        Implicit None
        INTEGER, Intent(In)    :: iind, nskip, ftype
        INTEGER, INTENT(OUT)   :: dxdtdt, dxdrdr, dxdrdt
        INTEGER, INTENT(OUT)   :: dxdrdp, dxdpdp, dxdtdp
        INTEGER, INTENT(IN)    :: x, dxdr, dxdt, dxdp
        INTEGER, INTENT(IN), Optional :: y, dydr, dydt, dydp
        dxdtdt = iind+nskip*3
        dxdrdr = iind+nskip*4
        dxdrdt = iind+nskip*5
        dxdrdp = iind+nskip*6
        dxdpdp = iind+nskip*7
        dxdtdp = iind+nskip*8

        dd_type(iind)  = ftype
        dd_src(1,iind) = x
        dd_src(2,iind) = dxdr
        dd_src(3,iind) = dxdt
        dd_src(4,iind) = dxdp
        If (ftype .ne. dd_scalar) Then
            dd_src(5,iind) = y
            dd_src(6,iind) = dydr
            dd_src(7,iind) = dydt
            dd_src(8,iind) = dydp
        Endif
        dd_sgn(iind) = 1.0d0
        If (ftype .eq. dd_phi) dd_sgn(iind) = -1.0d0

    End Subroutine Set_DD_Indices
End Module Diagnostics_Second_Derivatives
