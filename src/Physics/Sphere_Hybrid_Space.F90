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
    Use Math_Constants, Only : sqrt2pi, Ksp
    Use Controls
    Use Legendre_Polynomials, Only : iwp_lm, iwm_lm, wp_lm, wm_lm, ws0_lm, iws0_lm
    Use Spin_Conversions, Only : spin_state_is_slot, Convert_Slot_Pair_To_Spin, Convert_Spin_Pair_To_Slot
    Use ProblemSize
    Use Legendre_Polynomials, Only : p_lm_array
    Use Legendre_Transforms, Only : Legendre_Transform
    Use Spectral_Derivatives
    Use Fields
    Use Timers
    Use ClockInfo
    Use PDE_Coefficients

    Implicit None
    Real*8, Allocatable :: over_rhor(:), over_rhorsq(:), drho_term(:)

    Type(rmcontainer3D), Allocatable :: ftemp1(:), ftemp2(:),ftemp3(:), ftemp4(:)
    ! Spin-horizontal (q+/-) machinery: persistent explicit viscous force for the pair,
    ! assembled spectrally in rlm_spacea (constant-nu v1; dlnu terms TODO before
    ! variable-nu production) and added to the analyzed RHS in rlm_spaceb.
    Type(rmcontainer4D), Allocatable :: fviscq(:)
    Logical :: fviscq_alloc = .false.
    ! One-shot restart conversion: checkpoints store componentwise slot
    ! coefficients (format-compatible); on the first pass after a restart
    ! under spin_horizontal, convert the pair's slots to native q+/-.
Contains


    Subroutine Hybrid_Init()
        Integer :: r1, r2
        ! Allocate a few useful arrays that prevent extra mult/adds
        Allocate(over_rhor(my_r%min:my_r%max))
        Allocate(over_rhorsq(my_r%min:my_r%max))
        Allocate(drho_term(my_r%min:my_r%max))
        r1 = my_r%min
        r2 = my_r%max
        If (pseudo_incompressible) Then
            over_rhor(r1:r2)  = one_over_r(r1:r2)/(ref%density(r1:r2)*ref%exp_entropy(r1:r2))
            over_rhorsq(r1:r2)= OneOverRSquared(r1:r2)/(ref%density(r1:r2)*ref%exp_entropy(r1:r2))
            drho_term(r1:r2)  = ref%dlnrho(r1:r2)+one_over_r(r1:r2) + ref%dsdr_over_cp(r1:r2)
        Else
            over_rhor(r1:r2)  = one_over_r(r1:r2)/ref%density(r1:r2)
            over_rhorsq(r1:r2)= OneOverRSquared(r1:r2)/ref%density(r1:r2)
            drho_term(r1:r2)  = ref%dlnrho(r1:r2)+one_over_r(r1:r2)
        Endif

    End Subroutine Hybrid_Init

    Subroutine rlm_spacea()
        Implicit None
        Integer :: mp, i, r, imi, m, l, ind_top
        Call StopWatch(rlma_time)%startclock()

        ! Zero out l_max mode
        Do mp = my_mp%min, my_mp%max
            SBUFFA(l_max,:,:,:) = 0.0d0
        Enddo

        ! Allocate two work arrays
        Call Allocate_rlm_Field(ftemp1)
        Call Allocate_rlm_Field(ftemp2)


        If (output_iteration) Call Hybrid_Output_Initial()
        
        If (compressible) Then
           If (spin_horizontal .and. spin_state_is_slot) Then
              Call Convert_Slot_Pair_To_Spin(wsp%s2a, vtheta, vphi)
              spin_state_is_slot = .false.
              If (my_rank .eq. 0) Write(6,*) 'spin_horizontal: restart state converted slot -> q+/-'
           Endif
           If (.not. spin_horizontal) Call M0_Pole_Projection()
           ! First theta-derivatives...
           Call d_by_dtheta(wsp%s2a  ,     vr, dvrdt)
           If (.not. spin_horizontal) Then
              Call d_by_dtheta(wsp%s2a  , vtheta, dvtdt)
              Call d_by_dtheta(wsp%s2a  ,   vphi, dvpdt)
           Endif
           Call d_by_dtheta(wsp%s2a  ,   tvar, dtdt)
           Call d_by_dtheta(wsp%s2a  , rhovar, drhodt)

           Call d_by_dtheta(wsp%s2a  , dvrdr, d2vrdrdt)
           If (.not. spin_horizontal) Then
              Call d_by_dtheta(wsp%s2a  , dvtdr, d2vtdrdt)

              Call d_by_dtheta(wsp%s2a  , dvtdt, d2vtdt2)
           Else
              ! Spin representation: the (dvtdt,dvpdt) slots receive COPIES of the
              ! native q+/- coefficients; the Legendre synthesis applies the
              ! d(W+/-)/dtheta tables to them, delivering physical (dvth/dth,dvph/dth)
              ! pole-cleanly.  d2vtdt2 has no consumers under spin (its chain-rule
              ! and divu(3) consumers are all in .not.spin_horizontal regions);
              ! zero it explicitly so the slot never carries stale/uninitialized
              ! content through the transform.
              Do mp = my_mp%min, my_mp%max
                 wsp%s2a(mp)%data(:,:,:,dvtdt)    = wsp%s2a(mp)%data(:,:,:,vtheta)
                 wsp%s2a(mp)%data(:,:,:,dvpdt)    = wsp%s2a(mp)%data(:,:,:,vphi)
                 wsp%s2a(mp)%data(:,:,:,d2vtdt2)  = 0.0d0
              Enddo
              Call Assemble_Spin_Viscous()   ! fviscq for the pair; also supplies
                                             ! d_r(div_h) -> d2vtdrdt (divu(:,2))
                                             ! and zeros hvtheta/hvphi
           Endif

           ! Compute l_l_plus/r^2 for each term
           ! vr
           DO_IDX2
           SBUFFA(IDX2,hvr) = l_l_plus1(m:l_max)*SBUFFA(IDX2,vr)*OneOverRsquared(r)
           END_DO

           If (.not. spin_horizontal) Then
              ! vtheta
              DO_IDX2
              SBUFFA(IDX2,hvtheta) = l_l_plus1(m:l_max)*SBUFFA(IDX2,vtheta)*OneOverRsquared(r)
              END_DO

              ! vphi
              DO_IDX2
              SBUFFA(IDX2,hvphi) = l_l_plus1(m:l_max)*SBUFFA(IDX2,vphi)*OneOverRSquared(r)
              END_DO
           Else
              Continue   ! hvtheta/hvphi carry the grad(div_h) supply, loaded
                         ! in Assemble_Spin_Viscous (which runs above).
           Endif

           ! T
           DO_IDX2
           SBUFFA(IDX2,htvar) = l_l_plus1(m:l_max)*SBUFFA(IDX2,tvar)*OneOverRSquared(r)
           END_DO
        Else
           Call Velocity_Components()
           Call Velocity_Derivatives()
           Call d_by_dtheta(wsp%s2a,tvar,dtdt)
           do i = 1, n_active_scalars
              Call d_by_dtheta(wsp%s2a,chiavar(i),dchiadt(i))
           end do
           do i = 1, n_passive_scalars
              Call d_by_dtheta(wsp%s2a,chipvar(i),dchipdt(i))
           end do
        Endif 
        
        If (output_iteration) Call Hybrid_Output_Final()

        If (magnetism) Call compute_BandCurlB()

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
        If (compressible .and. spin_horizontal) Then
            Call Legendre_Transform(wsp%s2a,wsp%p2a, &
                spin_fa=(/vtheta,dvtdr,d2vtdr2,hvtheta/), spin_fb=(/vphi,dvpdr,d2vpdr2,hvphi/), &
                spin_da=(/dvtdt/), spin_db=(/dvpdt/))
        Else
            Call Legendre_Transform(wsp%s2a,wsp%p2a)
        Endif
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



    Subroutine Assemble_Spin_Viscous()
        ! Explicit spectral viscous force for the compressible horizontal pair
        ! (native q+/- representation), per the certified operator table:
        !   F+/- = nu(r) [ d2q/dr2 + (2/r) dq/dr - L2/r^2 q ]
        !          - nu(r) (2/r^2) sqrt(L2) vr_lm
        !          - (nu(r)/3) (1/r) sqrt(L2) div_lm ,
        !   div_lm = (sqrt(L2)/(2r)) (q+ - q-)
        ! Constant-nu v1: dlnu (viscosity-variation) terms are TODO before
        ! variable-nu production runs.  Signs certified in spin/battery_final.py;
        ! the code's CS-phased Pbar tables match the paper's Y0, so no phase
        ! factor appears in the vr/div couplings.
        !
        ! TIER-1 CONSISTENCY: when implicit_compressible_diffusion is on, the
        ! CN matrices apply nu[D2 + (2/r+dlnrho)D1 - (dlnrho/r)D0] to the pair
        ! rows (SLT loader) acting on the q+/- coefficients.  We subtract that
        ! operator SPECTRALLY here (same slots the matrices act on) so the
        ! explicit remainder is exact by construction.  The physical-space
        ! Subtract_Implicit_Diffusion_V pair blocks are gated off under spin:
        ! subtracting there and analyzing through the csc-prescaled transform
        ! applies the operator to (csc-roundtrip)-mapped coefficients, leaving
        ! a residual nu*k_r^2-scale operator that destabilizes driven runs
        ! (measured: blowup at t~150 at 48x64, dt-independent; clean with the
        ! spectral subtraction or with explicit diffusion).
        Implicit None
        Integer :: mp, m, r, imi, l
        Real*8  :: L2, sL, divl, oor, oor2, crossg

        If (.not. fviscq_alloc) Then
            Allocate(fviscq(my_mp%min:my_mp%max))
            Do mp = my_mp%min, my_mp%max
                m = m_values(mp)
                Allocate(fviscq(mp)%data(m:l_max,my_r%min:my_r%max,1:2,1:2))
            Enddo
            fviscq_alloc = .true.
        Endif

        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            Do imi = 1, 2
            Do r = my_r%min, my_r%max
                oor  = one_over_r(r)
                oor2 = OneOverRsquared(r)
                Do l = m, l_max
                    L2 = l*(l+1.0d0)
                    sL = sqrt(L2)
                    ! FULL divergence coefficient: horizontal part from the pair
                    ! + the vr radial part (dvrdr + 2 vr / r).  Using div_h alone
                    ! here (the pre-v7 bug) drops the vr->pair block of the
                    ! nu/3 grad(div) coupling while the vr equation keeps its
                    ! pair->vr block (divu(:,2) supply): the one-directional
                    ! coupling destabilizes the acoustic-dilatational family at
                    ! sigma ~ nu*kr*kh (deaths: 459 its at 128x192, 2657 at
                    ! 64x96 -- ratio 5.8 vs the coefficient ratio 5.6).
                    ! In-tree horizontal-divergence coefficient: 2*sqrt(L2)/r *
                    ! (q+ - q-).  MEASURED against the true physical divergence
                    ! (utest3: consistent injected state; div_l0 from the ungated
                    ! componentwise divu(:,1) vs the coefficient here): the paper
                    ! identity's (sqrt(L2)/2r) transcribes to 4x that in the
                    ! code's q+/- convention (ratio measured 3.94-4.08 across
                    ! radii, vr channel calibrated to 5%).  The same factor
                    ! applies to the d_r(div_h) supply below.
                    ! TRUE scalar-convention divergence.  The pair part's
                    ! constant is K*sqrt(L2)/r with K = 1/(2*sqrt(2*pi)) =
                    ! Ksp = 1/(2 sqrt(2 pi)): the synthesis-kernel constant
                    ! (multi-m plane scan, ratio flat 0.0996 vs the old 2*sL
                    ! form at every m from 0 to 31, i.e. true = K*sL exactly;
                    ! K pinned to 1/(2*sqrt(2*pi)) at 1e-5 via the exact
                    ! radius).  The earlier x4 (2*sL) calibration was a
                    ! normalization bug in the transposed-era analysis chain
                    ! and is retired.  This divl feeds BOTH the d_r(div_h)
                    ! supply (scalar consumers: must be true-convention) and
                    ! the fviscq div force below (which carries the exact
                    ! reverse bridge sqrt(2*pi), so the pair->pair round trip
                    ! is convention-free).
                    divl = Ksp*sL*oor*( wsp%s2a(mp)%data(l,r,imi,vtheta) &
                                        - wsp%s2a(mp)%data(l,r,imi,vphi) ) &
                         + wsp%s2a(mp)%data(l,r,imi,dvrdr) &
                         + 2.0d0*oor*wsp%s2a(mp)%data(l,r,imi,vr)
                    ! Cross terms (vr coupling and grad-div): these are
                    ! nabla-of-scalar-class forces, whose spin structure is
                    ! ANTISYMMETRIC across the q+/q- slots (the edth-sign),
                    ! NOT equal-slot.  Patterns and magnitudes MEASURED
                    ! end-to-end (utest planes, m=0 kernel: vth = K(q- - q+),
                    ! machine-clean single mode; true forces from the ungated
                    ! physical divu/fields):
                    !   vr term:  a± = ∓ 1.000 * nu (2/r^2) sqrt(L2) vr_lm
                    !             (measured coef 0.9965)
                    !   div term: a± = ∓ 0.250 * (nu/3)(sqrt(L2)/r) divl
                    !             (measured coef 0.2478; divl = TRUE scalar-
                    !              convention divergence, the x4 form below)
                    ! The equal-slot coding synthesized ZERO theta-force from
                    ! these terms (the (a+ + a-) channel), which both broke the
                    ! vr<->pair grad-div symmetry (the gate-G collapse family)
                    ! and dropped 28% of the true pair viscous force
                    ! (F_syn/F_true = 0.72 measured pre-fix).
                    If (implicit_compressible_acoustics) Then
                        ! The cross terms are implicit (coupled block);
                        ! nothing explicit remains of them.
                        crossg = 0.0d0
                    Else
                        crossg = nu(r)*2.0d0*oor2*sL*wsp%s2a(mp)%data(l,r,imi,vr) &
                               + sqrt2pi*(nu(r)/3.0d0)*oor*sL*divl
                    Endif
                    If (implicit_compressible_diffusion) Then
                        ! The T1 CN rows own the full
                        ! nu[D2+(2/r+dlnrho)D1+(Hlap-dlnrho/r)D0] pair operator,
                        ! so the explicit remainder is assembled directly as
                        !   -nu*dlnrho*(D1 - 1/r) q+/-   (l-independent).
                        fviscq(mp)%data(l,r,imi,1) = -nu(r)*ref%dlnrho(r)*( &
                              wsp%s2a(mp)%data(l,r,imi,dvtdr) &
                            - oor*wsp%s2a(mp)%data(l,r,imi,vtheta) ) &
                            - crossg
                        fviscq(mp)%data(l,r,imi,2) = -nu(r)*ref%dlnrho(r)*( &
                              wsp%s2a(mp)%data(l,r,imi,dvpdr) &
                            - oor*wsp%s2a(mp)%data(l,r,imi,vphi) ) &
                            + crossg
                    Else
                        fviscq(mp)%data(l,r,imi,1) = nu(r)*( wsp%s2a(mp)%data(l,r,imi,d2vtdr2) &
                              + 2.0d0*oor*wsp%s2a(mp)%data(l,r,imi,dvtdr) &
                              - L2*oor2*wsp%s2a(mp)%data(l,r,imi,vtheta) ) &
                              - crossg
                        fviscq(mp)%data(l,r,imi,2) = nu(r)*( wsp%s2a(mp)%data(l,r,imi,d2vpdr2) &
                              + 2.0d0*oor*wsp%s2a(mp)%data(l,r,imi,dvpdr) &
                              - L2*oor2*wsp%s2a(mp)%data(l,r,imi,vphi) ) &
                              + crossg
                    Endif
                    ! d_r(div_h) supply: the d2vtdrdt slot (scalar synthesis path)
                    ! feeds the spin-branch divu(:,2) = d_r(div u), whose consumer
                    ! is the ungated vr-equation nu/3 term:
                    If (implicit_compressible_acoustics) Then
                        ! The d_r(div_h) supply is implicit (vreq <- pair
                        ! columns); the explicit slot must carry zero.
                        wsp%s2a(mp)%data(l,r,imi,d2vtdrdt) = 0.0d0
                    Else
                        wsp%s2a(mp)%data(l,r,imi,d2vtdrdt) = Ksp*sL*( &
                              ( wsp%s2a(mp)%data(l,r,imi,dvtdr) &
                              - wsp%s2a(mp)%data(l,r,imi,dvpdr) )*oor &
                              - ( wsp%s2a(mp)%data(l,r,imi,vtheta) &
                              - wsp%s2a(mp)%data(l,r,imi,vphi) )*oor2 )
                    Endif
                    ! hvtheta/hvphi fed only the theta/phi nu/3 blocks, which are
                    ! gated off under spin (divu(3:4) consumers all in the OFF-only
                    ! region).  The prior fill here was a half-implemented gradient
                    ! route (wrong slots for the W+/- pair mechanism): honest zeros.
                    wsp%s2a(mp)%data(l,r,imi,hvtheta) = 0.0d0
                    wsp%s2a(mp)%data(l,r,imi,hvphi)   = 0.0d0
                Enddo
            Enddo
            Enddo
        Enddo
    End Subroutine Assemble_Spin_Viscous

    Subroutine M0_Pole_Projection()
        ! Remove the pole-singular m=0 sector from the sintheta-weighted
        ! vtheta and vphi slot fields.  Physical v = synth(slot)/sin(theta).
        ! For m=0, slot content with nonzero pole value gives v ~ 1/sin(theta):
        ! unphysical.  Enforce synth(slot)(costheta = +/-1) = 0 at every radius
        ! via the minimal-norm coefficient correction; the two pole constraints
        ! decouple by parity (even/odd l):
        !     c_l <- c_l - (sum_class c_l*w_l / sum_class w_l^2) * w_l
        ! with w_l = sqrt(2l+1).  The global normalization of Pbar_l0(+/-1)
        ! cancels in the ratio; only the relative l-weighting matters.
        ! Loop bound l_max-1: the l = l_max slot never reaches physical space
        ! (measured), so it is excluded from sums and corrections.
        ! Regular m=0 physics (e.g. meridional circulation) has zero pole
        ! sum already: the projection is a strict no-op on it.
        Implicit None
        Integer :: mp, r, l, imi, f, fld
        Real*8  :: pe, po, we, wo

        Do mp = my_mp%min, my_mp%max
            If (m_values(mp) .ne. 0) Cycle
            Do f = 1, 2
                If (f .eq. 1) fld = vtheta
                If (f .eq. 2) fld = vphi
                Do imi = 1, 2
                Do r = my_r%min, my_r%max
                    pe = 0.0d0; po = 0.0d0; we = 0.0d0; wo = 0.0d0
                    Do l = 0, l_max-1
                        If (mod(l,2) .eq. 0) Then
                            pe = pe + wsp%s2a(mp)%data(l,r,imi,fld)*sqrt(2.0d0*l+1.0d0)
                            we = we + (2.0d0*l+1.0d0)
                        Else
                            po = po + wsp%s2a(mp)%data(l,r,imi,fld)*sqrt(2.0d0*l+1.0d0)
                            wo = wo + (2.0d0*l+1.0d0)
                        Endif
                    Enddo
                    Do l = 0, l_max-1
                        If (mod(l,2) .eq. 0) Then
                            wsp%s2a(mp)%data(l,r,imi,fld) = wsp%s2a(mp)%data(l,r,imi,fld) &
                                                          - (pe/we)*sqrt(2.0d0*l+1.0d0)
                        Else
                            wsp%s2a(mp)%data(l,r,imi,fld) = wsp%s2a(mp)%data(l,r,imi,fld) &
                                                          - (po/wo)*sqrt(2.0d0*l+1.0d0)
                        Endif
                    Enddo
                Enddo
                Enddo
            Enddo
        Enddo
    End Subroutine M0_Pole_Projection
    
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
        If (compressible .and. spin_horizontal) Then
            ! The momentum RHS components (Ftheta,Fphi) analyze as a spin pair;
            ! the equation slots then hold native F+/- for the q+/- updates.
            ! spin_sindiv MUST stay .false.: the physical-space assembly
            ! delivers true unweighted (Ftheta,Fphi), and analysis with
            ! sindiv off is the exact inverse of the unweighted synthesis.
            ! With sindiv on, every force term is analyzed as F/sintheta,
            ! injecting a broad-l pole-band error each step (worst at m=1)
            ! that drives an l-tail instability.
            Call Legendre_Transform(wsp%p2b,wsp%s2b, spin_fa=(/vtheta/), spin_fb=(/vphi/), spin_sindiv=.false.)
        Else
            Call Legendre_Transform(wsp%p2b,wsp%s2b)
        Endif
        Call StopWatch(legendre_time)%increment()

        Call wsp%deconstruct('p2b')
        wsp%config = 's2b'
        Call StopWatch(rlmb_time)%startclock()


        If (.NOT. compressible) Then
           ! The NL RHS for W is r^2/(l(l+1)) * the NL RHS for Ur
           ! We already have the r^2 taken care of.  Now for the l(l+1)

           DO_IDX2
           SBUFFB(IDX2,wvar) = SBUFFB(IDX2,wvar)*over_l_l_plus1(m:l_max)
           END_DO
        Endif 

        ! Now for the Z RHS, formed from the radial component of the curl of u dot grad u
        If (compressible) Then
            If (spin_horizontal) Then
                ! Add the spectrally assembled explicit viscous force for the pair
                Do mp = my_mp%min, my_mp%max
                    m = m_values(mp)
                    Do imi = 1, 2
                    Do r = my_r%min, my_r%max
                        wsp%s2b(mp)%data(m:l_max,r,imi,vtheta) = wsp%s2b(mp)%data(m:l_max,r,imi,vtheta) &
                                                               + fviscq(mp)%data(m:l_max,r,imi,1)
                        wsp%s2b(mp)%data(m:l_max,r,imi,vphi)   = wsp%s2b(mp)%data(m:l_max,r,imi,vphi) &
                                                               + fviscq(mp)%data(m:l_max,r,imi,2)
                    Enddo
                    Enddo
                Enddo
            Endif
            If (magnetism) Then
                Call Allocate_rlm_Field(ftemp1)
                Call Allocate_rlm_Field(ftemp2)
                Call adjust_emf()
                Call DeAllocate_rlm_Field(ftemp1)
                Call DeAllocate_rlm_Field(ftemp2)
            Endif
        Else
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

        Endif 

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
        Integer :: r, m, mp, imi
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
        Integer :: r, m, mp, imi
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
        ! for efficiency later.
        ! The pseudo-incompressible correction has been implemented through drho_term.

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
        ! ...  pseudo-incompressible correction has been implemented through drho_term.
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
        If (pseudo_incompressible) Then
            DO_IDX2
                SBUFFA(IDX2,dvrdr) = SBUFFA(IDX2,dvrdr)- &
                    & SBUFFA(IDX2,vr)*ref%dsdr_over_cp(r)
            END_DO
        Endif        


        Call d_by_dtheta(wsp%s2a,vr,dvrdt)


        ! Convert Z to ell(ell+1) Z/r^2  (i.e. omega_r)
        ! Computing the radial component of the curl of the velocity
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
        Integer :: r, m, mp, imi

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
        Integer :: r, m, mp, imi
        If (compressible) Call cobuffer%construct("s2a")
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
        Integer :: mp
        Do mp = my_mp%min, my_mp%max
            ASBUFFA(l_max,:,:,:) = 0.0d0
        Enddo
        if (.NOT. compressible) Call Hydro_Output_Derivatives()
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
