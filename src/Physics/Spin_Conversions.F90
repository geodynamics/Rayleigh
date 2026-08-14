Module Spin_Conversions
    !////////////////////////////////////////////////////////////////////
    !  Slot <-> spin-pair conversions for the compressible horizontal
    !  velocity pair.  Pure table algebra (no transforms, no FFT): kept in
    !  a low-level module so Checkpointing and Initial_Conditions can use
    !  it without depending on Sphere_Hybrid_Space (build-order hygiene).
    !  The conv_alpha factor undoes the FFT-convention (PTS/STP) factors
    !  that the spin tables carry, since these composites are FFT-less;
    !  stp*pts = 1/(2*n_theta) for every m, so the correction is m-uniform.
    !////////////////////////////////////////////////////////////////////
    Use ProblemSize
    Use Load_Balance, Only : mp_lm_values, l_lm_values, my_num_lm, m_lm_values
    Use Legendre_Polynomials, Only : iwp_lm, iwm_lm, wp_lm, wm_lm, ws0_lm, iws0_lm
    Use Structures, Only : rmcontainer4D
    Implicit None
    Logical :: spin_state_is_slot = .false.
Contains

    Subroutine Convert_Slot_Pair_To_Spin(buf, fa, fb)
        ! One-shot after restart: (vtheta,vphi) slots hold componentwise
        ! sintheta-weighted coefficients (checkpoint format).  Convert to native
        ! q+/- using only certified in-tree machinery:
        !   synthesize with p_lm -> divide by sintheta -> rotate to Q+/- ->
        !   analyze with iwp/iwm.
        ! All-l-local per mp in s2a: the l-mixing is contained.
        Implicit None
        Type(rmcontainer4D), Intent(InOut) :: buf(my_mp%min:)
        Integer, Intent(In) :: fa, fb
        Integer :: mp, m, r, imi, i, nl, mloc
        Real*8, Allocatable :: slth(:,:), slph(:,:), qp(:,:), qm(:,:)
        Real*8 :: alpha, beta
        Real*8 :: conv_alpha, vthv, vphv

        alpha = 1.0d0; beta = 0.0d0
        conv_alpha = 2.0d0*n_theta
        Allocate(slth(1:n_theta,1:2), slph(1:n_theta,1:2))
        Allocate(qp(1:n_theta,1:2), qm(1:n_theta,1:2))
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            nl = l_max-m+1
            mloc = mp - my_mp%min + 1   ! Legendre tables are local-indexed (1:n_m)
            Do r = my_r%min, my_r%max
                ! synthesize both slots, both imi, at this radius
                CALL DGEMM('T','N',n_theta,2,nl, alpha, ws0_lm(mloc)%data(m:l_max,:), nl, &
                     buf(mp)%data(m:l_max,r,1:2,fa), nl, beta, slth, n_theta)
                CALL DGEMM('T','N',n_theta,2,nl, alpha, ws0_lm(mloc)%data(m:l_max,:), nl, &
                     buf(mp)%data(m:l_max,r,1:2,fb),   nl, beta, slph, n_theta)
                ! sintheta-divide (slot -> physical v), rotate to Q+/-
                Do i = 1, n_theta
                    vthv = slth(i,1)/sintheta(i); slth(i,1) = vthv
                    vthv = slth(i,2)/sintheta(i); slth(i,2) = vthv
                    vphv = slph(i,1)/sintheta(i); slph(i,1) = vphv
                    vphv = slph(i,2)/sintheta(i); slph(i,2) = vphv
                    qp(i,1) = slth(i,1) - slph(i,2)   ! Q+ re = vth_re - vph_im
                    qp(i,2) = slth(i,2) + slph(i,1)   ! Q+ im
                    qm(i,1) = slth(i,1) + slph(i,2)   ! Q- re
                    qm(i,2) = slth(i,2) - slph(i,1)   ! Q- im
                Enddo
                ! analyze against the spin tables: native q+/-
                ! conv_alpha undoes the FFT-convention factors (stp*pts = 1/(2 n_theta)
                ! for every m) that the tables now carry; the converters are FFT-less.
                CALL DGEMM('T','N',nl,2,n_theta, conv_alpha, iwp_lm(mloc)%data(:,m:l_max), n_theta, &
                     qp, n_theta, beta, buf(mp)%data(m:l_max,r,1:2,fa), nl)
                CALL DGEMM('T','N',nl,2,n_theta, conv_alpha, iwm_lm(mloc)%data(:,m:l_max), n_theta, &
                     qm, n_theta, beta, buf(mp)%data(m:l_max,r,1:2,fb), nl)
            Enddo
        Enddo
        DeAllocate(slth,slph,qp,qm)
    End Subroutine Convert_Slot_Pair_To_Spin

    Subroutine Convert_Spin_Pair_To_Slot(buf, fa, fb)
        ! Inverse of Convert_Slot_Pair_To_Spin, on an arbitrary s2a-layout
        ! container: (fa,fb) slots hold native q+/- coefficients; convert to
        ! the componentwise sintheta-weighted slot representation (checkpoint
        ! format).  Synthesize with wp/wm -> rotate to (vth,vph) [real parts by
        ! construction: real*8 arrays] -> multiply by sintheta -> analyze with
        ! ip_lm.  All-l-local per mp; used by the checkpoint writer.
        Implicit None
        Type(rmcontainer4D), Intent(InOut) :: buf(my_mp%min:)
        Integer, Intent(In) :: fa, fb
        Integer :: mp, m, r, i, nl, mloc
        Real*8, Allocatable :: tQp(:,:), tQm(:,:), slth(:,:), slph(:,:)
        Real*8 :: alpha, beta
        Real*8 :: conv_alpha, vthv, vphv

        alpha = 1.0d0; beta = 0.0d0
        conv_alpha = 2.0d0*n_theta
        Allocate(tQp(1:n_theta,1:2), tQm(1:n_theta,1:2))
        Allocate(slth(1:n_theta,1:2), slph(1:n_theta,1:2))
        Do mp = my_mp%min, my_mp%max
            m = m_values(mp)
            nl = l_max-m+1
            mloc = mp - my_mp%min + 1
            Do r = my_r%min, my_r%max
                CALL DGEMM('T','N',n_theta,2,nl, alpha, wp_lm(mloc)%data(m:l_max,:), nl, &
                     buf(mp)%data(m:l_max,r,1:2,fa), nl, beta, tQp, n_theta)
                CALL DGEMM('T','N',n_theta,2,nl, alpha, wm_lm(mloc)%data(m:l_max,:), nl, &
                     buf(mp)%data(m:l_max,r,1:2,fb), nl, beta, tQm, n_theta)
                Do i = 1, n_theta
                    vthv = 0.5d0*(tQp(i,1)+tQm(i,1))
                    slth(i,1) = sintheta(i)*vthv
                    vthv = 0.5d0*(tQp(i,2)+tQm(i,2))
                    slth(i,2) = sintheta(i)*vthv
                    vphv = 0.5d0*(tQp(i,2)-tQm(i,2))
                    slph(i,1) = sintheta(i)*vphv
                    vphv = 0.5d0*(tQm(i,1)-tQp(i,1))
                    slph(i,2) = sintheta(i)*vphv
                Enddo
                ! conv_alpha: see the reader converter -- undoes the FFT-convention
                ! table factors for this FFT-less composite.
                CALL DGEMM('T','N',nl,2,n_theta, conv_alpha, iws0_lm(mloc)%data(:,m:l_max), n_theta, &
                     slth, n_theta, beta, buf(mp)%data(m:l_max,r,1:2,fa), nl)
                CALL DGEMM('T','N',nl,2,n_theta, conv_alpha, iws0_lm(mloc)%data(:,m:l_max), n_theta, &
                     slph, n_theta, beta, buf(mp)%data(m:l_max,r,1:2,fb), nl)
            Enddo
        Enddo
        DeAllocate(tQp,tQm,slth,slph)
    End Subroutine Convert_Spin_Pair_To_Slot

End Module Spin_Conversions
