    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Vorticity Outputs
    Integer, Parameter :: vort_off = rhovoffset+100 ! :OFFSET CODE:

    Integer, Parameter :: vort_r      = vort_off+1  ! :tex: $\omega_r$
    Integer, Parameter :: vort_theta  = vort_off+2  ! :tex: $\omega_\theta$
    Integer, Parameter :: vort_phi    = vort_off+3  ! :tex: $\omega_\phi$

    Integer, Parameter :: vortp_r     = vort_off+4  ! :tex: $\omega_r^\prime$
    Integer, Parameter :: vortp_theta = vort_off+5  ! :tex: $\omega_\theta^\prime$
    Integer, Parameter :: vortp_phi   = vort_off+6  ! :tex: $\omega_\phi^\prime$

    Integer, Parameter :: vortm_r     = vort_off+7  ! :tex: $\overline{\omega_r}$
    Integer, Parameter :: vortm_theta = vort_off+8  ! :tex: $\overline{\omega_\theta}$
    Integer, Parameter :: vortm_phi   = vort_off+9  ! :tex: $\overline{\omega_\phi}$

    Integer, Parameter :: enstrophy    = vort_off+10 ! :tex: ${\boldsymbol \omega}\cdot{\boldsymbol\omega}$
    Integer, Parameter :: enstrophy_pm = vort_off+11 ! :tex: ${\boldsymbol \omega^\prime}\cdot \overline{{\boldsymbol\omega} } $
    Integer, Parameter :: enstrophy_mm = vort_off+12 ! :tex: $\overline{\boldsymbol \omega}\cdot \overline{{\boldsymbol\omega} } $
    Integer, Parameter :: enstrophy_pp = vort_off+13 ! :tex: ${\boldsymbol \omega^\prime}\cdot {\boldsymbol \omega^\prime} $

    Integer, Parameter :: vort_r_sq      = vort_off+14  ! :tex: $\omega_r^2$
    Integer, Parameter :: vort_theta_sq  = vort_off+15  ! :tex: $\omega_\theta^2$
    Integer, Parameter :: vort_phi_sq    = vort_off+16  ! :tex: $\omega_\phi^2$

    Integer, Parameter :: vortp_r_sq     = vort_off+17  ! :tex: $\omega_r^\prime^2$ 
    Integer, Parameter :: vortp_theta_sq = vort_off+18  ! :tex: $\omega_\theta^\prime^2$
    Integer, Parameter :: vortp_phi_sq   = vort_off+19  ! :tex: $\omega_\phi^\prime^2$

    Integer, Parameter :: vortm_r_sq     = vort_off+20  ! :tex: $\overline{\omega_r}^2$
    Integer, Parameter :: vortm_theta_sq = vort_off+21  ! :tex: $\overline{\omega_\theta}^2$
    Integer, Parameter :: vortm_phi_sq   = vort_off+22  ! :tex: $\overline{\omega_\phi}^2$

    Integer, Parameter :: zstream        = vort_off+23  ! :tex: $Z$

    Integer, Parameter :: kin_helicity_r     = vort_off+24 ! :tex: $v_r\omega_r$
    Integer, Parameter :: kin_helicity_theta = vort_off+25 ! :tex: $v_\theta\omega_\theta$
    Integer, Parameter :: kin_helicity_phi   = vort_off+26 ! :tex: $v_\phi\omega_\phi$

    Integer, Parameter :: kin_helicity_pp_r     = vort_off+27 ! :tex: $v_r^\prime\omega_r^\prime$
    Integer, Parameter :: kin_helicity_pp_theta = vort_off+28 ! :tex: $v_\theta^\prime\omega_\theta^\prime$
    Integer, Parameter :: kin_helicity_pp_phi   = vort_off+29 ! :tex: $v_\phi^\prime\omega_\phi^\prime$

    Integer, Parameter :: kin_helicity_mm_r     = vort_off+30 ! :tex: $\overline{v_r}\overline{\omega_r}$
    Integer, Parameter :: kin_helicity_mm_theta = vort_off+31 ! :tex: $\overline{v_\theta}\overline{\omega_\theta}$
    Integer, Parameter :: kin_helicity_mm_phi   = vort_off+32 ! :tex: $\overline{v_\phi}\overline{\omega_\phi}$

    Integer, Parameter :: kin_helicity_mp_r     = vort_off+33 ! :tex: $\overline{v_r}\omega_r^\prime$
    Integer, Parameter :: kin_helicity_mp_theta = vort_off+34 ! :tex: $\overline{v_\theta}\omega_\theta^\prime$
    Integer, Parameter :: kin_helicity_mp_phi   = vort_off+35 ! :tex: $\overline{v_\phi}\omega_\phi^\prime$

    Integer, Parameter :: kin_helicity_pm_r     = vort_off+36 ! :tex: $v_r^\prime\overline{\omega_r}$
    Integer, Parameter :: kin_helicity_pm_theta = vort_off+37 ! :tex: $v_\theta^\prime\overline{\omega_\theta}$
    Integer, Parameter :: kin_helicity_pm_phi   = vort_off+38 ! :tex: $v_\phi^\prime\overline{\omega_\phi}$

    Integer, Parameter :: kin_helicity          = vort_off+39 ! :tex: $v\cdot\omega$
    Integer, Parameter :: kin_helicity_pp       = vort_off+40 ! :tex: $v^\prime\cdot\omega^\prime$
    Integer, Parameter :: kin_helicity_mm       = vort_off+41 ! :tex: $\overline{v}\cdot\overline{\omega}$
    Integer, Parameter :: kin_helicity_mp       = vort_off+42 ! :tex: $\overline{v}\cdot\omega^\prime$
    Integer, Parameter :: kin_helicity_pm       = vort_off+43 ! :tex: $v^\prime\cdot\overline{\omega}$

