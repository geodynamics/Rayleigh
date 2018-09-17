    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Vorticity Outputs
    Integer, Parameter :: vort_off = rhovoffset+100 ! :OFFSET CODE: 

    Integer, Parameter :: vort_r      = vort_off+1  ! :tex: $\omega_r$
    Integer, Parameter :: vort_theta  = vort_off+2  ! :tex: $\omega_\theta$
    Integer, Parameter :: vort_phi    = vort_off+3  ! :tex: $\omega_\phi$

    Integer, Parameter :: vortp_r     = vort_off+4  ! :tex: $\omega_r'$    
    Integer, Parameter :: vortp_theta = vort_off+5  ! :tex: $\omega_\theta'$
    Integer, Parameter :: vortp_phi   = vort_off+6  ! :tex: $\omega_\phi'$

    Integer, Parameter :: vortm_r     = vort_off+7  ! :tex: $\overline{\omega_r}$
    Integer, Parameter :: vortm_theta = vort_off+8  ! :tex: $\overline{\omega_\theta}$
    Integer, Parameter :: vortm_phi   = vort_off+9  ! :tex: $\overline{\omega_\phi}$

    Integer, Parameter :: enstrophy    = vort_off+10 ! :tex: ${\boldsymbol \omega}\cdot{\boldsymbol\omega}$
    Integer, Parameter :: enstrophy_pm = vort_off+11 ! :tex: ${\boldsymbol \omega'}\cdot \overline{{\boldsymbol\omega} } $
    Integer, Parameter :: enstrophy_mm = vort_off+12 ! :tex: $\overline{\boldsymbol \omega}\cdot \overline{{\boldsymbol\omega} } $
    Integer, Parameter :: enstrophy_pp = vort_off+13 ! :tex: ${\boldsymbol \omega'}\cdot {\boldsymbol \omega'} $

    Integer, Parameter :: vort_r_sq      = vort_off+14  ! :tex: $\omega_r^2$
    Integer, Parameter :: vort_theta_sq  = vort_off+15  ! :tex: $\omega_\theta^2$
    Integer, Parameter :: vort_phi_sq    = vort_off+16  ! :tex: $\omega_\phi^2$

    Integer, Parameter :: vortp_r_sq     = vort_off+17  ! :tex: $\omega_r'^2$    
    Integer, Parameter :: vortp_theta_sq = vort_off+18  ! :tex: $\omega_\theta'^2$
    Integer, Parameter :: vortp_phi_sq   = vort_off+19  ! :tex: $\omega_\phi'^2$

    Integer, Parameter :: vortm_r_sq     = vort_off+20  ! :tex: $\overline{\omega_r}^2$
    Integer, Parameter :: vortm_theta_sq = vort_off+21  ! :tex: $\overline{\omega_\theta}^2$
    Integer, Parameter :: vortm_phi_sq   = vort_off+22  ! :tex: $\overline{\omega_\phi}^2$

    Integer, Parameter :: zstream        = vort_off+23  ! :tex: $Z$

    Integer, Parameter :: maria          = vort_off+34  ! :tex: $maria$
