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

    !//////////////////////////////////////////////////////////////////
    !  Composition/Scalar (Active)
    Integer, Parameter :: scalar_skip     = 200
    Integer, Parameter :: a_scalar_offset = 10000 ! :OFFSET CODE:
    Integer, Parameter :: p_scalar_offset = 20000 ! :OFFSET CODE:

    Integer, Parameter :: chi_off = 0 ! :OFFSET CODE:
    !------------ Fields ---------------------!
    ! Full
    Integer, Parameter :: chi   = chi_off + 1 !  :tex:  $\chi$
    ! Fluctuating
    Integer, Parameter :: chi_p = chi_off + 2 !  :tex:  $\chi'$
    ! Mean
    Integer, Parameter :: chi_m = chi_off + 3 !  :tex:  $\overline{\chi}$

    !------------ Radial Derivatives --------------!
    !  Full
    Integer, Parameter :: chi_dr   = chi_off + 4  ! :tex:  $\frac{\partial \chi}{\partial r}$
    ! Fluctuating
    Integer, Parameter :: chi_p_dr = chi_off + 5  ! :tex:  $\frac{\partial \chi'}{\partial r}$
    ! Mean
    Integer, Parameter :: chi_m_dr = chi_off + 6 ! :tex:  $\frac{\partial \overline{\chi}}{\partial r}$

    !------------ Theta Derivatives ---------------!
    Integer, Parameter :: chi_dtheta   = chi_off + 7  ! :tex:  $\frac{\partial \chi}{\partial \theta}$
    Integer, Parameter :: chi_p_dtheta = chi_off + 8  ! :tex:  $\frac{\partial \chi'}{\partial \theta}$
    Integer, Parameter :: chi_m_dtheta = chi_off + 9  ! :tex:  $\frac{\partial \overline{\chi}}{\partial \theta}$

    !------------ Phi Derivatives -----------------!
    Integer, Parameter :: chi_dphi     = chi_off + 10  ! :tex:  $\frac{\partial \chi}{\partial \phi}$
    Integer, Parameter :: chi_p_dphi   = chi_off + 11  ! :tex:  $\frac{\partial \chi'}{\partial \phi}$
    Integer, Parameter :: chi_m_dphi   = chi_off + 12  ! :tex:  $\frac{\partial \overline{\chi}}{\partial \phi}$


    ! FIXME: not implemented, placeholder
    ! !------------ (1/r) * Theta Derivatives --------!
    ! Integer, Parameter :: chi_dtr     = chi_off + 13 !  :tex:  $\frac{1}{r}\frac{\partial \chi}{\partial \theta}$
    ! Integer, Parameter :: chi_p_dtr   = chi_off + 14 !  :tex:  $\frac{1}{r}\frac{\partial \chi'}{\partial \theta}$
    ! Integer, Parameter :: chi_m_dtr   = chi_off + 15 !  :tex:  $\frac{1}{r}\frac{\partial \overline{\chi}}{\partial \theta}$
    ! !--- (1/{r sintheta}) * Phi Derivatives ---------!
    ! Integer, Parameter :: chi_dprs     = chi_off + 16 ! :tex: $\frac{1}{r\mathrm{sin}\theta} \frac{\partial \chi}{\partial \phi}$
    ! Integer, Parameter :: chi_p_dprs   = chi_off + 17 ! :tex: $\frac{1}{r\mathrm{sin}\theta} \frac{\partial \chi'}{\partial \phi}$
    ! Integer, Parameter :: chi_m_dprs   = chi_off + 18 ! :tex: $\frac{1}{r\mathrm{sin}\theta} \frac{\partial \overline{\chi}}{\partial \phi}$
    ! !------------ Radial Second Derivatives --------------!
    ! Integer, Parameter :: chi_d2r     = chi_off + 19  ! :tex:  $\frac{\partial^2 \chi}{\partial r^2}$
    ! Integer, Parameter :: chi_p_d2r   = chi_off + 20  ! :tex:  $\frac{\partial^2 \chi'}{\partial r^2}$
    ! Integer, Parameter :: chi_m_d2r   = chi_off + 21  !  :tex:  $\frac{\partial^2 \overline{\chi}}{\partial r^2}$
    ! !------------ Theta Second Derivatives --------------!
    ! Integer, Parameter :: chi_d2t     = chi_off + 22  ! :tex:  $\frac{\partial^2 \chi}{\partial \theta^2}$
    ! Integer, Parameter :: chi_p_d2t   = chi_off + 23  ! :tex:  $\frac{\partial^2 \chi'}{\partial \theta^2}$
    ! Integer, Parameter :: chi_m_d2t   = chi_off + 24  ! :tex:  $\frac{\partial^2 \overline{\chi}}{\partial \theta^2}$
    ! !------------ Phi Second Derivatives --------------!
    ! Integer, Parameter :: chi_d2p     = chi_off + 25 ! :tex:  $\frac{\partial^2 \chi}{\partial \phi^2}$
    ! Integer, Parameter :: chi_p_d2p   = chi_off + 26  ! :tex:  $\frac{\partial^2 \chi'}{\partial \phi^2}$
    ! Integer, Parameter :: chi_m_d2p   = chi_off + 27  ! :tex:  $\frac{\partial^2 \overline{\chi}}{\partial \phi^2}$
    ! !------------ Radial-Theta Second Derivatives --------------!
    ! Integer, Parameter :: chi_d2rt     = chi_off + 28  ! :tex:  $\frac{\partial^2 \chi}{\partial r \partial \theta}$
    ! Integer, Parameter :: chi_p_d2rt   = chi_off + 29  ! :tex:  $\frac{\partial^2 \chi'}{\partial r \partial \theta}$
    ! Integer, Parameter :: chi_m_d2rt   = chi_off + 30  ! :tex:  $\frac{\partial^2 \overline{\chi}}{\partial r \partial \theta}$
    ! !------------ Radial-Phi Derivatives --------------!
    ! Integer, Parameter :: chi_d2rp     = chi_off + 31  ! :tex:  $\frac{\partial^2 \chi}{\partial r \partial \phi}$
    ! Integer, Parameter :: chi_p_d2rp   = chi_off + 32  ! :tex:  $\frac{\partial^2 \chi'}{\partial r \partial \phi}$
    ! Integer, Parameter :: chi_m_d2rp   = chi_off + 33  ! :tex:  $\frac{\partial^2 \overline{\chi}}{\partial r \partial \phi}$
    ! !------------ Theta-Phi Derivatives --------------!
    ! Integer, Parameter :: chi_d2tp     = chi_off + 34  ! :tex:  $\frac{\partial^2 \chi}{\partial \theta \partial \phi}$
    ! Integer, Parameter :: chi_p_d2tp   = chi_off + 35  ! :tex:  $\frac{\partial^2 \chi'}{\partial \theta \partial \phi}$
    ! Integer, Parameter :: chi_m_d2tp   = chi_off + 36  ! :tex:  $\frac{\partial^2 \overline{\chi}}{\partial \theta \partial \phi}$


    ! FIXME: not implemented, placeholder - NOTE definitions are not correct/finalized
    Integer, Parameter :: chim_off = chi_off + 50 ! :OFFSET CODE:
    ! !---  Masses
    ! Integer, Parameter :: compositional_mass_full = chim_off + 1  ! rho_bar chi  :tex: $\mathrm{f}_1\mathrm{f}_4\Theta$
    ! Integer, Parameter :: compositional_mass_p    = chim_off + 2  ! rho_bar chi' :tex: $\mathrm{f}_1\mathrm{f}_4\Theta$
    ! Integer, Parameter :: compositional_mass_m    = chim_off + 3  ! rho_bar chi_bar <S> :tex: $\mathrm{f}_1\mathrm{f}_4\overline{\Theta}$

    ! !--- Quadratics
    ! Integer, Parameter :: compositional_mass_sq = chim_off + 4  ! rho_bar T_bar S  :tex: $\left(\mathrm{f}_1\mathrm{f}_4\Theta\right)^2$
    ! Integer, Parameter :: compositional_massp_sq    = chim_off + 5  ! rho_bar T_bar S' :tex: $\left(\mathrm{f}_1\mathrm{f}_4\Theta\right)^2$
    ! Integer, Parameter :: compositional_massm_sq    = chim_off + 6  ! rho_bar T_bar <S> :tex: $\left(\mathrm{f}_1\mathrm{f}_4\overline{\Theta}\right)^2$

    ! !--- Momentum equation Chi contributions
    Integer, Parameter :: chi_mom_eq_off = chim_off + 20 ! :OFFSET CODE:
    Integer, Parameter :: chi_buoyancy_force  =  chi_mom_eq_off+1  ! :tex: $d_{7_i}\mathrm{g}_{7_i}\chi$
    Integer, Parameter :: chi_buoyancy_pforce =  chi_mom_eq_off+2  ! :tex: $d_{7_i}\mathrm{g}_{7_i}\chi'$
    Integer, Parameter :: chi_buoyancy_mforce =  chi_mom_eq_off+3  ! :tex: $d_{7_i}\mathrm{g}_{7_i}\overline{\chi}$

    ! FIXME: not implemented, placeholder
    ! ! ell=0 pressure and buoyancy forces (r-direction only)
    ! ! These are substracted out from the radial terms above
    ! ! because they do not contribute to the dynamics.
    ! ! The ell=0 chi buoyancy force is determined by the mean chi profile.
    ! ! The ell=0 pressure force is DEFINED such that it exactly balances the ell=0 buoyancy
    Integer, Parameter :: chi_buoyancy_force_ell0 = chi_mom_eq_off+4    ! :tex: $d_{7_i}\mathrm{g}_{7_i}\chi_{00}$


    ! FIXME: not implemented, placeholder - NOTE definitions are not correct/finalized
    ! !--- Chi equations contributions
    Integer, Parameter :: chi_eq_off = chi_mom_eq_off + 10 ! :OFFSET CODE:
    ! ! Advection (Full)
    ! Integer, Parameter :: rhotv_grad_chi   = chi_eq_off+1 ! :tex: $\mathrm{f}_1\mathrm{f}_4\boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta$
    ! Integer, Parameter :: rhotvp_grad_chip = chi_eq_off+2 ! :tex: $\mathrm{f}_1\mathrm{f}_4\boldsymbol{v'}\cdot\boldsymbol{\nabla}\Theta'$
    ! Integer, Parameter :: rhotvp_grad_chim = chi_eq_off+3 ! :tex: $\mathrm{f}_1\mathrm{f}_4\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\Theta}$
    ! Integer, Parameter :: rhotvm_grad_chim = chi_eq_off+4 ! :tex: $\mathrm{f}_1\mathrm{f}_4\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\Theta}$
    ! Integer, Parameter :: rhotvm_grad_chip = chi_eq_off+5 ! :tex: $\mathrm{f}_1\mathrm{f}_4\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\Theta'$
    ! ! Advection (Pieces)
    ! Integer, Parameter :: rhotvr_grad_chi   = chi_eq_off+6 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r\frac{\partial s}{\partial r}$
    ! Integer, Parameter :: rhotvpr_grad_chip = chi_eq_off+7 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r'\frac{\partial \Theta'}{\partial r}$
    ! Integer, Parameter :: rhotvpr_grad_chim = chi_eq_off+8 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r'\frac{\partial \overline{\Theta}}{\partial r}$
    ! Integer, Parameter :: rhotvmr_grad_chim = chi_eq_off+9 ! :tex: $\mathrm{f}_1\mathrm{f}_4\overline{v_r}\frac{\partial \overline{\Theta}}{\partial r}$
    ! Integer, Parameter :: rhotvmr_grad_chip = chi_eq_off+10 ! :tex: $\mathrm{f}_1\mathrm{f}_4\overline{v_r}\frac{\partial \Theta'}{\partial r}$

    ! Integer, Parameter :: rhotvt_grad_chi   = chi_eq_off+11 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\theta}{r}\frac{\partial \Theta}{\partial \theta}$
    ! Integer, Parameter :: rhotvpt_grad_chip = chi_eq_off+12 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\theta'}{r}\frac{\partial \Theta'}{\partial \theta}$
    ! Integer, Parameter :: rhotvpt_grad_chim = chi_eq_off+13 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\theta'}{r}\frac{\partial \overline{\Theta}}{\partial \theta}$
    ! Integer, Parameter :: rhotvmt_grad_chim = chi_eq_off+14 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{\overline{v_\theta}}{r}\frac{\partial \overline{\Theta}}{\partial \theta}$
    ! Integer, Parameter :: rhotvmt_grad_chip = chi_eq_off+15 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{\overline{v_\theta}}{r}\frac{\partial \Theta'}{\partial \theta}$

    ! Integer, Parameter :: rhotvp_grad_chi   = chi_eq_off+16 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\phi}{r \mathrm{sin}\theta}\frac{\partial \Theta}{\partial \phi}$
    ! Integer, Parameter :: rhotvpp_grad_chip = chi_eq_off+17 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\phi'}{r \mathrm{sin}\theta}\frac{\partial \Theta'}{\partial \phi}$
    ! Integer, Parameter :: rhotvpp_grad_chim = chi_eq_off+18 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\phi'}{r \mathrm{sin}\theta}\frac{\partial \overline{\Theta}}{\partial \phi}$
    ! Integer, Parameter :: rhotvmp_grad_chim = chi_eq_off+19 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{\overline{v_\phi}}{r \mathrm{sin}\theta}\frac{\partial \overline{\Theta}}{\partial \phi}$
    ! Integer, Parameter :: rhotvmp_grad_chip = chi_eq_off+20 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{\overline{v_\phi}}{r \mathrm{sin}\theta}\frac{\partial \Theta'}{\partial \phi}$

    ! ! Diffusion
    ! Integer, Parameter :: chi_diff  = chi_eq_off+21 !  :tex:  $-c_6\boldsymbol{\nabla}\cdot\boldsymbol{F}_{cond}$
    ! Integer, Parameter :: chip_diff = chi_eq_off+22 !  :tex:  $-c_6\boldsymbol{\nabla}\cdot\boldsymbol{F'}_{cond}$
    ! Integer, Parameter :: chim_diff = chi_eq_off+23 !  :tex:  $-c_6\boldsymbol{\nabla}\cdot\boldsymbol{\overline{F}}_{cond}$

    ! Integer, Parameter :: chi_diff_r = chi_eq_off+24  !  :tex:  $c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\left(\frac{\partial^2 \Theta}{\partial r^2}+\frac{\partial \Theta}{\partial r}\left[\frac{2}{r}+\frac{\mathrm{d}}{\mathrm{dr}}\mathrm{ln}\left\{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\right\} \right]\right)$
    ! Integer, Parameter :: chip_diff_r = chi_eq_off+25 !  :tex:  $c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\left(\frac{\partial^2 \Theta'}{\partial r^2}+\frac{\partial \Theta'}{\partial r}\left[\frac{2}{r}+\frac{\mathrm{d}}{\mathrm{dr}}\mathrm{ln}\left\{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\right\} \right]\right)$
    ! Integer, Parameter :: chim_diff_r = chi_eq_off+26 !  :tex:  $c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\left(\frac{\partial^2 \overline{\Theta}}{\partial r^2}+\frac{\partial \overline{\Theta}}{\partial r}\left[\frac{2}{r}+\frac{\mathrm{d}}{\mathrm{dr}}\mathrm{ln}\left\{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\right\} \right]\right)$

    ! Integer, Parameter :: chi_diff_theta = chi_eq_off+27  !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2}\left(\frac{\partial^2 \Theta}{\partial \theta^2}+\mathrm{cot}\theta\frac{\partial s}{\partial\theta}\right)$
    ! Integer, Parameter :: chip_diff_theta = chi_eq_off+28 !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2}\left(\frac{\partial^2 \Theta'}{\partial \theta^2}+\mathrm{cot}\theta\frac{\partial \Theta'}{\partial\theta}\right)$
    ! Integer, Parameter :: chim_diff_theta = chi_eq_off+29 !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2}\left(\frac{\partial^2 \overline{\Theta}}{\partial \theta^2}+\mathrm{cot}\theta\frac{\partial \overline{\Theta}}{\partial\theta}\right)$

    ! Integer, Parameter :: chi_diff_phi = chi_eq_off+30  !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2\mathrm{sin^2}\theta}\frac{\partial^2 \Theta}{\partial \phi^2}$
    ! Integer, Parameter :: chip_diff_phi = chi_eq_off+31 !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2\mathrm{sin^2}\theta}\frac{\partial^2 \Theta'}{\partial \phi^2}$
    ! Integer, Parameter :: chim_diff_phi = chi_eq_off+32 !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2\mathrm{sin^2}\theta}\frac{\partial^2 \overline{\Theta}}{\partial \phi^2}$

    ! ! Mass sources...
    ! ! Q
    ! Integer, Parameter :: chi_vol_mass_flux = chi_eq_off+33 ! :tex: $F_Q(r)$
    ! Integer, Parameter :: chi_vol_source   = chi_eq_off+34 ! :tex: $\mathrm{f}_6(r)$


    ! ! Advective Fluxes (mass)
    ! Integer, Parameter :: rhot_vr_chi = chi_eq_off+40 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r\Theta$
    ! Integer, Parameter :: rhot_vrp_chip = chi_eq_off+41 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r'\Theta'$
    ! Integer, Parameter :: rhot_vrp_chim = chi_eq_off+42 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r'\overline{\Theta}$
    ! Integer, Parameter :: rhot_vrm_chip = chi_eq_off+43 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_r}\Theta'$
    ! Integer, Parameter :: rhot_vrm_chim = chi_eq_off+44 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_r}\,\overline{\Theta}$

    ! Integer, Parameter :: rhot_vt_chi   = chi_eq_off+45 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\theta \Theta$
    ! Integer, Parameter :: rhot_vtp_chip = chi_eq_off+46 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\theta'\Theta'$
    ! Integer, Parameter :: rhot_vtp_chim = chi_eq_off+47 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\theta'\overline{\Theta}$
    ! Integer, Parameter :: rhot_vtm_chip = chi_eq_off+48 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_\theta}\Theta'$
    ! Integer, Parameter :: rhot_vtm_chim = chi_eq_off+49 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_\theta}\,\overline{\Theta}$

    ! Integer, Parameter :: rhot_vp_chi   = chi_eq_off+50 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\phi \Theta$
    ! Integer, Parameter :: rhot_vpp_chip = chi_eq_off+51 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\phi'\Theta'$
    ! Integer, Parameter :: rhot_vpp_chim = chi_eq_off+52 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\phi'\overline{\Theta}$
    ! Integer, Parameter :: rhot_vpm_chip = chi_eq_off+53 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_\phi}\Theta'$
    ! Integer, Parameter :: rhot_vpm_chim = chi_eq_off+54 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_\phi}\,\overline{\Theta}$

    ! ! Diffusive Fluxes
    ! Integer, Parameter :: chi_cond_flux_r     = chi_eq_off+70    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{\partial \Theta}{\partial r}$
    ! Integer, Parameter :: chi_cond_flux_theta = chi_eq_off+71    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r}\frac{\partial \Theta}{\partial \theta}$
    ! Integer, Parameter :: chi_cond_flux_phi   = chi_eq_off+72    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r \mathrm{sin}\theta}\frac{\partial \Theta}{\partial \phi}$

    ! Integer, Parameter :: chi_cond_fluxp_r     = chi_eq_off+73    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{\partial \Theta'}{\partial r}$
    ! Integer, Parameter :: chicond_fluxp_theta = chi_eq_off+74    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r}\frac{\partial \Theta'}{\partial \theta}$
    ! Integer, Parameter :: chi_cond_fluxp_phi   = chi_eq_off+75    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r \mathrm{sin}\theta}\frac{\partial \Theta'}{\partial \phi}$

    ! Integer, Parameter :: chi_cond_fluxm_r     = chi_eq_off+76    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{\partial \overline{\Theta}}{\partial r}$
    ! Integer, Parameter :: chi_cond_fluxm_theta = chi_eq_off+77    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r}\frac{\partial \overline{\Theta}}{\partial \theta}$
    ! Integer, Parameter :: chi_cond_fluxm_phi   = chi_eq_off+78    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r \mathrm{sin}\theta}\frac{\partial \overline{\Theta}}{\partial \phi}$

    ! Kinetic energy chi contributions
    Integer, Parameter :: chi_keq_off = chi_eq_off + 90 ! :OFFSET CODE:

    Integer, Parameter :: chi_buoy_work    = chi_keq_off+1  ! :tex: $d_{7_i}v_r\mathrm{g}_{7_i}\chi$
    Integer, Parameter :: chi_buoy_work_pp = chi_keq_off+2  ! :tex: $d_{7_i}v'_r\mathrm{g}_{7_i}\chi'$
    Integer, Parameter :: chi_buoy_work_mm = chi_keq_off+3  ! :tex: $d_{7_i}\overline{v_r}\mathrm{g}_{7_i}\overline{\chi}$

    ! Turbulent kinetic energy chi contributions
    Integer, Parameter :: chi_turbke_offset = chi_keq_off + 10 ! :OFFSET CODE:
    Integer, Parameter :: chi_production_buoyant_pKE   = chi_turbke_offset + 1  ! :tex: $\mathrm{g}_{7_i}{\chi^\prime}{v^\prime}_r$