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

    Integer, Parameter :: te_eq_off = mom_eq_off+200 ! :OFFSET CODE:
    ! Advection (Full)
    Integer, Parameter :: rhotv_grad_s   = te_eq_off+1 ! :tex: $\mathrm{f}_1\mathrm{f}_4\boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta$
    Integer, Parameter :: rhotvp_grad_sp = te_eq_off+2 ! :tex: $\mathrm{f}_1\mathrm{f}_4\boldsymbol{v'}\cdot\boldsymbol{\nabla}\Theta'$
    Integer, Parameter :: rhotvp_grad_sm = te_eq_off+3 ! :tex: $\mathrm{f}_1\mathrm{f}_4\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\Theta}$
    Integer, Parameter :: rhotvm_grad_sm = te_eq_off+4 ! :tex: $\mathrm{f}_1\mathrm{f}_4\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\Theta}$
    Integer, Parameter :: rhotvm_grad_sp = te_eq_off+5 ! :tex: $\mathrm{f}_1\mathrm{f}_4\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\Theta'$
    ! Advection (Pieces)
    Integer, Parameter :: rhotvr_grad_s   = te_eq_off+6 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r\frac{\partial s}{\partial r}$
    Integer, Parameter :: rhotvpr_grad_sp = te_eq_off+7 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r'\frac{\partial \Theta'}{\partial r}$
    Integer, Parameter :: rhotvpr_grad_sm = te_eq_off+8 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r'\frac{\partial \overline{\Theta}}{\partial r}$
    Integer, Parameter :: rhotvmr_grad_sm = te_eq_off+9 ! :tex: $\mathrm{f}_1\mathrm{f}_4\overline{v_r}\frac{\partial \overline{\Theta}}{\partial r}$
    Integer, Parameter :: rhotvmr_grad_sp = te_eq_off+10 ! :tex: $\mathrm{f}_1\mathrm{f}_4\overline{v_r}\frac{\partial \Theta'}{\partial r}$

    Integer, Parameter :: rhotvt_grad_s   = te_eq_off+11 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\theta}{r}\frac{\partial \Theta}{\partial \theta}$
    Integer, Parameter :: rhotvpt_grad_sp = te_eq_off+12 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\theta'}{r}\frac{\partial \Theta'}{\partial \theta}$
    Integer, Parameter :: rhotvpt_grad_sm = te_eq_off+13 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\theta'}{r}\frac{\partial \overline{\Theta}}{\partial \theta}$
    Integer, Parameter :: rhotvmt_grad_sm = te_eq_off+14 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{\overline{v_\theta}}{r}\frac{\partial \overline{\Theta}}{\partial \theta}$
    Integer, Parameter :: rhotvmt_grad_sp = te_eq_off+15 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{\overline{v_\theta}}{r}\frac{\partial \Theta'}{\partial \theta}$

    Integer, Parameter :: rhotvp_grad_s   = te_eq_off+16 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\phi}{r \mathrm{sin}\theta}\frac{\partial \Theta}{\partial \phi}$
    Integer, Parameter :: rhotvpp_grad_sp = te_eq_off+17 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\phi'}{r \mathrm{sin}\theta}\frac{\partial \Theta'}{\partial \phi}$
    Integer, Parameter :: rhotvpp_grad_sm = te_eq_off+18 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{v_\phi'}{r \mathrm{sin}\theta}\frac{\partial \overline{\Theta}}{\partial \phi}$
    Integer, Parameter :: rhotvmp_grad_sm = te_eq_off+19 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{\overline{v_\phi}}{r \mathrm{sin}\theta}\frac{\partial \overline{\Theta}}{\partial \phi}$
    Integer, Parameter :: rhotvmp_grad_sp = te_eq_off+20 ! :tex: $\mathrm{f}_1\mathrm{f}_4\frac{\overline{v_\phi}}{r \mathrm{sin}\theta}\frac{\partial \Theta'}{\partial \phi}$

    ! Diffusion

    Integer, Parameter :: s_diff  = te_eq_off+21 !  :tex:  $-c_6\boldsymbol{\nabla}\cdot\boldsymbol{F}_{cond}$
    Integer, Parameter :: sp_diff = te_eq_off+22 !  :tex:  $-c_6\boldsymbol{\nabla}\cdot\boldsymbol{F'}_{cond}$
    Integer, Parameter :: sm_diff = te_eq_off+23 !  :tex:  $-c_6\boldsymbol{\nabla}\cdot\boldsymbol{\overline{F}}_{cond}$

    Integer, Parameter :: s_diff_r = te_eq_off+24  !  :tex:  $c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\left(\frac{\partial^2 \Theta}{\partial r^2}+\frac{\partial \Theta}{\partial r}\left[\frac{2}{r}+\frac{\mathrm{d}}{\mathrm{dr}}\mathrm{ln}\left\{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\right\} \right]\right)$
    Integer, Parameter :: sp_diff_r = te_eq_off+25 !  :tex:  $c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\left(\frac{\partial^2 \Theta'}{\partial r^2}+\frac{\partial \Theta'}{\partial r}\left[\frac{2}{r}+\frac{\mathrm{d}}{\mathrm{dr}}\mathrm{ln}\left\{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\right\} \right]\right)$
    Integer, Parameter :: sm_diff_r = te_eq_off+26 !  :tex:  $c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\left(\frac{\partial^2 \overline{\Theta}}{\partial r^2}+\frac{\partial \overline{\Theta}}{\partial r}\left[\frac{2}{r}+\frac{\mathrm{d}}{\mathrm{dr}}\mathrm{ln}\left\{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\right\} \right]\right)$

    Integer, Parameter :: s_diff_theta = te_eq_off+27  !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2}\left(\frac{\partial^2 \Theta}{\partial \theta^2}+\mathrm{cot}\theta\frac{\partial s}{\partial\theta}\right)$
    Integer, Parameter :: sp_diff_theta = te_eq_off+28 !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2}\left(\frac{\partial^2 \Theta'}{\partial \theta^2}+\mathrm{cot}\theta\frac{\partial \Theta'}{\partial\theta}\right)$
    Integer, Parameter :: sm_diff_theta = te_eq_off+29 !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2}\left(\frac{\partial^2 \overline{\Theta}}{\partial \theta^2}+\mathrm{cot}\theta\frac{\partial \overline{\Theta}}{\partial\theta}\right)$

    Integer, Parameter :: s_diff_phi = te_eq_off+30  !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2\mathrm{sin^2}\theta}\frac{\partial^2 \Theta}{\partial \phi^2}$
    Integer, Parameter :: sp_diff_phi = te_eq_off+31 !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2\mathrm{sin^2}\theta}\frac{\partial^2 \Theta'}{\partial \phi^2}$
    Integer, Parameter :: sm_diff_phi = te_eq_off+32 !  :tex:  $c_6\frac{\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5}{r^2\mathrm{sin^2}\theta}\frac{\partial^2 \overline{\Theta}}{\partial \phi^2}$

    ! Heat sources...
    ! Q
    Integer, Parameter :: vol_heat_flux = te_eq_off+33 ! :tex: $F_Q(r)$
    Integer, Parameter :: vol_heating   = te_eq_off+34 ! :tex: $\mathrm{f}_6(r)$


    ! Viscous Heating
    Integer, Parameter :: visc_heating = te_eq_off+35 ! :tex:  $c_5\Phi(r,\theta,\phi)$

    ! Ohmic Heating

    Integer, Parameter :: ohmic_heat    = te_eq_off+36 !  :tex:  $ \mathrm{f}_7c_4\left( {\boldsymbol{\mathcal{J}'}}\cdot{\boldsymbol{\mathcal{J}'}} \right)$
    Integer, Parameter :: ohmic_heat_pp = te_eq_off+37 !  :tex: $ \mathrm{f}_7c_4\left( {\boldsymbol{\mathcal{J}'}}\cdot{\boldsymbol{\mathcal{J}'}} \right)$
    Integer, Parameter :: ohmic_heat_pm = te_eq_off+38 !  :tex: $\mathrm{f}_7c_4\left( \overline{\boldsymbol{\mathcal{J}}}\cdot\overline{\boldsymbol{\mathcal{J}}} \right)$
    Integer, Parameter :: ohmic_heat_mm = te_eq_off+39 !  :tex: $\mathrm{f}_7c_4\left( \overline{\boldsymbol{\mathcal{J}}}\cdot{\boldsymbol{\mathcal{J}'}} \right)$



    ! Advective Fluxes (heat)
    Integer, Parameter :: rhot_vr_s = te_eq_off+40 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r\Theta$
    Integer, Parameter :: rhot_vrp_sp = te_eq_off+41 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r'\Theta'$
    Integer, Parameter :: rhot_vrp_sm = te_eq_off+42 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_r'\overline{\Theta}$
    Integer, Parameter :: rhot_vrm_sp = te_eq_off+43 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_r}\Theta'$
    Integer, Parameter :: rhot_vrm_sm = te_eq_off+44 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_r}\,\overline{\Theta}$

    Integer, Parameter :: rhot_vt_s   = te_eq_off+45 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\theta \Theta$
    Integer, Parameter :: rhot_vtp_sp = te_eq_off+46 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\theta'\Theta'$
    Integer, Parameter :: rhot_vtp_sm = te_eq_off+47 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\theta'\overline{\Theta}$
    Integer, Parameter :: rhot_vtm_sp = te_eq_off+48 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_\theta}\Theta'$
    Integer, Parameter :: rhot_vtm_sm = te_eq_off+49 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_\theta}\,\overline{\Theta}$

    Integer, Parameter :: rhot_vp_s   = te_eq_off+50 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\phi \Theta$
    Integer, Parameter :: rhot_vpp_sp = te_eq_off+51 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\phi'\Theta'$
    Integer, Parameter :: rhot_vpp_sm = te_eq_off+52 ! :tex: $\mathrm{f}_1\mathrm{f}_4v_\phi'\overline{\Theta}$
    Integer, Parameter :: rhot_vpm_sp = te_eq_off+53 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_\phi}\Theta'$
    Integer, Parameter :: rhot_vpm_sm = te_eq_off+54 ! :tex: $\mathrm{f}_1\mathrm{f}_4\,\overline{v_\phi}\,\overline{\Theta}$


    ! Advective Fluxes (enthalpy -- thermo-mechanical)
    Integer, Parameter :: enth_flux_r     = te_eq_off+55    ! :tex: $c_P \mathrm{f}_1 v_r T$
    Integer, Parameter :: enth_flux_theta = te_eq_off+56    ! :tex: $c_P \mathrm{f}_1 v_\theta T$
    Integer, Parameter :: enth_flux_phi   = te_eq_off+57    ! :tex: $c_P \mathrm{f}_1 v_\phi T$

    Integer, Parameter :: enth_flux_rpp     = te_eq_off+58    ! :tex: $c_P \mathrm{f}_1 v_r' T'$
    Integer, Parameter :: enth_flux_thetapp = te_eq_off+59    ! :tex: $c_P \mathrm{f}_1 v_\theta' T'$
    Integer, Parameter :: enth_flux_phipp   = te_eq_off+60    ! :tex: $c_P \mathrm{f}_1 v_\phi' T'$

    Integer, Parameter :: enth_flux_rpm     = te_eq_off+61    ! :tex: $c_P \mathrm{f}_1 v_r'\, \overline{T}$
    Integer, Parameter :: enth_flux_thetapm = te_eq_off+62    ! :tex: $c_P \mathrm{f}_1 v_\theta'\, \overline{T}$
    Integer, Parameter :: enth_flux_phipm   = te_eq_off+63    ! :tex: $c_P \mathrm{f}_1 v_\phi'\, \overline{T}$

    Integer, Parameter :: enth_flux_rmp     = te_eq_off+64    ! :tex: $c_P \mathrm{f}_1\, \overline{v_r}\, T'$
    Integer, Parameter :: enth_flux_thetamp = te_eq_off+65    ! :tex: $c_P \mathrm{f}_1\, \overline{v_\theta}\, T'$
    Integer, Parameter :: enth_flux_phimp   = te_eq_off+66    ! :tex: $c_P \mathrm{f}_1\, \overline{v_\phi}\, T'$

    Integer, Parameter :: enth_flux_rmm     = te_eq_off+67    ! :tex: $c_P \mathrm{f}_1\, \overline{v_r}\, \overline{T}$
    Integer, Parameter :: enth_flux_thetamm = te_eq_off+68    ! :tex: $c_P \mathrm{f}_1\, \overline{v_\theta}\, \overline{T}$
    Integer, Parameter :: enth_flux_phimm   = te_eq_off+69    ! :tex: $c_P \mathrm{f}_1\, \overline{v_\phi}\, \overline{T}$

    ! Diffusive Fluxes
    Integer, Parameter :: cond_flux_r     = te_eq_off+70    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{\partial \Theta}{\partial r}$
    Integer, Parameter :: cond_flux_theta = te_eq_off+71    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r}\frac{\partial \Theta}{\partial \theta}$
    Integer, Parameter :: cond_flux_phi   = te_eq_off+72    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r \mathrm{sin}\theta}\frac{\partial \Theta}{\partial \phi}$

    Integer, Parameter :: cond_fluxp_r     = te_eq_off+73    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{\partial \Theta'}{\partial r}$
    Integer, Parameter :: cond_fluxp_theta = te_eq_off+74    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r}\frac{\partial \Theta'}{\partial \theta}$
    Integer, Parameter :: cond_fluxp_phi   = te_eq_off+75    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r \mathrm{sin}\theta}\frac{\partial \Theta'}{\partial \phi}$

    Integer, Parameter :: cond_fluxm_r     = te_eq_off+76    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{\partial \overline{\Theta}}{\partial r}$
    Integer, Parameter :: cond_fluxm_theta = te_eq_off+77    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r}\frac{\partial \overline{\Theta}}{\partial \theta}$
    Integer, Parameter :: cond_fluxm_phi   = te_eq_off+78    !  :tex:  $-c_6\mathrm{f}_1\mathrm{f}_4\mathrm{f}_5\frac{1}{r \mathrm{sin}\theta}\frac{\partial \overline{\Theta}}{\partial \phi}$
