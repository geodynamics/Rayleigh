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

    !////////////////////////  Advection Terms ////////////////////
    ! Reynolds decomposition about the azimuthal mean may also be output
    ! NOTE:  ADVECTION TERMS ARE SCALED BY DENSITY (so that they represent a force density)

    Integer, Parameter :: curl_mom_eq_off = mom_eq_off+100 ! :OFFSET CODE:
    Integer, Parameter :: curl_v_grad_v_r       = curl_mom_eq_off+1 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right]_r$
    Integer, Parameter :: curl_v_grad_v_theta   = curl_mom_eq_off+2 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right]_\theta$
    Integer, Parameter :: curl_v_grad_v_phi     = curl_mom_eq_off+3 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right]_\phi$

    Integer, Parameter :: curl_v_grad_v_r_squared       = curl_mom_eq_off+4 ! :tex: $\left(\left[\nabla\times\left(\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right)\right]_r\right)^2$
    Integer, Parameter :: curl_v_grad_v_theta_squared   = curl_mom_eq_off+5 ! :tex: $\left(\left[\nabla\times\left(\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right)\right]_\theta\right)^2$
    Integer, Parameter :: curl_v_grad_v_phi_squared     = curl_mom_eq_off+6 ! :tex: $\left(\left[\nabla\times\left(\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]\right)\right]_\phi\right)^2$

!    Integer, Parameter :: curl_vp_grad_vm_r     = curl_mom_eq_off+7 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_r$
!    Integer, Parameter :: curl_vp_grad_vm_theta = curl_mom_eq_off+8 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\theta$
!    Integer, Parameter :: curl_vp_grad_vm_phi   = curl_mom_eq_off+9 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\phi$

!    Integer, Parameter :: curl_vm_grad_vp_r     = curl_mom_eq_off+10 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_r$
!    Integer, Parameter :: curl_vm_grad_vp_theta = curl_mom_eq_off+11 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\theta$
!    Integer, Parameter :: curl_vm_grad_vp_phi   = curl_mom_eq_off+12 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\phi$

!    Integer, Parameter :: curl_vp_grad_vp_r     = curl_mom_eq_off+13 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_r$
!    Integer, Parameter :: curl_vp_grad_vp_theta = curl_mom_eq_off+14 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\theta$
!    Integer, Parameter :: curl_vp_grad_vp_phi   = curl_mom_eq_off+15 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]\right]_\phi$

!    Integer, Parameter :: curl_vm_grad_vm_r     = curl_mom_eq_off+16 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_r$
!    Integer, Parameter :: curl_vm_grad_vm_theta = curl_mom_eq_off+17 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\theta$
!    Integer, Parameter :: curl_vm_grad_vm_phi   = curl_mom_eq_off+18 ! :tex: $\left[\nabla\times\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]\right]_\phi$

    !/////////////////////////////////////////////////////////////
    !  Linear forces
    !  Note:  the pressure gradient diagnostic codes are  above

    Integer, Parameter :: curl_buoyancy_force_theta  =  curl_mom_eq_off+19  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\theta$
    Integer, Parameter :: curl_buoyancy_force_phi  =  curl_mom_eq_off+20  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\phi$

    Integer, Parameter :: curl_buoyancy_force_theta_squared  =  curl_mom_eq_off+21  ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\theta\right)^2$
    Integer, Parameter :: curl_buoyancy_force_phi_squared  =  curl_mom_eq_off+22  ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta\right)\right]_\phi\right)^2$

!    Integer, Parameter :: curl_buoyancy_pforce_theta  =  curl_mom_eq_off+23  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\overline{\Theta}\right)\right]_\theta$
!    Integer, Parameter :: curl_buoyancy_pforce_phi  =  curl_mom_eq_off+24  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\overline{\Theta}\right)\right]_\phi$

!    Integer, Parameter :: curl_buoyancy_mforce_theta  =  curl_mom_eq_off+25  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta'\right)\right]_\theta$
!    Integer, Parameter :: curl_buoyancy_mforce_phi  =  curl_mom_eq_off+26  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_2\mathrm{f}_2\Theta'\right)\right]_\phi$


    Integer, Parameter :: curl_coriolis_force_r      = curl_mom_eq_off+27 ! :tex: $\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_r$
    Integer, Parameter :: curl_coriolis_force_theta  = curl_mom_eq_off+28 ! :tex: $\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\theta$
    Integer, Parameter :: curl_coriolis_force_phi    = curl_mom_eq_off+29 ! :tex: $\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\phi$

    Integer, Parameter :: curl_coriolis_force_r_squared      = curl_mom_eq_off+30 ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_r\right)^2$
    Integer, Parameter :: curl_coriolis_force_theta_squared  = curl_mom_eq_off+31 ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\theta\right)^2$
    Integer, Parameter :: curl_coriolis_force_phi_squared    = curl_mom_eq_off+32 ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v})\right)\right]_\phi\right)^2$

!    Integer, Parameter :: curl_coriolis_pforce_r     = curl_mom_eq_off+33 ! :tex: $\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_r$
!    Integer, Parameter :: curl_coriolis_pforce_theta = curl_mom_eq_off+34 ! :tex: $\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_\theta$
!    Integer, Parameter :: curl_coriolis_pforce_phi   = curl_mom_eq_off+35 ! :tex: $\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\boldsymbol{v}')\right)\right]_\phi$

!    Integer, Parameter :: curl_coriolis_mforce_r     = curl_mom_eq_off+36 ! :tex: $\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_r$
!    Integer, Parameter :: curl_coriolis_mforce_theta = curl_mom_eq_off+37 ! :tex: $\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_\theta$
!    Integer, Parameter :: curl_coriolis_mforce_phi   = curl_mom_eq_off+38 ! :tex: $\left[\boldsymbol{\nabla}\times\left(-c_1\mathrm{f}_1(\hat{\boldsymbol{z}}\times\overline{\boldsymbol{\hat{v}}})\right)\right]_\phi$

    ! Viscous forces
    Integer, Parameter :: curl_viscous_force_r       = curl_mom_eq_off+39 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_r$
    Integer, Parameter :: curl_viscous_force_theta   = curl_mom_eq_off+40 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\theta$
    Integer, Parameter :: curl_viscous_force_phi     = curl_mom_eq_off+41 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\phi$

    Integer, Parameter :: curl_viscous_force_r_squared       = curl_mom_eq_off+42 ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_r\right)^2$
    Integer, Parameter :: curl_viscous_force_theta_squared   = curl_mom_eq_off+43 ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\theta\right)^2$
    Integer, Parameter :: curl_viscous_force_phi_squared     = curl_mom_eq_off+44 ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D})\right)\right]_\phi\right)^2$

    Integer, Parameter :: curl_viscous_pforce_r      = curl_mom_eq_off+45 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_r$
    Integer, Parameter :: curl_viscous_pforce_theta  = curl_mom_eq_off+46 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_\theta$
    Integer, Parameter :: curl_viscous_pforce_phi    = curl_mom_eq_off+47 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal D'})\right)\right]_\phi$

    Integer, Parameter :: curl_viscous_mforce_r      = curl_mom_eq_off+51 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_r$
    Integer, Parameter :: curl_viscous_mforce_theta  = curl_mom_eq_off+52 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_\theta$
    Integer, Parameter :: curl_viscous_mforce_phi    = curl_mom_eq_off+53 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_5\,(\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal D}})\right)\right]_\phi$

    ! Pressure forces
    Integer, Parameter :: curl_pressure_force_theta   = curl_mom_eq_off+58 ! :tex: $\left[\nabla\times\nabla P\right]_\theta$
    Integer, Parameter :: curl_pressure_force_phi     = curl_mom_eq_off+59 ! :tex:  $\left[\nabla\times\nabla P\right]_\phi$

    Integer, Parameter :: curl_pressure_force_theta_squared   = curl_mom_eq_off+61 ! :tex: $\left[\nabla\times\nabla P\right]_\theta^2$
    Integer, Parameter :: curl_pressure_force_phi_squared     = curl_mom_eq_off+62 ! :tex: $\left[\nabla\times\nabla P\right]_\phi^2$

!    Integer, Parameter :: curl_pressure_pforce_theta  = curl_mom_eq_off+64 ! :tex: $\left[\nabla\times\nabla P\right]_\theta$
!    Integer, Parameter :: curl_pressure_pforce_phi    = curl_mom_eq_off+65 ! :tex: $\left[\nabla\times\nabla P\right]_\phi$

!    Integer, Parameter :: curl_pressure_mforce_theta  = curl_mom_eq_off+67 ! :tex: $\left[\nabla\times\nabla P\right]_\theta$
!    Integer, Parameter :: curl_pressure_mforce_phi    = curl_mom_eq_off+68 ! :tex: $\left[\nabla\times\nabla P\right]_\phi$


    !/////////////////////////// Lorentz forces ///////////////////////////////
    !  ref%Lorentz_Coeff * del  (del x B) x B
    !  ref%Lorentz_Coeff = 1/4pi when dimensional, Pr/(Pr_m E) when nondimesional
    !  j (below) is shorthand for ref%Lorentz_Coeff*delxB  (not quite the current density)

    Integer, Parameter :: curl_j_cross_b_r     = curl_mom_eq_off+69  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_r$
    Integer, Parameter :: curl_j_cross_b_theta = curl_mom_eq_off+70  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\theta$
    Integer, Parameter :: curl_j_cross_b_phi   = curl_mom_eq_off+71  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\phi$

    Integer, Parameter :: curl_j_cross_b_r_squared     = curl_mom_eq_off+72  ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_r\right)^2$
    Integer, Parameter :: curl_j_cross_b_theta_squared = curl_mom_eq_off+73  ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\theta\right)^2$
    Integer, Parameter :: curl_j_cross_b_phi_squared   = curl_mom_eq_off+74  ! :tex: $\left(\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B)\times\boldsymbol B\right)\right)\right]_\phi\right)^2$


!    Integer, Parameter :: curl_jp_cross_bm_r       = curl_mom_eq_off+75  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_r$
!    Integer, Parameter :: curl_jp_cross_bm_theta   = curl_mom_eq_off+76  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_\theta$
!    Integer, Parameter :: curl_jp_cross_bm_phi     = curl_mom_eq_off+77  ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\overline{\boldsymbol B}\right)\right)\right]_\phi$

!    Integer, Parameter :: curl_jm_cross_bp_r     = curl_mom_eq_off+78  ! :tex: $\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_r$
!    Integer, Parameter :: curl_jm_cross_bp_theta = curl_mom_eq_off+79  ! :tex: $\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_\theta$
!    Integer, Parameter :: curl_jm_cross_bp_phi   = curl_mom_eq_off+80  ! :tex: $\left[\boldsymbol{\nabla}\times c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\boldsymbol B'\right)\right]_\phi$

!    Integer, Parameter :: curl_jm_cross_bm_r     = curl_mom_eq_off+81 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_r$
!    Integer, Parameter :: curl_jm_cross_bm_theta = curl_mom_eq_off+82 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_\theta$
!    Integer, Parameter :: curl_jm_cross_bm_phi   = curl_mom_eq_off+83 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\boldsymbol B')\times\boldsymbol B'\right)\right)\right]_\phi$

!    Integer, Parameter :: curl_jp_cross_bp_r     = curl_mom_eq_off+84 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_r$
!    Integer, Parameter :: curl_jp_cross_bp_theta = curl_mom_eq_off+85 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_\theta$
!    Integer, Parameter :: curl_jp_cross_bp_phi   = curl_mom_eq_off+86 ! :tex: $\left[\boldsymbol{\nabla}\times\left(c_4\left((\boldsymbol{\nabla}\times\overline{\boldsymbol B})\times\overline{\boldsymbol B}\right)\right)\right]_\phi$

    Integer, Parameter :: curl_buoyancy_force_abs = curl_mom_eq_off+87
    Integer, Parameter :: curl_coriolis_force_abs = curl_mom_eq_off+88
    Integer, Parameter :: curl_v_grad_v_abs       = curl_mom_eq_off+89
    Integer, Parameter :: curl_j_cross_b_abs      = curl_mom_eq_off+90
    Integer, Parameter :: curl_viscous_force_abs  = curl_mom_eq_off+91
