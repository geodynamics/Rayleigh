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

    Integer, Parameter :: mom_eq_off = meoffset+100 ! :OFFSET CODE:
    Integer, Parameter :: v_grad_v_r       = mom_eq_off+1 ! :tex: $\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_r$
    Integer, Parameter :: v_grad_v_theta   = mom_eq_off+2 ! :tex: $\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\theta$
    Integer, Parameter :: v_grad_v_phi     = mom_eq_off+3 ! :tex: $\mathrm{f}_1\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\phi$

    Integer, Parameter :: vp_grad_vm_r     = mom_eq_off+4 ! :tex: $\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r$
    Integer, Parameter :: vp_grad_vm_theta = mom_eq_off+5 ! :tex: $\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta$
    Integer, Parameter :: vp_grad_vm_phi   = mom_eq_off+6 ! :tex: $\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi$

    Integer, Parameter :: vm_grad_vp_r     = mom_eq_off+7 ! :tex: $\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r$
    Integer, Parameter :: vm_grad_vp_theta = mom_eq_off+8 ! :tex: $\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta$
    Integer, Parameter :: vm_grad_vp_phi   = mom_eq_off+9 ! :tex: $\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi$

    Integer, Parameter :: vp_grad_vp_r     = mom_eq_off+10 ! :tex: $\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r$
    Integer, Parameter :: vp_grad_vp_theta = mom_eq_off+11 ! :tex: $\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta$
    Integer, Parameter :: vp_grad_vp_phi   = mom_eq_off+12 ! :tex: $\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi$

    Integer, Parameter :: vm_grad_vm_r     = mom_eq_off+13 ! :tex: $\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r$
    Integer, Parameter :: vm_grad_vm_theta = mom_eq_off+14 ! :tex: $\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta$
    Integer, Parameter :: vm_grad_vm_phi   = mom_eq_off+15 ! :tex: $\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi$

    !/////////////////////////////////////////////////////////////
    !  Linear Forces
    !  Note:  the pressure gradient diagnostic codes are  above

    Integer, Parameter :: buoyancy_force  =  mom_eq_off+16  ! :tex: $c_2\mathrm{f}_2\Theta$
    Integer, Parameter :: buoyancy_pforce =  mom_eq_off+17  ! :tex: $c_2\mathrm{f}_2\Theta'$
    Integer, Parameter :: buoyancy_mforce =  mom_eq_off+18  ! :tex: $c_2\mathrm{f}_2\overline{\Theta}$

    Integer, Parameter :: Coriolis_Force_r      = mom_eq_off+19 ! :tex: $-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v}\right]_r$
    Integer, Parameter :: Coriolis_Force_theta  = mom_eq_off+20 ! :tex: $-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v}\right]_\theta$
    Integer, Parameter :: Coriolis_Force_phi    = mom_eq_off+21 ! :tex: $-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v}\right]_\phi$

    Integer, Parameter :: Coriolis_pForce_r     = mom_eq_off+22 ! :tex: $-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v'}\right]_r$
    Integer, Parameter :: Coriolis_pForce_theta = mom_eq_off+23 ! :tex: $-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v'}\right]_\theta$
    Integer, Parameter :: Coriolis_pForce_phi   = mom_eq_off+24 ! :tex: $-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\boldsymbol{v'}\right]_\phi$

    Integer, Parameter :: Coriolis_mForce_r     = mom_eq_off+25 ! :tex: $-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\overline{\boldsymbol{v}}\right]_r$
    Integer, Parameter :: Coriolis_mForce_theta = mom_eq_off+26 ! :tex: $-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\overline{\boldsymbol{v}}\right]_\theta$
    Integer, Parameter :: Coriolis_mForce_phi   = mom_eq_off+27 ! :tex: $-c_1\mathrm{f}_1\left[\boldsymbol{\hat{z}}\times\overline{\boldsymbol{v}}\right]_\phi$

    ! Viscous forces
    Integer, Parameter :: viscous_Force_r       = mom_eq_off+28 ! :tex: $c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\right]_r$
    Integer, Parameter :: viscous_Force_theta   = mom_eq_off+29 ! :tex: $c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\right]_\theta$
    Integer, Parameter :: viscous_Force_phi     = mom_eq_off+30 ! :tex: $c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\right]_\phi$

    Integer, Parameter :: viscous_pForce_r      = mom_eq_off+31 ! :tex: $c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D'}}\right]_r$
    Integer, Parameter :: viscous_pForce_theta  = mom_eq_off+32 ! :tex: $c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D'}}\right]_\theta$
    Integer, Parameter :: viscous_pForce_phi    = mom_eq_off+33 ! :tex: $c_5\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D'}}\right]_\phi$

    Integer, Parameter :: viscous_mForce_r      = mom_eq_off+34 ! :tex: $c_5\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]_r$
    Integer, Parameter :: viscous_mForce_theta  = mom_eq_off+35 ! :tex: $c_5\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]_\theta$
    Integer, Parameter :: viscous_mForce_phi    = mom_eq_off+36 ! :tex: $c_5\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]_\phi$

    ! Pressure forces
    Integer, Parameter :: pressure_Force_r       = mom_eq_off+37 ! :tex: $ -c_3\mathrm{f}_1\frac{\partial}{\partial r}\left(\frac{P}{\mathrm{f}_1} \right)$
    Integer, Parameter :: pressure_Force_theta   = mom_eq_off+38 ! :tex: $ -c_3\frac{1}{r}\frac{\partial P}{\partial \theta}$
    Integer, Parameter :: pressure_Force_phi     = mom_eq_off+39 ! :tex: $ -c_3\frac{1}{r\mathrm{sin}\theta}\frac{\partial P}{\partial \phi}$

    Integer, Parameter :: pressure_pForce_r      = mom_eq_off+40 ! :tex: $ -c_3\mathrm{f}_1\frac{\partial}{\partial r}\left(\frac{P'}{\mathrm{f}_1} \right)$
    Integer, Parameter :: pressure_pForce_theta  = mom_eq_off+41 ! :tex: $ -c_3\frac{1}{r}\frac{\partial P'}{\partial \theta}$
    Integer, Parameter :: pressure_pForce_phi    = mom_eq_off+42 ! :tex: $ -c_3\frac{1}{r\mathrm{sin}\theta}\frac{\partial P'}{\partial \phi}$

    Integer, Parameter :: pressure_mForce_r      = mom_eq_off+43 ! :tex: $ -c_3\mathrm{f}_1\frac{\partial}{\partial r}\left(\frac{\overline{P}}{\mathrm{f}_1} \right)$
    Integer, Parameter :: pressure_mForce_theta  = mom_eq_off+44 ! :tex: $ -c_3\frac{1}{r}\frac{\partial \overline{P}}{\partial \theta}$
    Integer, Parameter :: pressure_mForce_phi    = mom_eq_off+45 ! :tex: $ -c_3\frac{1}{r\mathrm{sin}\theta}\frac{\partial \overline{P}}{\partial \phi}$

    ! ell=0 pressure and buoyancy forces (r-direction only)
    ! These are substracted out from the radial terms above
    ! because they do not contribute to the dynamics.
    ! The ell=0 buoyancy force is determined by the mean entropy profile.
    ! The ell=0 pressure force is DEFINED such that it exactly balances the ell=0 buoyancy
    Integer, Parameter :: buoyancy_force_ell0 = mom_eq_off+46    ! :tex: $c_2\mathrm{f}_2\Theta_{00}$
    Integer, Parameter :: pressure_force_ell0_r = mom_eq_off+47  ! :tex: $-c_3\mathrm{f}_1\frac{\partial}{\partial r}\left(\frac{P_{00}}{\mathrm{f}_1} \right)$

    !/////////////////////////// Lorentz Forces ///////////////////////////////
    !  ref%Lorentz_Coeff * (del x B) x B
    !  ref%Lorentz_Coeff = 1/4pi when dimensional, Pr/(Pr_m E) when nondimesional
    !  j (below) is shorthand for ref%Lorentz_Coeff*delxB  (not quite the current density)

    Integer, Parameter :: j_cross_b_r       = mom_eq_off+48  ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B}\right]_r$
    Integer, Parameter :: j_cross_b_theta   = mom_eq_off+49  ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B}\right]_\theta$
    Integer, Parameter :: j_cross_b_phi     = mom_eq_off+50  ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B}\right]_\phi$

    Integer, Parameter :: jp_cross_bm_r     = mom_eq_off+51  ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\overline{\boldsymbol{B}}\right]_r$
    Integer, Parameter :: jp_cross_bm_theta = mom_eq_off+52  ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\overline{\boldsymbol{B}}\right]_\theta$
    Integer, Parameter :: jp_cross_bm_phi   = mom_eq_off+53  ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\overline{\boldsymbol{B}}\right]_\phi$

    Integer, Parameter :: jm_cross_bp_r     = mom_eq_off+54  ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\boldsymbol{B'}\right]_r$
    Integer, Parameter :: jm_cross_bp_theta = mom_eq_off+55  ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\boldsymbol{B'}\right]_\theta$
    Integer, Parameter :: jm_cross_bp_phi   = mom_eq_off+56  ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\boldsymbol{B'}\right]_\phi$

    Integer, Parameter :: jm_cross_bm_r     = mom_eq_off+57 ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right]_r$
    Integer, Parameter :: jm_cross_bm_theta = mom_eq_off+58 ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right]_\theta$
    Integer, Parameter :: jm_cross_bm_phi   = mom_eq_off+59 ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right]_\phi$

    Integer, Parameter :: jp_cross_bp_r     = mom_eq_off+60 ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_r$
    Integer, Parameter :: jp_cross_bp_theta = mom_eq_off+61 ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_\theta$
    Integer, Parameter :: jp_cross_bp_phi   = mom_eq_off+62 ! :tex: $c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_\phi$
