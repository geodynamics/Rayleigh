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

    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Induction equation terms
    Integer, Parameter :: indoff = te_eq_off+200 ! :OFFSET CODE:

    !--------------- Terms involving v x B (full)
    Integer, Parameter :: induct_shear_r          = indoff+1  ! :tex: $\left[\boldsymbol{B}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_r$
    Integer, Parameter :: induct_comp_r           = indoff+2  ! :tex: $-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v} \right)B_r$
    Integer, Parameter :: induct_advec_r          = indoff+3  ! :tex: $-\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{B}\right]_r$
    Integer, Parameter :: induct_r                = indoff+4  ! :tex: $\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v}\times\boldsymbol{B}\right)\right]_r$
    Integer, Parameter :: induct_diff_r           = indoff+5  ! :tex: $-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B}\right)\right]_r$

    Integer, Parameter :: induct_shear_theta      = indoff+6  ! :tex: $\left[\boldsymbol{B}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\theta$
    Integer, Parameter :: induct_comp_theta       = indoff+7  ! :tex: $-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v} \right)B_\theta$
    Integer, Parameter :: induct_advec_theta      = indoff+8  ! :tex: $-\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{B}\right]_\theta$
    Integer, Parameter :: induct_theta            = indoff+9  ! :tex: $\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v}\times\boldsymbol{B}\right)\right]_\theta$
    Integer, Parameter :: induct_diff_theta       = indoff+10 ! :tex: $-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B}\right)\right]_\theta$

    Integer, Parameter :: induct_shear_phi        = indoff+11 ! :tex: $\left[\boldsymbol{B}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]_\phi$
    Integer, Parameter :: induct_comp_phi         = indoff+12 ! :tex: $-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v} \right)B_\phi$
    Integer, Parameter :: induct_advec_phi        = indoff+13 ! :tex: $-\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{B}\right]_\phi$
    Integer, Parameter :: induct_phi              = indoff+14 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v}\times\boldsymbol{B}\right)\right]_\phi$
    Integer, Parameter :: induct_diff_phi         = indoff+15 ! :tex: $-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B}\right)\right]_\phi$

    !--------------- Terms involving <v> x <B>
    Integer, Parameter :: induct_shear_vmbm_r     = indoff+16 ! :tex: $\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r$
    Integer, Parameter :: induct_comp_vmbm_r      = indoff+17 ! :tex: $-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)\overline{B_r}$
    Integer, Parameter :: induct_advec_vmbm_r     = indoff+18 ! :tex: $-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_r$
    Integer, Parameter :: induct_vmbm_r           = indoff+19 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\overline{\boldsymbol{B}}\right)\right]_r$
    Integer, Parameter :: induct_diff_bm_r        = indoff+20 ! :tex: $-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\right]_r$

    Integer, Parameter :: induct_shear_vmbm_theta = indoff+21 ! :tex: $\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta$
    Integer, Parameter :: induct_comp_vmbm_theta  = indoff+22 ! :tex: $-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)\overline{B_\theta}$
    Integer, Parameter :: induct_advec_vmbm_theta = indoff+23 ! :tex: $-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_\theta$
    Integer, Parameter :: induct_vmbm_theta       = indoff+24 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\overline{\boldsymbol{B}}\right)\right]_\theta$
    Integer, Parameter :: induct_diff_bm_theta    = indoff+25 ! :tex: $-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\right]_\theta$

    Integer, Parameter :: induct_shear_vmbm_phi   = indoff+26 ! :tex: $\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi$
    Integer, Parameter :: induct_comp_vmbm_phi    = indoff+27 ! :tex: $-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)\overline{B_\phi}$
    Integer, Parameter :: induct_advec_vmbm_phi   = indoff+28 ! :tex: $-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_\phi$
    Integer, Parameter :: induct_vmbm_phi         = indoff+29 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\overline{\boldsymbol{B}}\right)\right]_\phi$
    Integer, Parameter :: induct_diff_bm_phi      = indoff+30 ! :tex: $-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\right]_\phi$


    !--------------- Terms involving <v> x B'
    Integer, Parameter :: induct_shear_vmbp_r     = indoff+31 ! :tex: $\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_r$
    Integer, Parameter :: induct_comp_vmbp_r      = indoff+32 ! :tex: $-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)B'_r$
    Integer, Parameter :: induct_advec_vmbp_r     = indoff+33 ! :tex: $-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_r$
    Integer, Parameter :: induct_vmbp_r           = indoff+34 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\boldsymbol{B'}\right)\right]_r$
    Integer, Parameter :: induct_diff_bp_r        = indoff+35 ! :tex: $-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\right]_r$

    Integer, Parameter :: induct_shear_vmbp_theta = indoff+36 ! :tex: $\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\theta$
    Integer, Parameter :: induct_comp_vmbp_theta  = indoff+37 ! :tex: $-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)B'_\theta$
    Integer, Parameter :: induct_advec_vmbp_theta = indoff+38 ! :tex: $-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_\theta$
    Integer, Parameter :: induct_vmbp_theta       = indoff+39 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\boldsymbol{B'}\right)\right]_\theta$
    Integer, Parameter :: induct_diff_bp_theta    = indoff+40 ! :tex: $-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\right]_\theta$

    Integer, Parameter :: induct_shear_vmbp_phi   = indoff+41 ! :tex: $\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi$
    Integer, Parameter :: induct_comp_vmbp_phi    = indoff+42 ! :tex: $-\left(\overline{\boldsymbol{\nabla}}\cdot\overline{\boldsymbol{v}} \right)B'_\phi$
    Integer, Parameter :: induct_advec_vmbp_phi   = indoff+43 ! :tex: $-\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_\phi$
    Integer, Parameter :: induct_vmbp_phi         = indoff+44 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\overline{\boldsymbol{v}}\times\boldsymbol{B'}\right)\right]_\phi$
    Integer, Parameter :: induct_diff_bp_phi      = indoff+45 ! :tex: $-c_7\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\right]_\phi$

    !--------------- Terms involving v' x <B>
    Integer, Parameter :: induct_shear_vpbm_r     = indoff+46 ! :tex: $\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r$
    Integer, Parameter :: induct_comp_vpbm_r      = indoff+47 ! :tex: $-\left(\overline{\boldsymbol{\nabla}}\cdot\boldsymbol{v'} \right)\overline{B_r}$
    Integer, Parameter :: induct_advec_vpbm_r     = indoff+48 ! :tex: $-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_r$
    Integer, Parameter :: induct_vpbm_r           = indoff+49 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\overline{\boldsymbol{B}}\right)\right]_r$

    Integer, Parameter :: induct_shear_vpbm_theta = indoff+50 ! :tex: $\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta$
    Integer, Parameter :: induct_comp_vpbm_theta  = indoff+51 ! :tex: $-\left(\overline{\boldsymbol{\nabla}}\cdot\boldsymbol{v'} \right)\overline{B_\theta}$
    Integer, Parameter :: induct_advec_vpbm_theta = indoff+52 ! :tex: $-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_\theta$
    Integer, Parameter :: induct_vpbm_theta       = indoff+53 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\overline{\boldsymbol{B}}\right)\right]_\theta$

    Integer, Parameter :: induct_shear_vpbm_phi   = indoff+54 ! :tex: $\left[\overline{\boldsymbol{B}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi$
    Integer, Parameter :: induct_comp_vpbm_phi    = indoff+55 ! :tex: $-\left(\overline{\boldsymbol{\nabla}}\cdot\boldsymbol{v'} \right)\overline{B_\phi}$
    Integer, Parameter :: induct_advec_vpbm_phi   = indoff+56 ! :tex: $-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{B}}\right]_\phi$
    Integer, Parameter :: induct_vpbm_phi         = indoff+57 ! :tex: $\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\overline{\boldsymbol{B}}\right)\right]_\phi$

    !--------------- Terms involving v' x B'
    Integer, Parameter :: induct_shear_vpbp_r     = indoff+58 !  :tex: $\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_r$
    Integer, Parameter :: induct_comp_vpbp_r      = indoff+59 !  :tex: $-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v'} \right)B'_r$
    Integer, Parameter :: induct_advec_vpbp_r     = indoff+60 !  :tex: $-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_r$
    Integer, Parameter :: induct_vpbp_r           = indoff+61 !  :tex: $\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\right]_r$

    Integer, Parameter :: induct_shear_vpbp_theta = indoff+62 !  :tex: $\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\theta$
    Integer, Parameter :: induct_comp_vpbp_theta  = indoff+63 !  :tex: $-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v'} \right)B'_\theta$
    Integer, Parameter :: induct_advec_vpbp_theta = indoff+64 !  :tex: $-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_\theta$
    Integer, Parameter :: induct_vpbp_theta       = indoff+65 !  :tex: $\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\right]_\theta$

    Integer, Parameter :: induct_shear_vpbp_phi   = indoff+66 !  :tex: $\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi$
    Integer, Parameter :: induct_comp_vpbp_phi    = indoff+67 !  :tex: $-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v'} \right)B'_\phi$
    Integer, Parameter :: induct_advec_vpbp_phi   = indoff+68 !  :tex: $-\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right]_\phi$
    Integer, Parameter :: induct_vpbp_phi         = indoff+69 !  :tex: $\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\right]_\phi$

