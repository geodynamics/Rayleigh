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
    !       Magnetic Energy Equation Terms
    !       Source terms and fluxes
    !       As far as the Reynold decomposition goes,
    !       Only those terms which contribute to a zonal mean are included here
    Integer, Parameter :: meq_off = keq_off +100 ! :OFFSET CODE:

    ! This one's a little different in that we place the fluxes first
    Integer, Parameter :: ecrossb_r     = meq_off+1      ! :tex: $\left[\left(\boldsymbol{v}\times\boldsymbol{B}-\eta\boldsymbol{\mathcal{J}}\right)\times\boldsymbol{B}\right]_r$
    Integer, Parameter :: ecrossb_theta = meq_off+2      ! :tex: $\left[\left(\boldsymbol{v}\times\boldsymbol{B}-\eta\boldsymbol{\mathcal{J}}\right)\times\boldsymbol{B}\right]_\theta$
    Integer, Parameter :: ecrossb_phi   = meq_off+3      ! :tex: $\left[\left(\boldsymbol{v}\times\boldsymbol{B}-\eta\boldsymbol{\mathcal{J}}\right)\times\boldsymbol{B}\right]_\phi $

    Integer, Parameter :: ecrossb_ppp_r     = meq_off+4  ! :tex: $\left[\left(\boldsymbol{v'}\times\boldsymbol{B'}-\eta\boldsymbol{\mathcal{J'}}\right)\times\boldsymbol{B'}\right]_r $
    Integer, Parameter :: ecrossb_ppp_theta = meq_off+5  ! :tex: $\left[\left(\boldsymbol{v'}\times\boldsymbol{B'}-\eta\boldsymbol{\mathcal{J'}}\right)\times\boldsymbol{B'}\right]_\theta $
    Integer, Parameter :: ecrossb_ppp_phi   = meq_off+6  ! :tex: $\left[\left(\boldsymbol{v'}\times\boldsymbol{B'}-\eta\boldsymbol{\mathcal{J'}}\right)\times\boldsymbol{B'}\right]_\phi$

    Integer, Parameter :: ecrossb_mmm_r     = meq_off+7  ! :tex: $\left[\left(\boldsymbol{\overline{v}}\times\boldsymbol{\overline{B}}-\eta\boldsymbol{\mathcal{\overline{J}}}\right)\times\boldsymbol{\overline{B}}\right]_r $
    Integer, Parameter :: ecrossb_mmm_theta = meq_off+8  ! :tex: $\left[\left(\boldsymbol{\overline{v}}\times\boldsymbol{\overline{B}}-\eta\boldsymbol{\mathcal{\overline{J}}}\right)\times\boldsymbol{\overline{B}}\right]_\theta $
    Integer, Parameter :: ecrossb_mmm_phi = meq_off+9  ! :tex: $\left[\left(\boldsymbol{\overline{v}}\times\boldsymbol{\overline{B}}-\eta\boldsymbol{\mathcal{\overline{J}}}\right)\times\boldsymbol{\overline{B}}\right]_\phi $

    Integer, Parameter :: ecrossb_ppm_r     = meq_off+10  ! :tex: $\left[\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\times\boldsymbol{\overline{B}}\right]_r  $
    Integer, Parameter :: ecrossb_ppm_theta = meq_off+11  ! :tex: $\left[\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\times\boldsymbol{\overline{B}}\right]_\theta $
    Integer, Parameter :: ecrossb_ppm_phi   = meq_off+12  ! :tex: $\left[\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\times\boldsymbol{\overline{B}}\right]_\phi $

    Integer, Parameter :: ecrossb_pmp_r     = meq_off+13  ! :tex: $\left[\left(\boldsymbol{v'}\times\boldsymbol{\overline{B}}\right)\times\boldsymbol{B'}\right]_r $
    Integer, Parameter :: ecrossb_pmp_theta = meq_off+14  ! :tex: $\left[\left(\boldsymbol{v'}\times\boldsymbol{\overline{B}}\right)\times\boldsymbol{B'}\right]_\theta $
    Integer, Parameter :: ecrossb_pmp_phi   = meq_off+15  ! :tex: $\left[\left(\boldsymbol{v'}\times\boldsymbol{\overline{B}}\right)\times\boldsymbol{B'}\right]_\phi $

    Integer, Parameter :: ecrossb_mpp_r     = meq_off+16  ! :tex: $\left[\left(\boldsymbol{\overline{v}}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_r $
    Integer, Parameter :: ecrossb_mpp_theta = meq_off+17  ! :tex: $\left[\left(\boldsymbol{\overline{v}}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_\theta $
    Integer, Parameter :: ecrossb_mpp_phi   = meq_off+18  ! :tex: $\left[\left(\boldsymbol{\overline{v}}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_\phi $

    !Next, the source terms (ppp etc. order =>  b dot del cross (v x b)
    Integer, Parameter :: induct_work     = meq_off+19 ! :tex: $\boldsymbol{B}\cdot\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v}\times\boldsymbol{B}\right)\right] $
    Integer, Parameter :: induct_work_ppp = meq_off+20 ! :tex: $\boldsymbol{B'}\cdot\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\right] $
    Integer, Parameter :: induct_work_ppm = meq_off+21 ! :tex: $\boldsymbol{B'}\cdot\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\boldsymbol{\overline{B}}\right)\right] $
    Integer, Parameter :: induct_work_pmp = meq_off+22 ! :tex: $\boldsymbol{B'}\cdot\left[\boldsymbol{\nabla}\times\left(\boldsymbol{\overline{v}}\times\boldsymbol{B'}\right)\right] $
    Integer, Parameter :: induct_work_mpp = meq_off+23 ! :tex: $\boldsymbol{\overline{B}}\cdot\left[\boldsymbol{\nabla}\times\left(\boldsymbol{v'}\times\boldsymbol{B'}\right)\right] $
    Integer, Parameter :: induct_work_mmm = meq_off+24 ! :tex: $\boldsymbol{\overline{B}}\cdot\left[\boldsymbol{\nabla}\times\left(\boldsymbol{\overline{v}}\times\boldsymbol{\overline{B}}\right)\right] $


    Integer, Parameter :: ishear_work = meq_off+25 ! :tex: $\boldsymbol{B}\cdot\left[\boldsymbol{B}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right] $
    Integer, Parameter :: iadvec_work = meq_off+26 ! :tex: $-\boldsymbol{B}\cdot\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{B}\right] $
    Integer, Parameter :: icomp_work  = meq_off+27 ! :tex: $-\boldsymbol{B}\cdot\left(\boldsymbol{\nabla}\cdot\boldsymbol{v} \right)\boldsymbol{B}$

    Integer, Parameter :: ishear_work_pmp = meq_off+28 ! :tex: $\boldsymbol{B'}\cdot\left[\boldsymbol{\overline{B}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right] $
    Integer, Parameter :: iadvec_work_pmp = meq_off+29 ! :tex: $-\boldsymbol{B'}\cdot\left[\boldsymbol{\overline{v}}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right] $
    Integer, Parameter :: icomp_work_pmp  = meq_off+30 ! :tex: $-\boldsymbol{B'}\cdot\left(\boldsymbol{\nabla}\cdot\boldsymbol{\overline{v}} \right)\boldsymbol{B'} $

    Integer, Parameter :: ishear_work_ppm = meq_off+31 ! :tex: $\boldsymbol{B'}\cdot\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\boldsymbol{\overline{v}}\right] $
    Integer, Parameter :: iadvec_work_ppm = meq_off+32 ! :tex: $-\boldsymbol{B'}\cdot\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{\overline{B}}\right] $
    Integer, Parameter :: icomp_work_ppm  = meq_off+33 ! :tex: $-\boldsymbol{B'}\cdot\left(\boldsymbol{\nabla}\cdot\boldsymbol{v'} \right)\boldsymbol{\overline{B}}  $

    Integer, Parameter :: ishear_work_mmm = meq_off+34 ! :tex: $\boldsymbol{\overline{B}}\cdot\left[\boldsymbol{\overline{B}}\cdot\boldsymbol{\nabla}\,\boldsymbol{\overline{v}}\right] $
    Integer, Parameter :: iadvec_work_mmm = meq_off+35 ! :tex: $-\boldsymbol{\overline{B}}\cdot\left[\boldsymbol{\overline{v}}\cdot\boldsymbol{\nabla}\boldsymbol{\overline{B}}\right] $
    Integer, Parameter :: icomp_work_mmm  = meq_off+36 ! :tex: $-\boldsymbol{\overline{B}}\cdot\left(\boldsymbol{\nabla}\cdot\boldsymbol{\overline{v}} \right)\boldsymbol{\overline{B}} $

    Integer, Parameter :: ishear_work_mpp = meq_off+37 ! :tex: $\boldsymbol{\overline{B}}\cdot\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right] $
    Integer, Parameter :: iadvec_work_mpp = meq_off+38 ! :tex: $-\boldsymbol{\overline{B}}\cdot\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right] $
    Integer, Parameter :: icomp_work_mpp  = meq_off+39 ! :tex: $-\boldsymbol{\overline{B}}\cdot\left(\boldsymbol{\nabla}\cdot\boldsymbol{v'} \right)\boldsymbol{B'} $

    Integer, Parameter :: ishear_work_ppp = meq_off+40 ! :tex: $\boldsymbol{B'}\cdot\left[\boldsymbol{B'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]  $
    Integer, Parameter :: iadvec_work_ppp = meq_off+41 ! :tex: $-\boldsymbol{B'}\cdot\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{B'}\right] $
    Integer, Parameter :: icomp_work_ppp  = meq_off+42 ! :tex: $-\boldsymbol{B'}\cdot\left(\boldsymbol{\nabla}\cdot\boldsymbol{v'} \right)\boldsymbol{B'} $

    Integer, Parameter :: idiff_work    = meq_off+43 ! :tex: $-c_7\boldsymbol{B}\cdot\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B}\right)\right] $
    Integer, Parameter :: idiff_work_pp = meq_off+44 ! :tex: $-c_7\boldsymbol{B'}\cdot\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\right] $
    Integer, Parameter :: idiff_work_mm = meq_off+45 ! :tex: $-c_7\boldsymbol{\overline{B}}\cdot\left[ \boldsymbol{\nabla}\times\left(\mathrm{f}_7\boldsymbol{\nabla}\times\boldsymbol{\overline{B}}\right)\right] $
