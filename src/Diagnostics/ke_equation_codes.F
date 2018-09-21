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
    !       Kinetic Equation Terms
    !       Source terms and fluxes
    !       Only those terms which contribute to a zonal mean are included here
    Integer, Parameter :: keq_off = amoff +100 ! :OFFSET CODE:

    ! First the source terms
    Integer, Parameter :: press_work    = keq_off+1 ! :tex: $-c_3\mathrm{f}_1\boldsymbol{v}\cdot\boldsymbol{\nabla}\left(\frac{P}{\mathrm{f}_1}\right)$
    Integer, Parameter :: press_work_pp = keq_off+2 ! :tex: $-c_3\mathrm{f}_1\boldsymbol{v'}\cdot\boldsymbol{\nabla}\left(\frac{P'}{\mathrm{f}_1}\right)$
    Integer, Parameter :: press_work_mm = keq_off+3 ! :tex: $-c_3\mathrm{f}_1\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\left(\frac{\overline{P}}{\mathrm{f}_1}\right)$

    Integer, Parameter :: buoy_work    = keq_off+4  ! :tex: $c_2v_r\mathrm{f}_2\Theta$
    Integer, Parameter :: buoy_work_pp = keq_off+5  ! :tex: $c_2v'_r\mathrm{f}_2\Theta'$
    Integer, Parameter :: buoy_work_mm = keq_off+6  ! :tex: $c_2\overline{v_r}\mathrm{f}_2\overline{\Theta}$

    Integer, Parameter :: visc_work    = keq_off+7  ! :tex: $c_5\boldsymbol{v}\cdot\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\right]$
    Integer, Parameter :: visc_work_pp = keq_off+8  ! :tex: $c_5\boldsymbol{v'}\cdot\left[\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D'}}\right]$
    Integer, Parameter :: visc_work_mm = keq_off+9  ! :tex: $c_5\overline{\boldsymbol{v}}\cdot\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]$

    Integer, Parameter :: advec_work     = keq_off+10 ! :tex: $\mathrm{f}_1\boldsymbol{v}\cdot\left[\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}\right]$
    Integer, Parameter :: advec_work_ppp = keq_off+11 ! :tex: $\mathrm{f}_1\boldsymbol{v'}\cdot\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]$
    Integer, Parameter :: advec_work_mpp = keq_off+12 ! :tex: $\mathrm{f}_1\overline{\boldsymbol{v}}\cdot\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]$
    Integer, Parameter :: advec_work_pmp = keq_off+13 ! :tex: $\mathrm{f}_1\boldsymbol{v'}\cdot\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]$
    Integer, Parameter :: advec_work_ppm = keq_off+14 ! :tex: $\mathrm{f}_1\boldsymbol{v'}\cdot\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]$
    Integer, Parameter :: advec_work_mmm = keq_off+15 ! :tex: $\mathrm{f}_1\overline{\boldsymbol{v}}\cdot\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]$


    Integer, Parameter :: mag_work     = keq_off+16 ! :tex: $c_4\boldsymbol{v}\cdot\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B}\right]$
    Integer, Parameter :: mag_work_ppp = keq_off+18 ! :tex: $c_4\boldsymbol{v'}\cdot\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]$
    Integer, Parameter :: mag_work_mpp = keq_off+19 ! :tex: $c_4\overline{\boldsymbol{v}}\cdot\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right] $
    Integer, Parameter :: mag_work_pmp = keq_off+20 ! :tex: $c_4\boldsymbol{v'}\cdot\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\boldsymbol{B'}\right] $
    Integer, Parameter :: mag_work_ppm = keq_off+21 ! :tex: $c_4\boldsymbol{v'}\cdot\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\overline{\boldsymbol{B}}\right]$
    Integer, Parameter :: mag_work_mmm = keq_off+22 ! :tex: $c_4\overline{\boldsymbol{v}}\cdot\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right] $

    ! Now the fluxes
    Integer, Parameter :: ke_flux_radial = keq_off+23    ! :tex: $\frac{1}{2}\mathrm{f}_1v_r\,v^2$
    Integer, Parameter :: ke_flux_theta  = keq_off+24    ! :tex: $\frac{1}{2}\mathrm{f}_1v_\theta\,v^2 $
    Integer, Parameter :: ke_flux_phi    = keq_off+25    ! :tex: $\frac{1}{2}\mathrm{f}_1v_\phi\,v^2$

    Integer, Parameter :: mke_mflux_radial = keq_off+26  ! :tex: $\frac{1}{2}\mathrm{f}_1\,\overline{v_r}\,\overline{v}^2 $
    Integer, Parameter :: mke_mflux_theta  = keq_off+27  ! :tex: $\frac{1}{2}\mathrm{f}_1\,\overline{v_\theta}\,\overline{v}^2$
    Integer, Parameter :: mke_mflux_phi    = keq_off+28  ! :tex: $\frac{1}{2}\mathrm{f}_1\,\overline{v_\phi}\,\overline{v}^2$

    Integer, Parameter :: pke_mflux_radial  = keq_off+29 ! :tex: $\frac{1}{2}\mathrm{f}_1\,\overline{v_r}\,{v'}^2 $
    Integer, Parameter :: pke_mflux_theta   = keq_off+30 ! :tex: $\frac{1}{2}\mathrm{f}_1\,\overline{v_\theta}\,{v'}^2  $
    Integer, Parameter :: pke_mflux_phi     = keq_off+31 ! :tex: $\frac{1}{2}\mathrm{f}_1\,\overline{v_\phi}\,{v'}^2  $

    Integer, Parameter :: pke_pflux_radial  = keq_off+32 ! :tex: $\frac{1}{2}\mathrm{f}_1\,v'_r\,{v'}^2 $
    Integer, Parameter :: pke_pflux_theta   = keq_off+33 ! :tex: $\frac{1}{2}\mathrm{f}_1\,v'_\theta\,{v'}^2 $
    Integer, Parameter :: pke_pflux_phi     = keq_off+34 ! :tex: $\frac{1}{2}\mathrm{f}_1\,v'_\phi\,{v'}^2 $

    Integer, Parameter :: visc_flux_r     = keq_off+35 ! :tex: $c_5\left[\boldsymbol{v}\cdot\boldsymbol{\mathcal{D}} \right]_r$
    Integer, Parameter :: visc_flux_theta = keq_off+36 ! :tex: $c_5\left[\boldsymbol{v}\cdot\boldsymbol{\mathcal{D}} \right]_\theta $
    Integer, Parameter :: visc_flux_phi   = keq_off+37 ! :tex: $c_5\left[\boldsymbol{v}\cdot\boldsymbol{\mathcal{D}} \right]_\phi$

    Integer, Parameter :: visc_fluxpp_r     = keq_off+38 ! :tex: $c_5\left[\boldsymbol{v'}\cdot\boldsymbol{\mathcal{D'}} \right]_r$
    Integer, Parameter :: visc_fluxpp_theta = keq_off+39 ! :tex: $c_5\left[\boldsymbol{v'}\cdot\boldsymbol{\mathcal{D'}} \right]_\theta$
    Integer, Parameter :: visc_fluxpp_phi   = keq_off+40 ! :tex: $c_5\left[\boldsymbol{v'}\cdot\boldsymbol{\mathcal{D'}} \right]_\phi$

    Integer, Parameter :: visc_fluxmm_r     = keq_off+41 ! :tex: $c_5\left[\boldsymbol{\overline{v}}\cdot\boldsymbol{\overline{\mathcal{D}}} \right]_r$
    Integer, Parameter :: visc_fluxmm_theta = keq_off+42 ! :tex: $c_5\left[\boldsymbol{\overline{v}}\cdot\boldsymbol{\overline{\mathcal{D}}} \right]_\theta$
    Integer, Parameter :: visc_fluxmm_phi   = keq_off+43 ! :tex: $c_5\left[\boldsymbol{\overline{v}}\cdot\boldsymbol{\overline{\mathcal{D}}} \right]_\phi$

    Integer, Parameter :: press_flux_r     = keq_off+44  ! :tex: $-c_3v_r P$
    Integer, Parameter :: press_flux_theta = keq_off+45  ! :tex: $-c_3v_\theta P$
    Integer, Parameter :: press_flux_phi   = keq_off+46  ! :tex: $-c_3v_\phi P$

    Integer, Parameter :: press_fluxpp_r     = keq_off+47  ! :tex: $-c_3v'_r P'$
    Integer, Parameter :: press_fluxpp_theta = keq_off+48  ! :tex: $-c_3v'_\theta P'$
    Integer, Parameter :: press_fluxpp_phi   = keq_off+49  ! :tex: $-c_3v'_\phi P'$

    Integer, Parameter :: press_fluxmm_r     = keq_off+50  ! :tex: $-c_3\overline{v_r}\, \overline{P}$
    Integer, Parameter :: press_fluxmm_theta = keq_off+51  ! :tex: $-c_3\overline{v_\theta}\, \overline{P}$
    Integer, Parameter :: press_fluxmm_phi   = keq_off+52  ! :tex: $-c_3\overline{v_\phi}\, \overline{P}$

    Integer, Parameter :: production_shear_ke  = keq_off+53 ! :tex: $ -- $
    Integer, Parameter :: production_shear_pke = keq_off+54 ! :tex: $ -- $
    Integer, Parameter :: production_shear_mke = keq_off+55 ! :tex: $ -- $

