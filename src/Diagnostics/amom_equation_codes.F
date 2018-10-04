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
    !       Mean Angular Momentum Equation Terms
    !       i.e., terms that govern the evolution of d_by_dt{r sin\theta rho_hat <v>_\phi}
    Integer, Parameter :: amoff = indoff+200 ! :OFFSET CODE:
    ! Source terms
    !/////////////////////////// Lorentz Forces ///////////////////////////////


    Integer, Parameter :: samom_advec_pp = amoff+1 ! :tex: $r\,\mathrm{sin}\theta\mathrm{f}_1\left[\boldsymbol{v'}\cdot\boldsymbol{\nabla}\boldsymbol{v'}\right]_\phi$
    Integer, Parameter :: samom_advec_mm = amoff+2 ! :tex: $r\,\mathrm{sin}\theta\mathrm{f}_1\left[\overline{\boldsymbol{v}}\cdot\boldsymbol{\nabla}\overline{\boldsymbol{v}}\right]_\phi$

    Integer, Parameter :: samom_coriolis = amoff+3 ! :tex: $-c_1\mathrm{f}_1r\mathrm{sin}\theta \left(\mathrm{cos}\theta\, v_\theta + \mathrm{sin}\theta\, v_r\right)$

    Integer, Parameter :: samom_diffusion = amoff+4 ! :tex: $r\,\mathrm{sin}\theta\left[\boldsymbol{\nabla}\cdot\overline{\boldsymbol{\mathcal{D}}}\right]_\phi$

    Integer, Parameter :: samom_lorentz_mm   = amoff+5 ! :tex: $r\,\mathrm{sin}\theta c_4\left[\left(\boldsymbol{\nabla}\times\overline{\boldsymbol{B}}\right)\times\overline{\boldsymbol{B}}\right]_\phi$
    Integer, Parameter :: samom_lorentz_pp   = amoff+6 ! :tex: $r\,\mathrm{sin}\theta c_4\left[\left(\boldsymbol{\nabla}\times\boldsymbol{B'}\right)\times\boldsymbol{B'}\right]_\phi$


    ! Fluxes
    Integer, Parameter :: famom_fluct_r     = amoff+7 ! :tex:  $\mathrm{f}_1r\mathrm{sin}\theta v'_r v'_\phi $
    Integer, Parameter :: famom_fluct_theta = amoff+8 ! :tex:  $\mathrm{f}_1r\mathrm{sin}\theta v'_\theta v'_\phi $

    Integer, Parameter :: famom_dr_r        = amoff+9  ! :tex:  $\mathrm{f}_1r\mathrm{sin}\theta\, \overline{v_r}\,\, \overline{v_\phi} $
    Integer, Parameter :: famom_dr_theta    = amoff+10 ! :tex:  $\mathrm{f}_1r\mathrm{sin}\theta\, \overline{v_\theta}\,\, \overline{v_\phi} $

    Integer, Parameter :: famom_mean_r      = amoff+11 ! :tex:  $\frac{c_1}{2}\mathrm{f}_1r^2\mathrm{sin}^2\theta\, \overline{v_r}$
    Integer, Parameter :: famom_mean_theta  = amoff+12 ! :tex:  $\frac{c_1}{2}\mathrm{f}_1r^2\mathrm{sin}^2\theta\, \overline{v_\theta}$

    Integer, Parameter :: famom_diff_r       = amoff+13 ! :tex: $\mathrm{f}_1\nu\mathrm{sin}\theta\left(v_\phi-r\frac{\partial\,\overline{v_\phi}}{\partial r}   \right)$
    Integer, Parameter :: famom_diff_theta   = amoff+14 ! :tex: $\mathrm{f}_1\nu\left(\mathrm{cos}\theta\, \overline{v_\phi}-\mathrm{sin}\theta\frac{\partial\, \overline{v_\phi}}{\partial \theta}\right)$

    Integer, Parameter :: famom_maxstr_r     = amoff+15 ! :tex:  $-r\mathrm{sin}\theta c_4\,B'_r\,B'_\phi$
    Integer, Parameter :: famom_maxstr_theta = amoff+16 ! :tex:  $-r\mathrm{sin}\theta c_4\,B'_\theta\,B'_\phi$

    Integer, Parameter :: famom_magtor_r     = amoff+17 !:tex:  $-r\mathrm{sin}\theta c_4\,\overline{B_r}\,\overline{B_\phi}$
    Integer, Parameter :: famom_magtor_theta = amoff+18 ! :tex:  $-r\mathrm{sin}\theta c_4\,\overline{B_\theta}\,\overline{B_\phi}$


