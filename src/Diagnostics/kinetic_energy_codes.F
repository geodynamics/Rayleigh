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

    !///////////////////////////////////////////////////////////
    !           Kinetic Energies
    Integer, Parameter :: keoffset = vort_off+100 ! :OFFSET CODE:
    Integer, Parameter :: kinetic_energy = keoffset+1   ! :tex: $\frac{1}{2}\mathrm{f}_1{\boldsymbol v}^2$
    Integer, Parameter :: radial_ke      = keoffset+2   ! :tex: $\frac{1}{2}\mathrm{f}_1v_r^2$
    Integer, Parameter :: theta_ke       = keoffset+3   ! :tex: $\frac{1}{2}\mathrm{f}_1v_\theta^2$
    Integer, Parameter :: phi_ke         = keoffset+4   ! :tex: $\frac{1}{2}\mathrm{f}_1v_\phi^2$

    Integer, Parameter :: mkinetic_energy = keoffset+5 ! :tex: $\frac{1}{2}\mathrm{f}_1\overline{\boldsymbol v}^2$
    Integer, Parameter :: radial_mke      = keoffset+6 ! :tex: $\frac{1}{2}\mathrm{f}_1\overline{v_r}^2$
    Integer, Parameter :: theta_mke       = keoffset+7 ! :tex: $\frac{1}{2}\mathrm{f}_1\overline{v_\theta}^2$
    Integer, Parameter :: phi_mke         = keoffset+8 ! :tex: $\frac{1}{2}\mathrm{f}_1\overline{v_\phi}^2$

    Integer, Parameter :: pkinetic_energy = keoffset+9  ! :tex: $\frac{1}{2}\mathrm{f}_1{\boldsymbol v'}^2$
    Integer, Parameter :: radial_pke      = keoffset+10 ! :tex: $\frac{1}{2}\mathrm{f}_1{v_r'}^2$
    Integer, Parameter :: theta_pke       = keoffset+11 ! :tex: $\frac{1}{2}\mathrm{f}_1{v_\theta'}^2$
    Integer, Parameter :: phi_pke         = keoffset+12 ! :tex: $\frac{1}{2}\mathrm{f}_1{v_\phi'}^2$

    !--- Since density varies with radius, it can be useful to output the
    !--- squared fields, sans density, as well.
    Integer, Parameter :: vsq         = keoffset+13   !  :tex: ${\boldsymbol v}^2$
    Integer, Parameter :: radial_vsq  = keoffset+14   !  :tex: ${v_r}^2$
    Integer, Parameter :: theta_vsq   = keoffset+15   !  :tex: ${v_\theta}^2$
    Integer, Parameter :: phi_vsq     = keoffset+16   !  :tex: ${v_\phi}^2$

    Integer, Parameter :: mvsq        = keoffset+17   ! :tex: $\overline{\boldsymbol v}^2$
    Integer, Parameter :: radial_mvsq = keoffset+18   ! :tex: $\overline{v_r}^2$
    Integer, Parameter :: theta_mvsq  = keoffset+19   ! :tex: $\overline{v_\theta}^2$
    Integer, Parameter :: phi_mvsq    = keoffset+20   ! :tex: $\overline{v_\phi}^2$

    Integer, Parameter :: pvsq        = keoffset+21   ! :tex: ${\boldsymbol v'}^2$
    Integer, Parameter :: radial_pvsq = keoffset+22   ! :tex: ${v_r'}^2$
    Integer, Parameter :: theta_pvsq  = keoffset+23   ! :tex: ${v_\theta'}^2$
    Integer, Parameter :: phi_pvsq    = keoffset+24   ! :tex: ${v_\phi'}^2$
