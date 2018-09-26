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

    Integer, Parameter :: rhovoffset = voffset+200 ! :OFFSET CODE:
    !------------ Mass Flux ---------------------!
    ! Full
    Integer, Parameter :: rhov_r      = rhovoffset+1 ! :tex: $\mathrm{f}_1v_r$
    Integer, Parameter :: rhov_theta  = rhovoffset+2 ! :tex: $\mathrm{f}_1v_\theta$
    Integer, Parameter :: rhov_phi    = rhovoffset+3 ! :tex: $\mathrm{f}_1v_\phi$
    ! Fluctuating
    Integer, Parameter :: rhovp_r     = rhovoffset+4 ! :tex: $\mathrm{f}_1v_r'$
    Integer, Parameter :: rhovp_theta = rhovoffset+5 ! :tex: $\mathrm{f}_1v_\theta'$
    Integer, Parameter :: rhovp_phi   = rhovoffset+6 ! :tex: $\mathrm{f}_1v_\phi'$
    ! Mean
    Integer, Parameter :: rhovm_r     = rhovoffset+7 ! :tex: $\mathrm{f}_1\overline{v_r}$
    Integer, Parameter :: rhovm_theta = rhovoffset+8 ! :tex: $\mathrm{f}_1\overline{v_\theta}$
    Integer, Parameter :: rhovm_phi   = rhovoffset+9 ! :tex: $\mathrm{f}_1\overline{v_\phi}$
