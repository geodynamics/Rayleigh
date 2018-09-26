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

    !///////////////////////////////////////////////////
    !       Current Density Outputs (Including Ohmic Heating)
    !       This is Curl B -- rename accordingly
    Integer, Parameter :: joffset = boffset+200 ! :OFFSET CODE:

    ! Radial components (full, fluctuating, mean)
    Integer, Parameter :: j_r  = joffset+1      ! :tex: $\mathcal{J}_r$
    Integer, Parameter :: jp_r = joffset+2      ! :tex: $\mathcal{J}_r'$
    Integer, Parameter :: jm_r = joffset+3      ! :tex: $\overline{\mathcal{J}}_r$

    ! Theta components (full, fluctuating, mean)
    Integer, Parameter :: j_theta  = joffset+4      ! :tex: $\mathcal{J}_\theta$
    Integer, Parameter :: jp_theta = joffset+5      ! :tex: $\mathcal{J}_\theta'$
    Integer, Parameter :: jm_theta = joffset+6      ! :tex: $\overline{\mathcal{J}}_\theta$

    ! Phi components (full, fluctuating, mean)
    Integer, Parameter :: j_phi  = joffset+7    ! :tex: $\mathcal{J}_\phi$
    Integer, Parameter :: jp_phi = joffset+8    ! :tex: $\mathcal{J}_\phi'$
    Integer, Parameter :: jm_phi = joffset+9    ! :tex: $\overline{\mathcal{J}}_\phi$

    ! J-squared terms
    Integer, Parameter :: j_sq   = joffset+10 ! :tex: ${\boldsymbol{\mathcal{J}}}\cdot{\boldsymbol{\mathcal{J}}}$
    Integer, Parameter :: jp_sq  = joffset+11 ! :tex: ${\boldsymbol{\mathcal{J}'}}\cdot{\boldsymbol{\mathcal{J}'}}$
    Integer, Parameter :: jm_sq  = joffset+12 ! :tex: $\overline{\boldsymbol{\mathcal{J}}}\cdot\overline{\boldsymbol{\mathcal{J}}}$
    Integer, Parameter :: jpm_sq = joffset+13 ! :tex: $\overline{\boldsymbol{\mathcal{J}}}\cdot{\boldsymbol{\mathcal{J}'}}$

    ! Radial components (full, fluctuating, mean)
    Integer, Parameter :: j_r_sq  = joffset+14      ! :tex: $\left(\mathcal{J}_r\right)^2$
    Integer, Parameter :: jp_r_sq = joffset+15      ! :tex: $\left(\mathcal{J}_r'\right)^2$
    Integer, Parameter :: jm_r_sq = joffset+16      ! :tex: $\left(\overline{\mathcal{J}}_r\right)^2$

    ! Theta components (full, fluctuating, mean)
    Integer, Parameter :: j_theta_sq  = joffset+17      ! :tex: $\left(\mathcal{J}_\theta\right)^2$
    Integer, Parameter :: jp_theta_sq = joffset+18      ! :tex: $\left(\mathcal{J}_\theta'\right)^2$
    Integer, Parameter :: jm_theta_sq = joffset+19      ! :tex: $\left(\overline{\mathcal{J}}_\theta\right)^2$

    ! Phi components (full, fluctuating, mean)
    Integer, Parameter :: j_phi_sq  = joffset+20    ! :tex: $\left(\mathcal{J}_\phi\right)^2$
    Integer, Parameter :: jp_phi_sq = joffset+21    ! :tex: $\left(\mathcal{J}_\phi'\right)^2$
    Integer, Parameter :: jm_phi_sq = joffset+22    ! :tex: $\left(\overline{\mathcal{J}}_\phi\right)^2$
