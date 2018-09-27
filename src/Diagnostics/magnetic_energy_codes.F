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
    !           Magnetic Energies
    Integer, Parameter :: meoffset = joffset+100 ! :OFFSET CODE:

    Integer, Parameter :: magnetic_energy = meoffset+1 ! :tex: $\frac{1}{2}c_4\boldsymbol{B}^2$
    Integer, Parameter :: radial_me       = meoffset+2 ! :tex: $\frac{1}{2}c_4B_r^2$
    Integer, Parameter :: theta_me        = meoffset+3 ! :tex: $\frac{1}{2}c_4B_\theta^2$
    Integer, Parameter :: phi_me        = meoffset+4 ! :tex: $\frac{1}{2}c_4B_\phi^2$

    Integer, Parameter :: mmagnetic_energy = meoffset+5 ! :tex: $\frac{1}{2}c_4\overline{\boldsymbol{B}}^2$
    Integer, Parameter :: radial_mme       = meoffset+6 ! :tex: $\frac{1}{2}c_4\overline{B_r}^2$
    Integer, Parameter :: theta_mme        = meoffset+7 ! :tex: $\frac{1}{2}c_4\overline{B_\theta}^2$
    Integer, Parameter :: phi_mme        = meoffset+8 ! :tex: $\frac{1}{2}c_4\overline{B_\phi}^2$

    Integer, Parameter :: pmagnetic_energy = meoffset+9  ! :tex: $\frac{1}{2}c_4\boldsymbol{B'}^2$
    Integer, Parameter :: radial_pme       = meoffset+10 ! :tex: $\frac{1}{2}c_4{B_r'}^2$
    Integer, Parameter :: theta_pme        = meoffset+11 ! :tex: $\frac{1}{2}c_4{B_\theta'}^2$
    Integer, Parameter :: phi_pme        = meoffset+12 ! :tex: $\frac{1}{2}c_4{B_\phi'}^2$

