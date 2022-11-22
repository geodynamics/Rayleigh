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

    !//////////////////////////////////////////////////////////////////
    !  Composition/Scalar (Active)
    Integer, Parameter :: scalar_skip     = 200
    Integer, Parameter :: a_scalar_offset = 10000 ! :OFFSET CODE:
    Integer, Parameter :: p_scalar_offset = 20000 ! :OFFSET CODE:

    !------------ Fields ---------------------!
    ! Full
    Integer, Parameter :: chi   = 1 !  :tex:  $\chi$
    ! Fluctuating
    Integer, Parameter :: chi_p = 2 !  :tex:  $\chi'$
    ! Mean
    Integer, Parameter :: chi_m = 3 !  :tex:  $\overline{\chi}$

    !------------ Radial Derivatives --------------!
    !  Full
    Integer, Parameter :: chi_dr   = 4  ! :tex:  $\frac{\partial \chi}{\partial r}$
    ! Fluctuating
    Integer, Parameter :: chi_p_dr = 5  ! :tex:  $\frac{\partial \chi'}{\partial r}$
    ! Mean
    Integer, Parameter :: chi_m_dr = 6 ! :tex:  $\frac{\partial \overline{\chi}}{\partial r}$

    !------------ Theta Derivatives ---------------!
    Integer, Parameter :: chi_dtheta   = 7  ! :tex:  $\frac{\partial \chi}{\partial \theta}$

    Integer, Parameter :: chi_p_dtheta = 8  ! :tex:  $\frac{\partial \chi'}{\partial \theta}$

    Integer, Parameter :: chi_m_dtheta = 9  ! :tex:  $\frac{\partial \overline{\chi}}{\partial \theta}$

    !------------ Phi Derivatives -----------------!
    Integer, Parameter :: chi_dphi     = 10  ! :tex:  $\frac{\partial \chi}{\partial \phi}$

    Integer, Parameter :: chi_p_dphi   = 11  ! :tex:  $\frac{\partial \chi'}{\partial \phi}$

    Integer, Parameter :: chi_m_dphi   = 12  ! :tex:  $\frac{\partial \overline{\chi}}{\partial \phi}$

