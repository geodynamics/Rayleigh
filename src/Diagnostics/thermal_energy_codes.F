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

    !--- Thermal Energies
    Integer, Parameter :: teoffset = pt_off+200 ! :OFFSET CODE:
    Integer, Parameter :: thermal_energy_full = teoffset+1  ! rho_bar T_bar S  :tex: $\mathrm{f}_1\mathrm{f}_4\Theta$
    Integer, Parameter :: thermal_energy_p    = teoffset+2  ! rho_bar T_bar S' :tex: $\mathrm{f}_1\mathrm{f}_4\Theta$
    Integer, Parameter :: thermal_energy_m    = teoffset+3  ! rho_bar T_bar <S> :tex: $\mathrm{f}_1\mathrm{f}_4\overline{\Theta}$

    !--- Enthalpy
    Integer, Parameter :: enthalpy_full = teoffset+4 ! :tex: $c_P \hat{\rho}T$
    Integer, Parameter :: enthalpy_p = teoffset+5    ! :tex: $c_P \hat{\rho}T'$
    Integer, Parameter :: enthalpy_m = teoffset+6    ! :tex: $c_P \hat{\rho}\overline{T}$

    !--- Quadratics
    Integer, Parameter :: thermal_energy_sq = teoffset+7  ! rho_bar T_bar S  :tex: $\left(\mathrm{f}_1\mathrm{f}_4\Theta\right)^2$
    Integer, Parameter :: thermal_energyp_sq    = teoffset+8  ! rho_bar T_bar S' :tex: $\left(\mathrm{f}_1\mathrm{f}_4\Theta\right)^2$
    Integer, Parameter :: thermal_energym_sq    = teoffset+9  ! rho_bar T_bar <S> :tex: $\left(\mathrm{f}_1\mathrm{f}_4\overline{\Theta}\right)^2$

    Integer, Parameter :: enthalpy_sq = teoffset+10 ! :tex: $\left(c_P \hat{\rho}T\right)^2$
    Integer, Parameter :: enthalpyp_sq = teoffset+11    ! :tex: $\left(c_P \hat{\rho}T'\right)^2$
    Integer, Parameter :: enthalpym_sq = teoffset+12    ! :tex: $\left(c_P \hat{\rho}\overline{T}\right)^2$
