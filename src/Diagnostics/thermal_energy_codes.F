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
