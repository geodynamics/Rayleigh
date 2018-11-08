    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Turbulent KE Outputs
    Integer, Parameter :: turbke_offset = custom_offset+500

    Integer, Parameter :: production_buoyant_pKE   = turbke_offset + 1  ! :tex: $\mathrm{f}_2{\Theta^\prime}{v^\prime}_r$
    Integer, Parameter :: production_shear2_pKE    = turbke_offset + 2  ! :tex: $-\mathrm{f}_1{v^\prime_i}{v^\prime_j}\overline{e}_{ij}$
    Integer, Parameter :: dissipation_viscous_pKE  = turbke_offset + 3  ! :tex: $2\mathrm{f}_1\mathrm{f}_3\left[e^\prime_{ij}e^\prime_{ij}-\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}^\prime\right)^2/3\right]$

    Integer, Parameter :: transport_pressure_pKE   = turbke_offset + 4  ! :tex: $-\boldsymbol{\nabla}\cdot\left(P^\prime\boldsymbol{v}^\prime\right)$
    Integer, Parameter :: transport_viscous_pKE    = turbke_offset + 5  ! :tex: $\boldsymbol{\nabla}\cdot\left(\boldsymbol{\mathcal{D}}^\prime\cdot\boldsymbol{v}^\prime\right)$
    Integer, Parameter :: transport_turbadvect_pKE = turbke_offset + 6  ! :tex: $-\boldsymbol{\nabla}\cdot\left(\frac{1}{2}\mathrm{f}_1{v^\prime}^2\boldsymbol{v}^\prime\right)$
    Integer, Parameter :: transport_meanadvect_pKE = turbke_offset + 7  ! :tex: $-\boldsymbol{\nabla}\cdot\left(\frac{1}{2}\mathrm{f}_1{v^\prime}^2\overline{\boldsymbol{v}}\right)$

    Integer, Parameter :: rflux_pressure_pKE       = turbke_offset + 8  ! :tex: $P^\prime{v_r^\prime}$
    Integer, Parameter :: rflux_viscous_pKE        = turbke_offset + 9  ! :tex: $-\left[\boldsymbol{\mathcal{D}}^\prime\cdot\boldsymbol{v}^\prime\right]_r$
    Integer, Parameter :: rflux_turbadvect_pKE     = turbke_offset + 10 ! :tex: $\frac{1}{2}\mathrm{f}_1{v^\prime}^2v_r^\prime$
    Integer, Parameter :: rflux_meanadvect_pKE     = turbke_offset + 11 ! :tex: $\frac{1}{2}\mathrm{f}_1{v^\prime}^2\overline{v}_r$

    Integer, Parameter :: thetaflux_pressure_pKE   = turbke_offset + 12 ! :tex: $P^\prime{v_\theta^\prime}$
    Integer, Parameter :: thetaflux_viscous_pKE    = turbke_offset + 13 ! :tex: $-\left[\boldsymbol{\mathcal{D}}^\prime\cdot\boldsymbol{v}^\prime\right]_\theta$
    Integer, Parameter :: thetaflux_turbadvect_pKE = turbke_offset + 14 ! :tex: $\frac{1}{2}\mathrm{f}_1{v^\prime}^2v_\theta^\prime$
    Integer, Parameter :: thetaflux_meanadvect_pKE = turbke_offset + 15 ! :tex: $\frac{1}{2}\mathrm{f}_1{v^\prime}^2\overline{v}_\theta$

    Integer, Parameter :: phiflux_pressure_pKE     = turbke_offset + 16 ! :tex: $P^\prime{v_\phi^\prime}$
    Integer, Parameter :: phiflux_viscous_pKE      = turbke_offset + 17 ! :tex: $-\left[\boldsymbol{\mathcal{D}}^\prime\cdot\boldsymbol{v}^\prime\right]_\phi$
    Integer, Parameter :: phiflux_turbadvect_pKE   = turbke_offset + 18 ! :tex: $\frac{1}{2}\mathrm{f}_1{v^\prime}^2v_\phi^\prime$
    Integer, Parameter :: phiflux_meanadvect_pKE   = turbke_offset + 19 ! :tex: $\frac{1}{2}\mathrm{f}_1{v^\prime}^2\overline{v}_\phi$

