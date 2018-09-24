    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Turbulent KE Outputs
    Integer, Parameter :: turbke_offset = custom_offset+500

    Integer, Parameter :: production_buoyant_pKE   = turbke_offset + 1  ! Buoyant Production of turbulent kinetic energy
    !Integer, Parameter :: production_shear_pKE     = turbke_offset + 2     ! Shear Production of turbulent kinetic energy
    Integer, Parameter :: dissipation_viscous_pKE  = turbke_offset + 3  ! Viscous Dissipation of turbulent kinetic energy
    Integer, Parameter :: transport_pressure_pKE   = turbke_offset + 4  ! Pressure Transport of turbulent kinetic energy
    Integer, Parameter :: transport_viscous_pKE    = turbke_offset + 5      ! Viscous Transport of turbulent kinetic energy
    Integer, Parameter :: transport_turbadvect_pKE = turbke_offset + 6  ! Turbulent Advective Transport of turbulent kinetic energy
    Integer, Parameter :: transport_meanadvect_pKE = turbke_offset + 7  ! Mean Advective Transport of turbulent kinetic energy
    Integer, Parameter :: rflux_pressure_pKE       = turbke_offset + 8          ! Radial Pressure Flux of turbulent kinetic energy
    Integer, Parameter :: rflux_viscous_pKE        = turbke_offset + 9              ! Radial Viscous Flux of turbulent kinetic energy  $
    Integer, Parameter :: rflux_turbadvect_pKE     = turbke_offset + 10     ! Radial Turbulent Advective Flux of turbulent kinetic ener$
    Integer, Parameter :: rflux_meanadvect_pKE     = turbke_offset + 11     ! Radial Mean Advective Flux of turbulent kinetic energy
    Integer, Parameter :: thetaflux_pressure_pKE   = turbke_offset + 12 ! Colatitudinal Pressure Flux of turbulent kinetic energy
    Integer, Parameter :: thetaflux_viscous_pKE    = turbke_offset + 13 ! Colatitudinal Viscous Flux of turbulent kinetic energy
    Integer, Parameter :: thetaflux_turbadvect_pKE = turbke_offset + 14 ! Colatitudinal Turbulent Advective Flux of turbulent kinetic e$
    Integer, Parameter :: thetaflux_meanadvect_pKE = turbke_offset + 15 ! Colatitudinal Mean Advective Flux of turbulent kinetic energy


