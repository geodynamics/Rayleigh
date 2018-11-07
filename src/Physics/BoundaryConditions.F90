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

Module BoundaryConditions
    Use Math_Constants
    Use ProblemSize
    !Use ReferenceState
    !Use TransportCoefficients
    Implicit None

    Logical :: Fix_Tvar_Top    = .True.
    Logical :: Fix_Tvar_Bottom = .True.
    Logical :: Fix_dTdr_Top    = .False.
    Logical :: Fix_dTdr_Bottom = .False.
    Logical :: Fix_divrt_top = .False.
    Logical :: Fix_divt_top = .False.
    Logical :: Fix_divrfc_top = .False.
    Logical :: Fix_divfc_top = .False.
    Logical :: Fix_poloidalfield_top = .False.
    Logical :: Fix_poloidalfield_bottom = .False.
    Logical :: Impose_Dipole_Field = .False.
    Logical :: fix_tdt_bottom = .false.

    Real*8  :: T_Bottom     = 1.0d0
    Real*8  :: T_Top        = 0.0d0
    Real*8  :: dTdr_Top     = 0.0d0
    Real*8  :: dTdr_Bottom  = 0.0d0
    Real*8  :: C10_bottom = 0.0d0
    Real*8  :: C10_top = 0.0d0
    Real*8  :: C11_bottom = 0.0d0
    Real*8  :: C11_top = 0.0d0
    Real*8  :: C1m1_bottom = 0.0d0
    Real*8  :: C1m1_top = 0.0d0
    Real*8  :: Br_bottom = 0.0d0
    Real*8  :: Dipole_Tilt_Degrees = 0.0d0

    Logical :: Strict_L_Conservation = .false.         ! (In-Progress) Turn on to enforce angular momentum conservation abous x,y, and z-axes
    Logical :: no_slip_boundaries = .false. ! Set to true to use no-slip boundaries.  Stree-free boundaries are the default.
    Logical :: stress_free_top = .true., stress_free_bottom = .true.
    Logical :: no_slip_top = .false., no_slip_bottom = .false.


    Namelist /Boundary_Conditions_Namelist/ Fix_Tvar_Top, Fix_Tvar_Bottom, T_Bottom, T_Top, dTdr_top, dTdr_bottom, &
        fix_dtdr_bottom, fix_dtdr_top, fix_divrt_top, fix_divt_top, fix_divrfc_top, fix_divfc_top, &
        no_slip_boundaries, strict_L_Conservation, fix_poloidalfield_top, fix_poloidalfield_bottom, &
        C10_bottom, C10_top, C11_bottom, C11_top, C1m1_bottom, C1m1_top, Br_bottom, &
        dipole_tilt_degrees, impose_dipole_field, fix_tdt_bottom, no_slip_top, no_slip_bottom, &
        stress_free_top, stress_free_bottom

Contains

    Subroutine Initialize_Boundary_Conditions()
        Implicit None
        Real*8 :: tilt_angle_radians,a,b
        Real*8 :: fsun


        fix_tvar_top = .not. fix_dtdr_top
        fix_tvar_bottom = .not. fix_dtdr_bottom


        If (no_slip_boundaries) Then
            no_slip_top = .true.
            no_slip_bottom = .true.
        Endif
        stress_free_top = .not. no_slip_top
        stress_free_bottom = .not. no_slip_bottom


        If (impose_dipole_field) Then
            fix_poloidalfield_top = .true.
            fix_poloidalfield_bottom = .true.
            tilt_angle_radians = pi/180.0*dipole_tilt_degrees
            a = cos(tilt_angle_radians)*sqrt(4.0d0*Pi/3.0d0)
            b = sin(tilt_angle_radians)*sqrt(4.0d0*Pi/3.0d0)

            ! We use the bottom values to derive the top values
            C10_bottom = a*Br_bottom/2.0d0*(radius(N_r)**3)
            C10_top = C10_bottom*(radius(N_R)/radius(1))

            Write(6,*)'C10_bottom is: ', c10_bottom, c10_top

            C11_bottom = b*Br_bottom/2.0d0*(radius(N_r)**3)
            C11_top = C11_bottom*(radius(N_R)/radius(1))

            C1m1_bottom = 0.0d0*Br_bottom/2.0d0*(radius(N_r)**3)
            C1m1_top = C1m1_bottom*(radius(N_R)/radius(1))

        Endif

    End Subroutine Initialize_Boundary_Conditions

    Subroutine Restore_BoundaryCondition_Defaults()
        Implicit None
        Fix_Tvar_Top    = .True.
        Fix_Tvar_Bottom = .True.
        Fix_dTdr_Top    = .False.
        Fix_dTdr_Bottom = .False.
        Fix_divrt_top   = .False.
        Fix_divt_top    = .False.
        Fix_divrfc_top  = .False.
        Fix_divfc_top   = .False.
        Fix_poloidalfield_top    = .False.
        Fix_poloidalfield_bottom = .False.
        Impose_Dipole_Field      = .False.
        fix_tdt_bottom           = .False.

        T_Bottom     = 1.0d0
        T_Top        = 0.0d0
        dTdr_Top     = 0.0d0
        dTdr_Bottom  = 0.0d0
        C10_bottom = 0.0d0
        C10_top = 0.0d0
        C11_bottom = 0.0d0
        C11_top = 0.0d0
        C1m1_bottom = 0.0d0
        C1m1_top = 0.0d0
        Br_bottom = 0.0d0
        Dipole_Tilt_Degrees = 0.0d0

        Strict_L_Conservation = .false.
        no_slip_boundaries = .false.
    End Subroutine Restore_BoundaryCondition_Defaults
End Module BoundaryConditions
