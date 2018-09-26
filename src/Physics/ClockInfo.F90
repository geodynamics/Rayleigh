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

Module ClockInfo
    Use Controls, Only : max_time_step
    Implicit None
    Logical :: new_timestep = .true.
    Logical :: euler_step = .false.
    Real*8  :: new_deltat, deltat, old_deltat
    Real*8  :: min_dt_change = 0.1d0
    !Real*8  :: max_time_step = 5.0d-4
    !Real*8  :: min_time_step = 1.0d-13
    Real*8  :: old_ab_factor = 1.0d0, new_ab_factor = 1.0d0
    Real*8  :: simulation_time

    Integer :: iteration
    Character*8 :: t_ofmt = '(ES12.5)'    ! For formatted timestep output
    Logical :: output_iteration = .false.
End Module ClockInfo
