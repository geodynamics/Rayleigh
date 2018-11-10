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

Module Timing
USE RA_MPI_BASE, Only : MPI_WTIME, MPI_WTICK
Type, Public :: Timer
    Real*8 :: delta, elapsed
    Real*8 :: t1

    Contains
    Procedure :: Init  => Initialize_Timer
    Procedure :: startclock
    Procedure :: stopclock
    Procedure :: increment
End Type Timer


Contains

Subroutine get_ticklength(tl)
    Implicit None
    Real*8, Intent(InOUt) :: tl
    tl = mpi_wtick()
End Subroutine get_ticklength


Subroutine Initialize_Timer(self)
    Implicit None
    Class(Timer) :: self
    self%elapsed = 0.0d0
    self%t1 = 0.0d0
    self%delta = 0.0d0
End Subroutine Initialize_Timer

Subroutine Startclock(self)
    Implicit None
    Class(Timer) :: self
    self%t1 = MPI_WTIME()
End Subroutine Startclock

Subroutine Stopclock(self)
    Implicit None
    Real*8 :: t2
    Class(Timer) :: self
    t2 = MPI_WTIME()
    self%delta = t2-self%t1
End Subroutine Stopclock

Subroutine Increment(self)
    Implicit None
    Class(Timer) :: self
    Call self%stopclock()
    self%elapsed = self%elapsed+self%delta
End Subroutine Increment
End Module Timing
