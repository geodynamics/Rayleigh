Module Timing
USE MPI_BASE, Only : MPI_WTIME, MPI_WTICK
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
