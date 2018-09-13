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
	Character*8 :: t_ofmt = '(ES12.5)'	! For formatted timestep output
    Logical :: output_iteration = .false.
End Module ClockInfo
