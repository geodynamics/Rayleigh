Module Sphere_Driver
	Use ClockInfo
	Use Sphere_Hybrid_Space,   Only : rlm_spacea, rlm_spaceb, hybrid_init
	Use Sphere_Physical_Space, Only : physical_space, ohmic_heating_coeff
	Use Sphere_Spectral_Space, Only : post_solve, advancetime, ctemp
    Use Diagnostics_Interface, Only : Reboot_Diagnostics
    Use Spherical_IO, Only : time_to_output
	Use Checkpointing
	Use Controls
	Use Timers
	Use Fields
Contains

	Subroutine Initialize_TimeStepping(iter)
		Implicit None
		Integer, Intent(In) :: iter
		Real*8 :: dr, tdiff

		dr = radius(1)-radius(2)  ! assume uniform grid for now
		tdiff = dr*dr				  ! Viscous diffusion time across one grid point		
		! This may seem stringent, but it is in-line with our max time step in ASH
		!max_time_step = tdiff *10.0d0 
		!max_time_step = max_time_step*4.0d0
		!min_time_step = max_time_step*1.0d-4

		If (iter .eq. 1) Then
			new_deltat   = max_time_step
			deltat   = 0.0d0
			old_deltat   = 0.0d0
			simulation_time = 0.0d0
		Else
			! We have restarted from a checkpoint
			! Change new_deltat and deltat appropriately
			new_deltat = checkpoint_newdt
			deltat = checkpoint_dt
			old_deltat = 0.0d0
            simulation_time = checkpoint_time
		Endif
		new_timestep = .true.

	End Subroutine Initialize_TimeStepping

	Subroutine Main_Loop_Sphere()
		Implicit None
		Integer ::  last_iteration, first_iteration,i
		Real*8  :: captured_time, max_time_seconds	
        Character*14 :: tmstr
        Character*8 :: istr, dtfmt ='(ES10.4)'
        Character*7 :: fmtstr = '(F14.4)', ifmtstr = '(i8.8)' 
		! We enter the main loop assuming that the solve has just been performed
		! and that the equation set structure's RHS contains our primary fields with 
		! radial dimension in-processor.
		! Care needs to be taken at init to ensure fields (W,Z,P,T) are stored
		! in the RHS (they are copied out upon entry into the loop).



		first_iteration = 1+checkpoint_iter ! checkpoint_iter is 0 by default
        If (first_iteration .eq. 1) Then
            Euler_step = .true.
        Endif
        If ((new_iteration .gt. 0) .and. (checkpoint_iter .gt. 0) ) Then
            ! Reset the time step # provided by the checkpoint file.
            first_iteration = new_iteration  
            If (my_rank .eq. 0) Then

                Call stdout%print(' ')
                Call stdout%print('///////////////////////////////////////////////////////////////')
                Call stdout%print(' WARNING:  Time-step counter has been manually reset.')
                Write(istr,ifmtstr)checkpoint_iter
                Call stdout%print('           Checkpoint time-step ID: '//istr//'.')
                Write(istr,ifmtstr)first_iteration
                Call stdout%print('                  New time-step ID: '//istr//'.')
                Call stdout%print('           Revise main_input before the next restart!')
                Call stdout%print('///////////////////////////////////////////////////////////////')
                Call stdout%print(' ')
            Endif
        Endif
		last_iteration = first_iteration + max_iterations-1
		Call Initialize_TimeStepping(first_iteration)
		If ((chebyshev .or. magnetism)) Then
			! work structure for post_solve_cheby
			Call ctemp%init(field_count = wsfcount, config = 'p1b')
		Endif


		Call Hybrid_Init()
		Call StopWatch(loop_time)%StartClock()
        max_time_seconds = 60*max_time_minutes

        iteration = first_iteration


        !//////////////   BEGIN MAIN LOOP
        Do while (iteration .le. last_iteration)    
            !Check here to see if this is an output iteration.
            !If so, we will want to transfer additional information within
            !The transpose buffers
            output_iteration = time_to_output(iteration)

			Call Post_Solve() ! Linear Solve Configuration


			If (my_rank .eq. 0) Then
                Write(istr,ifmtstr)iteration
                Write(tmstr,dtfmt)deltat
                Call stdout%print(' On iteration : '//istr//'    DeltaT :   '//tmstr)
            Endif
			Call rlm_spacea()

			Call Physical_Space()

			Call rlm_spaceb()

			Call AdvanceTime()

            
            ! Disabling this for the time being.  It needs to be brought up-to-date with the
            ! rest of the code.
            !Call Reboot_Diagnostics(iteration)
            !Note:  The above only reboots the diagnostics if: 
            !       1.) It's time to check for the reboot file (based on diagnostic_reboot_interval)
            !       2.) The appropriately named reboot file exists (reboot_diagnostics0, reboot_diagnostics1, etc.)


            Call StopWatch(walltime)%increment() ! Keep track of the walltime
            Call StopWatch(walltime)%startclock()


            !////////////////////////////////////////////////////////////////////
            !   The final part of the loop just deals with cleaning up if it's 
            !      time to end the run.
            global_msgs(2) = stopwatch(walltime)%elapsed !/timer_ticklength
            If (global_msgs(2) .gt. max_time_seconds) Then
                If (my_rank .eq. 0) Then
                    Call stdout%print(' User-specified maximum walltime exceeded.  Cleaning up.')
                Endif
                last_iteration = iteration !force loop to end
            Endif

            If(iteration .eq. last_iteration) Then
                If (save_last_timestep) Then
                    checkpoint_interval = iteration ! force a checkpoint on final iteration
                Endif
            Endif

            Call IsItTimeForACheckpoint(iteration)
            If (ItIsTimeForACheckpoint) Then
                Call StopWatch(cwrite_time)%StartClock()
                If (chk_type .ne. 2) Then
                    Call Write_Checkpoint(wsp%p1b,iteration, deltat,new_deltat,simulation_time)                    
					 
                Else
                    Call Write_Checkpoint_Alt(wsp%p1b,iteration, deltat,new_deltat,simulation_time)

                Endif
				Call StopWatch(cwrite_time)%Increment()
            Endif

            iteration = iteration+1
		Enddo
        !///////////////////// END MAIN LOOP

        Call StopWatch(walltime)%increment()
		Call StopWatch(loop_time)%Increment()
		if (my_rank .eq. 0) Then
            Call stdout%print('  ')
            Call stdout%print('  ')
			Call stdout%print('//////////////////////////////////////////////')
			Call stdout%print('   Measured Timings for Process 0  (seconds)  ')
            Call stdout%print('  ')
            Write(tmstr,fmtstr)StopWatch(loop_time)%elapsed
            Call stdout%print(' Elapsed time: '//tmstr )


            Write(tmstr,fmtstr)StopWatch(ctranspose_time)%elapsed
            Call stdout%print('  Column time: '//tmstr)



            Write(tmstr,fmtstr)StopWatch(rtranspose_time)%elapsed
            Call stdout%print('     Row time: '//tmstr)

            Write(tmstr,fmtstr)StopWatch(legendre_time)%elapsed
            Call stdout%print('Legendre time: '//tmstr)
            Write(tmstr,fmtstr)StopWatch(fft_time)%elapsed
            Call stdout%print('     FFT time: '//tmstr)
            Write(tmstr,fmtstr)StopWatch(solve_time)%elapsed
            Call stdout%print('   Solve time: '//tmstr)
            Write(tmstr,fmtstr)StopWatch(rlma_time)%elapsed
            Call stdout%print('    rlma time: '//tmstr)
            Write(tmstr,fmtstr)StopWatch(rlmb_time)%elapsed
            Call stdout%print('    rlmb time: '//tmstr)
            Write(tmstr,fmtstr)StopWatch(pspace_time)%elapsed
            Call stdout%print('  pspace time: '//tmstr)
            Write(tmstr,fmtstr)StopWatch(psolve_time)%elapsed
            Call stdout%print('  psolve time: '//tmstr)
            Write(tmstr,fmtstr)StopWatch(dphi_time)%elapsed
            Call stdout%print('    dphi time: '//tmstr)

			captured_time = 0.0d0
			Do i = 2, 11
				captured_time = captured_time + StopWatch(i)%elapsed
			Enddo
            Write(tmstr,fmtstr)captured_time
            Call stdout%print('captured time: '//tmstr)

            Write(tmstr,fmtstr)(last_iteration-first_iteration)/StopWatch(loop_time)%elapsed
            Call stdout%print('   ')
            Call stdout%print('     iter/sec: '//tmstr)

			Call stdout%print('//////////////////////////////////////////////')

		Endif
		Call Finalize_Timing(n_r,l_max,max_iterations)
	End Subroutine Main_Loop_Sphere

End Module Sphere_Driver
