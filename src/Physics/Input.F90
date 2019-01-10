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

Module Input
    Use ProblemSize,  Only : problemsize_namelist, nprow, npcol, n_r,n_theta, npout, global_rank, &
                             & ncpu_global
    Use Controls,     Only : temporal_controls_namelist, numerical_controls_namelist, &
                            & physical_controls_namelist, max_iterations, pad_alltoall, &
                            & multi_run_mode, nruns, rundirs, my_path, run_cpus, &
                            & io_controls_namelist, new_iteration, jobinfo_file
    Use Spherical_IO, Only : output_namelist
    Use BoundaryConditions, Only : boundary_conditions_namelist
    Use Initial_Conditions, Only : initial_conditions_namelist, alt_check, init_type, magnetic_init_type
    Use TestSuite, Only : test_namelist
    Use ReferenceState, Only : reference_namelist, Prandtl_Number, Rayleigh_Number, &
                               Magnetic_Prandtl_Number, Ekman_Number
    Use TransportCoefficients, Only : Transport_Namelist
    Use Parallel_Framework, Only : pfi
    Use Stable_Plugin, Only : stable_namelist
    Use Run_Parameters, Only : write_run_parameters
    Implicit None

Contains

    Subroutine Main_Input()
        Implicit None
        Character*120 :: input_file
        input_file = Trim(my_path)//'main_input'

        ! First read the main input file
        Open(unit=20, file=input_file, status="old", position="rewind")
        Read(unit=20, nml=problemsize_namelist)
        Read(unit=20, nml=numerical_controls_namelist)
        Read(unit=20, nml=physical_controls_namelist)
        Read(unit=20, nml=temporal_controls_namelist)
        Read(unit=20, nml=io_controls_namelist)
        Read(unit=20, nml=output_namelist)
        Read(unit=20, nml=boundary_conditions_namelist)
        Read(unit=20, nml=initial_conditions_namelist)
        Read(unit=20, nml=test_namelist)
        Read(unit=20, nml=reference_namelist)
        Read(unit=20, nml=Transport_Namelist)
        !Read(unit=20, nml=Stable_Namelist)
        Close(20)

        ! Check the command line to see if any arguments were passed explicitly
        Call CheckArgs()

        ! write input parameters and other build information
        Call Write_Run_Parameters()

    End Subroutine Main_Input

    Subroutine Check_Run_Mode()
        ! Checks the command line for multiple run flag
        Implicit None
        Character*10 :: arg, arg2
        Integer :: i, itemp
        Logical :: eof_err
        Character*120 :: aline, ifile
        Integer :: errcheck,  rcount, cpu_count, cmin, mx_rank, mn_rank



        ncpu_global = pfi%wcomm%np
        global_rank = pfi%wcomm%rank

        i = 1
        nruns = 0
        DO
            CALL get_command_argument(i, arg)
            IF (LEN_TRIM(arg) == 0) EXIT
            arg2 = TRIM(AdjustL(arg))
            If (arg2 .eq. '-nruns') then
                multi_run_mode = .true.
                CALL get_command_argument(i+1, arg)
                arg2 = TRIM(AdjustL(arg))
                Read (arg2,*) nruns
            Endif
            i = i+1
        ENDDO



        If (multi_run_mode) Then
            If (nruns .gt. 0) Then
                if (global_rank .eq. 0) Write(6,*)'    Multi-run mode enabled.'
                Allocate(rundirs(1:nruns))

                Open(unit=20,file='run_list',form='formatted', status='old',access='stream',iostat = errcheck)

                If (errcheck .eq. 0) Then
                    rcount = 0
                    Do i = 1, nruns
                        Read(20,'(a)',end=314) aline
                        rundirs(i) = trim(adjustl(aline))
                        rcount = rcount+1
                    Enddo
                    eof_err = .false.
                    GOTO 315
        314         eof_err = .true.

        315         If (eof_err) Then
                        If (global_rank .eq. 0) Then
                            Write(6,*)'Multi-run error:  Dirs expected, found: ',nruns, rcount
                Write(6,*)'dirs found: '
                Write(6,*)rundirs
                        Endif
                        multi_run_mode = .false.
                    Endif

                    Close(20)
                Else
                    if (global_rank .eq. 0) Write(6,*)'Multi-run error: cannot open runfile!'
                    multi_run_mode = .false.
                Endif
            Endif
        Endif

        If (multi_run_mode) Then
            Allocate(run_cpus(1:nruns))
            If (global_rank .eq. 0) Then
                cpu_count = 0
                Do i = 1, nruns
                    ifile = rundirs(i)
                    ifile = Trim(ifile)//'/main_input'

                    Open(unit=20, file=ifile, status="old", position="rewind")
                    Read(unit=20, nml=problemsize_namelist)
                    close(20)
                    run_cpus(i) = nprow*npcol

                Enddo
                cpu_count = Sum(run_cpus)
                If (cpu_count .ne. ncpu_global) Then
                    Write(6,*)'////////////////////////////////////////////////////////////////'
                    Write(6,*)'//    Error! Total required CPU count does not equal available cpus.'
                    Write(6,*)'//    Available CPUs (from mpiexec): ', ncpu_global
                    Write(6,*)'//    Required CPUs  (for multiple runs): ', cpu_count
                    Write(6,*)'//    Run will now terminate!'
                    Write(6,*)'////////////////////////////////////////////////////////////////'
                    !run_cpus(:) = -1  ! Will use this as a way of informing other processors of error
                Endif
            Endif
            Call pfi%broadcast_intarr(run_cpus)
            cmin = minval(run_cpus)
            If (cmin .gt. 0) Then
                mn_rank = 0
                Do i = 1, nruns
                    mx_rank = run_cpus(i)-1 +mn_rank
                    If ( (global_rank .le. mx_rank) .and. (global_rank .ge. mn_rank) ) Then
                        my_path = rundirs(i)
                        my_path = Trim(my_path)//'/'
                    Endif
                    mn_rank = mx_rank+1
                Enddo
            Endif


        Endif

    End Subroutine Check_Run_Mode

    Subroutine CheckArgs()
            ! Checks the command line for acceptable arguments.
            ! Specified values overwrite namelist inputs.
            Implicit None
            Character*120 :: arg, arg2
            Integer :: i, itemp
            i = 1
            DO
              CALL get_command_argument(i, arg)
             IF (LEN_TRIM(arg) == 0) EXIT

                arg2 = TRIM(AdjustL(arg))
                If (arg .eq. '-new_iter') then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                    Read (arg2,*) new_iteration
                Endif

                If (arg .eq. '-nprow') then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) nprow
                Endif
                If (arg .eq. '-npcol') Then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) npcol
                Endif
                If (arg .eq. '-npout') Then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) npout
                Endif
                If (arg .eq. '-niter') Then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) max_iterations
                Endif
                If (arg .eq. '-nr') Then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) n_r
                Endif
                If (arg .eq. '-ntheta') Then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) n_theta
                Endif
                If (arg .eq. '-pata') Then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) itemp
                    if (itemp .eq. 1) then
                        pad_alltoall = .true.
                    else
                        pad_alltoall = .false.
                    endif
                Endif
                If (arg .eq. '-altc') Then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) itemp
                    if (itemp .eq. 1) then
                        alt_check = .true.
                    else
                        alt_check = .false.
                    endif
                Endif
                If (arg .eq. '-jobinfo') then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) jobinfo_file
                Endif

                ! Physical control parameters
                If (arg .eq. '-Ra') then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) Rayleigh_Number
                Endif

                If (arg .eq. '-E') then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) Ekman_Number
                Endif

                If (arg .eq. '-Pr') then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) Prandtl_Number
                Endif

                If (arg .eq. '-Pm') then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) Magnetic_Prandtl_Number
                Endif

                If (arg .eq. '-init') then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) init_type
                Endif

                If (arg .eq. '-mag-init') then
                    CALL get_command_argument(i+1, arg)
                    arg2 = TRIM(AdjustL(arg))
                  Read (arg2,*) magnetic_init_type
                Endif

              i = i+1

          END DO
    End Subroutine CheckArgs

End Module Input
