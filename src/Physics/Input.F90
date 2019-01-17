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
                             & ncpu_global, aspect_ratio, l_max
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

    Implicit None

    Interface Read_CMD_Line
        Module Procedure Read_CMD_Integer, Read_CMD_Double, Read_CMD_Logical
        Module Procedure Read_CMD_String
    End Interface

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
        Close(20)

        ! Check the command line to see if any arguments were passed explicitly
        Call CheckArgs()

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

            Call Read_CMD_Line('-new_iter' , new_iteration)
            Call Read_CMD_Line('-nprow'    , nprow)
            Call Read_CMD_Line('-npcol'    , npcol)

            Call Read_CMD_Line('-npout'  , npout)
            Call Read_CMD_Line('-niter'  , max_iterations)
            Call Read_CMD_Line('-nr'     , n_r)
            Call Read_CMD_Line('-ntheta' , n_theta)
            Call Read_CMD_Line('-aspect' , aspect_ratio)
            Call Read_CMD_Line('-lmax' , l_max)

            Call Read_CMD_Line('-pata' , pad_alltoall)
            Call Read_CMD_Line('-altc' , alt_check)

            Call Read_CMD_Line('-jobinfo' , jobinfo_file)               

            Call Read_CMD_Line('-Ra' , Rayleigh_Number)
            Call Read_CMD_Line('-E'  , Ekman_Number)
            Call Read_CMD_Line('-Pr' , Prandtl_Number)
            Call Read_CMD_Line('-Pm' , Magnetic_Prandtl_Number)

            Call Read_CMD_Line('-init'     , init_type)
            Call Read_CMD_Line('-mag-init' , magnetic_init_type)

    End Subroutine CheckArgs

    !////////////////////////////////////////////////////////////
    ! The following subroutines are members of the Read_CMD_Line
    ! interface. Each routine sets the value of ivar to the value 
    ! following istring if it is specified at the command line.  
    ! The one exception is Read_CMD_Logical, which sets ivar to .true.
    ! if the value following istring is 1.  Ivar is set to .false.
    ! otherwise.

    Subroutine Read_CMD_String(istring, ivar)
        Implicit None
        Character(*),Intent(In) :: istring
        Character(*), Intent(InOut) :: ivar
        Integer :: i, n
        Character*1024 :: argname, argval, argshift
        n = command_argument_count()
        Do i = 1, n,2
            Call get_command_argument(i,argname)
            Call get_command_argument(i+1,argval)
            If (istring .eq. argname) Then
                argshift = TRIM(AdjustL(argval))
                Read(argshift,*) ivar
            Endif
        Enddo
    End Subroutine Read_CMD_String

    Subroutine Read_CMD_Logical(istring, ivar)
        Implicit None
        Character(*),Intent(In) :: istring
        Logical, Intent(InOut) :: ivar
        Integer :: i, n, itemp
        Character*1024 :: argname, argval, argshift
        n = command_argument_count()
        Do i = 1, n,2
            Call get_command_argument(i,argname)
            Call get_command_argument(i+1,argval)
            If (istring .eq. argname) Then
                argshift = TRIM(AdjustL(argval))
                Read(argshift,*) itemp
                If (itemp .eq. 1) Then
                    ivar = .true.
                Else
                    ivar = .false.
                Endif
            Endif
        Enddo
    End Subroutine Read_CMD_Logical

    Subroutine Read_CMD_Integer(istring, ivar)
        Implicit None
        Character(*),Intent(In) :: istring
        Integer, Intent(InOut) :: ivar
        Integer :: i, n
        Character*1024 :: argname, argval, argshift
        n = command_argument_count()
        Do i = 1, n,2
            Call get_command_argument(i,argname)
            Call get_command_argument(i+1,argval)
            If (istring .eq. argname) Then
                argshift = TRIM(AdjustL(argval))
                Read(argshift,*) ivar
            Endif
        Enddo
    End Subroutine Read_CMD_Integer

    Subroutine Read_CMD_Double(istring, ivar)
        Implicit None
        Character(*),Intent(In) :: istring
        Real*8, Intent(InOut) :: ivar
        Integer :: i, n
        Character*1024 :: argname, argval, argshift
        n = command_argument_count()
        Do i = 1, n,2
            Call get_command_argument(i,argname)
            Call get_command_argument(i+1,argval)
            If (istring .eq. argname) Then
                argshift = TRIM(AdjustL(argval))
                Read(argshift,*) ivar
            Endif
        Enddo
    End Subroutine Read_CMD_Double

End Module Input
