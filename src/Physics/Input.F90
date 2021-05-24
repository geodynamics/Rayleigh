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
                             ncpu_global, aspect_ratio, l_max
    Use Controls,     Only : temporal_controls_namelist, numerical_controls_namelist, &
                             physical_controls_namelist, max_iterations, pad_alltoall, &
                             multi_run_mode, nruns, rundirs, my_path, run_cpus, &
                             io_controls_namelist, new_iteration, jobinfo_file, my_sim_id, &
                             integer_input_digits, int_in_fmt, initialize_io_format_codes, &
                             full_restart
    Use Spherical_IO, Only : output_namelist
    Use BoundaryConditions, Only : boundary_conditions_namelist, T_top_file, T_bottom_file, &
                                   dTdr_top_file, dTdr_bottom_file, C_Top_File, C_Bottom_File
    Use Initial_Conditions, Only : initial_conditions_namelist, alt_check, init_type, &
                                   magnetic_init_type, restart_iter
    Use TestSuite, Only : test_namelist
    Use PDE_Coefficients, Only : reference_namelist, Prandtl_Number, Rayleigh_Number, &
                               Magnetic_Prandtl_Number, Ekman_Number, Transport_Namelist, &
                               custom_reference_file
    Use Parallel_Framework, Only : pfi
    Use Checkpointing, Only : auto_fmt

    Implicit None

    Interface Read_CMD_Line
        Module Procedure Read_CMD_Integer, Read_CMD_Double, Read_CMD_Logical
        Module Procedure Read_CMD_String
    End Interface

    Integer, Private :: my_sim_rank
    Integer, Private :: input_error_code=13

Contains

    Subroutine Main_Input()
        Use RA_MPI_Base
        Use MPI_Layer
        Implicit None
        Integer :: nlines, line_len,  ierr, pars(2), isave
        Character*120 :: input_file, iter_string, input_prefix, res_eq_file
        Character(len=:), Allocatable :: input_as_string(:)
        Type(Communicator) :: sim_comm
        Logical :: read_complete
        Integer :: full_restart_iter=0
        input_file = Trim(my_path)//'main_input'

        ! Check command line to see if this is a full restart,
        ! in which case all input comes from checkpoint directory.    
        restart_iter = 0

        full_restart=.false.
        Call Read_CMD_Line('-full_restart' , full_restart_iter)

        If (full_restart_iter .ne. 0) Then
            full_restart = .true.
            If (full_restart_iter .lt. 0) Then
                Write(iter_string,auto_fmt) -full_restart_iter
                input_prefix= Trim(my_path)//'Checkpoints/quicksave_'//TRIM(iter_string)
            Else
                ! We may be restarting from a checkpoint written with 
                ! more/less than the default 8 digits.
                Call Read_CMD_Line('-input_digits',integer_input_digits)
                isave = integer_input_digits
                Call Initialize_IO_Format_Codes()
                Write(iter_string,int_in_fmt) full_restart_iter
                input_prefix= Trim(my_path)//'Checkpoints/'//TRIM(iter_string)
            Endif
            res_eq_file = TRIM(input_prefix)//'/equation_coefficients'
            input_file = TRIM(input_prefix)//'/main_input'
            

        Endif

        ! For multi-run mode, we need a unique communicator for each run before the read/broadcast.
        ! Simulation-specific communicators have yet to be set up at this point in the program flow.
        ! We create some temporary communicators here so that run parameters can be broadcast.
        Call MPI_Comm_Split(MPI_COMM_WORLD, my_sim_id, global_rank, sim_comm%comm, ierr)
        Call MPI_Comm_Rank(sim_comm%comm, sim_comm%rank, ierr)
        my_sim_rank = sim_comm%rank

        ! Rank zero read the file and broadcasts the file size to all other ranks.
        If (my_sim_rank .eq. 0)  Call File_to_String(input_file,input_as_string,nlines,line_len)
        pars(1:2) = (/ nlines, line_len /)
        Call MPI_Bcast(pars, 2, MPI_INTEGER, 0, sim_comm%comm,ierr)

        nlines = pars(1)
        line_len = pars(2)
        If (my_sim_rank .gt. 0)  Allocate(  character(len=line_len) :: input_as_string(nlines) )

        ! Rank 0 broadcasts the file contents
        Call MPI_Bcast(input_as_string, nlines*line_len, MPI_CHARACTER, 0, sim_comm%comm,ierr)       
        
        ! Each rank then reads the namelists from the character array
        read_complete = .false.
        Read(input_as_string, nml=problemsize_namelist        , err=314)
        Read(input_as_string, nml=numerical_controls_namelist , err=314)
        Read(input_as_string, nml=physical_controls_namelist  , err=314)
        Read(input_as_string, nml=temporal_controls_namelist  , err=314)
        Read(input_as_string, nml=io_controls_namelist        , err=314)
        Read(input_as_string, nml=output_namelist             , err=314)
        Read(input_as_string, nml=boundary_conditions_namelist, err=314)
        Read(input_as_string, nml=initial_conditions_namelist , err=314)
        Read(input_as_string, nml=test_namelist               , err=314)
        Read(input_as_string, nml=reference_namelist          , err=314)
        Read(input_as_string, nml=transport_namelist          , err=314)

        read_complete = .true.

314     If (.not. read_complete) Then
            If (my_sim_rank .eq. 0) Then
                Write(6,*)' '
                Write(6,*)' Error reading main_input file as broadcast by rank 0.'
                Write(6,*)' Each rank will now attempt to read main_input independently.'
                Write(6,*)' '
            Endif
            Call Main_Input_All_Read()
        Endif

        ! Check the command line to see if any arguments were passed explicitly
        Call CheckArgs()

        ! Finally, if this is a full restart, make sure we actually restart 
        ! using only data from the desired checkpoint directory.
        ! This entails repointing a few file names and input_flags.
        If (full_restart) Then
            init_type=-1
            magnetic_init_type=-1
            restart_iter = full_restart_iter
            custom_reference_file = res_eq_file
            integer_input_digits = isave
            If (global_rank .eq. 0) Then
                Write(6,*)'Initializing in full-restart mode.'
                Write(6,*)'Input file: ', input_file
                Write(6,*)'Equation coefficients file: ' , custom_reference_file
            Endif
        Endif

        DeAllocate(input_as_string)

    End Subroutine Main_Input

    Subroutine Main_Input_All_Read()
        Use RA_MPI_Base
        Implicit None
        Character*120 :: input_file
        Character*256 :: emsg
        Logical :: read_complete
        Integer :: ierr

        input_file = Trim(my_path)//'main_input'

        ! If we get here, the broadcast failed.
        ! All processes now attempt to read main_input.
        read_complete = .false.
        Open(unit=20, file=input_file, status="old", position="rewind")
        Read(unit=20, nml=problemsize_namelist        , err=315, iomsg=emsg)
        Read(unit=20, nml=numerical_controls_namelist , err=315, iomsg=emsg)
        Read(unit=20, nml=physical_controls_namelist  , err=315, iomsg=emsg)
        Read(unit=20, nml=temporal_controls_namelist  , err=315, iomsg=emsg)
        Read(unit=20, nml=io_controls_namelist        , err=315, iomsg=emsg)
        Read(unit=20, nml=output_namelist             , err=315, iomsg=emsg)
        Read(unit=20, nml=boundary_conditions_namelist, err=315, iomsg=emsg)
        Read(unit=20, nml=initial_conditions_namelist , err=315, iomsg=emsg)
        Read(unit=20, nml=test_namelist               , err=315, iomsg=emsg)
        Read(unit=20, nml=reference_namelist          , err=315, iomsg=emsg)
        Read(unit=20, nml=Transport_Namelist          , err=315, iomsg=emsg)
        Close(20)
        read_complete = .true.

315     If (.not. read_complete) Then
            If (my_sim_rank .eq. 0) Then
                Write(6,*)' '
                Write(6,*)' Error:  Multi-process read of main_input also failed.'
                Write(6,*)' Check the contents of your main_input file. '
                Write(6,*)' Compiler error message:  '
                Write(6,*)'   ', TRIM(emsg)
                Write(6,*)'  '
            Endif

            Call MPI_Finalize(ierr)
            Call Exit(input_error_code)
        Else
            If (my_sim_rank .eq. 0) Then
                Write(6,*)' '
                Write(6,*)'...  Multi-process read of main_input succeeded.'
                Write(6,*)' '
            Endif
        Endif       

    End Subroutine Main_Input_All_Read


    Subroutine File_to_String(filename, lines, nlines,max_len)
        Implicit None
        Character(len=*),Intent(in) :: filename
        Character(len=:), Allocatable, Intent(out) :: lines(:)
        Integer, Intent(Out) :: nlines, max_len
        Integer :: iunit,istat
        Character(2048):: line, line2
        Integer :: com_check, line_len, i

        ! This routine reads filename into the output character array 'lines.'
        ! For each line, it filters out any comments( !-marks and ensuing text).

        ! First pass: Count lines and record longest line length.
        max_len = 0
        Open(NewUnit=iunit,File=filename,Status='OLD' ,IOStat=istat)
        If (istat .eq. 0) Then

            nlines = 0 
            Do
                Read(iunit,'(a)',IOStat=istat) line
                If (istat .ne. 0) Exit

                com_check = index(line,'!') ! Line length ends at comment symbol
                If (com_check .ne. 0) Then
                    line_len = com_check-1
                Else
                    line_len = LEN(TRIM(line))
                Endif
                max_len = MAX(max_len,line_len)
                nlines = nlines + 1
            Enddo
            Close (iunit)

        Endif

        Allocate(character(len=max_len) :: lines(nlines))

        ! Second pass: Read the lines + filter out comments
        Open(NewUnit=iunit,File=filename,Status='OLD' , IOStat=istat)
        If (istat .eq. 0) Then

            Do i = 1, nlines

                Read(iunit,'(a)',iostat=istat) line
                com_check = index(line,'!')
                If (com_check .gt. 1) Then
                    lines(i)(1:com_check-1) = line(1:com_check-1)
                    lines(i)(com_check:max_len)=''  ! Failure to initialize remainder of string can corrupt input.
                Else If (com_check .eq. 1) Then
                    lines(i) = ''
                Else
                    lines(i) = line
                Endif

            Enddo
            Close(iunit)

        Endif        

    End Subroutine File_to_String

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
                        my_sim_id = i
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
