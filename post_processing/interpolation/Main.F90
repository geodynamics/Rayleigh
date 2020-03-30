Program Driver

    Use Interpolation
    Use Read_CMD

    Implicit None
    Logical :: display_help = .false.

    Call Initialize()
    If (display_help) Then
        Call Print_Help_Message()
    Else
        Call Interpolate()
        Call Finalize()
    Endif
Contains

    Subroutine Print_Help_Message()
        Implicit None
        Write(6,*)'This is a help message.'
        Write(6,*)'Exiting...'
    End Subroutine Print_Help_Message

    Subroutine Initialize()
        Call Read_Input()
        Call omp_set_num_threads(nthrd)
    End Subroutine Initialize

    Subroutine Finalize()
        Call Finalize_Interp()
        Print*, 'Finished!'
        Stop
    End Subroutine Finalize

    Subroutine Read_Input()
        Implicit None
        Logical :: exit_program = .false.

        ! Read from command line
        Call Read_CMD_Line('-N'     , ncube)   
        Call Read_CMD_Line('-i'     , input_file)
        Call Read_CMD_Line('-o'     , output_file)
        Call Read_CMD_Line('-d'     , double_precision_output)
        Call Read_CMD_Line('-g'     , grid_file)  ! for legacy support
        Call Read_CMD_Line('-nthread', nthrd)
        Call Read_CMD_Line('-s'     , single_precision_input)
        Call Read_CMD_Line('-h', display_help)
        Call Read_CMD_Line('--help', display_help)
        Call Read_CMD_Line('-spherical', to_spherical)

        to_cartesian = .not. to_spherical

        If (ncube .le. 0) Then
            Write(6,*)' Error:  Resolution of data cube must be set using -N X option, with X > 0.'
            exit_program = .true.
        Endif

        If (input_file .eq. 'None') Then
            Write(6,*)' Error: A Rayleigh-format 3-D input file must be specified by using the -i flag as: -i input_filename'
            exit_program = .true.
        Endif

        If (output_file .eq. 'None') Then
            Write(6,*)' Error: An output file must be specified by using -o flag as: -o output_filename'
            exit_program = .true.
        Endif

        If (exit_program) Then
            Write(6,*)'Exiting...'
            Stop
        Endif

        Write(6,*)'nthrd is: ', nthrd
    End Subroutine Read_Input

End Program Driver
