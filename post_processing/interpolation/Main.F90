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
    !Namelist /ProblemSpec_Namelist/ nr, ncube, nthrd, nfloat, cartesian, perm_dir, &
    !        initial_iteration, final_iteration, iteration_step, quantity

    ! Read from command line
    Call Read_CMD_Line('-N'     , ncube)   
    Call Read_CMD_Line('-i'     , input_file)
    Call Read_CMD_Line('-o'     , output_file)
    Call Read_CMD_Line('-do'     , double_precision_output)
    Call Read_CMD_Line('-g'     , grid_file)  ! for legacy support
    Call Read_CMD_Line('-nthread', nthrd)
    Call Read_CMD_Line('-si'     , single_precision_input)
    Call Read_CMD_Line('-h', display_help)
    Call Read_CMD_Line('--help', display_help)

  End Subroutine Read_Input

End Program Driver
