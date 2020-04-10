Program Driver

    Use Interpolation
    Use Read_CMD

    Implicit None
    Logical :: display_help = .false.
    Integer :: t1, t2, count_rate, count_max

    Call System_Clock(t1,count_rate,count_max)
    Call Initialize()
    Call Interpolate()
    Call Finalize()

Contains

    Subroutine Print_Help_Message()
        Implicit None
        Write(6,*)''
        Write(6,*)'    Interp3d'
        Write(6,*)' '
        Write(6,*)'    Interpolates Rayleigh Spherical_3D output to a uniform Cartesian grid.'
        Write(6,*)' '
        Write(6,*)'    Calling syntax:'
        Write(6,*)''
        Write(6,*)'      (1)   Scalar mode:'
        Write(6,*)''
        Write(6,*)'            interp3d -i input_file -o output_file -N X'
        Write(6,*)''
        Write(6,*)'      (2)   Vector mode (spherical vector input; Cartesian vector output)'
        Write(6,*)''
        Write(6,*)'            interp3d -ir r_file -it t_file -ip p_file -ox x_file -oy y_file -oz z_file -N X'
        Write(6,*)''
        Write(6,*)'    Required Flags: '
        Write(6,*)'    '
        Write(6,*)"            -N {X}         :  number of gridpoints along each axis of output cube (X must be an integer)"
        Write(6,*)""
        Write(6,*)"      (1)  Scalar mode: "
        Write(6,*)""
        Write(6,*)"            -i {file name} :  specifies the input file"
        Write(6,*)"            -o {file name} :  specifies output file"
        Write(6,*)""
        Write(6,*)"      (2)  Vector mode: "
        Write(6,*)""
        Write(6,*)"            -ir {file name} :  specifies the radial-component input file in vector mode"
        Write(6,*)"            -it {file name} :  specifies the theta-component input file in vector"
        Write(6,*)"            -ip {file name} :  specifies the phi-component input file in vector mode"
        Write(6,*)"            -ox {file name} :  specifies x-componet output file in vector mode"
        Write(6,*)"            -oy {file name} :  specifies y-componet output file in vector mode"
        Write(6,*)"            -oz {file name} :  specifies z-componet output file in vector mode"
        Write(6,*)" "
        Write(6,*)"    Optional Flags: "
        Write(6,*)" "
        Write(6,*)"            -d              :  double-precision output (default is single precision)"
        Write(6,*)"            -g {file name}  :  grid file (for legacy Spherical_3D format)"
        Write(6,*)"            -h , --help     :  display this help message"
        Write(6,*)""
        Write(6,*)"            -nthread X      :  Specifies the number (X) of OpenMP threads to use."
        Write(6,*)"                               By default, the thread count is determined via"
        Write(6,*)"                               the OMP_NUM_THREADS environment variable."
        Write(6,*)""
        Write(6,*)"            -om {file name} :  (Vector mode only)  If specified, the magnitude of "
        Write(6,*)"                               the Cartesian vector generated in vector mode is written "
        Write(6,*)"                               to {file name}."
        Write(6,*)""
        !Write(6,*)"            -s             :  single-precision input (default is double precision)"
        Write(6,*)""
        Write(6,*)"            -rmin X        :  Set input data to zero where r < X."
        Write(6,*)"            -rmax X        :  Set input data to zero where r > X."
        Write(6,*)"            -rpm           :  Remove the phi (longitudinal) mean from input data."
        Write(6,*)"            -rsm           :  Remove the full spherical mean from input data."
        Write(6,*)"            -v             :  Verbose mode (produce status output)."
        Write(6,*)" "
        Stop
    End Subroutine Print_Help_Message

    Subroutine Initialize()
        Call Read_Input()
        If (nthrd .gt. 0) Call omp_set_num_threads(nthrd)
    End Subroutine Initialize

    Subroutine Finalize()
        Call Finalize_Interp()
        If (verbose) Then
            Call System_Clock(t2,count_rate,count_max)
            Write(6,*)'Complete.  Elapsed time (s): ', real(t2-t1)/real(count_rate)
            Write(6,*)''
        Endif
        Stop
    End Subroutine Finalize

    Subroutine Read_Input()
        Implicit None
        Logical :: exit_program = .false.

        ! Read from command line
        Call Read_CMD_Line('-N'     , ncube)   
        Call Read_CMD_Line('-i'     , input_file)
        Call Read_CMD_Line('-ir'    , input_rfile)
        Call Read_CMD_Line('-it'    , input_tfile)
        Call Read_CMD_Line('-ip'    , input_pfile)
        Call Read_CMD_Line('-o'     , output_file)
        Call Read_CMD_Line('-ox'    , output_xfile)
        Call Read_CMD_Line('-oy'    , output_yfile)
        Call Read_CMD_Line('-oz'    , output_zfile)
        Call Read_CMD_Line('-d'     , double_precision_output)
        Call Read_CMD_Line('-g'     , grid_file)  ! for legacy support
        Call Read_CMD_Line('-nthread', nthrd)
        Call Read_CMD_Line('-s'     , single_precision_input)
        Call Read_CMD_Line('-h', display_help)
        Call Read_CMD_Line('--help', display_help)
        Call Read_CMD_Line('-spherical', to_spherical)
        Call Read_CMD_Line('-rsm', sphere_mean)
        Call Read_CMD_Line('-rpm', phi_mean)
        Call Read_CMD_Line('-rmax',rmax_zero)
        Call Read_CMD_Line('-rmin',rmin_zero)
        Call Read_CMD_Line('-v',verbose)
        Call Read_CMD_Line('-om', mag_file)

        If (mag_file .ne. 'None') output_mag=.true.

        If (display_help) Then 
            Call Print_Help_Message()
            exit_program = .true.
        Endif

        to_cartesian = .not. to_spherical

        If (ncube .le. 0) Then
            Write(6,*)' Error:  Resolution of data cube must be set using -N X option, with X > 0.'
            exit_program = .true.
        Endif

        If (input_file .eq. 'None') Then
            If ( (input_rfile .eq. 'None') .and. (input_tfile .eq. 'None') &
                  .and. (input_pfile .eq. 'None') ) Then 

                Write(6,*)' Error: A Rayleigh-format 3-D input file must be specified by using the -i flag as: -i filename.'
                Write(6,*)'        Optionally, specify three input files using the -ir, -it, -ip flags for vector mode.'
            
                exit_program = .true.
            Else
                vector_mode = .true.
            Endif
        Endif

        If (output_file .eq. 'None') Then

            If ( (output_xfile .eq. 'None') .and. (output_yfile .eq. 'None') &
                  .and. (output_zfile .eq. 'None') ) Then 

                Write(6,*)' Error: An output file must be specified using the -o flag as: -o filename.'
                Write(6,*)'        Optionally, specify three output files using the -ox, -oy, -oz flags for vector mode.'
            
                exit_program = .true.
            Else
                vector_mode = .true.
            Endif

        Endif

        If (exit_program) Then
            Write(6,*)'Exiting...'
            Stop
        Endif

    End Subroutine Read_Input

End Program Driver
