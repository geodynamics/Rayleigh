Program Driver

  Use Interpolation

  Implicit None

  Integer :: namelist_unit=20, kk, nk
  Real*8 :: tbeg, tend, tavg
  Character(len=64) :: arg

  Call getarg(1,arg)
  Print*, 'Have read file name part:', arg

  Call Initialize()

  tavg = 0d0
  kk = 0
  nk = (final_iteration-initial_iteration)/iteration_step
  Do iteration=initial_iteration, final_iteration, iteration_step
     tbeg = 0d0
     Call cpu_time(tbeg)
     Call Make_Iteration_Filename()
     Call Interpolate()
     tend = 0d0
     Call cpu_time(tend)
     tend = tend-tbeg
     If (kk .eq. 0) Then
        tavg = tend
     Else
        tavg = 0.5d0*(tavg + tend)
     EndIf

     If (kk .gt. 1) Then
        Print*, 'Time remaining:', Dble(nk-kk)*tavg
     EndIf
     kk = kk + 1
  EndDo

  Call Finalize()

Contains

  Subroutine Initialize()
    Call Read_Input()
    Call omp_set_num_threads(nthrd)
    Call Make_Iteration_Filename()
    Call Read_Grid()
  End Subroutine Initialize

  Subroutine Finalize()
    Call Finalize_Interp()
    Print*, 'Finished!'
    Stop
  End Subroutine Finalize

  Subroutine Read_Input()

    Namelist /ProblemSpec_Namelist/ nr, ncube, nthrd, nfloat, cartesian, perm_dir, &
            initial_iteration, final_iteration, iteration_step, quantity

    !  read in namelists
    Open(namelist_unit,file='input',status='old')
    Read(namelist_unit,nml=ProblemSpec_Namelist)
    Close(namelist_unit)

    iteration = initial_iteration

    quantity = '_'//trim(arg)

  End Subroutine Read_Input

End Program Driver
