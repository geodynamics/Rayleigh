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

!/////////////////////////////////////////////////////////////////////////////////////////
!
!
!
!                                      Rayleigh
!
!
!/////////////////////////////////////////////////////////////////////////////////////////
!Once, there was the
Program Main!
    Use MakeDir
    Use Controls
    Use Fields
    Use Initial_Conditions
    Use Parallel_Framework
    Use ProblemSize
    Use Input
    Use Diagnostics_Interface, Only : Initialize_Diagnostics
    Use TestSuite
    Use Checkpointing
    Use Sphere_Linear_Terms
    Use Sphere_Driver, Only : Main_Loop_Sphere
    Use Timers
    Use Fourier_Transform, Only : Initialize_FFTs
    Use Benchmarking, Only : Initialize_Benchmarking, Benchmark_Input_Reset
    Use Run_Parameters, Only : write_run_parameters

    Implicit None

    Call Main_MPI_Init(global_rank)   !Initialize MPI

    Call Check_Run_Mode()   !This needs to be done before ever reading main input (handles multiple runs)


    Call Main_Input()
    Call Benchmark_Input_Reset() ! Sets run parameters to benchmark parameters if benchmark_mode .ge. 0

    If (test_mode) Then
        Call Initialize_Controls()
        Call Set_Math_Constants()
        Call Init_ProblemSize()
        Call Test_Lib()
    Else
        Call Main_Initialization()
        Call Write_Run_Parameters()  ! write input parameters and other build information
        Call Main_Loop_Sphere()
    Endif
    Call Finalization()

Contains
    Subroutine Main_Initialization()
        Implicit None


        Call Initialize_Controls()

        Call Set_Math_Constants()
        Call Init_ProblemSize()

        Call Initialize_Directory_Structure()

        Call Initialize_Benchmarking()

        Call Initialize_FFts()
        Call Initialize_Reference()

        Call Initialize_Field_Structure()
        Call Initialize_Checkpointing()
        Call Initialize_Transport_Coefficients()
        Call Initialize_Boundary_Conditions()


        Call Write_Equation_Coefficients_File()

        Call Initialize_Diagnostics()

        Call Full_Barrier()

        Call Linear_Init()

        Call Initialize_Fields()
       
        Call Finalize_Boundary_Conditions()

        Call StopWatch(init_time)%increment() ! started in Init_Problemsize just after MPI is started up

        If (my_rank .eq. 0) Then
            Call stdout%print(" Initialization Complete.")
            Call stdout%print(" //////////////////////////////////////")
            Call stdout%print(" ")
            Call stdout%print(" Remember to cite Rayleigh in your publications:")
            Call stdout%print(" Check the ACKNOWLEDGE file for more information.")
            Call stdout%print(" ")
            Call stdout%print(" //////////////////////////////////////")
            Call stdout%print(" ")
        Endif
    End Subroutine Main_Initialization

    Subroutine Initialize_Directory_Structure()
        Implicit None
        Integer :: ecode
        Logical :: file_exists
        If (my_rank .eq. 0) Then
            Call Make_Directory(Trim(my_path)//'Checkpoints',ecode)
            Call Make_Directory(Trim(my_path)//'Timings',ecode)
            Call Make_Directory(Trim(my_path)//'Benchmark_Reports',ecode)

            !Delete possibly existing terminate file from last run
            Inquire(file=Trim(my_path)//Trim(terminate_file),exist=file_exists)
            If (file_exists) Then
               Open(unit=25,file=Trim(my_path)//Trim(terminate_file))
               close(unit=25,status='delete')
            Endif
        Endif
    End Subroutine Initialize_Directory_Structure

    Subroutine Finalization()
        If (.not. test_mode) Then
            If (my_rank .eq. 0) Call stdout%finalize()
        Endif
        Call pfi%exit()
    End Subroutine Finalization
End Program Main
