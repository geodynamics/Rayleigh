Module Read_CMD
    Implicit None

    Interface Read_CMD_Line
        Module Procedure Read_CMD_Integer, Read_CMD_Double, Read_CMD_Logical
        Module Procedure Read_CMD_String
    End Interface

Contains

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
        Do i = 1, n-1
            Call get_command_argument(i,argname)
            Call get_command_argument(i+1,argval)
            If (istring .eq. argname) Then
                argshift = TRIM(AdjustL(argval))
                ivar=argshift
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
        Do i = 1, n
            Call get_command_argument(i,argname)
            !Call get_command_argument(i+1,argval)
            If (istring .eq. argname) ivar = .true.
        Enddo
    End Subroutine Read_CMD_Logical

    Subroutine Read_CMD_Integer(istring, ivar)
        Implicit None
        Character(*),Intent(In) :: istring
        Integer, Intent(InOut) :: ivar
        Integer :: i, n
        Character*1024 :: argname, argval, argshift
        n = command_argument_count()
        Do i = 1, n-1
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
        Do i = 1, n-1
            Call get_command_argument(i,argname)
            Call get_command_argument(i+1,argval)
            If (istring .eq. argname) Then
                argshift = TRIM(AdjustL(argval))
                Read(argshift,*) ivar
            Endif
        Enddo
    End Subroutine Read_CMD_Double

End Module Read_CMD
