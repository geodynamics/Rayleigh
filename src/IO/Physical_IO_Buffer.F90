Module Physical_IO_Buffer
    Use RA_MPI_Base
    Type, Public :: IO_Buffer_Physical
        Type(communicator) :: ocomm ! output communicator for MPI File I/O
        

    End Type IO_Buffer_Physical

Contains

    Subroutine Initialize_Self(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self

    End Subroutine Initialize_Self

    Subroutine Cache_Data(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self

    End Subroutine Cache_Data

    Subroutine Cascade_Down_Row(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self

    End Subroutine Cascade_Down_Row

    Subroutine Write_Data(self)
        Implicit None
        Class(IO_Buffer_Physical) :: self

    End Subroutine Write_Data


End Module Physical_IO_Buffer
