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

Module BufferedOutput
    Implicit None
    Type, Public :: OutputBuffer
        Integer :: current_index =1
        Integer :: nlines = 1
        Integer :: file_unit
        Logical :: file_created = .false.
        Character*120 :: file_name
        Character*120, Allocatable :: lines(:)
        Contains
        Procedure :: init => Initialize_OutputBuffer
        Procedure :: print => Write_to_Buffer
        Procedure :: flush_buffer
        Procedure :: partial_flush
        Procedure :: finalize => Reset_OutputBuffer
    End Type OutputBuffer
    Type(OutputBuffer) :: stdout
Contains

    Subroutine Initialize_OutputBuffer(self,file_unit,line_count,filename)
        Class(OutputBuffer) :: self
        Integer, Intent(In) :: file_unit
        Integer, Intent(In), Optional :: line_count
        Character*120, Intent(In), Optional :: filename
        self%nlines = 1
        If (present(line_count)) Then
            If (line_count .gt. 0) Then
                self%nlines = line_count
            Endif
        Endif
        If (present(filename)) Then
            self%file_name = filename
        Endif

        Allocate(self%lines(1:self%nlines))
        self%file_unit = file_unit

    End Subroutine Initialize_OutputBuffer

    Subroutine Write_to_Buffer(self,msg)
        Class (OutputBuffer) :: self
        Character(LEN=*), Intent(In) :: msg
        self%lines(self%current_index) = msg
        If (self%current_index .eq. self%nlines) Then
            Call self%flush_buffer()
        Else
            self%current_index = self%current_index+1
        Endif
    End Subroutine Write_to_Buffer

    Subroutine flush_buffer(self,istop)
        Class (OutputBuffer) :: self
        Integer :: i, imax
        Integer, Intent(In), optional :: istop
        !Flushes current contents of buffer
        !Resets buffer index to 1

        If (self%file_unit .ne. 6) Then
            If (self%file_created) Then
                Open(unit = self%file_unit, file = self%file_name,action="write", status="OLD", POSITION = "APPEND", &
                    & FORM = 'FORMATTED')
            Else
                Open(unit = self%file_unit, file = self%file_name,action="write", status="REPLACE", FORM = 'FORMATTED')
                self%file_created = .true.
            Endif
        Endif

        if (present(istop)) Then
            imax = istop
        Else
            imax = self%current_index
        Endif
        Do i = 1, imax
            Write(self%file_unit,*)Trim(self%lines(i))
        Enddo

        If (self%file_unit .ne. 6) Close(self%file_unit)

        self%current_index = 1
    End Subroutine flush_buffer

    Subroutine partial_flush(self)
        Class(OutputBuffer) :: self
        Integer :: iremain
        ! Forces a flush of the buffer, even if it isn't full
        iremain = self%current_index-1
        Call self%flush_buffer(istop = iremain)
    End Subroutine partial_flush

    Subroutine Reset_OutputBuffer(self)
        !Flushes remaining contents of buffer and resets to initial state
        Class(OutputBuffer) :: self



        Call self%partial_flush()

        self%current_index =1
        self%nlines = 1
        self%file_created = .false.
        !Call self%flush_buffer()
        If (Allocated(self%lines)) DeAllocate(self%lines)
    End Subroutine Reset_OutputBuffer

End Module BufferedOutput
