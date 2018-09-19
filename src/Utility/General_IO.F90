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

Module General_IO
    Use Parallel_Framework, Only : pfi
    Use General_MPI, Only : BCAST2D, BCAST1D
        ! NOTE:  need to modify filename before calling this full_path = Trim(my_path)//filename

    Implicit None
Contains
    Subroutine Read_Rayleigh_Array(filename,array_data, dims)
        Character*120, Intent(In) :: filename
        Real*8, Intent(InOut) :: array_data(1:,1:)
        Integer :: pi_integer,nrows, ncolumns
        Integer, Optional, Intent(InOut) :: dims(2)
        Integer :: rcheck, ccheck, arrdims(2)

        arrdims = shape(array_data)
        rcheck = arrdims(1)
        ccheck = arrdims(2)

        If (pfi%gcomm%rank .eq. 0) Then
            !Only one processes actually opens the file
            !After that, the contents of the array are broadcast across columns and rows
            Open(unit=15,file=filename,form='unformatted', status='old',access='stream')
            Read(15)pi_integer

            If (pi_integer .ne. 314) Then
                close(15)
                Open(unit=15,file=filename,form='unformatted', status='old', &
                     CONVERT = 'BIG_ENDIAN' , access='stream')
                Read(15)pi_integer
                If (pi_integer .ne. 314) Then
                    Close(15)
                    Open(unit=15,file=filename,form='unformatted', status='old', &
                     CONVERT = 'LITTLE_ENDIAN' , access='stream')
                    Read(15)pi_integer
                Endif
            Endif

            If (pi_integer .eq. 314) Then
                Read(15)nrows
                Read(15)ncolumns
                If ( (nrows*ncolumns) .gt. (rcheck*ccheck) ) Then
                    Write(6,*)'Error:  file contents exceed array capacity'
                    Write(6,*)'Filename       : ', filename
                    Write(6,*)'File ncolumns  : ', ncolumns
                    Write(6,*)'File nrows     : ', nrows
                    Write(6,*)'Array ncolumns : ', ccheck
                    Write(6,*)'Array nrows    : ', rcheck
                    Write(6,*)'Setting array contents to zero.'
                    array_data(:,:) = 0.0d0
                Else
                    Read(15)array_data(1:nrows,1:ncolumns)
                Endif
                Close(15)
                If (present(dims)) Then
                    dims(1) = nrows
                    dims(2) = ncolumns
                Endif
            Else
                Write(6,*)'Could not resolve endianness.  Setting array elements to zero.'
            Endif
        Endif

        If (present(dims)) Then
            If (pfi%rcomm%rank .eq. 0) Then
                ! Broadcast along the column
                Call BCAST1D(dims,grp = pfi%ccomm)
            Endif
            ! And then along rows
            Call BCAST1D(dims,grp = pfi%rcomm)
        Endif


        If (pfi%rcomm%rank .eq. 0) Then
            ! Broadcast along the column
            Call BCAST2D(array_data,grp = pfi%ccomm)
        Endif
        ! And then along rows
        Call BCAST2D(array_data,grp = pfi%rcomm)


    End Subroutine Read_Rayleigh_Array
    Subroutine Write_Rayleigh_Array(arr,filename)
        Implicit None
        Character*120, Optional, Intent(In) :: filename
        Character*120 :: ref_file
        Integer :: i,j,nq,sig = 314, nx
        Real*8, Intent(In) :: arr(1:,1:)
        nx = size(arr,1)
        nq = size(arr,2)
        If (pfi%gcomm%rank .eq. 0) Then
            Open(unit=15,file=filename,form='unformatted', status='replace',access='stream')
            Write(15)sig
            Write(15)nx
            Write(15)nq
            Write(15)((arr(i,j),i=1,nx),j = 1, nq)

            Close(15)
        Endif
    End Subroutine Write_Rayleigh_Array
End Module General_IO
