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

Module Fourier_Derivatives

Interface d_by_dphi
    Module Procedure d_by_dphi3D
End Interface

Contains
    !/////////////////////////////////////////////
    !  d_by_dm3D
    !  Computes derivative of arrin along first dimension
    !  Arrin is assumed to be in spectral space in dimension 1
    !  Storage of m values is assumed to be in FFTW's r2c in place format
    !    Note that this has no spatial scale factor (d_by_dphi vs d_by_dx)
    Subroutine d_by_dphi3D(arrin,arrout)
        Implicit None
        Integer :: i, j, k, m, ashape(1:3),ni,nj,nk
        Real*8, Intent(InOut) :: arrin(:,1:,1:)
        Real*8, Intent(InOut) :: arrout(:,1:,1:)
        ashape = shape(arrin)
        ni = ashape(2)
        nj = ashape(3)
        nk = ashape(1)
        Do j = 1, nj
            Do i = 1, ni
                Do k = 1, nk,2
                    m = (k-1)/2
                    arrout(k,i,j) = -m*arrin(k+1,i,j)
                    arrout(k+1,i,j) = m*arrin(k,i,j)
                Enddo
            Enddo
        Enddo
    End Subroutine d_by_dphi3D



End Module Fourier_Derivatives
