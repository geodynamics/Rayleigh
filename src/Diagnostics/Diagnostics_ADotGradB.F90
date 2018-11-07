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

#include "indices.F"

Module Diagnostics_ADotGradB
    Use Diagnostics_Base

    Interface ADotGradB
        Module Procedure ADotGradB_3D3D, ADotGradB_3D2D
        Module Procedure ADotGradB_2D3D, ADotGradB_2D2D
    End Interface

Contains
    Subroutine ADotGradB_3D3D(abuff,bbuff,cbuff,aindices, bindices, cindices)
        Implicit None
        !Computes C = A dot grad B in spherical coordinates
        !             A,B, and C are assumed to be functions of r,theta, and phi
        !
        !Inputs:
        !   abuff -- 4-D array containing r,theta,phi components of A on r,theta,phi grid
        !   bbuff -- 4-D array containing r,theta,phi, components of B and their derivatives
        !               on r,theta,phi grid
        !
        !Optional Inputs:
        !   aindices,bindices, cindices --  These are optional arrays describing where the
        !                                   quantities above are located within their respective
        !                                   buffers.  If these arrays are provided, quantities
        !                                   are assumed to reside at default locations within
        !                                   each buffer (see below)
        !
        !Outputs
        !   cbuff -- buffer containing r,theta,phi components of C on an r,theta,phi grid

        Real*8, Intent(InOut) :: abuff(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Intent(InOut) :: bbuff(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Intent(InOut) :: cbuff(1:,my_r%min:,my_theta%min:,1:)
        Integer, Intent(In), Optional :: aindices(1:), bindices(1:), cindices(1:)

        Integer :: r,k, t   !indices for iteration

        Integer :: iar    , iatheta , iaphi    ! Abuff indices for a_r, a_theta, and a_phi
        Integer :: ibr    , ibtheta , ibphi    ! Bbuff indices for b_r, b_theta, and b_phi
        Integer :: idbrdr , idbtdr  , idbpdr   ! Bbuff indices for d[b_r]/dr     , d[b_theta]/dr     , d[b_phi]/dr
        Integer :: idbrdt , idbtdt  , idbpdt   ! Bbuff indices for d[b_r]/dtheta , d[b_theta]/dtheta , d[b_phi]/dtheta
        Integer :: idbrdp , idbtdp  , idbpdp   ! Bbuff indices for d[b_r]/dphi   , d[b_theta]/dphi   , d[b_phi]/dphi
        Integer :: icr    , ictheta , icphi    ! Cbuff indices for c_r, c_theta, and c_phi



        If (present(aindices)) Then
            iar     = aindices(1)
            iatheta = aindices(2)
            iaphi   = aindices(3)
        Else
            iar     = 1
            iatheta = 2
            iaphi   = 3
        Endif

        If (present(bindices)) Then
            ibr     = bindices(1)
            ibtheta = bindices(2)
            ibphi   = bindices(3)

            idbrdr  = bindices(4)
            idbrdt  = bindices(5)
            idbrdp  = bindices(6)

            idbtdr  = bindices(7)
            idbtdt  = bindices(8)
            idbtdp  = bindices(9)

            idbpdr  = bindices(10)
            idbpdt  = bindices(11)
            idbpdp  = bindices(12)
        Else
            ibr     = 1
            ibtheta = 2
            ibphi   = 3

            idbrdr  = 4
            idbrdt  = 5
            idbrdp  = 6

            idbtdr  = 7
            idbtdt  = 8
            idbtdp  = 9

            idbpdr  = 10
            idbpdt  = 11
            idbpdp  = 12
        Endif

        If (present(cindices)) Then
            icr     = cindices(1)
            ictheta = cindices(2)
            icphi   = cindices(3)
        Else
            icr     = 1
            ictheta = 2
            icphi   = 3
        Endif

        !//////////////////////////////////////////////////////////////////////////////////

        !--- Radial Component

        DO_PSI
            cbuff(PSI,icr) = abuff(PSI,iar) * bbuff(PSI,idbrdr)                    &

                & + ( abuff(PSI,iatheta) * ( bbuff(PSI,idbrdt)-bbuff(PSI,ibtheta)) &

                & +   abuff(PSI,iaphi)   * ( bbuff(PSI,idbrdp)*csctheta(t)         &
                &                           -bbuff(PSI,ibphi) ) )* one_over_r(r)
        END_DO


        !--- Theta Component

        DO_PSI
            cbuff(PSI,ictheta) = abuff(PSI,iar  ) *  bbuff(PSI,idbtdr)          &

                & + ( abuff(PSI,iatheta) * (bbuff(PSI,idbtdt) + bbuff(PSI,ibr)) &

                & +   abuff(PSI,iaphi)   * ( bbuff(PSI,idbtdp)*csctheta(t)      &
                                            -bbuff(PSI,ibphi )*cottheta(t) ) )  &
                & *  one_over_r(r)
        END_DO

        !--- Phi Component

        DO_PSI
            cbuff(PSI,icphi) =  abuff(PSI,iar  ) *  bbuff(PSI,idbpdr)                       &

                & + ( abuff(PSI,iatheta)*bbuff(PSI,idbpdt)  &

                & +  abuff(PSI,iaphi)*(bbuff(PSI,idbpdp)*csctheta(t) + bbuff(PSI,ibr) &
                & +  bbuff(PSI,ibtheta)*cottheta(t) ) ) &
                & *  one_over_r(r)
        END_DO

    End Subroutine ADotGradB_3D3D

    Subroutine ADotGradB_3D2D(abuff,bbuff,cbuff,aindices, bindices, cindices)
        Implicit None
        !Computes C = A dot grad B in spherical coordinates
        !             A and C are assumed to be functions of r,theta, and phi
        !             B is assumed to be a 2-D function of r,theta
        !
        !Inputs:
        !   abuff -- 4-D array containing r,theta,phi components of A on r,theta,phi grid
        !   bbuff -- 3-D array containing r,theta,phi, components of B and their derivatives
        !               on r,theta grid
        !
        !Optional Inputs:
        !   aindices,bindices, cindices --  These are optional arrays describing where the
        !                                   quantities above are located within their respective
        !                                   buffers.  If these arrays are provided, quantities
        !                                   are assumed to reside at default locations within
        !                                   each buffer (see below)
        !
        !Outputs
        !   cbuff -- buffer containing r,theta,phi components of C on an r,theta,phi grid

        Real*8, Intent(InOut) :: abuff(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Intent(InOut) :: bbuff(my_r%min:,my_theta%min:,1:)
        Real*8, Intent(InOut) :: cbuff(1:,my_r%min:,my_theta%min:,1:)
        Integer, Intent(In), Optional :: aindices(1:), bindices(1:), cindices(1:)

        Integer :: r,k, t   !indices for iteration

        Integer :: iar    , iatheta , iaphi    ! Abuff indices for a_r, a_theta, and a_phi
        Integer :: ibr    , ibtheta , ibphi    ! Bbuff indices for b_r, b_theta, and b_phi
        Integer :: idbrdr , idbtdr  , idbpdr   ! Bbuff indices for d[b_r]/dr     , d[b_theta]/dr     , d[b_phi]/dr
        Integer :: idbrdt , idbtdt  , idbpdt   ! Bbuff indices for d[b_r]/dtheta , d[b_theta]/dtheta , d[b_phi]/dtheta
        Integer :: idbrdp , idbtdp  , idbpdp   ! Bbuff indices for d[b_r]/dphi   , d[b_theta]/dphi   , d[b_phi]/dphi
        Integer :: icr    , ictheta , icphi    ! Cbuff indices for c_r, c_theta, and c_phi



        If (present(aindices)) Then
            iar     = aindices(1)
            iatheta = aindices(2)
            iaphi   = aindices(3)
        Else
            iar     = 1
            iatheta = 2
            iaphi   = 3
        Endif

        If (present(bindices)) Then
            ibr     = bindices(1)
            ibtheta = bindices(2)
            ibphi   = bindices(3)

            idbrdr  = bindices(4)
            idbrdt  = bindices(5)
            idbrdp  = bindices(6)

            idbtdr  = bindices(7)
            idbtdt  = bindices(8)
            idbtdp  = bindices(9)

            idbpdr  = bindices(10)
            idbpdt  = bindices(11)
            idbpdp  = bindices(12)
        Else
            ibr     = 1
            ibtheta = 2
            ibphi   = 3

            idbrdr  = 4
            idbrdt  = 5
            idbrdp  = 6

            idbtdr  = 7
            idbtdt  = 8
            idbtdp  = 9

            idbpdr  = 10
            idbpdt  = 11
            idbpdp  = 12
        Endif

        If (present(cindices)) Then
            icr     = cindices(1)
            ictheta = cindices(2)
            icphi   = cindices(3)
        Else
            icr     = 1
            ictheta = 2
            icphi   = 3
        Endif

        !//////////////////////////////////////////////////////////////////////////////////

        !--- Radial Component

        DO_PSI
            cbuff(PSI,icr) = abuff(PSI,iar) * bbuff(PSI2,idbrdr)                    &

                & + ( abuff(PSI,iatheta) * ( bbuff(PSI2,idbrdt)-bbuff(PSI2,ibtheta)) &

                & +   abuff(PSI,iaphi)   * ( bbuff(PSI2,idbrdp)*csctheta(t)         &
                &                           -bbuff(PSI2,ibphi) ) )* one_over_r(r)
        END_DO


        !--- Theta Component

        DO_PSI
            cbuff(PSI,ictheta) = abuff(PSI,iar  ) *  bbuff(PSI2,idbtdr)          &

                & + ( abuff(PSI,iatheta) * (bbuff(PSI2,idbtdt) + bbuff(PSI2,ibr)) &

                & +   abuff(PSI,iaphi)   * ( bbuff(PSI2,idbtdp)*csctheta(t)      &
                                            -bbuff(PSI2,ibphi )*cottheta(t) ) )  &
                & *  one_over_r(r)
        END_DO

        !--- Phi Component

        DO_PSI
            cbuff(PSI,icphi) =  abuff(PSI,iar  ) *  bbuff(PSI2,idbpdr)                       &

                & + (abuff(PSI,iatheta)*bbuff(PSI2,idbpdt)  &

                & +  abuff(PSI,iaphi)  *(bbuff(PSI2,idbpdp)*csctheta(t) + bbuff(PSI2,ibr ) &
                & +  bbuff(PSI2,ibtheta)*cottheta(t)) ) &
                & *  one_over_r(r)
        END_DO

    End Subroutine ADotGradB_3D2D


    Subroutine ADotGradB_2D3D(abuff,bbuff,cbuff,aindices, bindices, cindices)
        Implicit None
        !Computes C = A dot grad B in spherical coordinates
        !             A is assumed to be a function or r,theta
        !             B and C are assumed to be functions of r,theta, and phi
        !
        !Inputs:
        !   abuff -- 3-D array containing r,theta,phi components of A on r,theta grid
        !   bbuff -- 4-D array containing r,theta,phi, components of B and their derivatives
        !               on r,theta,phi grid
        !
        !Optional Inputs:
        !   aindices,bindices, cindices --  These are optional arrays describing where the
        !                                   quantities above are located within their respective
        !                                   buffers.  If these arrays are provided, quantities
        !                                   are assumed to reside at default locations within
        !                                   each buffer (see below)
        !
        !Outputs
        !   cbuff -- buffer containing r,theta,phi components of C on an r,theta,phi grid

        Real*8, Intent(InOut) :: abuff(my_r%min:,my_theta%min:,1:)
        Real*8, Intent(InOut) :: bbuff(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Intent(InOut) :: cbuff(1:,my_r%min:,my_theta%min:,1:)
        Integer, Intent(In), Optional :: aindices(1:), bindices(1:), cindices(1:)

        Integer :: r,k, t   !indices for iteration

        Integer :: iar    , iatheta , iaphi    ! Abuff indices for a_r, a_theta, and a_phi
        Integer :: ibr    , ibtheta , ibphi    ! Bbuff indices for b_r, b_theta, and b_phi
        Integer :: idbrdr , idbtdr  , idbpdr   ! Bbuff indices for d[b_r]/dr     , d[b_theta]/dr     , d[b_phi]/dr
        Integer :: idbrdt , idbtdt  , idbpdt   ! Bbuff indices for d[b_r]/dtheta , d[b_theta]/dtheta , d[b_phi]/dtheta
        Integer :: idbrdp , idbtdp  , idbpdp   ! Bbuff indices for d[b_r]/dphi   , d[b_theta]/dphi   , d[b_phi]/dphi
        Integer :: icr    , ictheta , icphi    ! Cbuff indices for c_r, c_theta, and c_phi



        If (present(aindices)) Then
            iar     = aindices(1)
            iatheta = aindices(2)
            iaphi   = aindices(3)
        Else
            iar     = 1
            iatheta = 2
            iaphi   = 3
        Endif

        If (present(bindices)) Then
            ibr     = bindices(1)
            ibtheta = bindices(2)
            ibphi   = bindices(3)

            idbrdr  = bindices(4)
            idbrdt  = bindices(5)
            idbrdp  = bindices(6)

            idbtdr  = bindices(7)
            idbtdt  = bindices(8)
            idbtdp  = bindices(9)

            idbpdr  = bindices(10)
            idbpdt  = bindices(11)
            idbpdp  = bindices(12)
        Else
            ibr     = 1
            ibtheta = 2
            ibphi   = 3

            idbrdr  = 4
            idbrdt  = 5
            idbrdp  = 6

            idbtdr  = 7
            idbtdt  = 8
            idbtdp  = 9

            idbpdr  = 10
            idbpdt  = 11
            idbpdp  = 12
        Endif

        If (present(cindices)) Then
            icr     = cindices(1)
            ictheta = cindices(2)
            icphi   = cindices(3)
        Else
            icr     = 1
            ictheta = 2
            icphi   = 3
        Endif

        !//////////////////////////////////////////////////////////////////////////////////

        !--- Radial Component

        DO_PSI
            cbuff(PSI,icr) = abuff(PSI2,iar) * bbuff(PSI,idbrdr)                    &

                & + ( abuff(PSI2,iatheta) * ( bbuff(PSI,idbrdt)-bbuff(PSI,ibtheta)) &

                & +   abuff(PSI2,iaphi)   * ( bbuff(PSI,idbrdp)*csctheta(t)         &
                &                           -bbuff(PSI,ibphi) ) )* one_over_r(r)
        END_DO


        !--- Theta Component

        DO_PSI
            cbuff(PSI,ictheta) = abuff(PSI2,iar  ) *  bbuff(PSI,idbtdr)          &

                & + ( abuff(PSI2,iatheta) * (bbuff(PSI,idbtdt) + bbuff(PSI,ibr)) &

                & +   abuff(PSI2,iaphi)   * ( bbuff(PSI,idbtdp)*csctheta(t)      &
                                            -bbuff(PSI,ibphi )*cottheta(t) ) )  &
                & *  one_over_r(r)
        END_DO

        !--- Phi Component

        DO_PSI
            cbuff(PSI,icphi) =  abuff(PSI2,iar  ) *  bbuff(PSI,idbpdr)                       &

                & + (abuff(PSI2,iatheta)*bbuff(PSI,idbpdt)  &

                & +  abuff(PSI2,iaphi)  *(bbuff(PSI,idbpdp)*csctheta(t) + bbuff(PSI,ibr ) &
                & +  bbuff(PSI,ibtheta)*cottheta(t)) ) &
                & *  one_over_r(r)
        END_DO

    End Subroutine ADotGradB_2D3D


    Subroutine ADotGradB_2D2D(abuff,bbuff,cbuff,aindices, bindices, cindices)
        Implicit None
        !Computes C = A dot grad B in spherical coordinates
        !             A and B are assumed to be 2-D functions of r and theta
        !             C is assumed to be a 3-D function of r,theta, and phi
        !
        !Inputs:
        !   abuff -- 3-D array containing r,theta,phi components of A on r,theta grid
        !   bbuff -- 3-D array containing r,theta,phi, components of B and their derivatives
        !               on r,theta grid
        !
        !Optional Inputs:
        !   aindices,bindices, cindices --  These are optional arrays describing where the
        !                                   quantities above are located within their respective
        !                                   buffers.  If these arrays are provided, quantities
        !                                   are assumed to reside at default locations within
        !                                   each buffer (see below)
        !
        !Outputs
        !   cbuff -- buffer containing r,theta,phi components of C on an r,theta,phi grid

        Real*8, Intent(InOut) :: abuff(my_r%min:,my_theta%min:,1:)
        Real*8, Intent(InOut) :: bbuff(my_r%min:,my_theta%min:,1:)
        Real*8, Intent(InOut) :: cbuff(1:,my_r%min:,my_theta%min:,1:)
        Integer, Intent(In), Optional :: aindices(1:), bindices(1:), cindices(1:)

        Integer :: r,k, t   !indices for iteration

        Integer :: iar    , iatheta , iaphi    ! Abuff indices for a_r, a_theta, and a_phi
        Integer :: ibr    , ibtheta , ibphi    ! Bbuff indices for b_r, b_theta, and b_phi
        Integer :: idbrdr , idbtdr  , idbpdr   ! Bbuff indices for d[b_r]/dr     , d[b_theta]/dr     , d[b_phi]/dr
        Integer :: idbrdt , idbtdt  , idbpdt   ! Bbuff indices for d[b_r]/dtheta , d[b_theta]/dtheta , d[b_phi]/dtheta
        Integer :: idbrdp , idbtdp  , idbpdp   ! Bbuff indices for d[b_r]/dphi   , d[b_theta]/dphi   , d[b_phi]/dphi
        Integer :: icr    , ictheta , icphi    ! Cbuff indices for c_r, c_theta, and c_phi



        If (present(aindices)) Then
            iar     = aindices(1)
            iatheta = aindices(2)
            iaphi   = aindices(3)
        Else
            iar     = 1
            iatheta = 2
            iaphi   = 3
        Endif

        If (present(bindices)) Then
            ibr     = bindices(1)
            ibtheta = bindices(2)
            ibphi   = bindices(3)

            idbrdr  = bindices(4)
            idbrdt  = bindices(5)
            idbrdp  = bindices(6)

            idbtdr  = bindices(7)
            idbtdt  = bindices(8)
            idbtdp  = bindices(9)

            idbpdr  = bindices(10)
            idbpdt  = bindices(11)
            idbpdp  = bindices(12)
        Else
            ibr     = 1
            ibtheta = 2
            ibphi   = 3

            idbrdr  = 4
            idbrdt  = 5
            idbrdp  = 6

            idbtdr  = 7
            idbtdt  = 8
            idbtdp  = 9

            idbpdr  = 10
            idbpdt  = 11
            idbpdp  = 12
        Endif

        If (present(cindices)) Then
            icr     = cindices(1)
            ictheta = cindices(2)
            icphi   = cindices(3)
        Else
            icr     = 1
            ictheta = 2
            icphi   = 3
        Endif

        !//////////////////////////////////////////////////////////////////////////////////

        !--- Radial Component

        DO_PSI2
            cbuff(:,PSI2,icr) = abuff(PSI2,iar) * bbuff(PSI2,idbrdr)                    &

                & + ( abuff(PSI2,iatheta) * ( bbuff(PSI2,idbrdt)-bbuff(PSI2,ibtheta)) &

                & +   abuff(PSI2,iaphi)   * ( bbuff(PSI2,idbrdp)*csctheta(t)         &
                &                           -bbuff(PSI2,ibphi) ) )* one_over_r(r)
        END_DO2


        !--- Theta Component

        DO_PSI2
            cbuff(:,PSI2,ictheta) = abuff(PSI2,iar  ) *  bbuff(PSI2,idbtdr)          &

                & + ( abuff(PSI2,iatheta) * (bbuff(PSI2,idbtdt) + bbuff(PSI2,ibr)) &

                & +   abuff(PSI2,iaphi)   * ( bbuff(PSI2,idbtdp)*csctheta(t)      &
                                            -bbuff(PSI2,ibphi )*cottheta(t) ) )  &
                & *  one_over_r(r)
        END_DO2

        !--- Phi Component

        DO_PSI2
            cbuff(:,PSI2,icphi) =  abuff(PSI2,iar  ) *  bbuff(PSI2,idbpdr)                       &

                & + (abuff(PSI2,iatheta)*bbuff(PSI2,idbpdt)   &

                & +  abuff(PSI2,iaphi)  *(bbuff(PSI2,idbpdp)*csctheta(t) + bbuff(PSI2,ibr )  &
                & +  bbuff(PSI2,ibtheta)*cottheta(t) ) ) &
                & *  one_over_r(r)
        END_DO2

    End Subroutine ADotGradB_2D2D

End Module Diagnostics_ADotGradB
