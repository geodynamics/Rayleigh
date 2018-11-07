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
!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_CUSTOM
!   This is where users can implement their own diagnostics
!   Custom diagnostics are most easily formed using the contents of the buffer
!
!   That array is dimensioned as:
!   buffer(1:n_phi+2, my_r%min:my_r%max, my_theta%min:my_theta%max,1:nvariables)
!
!   The extra 2 in the first index is needed for the in-place FFTs.  Care should be taken
!   to only loop over 1 to n_phi.
!
!   Each index along the 4th dimension of buffer corresponds to a different variable.
!   These indices (left) and the variables they correspond to (right) are given below.
!
! Field variables:
!   vr      -- radial velocity
!   vtheta  -- theta velocity
!   vphi    -- phi velocity
!   tvar    -- temperature or entropy
!   pvar    -- pressure
!   zvar    -- l(l+1)*Z/r^2  where Z is the toroidal streamfunction
!
! Radial Derivatives:
!   dvrdr   -- d(v_r)/dr
!   dvtdr   -- d(v_theta)/dr
!   dvpdr   -- d(v_phi)/dr
!   dtdr    -- d(temperature or entropy)/dr
!
!
! Theta Derivatives:
!   dvrdt   -- d(v_r)/dtheta
!   dvtdt   -- d(v_theta)/dtheta
!   dvpdt   -- d(v_phi)/dtheta
!   dtdt    -- (1/r)*d(temperature or entropy)/dtheta (<--- Note 1/r)
!
!
! Phi Derivatives:
!   dvrdp   --  d(v_r)/dphi
!   dvtdp   --  d(v_theta)/dphi
!   dvpdp   --  d(v_phi)/dphi
!   dtdp    --  (1/r)*d(temperature or entropy)/dphi   (<--- Note 1/r)
!
!
! If Magnetism is On, six additional variables are present:
!   br      -- radial magnetic field
!   btheta  -- theta magnetic field
!   bphi    -- phi magnetic field
!   curlbr      -- [Del x B]_r
!   curlbtheta  -- [Del x B]_theta
!   curlbphi    -- [Del x B]_phi
!
!
! If Induction Output is needed for this iteration,
!   the buffer also holds the derivatives of each
!   component of B.
!
! Radial Derivatives:
!   dbrdr   -- d(b_r)/dr
!   dbtdr   -- d(b_theta)/dr
!   dbpdr   -- d(b_phi)/dr
!
!
! Theta Derivatives:
!   dbrdt   -- d(b_r)/dtheta
!   dbtdt   -- d(b_theta)/dtheta
!   dbpdt   -- d(b_phi)/dtheta


! Phi Derivatives:
!   dbrdp   --  d(b_r)/dphi
!   dbtdp   --  d(b_theta)/dphi
!   dbpdp   --  d(b_phi)/dphi
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Custom
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Custom_MHD_Diagnostics(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        !=============================================
        ! Edit Below This Line (you may define your own variables below)
        Real*8 :: mtmp1, mtmp2, mtmp3   ! temporary variables for use as needed


        ! TUTORIAL EXAMPLE 1:
        ! We begin with an example of cross helicity
        ! Note:  qty is defined and allocated elsewhere
        !        it is dimensioned as qty(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max)

        If (compute_quantity(cross_helicity)) Then
            !Note the call to compute_quantity. This must always be done before adding a diagnostic.
            !The function compute_quantity peforms two functions:
            !   1.  It returns True or False based on whether or not it's time to
            !       output this particular diagnostic.
            !   2.  It sets internal flags in the IO module that allow us to add the diagnostic.
            DO_PSI
                qty(PSI) = buffer(PSI,vr)*buffer(PSI,br)
                qty(PSI) = qty(PSI)+buffer(PSI,vtheta)*buffer(PSI,btheta)
                qty(PSI) = qty(PSI)+buffer(PSI,vphi)*buffer(PSI,vphi)
            END_DO
            ! After we have assigned values to qty, we have to send that quantity to the IO
            ! module for processing (slicing, averaging, etc.)
            ! This is done with a call to Add_quantity as below.
            Call Add_Quantity(qty)
        Endif

        !  TUTORIAL EXAMPLE 2:
        !   Quantities in fbuffer are identical to those in buffer except
        !   that their m = 0 component has been subtracted off.  fbuffer
        !   is defined and allocated elsewhere.  We might be concerned with
        !   the turbulent cross helicity:  v' dot B '.  Calculating this is
        !   a simple modification to the code above:
        If (compute_quantity(turb_cross_helicity)) Then
            !Note the call to compute_quantity. This must always be done before adding a diagnostic.
            !The function compute_quantity peforms two functions:
            !   1.  It returns True or False based on whether or not it's time to
            !       output this particular diagnostic.
            !   2.  It sets internal flags in the IO module that allow us to add the diagnostic.
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*fbuffer(PSI,br)
                qty(PSI) = qty(PSI)+fbuffer(PSI,vtheta)*fbuffer(PSI,btheta)
                qty(PSI) = qty(PSI)+fbuffer(PSI,vphi)*fbuffer(PSI,vphi)
            END_DO
            ! After we have assigned values to qty, we have to send that quantity to the IO
            ! module for processing (slicing, averaging, etc.)
            ! This is done with a call to Add_quantity as below.
            Call Add_Quantity(qty)
        Endif


        !NOTE: we have made use of the macros defined at the top of this file.
        !The DO_PSI loop above is shorthand for the following:
        !
        ! Do t = my_theta%min, my_theta%max
        !    Do r = my_r%min, my_r%max
        !        Do k = 1, n_phi
        !            qty(k,r,t) = buffer(k,r,t,vr)*buffer(k,r,t,br)
        !            qty(k,r,t) = qty(k,r,t)+buffer(k,r,t,vtheta)*buffer(k,r,t,btheta)
        !            qty(k,r,t) = qty(k,r,t)+buffer(k,r,t,vphi)*buffer(k,r,t,vphi)
        !        Enddo
        !    Enddo
        ! Enddo
        !
        ! We highly encourage you to use the PSI macros.
        ! Doing so tends to lead to fewer bugs.


        !  Tutorial Exercise 1:
        !    After examining the example above, try adding your own code to compute
        !    the vb_angle diagnostic (the cosine of the angle between v and B)

        ! YOUR CODE GOES HERE


        ! Edit Above This Line
        !=============================================


    End Subroutine Custom_MHD_Diagnostics

    Subroutine Custom_Hydro_Diagnostics(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        !=============================================
        ! Edit Below This Line (you may define your own variables below)
        Real*8 :: htmp1, htmp2, htmp3   ! temporary variables for use if needed

        ! Tutorial Exercise 2:
        !   Uncomment and modify the code below to assign
        !   v dot grad {T or S} to qty and add it to the outputs.
        !
        !   Note that dtdp and dtdt contain {1/r} d{T or S}/dphi
        !   and {1/r} d {T or S}/dtheta respectively (see comments at the top).
        !
        !   We begin by defining the phi-advection piece,
        !   but note that this isn't added to the outputs yet.
        !
        !   You will need to mimic the compute_quantity/add_quantity
        !   logic from custom_mhd_diagnostics to make this work.

        !DO_PSI
        !    qty(PSI) = wsp%p3a(PSI,vphi)*wsp%p3a(PSI,dtdp)*csctheta(t)  ! note
        !END_DO

        If (compute_quantity(ell0_vr)) Then
            Write(6,*)'Ell0_vr', ell0_vr
            DO_PSI
                qty(PSI) = ell0_values(r,vr)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ell0_tvar)) Then
            Write(6,*)'Ell0 tvar', ell0_tvar
            DO_PSI
                qty(PSI) = ell0_values(r,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ell0_dpdr)) Then
            Write(6,*)'Ell0 tvar', ell0_dpdr
            DO_PSI
                qty(PSI) = ell0_values(r,dpdr)
            END_DO
            Call Add_Quantity(qty)
        Endif



        ! Edit Above This Line
        !=============================================

    End Subroutine Custom_Hydro_Diagnostics

End Module Diagnostics_Custom

