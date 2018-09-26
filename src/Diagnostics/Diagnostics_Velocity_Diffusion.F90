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
Module Diagnostics_Velocity_Diffusion
    Use Diagnostics_Base
    Implicit None


Contains
    !//////////////////////////////////////////////////////////////////////////////
    ! Note:  These diagnostic quantities were programmed shortly before
    !        Rayleigh's release.  They were programmed with readability, and
    !        not efficiency or vectorization, in mind.

    Subroutine Compute_Velocity_Diffusion(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: mypi, amp,del2u, estress
        Real*8, Allocatable :: mu_visc(:), dmudr(:), ovstheta(:), ovs2theta(:)

        Allocate(ovstheta(1:N_theta), ovs2theta(1:N_theta)) ! 1/sin; 1/sin^2
        ovstheta = 1.0d0/sintheta
        ovs2theta = 1.0d0/sin2theta

        !Compute the dynamic viscosity mu=rho*nu (nu is kinematic viscosity)
        Allocate(mu_visc(1:N_R), dmudr(1:N_R))

        mu_visc = ref%density*nu
        dmudr = mu_visc*(ref%dlnrho+dlnu)


        !////////////////////////////////////////////////////////
        ! Full viscous forces
        ! r-direction
        If (compute_quantity(viscous_force_r)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_r}
                del2u = DDBUFF(PSI,dvrdrdr)+Two_Over_R(r)*buffer(PSI,dvrdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvrdtdt)+cottheta(t)*buffer(PSI,dvrdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u-2.0d0*OneOverRsquared(r)*( &
                        buffer(PSI,vr) + &
                        buffer(PSI,dvtdt)+buffer(PSI,vtheta)*cottheta(t) + &
                        ovstheta(t)*buffer(PSI,dvpdp) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*(buffer(PSI,dvrdr)*ref%dlnrho(r)+ &
                        buffer(PSI,vr)*ref%d2lnrho(r) )



                ! Finally, add the piece due to the gradient of mu
                estress = buffer(PSI,dvrdr)-One_Third*buffer(PSI,vr)*ref%dlnrho(r)

                qty(PSI) = 2.0d0*dmudr(r)*estress + mu_visc(r)*del2u


            END_DO

            Call Add_Quantity(qty)
        Endif

        !Theta-direction
        If (compute_quantity(viscous_force_theta)) Then

            DO_PSI
                ! first, compute all the terms multiplied by mu
                ! Del^2 {u_theta}
                del2u = DDBUFF(PSI,dvtdrdr)+Two_Over_R(r)*buffer(PSI,dvtdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvtdtdt)+cottheta(t)*buffer(PSI,dvtdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{u} }_r
                del2u = del2u +OneOverRSquared(r)*( -2.0d0*buffer(PSI,dvrdt) + &
                        ovs2theta(t)*(   buffer(PSI,vtheta) + &
                        2.0d0*costheta(t)*buffer(PSI,dvpdp) ) )

                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*buffer(PSI,dvrdt)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(buffer(PSI,dvrdt)-buffer(PSI,vtheta) ) +buffer(PSI,dvtdr)

                qty(PSI) = dmudr(r)*estress +mu_visc(r)*del2u
            END_DO

            Call Add_Quantity(qty)
        Endif


        !Phi-direction
        If (compute_quantity(viscous_force_phi)) Then

            DO_PSI
                del2u = DDBUFF(PSI,dvpdrdr)+Two_Over_R(r)*buffer(PSI,dvpdr)
                del2u = del2u+OneOverRSquared(r)*(DDBUFF(PSI,dvpdtdt)+cottheta(t)*buffer(PSI,dvpdt))
                del2u = del2u+OneOverRSquared(r)*DDBUFF(PSI,dvpdpdp)*ovs2theta(t)


                ! Move onto the non del squared terms (compressibility term)
                del2u = del2u-One_Third*One_Over_R(r)*ovstheta(t)*buffer(PSI,dvrdp)*ref%dlnrho(r)

                ! Finally, add the piece due to the gradient of mu
                estress = One_Over_R(r)*(ovstheta(t)*buffer(PSI,dvrdp)-buffer(PSI,vphi) )+ &
                          buffer(PSI,dvpdr)

                qty(PSI) =dmudr(r)*estress + mu_visc(r)*del2u
            END_DO

            Call Add_Quantity(qty)
        Endif


        DeAllocate(mu_visc, dmudr)
        DeAllocate(ovstheta,ovs2theta)
    End Subroutine Compute_Velocity_Diffusion

End Module Diagnostics_Velocity_Diffusion
