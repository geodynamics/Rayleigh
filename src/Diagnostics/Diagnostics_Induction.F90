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
!               DIAGNOSTICS_INDUCTION
!               This module computes del x (vxB), its
!               constituent terms (i.e., B dot grad v),
!               and their Reynolds decomposition
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Induction
    Use Diagnostics_Base
    Use Diagnostics_ADotGradB
    Implicit None
    Logical :: allocate_indr = .false.
    Logical :: allocate_indt = .false.
    Logical :: allocate_indp = .false.
    Logical :: allocate_indw = .false.
    Logical :: compute_shear = .false.
    Logical :: compute_advec = .false.
    Logical :: compute_vmbm_shear = .false.
    Logical :: compute_vmbm_advec = .false.
    Logical :: compute_vmbp_shear = .false.
    Logical :: compute_vmbp_advec = .false.
    Logical :: compute_vpbm_shear = .false.
    Logical :: compute_vpbm_advec = .false.
    Logical :: compute_vpbp_shear = .false.
    Logical :: compute_vpbp_advec = .false.

Contains

    Subroutine Compute_Induction_Terms(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Allocatable :: ind_r(:,:,:), ind_theta(:,:,:), ind_phi(:,:,:)
        Real*8, Allocatable :: ind_work(:,:,:)
        Real*8, Allocatable :: cbuffer(:,:,:,:)
        Integer :: r,k, t
        Allocate(cbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))
        Call Reset_Induction_Flags()

        If (allocate_indr) Then
            Allocate(    ind_r(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif
        If (allocate_indt) Then
            Allocate(ind_theta(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif
        If (allocate_indp) Then
            Allocate(  ind_phi(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif
        If (allocate_indw) Then
            Allocate(  ind_work(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max))
        Endif
        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 1.    Terms resulting full v cross full B.
        !
        !////////////////////////////////////////////////////////////////////////
        !1a.  B dot grad v
        If (compute_shear) Then

            Call ADotGradB(buffer,buffer,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induct_shear_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(ishear_work)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*buffer(PSI,br)+ &
                               cbuffer(PSI,2)*buffer(PSI,btheta)+ &
                               cbuffer(PSI,3)*buffer(PSI,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)
            Endif

            If (compute_quantity(induct_work)) Then
                DO_PSI
                    ind_work(PSI) = cbuffer(PSI,1)*buffer(PSI,br)+&
                                    cbuffer(PSI,2)*buffer(PSI,btheta)+&
                                    cbuffer(PSI,3)*buffer(PSI,bphi)
                END_DO
            Endif
        Endif

        !1b.  -v dot grad B
        If (compute_advec) Then

            Call ADotGradB(buffer,buffer,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induct_advec_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_advec_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_advec_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(iadvec_work)) Then
                DO_PSI
                    qty(PSI) = -cbuffer(PSI,1)*buffer(PSI,br)     &
                               -cbuffer(PSI,2)*buffer(PSI,btheta) &
                               -cbuffer(PSI,3)*buffer(PSI,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)
            Endif
            If (compute_quantity(induct_work)) Then
                DO_PSI
                    ind_work(PSI) = ind_work(PSI)-&
                                    cbuffer(PSI,1)*buffer(PSI,br)     - &
                                    cbuffer(PSI,2)*buffer(PSI,btheta) - &
                                    cbuffer(PSI,3)*buffer(PSI,bphi)
                END_DO
            Endif
        Endif

        !1c.  -B (div dot v)
        ! Take care with the logic here...

        If (compute_quantity(induct_comp_r) .or. compute_quantity(induct_r)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,br)*buffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_r)) Call Add_Quantity(qty)
            If (compute_quantity(induct_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induct_comp_theta) .or. compute_quantity(induct_theta)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,btheta)*buffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induct_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induct_comp_phi) .or. compute_quantity(induct_phi)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,bphi)*buffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induct_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_phi)
            Endif
        Endif
        If (compute_quantity(induct_work)) Then
            DO_PSI
                ind_work(PSI) = ind_work(PSI)+buffer(PSI,vr)*ref%dlnrho(r)*&
                                ( buffer(PSI,br)**2     + &
                                  buffer(PSI,btheta)**2 + &
                                  buffer(PSI,bphi)**2 )
            END_DO
            Call Add_Quantity(ind_work)
        Endif
        If (compute_quantity(icomp_work)) Then
            DO_PSI
                qty(PSI) = buffer(PSI,vr)*ref%dlnrho(r)*&
                                ( buffer(PSI,br)**2     + &
                                  buffer(PSI,btheta)**2 + &
                                  buffer(PSI,bphi)**2 )
            END_DO
            Call Add_Quantity(qty)
        Endif

        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 2.    Terms resulting from <v> x B'.
        !
        !////////////////////////////////////////////////////////////////////////
        !2a.  B' dot grad <v>
        If (compute_vmbp_shear) Then

            Call ADotGradB(fbuffer,m0_values,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induct_shear_vmbp_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_vmbp_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_vmbp_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ishear_work_pmp)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)+ &
                               cbuffer(PSI,2)*fbuffer(PSI,btheta)+ &
                               cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif



            If (compute_quantity(induct_vmbp_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_vmbp_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_vmbp_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)
            Endif
            If (compute_quantity(induct_work_pmp)) Then
                DO_PSI
                    ind_work(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)+&
                                    cbuffer(PSI,2)*fbuffer(PSI,btheta)+&
                                    cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
            Endif
        Endif

        !2b.  -<v> dot grad B'
        If (compute_vmbp_advec) Then

            Call ADotGradB(m0_values,fbuffer,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induct_advec_vmbp_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_advec_vmbp_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_advec_vmbp_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(iadvec_work_pmp)) Then
                DO_PSI
                    qty(PSI) = -cbuffer(PSI,1)*fbuffer(PSI,br)     &
                               -cbuffer(PSI,2)*fbuffer(PSI,btheta) &
                               -cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_vmbp_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_vmbp_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_vmbp_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)
            Endif
            If (compute_quantity(induct_work_pmp)) Then
                DO_PSI
                    ind_work(PSI) = ind_work(PSI)-&
                                    cbuffer(PSI,1)*fbuffer(PSI,br)     - &
                                    cbuffer(PSI,2)*fbuffer(PSI,btheta) - &
                                    cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
            Endif

        Endif

        !2c.  -B' (div dot <v>)
        ! Take care with the logic here...

        If (compute_quantity(induct_comp_vmbp_r) .or. compute_quantity(induct_vmbp_r)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,br)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vmbp_r)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vmbp_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induct_comp_vmbp_theta) .or. compute_quantity(induct_vmbp_theta)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,btheta)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vmbp_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vmbp_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induct_comp_vmbp_phi) .or. compute_quantity(induct_vmbp_phi)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,bphi)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vmbp_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vmbp_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_phi)
            Endif
        Endif
        If (compute_quantity(induct_work_pmp)) Then
            DO_PSI
                ind_work(PSI) = ind_work(PSI)+m0_values(PSI2,vr)*ref%dlnrho(r)*&
                                ( fbuffer(PSI,br)**2     + &
                                  fbuffer(PSI,btheta)**2 + &
                                  fbuffer(PSI,bphi)**2 )
            END_DO
            Call Add_Quantity(ind_work)
        Endif

        If (compute_quantity(icomp_work_pmp)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)*ref%dlnrho(r)*&
                                ( fbuffer(PSI,br)**2     + &
                                  fbuffer(PSI,btheta)**2 + &
                                  fbuffer(PSI,bphi)**2 )
            END_DO
            Call Add_Quantity(qty)
        Endif
        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 3.    Terms resulting from v' x <B>.
        !
        !////////////////////////////////////////////////////////////////////////
        !3a.  <B> dot grad v'
        If (compute_vpbm_shear) Then

            Call ADotGradB(m0_values,fbuffer,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induct_shear_vpbm_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_vpbm_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_vpbm_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ishear_work_ppm)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)+ &
                               cbuffer(PSI,2)*fbuffer(PSI,btheta)+ &
                               cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif


            If (compute_quantity(induct_vpbm_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_vpbm_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_vpbm_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)
            Endif
            If (compute_quantity(induct_work_ppm)) Then
                DO_PSI
                    ind_work(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)+&
                                    cbuffer(PSI,2)*fbuffer(PSI,btheta)+&
                                    cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
            Endif
        Endif

        !3b.  -v' dot grad <B>
        If (compute_vpbm_advec) Then

            Call ADotGradB(fbuffer,m0_values,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induct_advec_vpbm_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_advec_vpbm_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_advec_vpbm_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(iadvec_work_ppm)) Then
                DO_PSI
                    qty(PSI) = -cbuffer(PSI,1)*fbuffer(PSI,br)     &
                               -cbuffer(PSI,2)*fbuffer(PSI,btheta) &
                               -cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_vpbm_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_vpbm_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_vpbm_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)
            Endif
            If (compute_quantity(induct_work_ppm)) Then
                DO_PSI
                    ind_work(PSI) = ind_work(PSI)-&
                                    cbuffer(PSI,1)*fbuffer(PSI,br)     - &
                                    cbuffer(PSI,2)*fbuffer(PSI,btheta) - &
                                    cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
            Endif




        Endif

        !3c.  -<B> (div dot v')
        ! Take care with the logic here...

        If (compute_quantity(induct_comp_vpbm_r) .or. compute_quantity(induct_vpbm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,br)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vpbm_r)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vpbm_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induct_comp_vpbm_theta) .or. compute_quantity(induct_vpbm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,btheta)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vpbm_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vpbm_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induct_comp_vpbm_phi) .or. compute_quantity(induct_vpbm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,bphi)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vpbm_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vpbm_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_phi)
            Endif
        Endif
        If (compute_quantity(induct_work_ppm)) Then
            DO_PSI
                ind_work(PSI) = ind_work(PSI)+fbuffer(PSI,vr)*ref%dlnrho(r)*&
                                ( fbuffer(PSI,br)*m0_values(PSI2,br)     + &
                                  fbuffer(PSI,btheta)*m0_values(PSI2,btheta) + &
                                  fbuffer(PSI,bphi)*m0_values(PSI2,bphi) )
            END_DO
            Call Add_Quantity(ind_work)
        Endif
        If (compute_quantity(icomp_work_ppm)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%dlnrho(r)*&
                                ( fbuffer(PSI,br)*m0_values(PSI2,br)         + &
                                  fbuffer(PSI,btheta)*m0_values(PSI2,btheta) + &
                                  fbuffer(PSI,bphi)*m0_values(PSI2,bphi) )
            END_DO
            Call Add_Quantity(qty)
        Endif
        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 4.    Terms resulting from <v> x <B>.
        !
        !////////////////////////////////////////////////////////////////////////
        !4a.  <B> dot grad <v>
        If (compute_vmbm_shear) Then

            Call ADotGradB(m0_values,m0_values,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induct_shear_vmbm_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_vmbm_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_vmbm_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(ishear_work_mmm)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*m0_values(PSI2,br)+ &
                               cbuffer(PSI,2)*m0_values(PSI2,btheta)+ &
                               cbuffer(PSI,3)*m0_values(PSI2,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif


            If (compute_quantity(induct_vmbm_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_vmbm_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_vmbm_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)
            Endif

            If (compute_quantity(induct_work_mmm)) Then
                DO_PSI
                    ind_work(PSI) = cbuffer(PSI,1)*m0_values(PSI2,br)+&
                                    cbuffer(PSI,2)*m0_values(PSI2,btheta)+&
                                    cbuffer(PSI,3)*m0_values(PSI2,bphi)
                END_DO
            Endif
        Endif

        !4b.  -<v> dot grad <B>
        If (compute_vmbm_advec) Then

            Call ADotGradB(m0_values,m0_values,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induct_advec_vmbm_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_advec_vmbm_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_advec_vmbm_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif


            If (compute_quantity(iadvec_work_mmm)) Then
                DO_PSI
                    qty(PSI) = -cbuffer(PSI,1)*m0_values(PSI2,br)     &
                               -cbuffer(PSI,2)*m0_values(PSI2,btheta) &
                               -cbuffer(PSI,3)*m0_values(PSI2,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_vmbm_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_vmbm_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_vmbm_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)
            Endif
            If (compute_quantity(induct_work_mmm)) Then
                DO_PSI
                    ind_work(PSI) = ind_work(PSI)-&
                                    cbuffer(PSI,1)*m0_values(PSI2,br)     - &
                                    cbuffer(PSI,2)*m0_values(PSI2,btheta) - &
                                    cbuffer(PSI,3)*m0_values(PSI2,bphi)
                END_DO
            Endif
        Endif

        !4c.  -<B> (div dot <v>)
        ! Take care with the logic here...

        If (compute_quantity(induct_comp_vmbm_r) .or. compute_quantity(induct_vmbm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,br)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vmbm_r)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vmbm_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induct_comp_vmbm_theta) .or. compute_quantity(induct_vmbm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,btheta)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vmbm_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vmbm_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induct_comp_vmbm_phi) .or. compute_quantity(induct_vmbm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,bphi)*m0_values(PSI2,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vmbm_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vmbm_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_phi)
            Endif
        Endif
        If (compute_quantity(induct_work_mmm)) Then
            DO_PSI
                ind_work(PSI) = ind_work(PSI)+m0_values(PSI2,vr)*ref%dlnrho(r)*&
                                ( m0_values(PSI2,br)**2     + &
                                  m0_values(PSI2,btheta)**2 + &
                                  m0_values(PSI2,bphi)**2 )
            END_DO
            Call Add_Quantity(ind_work)
        Endif
        If (compute_quantity(icomp_work_mmm)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,vr)*ref%dlnrho(r)*&
                                ( m0_values(PSI2,br)**2     + &
                                  m0_values(PSI2,btheta)**2 + &
                                  m0_values(PSI2,bphi)**2 )
            END_DO
            Call Add_Quantity(qty)
        Endif
        !////////////////////////////////////////////////////////////////////////
        !
        !   Part 5.    Terms resulting from v' x B'.
        !
        !////////////////////////////////////////////////////////////////////////
        !5a.  B' dot grad v'
        If (compute_vpbp_shear) Then

            Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices = bindex, bindices=vindex)

            If (compute_quantity(induct_shear_vpbp_r)) Then
                qty(:,:,:) = cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_vpbp_theta)) Then
                qty(:,:,:) = cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_shear_vpbp_phi)) Then
                qty(:,:,:) = cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(ishear_work_mpp)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*m0_values(PSI2,br)+ &
                               cbuffer(PSI,2)*m0_values(PSI2,btheta)+ &
                               cbuffer(PSI,3)*m0_values(PSI2,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(ishear_work_ppp)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)+ &
                               cbuffer(PSI,2)*fbuffer(PSI,btheta)+ &
                               cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_vpbp_r)) Then
                 ind_r(:,:,:) = cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_vpbp_theta)) Then
                 ind_theta(:,:,:) = cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_vpbp_phi)) Then
                 ind_phi(:,:,:) = cbuffer(:,:,:,3)
            Endif
            If (compute_quantity(induct_work_ppp)) Then
                DO_PSI
                    ind_work(PSI) = cbuffer(PSI,1)*fbuffer(PSI,br)+&
                                    cbuffer(PSI,2)*fbuffer(PSI,btheta)+&
                                    cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
            Endif
            If (compute_quantity(induct_work_mpp)) Then
                DO_PSI
                    tmp1(PSI) =     cbuffer(PSI,1)*m0_values(PSI2,br)+&
                                    cbuffer(PSI,2)*m0_values(PSI2,btheta)+&
                                    cbuffer(PSI,3)*m0_values(PSI2,bphi)
                END_DO
            Endif

        Endif

        !5b.  -v' dot grad B'
        If (compute_vpbp_advec) Then

            Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices = vindex, bindices=bindex)

            If (compute_quantity(induct_advec_vpbp_r)) Then
                qty(:,:,:) = -cbuffer(:,:,:,1)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_advec_vpbp_theta)) Then
                qty(:,:,:) = -cbuffer(:,:,:,2)
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(induct_advec_vpbp_phi)) Then
                qty(:,:,:) = -cbuffer(:,:,:,3)
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(iadvec_work_mpp)) Then
                DO_PSI
                    qty(PSI) = -cbuffer(PSI,1)*m0_values(PSI2,br)     &
                               -cbuffer(PSI,2)*m0_values(PSI2,btheta) &
                               -cbuffer(PSI,3)*m0_values(PSI2,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(iadvec_work_ppp)) Then
                DO_PSI
                    qty(PSI) = -cbuffer(PSI,1)*fbuffer(PSI,br)     &
                               -cbuffer(PSI,2)*fbuffer(PSI,btheta) &
                               -cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(induct_vpbp_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)-cbuffer(:,:,:,1)
            Endif
            If (compute_quantity(induct_vpbp_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)-cbuffer(:,:,:,2)
            Endif
            If (compute_quantity(induct_vpbp_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)-cbuffer(:,:,:,3)
            Endif
            If (compute_quantity(induct_work_ppp)) Then
                DO_PSI
                    ind_work(PSI) = ind_work(PSI)-&
                                    cbuffer(PSI,1)*fbuffer(PSI,br)     - &
                                    cbuffer(PSI,2)*fbuffer(PSI,btheta) - &
                                    cbuffer(PSI,3)*fbuffer(PSI,bphi)
                END_DO
            Endif
            If (compute_quantity(induct_work_mpp)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)-&
                                    cbuffer(PSI,1)*m0_values(PSI2,br)     - &
                                    cbuffer(PSI,2)*m0_values(PSI2,btheta) - &
                                    cbuffer(PSI,3)*m0_values(PSI2,bphi)
                END_DO
            Endif

        Endif

        !5c.  -B' (div dot v')
        ! Take care with the logic here...

        If (compute_quantity(induct_comp_vpbp_r) .or. compute_quantity(induct_vpbp_r)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,br)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vpbp_r)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vpbp_r)) Then
                 ind_r(:,:,:) = ind_r(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_r)
            Endif
        Endif

        If (compute_quantity(induct_comp_vpbp_theta) .or. compute_quantity(induct_vpbp_theta)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,btheta)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vpbp_theta)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vpbp_theta)) Then
                 ind_theta(:,:,:) = ind_theta(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_theta)
            Endif
        Endif

        If (compute_quantity(induct_comp_vpbp_phi) .or. compute_quantity(induct_vpbp_phi)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,bphi)*fbuffer(PSI,vr)*ref%dlnrho(r)
            END_DO
            If (compute_quantity(induct_comp_vpbp_phi)) Call Add_Quantity(qty)
            If (compute_quantity(induct_vpbp_phi)) Then
                 ind_phi(:,:,:) = ind_phi(:,:,:)+qty(:,:,:)
                 Call Add_Quantity(ind_phi)
            Endif
        Endif
        If (compute_quantity(induct_work_ppp)) Then
            DO_PSI
                ind_work(PSI) = ind_work(PSI)+fbuffer(PSI,vr)*ref%dlnrho(r)*&
                                ( fbuffer(PSI,br    )**2  + &
                                  fbuffer(PSI,btheta)**2  + &
                                  fbuffer(PSI,bphi  )**2 )
            END_DO
            Call Add_Quantity(ind_work)
        Endif

        If (compute_quantity(induct_work_mpp)) Then
            DO_PSI
                tmp1(PSI) = tmp1(PSI)+fbuffer(PSI,vr)*ref%dlnrho(r)*&
                                ( fbuffer(PSI,br    )*m0_values(PSI2,br    ) + &
                                  fbuffer(PSI,btheta)*m0_values(PSI2,btheta) + &
                                  fbuffer(PSI,bphi  )*m0_values(PSI2,bphi  ))
            END_DO
            Call Add_Quantity(tmp1)
        Endif

        If (compute_quantity(icomp_work_mpp)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%dlnrho(r)*&
                                ( fbuffer(PSI,br)*m0_values(PSI2,br)     + &
                                  fbuffer(PSI,btheta)*m0_values(PSI2,btheta) + &
                                  fbuffer(PSI,bphi)*m0_values(PSI2,bphi) )
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(icomp_work_ppp)) Then
            DO_PSI
                qty(PSI) = fbuffer(PSI,vr)*ref%dlnrho(r)*&
                                ( fbuffer(PSI,br)**2     + &
                                  fbuffer(PSI,btheta)**2 + &
                                  fbuffer(PSI,bphi)**2 )
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (allocated(ind_r    )) DeAllocate(ind_r)
        If (allocated(ind_theta)) DeAllocate(ind_theta)
        If (allocated(ind_phi  )) DeAllocate(ind_phi)
        If (allocated(ind_work )) DeAllocate(ind_work)
        DeAllocate(cbuffer)
    End Subroutine Compute_Induction_Terms



    Subroutine Compute_Magnetic_Diffusion(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8 :: del2b
        Real*8, Allocatable :: ovstheta(:), ovs2theta(:)

        Allocate(ovstheta(1:N_theta), ovs2theta(1:N_theta)) ! 1/sin; 1/sin^2
        ovstheta = 1.0d0/sintheta
        ovs2theta = 1.0d0/sin2theta



        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Terms involving B
        ! r-direction
        If (compute_quantity(induct_diff_r) .or. compute_quantity(idiff_work)) Then

            DO_PSI
                ! Del^2 {B_r}
                del2b = DDBUFF(PSI,dbrdrdr)+Two_Over_R(r)*buffer(PSI,dbrdr)
                del2b = del2b+OneOverRSquared(r)*(DDBUFF(PSI,dbrdtdt)+cottheta(t)*buffer(PSI,dbrdt))
                del2b = del2b+OneOverRSquared(r)*DDBUFF(PSI,dbrdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{B} }_r
                del2b = del2b-2.0d0*OneOverRsquared(r)*( &
                        buffer(PSI,br) + &
                        buffer(PSI,dbtdt)+buffer(PSI,btheta)*cottheta(t) + &
                        ovstheta(t)*buffer(PSI,dbpdp) )

                qty(PSI) = eta(r)*del2b


            END_DO
            If (compute_quantity(induct_diff_r)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(idiff_work)) Then
                DO_PSI
                    tmp1(PSI) = qty(PSI)*buffer(PSI,br)
                END_DO
            Endif
        Endif

        !Theta-direction; Full
        If (compute_quantity(induct_diff_theta) .or. compute_quantity(idiff_work)) Then

            DO_PSI

                ! Del^2 {B_theta}
                del2b = DDBUFF(PSI,dbtdrdr)+Two_Over_R(r)*buffer(PSI,dbtdr)
                del2b = del2b+OneOverRSquared(r)*(DDBUFF(PSI,dbtdtdt)+cottheta(t)*buffer(PSI,dbtdt))
                del2b = del2b+OneOverRSquared(r)*DDBUFF(PSI,dbtdpdp)*ovs2theta(t)

                !Add geometric terms to make this { Del^2{B} }_theta
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*buffer(PSI,dbrdt) - &
                        ovs2theta(t)*(   buffer(PSI,btheta) + &
                        2.0d0*costheta(t)*buffer(PSI,dbpdp) ) )

                ! Add the contribution from a gradient in eta
                qty(PSI) = eta(r)*(del2b-buffer(PSI,curlbphi)*dlneta(r))

            END_DO

            If (compute_quantity(induct_diff_theta)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(idiff_work)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*buffer(PSI,btheta)
                END_DO
            Endif

        Endif

        !Phi-direction
        If (compute_quantity(induct_diff_phi) .or. compute_quantity(idiff_work)) Then

            DO_PSI
                ! build Del^2{B_phi}
                del2b = DDBUFF(PSI,dbpdrdr)+Two_Over_R(r)*buffer(PSI,dbpdr)
                del2b = del2b+OneOverRSquared(r)*(DDBUFF(PSI,dbpdtdt)+cottheta(t)*buffer(PSI,dbpdt))
                del2b = del2b+OneOverRSquared(r)*DDBUFF(PSI,dbpdpdp)*ovs2theta(t)


                !Add geometric terms to make this { Del^2{u} }_phi
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*buffer(PSI,dbrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   buffer(PSI,bphi) - &
                        2.0d0*costheta(t)*buffer(PSI,dbtdp) ) )

                qty(PSI) = eta(r)*(del2b+buffer(PSI,curlbtheta)*dlneta(r))

            END_DO


            If (compute_quantity(induct_diff_phi)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(idiff_work)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*buffer(PSI,bphi)
                END_DO
                Call Add_Quantity(tmp1)
            Endif

        Endif


        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        ! Next, terms involving B'
        ! r-direction; fluctuating
        If (compute_quantity(induct_diff_bp_r) .or. compute_quantity(idiff_work_pp)) Then

            DO_PSI

                del2b = d2_fbuffer(PSI,dbrdrdr)+Two_Over_R(r)*fbuffer(PSI,dbrdr)
                del2b = del2b+OneOverRSquared(r)*(d2_fbuffer(PSI,dbrdtdt)+cottheta(t)*fbuffer(PSI,dbrdt))
                del2b = del2b+OneOverRSquared(r)*d2_fbuffer(PSI,dbrdpdp)*ovs2theta(t)


                del2b = del2b-2.0d0*OneOverRsquared(r)*( &
                        fbuffer(PSI,br) + &
                        fbuffer(PSI,dbtdt)+fbuffer(PSI,btheta)*cottheta(t) + &
                        ovstheta(t)*fbuffer(PSI,dbpdp) )


                qty(PSI) = eta(r)*del2b


            END_DO
            If (compute_quantity(induct_diff_bp_r)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(idiff_work_pp)) Then
                DO_PSI
                    tmp1(PSI) = qty(PSI)*buffer(PSI,br)
                END_DO
            Endif
        Endif

        !Theta-direction
        If (compute_quantity(induct_diff_bp_theta) .or. compute_quantity(idiff_work_pp)) Then

            DO_PSI

                del2b = d2_fbuffer(PSI,dbtdrdr)+Two_Over_R(r)*fbuffer(PSI,dbtdr)
                del2b = del2b+OneOverRSquared(r)*(d2_fbuffer(PSI,dbtdtdt)+cottheta(t)*fbuffer(PSI,dbtdt))
                del2b = del2b+OneOverRSquared(r)*d2_fbuffer(PSI,dbtdpdp)*ovs2theta(t)


                del2b = del2b +OneOverRSquared(r)*( 2.0d0*fbuffer(PSI,dbrdt) - &
                        ovs2theta(t)*(   fbuffer(PSI,btheta) + &
                        2.0d0*costheta(t)*fbuffer(PSI,dbpdp) ) )


                qty(PSI) = eta(r)*(del2b-fbuffer(PSI,curlbphi)*dlneta(r))

            END_DO

            If (compute_quantity(induct_diff_bp_theta)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(idiff_work_pp)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*buffer(PSI,btheta)
                END_DO
            Endif
        Endif

        !Phi-direction
        If (compute_quantity(induct_diff_bp_phi) .or. compute_quantity(idiff_work_pp)) Then

            DO_PSI
                del2b = d2_fbuffer(PSI,dbpdrdr)+Two_Over_R(r)*fbuffer(PSI,dbpdr)
                del2b = del2b+OneOverRSquared(r)*(d2_fbuffer(PSI,dbpdtdt)+cottheta(t)*fbuffer(PSI,dbpdt))
                del2b = del2b+OneOverRSquared(r)*d2_fbuffer(PSI,dbpdpdp)*ovs2theta(t)


                !Add geometric terms to make this { Del^2{u} }_phi
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*fbuffer(PSI,dbrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   fbuffer(PSI,bphi) - &
                        2.0d0*costheta(t)*fbuffer(PSI,dbtdp) ) )

                qty(PSI) = eta(r)*(del2b+fbuffer(PSI,curlbtheta)*dlneta(r))
            END_DO


            If (compute_quantity(induct_diff_bp_phi)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(idiff_work_pp)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*buffer(PSI,bphi)
                END_DO
                Call Add_Quantity(tmp1)
            Endif
        Endif

        !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        ! Finally, terms involving <B>
        ! r-direction; mean
        If (compute_quantity(induct_diff_bm_r) .or. compute_quantity(idiff_work_mm)) Then

            DO_PSI

                del2b = d2_m0(PSI2,dbrdrdr)+Two_Over_R(r)*m0_values(PSI2,dbrdr)
                del2b = del2b+OneOverRSquared(r)*(d2_m0(PSI2,dbrdtdt)+cottheta(t)*m0_values(PSI2,dbrdt))
                del2b = del2b+OneOverRSquared(r)*d2_m0(PSI2,dbrdpdp)*ovs2theta(t)


                del2b = del2b-2.0d0*OneOverRsquared(r)*( &
                        m0_values(PSI2,br) + &
                        m0_values(PSI2,dbtdt)+m0_values(PSI2,btheta)*cottheta(t) + &
                        ovstheta(t)*m0_values(PSI2,dbpdp) )

                qty(PSI) = eta(r)*del2b


            END_DO

            If (compute_quantity(induct_diff_bm_r)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(idiff_work_mm)) Then
                DO_PSI
                    tmp1(PSI) = qty(PSI)*m0_values(PSI2,br)
                END_DO
            Endif
        Endif

        !Theta-direction; Mean
        If (compute_quantity(induct_diff_bm_theta) .or. compute_quantity(idiff_work_mm)) Then

            DO_PSI

                del2b = d2_m0(PSI2,dbtdrdr)+Two_Over_R(r)*m0_values(PSI2,dbtdr)
                del2b = del2b+OneOverRSquared(r)*(d2_m0(PSI2,dbtdtdt)+cottheta(t)*m0_values(PSI2,dbtdt))
                del2b = del2b+OneOverRSquared(r)*d2_m0(PSI2,dbtdpdp)*ovs2theta(t)


                del2b = del2b +OneOverRSquared(r)*( 2.0d0*m0_values(PSI2,dbrdt) - &
                        ovs2theta(t)*(   m0_values(PSI2,btheta) + &
                        2.0d0*costheta(t)*m0_values(PSI2,dbpdp) ) )

                qty(PSI) = eta(r)*(del2b-m0_values(PSI2,curlbphi)*dlneta(r))
            END_DO


            If (compute_quantity(induct_diff_bm_theta)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(idiff_work_mm)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*m0_values(PSI2,btheta)
                END_DO
            Endif
        Endif

        !Phi-direction (mean)
        If (compute_quantity(induct_diff_bm_phi) .or. compute_quantity(idiff_work_mm)) Then

            DO_PSI
                del2b = d2_m0(PSI2,dbpdrdr)+Two_Over_R(r)*m0_values(PSI2,dbpdr)
                del2b = del2b+OneOverRSquared(r)*(d2_m0(PSI2,dbpdtdt)+cottheta(t)*m0_values(PSI2,dbpdt))
                del2b = del2b+OneOverRSquared(r)*d2_m0(PSI2,dbpdpdp)*ovs2theta(t)


                !Add geometric terms to make this { Del^2{u} }_phi
                del2b = del2b +OneOverRSquared(r)*( 2.0d0*m0_values(PSI2,dbrdp)*ovstheta(t) - &
                        ovs2theta(t)*(   m0_values(PSI2,bphi) - &
                        2.0d0*costheta(t)*m0_values(PSI2,dbtdp) ) )
                qty(PSI) = eta(r)*(del2b+m0_values(PSI2,curlbtheta)*dlneta(r))



            END_DO
            If (compute_quantity(induct_diff_bm_phi)) Then
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(idiff_work_mm)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)*m0_values(PSI2,bphi)
                END_DO
                Call Add_Quantity(tmp1)
            Endif
        Endif

        DeAllocate(ovstheta,ovs2theta)
    End Subroutine Compute_Magnetic_Diffusion

    Subroutine Reset_Induction_Flags()
        Implicit None
        !/////////////////////////////////////////////////
        !   Compute_Induction_Terms makes use of several logical flags.
        !   Before carrying out the computations in that subroutine,
        !   we set these flags to the appropriate value.

        !  First, reset ALL flags to false
        !  (This routine should be called at each output iteration)
        allocate_indr = .false.
        allocate_indt = .false.
        allocate_indp = .false.
        allocate_indw = .false.

        compute_shear = .false.
        compute_advec = .false.
        compute_vmbm_shear = .false.
        compute_vmbm_advec = .false.
        compute_vmbp_shear = .false.
        compute_vmbp_advec = .false.
        compute_vpbp_shear = .false.
        compute_vpbp_advec = .false.
        compute_vpbm_shear = .false.
        compute_vpbm_advec = .false.

        ! 1.  Decide if we need to allocate the induction arrays
        If (compute_quantity(induct_r     )) allocate_indr = .true.
        If (compute_quantity(induct_vmbm_r)) allocate_indr = .true.
        If (compute_quantity(induct_vmbp_r)) allocate_indr = .true.
        If (compute_quantity(induct_vpbm_r)) allocate_indr = .true.
        If (compute_quantity(induct_vpbp_r)) allocate_indr = .true.

        If (compute_quantity(induct_theta     )) allocate_indt = .true.
        If (compute_quantity(induct_vmbm_theta)) allocate_indt = .true.
        If (compute_quantity(induct_vmbp_theta)) allocate_indt = .true.
        If (compute_quantity(induct_vpbm_theta)) allocate_indt = .true.
        If (compute_quantity(induct_vpbp_theta)) allocate_indt = .true.

        If (compute_quantity(induct_phi     )) allocate_indp = .true.
        If (compute_quantity(induct_vmbm_phi)) allocate_indp = .true.
        If (compute_quantity(induct_vmbp_phi)) allocate_indp = .true.
        If (compute_quantity(induct_vpbm_phi)) allocate_indp = .true.
        If (compute_quantity(induct_vpbp_phi)) allocate_indp = .true.

        If (compute_quantity(induct_work))     allocate_indw = .true.
        If (compute_quantity(induct_work_ppp)) allocate_indw = .true.
        If (compute_quantity(induct_work_mpp)) allocate_indw = .true.
        If (compute_quantity(induct_work_pmp)) allocate_indw = .true.
        If (compute_quantity(induct_work_ppm)) allocate_indw = .true.
        If (compute_quantity(induct_work_mmm)) allocate_indw = .true.

        !2.  Set flags related to full v x B
        If (compute_quantity(induct_work)) Then
            compute_shear = .true.
            compute_advec = .true.
        Endif
        If (compute_quantity(ishear_work)) compute_shear = .true.
        If (compute_quantity(iadvec_work)) compute_advec = .true.

        If (compute_quantity(induct_r)) Then
            compute_shear = .true.
            compute_advec = .true.
        Endif

        If (compute_quantity(induct_shear_r)) compute_shear = .true.
        If (compute_quantity(induct_advec_r)) compute_advec = .true.

        If (compute_quantity(induct_theta)) Then
            compute_shear = .true.
            compute_advec = .true.
        Endif

        If (compute_quantity(induct_shear_theta)) compute_shear = .true.
        If (compute_quantity(induct_advec_theta)) compute_advec = .true.

        If (compute_quantity(induct_phi)) Then
            compute_shear = .true.
            compute_advec = .true.
        Endif

        If (compute_quantity(induct_shear_phi)) compute_shear = .true.
        If (compute_quantity(induct_advec_phi)) compute_advec = .true.

        !3.  Set flags related to <v> x <B>
        If (compute_quantity(induct_work_mmm)) Then
            compute_vmbm_shear = .true.
            compute_vmbm_advec = .true.
        Endif

        If (compute_quantity(ishear_work_mmm)) compute_vmbm_shear = .true.
        If (compute_quantity(iadvec_work_mmm)) compute_vmbm_advec = .true.

        If (compute_quantity(induct_vmbm_r)) Then
            compute_vmbm_shear = .true.
            compute_vmbm_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vmbm_r)) compute_vmbm_shear = .true.
        If (compute_quantity(induct_advec_vmbm_r)) compute_vmbm_advec = .true.

        If (compute_quantity(induct_vmbm_theta)) Then
            compute_vmbm_shear = .true.
            compute_vmbm_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vmbm_theta)) compute_vmbm_shear = .true.
        If (compute_quantity(induct_advec_vmbm_theta)) compute_vmbm_advec = .true.

        If (compute_quantity(induct_vmbm_phi)) Then
            compute_vmbm_shear = .true.
            compute_vmbm_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vmbm_phi)) compute_vmbm_shear = .true.
        If (compute_quantity(induct_advec_vmbm_phi)) compute_vmbm_advec = .true.


        !4. Set flags related to v' x <B>
        If (compute_quantity(induct_work_ppm)) Then
            compute_vpbm_shear = .true.
            compute_vpbm_advec = .true.
        Endif

        If (compute_quantity(ishear_work_ppm)) compute_vpbm_shear = .true.
        If (compute_quantity(iadvec_work_ppm)) compute_vpbm_advec = .true.


        If (compute_quantity(induct_vpbm_r)) Then
            compute_vpbm_shear = .true.
            compute_vpbm_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vpbm_r)) compute_vpbm_shear = .true.
        If (compute_quantity(induct_advec_vpbm_r)) compute_vpbm_advec = .true.

        If (compute_quantity(induct_vpbm_theta)) Then
            compute_vpbm_shear = .true.
            compute_vpbm_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vpbm_theta)) compute_vpbm_shear = .true.
        If (compute_quantity(induct_advec_vpbm_theta)) compute_vpbm_advec = .true.

        If (compute_quantity(induct_vpbm_phi)) Then
            compute_vpbm_shear = .true.
            compute_vpbm_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vpbm_phi)) compute_vpbm_shear = .true.
        If (compute_quantity(induct_advec_vpbm_phi)) compute_vpbm_advec = .true.

        !5. Set flags related to v' x 'B'
        If (compute_quantity(induct_work_ppp) .or. compute_quantity(induct_work_mpp)) Then
            compute_vpbp_shear = .true.
            compute_vpbp_advec = .true.
        Endif

        If (compute_quantity(ishear_work_ppp)) compute_vpbp_shear = .true.
        If (compute_quantity(iadvec_work_ppp)) compute_vpbp_advec = .true.

        If (compute_quantity(ishear_work_mpp)) compute_vpbp_shear = .true.
        If (compute_quantity(iadvec_work_mpp)) compute_vpbp_advec = .true.

        If (compute_quantity(induct_vpbp_r)) Then
            compute_vpbp_shear = .true.
            compute_vpbp_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vpbp_r)) compute_vpbp_shear = .true.
        If (compute_quantity(induct_advec_vpbp_r)) compute_vpbp_advec = .true.

        If (compute_quantity(induct_vpbp_theta)) Then
            compute_vpbp_shear = .true.
            compute_vpbp_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vpbp_theta)) compute_vpbp_shear = .true.
        If (compute_quantity(induct_advec_vpbp_theta)) compute_vpbp_advec = .true.

        If (compute_quantity(induct_vpbp_phi)) Then
            compute_vpbp_shear = .true.
            compute_vpbp_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vpbp_phi)) compute_vpbp_shear = .true.
        If (compute_quantity(induct_advec_vpbp_phi)) compute_vpbp_advec = .true.

        !6. Set flags related to <v> x 'B'
        If (compute_quantity(induct_work_pmp)) Then
            compute_vmbp_shear = .true.
            compute_vmbp_advec = .true.
        Endif

        If (compute_quantity(ishear_work_pmp)) compute_vmbp_shear = .true.
        If (compute_quantity(iadvec_work_pmp)) compute_vmbp_advec = .true.


        If (compute_quantity(induct_vmbp_r)) Then
            compute_vmbp_shear = .true.
            compute_vmbp_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vmbp_r)) compute_vmbp_shear = .true.
        If (compute_quantity(induct_advec_vmbp_r)) compute_vmbp_advec = .true.

        If (compute_quantity(induct_vmbp_theta)) Then
            compute_vmbp_shear = .true.
            compute_vmbp_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vmbp_theta)) compute_vmbp_shear = .true.
        If (compute_quantity(induct_advec_vmbp_theta)) compute_vmbp_advec = .true.

        If (compute_quantity(induct_vmbp_phi)) Then
            compute_vmbp_shear = .true.
            compute_vmbp_advec = .true.
        Endif

        If (compute_quantity(induct_shear_vmbp_phi)) compute_vmbp_shear = .true.
        If (compute_quantity(induct_advec_vmbp_phi)) compute_vmbp_advec = .true.

    End Subroutine Reset_Induction_Flags

End Module Diagnostics_Induction
