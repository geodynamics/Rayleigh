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

Module Diagnostics_Inertial_Forces

    Use Diagnostics_Base
    Use Diagnostics_ADotGradB
Contains
    Subroutine Compute_Inertial_Terms(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8, Allocatable :: cbuffer(:,:,:,:)



        Logical :: compute_fluctuations = .false.
        Logical :: compute_full_full    = .false.
        Logical :: compute_fluct_mean   = .false.
        Logical :: compute_mean_fluct   = .false.
        Logical :: compute_fluct_fluct  = .false.
        Logical :: compute_mean_mean    = .false.

        !First, figure out what we need to bother computing

        !///-- Flags related to Angular Momentum
        If (compute_quantity(samom_advec_pp)) Then
            compute_fluct_fluct  = .true.
        Endif
        If (compute_quantity(samom_advec_mm)) Then
            compute_mean_mean  = .true.
        Endif

        !///-- Flags for full dot grad full terms
        If (compute_quantity(v_grad_v_r)) Then
            compute_full_full = .true.
        Endif
        If (compute_quantity(v_grad_v_theta)) Then
            compute_full_full = .true.
        Endif
        If (compute_quantity(v_grad_v_phi)) Then
            compute_full_full = .true.
        Endif

        !///-- Flags for mean dot grad mean terms
        If (compute_quantity(vm_grad_vm_r)) Then
            compute_mean_mean   = .true.
        Endif
        If (compute_quantity(vm_grad_vm_theta)) Then
            compute_mean_mean   = .true.
        Endif
        If (compute_quantity(vm_grad_vm_phi)) Then
            compute_mean_mean   = .true.
        Endif

        !///-- Flags for fluctuating dot grad mean terms
        If (compute_quantity(vp_grad_vm_r)) Then
            compute_fluctuations = .true.
            compute_fluct_mean   = .true.
        Endif
        If (compute_quantity(vp_grad_vm_theta)) Then
            compute_fluctuations = .true.
            compute_fluct_mean   = .true.
        Endif
        If (compute_quantity(vp_grad_vm_phi)) Then
            compute_fluctuations = .true.
            compute_fluct_mean   = .true.
        Endif

        !///-- Flags for mean dot grad fluctuating terms
        If (compute_quantity(vm_grad_vp_r)) Then
            compute_fluctuations = .true.
            compute_mean_fluct   = .true.
        Endif
        If (compute_quantity(vm_grad_vp_theta)) Then
            compute_fluctuations = .true.
            compute_mean_fluct   = .true.
        Endif
        If (compute_quantity(vm_grad_vp_phi)) Then
            compute_fluctuations = .true.
            compute_mean_fluct   = .true.
        Endif

        !///-- Flags for fluctuating dot grad fluctuating terms
        If (compute_quantity(vp_grad_vp_r)) Then
            compute_fluctuations = .true.
            compute_fluct_fluct  = .true.
        Endif
        If (compute_quantity(vp_grad_vp_theta)) Then
            compute_fluctuations = .true.
            compute_fluct_fluct  = .true.
        Endif
        If (compute_quantity(vp_grad_vp_phi)) Then
            compute_fluctuations = .true.
            compute_fluct_fluct  = .true.
        Endif

        Allocate(cbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:3))


        !//////////////////////////////////////////////////////////////////////////////////
        !/////////////// v dot grad v (full)//////////////////
        If (compute_full_full .or. compute_quantity(advec_work)) Then
            Call ADotGradB(buffer,buffer,cbuffer,aindices=vindex,bindices=vindex)

            If (compute_quantity(v_grad_v_r) .or. compute_quantity(advec_work)) Then
                DO_PSI
                    qty(PSI) = mean_3dbuffer(PSI,aforce_r)-mean_ell0buffer(r,aforce_r)
                END_DO
                If (compute_quantity(v_grad_v_r)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work)) Then
                    DO_PSI
                        tmp1(PSI) = -qty(PSI)*buffer(PSI,vr)
                    END_DO
                Endif
            Endif

            If (compute_quantity(v_grad_v_theta) .or. compute_quantity(advec_work)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO
                If (compute_quantity(v_grad_v_theta)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI)-qty(PSI)*buffer(PSI,vtheta)
                    END_DO
                Endif
            Endif

            If (compute_quantity(v_grad_v_phi) .or. compute_quantity(advec_work)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO
                If (compute_quantity(v_grad_v_phi)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI)-qty(PSI)*buffer(PSI,vphi)
                    END_DO
                    Call Add_Quantity(tmp1)
                Endif
            Endif

        Endif

        !/////////////// v' dot grad v' //////////////////
        If (compute_fluct_fluct .or. compute_quantity(advec_work_ppp) .or. &
            compute_quantity(advec_work_mpp)) Then
            Call ADotGradB(fbuffer,fbuffer,cbuffer,aindices=vindex,bindices=vindex)
            If (compute_quantity(vp_grad_vp_r) .or. compute_quantity(advec_work_ppp) &
                .or. compute_quantity(advec_work_mpp)) Then

                DO_PSI
                    qty(PSI) = mean_3dbuffer(PSI,aforcepp_r)-mean_ell0buffer(r,aforcepp_r) !cbuffer(PSI,1)*ref%density(r)
                END_DO
                If (compute_quantity(vp_grad_vp_r)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_ppp)) Then
                    DO_PSI
                        tmp1(PSI) = -qty(PSI)*fbuffer(PSI,vr)
                    END_DO
                Endif
                If (compute_quantity(advec_work_mpp)) Then
                    DO_PSI
                        tmp4(PSI) = -qty(PSI)*m0_values(PSI2,vr)
                    END_DO
                Endif
            Endif
            If (compute_quantity(vp_grad_vp_theta) .or. compute_quantity(advec_work_ppp) &
                .or. compute_quantity(advec_work_mpp)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO
                If (compute_quantity(vp_grad_vp_theta)) Call Add_Quantity(qty)

                If (compute_quantity(advec_work_ppp)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI) - qty(PSI)*fbuffer(PSI,vtheta)
                    END_DO
                Endif

                If (compute_quantity(advec_work_mpp)) Then
                    DO_PSI
                        tmp4(PSI) = tmp4(PSI) - qty(PSI)*m0_values(PSI2,vtheta)
                    END_DO
                Endif

            Endif
            If (compute_quantity(vp_grad_vp_phi) .or. compute_quantity(samom_advec_pp) &
                .or. compute_quantity(advec_work_ppp) .or. compute_quantity(advec_work_mpp)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO
                If (compute_quantity(vp_grad_vp_phi)) Call Add_Quantity(qty)

                If (compute_quantity(advec_work_ppp)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI) - qty(PSI)*fbuffer(PSI,vphi)
                    END_DO
                    Call Add_Quantity(tmp1)
                Endif

                If (compute_quantity(advec_work_mpp)) Then
                    DO_PSI
                        tmp4(PSI) = tmp4(PSI) - qty(PSI)*m0_values(PSI2,vphi)
                    END_DO
                    Call Add_Quantity(tmp4)
                Endif
                If (compute_quantity(samom_advec_pp)) Then
                    DO_PSI
                        qty(PSI) = qty(PSI)*radius(r)*sintheta(t)
                    END_DO
                    Call Add_Quantity(qty)
                Endif
            Endif
        Endif

        !/////////////// <v> dot grad <v> //////////////////
        If (compute_mean_mean .or. compute_quantity(advec_work_mmm)) Then
            Call ADotGradB(m0_values,m0_values,cbuffer,aindices=vindex,bindices=vindex)

            If (compute_quantity(vm_grad_vm_r) .or. compute_quantity(advec_work_mmm)) Then
                DO_PSI
                    qty(PSI) = mean_3dbuffer(PSI,aforcemm_r)-mean_ell0buffer(r,aforcemm_r) ! cbuffer(PSI,1)*ref%density(r)
                END_DO

                If (compute_quantity(vm_grad_vm_r)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_mmm)) Then
                    DO_PSI
                        tmp1(PSI) = -qty(PSI)*m0_values(PSI2,vr)
                    END_DO
                Endif

            Endif

            If (compute_quantity(vm_grad_vm_theta) .or. compute_quantity(advec_work_mmm)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO
                If (compute_quantity(vm_grad_vm_theta)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_mmm)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI)-qty(PSI)*m0_values(PSI2,vtheta)
                    END_DO
                Endif
            Endif

            If (compute_quantity(vm_grad_vm_phi) .or. compute_quantity(samom_advec_mm) &
                .or. compute_quantity(advec_work_mmm)) Then

                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO

                If (compute_quantity(vm_grad_vm_phi)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_mmm)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI)-qty(PSI)*m0_values(PSI2,vphi)
                    END_DO
                    Call Add_Quantity(tmp1)
                Endif
                If (compute_quantity(samom_advec_mm)) Then
                    DO_PSI
                        qty(PSI) = qty(PSI)*radius(r)*sintheta(t)
                    END_DO
                    Call Add_Quantity(qty)
                Endif
            Endif
        Endif

        !/////////////// v' dot grad <v> //////////////////
        If (compute_fluct_mean .or. compute_quantity(advec_work_ppm)) Then

            Call ADotGradB(fbuffer,m0_values,cbuffer,aindices = vindex, bindices=vindex)

            If (compute_quantity(vp_grad_vm_r) .or. compute_quantity(advec_work_ppm)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*ref%density(r)
                END_DO
                If (compute_quantity(vp_grad_vm_r)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_ppm)) Then
                    DO_PSI
                        tmp1(PSI) = -qty(PSI)*fbuffer(PSI,vr)
                    END_DO
                Endif

            Endif

            If (compute_quantity(vp_grad_vm_theta) .or. compute_quantity(advec_work_ppm)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO
                If (compute_quantity(vp_grad_vm_theta)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_ppm)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI)-qty(PSI)*fbuffer(PSI,vtheta)
                    END_DO
                Endif
            Endif

            If (compute_quantity(vp_grad_vm_phi) .or. compute_quantity(advec_work_ppm)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO

                If (compute_quantity(vp_grad_vm_phi)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_ppm)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI)-qty(PSI)*fbuffer(PSI,vphi)
                    END_DO
                    Call Add_Quantity(tmp1)
                Endif
            Endif
        Endif

        !/////////////// <v> dot grad v' //////////////////
        If (compute_mean_fluct .or. compute_quantity(advec_work_pmp)) Then

            Call ADotGradB(m0_values,fbuffer,cbuffer,aindices=vindex,bindices=vindex)

            If (compute_quantity(vm_grad_vp_r) .or. compute_quantity(advec_work_pmp)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,1)*ref%density(r)
                END_DO
                If (compute_quantity(vm_grad_vp_r)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_pmp)) Then
                    DO_PSI
                        tmp1(PSI) = -qty(PSI)*fbuffer(PSI,vphi)
                    END_DO
                Endif
            Endif

            If (compute_quantity(vm_grad_vp_theta) .or. compute_quantity(advec_work_pmp)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,2)*ref%density(r)
                END_DO

                If (compute_quantity(vm_grad_vp_theta)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_pmp)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI)-qty(PSI)*fbuffer(PSI,vphi)
                    END_DO
                Endif
            Endif

            If (compute_quantity(vm_grad_vp_phi) .or. compute_quantity(advec_work_pmp)) Then
                DO_PSI
                    qty(PSI) = cbuffer(PSI,3)*ref%density(r)
                END_DO
                If (compute_quantity(vm_grad_vp_phi)) Call Add_Quantity(qty)
                If (compute_quantity(advec_work_pmp)) Then
                    DO_PSI
                        tmp1(PSI) = tmp1(PSI)-qty(PSI)*fbuffer(PSI,vphi)
                    END_DO
                    Call Add_Quantity(tmp1)
                Endif
            Endif
        Endif

        DeAllocate(cbuffer)

    End Subroutine Compute_Inertial_Terms
End Module Diagnostics_Inertial_Forces
