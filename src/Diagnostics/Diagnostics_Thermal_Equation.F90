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

Module Diagnostics_Thermal_Equation
    Use Diagnostics_Base
    Use Diagnostics_ADotGradB
    Implicit None
Contains
    Subroutine Compute_Thermal_Equation_Terms(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Call Compute_Thermal_Advective_Terms(buffer)
        Call Compute_Thermal_Diffusion_Terms(buffer)
        Call Compute_Thermal_HeatSource(buffer)

        If (compute_quantity(visc_heating)) Then
            !Fairly complicated expression -- gets its own routine
            If (viscous_heating) Then
                Call Viscous_Heating_Diagnostics(buffer)
            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif
        If (magnetism) Call Compute_Ohmic_Heating_Diag(buffer)
    End Subroutine Compute_Thermal_Equation_Terms

    Subroutine Compute_Thermal_Diffusion_Terms(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Allocatable :: rhot(:), rhotk(:)
        Integer :: r,k, t

        Allocate(rhot(1:N_R), rhotk(1:N_R))
        rhot = ref%density*ref%temperature
        rhotk = rhot*kappa

        !////////////////////////////////////////////////////////////////////////
        ! Diffusive terms deriving from the divergence of the conductive flux

        ! First, the full-s terms
        If (compute_quantity(s_diff_phi) .or. compute_quantity(s_diff)) Then
            DO_PSI
                qty(PSI) = DDBUFF(PSI,dtdpdp)*rhotk(r)*(One_Over_R(r)*csctheta(t))**2
            END_DO
            If (compute_quantity(s_diff_phi)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(s_diff)) Then
                tmp1(:,:,:) = qty(:,:,:)
            Endif
        Endif

        If (compute_quantity(s_diff_theta) .or. compute_quantity(s_diff)) Then

            DO_PSI
                qty(PSI) = rhotk(r)*OneOverRSquared(r)*(DDBUFF(PSI,dtdtdt)+cottheta(t)*buffer(PSI,dtdt))
            END_DO

            If (compute_quantity(s_diff_theta)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(s_diff)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)
                END_DO
            Endif

        Endif

        If (compute_quantity(s_diff_r) .or. compute_quantity(s_diff)) Then

            DO_PSI
                qty(PSI)=Two_Over_R(r)+ref%dlnrho(r)+ref%dlnt(r)+dlnkappa(r)
            END_DO

            DO_PSI
                qty(PSI)=qty(PSI)*buffer(PSI,dtdr) + DDBUFF(PSI,dtdrdr)
            END_DO

            DO_PSI
                qty(PSI)=rhotk(r)*qty(PSI)
            END_DO

            If (compute_quantity(s_diff_r)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(s_diff)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif


        ! Next, fluctuating-s terms
        If (compute_quantity(sp_diff_phi) .or. compute_quantity(sp_diff)) Then
            DO_PSI
                qty(PSI) = d2_fbuffer(PSI,dtdpdp)*rhotk(r)*(One_Over_R(r)*csctheta(t))**2
            END_DO
            If (compute_quantity(sp_diff_phi)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(sp_diff)) Then
                tmp1(:,:,:) = qty(:,:,:)
            Endif
        Endif

        If (compute_quantity(sp_diff_theta) .or. compute_quantity(sp_diff)) Then

            DO_PSI
                qty(PSI) = rhotk(r)*OneOverRSquared(r)*(d2_fbuffer(PSI,dtdtdt)+cottheta(t)*fbuffer(PSI,dtdt))
            END_DO

            If (compute_quantity(sp_diff_theta)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(sp_diff)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)
                END_DO
            Endif

        Endif

        If (compute_quantity(sp_diff_r) .or. compute_quantity(sp_diff)) Then

            DO_PSI
                qty(PSI)=Two_Over_R(r)+ref%dlnrho(r)+ref%dlnt(r)+dlnkappa(r)
            END_DO

            DO_PSI
                qty(PSI)=qty(PSI)*fbuffer(PSI,dtdr) + d2_fbuffer(PSI,dtdrdr)
            END_DO

            DO_PSI
                qty(PSI)=rhotk(r)*qty(PSI)
            END_DO

            If (compute_quantity(sp_diff_r)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(sp_diff)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        ! Finally, mean s terms
        If (compute_quantity(sm_diff_phi) .or. compute_quantity(sm_diff)) Then
            DO_PSI
                qty(PSI) = d2_m0(PSI2,dtdpdp)*rhotk(r)*(One_Over_R(r)*csctheta(t))**2
            END_DO
            If (compute_quantity(sm_diff_phi)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(sm_diff)) Then
                tmp1(:,:,:) = qty(:,:,:)
            Endif
        Endif

        If (compute_quantity(sm_diff_theta) .or. compute_quantity(sm_diff)) Then

            DO_PSI
                qty(PSI) = rhotk(r)*OneOverRSquared(r)*(d2_m0(PSI2,dtdtdt)+cottheta(t)*m0_values(PSI2,dtdt))
            END_DO

            If (compute_quantity(sm_diff_theta)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(sm_diff)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)
                END_DO
            Endif

        Endif

        If (compute_quantity(sm_diff_r) .or. compute_quantity(sm_diff)) Then

            DO_PSI
                qty(PSI)=Two_Over_R(r)+ref%dlnrho(r)+ref%dlnt(r)+dlnkappa(r)
            END_DO

            DO_PSI
                qty(PSI)=qty(PSI)*m0_values(PSI2,dtdr) + d2_m0(PSI2,dtdrdr)
            END_DO

            DO_PSI
                qty(PSI)=rhotk(r)*qty(PSI)
            END_DO

            If (compute_quantity(sm_diff_r)) Then
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(sm_diff)) Then
                DO_PSI
                    tmp1(PSI) = tmp1(PSI)+qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif


        !//////////////////////////////////////////////////////////////////////////
        !  Fluxes

        !Conductive Heat Flux (full)
        If (compute_quantity(cond_flux_r)) Then
            DO_PSI
                qty(PSI) = -rhot(r)*kappa(r)*buffer(PSI,dtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(cond_flux_theta)) Then
            DO_PSI
                qty(PSI) = -rhot(r)*one_over_r(r)* &
                  & kappa(r)*buffer(PSI,dtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(cond_flux_phi)) Then
            DO_PSI
                qty(PSI) = -rhot(r)*one_over_r(r)*csctheta(t) &
                  & *kappa(r)*buffer(PSI,dtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !Conductive Heat Flux (fluctuating)
        If (compute_quantity(cond_fluxp_r)) Then
            DO_PSI
                qty(PSI) = -rhot(r)*kappa(r)*fbuffer(PSI,dtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(cond_fluxp_theta)) Then
            DO_PSI
                qty(PSI) = -rhot(r)*one_over_r(r)* &
                  & kappa(r)*fbuffer(PSI,dtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(cond_fluxp_phi)) Then
            DO_PSI
                qty(PSI) = -rhot(r)*one_over_r(r)*csctheta(t) &
                  & *kappa(r)*fbuffer(PSI,dtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !Conductive Flux (mean)
        If (compute_quantity(cond_fluxm_r)) Then
            DO_PSI
                qty(PSI) = -rhot(r)*kappa(r)*m0_values(PSI2,dtdr)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(cond_fluxm_theta)) Then
            DO_PSI
                qty(PSI) = -rhot(r)*one_over_r(r)* &
                  & kappa(r)*m0_values(PSI2,dtdt)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(cond_fluxm_phi)) Then
            DO_PSI
                qty(PSI) = -rhot(r)*one_over_r(r)*csctheta(t) &
                  & *kappa(r)*m0_values(PSI2,dtdp)
            END_DO
            Call Add_Quantity(qty)
        Endif

        DeAllocate(rhot, rhotk)
    End Subroutine Compute_Thermal_Diffusion_Terms

    Subroutine Compute_Thermal_Advective_Terms(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Allocatable :: rhot(:)
        Real*8 :: dt_by_ds, dt_by_dp
        Integer :: r,k, t

        Allocate(rhot(1:N_R))
        rhot = ref%density*ref%temperature


        ! Advective terms (full)
        If (Compute_Quantity(rhotv_grad_s)) Then
            DO_PSI
                qty(PSI)=buffer(PSI,vr)*buffer(PSI,dtdr)+ &
                         one_over_r(r)*(buffer(PSI,vtheta)*buffer(PSI,dtdt) + &
                         buffer(PSI,vphi)*buffer(PSI,dtdp)*csctheta(t))
                qty(PSI)=qty(PSI)*rhot(r)

            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvp_grad_sp)) Then
            DO_PSI
                qty(PSI)=fbuffer(PSI,vr)*fbuffer(PSI,dtdr)+ &
                         one_over_r(r)*(fbuffer(PSI,vtheta)*fbuffer(PSI,dtdt) + &
                         fbuffer(PSI,vphi)*fbuffer(PSI,dtdp)*csctheta(t))

                qty(PSI)=qty(PSI)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhotvp_grad_sm)) Then
            DO_PSI
                qty(PSI)=fbuffer(PSI,vr)*m0_values(PSI2,dtdr)+ &
                         one_over_r(r)*(fbuffer(PSI,vtheta)*m0_values(PSI2,dtdt) + &
                         fbuffer(PSI,vphi)*m0_values(PSI2,dtdp)*csctheta(t))
                qty(PSI)=qty(PSI)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvm_grad_sp)) Then
            DO_PSI
                qty(PSI)=m0_values(PSI2,vr)*fbuffer(PSI,dtdr)+ &
                         one_over_r(r)*(m0_values(PSI2,vtheta)*fbuffer(PSI,dtdt) + &
                         m0_values(PSI2,vphi)*fbuffer(PSI,dtdp)*csctheta(t))
                qty(PSI)=qty(PSI)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvm_grad_sm)) Then
            DO_PSI
                qty(PSI)=m0_values(PSI2,vr)*m0_values(PSI2,dtdr)+ &
                         one_over_r(r)*(m0_values(PSI2,vtheta)*m0_values(PSI2,dtdt) + &
                         m0_values(PSI2,vphi)*m0_values(PSI2,dtdp)*csctheta(t))
                qty(PSI)=qty(PSI)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Advective terms (pieces)
        ! Radial Contributions
        If (Compute_Quantity(rhotvr_grad_s)) Then
            DO_PSI
                qty(PSI)=buffer(PSI,vr)*buffer(PSI,dtdr)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvpr_grad_sp)) Then
            DO_PSI
                qty(PSI)=fbuffer(PSI,vr)*fbuffer(PSI,dtdr)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhotvpr_grad_sm)) Then
            DO_PSI
                qty(PSI)=fbuffer(PSI,vr)*m0_values(PSI2,dtdr)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvmr_grad_sp)) Then
            DO_PSI
                qty(PSI)=m0_values(PSI2,vr)*fbuffer(PSI,dtdr)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvmr_grad_sm)) Then
            DO_PSI
                qty(PSI)=m0_values(PSI2,vr)*m0_values(PSI2,dtdr)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Theta contributions
        If (Compute_Quantity(rhotvt_grad_s)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*buffer(PSI,vtheta)*buffer(PSI,dtdt)*rhot(r)

            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvpt_grad_sp)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*fbuffer(PSI,vtheta)*fbuffer(PSI,dtdt)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhotvpt_grad_sm)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*fbuffer(PSI,vtheta)*m0_values(PSI2,dtdt)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvmt_grad_sp)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*m0_values(PSI2,vtheta)*fbuffer(PSI,dtdt)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvmt_grad_sm)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*m0_values(PSI2,vtheta)*m0_values(PSI2,dtdt)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        ! Phi Contributions
        If (Compute_Quantity(rhotvp_grad_s)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*buffer(PSI,vphi)*buffer(PSI,dtdp)*csctheta(t)*rhot(r)

            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvpp_grad_sp)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*fbuffer(PSI,vphi)*fbuffer(PSI,dtdp)*csctheta(t)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhotvpp_grad_sm)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*fbuffer(PSI,vphi)*m0_values(PSI2,dtdp)*csctheta(t)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvmp_grad_sp)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*m0_values(PSI2,vphi)*fbuffer(PSI,dtdp)*csctheta(t)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        If (Compute_Quantity(rhotvmp_grad_sm)) Then
            DO_PSI
                qty(PSI)=one_over_r(r)*m0_values(PSI2,vphi)*m0_values(PSI2,dtdp)*csctheta(t)*rhot(r)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !/////////////////////////////////////////
        ! Compute advective heat FLUXES while here as well...

        !Radial
        If (Compute_Quantity(rhot_vr_s)) Then
            DO_PSI
                qty(PSI)=rhot(r)*buffer(PSI,vr)*buffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vrp_sp)) Then
            DO_PSI
                qty(PSI)=rhot(r)*fbuffer(PSI,vr)*fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vrp_sm)) Then
            DO_PSI
                qty(PSI)=rhot(r)*fbuffer(PSI,vr)*m0_values(PSI2,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vrm_sp)) Then
            DO_PSI
                qty(PSI)=rhot(r)*m0_values(PSI2,vr)*fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vrm_sm)) Then
            DO_PSI
                qty(PSI)=rhot(r)*m0_values(PSI2,vr)*m0_values(PSI2,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !Theta
        If (Compute_Quantity(rhot_vt_s)) Then
            DO_PSI
                qty(PSI)=rhot(r)*buffer(PSI,vtheta)*buffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vtp_sp)) Then
            DO_PSI
                qty(PSI)=rhot(r)*fbuffer(PSI,vtheta)*fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vtp_sm)) Then
            DO_PSI
                qty(PSI)=rhot(r)*fbuffer(PSI,vtheta)*m0_values(PSI2,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vtm_sp)) Then
            DO_PSI
                qty(PSI)=rhot(r)*m0_values(PSI2,vtheta)*fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vtm_sm)) Then
            DO_PSI
                qty(PSI)=rhot(r)*m0_values(PSI2,vtheta)*m0_values(PSI2,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif

        !Phi
        If (Compute_Quantity(rhot_vp_s)) Then
            DO_PSI
                qty(PSI)=rhot(r)*buffer(PSI,vphi)*buffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vtp_sp)) Then
            DO_PSI
                qty(PSI)=rhot(r)*fbuffer(PSI,vphi)*fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vtp_sm)) Then
            DO_PSI
                qty(PSI)=rhot(r)*fbuffer(PSI,vphi)*m0_values(PSI2,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vtm_sp)) Then
            DO_PSI
                qty(PSI)=rhot(r)*m0_values(PSI2,vphi)*fbuffer(PSI,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif
        If (Compute_Quantity(rhot_vtm_sm)) Then
            DO_PSI
                qty(PSI)=rhot(r)*m0_values(PSI2,vphi)*m0_values(PSI2,tvar)
            END_DO
            Call Add_Quantity(qty)
        Endif


        !////////////////////////////////////////////////
        ! Finally, compute the enthalpy fluxes
        If (compute_quantity(enth_flux_r) .or. compute_quantity(enth_flux_theta) .or. &
            compute_quantity(enth_flux_phi)) Then

            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    dt_by_ds = ref%temperature(r)/pressure_specific_heat
                    dt_by_dp = 1.0d0/pressure_specific_heat/ref%density(r)
                    Do k = 1, n_phi
                        tmp1(PSI) = pressure_specific_heat*ref%density(r) * &
                            (dt_by_ds*buffer(PSI,tvar) + dt_by_dp*buffer(PSI,pvar))
                    Enddo
                Enddo
            Enddo

            If (compute_quantity(enth_flux_r)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,vr)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_theta)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,vtheta)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_phi)) Then
                DO_PSI
                    qty(PSI) = buffer(PSI,vphi)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(enth_flux_rpp) .or. compute_quantity(enth_flux_thetapp) .or. &
            compute_quantity(enth_flux_phipp)) Then

            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    dt_by_ds = ref%temperature(r)/pressure_specific_heat
                    dt_by_dp = 1.0d0/pressure_specific_heat/ref%density(r)
                    Do k = 1, n_phi
                        tmp1(PSI) = pressure_specific_heat*ref%density(r) * &
                            (dt_by_ds*fbuffer(PSI,tvar) + dt_by_dp*fbuffer(PSI,pvar))
                    Enddo
                Enddo
            Enddo

            If (compute_quantity(enth_flux_rpp)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,vr)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_thetapp)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,vtheta)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_phipp)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,vphi)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(enth_flux_rpm) .or. compute_quantity(enth_flux_thetapm) .or. &
            compute_quantity(enth_flux_phipm)) Then

            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    dt_by_ds = ref%temperature(r)/pressure_specific_heat
                    dt_by_dp = 1.0d0/pressure_specific_heat/ref%density(r)
                    Do k = 1, n_phi
                        tmp1(PSI) = pressure_specific_heat*ref%density(r) * &
                            (dt_by_ds*m0_values(PSI2,tvar) + dt_by_dp*m0_values(PSI2,pvar))
                    Enddo
                Enddo
            Enddo

            If (compute_quantity(enth_flux_rpm)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,vr)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_thetapm)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,vtheta)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_phipm)) Then
                DO_PSI
                    qty(PSI) = fbuffer(PSI,vphi)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(enth_flux_rmp) .or. compute_quantity(enth_flux_thetamp) .or. &
            compute_quantity(enth_flux_phimp)) Then

            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    dt_by_ds = ref%temperature(r)/pressure_specific_heat
                    dt_by_dp = 1.0d0/pressure_specific_heat/ref%density(r)
                    Do k = 1, n_phi
                        tmp1(PSI) = pressure_specific_heat*ref%density(r) * &
                            (dt_by_ds*fbuffer(PSI,tvar) + dt_by_dp*fbuffer(PSI,pvar))
                    Enddo
                Enddo
            Enddo

            If (compute_quantity(enth_flux_rmp)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vr)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_thetamp)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vtheta)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_phimp)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vphi)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(enth_flux_rmm) .or. compute_quantity(enth_flux_thetamm) .or. &
            compute_quantity(enth_flux_phimm)) Then

            Do t = my_theta%min, my_theta%max
                Do r = my_r%min, my_r%max
                    dt_by_ds = ref%temperature(r)/pressure_specific_heat
                    dt_by_dp = 1.0d0/pressure_specific_heat/ref%density(r)
                    Do k = 1, n_phi
                        tmp1(PSI) = pressure_specific_heat*ref%density(r) * &
                            (dt_by_ds*m0_values(PSI2,tvar) + dt_by_dp*m0_values(PSI2,pvar))
                    Enddo
                Enddo
            Enddo

            If (compute_quantity(enth_flux_rmm)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vr)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_thetamm)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vtheta)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

            If (compute_quantity(enth_flux_phimm)) Then
                DO_PSI
                    qty(PSI) = m0_values(PSI2,vphi)*tmp1(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        DeAllocate(rhot)
    End Subroutine Compute_Thermal_Advective_Terms

    Subroutine Compute_Ohmic_Heating_Diag(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8 :: dt_by_ds, dt_by_dp
        Integer :: r,k, t
        ! The "Diag" in the name refers to "diagnostics," and not "diagonal"

        If (compute_quantity(ohmic_heat)) Then
            If (ohmic_heating) Then
                DO_PSI
                    qty(PSI) = ohmic_heating_coeff(r)*(buffer(PSI,curlbr)**2 + &
                               &   buffer(PSI,curlbtheta)**2 + &
                               &   buffer(PSI,curlbphi)**2)
                END_DO
            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif
        If (compute_quantity(ohmic_heat_pp)) Then
            If (ohmic_heating) Then
                DO_PSI
                    qty(PSI) = ohmic_heating_coeff(r)*(fbuffer(PSI,curlbr)**2 + &
                               &   fbuffer(PSI,curlbtheta)**2 + &
                               &   fbuffer(PSI,curlbphi)**2)
                END_DO
            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif

        If (compute_quantity(ohmic_heat_pm)) Then
            If (ohmic_heating) Then
                DO_PSI
                    qty(PSI) = ohmic_heating_coeff(r)*(fbuffer(PSI,curlbr)*m0_values(PSI2,curlbr) + &
                               &   fbuffer(PSI,curlbtheta)*m0_values(PSI2,curlbtheta) + &
                               &   fbuffer(PSI,curlbphi)*m0_values(PSI2,curlbphi))
                END_DO
            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif


        If (compute_quantity(ohmic_heat_mm)) Then
            If (ohmic_heating) Then
                DO_PSI
                    qty(PSI) = ohmic_heating_coeff(r)*(m0_values(PSI2,curlbr)**2 + &
                               &   m0_values(PSI2,curlbtheta)**2 + &
                               &   m0_values(PSI2,curlbphi)**2)
                END_DO

            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Ohmic_Heating_Diag


    Subroutine Viscous_Heating_Diagnostics(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: t,r,k
        Real*8 :: tmp, tmp2
        ! Computes the viscous heating term that appears in the
        ! thermal energy equation

        ! Computing viscous heating means contracting a 3x3 tensor
        ! To avoid code-bloat, we do not performance a Reynolds-like
        ! decomposition on this quantity.  We comput only the full viscous heating.

        !Contributions from E_rr, E_theta_theta     & E_phi_phi

        DO_PSI
            tmp = (buffer(PSI,dvpdp)*csctheta(t) +buffer(PSI,vr) &
                    +buffer(PSI,vtheta)*cottheta(t))*one_over_r(r)    !e_phi_phi
            tmp2 = (buffer(PSI,dvtdt)+buffer(PSI,vr))*one_over_r(r) ! e_theta_theta
            qty(PSI) = buffer(PSI,dvrdr)*buffer(PSI,dvrdr)+tmp*tmp +tmp2*tmp2
        END_DO


        !E_r_phi
        DO_PSI
            tmp = (buffer(PSI,dvrdp)*csctheta(t)- buffer(PSI,vphi))*one_over_r(r) &
                    +buffer(PSI,dvpdr) ! 2*e_r_phi
            qty(PSI) = qty(PSI)+tmp*tmp*Half  ! +2 e_r_phi**2
        END_DO



        !E_r_theta

        DO_PSI
            tmp = (buffer(PSI,dvrdt)-buffer(PSI,vtheta))*one_over_r(r) &
                    +buffer(PSI,dvtdr) ! 2*e_r_theta

            qty(PSI) = qty(PSI)+tmp*tmp*Half   ! + 2+e_r_theta**2

        END_DO


        !E_phi_theta
        DO_PSI
            tmp = (buffer(PSI,dvpdt) &
                +buffer(PSI,dvtdp)*csctheta(t) &
                -buffer(PSI,vphi)*cottheta(t) )*one_over_r(r)        ! 2*e_phi_theta

            qty(PSI) = qty(PSI)+tmp*tmp*Half   ! + 2*e_phi_theta**2
        END_DO



        ! -1/3 (div dot v )**2
        DO_PSI
            tmp = -buffer(PSI,vr)*ref%dlnrho(r)
            qty(PSI) = qty(PSI)-tmp*tmp*one_third   ! + 2*e_phi_theta**2
        END_DO
        DO_PSI
            qty(PSI) = viscous_heating_coeff(r)*qty(PSI)
        END_DO


    End Subroutine Viscous_Heating_Diagnostics




    Subroutine Compute_Thermal_HeatSource(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Real*8, Allocatable :: rhot(:)
        Real*8 :: mean_rho, mean_t, mean_r2, mean_q, mean_dr
        Real*8 :: qadd, fpr2dr
        Integer :: r,k, t
        ! Volume Heating
        If (compute_quantity(vol_heating)) Then
            If (allocated(ref%heating)) Then
                DO_PSI
                    qty(PSI) = ref%heating(r)* &
                        & ref%density(r)*ref%temperature(r)
                END_DO
            Else
                qty(:,:,:) = 0.0d0
            Endif
            Call Add_Quantity(qty)
        Endif

        !The "Flux" associated with the volume heating
        If (compute_quantity(vol_heat_flux)) Then
            tmp1d(:) = 0.0d0
            If (heating_type .gt. 0) Then

                ! Note that radial_integral_weights give int{f r^2}/int(r^2}
                Do r = N_R-1, 1,-1
                    mean_rho = half*(ref%density(r)     +  ref%density(r+1)    )
                    mean_t   = half*(ref%temperature(r) +  ref%temperature(r+1))
                    mean_r2  = half*(r_squared(r)       +  r_squared(r+1)      )
                    mean_q   = half*(ref%heating(r)     +  ref%heating(r+1)    )
                    mean_dr  = radius(r)-radius(r+1)

                    qadd = mean_q*mean_rho*mean_t
                    fpr2dr = mean_r2*four_pi*mean_dr

                    tmp1d(r) = tmp1d(r+1)+qadd*fpr2dr

                Enddo
                tmp1d = tmp1d(1)-tmp1d
                tmp1d = tmp1d/four_pi/r_squared
            Endif
            DO_PSI
                qty(PSI) = tmp1d(r)
            END_DO
            Call Add_Quantity(qty)
        Endif
    End Subroutine Compute_Thermal_HeatSource

End Module Diagnostics_Thermal_Equation
