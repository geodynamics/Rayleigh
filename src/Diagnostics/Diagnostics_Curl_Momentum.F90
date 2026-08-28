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

Module Diagnostics_Curl_Momentum
    Use Diagnostics_Base
    Use Spectral_Derivatives
    Use Finite_Difference, Only : d_by_dx3d3
    Implicit None

    Integer, Allocatable :: vfdindmap(:,:)
    Integer :: nvffields

Contains

    Subroutine Compute_Curl_Momentum_Forces(buffer)
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Call Compute_Curl_Advection_Force(buffer)
        Call Compute_Curl_Buoyancy_Force(buffer)
        Call Compute_Curl_Magnetic_Force(buffer)
        Call Compute_Curl_Coriolis_Force(buffer)
        Call Compute_Curl_Pressure_Force(buffer)
        Call Compute_Curl_Viscous_Force(buffer)
    End Subroutine

    Subroutine Compute_Curl_Advection_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t

        Real*8  :: pfactor(my_r%min:my_r%max)
        pfactor(my_r%min:my_r%max) = ref%dpdr_w_term(my_r%min:my_r%max) &
                                        /ref%density(my_r%min:my_r%max)
        
        !!!!!!!!!!!!!!!!!!!!!!!!! Advection Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        If (compute_quantity(curl_v_grad_v_r) .or. compute_quantity(curl_v_grad_v_r_squared)) Then
            DO_PSI
                qty(PSI) = one_over_r(r) * ref%density(r) * ( &
                                buffer(PSI,vr) * one_over_r(r) * DDBUFF(PSI,dvpdrdt) + &
                                buffer(PSI,dvrdt) * one_over_r(r) * buffer(PSI,dvpdr) + &
                                buffer(PSI,vtheta) * one_over_r(r) * DDBUFF(PSI,dvpdtdt) + &
                                buffer(PSI,dvtdt) * one_over_r(r) * buffer(PSI,dvpdt) + &
                                buffer(PSI,vphi) * one_over_r(r) * (csctheta(t) * DDBUFF(PSI,dvpdtdp) - &
                                csctheta(t) * csctheta(t) * costheta(t) * buffer(PSI,dvpdp) + &
                                buffer(PSI,dvrdt) + &
                                cottheta(t) * buffer(PSI,dvtdt) - &
                                csctheta(t) * csctheta(t) * buffer(PSI,vtheta)) + &
                                buffer(PSI,dvpdt) * one_over_r(r) * (csctheta(t) * buffer(PSI,dvpdp) + &
                                buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) + &
                            one_over_r(r) * costheta(t) * csctheta(t) * ref%density(r) * (buffer(PSI,vr) * &
                            one_over_r(r) * buffer(PSI,dvpdr) + buffer(PSI,vtheta) * one_over_r(r) * &
                            buffer(PSI,dvpdt) + buffer(PSI,vphi) * one_over_r(r) * (csctheta(t) * &
                            buffer(PSI,dvpdp) + buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) - &
                            one_over_r(r) * csctheta(t) * ref%density(r) * (buffer(PSI,vr) * DDBUFF(PSI,dvtdrdp) + &
                            buffer(PSI,dvrdp) * buffer(PSI,dvtdr) + buffer(PSI,vtheta) * one_over_r(r) * &
                            (DDBUFF(PSI,dvtdtdp) + buffer(PSI,dvrdp)) + buffer(PSI,dvtdp) * one_over_r(r) * &
                            (buffer(PSI,dvtdt) + buffer(PSI,vr)) + buffer(PSI,vphi) * one_over_r(r) * &
                            (csctheta(t) * DDBUFF(PSI,dvtdpdp) - cottheta(t) * buffer(PSI,dvpdp)) + &
                            buffer(PSI,dvpdp) * one_over_r(r) * (csctheta(t) * buffer(PSI,dvtdp) - &
                            cottheta(t) * buffer(PSI,vphi)))                                       
            END_DO
            If (compute_quantity(curl_v_grad_v_r)) Call Add_Quantity(qty)
            If (compute_quantity(curl_v_grad_v_r_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif            
            
        Endif
    
        If (compute_quantity(curl_v_grad_v_theta) .or. compute_quantity(curl_v_grad_v_theta_squared)) Then
            DO_PSI
                qty(PSI) = one_over_r(r) * csctheta(t) * ref%density(r) * ( buffer(PSI,vr) * DDBUFF(PSI,dvrdrdp) + &
                            buffer(PSI,dvrdp) * buffer(PSI,dvrdr) + buffer(PSI,vtheta) * one_over_r(r) * &
                            (DDBUFF(PSI,dvrdtdp) - buffer(PSI,dvtdp)) + buffer(PSI,dvtdp) * one_over_r(r) * &
                            (buffer(PSI,dvrdt) - buffer(PSI,vtheta)) + buffer(PSI,vphi) * one_over_r(r) * &
                            (csctheta(t) * DDBUFF(PSI,dvrdpdp) - buffer(PSI,dvpdp)) + buffer(PSI,dvpdp) * &
                            one_over_r(r) * ( csctheta(t) * buffer(PSI,dvrdp) - buffer(PSI,vphi))) - &
                            one_over_r(r) * ref%density(r) * (buffer(PSI,vr) * one_over_r(r) * buffer(PSI,dvpdr) + &
                            buffer(PSI,vtheta) * one_over_r(r) * buffer(PSI,dvpdt) + buffer(PSI,vphi) * one_over_r(r) * &
                            (csctheta(t) * buffer(PSI,dvpdp) + buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) - &
                            ref%density(r) * (buffer(PSI,vr) * one_over_r(r) * DDBUFF(PSI,dvpdrdr) + &
                                buffer(PSI,dvrdr) * one_over_r(r) * buffer(PSI,dvpdr) - &
                                buffer(PSI,vr) * one_over_r(r) * one_over_r(r) * buffer(PSI,dvpdr) + &
                                buffer(PSI,vtheta) * one_over_r(r) * DDBUFF(PSI,dvpdrdt) + &
                                buffer(PSI,dvtdr) * one_over_r(r) * buffer(PSI,dvpdt) - &
                                buffer(PSI,vtheta) * one_over_r(r) * one_over_r(r) * buffer(PSI,dvpdt) + &
                                buffer(PSI,vtheta) * one_over_r(r) * (csctheta(t) * DDBUFF(PSI,dvpdrdp) + &
                                buffer(PSI,dvrdr) + &
                                cottheta(t) * buffer(PSI,dvtdr)) + &
                                buffer(PSI,dvpdr) * one_over_r(r) * (csctheta(t) * buffer(PSI,dvpdp) + &
                                buffer(PSI,vr) + &
                                cottheta(t) * buffer(PSI,vtheta)) - &
                                buffer(PSI,vphi) *  one_over_r(r) * one_over_r(r) * (csctheta(t) *  buffer(PSI,dvpdp) + &
                                buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) - &
                            ref%dlnrho(r) * &
                            ref%density(r) * (buffer(PSI,vr) * one_over_r(r) * buffer(PSI,dvpdr) + &
                            buffer(PSI,vtheta) * one_over_r(r) * buffer(PSI,dvpdt) + buffer(PSI,vphi) * one_over_r(r) * &
                            (csctheta(t) * buffer(PSI,dvpdp) + buffer(PSI,vr) + cottheta(t) * buffer(PSI,vtheta))) 
            END_DO
            If (compute_quantity(curl_v_grad_v_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_v_grad_v_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif
        
        If (compute_quantity(curl_v_grad_v_phi) .or. compute_quantity(curl_v_grad_v_phi_squared)) Then
            DO_PSI
                qty(PSI) = ref%density(r) * (buffer(PSI,vr) * DDBUFF(PSI,dvtdrdr) + &
                                buffer(PSI,dvrdr) * buffer(PSI,dvtdr) + &
                                buffer(PSI,vtheta) * one_over_r(r) * (DDBUFF(PSI,dvtdrdt) + &
                                buffer(PSI,dvrdr)) + &
                                buffer(PSI,dvtdr) * one_over_r(r) * (buffer(PSI,dvtdt) + &
                                buffer(PSI,vr)) - &
                                buffer(PSI,vtheta) * one_over_r(r) * one_over_r(r) * (buffer(PSI,dvtdt) + &
                                buffer(PSI,vr)) + &
                                buffer(PSI,vtheta) * one_over_r(r) * (csctheta(t) * DDBUFF(PSI,dvtdrdp) - &
                                cottheta(t) * buffer(PSI,dvpdr)) + &
                                buffer(PSI,dvpdr) * one_over_r(r) * (csctheta(t) * buffer(PSI,dvtdp) - &
                                cottheta(t) * buffer(PSI,vphi)) - &
                                buffer(PSI,vphi) * one_over_r(r) * one_over_r(r) * (csctheta(t) * buffer(PSI,dvtdp) - &
                                cottheta(t) * buffer(PSI,vphi))) + &   
                            ref%dlnrho(r) * ref%density(r) * (buffer(PSI,vr) * buffer(PSI,dvtdr) + buffer(PSI,vtheta) * &
                            one_over_r(r) * (buffer(PSI,dvtdt) + buffer(PSI,vr)) + buffer(PSI,vphi) * one_over_r(r) * &
                            (csctheta(t) * buffer(PSI,dvtdp) - cottheta(t) * buffer(PSI,vphi))) + &
                            one_over_r(r) * ref%density(r) * (buffer(PSI,vr) * buffer(PSI,dvtdr) + buffer(PSI,vtheta) * &
                            one_over_r(r) * (buffer(PSI,dvtdt) + buffer(PSI,vr)) + buffer(PSI,vphi) * one_over_r(r) * &
                            (csctheta(t) * buffer(PSI,dvtdp) - cottheta(t) * buffer(PSI,vphi))) - &
                            one_over_r(r) * ref%density(r) * (buffer(PSI,vr) * DDBUFF(PSI,dvrdrdt) + &
                                buffer(PSI,dvrdt) * buffer(PSI,dvrdr) + buffer(PSI,vtheta) * one_over_r(r) * &
                                (DDBUFF(PSI,dvrdtdt) + buffer(PSI,dvtdt)) + &
                            buffer(PSI,dvtdt) * one_over_r(r) * (buffer(PSI,dvrdt) - buffer(PSI,vtheta)) + &
                            buffer(PSI,vphi) * one_over_r(r) * (csctheta(t) * DDBUFF(PSI,dvrdtdp) - &
                                csctheta(t) * csctheta(t) * costheta(t) * buffer(PSI,dvrdp) - &
                                buffer(PSI,dvpdt)) + &
                            buffer(PSI,dvpdt) * one_over_r(r) * &
                            (csctheta(t) * buffer(PSI,dvrdp) - buffer(PSI,vphi)))
            END_DO
            If (compute_quantity(curl_v_grad_v_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_v_grad_v_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        If (compute_quantity(curl_v_grad_v_abs)) Then
            DO_PSI
                qty(PSI) =((one_over_r(r) * one_over_r(r) * buffer(PSI,dvtdt) * buffer(PSI,dvpdt) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,vtheta) * DDBUFF(PSI,dvpdtdt) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * buffer(PSI,vtheta) * buffer(PSI,vphi) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,dvtdt) * buffer(PSI,vphi) + &                         
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,vtheta) * &
                           buffer(PSI,dvpdt) - &
                           one_over_r(r) * buffer(PSI,dvtdt) * buffer(PSI,dvtdp) - &
                           one_over_r(r) * buffer(PSI,vtheta) * DDBUFF(PSI,dvtdtdp) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                           buffer(PSI,vr) * buffer(PSI,dvrdp) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * &
                           buffer(PSI,dvrdt) * buffer(PSI,dvrdp) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * &
                           buffer(PSI,vr) * DDBUFF(PSI,dvrdtdp) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,dvrdt) * buffer(PSI,vphi) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,vr) * buffer(PSI,dvpdt) + &
                           one_over_r(r) * buffer(PSI,dvrdt) * buffer(PSI,dvpdr) + &
                           one_over_r(r) * buffer(PSI,vr) * DDBUFF(PSI,dvpdrdt) + &                             !(first part ended)
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,vtheta) * buffer(PSI,dvpdt) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                           costheta(t) * costheta(t) * buffer(PSI,vtheta) * buffer(PSI,vphi) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,vtheta) * buffer(PSI,dvtdp) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                           costheta(t) * buffer(PSI,vr) * buffer(PSI,dvrdp) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,vr) * buffer(PSI,vphi) + &
                           one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,vr) * &
                           buffer(PSI,dvpdr) - &    !(second part ended)
                           one_over_r(r) * buffer(PSI,vr) * buffer(PSI,dvtdp) - &
                           one_over_r(r) * buffer(PSI,vtheta) * buffer(PSI,dvrdp) - &
                           buffer(PSI,dvrdp) * buffer(PSI,dvtdr) - &
                           buffer(PSI,vr) * DDBUFF(PSI,dvtdrdp) + &
                           buffer(PSI,dvrdp) * buffer(PSI,dvrdt) + & 
                           buffer(PSI,vr) * DDBUFF(PSI,dvrdtdp) + &
                           one_over_r(r) * buffer(PSI,dvpdp) * buffer(PSI,dvpdt) + &
                           one_over_r(r) * buffer(PSI,vphi) * DDBUFF(PSI,dvpdtdp) + &
                           2 * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,vphi) * buffer(PSI,dvpdp) - &
                           buffer(PSI,dvpdp) * buffer(PSI,dvtdp) - &               
                           buffer(PSI,vphi) * DDBUFF(PSI,dvtdpdp))**2 + &
                           (one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                           buffer(PSI,dvpdp) * buffer(PSI,dvrdp) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                           buffer(PSI,vphi) * DDBUFF(PSI,dvrdpdp) - &          
                           2 * one_over_r(r) * one_over_r(r) * csctheta(t) * &
                           buffer(PSI,vphi) * buffer(PSI,dvpdp) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI,dvpdp) * buffer(PSI,dvpdr) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI,vphi) * DDBUFF(PSI,dvpdrdp) - &
                           2 * one_over_r(r) * one_over_r(r) * csctheta(t) * &
                           buffer(PSI,vtheta) * buffer(PSI,dvtdp) - &    
                           one_over_r(r) * csctheta(t) * buffer(PSI,dvtdp) * buffer(PSI,dvtdr) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI,vtheta) * DDBUFF(PSI,dvtdrdp) + &
                           one_over_r(r) * csctheta(t) * buffer(PSI,dvtdp) * buffer(PSI,dvrdt) + &
                           one_over_r(r) * csctheta(t) * buffer(PSI,vtheta) * DDBUFF(PSI,dvrdtdp) - &    !(first part ended)
                           one_over_r(r) * one_over_r(r) * buffer(PSI,vtheta) * buffer(PSI,dvpdt) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,vtheta) * buffer(PSI,vphi) + &
                           one_over_r(r) * buffer(PSI,vtheta) * buffer(PSI,dvtdp) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,vr) * &
                           buffer(PSI,dvrdp) - &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,vr) * buffer(PSI,vphi) - &
                           one_over_r(r) * buffer(PSI,vr) * buffer(PSI,dvpdr) + &                            !(second part ended)
                           one_over_r(r) * one_over_r(r) * buffer(PSI,vtheta) * buffer(PSI,dvpdt) - &
                           one_over_r(r) * buffer(PSI,dvtdr) * buffer(PSI,dvpdt) - &
                           one_over_r(r) * buffer(PSI,vtheta) * DDBUFF(PSI,dvpdrdt) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,vtheta) * buffer(PSI,vphi) - &
                           one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,dvtdr) * buffer(PSI,vphi) - &
                           one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,vtheta) * buffer(PSI,dvpdr) + &
                           buffer(PSI,dvtdr) * buffer(PSI,dvtdp) + &
                           buffer(PSI,vtheta) * DDBUFF(PSI,dvtdrdp) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * &
                           buffer(PSI,vr) * buffer(PSI,dvrdp) + &
                           one_over_r(r) * csctheta(t) * buffer(PSI,dvrdr) * &
                           buffer(PSI,dvrdp) + &
                           one_over_r(r) * csctheta(t) * buffer(PSI,vr) * &
                           DDBUFF(PSI,dvrdrdp) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,vr) * buffer(PSI,vphi) - &
                           one_over_r(r) * buffer(PSI,dvrdr) * buffer(PSI,vphi) - &
                           one_over_r(r) * buffer(PSI,vr) * buffer(PSI,dvpdr) - &
                           buffer(PSI,dvrdr) * buffer(PSI,dvpdr) - &
                           buffer(PSI,vr) * DDBUFF(PSI,dvpdrdr))**2 + &
                           (one_over_r(r) * one_over_r(r) * buffer(PSI,vtheta) * buffer(PSI,vr) + &
                           one_over_r(r) * buffer(PSI,vr) * buffer(PSI,dvtdr) - &
                           one_over_r(r) * buffer(PSI,vr) * buffer(PSI,dvrdt) - &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,vphi) * buffer(PSI,dvpdt) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           (buffer(PSI,vphi))**2 + &
                           one_over_r(r) * buffer(PSI,vphi) * buffer(PSI,dvtdp) - &     !(first part ended)
                           one_over_r(r) * one_over_r(r) * buffer(PSI,vtheta) * &
                           buffer(PSI,vr) + &
                           one_over_r(r) * buffer(PSI,dvtdr) * buffer(PSI,vr) + &
                           one_over_r(r) * buffer(PSI,vtheta) * buffer(PSI,dvrdr) + &
                           buffer(PSI,dvrdr) * buffer(PSI,dvtdr) + &
                           buffer(PSI,vr) * DDBUFF(PSI,dvtdrdr) - &
                           buffer(PSI,dvrdr) * buffer(PSI,dvrdt) - &
                           buffer(PSI,vr) * DDBUFF(PSI,dvrdrdt) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,vphi) * buffer(PSI,dvpdt) - &
                           one_over_r(r) * buffer(PSI,dvpdr) * buffer(PSI,dvpdt) - &
                           one_over_r(r) * buffer(PSI,vphi) * DDBUFF(PSI,dvpdrdt) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           (buffer(PSI,vphi))**2 - &
                           2 * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI, vphi) * &
                           buffer(PSI,dvpdr) + &
                           buffer(PSI, dvpdr) * buffer(PSI, dvtdp) + &
                           buffer(PSI, vphi) * DDBUFF(PSI,dvtdrdp) + &                 !(second part ended)
                           one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                           buffer(PSI, vphi) * buffer(PSI,dvrdp) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI,dvpdt) * buffer(PSI,dvrdp) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI,vphi) * &
                           DDBUFF(PSI,dvrdtdp) + & 
                           2 * one_over_r(r) * buffer(PSI,vphi) * buffer(PSI,dvpdt) + &
                           buffer(PSI,dvpdt) * buffer(PSI,dvpdr) + &
                           buffer(PSI,vphi) * DDBUFF(PSI,dvpdrdt) + &
                           2 * one_over_r(r) * buffer(PSI,vtheta) * buffer(PSI,dvtdt) + &
                           buffer(PSI,dvtdt) * buffer(PSI,dvtdr) + &
                           buffer(PSI,vtheta) * DDBUFF(PSI,dvtdrdt) - &
                           buffer(PSI,dvtdt) * buffer(PSI,dvrdt) - &
                           buffer(PSI,vtheta) * DDBUFF(PSI,dvrdtdt))**2)**(0.5)
            END_DO
            Call Add_Quantity(qty)
        Endif


    End Subroutine Compute_Curl_Advection_Force

    Subroutine Compute_Curl_Buoyancy_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t
        
        !!!!!!!!!!!!!!!!!!!!!!!!!! Buoyancy Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  

        If (compute_quantity(curl_buoyancy_force_theta) .or. compute_quantity(curl_buoyancy_force_theta_squared)) Then
            DO_PSI
                qty(PSI) = ref%Buoyancy_Coeff(r) * (csctheta(t) * &
                            radius(r) * buffer(PSI,dtdp))  ! since dtdp = (1/r)*d(temperature or entropy)/dphi
            END_DO
            If (compute_quantity(curl_buoyancy_force_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_buoyancy_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif

        If (compute_quantity(curl_buoyancy_force_phi) .or. compute_quantity(curl_buoyancy_force_phi_squared)) Then
            DO_PSI
                qty(PSI) = -ref%Buoyancy_Coeff(r) * ( radius(r) * buffer(PSI,dtdt))   ! since dtdt = (1/r)*d(temperature or $
            END_DO
            If (compute_quantity(curl_buoyancy_force_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_buoyancy_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif

        Endif
        
        If (compute_quantity(curl_buoyancy_force_abs)) Then
            DO_PSI
                qty(PSI) = ((ref%Buoyancy_Coeff(r) * csctheta(t) * buffer(PSI,dtdp))**2 + &
                           (-ref%Buoyancy_Coeff(r) * buffer(PSI,dtdt))**2)**(0.5)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Curl_Buoyancy_Force

    Subroutine Compute_Curl_Magnetic_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t         

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Magnetic Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        If (compute_quantity(curl_j_cross_b_r) .or. compute_quantity(curl_j_cross_b_r_squared)) Then
            DO_PSI
                qty(PSI) = ref%Lorentz_Coeff*(one_over_r(r) * one_over_r(r) * buffer(PSI,dbtdt) * buffer(PSI,dbpdt) + &
                            one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * DDBUFF(PSI,dbpdtdt) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * buffer(PSI,btheta) * &
                            buffer(PSI,bphi) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,dbtdt) * &
                            buffer(PSI,bphi) + &                         
                            one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,btheta) * &
                            buffer(PSI,dbpdt) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                                buffer(PSI,dbtdt) * buffer(PSI,dbtdp) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                                buffer(PSI,btheta) * DDBUFF(PSI,dbtdtdp) + & 
                            one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * buffer(PSI,btheta) * &
                            buffer(PSI,dbtdp) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                            buffer(PSI,br) * buffer(PSI,dbrdp) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,dbrdt) * buffer(PSI,dbrdp) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,br) * DDBUFF(PSI,dbrdtdp) + &
                            one_over_r(r) * one_over_r(r) * buffer(PSI,dbrdt) * buffer(PSI,bphi) + &
                            one_over_r(r) * one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbpdt) + &
                            one_over_r(r) * buffer(PSI,dbrdt) * buffer(PSI,dbpdr) + &
                            one_over_r(r) * buffer(PSI,br) * DDBUFF(PSI,dbpdrdt) + &                             !(first part ended)
                            one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                            buffer(PSI,btheta) * buffer(PSI,dbpdt) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * costheta(t) * &
                            buffer(PSI,btheta) * buffer(PSI,bphi) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                            buffer(PSI,btheta) * buffer(PSI,dbtdp) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                            buffer(PSI,br) * buffer(PSI,dbrdp) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                            buffer(PSI,br) * buffer(PSI,bphi) + &
                            one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,br) * buffer(PSI,dbpdr) - &    !(second part ended)
                            one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbtdp) - &
                            one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbrdp) - &
                            buffer(PSI,dbrdp) * buffer(PSI,dbtdr) - &
                            buffer(PSI,br) * DDBUFF(PSI,dbtdrdp) + &
                            buffer(PSI,dbrdp) * one_over_r(r) * buffer(PSI,dbrdt) + & 
                            buffer(PSI,br) * one_over_r(r) * DDBUFF(PSI,dbrdtdp) + &
                            one_over_r(r) * buffer(PSI,dbpdp) * buffer(PSI,dbpdt) + &
                            one_over_r(r) * buffer(PSI,bphi) * DDBUFF(PSI,dbpdtdp) + &
                            2 * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,bphi) * buffer(PSI,dbpdp) - &
                            buffer(PSI,dbpdp) * buffer(PSI,dbtdp) - &               
                            buffer(PSI,bphi) * DDBUFF(PSI,dbtdpdp))
            END_DO
            If (compute_quantity(curl_j_cross_b_r)) Call Add_Quantity(qty)
            If (compute_quantity(curl_j_cross_b_r_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif            
        
        Endif
        
        
        If (compute_quantity(curl_j_cross_b_theta) .or. compute_quantity(curl_j_cross_b_theta_squared)) Then
            DO_PSI
                qty(PSI) = ref%Lorentz_Coeff*(one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                            buffer(PSI,dbpdp) * buffer(PSI,dbrdp) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                            buffer(PSI,bphi) * DDBUFF(PSI,dbrdpdp) - &          
                            2 * one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,bphi) * buffer(PSI,dbpdp) - &
                            one_over_r(r) * csctheta(t) * buffer(PSI,dbpdp) * buffer(PSI,dbpdr) - &
                            one_over_r(r) * csctheta(t) * buffer(PSI,bphi) * DDBUFF(PSI,dbpdrdp) - &
                            2 * one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * buffer(PSI,dbtdp) - &    
                            one_over_r(r) * csctheta(t) * buffer(PSI,dbtdp) * buffer(PSI,dbtdr) - &
                            one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * DDBUFF(PSI,dbtdrdp) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,dbtdp) * buffer(PSI,dbrdt) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * DDBUFF(PSI,dbrdtdp) - &    !(first part ended)
                            one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbpdt) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                            buffer(PSI,btheta) * buffer(PSI,bphi) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * buffer(PSI,dbtdp) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * &
                            buffer(PSI,br) * buffer(PSI,dbrdp) - &
                            one_over_r(r) * one_over_r(r) * buffer(PSI,br) * buffer(PSI,bphi) - &
                            one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbpdr) + &                            !(second part ended)
                            one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbpdt) - &
                            one_over_r(r) * buffer(PSI,dbtdr) * buffer(PSI,dbpdt) - &
                            one_over_r(r) * buffer(PSI,btheta) * DDBUFF(PSI,dbpdrdt) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                            buffer(PSI,btheta) * buffer(PSI,bphi) - &
                            one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,dbtdr) * buffer(PSI,bphi) - &
                            one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,btheta) * buffer(PSI,dbpdr) + &
                            one_over_r(r) * csctheta(t) * buffer(PSI,dbtdr) * buffer(PSI,dbtdp) - &
                            buffer(PSI,btheta) * one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,dbtdr) + &
                            one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * DDBUFF(PSI,dbtdrdp) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,br) * buffer(PSI,dbrdp) + &
                            one_over_r(r) * csctheta(t) * buffer(PSI,dbrdr) * buffer(PSI,dbrdp) + &
                            one_over_r(r) * csctheta(t) * buffer(PSI,br) * DDBUFF(PSI,dbrdrdp) + &
                            one_over_r(r) * one_over_r(r) * buffer(PSI,br) * buffer(PSI,bphi) - &
                            one_over_r(r) * buffer(PSI,dbrdr) * buffer(PSI,bphi) - &
                            one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbpdr) - &
                            buffer(PSI,dbrdr) * buffer(PSI,dbpdr) - &
                            buffer(PSI,br) * DDBUFF(PSI,dbpdrdr))
            END_DO
            If (compute_quantity(curl_j_cross_b_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_j_cross_b_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif            

        Endif
        
        If (compute_quantity(curl_j_cross_b_phi) .or. compute_quantity(curl_j_cross_b_phi_squared)) Then
            DO_PSI
                qty(PSI) = ref%Lorentz_Coeff*(one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,br) + &
                            one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbtdr) - &
                            one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbrdt) - &
                            one_over_r(r) * one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbpdt) - &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * (buffer(PSI, bphi))**2 + &
                            one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbtdp) - &     !(first part ended)
                            one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,br) + &
                            one_over_r(r) * buffer(PSI,dbtdr) * buffer(PSI,br) + &
                            one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbrdr) + &
                            buffer(PSI,dbrdr) * buffer(PSI,dbtdr) + &
                            buffer(PSI,br) * DDBUFF(PSI,dbtdrdr) - &
                            buffer(PSI,dbrdr) * buffer(PSI,dbrdt) - &
                            one_over_r(r) * buffer(PSI,br) * DDBUFF(PSI,dbrdrdt) + &
                            one_over_r(r) * one_over_r(r) * buffer(PSI,dbrdt) + &
                            one_over_r(r) * one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbpdt) - &
                            one_over_r(r) * buffer(PSI,dbpdr) * buffer(PSI,dbpdt) - &
                            one_over_r(r) * buffer(PSI,bphi) * DDBUFF(PSI,dbpdrdt) + &
                            one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                            (buffer(PSI, bphi))**2 - &
                            2 * one_over_r(r) * csctheta(t) * costheta(t) * &
                            buffer(PSI, bphi) * buffer(PSI,dbpdr) + &
                            one_over_r(r) * csctheta(t) * buffer(PSI, dbpdr) * buffer(PSI, dbtdp) - &
                            buffer(PSI, bphi) * one_over_r(r) * one_over_r(r) * csctheta(t)  * buffer(PSI, dbtdp) + &
                            one_over_r(r) * csctheta(t) * buffer(PSI, bphi) * DDBUFF(PSI,dbtdrdp) + &                 !(second part ended)
                            one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                            buffer(PSI, bphi) * buffer(PSI,dbrdp) - &
                            one_over_r(r) * csctheta(t) * buffer(PSI,dbpdt) * buffer(PSI,dbrdp) - &
                            one_over_r(r) * csctheta(t) * buffer(PSI,bphi) * DDBUFF(PSI,dbrdtdp) + & 
                            2 * one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbpdt) + &
                            buffer(PSI,dbpdt) * buffer(PSI,dbpdr) + &
                            buffer(PSI,bphi) * DDBUFF(PSI,dbpdrdt) + &
                            2 * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbtdt) + &
                            buffer(PSI,dbtdt) * buffer(PSI,dbtdr) + &
                            buffer(PSI,btheta) * DDBUFF(PSI,dbtdrdt) - &
                            one_over_r(r) * buffer(PSI,dbtdt) * buffer(PSI,dbrdt) - &
                            one_over_r(r) * buffer(PSI,btheta) * DDBUFF(PSI,dbrdtdt))
            END_DO
            If (compute_quantity(curl_j_cross_b_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_j_cross_b_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif            

        Endif
        
        If (compute_quantity(curl_j_cross_b_abs)) Then
            DO_PSI
                qty(PSI) = ((ref%Lorentz_Coeff * (one_over_r(r) * one_over_r(r) * &
                           buffer(PSI,dbtdt) * buffer(PSI,dbpdt) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * DDBUFF(PSI,dbpdtdt) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                           buffer(PSI,btheta) * buffer(PSI,bphi) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,dbtdt) * buffer(PSI,bphi) + &                         
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,btheta) * buffer(PSI,dbpdt) - &
                           one_over_r(r) * buffer(PSI,dbtdt) * buffer(PSI,dbtdp) - &
                           one_over_r(r) * buffer(PSI,btheta) * DDBUFF(PSI,dbtdtdp) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                           buffer(PSI,br) * buffer(PSI,dbrdp) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,dbrdt) * buffer(PSI,dbrdp) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,br) * DDBUFF(PSI,dbrdtdp) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,dbrdt) * buffer(PSI,bphi) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbpdt) + &
                           one_over_r(r) * buffer(PSI,dbrdt) * buffer(PSI,dbpdr) + &
                           one_over_r(r) * buffer(PSI,br) * DDBUFF(PSI,dbpdrdt) + &                             !(first part ended)
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,btheta) * buffer(PSI,dbpdt) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * costheta(t) * &
                           buffer(PSI,btheta) * buffer(PSI,bphi) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,btheta) * buffer(PSI,dbtdp) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                           buffer(PSI,br) * buffer(PSI,dbrdp) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,br) * buffer(PSI,bphi) + &
                           one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,br) * buffer(PSI,dbpdr) - &    !(second part ended)
                           one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbtdp) - &
                           one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbrdp) - &
                           buffer(PSI,dbrdp) * buffer(PSI,dbtdr) - &
                           buffer(PSI,br) * DDBUFF(PSI,dbtdrdp) + &
                           buffer(PSI,dbrdp) * buffer(PSI,dbrdt) + & 
                           buffer(PSI,br) * DDBUFF(PSI,dbrdtdp) + &
                           one_over_r(r) * buffer(PSI,dbpdp) * buffer(PSI,dbpdt) + &
                           one_over_r(r) * buffer(PSI,bphi) * DDBUFF(PSI,dbpdtdp) + &
                           2 * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,bphi) * buffer(PSI,dbpdp) - &
                           buffer(PSI,dbpdp) * buffer(PSI,dbtdp) - &               
                           buffer(PSI,bphi) * DDBUFF(PSI,dbtdpdp)))**2 + &
                           (ref%Lorentz_Coeff * (one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                           buffer(PSI,dbpdp) * buffer(PSI,dbrdp) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * csctheta(t) * &
                           buffer(PSI,bphi) * DDBUFF(PSI,dbrdpdp) - &          
                           2 * one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,bphi) * buffer(PSI,dbpdp) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI, dbpdp) * buffer(PSI,dbpdr) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI,bphi) * DDBUFF(PSI,dbpdrdp) - &
                           2 * one_over_r(r) * one_over_r(r) * csctheta(t) * &
                           buffer(PSI,btheta) * buffer(PSI,dbtdp) - &    
                           one_over_r(r) * csctheta(t) * buffer(PSI,dbtdp) * buffer(PSI,dbtdr) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * DDBUFF(PSI,dbtdrdp) + &
                           one_over_r(r) * csctheta(t) * buffer(PSI,dbtdp) * buffer(PSI,dbrdt) + &
                           one_over_r(r) * csctheta(t) * buffer(PSI,btheta) * DDBUFF(PSI,dbrdtdp) - &    !(first part ended)
                           one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbpdt) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,btheta) * buffer(PSI,bphi) + &
                           one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbtdp) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * buffer(PSI,br) * buffer(PSI,dbrdp) - &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,br) * buffer(PSI,bphi) - &
                           one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbpdr) + &                            !(second part ended)
                           one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbpdt) - &
                           one_over_r(r) * buffer(PSI,dbtdr) * buffer(PSI,dbpdt) - &
                           one_over_r(r) * buffer(PSI,btheta) * DDBUFF(PSI,dbpdrdt) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,btheta) * buffer(PSI,bphi) - &
                           one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI,dbtdr) * buffer(PSI,bphi) - &
                           one_over_r(r) * csctheta(t) * costheta(t) * &
                           buffer(PSI,btheta) * buffer(PSI,dbpdr) + &
                           buffer(PSI,dbtdr) * buffer(PSI,dbtdp) + &
                           buffer(PSI,btheta) * DDBUFF(PSI,dbtdrdp) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * &
                           buffer(PSI,br) * buffer(PSI,dbrdp) + &
                           one_over_r(r) * csctheta(t) * buffer(PSI,dbrdr) * buffer(PSI,dbrdp) + &
                           one_over_r(r) * csctheta(t) * buffer(PSI,br) * DDBUFF(PSI,dbrdrdp) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,br) * buffer(PSI,bphi) - &
                           one_over_r(r) * buffer(PSI,dbrdr) * buffer(PSI,bphi) - &
                           one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbpdr) - &
                           buffer(PSI,dbrdr) * buffer(PSI,dbpdr) - &
                           buffer(PSI,br) * DDBUFF(PSI,dbpdrdr)))**2 + &
                           (ref%Lorentz_Coeff * (one_over_r(r) * one_over_r(r) * &
                           buffer(PSI,btheta) * buffer(PSI,br) + &
                           one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbtdr) - &
                           one_over_r(r) * buffer(PSI,br) * buffer(PSI,dbrdt) - &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbpdt) - &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * &
                           (buffer(PSI, bphi))**2 + &
                           one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbtdp) - &     !(first part ended)
                           one_over_r(r) * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,br) + &
                           one_over_r(r) * buffer(PSI,dbtdr) * buffer(PSI,br) + &
                           one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbrdr) + &
                           buffer(PSI,dbrdr) * buffer(PSI,dbtdr) + &
                           buffer(PSI,br) * DDBUFF(PSI,dbtdrdr) - &
                           buffer(PSI,dbrdr) * buffer(PSI,dbrdt) - &
                           buffer(PSI,br) * DDBUFF(PSI,dbrdrdt) + &
                           one_over_r(r) * one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbpdt) - &
                           one_over_r(r) * buffer(PSI,dbpdr) * buffer(PSI,dbpdt) - &
                           one_over_r(r) * buffer(PSI,bphi) * DDBUFF(PSI,dbpdrdt) + &
                           one_over_r(r) * one_over_r(r) * csctheta(t) * costheta(t) * (buffer(PSI, bphi))**2 - &
                           2 * one_over_r(r) * csctheta(t) * costheta(t) * buffer(PSI, bphi) * buffer(PSI,dbpdr) + &
                           buffer(PSI, dbpdr) * buffer(PSI, dbtdp) + &
                           buffer(PSI, bphi) * DDBUFF(PSI,dbtdrdp) + &                 !(second part ended)
                           one_over_r(r) * csctheta(t) * csctheta(t) * costheta(t) * &
                           buffer(PSI, bphi) * buffer(PSI,dbrdp) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI,dbpdt) * buffer(PSI,dbrdp) - &
                           one_over_r(r) * csctheta(t) * buffer(PSI,bphi) * DDBUFF(PSI,dbrdtdp) + & 
                           2 * one_over_r(r) * buffer(PSI,bphi) * buffer(PSI,dbpdt) + &
                           buffer(PSI,dbpdt) * buffer(PSI,dbpdr) + &
                           buffer(PSI,bphi) * DDBUFF(PSI,dbpdrdt) + &
                           2 * one_over_r(r) * buffer(PSI,btheta) * buffer(PSI,dbtdt) + &
                           buffer(PSI,dbtdt) * buffer(PSI,dbtdr) + &
                           buffer(PSI,btheta) * DDBUFF(PSI,dbtdrdt) - &
                           buffer(PSI,dbtdt) * buffer(PSI,dbrdt) - &
                           buffer(PSI,btheta) * DDBUFF(PSI,dbrdtdt)))**2)**(0.5)
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Curl_Magnetic_Force

    Subroutine Compute_Curl_Coriolis_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t

        !!!!!!!!!!!!!!!!!!!!!!!!!!!! Coriolis Force !!!!!!!!!!!!!!!!!!!!!!!!!!

        If (compute_quantity(curl_coriolis_force_r) .or. compute_quantity(curl_coriolis_force_r_squared)) Then
            DO_PSI
                qty(PSI) = - ref%Coriolis_Coeff * ref%density(r) * one_over_r(r) * &
                            (-sintheta(t) * buffer(PSI,vtheta) + cottheta(t) * csctheta(t) * buffer(PSI,vtheta) + &
                            costheta(t) *  buffer(PSI,dvtdt) + 2 * costheta(t) *  buffer(PSI,vr) + sintheta(t) * &
                            buffer(PSI,dvrdt) + cottheta(t) *  buffer(PSI,dvpdt))    
            END_DO
            If (compute_quantity(curl_coriolis_force_r)) Call Add_Quantity(qty)
            If (compute_quantity(curl_coriolis_force_r_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        
        Endif

        If (compute_quantity(curl_coriolis_force_theta) .or. compute_quantity(curl_coriolis_force_theta_squared)) Then
            DO_PSI
                qty(PSI) = ref%Coriolis_Coeff * ref%density(r) * (one_over_r(r) * (buffer(PSI,dvpdp) + &
                            costheta(t) *  buffer(PSI,vtheta) + sintheta(t) * buffer(PSI,vr)) + costheta(t) * &
                            buffer(PSI,dvtdr) + sintheta(t) * buffer(PSI,dvrdr) + ref%dlnrho(r) * costheta(t) * &
                            buffer(PSI,vtheta) + ref%dlnrho(r) * sintheta(t) * buffer(PSI,vr))
            END_DO
            If (compute_quantity(curl_coriolis_force_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_coriolis_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
                
        Endif

        If (compute_quantity(curl_coriolis_force_phi) .or. compute_quantity(curl_coriolis_force_phi)) Then
            DO_PSI
                qty(PSI) = ref%Coriolis_Coeff * ref%density(r) * (ref%dlnrho(r) * costheta(t) * buffer(PSI,vphi) + &
                            costheta(t) * buffer(PSI,dvpdr) - one_over_r(r) * sintheta(t) * buffer(PSI,dvpdr)) 

            END_DO
            If (compute_quantity(curl_coriolis_force_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_coriolis_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        
        Endif

        If (compute_quantity(curl_coriolis_force_abs)) Then
            DO_PSI
                qty(PSI) =((- ref%Coriolis_Coeff * ref%density(r) * one_over_r(r) * &
                            (-sintheta(t) * buffer(PSI,vtheta) + cottheta(t) * csctheta(t) * buffer(PSI,vtheta) + &
                            costheta(t) *  buffer(PSI,dvtdt) + 2 * costheta(t) *  buffer(PSI,vr) + sintheta(t) * &
                            buffer(PSI,dvrdt) + cottheta(t) *  buffer(PSI,dvpdt)))**2 + &
                           (ref%Coriolis_Coeff * ref%density(r) * (one_over_r(r) * (buffer(PSI,dvpdp) + &
                            costheta(t) *  buffer(PSI,vtheta) + sintheta(t) * buffer(PSI,vr)) + costheta(t) * &
                            buffer(PSI,dvtdr) + sintheta(t) * buffer(PSI,dvrdr) + ref%dlnrho(r) * costheta(t) * &
                            buffer(PSI,vtheta) + ref%dlnrho(r) * sintheta(t) * buffer(PSI,vr)))**2 + &
                           (ref%Coriolis_Coeff * ref%density(r) * (ref%dlnrho(r) * costheta(t) * buffer(PSI,vphi) + &
                            costheta(t) * buffer(PSI,dvpdr) - one_over_r(r) * sintheta(t) * buffer(PSI,dvpdr)))**2)**(0.5)   
            END_DO
            Call Add_Quantity(qty)
        Endif

    End Subroutine Compute_Curl_Coriolis_Force

    Subroutine Compute_Curl_Pressure_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t
        ! NOTE: pfactor is assumed constant in the derivation below
        Real*8  :: pfactor(my_r%min:my_r%max)
        pfactor(my_r%min:my_r%max) = ref%dpdr_w_term(my_r%min:my_r%max) &
                                        /ref%density(my_r%min:my_r%max)

        !!!!!!!!!!!!!!!!!!!!!!!!!!! Pressure Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       ! curl_pressure_force_r is always 0
         
        If (compute_quantity(curl_pressure_force_theta) .or. compute_quantity(curl_pressure_force_theta_squared)) Then
            DO_PSI
                qty(PSI) = pfactor(r) * &
                            One_Over_R(r) * csctheta(t) * ref%dlnrho(r) * buffer(PSI,dpdt) 
            END_DO
            If (compute_quantity(curl_pressure_force_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_pressure_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

        If (compute_quantity(curl_pressure_force_phi) .or. compute_quantity(curl_pressure_force_phi_squared)) Then
            DO_PSI
                qty(PSI) = - pfactor(r) * &
                              One_Over_R(r) * ref%dlnrho(r) * buffer(PSI,dpdt)
            END_DO
            If (compute_quantity(curl_pressure_force_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_pressure_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)
            Endif
        Endif

    End Subroutine Compute_Curl_Pressure_Force

    Subroutine Compute_Curl_Viscous_Force(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r, k, t

        !!!!!!!!!!!!!!!!!!!!!!!!!!! Viscous Force !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
        If (compute_quantity(curl_viscous_force_r) .or. compute_quantity(curl_viscous_force_r_squared)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*(VFDBUFF(PSI,dvf_p_dt) + &
                                cottheta(t)*vforce_buffer(PSI,vf_p) - &
                                csctheta(t)*VFDBUFF(PSI,dvf_t_dp))
            END_DO
            If (compute_quantity(curl_viscous_force_r)) Call Add_Quantity(qty)
            If (compute_quantity(curl_viscous_force_r_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)   
            Endif
        Endif    
            
        If (compute_quantity(curl_viscous_force_theta) .or. compute_quantity(curl_viscous_force_theta_squared)) Then
            DO_PSI
                qty(PSI) = One_Over_R(r)*(csctheta(t)*VFDBUFF(PSI,dvf_r_dp) - &
                                vforce_buffer(PSI,vf_p)) - &
                            VFDBUFF(PSI,dvf_p_dr)
            END_DO
            If (compute_quantity(curl_viscous_force_theta)) Call Add_Quantity(qty)
            If (compute_quantity(curl_viscous_force_theta_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)   
            Endif
        Endif    

        If (compute_quantity(curl_viscous_force_phi) .or. compute_quantity(curl_viscous_force_phi_squared)) Then
            DO_PSI
                qty(PSI) = VFDBUFF(PSI,dvf_t_dr) + &
                            One_Over_R(r)*(vforce_buffer(PSI,vf_t) - &
                                VFDBUFF(PSI,dvf_r_dt))
            END_DO
            If (compute_quantity(curl_viscous_force_phi)) Call Add_Quantity(qty)
            If (compute_quantity(curl_viscous_force_phi_squared)) Then
                DO_PSI
                    qty(PSI) = qty(PSI)*qty(PSI)
                END_DO
                Call Add_Quantity(qty)   
            Endif
        Endif    
        If (compute_quantity(curl_viscous_force_abs)) Then
            DO_PSI
                qty(PSI) = sqrt((One_Over_R(r)*(VFDBUFF(PSI,dvf_p_dt) + &
                                cottheta(t)*vforce_buffer(PSI,vf_p) - &
                                csctheta(t)*VFDBUFF(PSI,dvf_t_dp)))**2 + &
                                (One_Over_R(r)*(csctheta(t)*VFDBUFF(PSI,dvf_r_dp) - &
                                vforce_buffer(PSI,vf_p)) - &
                                VFDBUFF(PSI,dvf_p_dr))**2 + &
                                (VFDBUFF(PSI,dvf_t_dr) + &
                                One_Over_R(r)*(vforce_buffer(PSI,vf_t) - &
                                VFDBUFF(PSI,dvf_r_dt)))**2)
            END_DO
            Call Add_Quantity(qty)   
        Endif    

    End Subroutine Compute_Curl_Viscous_Force

    Subroutine Vforce_Derivative_Logic(check, compute_vforce_i_dj)
        ! Trigger-code logic shared between the once-at-startup buffer-sizing
        ! pass (check => Sometimes_Compute, decides which vforce fields ever
        ! need derivatives, for vfdindmap/buffer sizing) and the per-iteration
        ! recheck of whether Grad_Viscous_Force needs to run this iteration
        ! (check => Compute_Quantity).
        Implicit None
        Procedure(Quantity_Check_If) :: check
        Logical, Intent(Out) :: compute_vforce_i_dj(9,3)
        Integer :: vfoff

        compute_vforce_i_dj = .false.

        If (check(curl_viscous_force_r) .or. &
            check(curl_viscous_force_r_squared) .or. &
            check(curl_viscous_force_abs)) Then
            compute_vforce_i_dj(2,3) = .true.
            compute_vforce_i_dj(3,2) = .true.
        Endif

        If (check(curl_viscous_force_theta) .or. &
            check(curl_viscous_force_theta_squared) .or. &
            check(curl_viscous_force_abs)) Then
            compute_vforce_i_dj(1,3) = .true.
            compute_vforce_i_dj(3,1) = .true.
        Endif

        If (check(curl_viscous_force_phi) .or. &
            check(curl_viscous_force_phi_squared) .or. &
            check(curl_viscous_force_abs)) Then
            compute_vforce_i_dj(1,2) = .true.
            compute_vforce_i_dj(2,1) = .true.
        Endif

        vfoff = 3
        If (check(curl_viscous_pforce_r)) Then
            compute_vforce_i_dj(vfoff+2,3) = .true.
            compute_vforce_i_dj(vfoff+3,2) = .true.
        Endif

        If (check(curl_viscous_pforce_theta)) Then
            compute_vforce_i_dj(vfoff+1,3) = .true.
            compute_vforce_i_dj(vfoff+3,1) = .true.
        Endif

        If (check(curl_viscous_pforce_phi)) Then
            compute_vforce_i_dj(vfoff+1,2) = .true.
            compute_vforce_i_dj(vfoff+2,1) = .true.
        Endif

        vfoff = 6
        If (check(curl_viscous_mforce_r)) Then
            compute_vforce_i_dj(vfoff+2,3) = .true.
            compute_vforce_i_dj(vfoff+3,2) = .true.
        Endif

        If (check(curl_viscous_mforce_theta)) Then
            compute_vforce_i_dj(vfoff+1,3) = .true.
            compute_vforce_i_dj(vfoff+3,1) = .true.
        Endif

        If (check(curl_viscous_mforce_phi)) Then
            compute_vforce_i_dj(vfoff+1,2) = .true.
            compute_vforce_i_dj(vfoff+2,1) = .true.
        Endif

    End Subroutine Vforce_Derivative_Logic

    Function Vforce_Derivatives_Needed() result(needed)
        ! Per-iteration determination of whether Grad_Viscous_Force needs to
        ! run, via compute_quantity (this iteration's menu) rather than
        ! sometimes_compute, using the same trigger-code logic as
        ! Initialize_Grad_Viscous_Force (Vforce_Derivative_Logic).
        Implicit None
        Logical :: needed
        Logical :: compute_vforce_i_dj(9,3)

        Call Vforce_Derivative_Logic(Compute_Quantity, compute_vforce_i_dj)
        needed = any(compute_vforce_i_dj)

    End Function Vforce_Derivatives_Needed

    Subroutine Initialize_Grad_Viscous_Force()
        Implicit None
        integer :: nvfind, nvfdind, vfoff
        integer :: nvfdrfields, nvfdtfields, nvfdpfields
        integer :: vfdfcount(3,2) ! buffer sizes
        Integer :: vf_i(9) ! indices to vforce_buffer
        Logical :: compute_vforce_i_dj(9,3)
        Integer :: i, j

        dvf_r_dt = -1
        dvf_r_dp = -1
        dvf_t_dr = -1
        dvf_t_dp = -1
        dvf_p_dr = -1
        dvf_p_dt = -1
        dvfp_r_dt = -1
        dvfp_r_dp = -1
        dvfp_t_dr = -1
        dvfp_t_dp = -1
        dvfp_p_dr = -1
        dvfp_p_dt = -1
        dvfm_r_dt = -1
        dvfm_r_dp = -1
        dvfm_t_dr = -1
        dvfm_t_dp = -1
        dvfm_p_dr = -1
        dvfm_p_dt = -1

        vf_i = [vf_r, vf_t, vf_p, vfp_r, vfp_t, vfp_p, vfm_r, vfm_t, vfm_p]

        Call Vforce_Derivative_Logic(Sometimes_Compute, compute_vforce_i_dj)

        ! work out how many vf fields we'll be taking the derivative of
        nvffields = count(count(compute_vforce_i_dj, dim=2) .gt. 0)
        Allocate(vfdindmap(nvffields,4))
        vfdindmap(:,:) = -1


        ! next assign indices to vf_i and vforce derivative entries
        ! this loop is designed to make sure the new derivative fields are indexed
        ! in r, t, p order so that we can grow the buffers appropriately
        nvfdind = nvffields
        do j = 1, 3
            nvfind = 0
            do i = 1, 9
                if (count(compute_vforce_i_dj(i,:)) .gt. 0) then
                    nvfind = nvfind + 1
                    ! assign indices to vforce_buffer in first column (on first outer loop)
                    if (j .eq. 1) vfdindmap(nvfind, 1) = vf_i(i)
                    if (compute_vforce_i_dj(i,j)) then
                        nvfdind = nvfdind + 1
                        if ((i .eq. 1) .and. (j .eq. 2)) dvf_r_dt = nvfdind
                        if ((i .eq. 1) .and. (j .eq. 3)) dvf_r_dp = nvfdind
                        if ((i .eq. 2) .and. (j .eq. 1)) dvf_t_dr = nvfdind
                        if ((i .eq. 2) .and. (j .eq. 3)) dvf_t_dp = nvfdind
                        if ((i .eq. 3) .and. (j .eq. 1)) dvf_p_dr = nvfdind
                        if ((i .eq. 3) .and. (j .eq. 2)) dvf_p_dt = nvfdind
                        vfoff = 3
                        if ((i .eq. vfoff+1) .and. (j .eq. 2)) dvfp_r_dt = nvfdind
                        if ((i .eq. vfoff+1) .and. (j .eq. 3)) dvfp_r_dp = nvfdind
                        if ((i .eq. vfoff+2) .and. (j .eq. 1)) dvfp_t_dr = nvfdind
                        if ((i .eq. vfoff+2) .and. (j .eq. 3)) dvfp_t_dp = nvfdind
                        if ((i .eq. vfoff+3) .and. (j .eq. 1)) dvfp_p_dr = nvfdind
                        if ((i .eq. vfoff+3) .and. (j .eq. 2)) dvfp_p_dt = nvfdind
                        vfoff = 6
                        if ((i .eq. vfoff+1) .and. (j .eq. 2)) dvfm_r_dt = nvfdind
                        if ((i .eq. vfoff+1) .and. (j .eq. 3)) dvfm_r_dp = nvfdind
                        if ((i .eq. vfoff+2) .and. (j .eq. 1)) dvfm_t_dr = nvfdind
                        if ((i .eq. vfoff+2) .and. (j .eq. 3)) dvfm_t_dp = nvfdind
                        if ((i .eq. vfoff+3) .and. (j .eq. 1)) dvfm_p_dr = nvfdind
                        if ((i .eq. vfoff+3) .and. (j .eq. 2)) dvfm_p_dt = nvfdind
                        vfdindmap(nvfind, j+1) = nvfdind
                    endif
                endif
            enddo
        enddo

        ! work out how many of each type of derivative we're taking
        nvfdrfields = count(compute_vforce_i_dj(:,1))
        nvfdtfields = count(compute_vforce_i_dj(:,2))
        nvfdpfields = count(compute_vforce_i_dj(:,3))

        ! size the buffers at each config stage
        vfdfcount(1,1) = nvffields + nvfdrfields ! config 1a
        vfdfcount(2,1) = nvffields + nvfdrfields + nvfdtfields ! 2a
        vfdfcount(3,1) = nvffields + nvfdrfields + nvfdtfields + nvfdpfields ! 3a
        vfdfcount(3,2) = nvffields ! 3b
        vfdfcount(2,2) = nvffields ! 2b
        vfdfcount(1,2) = nvffields + nvfdrfields ! 1b

        Call d_vforce_buffer%init(field_count = vfdfcount, config = 'p3b')

        Call d_vforce_buffer%construct('p3a')

        Call d_vforce_buffer%deconstruct('p3a')

    End Subroutine Initialize_Grad_Viscous_Force

    Subroutine Grad_Viscous_Force()
        Implicit None
        Integer :: i, r, k, t

        call d_vforce_buffer%construct('p3b')
        d_vforce_buffer%config = 'p3b'
        d_vforce_buffer%p3b = 0.0d0

        ! load the fields we want to take derivatives of
        do i = 1, nvffields
            d_vforce_buffer%p3b(1:n_phi,:,:,i) = vforce_buffer(:,:,:,vfdindmap(i, 1))
        enddo

        ! transform to Fourier m space
        Call fft_to_spectral(d_vforce_buffer%p3b, rsc = .true.)

        ! reform to hybrid rlm space
        call d_vforce_buffer%reform() ! move to p2b

        ! allocate spectral buffer and transform
        call d_vforce_buffer%construct('s2b')
        call Legendre_Transform(d_vforce_buffer%p2b, d_vforce_buffer%s2b)

        ! deallocate p2b
        call d_vforce_buffer%deconstruct('p2b')
        d_vforce_buffer%config = 's2b'

        ! reform
        call d_vforce_buffer%reform() ! move to p1b

        ! do a little gymnastics with p1a and p1b
        call d_vforce_buffer%construct('p1a')
        if (chebyshev) then
            ! store chebyshev coefficients in p1a and dealias
            call gridcp%to_Spectral(d_vforce_buffer%p1b, d_vforce_buffer%p1a)
            call gridcp%dealias_buffer(d_vforce_buffer%p1a)
        else
            d_vforce_buffer%p1a = d_vforce_buffer%p1b
        end if

        d_vforce_buffer%p1b = 0.0
        d_vforce_buffer%config = 'p1a'

        ! take d_by_dr
        ! (and transform back to grid space if in Chebyshev)
        if (chebyshev) then
            do i = 1, nvffields
                if (vfdindmap(i,2) .gt. 0) then
                    call gridcp%d_by_dr_cp(i, vfdindmap(i,2), d_vforce_buffer%p1a, 1)
                endif
            enddo
            call gridcp%from_spectral(d_vforce_buffer%p1a, d_vforce_buffer%p1b)
            d_vforce_buffer%p1a = d_vforce_buffer%p1b
        else
            do i = 1, nvffields
                if (vfdindmap(i,2) .gt. 0) then
                    call d_by_dx3d3(i, vfdindmap(i,2), d_vforce_buffer%p1a,1)
                endif
            enddo
        end if

        ! moving back
        call d_vforce_buffer%deconstruct('p1b')

        ! reform and start moving back
        call d_vforce_buffer%reform() ! now in s2a

        ! take theta derivatives
        do i = 1, nvffields
            if (vfdindmap(i,3) .gt. 0) then
                call d_by_dtheta(d_vforce_buffer%s2a, i, vfdindmap(i,3))
            endif
        enddo

        call d_vforce_buffer%construct('p2a')
        call Legendre_Transform(d_vforce_buffer%s2a, d_vforce_buffer%p2a)
        call d_vforce_buffer%deconstruct('s2a')
        d_vforce_buffer%config = 'p2a'

        ! reform
        call d_vforce_buffer%reform() ! move to p3a

        ! take d_by_dphi derivatives
        do i = 1, nvffields
            if (vfdindmap(i,4) .gt. 0) then
                call d_by_dphi(d_vforce_buffer%p3a, i, vfdindmap(i,4))
            endif
        enddo

        ! transform to grid space
        call FFT_To_Physical(d_vforce_buffer%p3a, rsc=.true.)

        ! Convert sintheta*{dxdt} to dxdt
        do i = 1, nvffields
            if (vfdindmap(i,3) .gt. 0) then
                DO_PSI
                    d_vforce_buffer%p3a(PSI,vfdindmap(i,3)) = d_vforce_buffer%p3a(PSI,vfdindmap(i,3))*csctheta(t)
                END_DO
            end if
        enddo

    End Subroutine Grad_Viscous_Force
 
End Module Diagnostics_Curl_Momentum
