#include "indices.F"
!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_CURRENT_DENSITY
!               This module computes the components of del x B. 
!               Zonal means and fluctuations about those means are 
!               also computed (if desired).
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Current_Density
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_J_Components(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t

        !/////////////////////////////////////////
        ! 1. terms involving J_r
        If (compute_quantity(j_r)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,curlbr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(j_r_sq) .or. compute_quantity(j_sq)) Then 
            DO_PSI
                qty(PSI) = buffer(PSI,curlbr)**2
            END_DO
            If (compute_quantity(j_r_sq)) Call Add_Quantity(qty)
            If (compute_quantity(j_sq)) tmp1 = qty
        Endif


        !/////////////////////////////////////////
        ! 2. terms involving J_theta
        If (compute_quantity(j_theta)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,curlbtheta)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(j_theta_sq) .or. compute_quantity(j_sq)) Then 
            DO_PSI
                qty(PSI) = buffer(PSI,curlbtheta)**2
            END_DO
            If (compute_quantity(j_theta_sq)) Call Add_Quantity(qty)
            If (compute_quantity(j_sq)) tmp1 = tmp1+qty
        Endif




        !/////////////////////////////////////////
        ! 3. terms involving J_phi
        If (compute_quantity(j_phi)) Then
            qty(1:n_phi,:,:) = buffer(1:n_phi,:,:,curlbphi)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(j_phi_sq) .or. compute_quantity(j_sq)) Then 
            DO_PSI
                qty(PSI) = buffer(PSI,curlbphi)**2
            END_DO
            If (compute_quantity(j_phi_sq)) Call Add_Quantity(qty)
            If (compute_quantity(j_sq)) Then
                tmp1 = tmp1+qty
                Call Add_Quantity(tmp1)
            Endif 
        Endif


        !////////////////////////////////////////////////////
        !       Now the perturbation  terms involving J'

        !1.)  J'_r
        If (compute_quantity(jp_r)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,curlbr)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jp_r_sq) .or. compute_quantity(jp_sq)) Then 
            DO_PSI
                qty(PSI) = fbuffer(PSI,curlbr)**2
            END_DO
            If (compute_quantity(jp_r_sq)) Call Add_Quantity(qty)
            If (compute_quantity(jp_sq)) tmp1= qty
        Endif

        !2.) J'_theta
        If (compute_quantity(jp_theta)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,curlbtheta)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jp_theta_sq) .or. compute_quantity(jp_sq)) Then 
            DO_PSI
                qty(PSI) = fbuffer(PSI,curlbtheta)**2
            END_DO
            If (compute_quantity(jp_theta_sq)) Call Add_Quantity(qty)
            If (compute_quantity(jp_sq)) tmp1= tmp1+qty
        Endif

        !3.) J'_phi
        If (compute_quantity(jp_phi)) Then
            qty(1:n_phi,:,:) = fbuffer(1:n_phi,:,:,curlbphi)
            Call Add_Quantity(qty)
        Endif		
        If (compute_quantity(jp_phi_sq) .or. compute_quantity(jp_sq)) Then 
            DO_PSI
                qty(PSI) = fbuffer(PSI,curlbphi)**2
            END_DO
            If (compute_quantity(jp_phi_sq)) Call Add_Quantity(qty)
            If (compute_quantity(jp_sq)) Then
                tmp1= tmp1+qty
                Call Add_Quantity(tmp1)
            Endif 
        Endif

        !////////////////////////////////////////////////////
        !       Finally, terms involving < J >

        !1.) < J >_r
        If (compute_quantity(jm_r)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,curlbr)
            END_DO
            Call Add_Quantity(qty)
        Endif	
        If (compute_quantity(jm_r_sq) .or. compute_quantity(jm_sq)) Then 
            DO_PSI
                qty(PSI) = m0_values(PSI2,curlbr)**2
            END_DO
            If (compute_quantity(jm_r_sq)) Call Add_Quantity(qty)
            If (compute_quantity(jm_sq)) tmp1= qty
        Endif


        !2.) < J >_theta
        If (compute_quantity(jm_theta)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,curlbtheta)
            END_DO
            Call Add_Quantity(qty)
        Endif	
        If (compute_quantity(jm_theta_sq) .or. compute_quantity(jm_sq)) Then 
            DO_PSI
                qty(PSI) = m0_values(PSI2,curlbtheta)**2
            END_DO
            If (compute_quantity(jm_theta_sq)) Call Add_Quantity(qty)
            If (compute_quantity(jm_sq)) tmp1= tmp1+qty
        Endif

        !3.) < J >_phi
        If (compute_quantity(jm_phi)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,curlbphi)
            END_DO
            Call Add_Quantity(qty)
        Endif	
        If (compute_quantity(jm_phi_sq) .or. compute_quantity(jm_sq)) Then 
            DO_PSI
                qty(PSI) = m0_values(PSI2,curlbphi)**2
            END_DO
            If (compute_quantity(jm_phi_sq)) Call Add_Quantity(qty)
            If (compute_quantity(jm_sq)) Then
                tmp1= tmp1+qty
                Call Add_Quantity(tmp1)
            Endif 
        Endif

        !////////////////////////////////////////////
        ! J' dot J_mean
        If (compute_quantity(jpm_sq)) Then
            DO_PSI
                qty(PSI) = m0_values(PSI2,curlbphi)*fbuffer(PSI,curlbphi)
            END_DO
            DO_PSI
                qty(PSI) = qty(PSI)+m0_values(PSI2,curlbtheta)*fbuffer(PSI,curlbtheta)
            END_DO
            DO_PSI
                qty(PSI) = qty(PSI)+m0_values(PSI2,curlbr)*fbuffer(PSI,curlbr)
            END_DO

            Call Add_Quantity(qty)
        Endif


    End Subroutine Compute_J_Components

End Module Diagnostics_Current_Density
