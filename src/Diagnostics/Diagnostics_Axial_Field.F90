#include "indices.F"

!///////////////////////////////////////////////////////////////////
!               DIAGNOSTICS_AXIAL_FIELD
!               This module computes the axial components of fields. 
!               Zonal means and fluctuations about those means are 
!               also computed (if desired).
!///////////////////////////////////////////////////////////////////

Module Diagnostics_Axial_Field
    Use Diagnostics_Base
    Implicit None
Contains

    Subroutine Compute_Axial_Field(buffer)
        Implicit None
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        Integer :: r,k, t
        Real*8  :: vort_r, vort_theta
        Real*8  :: tantheta(1:n_theta)
        Real*8, Allocatable  :: tmp2(:,:,:), tmp3(:,:,:)
        tantheta(:) = 1.0d0/cottheta(:)
        Allocate(tmp2(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
        Allocate(tmp3(1:n_phi, my_r%min:my_r%max, my_theta%min:my_theta%max))
       
        !/////////////////////////////////////////
        ! 1. terms involving Axial Vorticity and Velocity
        If (compute_quantity(v_z) .or. compute_quantity(vort_z) .or. &
             compute_quantity(kin_helicity_z)) Then
            DO_PSI
                tmp1(PSI)=costheta(t)*buffer(PSI,vr)-sintheta(t)*buffer(PSI,vtheta) !v_z
                vort_r=One_Over_R(r)*( buffer(PSI,dvpdt) + &
                           cottheta(t)*buffer(PSI,vphi) - &
                           csctheta(t)*buffer(PSI,dvtdp) )                          !vort_r 
                vort_theta= One_Over_R(r)*( csctheta(t)*buffer(PSI,dvrdp) - &
                           buffer(PSI,vphi) )-buffer(PSI,dvpdr)                     !vort_theta
                qty(PSI)=costheta(t)*vort_r-sintheta(t)*vort_theta !vort_z          !vort_z
            END_DO
            
            If (compute_quantity(v_z)) Then
                Call Add_Quantity(tmp1)
            Endif 
            If (compute_quantity(vort_z)) Then
                Call Add_Quantity(qty)
            Endif 
            If (compute_quantity(kin_helicity_z)) Then
                qty = tmp1*qty
                Call Add_Quantity(qty)
            Endif
        Endif	

        If (compute_quantity(vm_z) .or. compute_quantity(vortm_z) .or. compute_quantity(kin_helicity_z_mm) .or.&
            compute_quantity(kin_helicity_z_mp) .or. compute_quantity(kin_helicity_z_pm)) Then
            DO_PSI
                tmp1(PSI)=costheta(t)*m0_values(PSI2,vr)-sintheta(t)*m0_values(PSI2,vtheta) !vm_z
                vort_r=One_Over_R(r)*( m0_values(PSI2,dvpdt) + &
                           cottheta(t)*m0_values(PSI2,vphi) )                               !vortm_r
                vort_theta= -One_Over_R(r)*m0_values(PSI2,vphi)-m0_values(PSI2,dvpdr)       !vortm_theta 
                tmp2(PSI)=costheta(t)*vort_r-sintheta(t)*vort_theta                         !vortm_z
            END_DO
            
            If (compute_quantity(vm_z)) Then
                Call Add_Quantity(tmp1)
            Endif 
            If (compute_quantity(vortm_z)) Then
                Call Add_Quantity(tmp2)
            Endif 
            If (compute_quantity(kin_helicity_z_mm)) Then
                qty = tmp1*tmp2
                Call Add_Quantity(qty)
            Endif            
            !Remember tmp1/tmp2 is currently vm_z/vortm_z
        Endif
        If (compute_quantity(vp_z) .or. compute_quantity(vortp_z) .or. compute_quantity(kin_helicity_z_pp) .or.&
            compute_quantity(kin_helicity_z_mp) .or. compute_quantity(kin_helicity_z_pm)) Then
            DO_PSI
                tmp3(PSI)=costheta(t)*fbuffer(PSI,vr)-sintheta(t)*fbuffer(PSI,vtheta) !vp_z
                vort_r=One_Over_R(r)*( fbuffer(PSI,dvpdt) + &
                           cottheta(t)*fbuffer(PSI,vphi) - &
                           csctheta(t)*fbuffer(PSI,dvtdp) )                           !vortp_r 
                vort_theta= One_Over_R(r)*( csctheta(t)*fbuffer(PSI,dvrdp) - &
                           fbuffer(PSI,vphi) )-fbuffer(PSI,dvpdr)                     !vortp_theta
                tmp4(PSI)=costheta(t)*vort_r-sintheta(t)*vort_theta                   !vortp_z
            END_DO
            
            If (compute_quantity(vp_z)) Then
                Call Add_Quantity(tmp3)
            Endif 
            If (compute_quantity(vortp_z)) Then
                Call Add_Quantity(tmp4)
            Endif 
            If (compute_quantity(kin_helicity_z_pm)) Then
                qty = tmp3*tmp2
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(kin_helicity_z_mp)) Then
                qty = tmp1*tmp4
                Call Add_Quantity(qty)
            Endif
            If (compute_quantity(kin_helicity_z_pp)) Then
                qty = tmp3*tmp4
                Call Add_Quantity(qty)
            Endif
        Endif	

        !/////////////////////////////////////////////////
        ! 2. The term involving dvzdz
        If (compute_quantity(dvzdz)) Then
            DO_PSI
                qty(PSI)=buffer(PSI,dvrdr)-buffer(PSI,dvtdr)*tantheta(t)-cottheta(t)*one_over_r(r)*&
                (buffer(PSI,dvrdt)-buffer(PSI,vtheta))+one_over_r(r)*(buffer(PSI,vr)+buffer(PSI,dvtdt)) !dvzdz
            END_DO
            Call Add_Quantity(qty)
        Endif	

        If (compute_quantity(dvzdz_m)) Then
            DO_PSI
                qty(PSI)=m0_values(PSI2,dvrdr)-m0_values(PSI2,dvtdr)*tantheta(t)-cottheta(t)*one_over_r(r)*&
                (m0_values(PSI2,dvrdt)-m0_values(PSI2,vtheta))+one_over_r(r)*(m0_values(PSI2,vr)+&
                m0_values(PSI2,dvtdt)) !dvzdz_m
            END_DO
            Call Add_Quantity(qty)
        Endif	
        If (compute_quantity(dvzdz_p)) Then
            DO_PSI
                qty(PSI)=fbuffer(PSI,dvrdr)-fbuffer(PSI,dvtdr)*tantheta(t)-cottheta(t)*one_over_r(r)*&
                (fbuffer(PSI,dvrdt)-fbuffer(PSI,vtheta))+one_over_r(r)*(fbuffer(PSI,vr)+fbuffer(PSI,dvtdt)) !dvzdz_p
            END_DO
            Call Add_Quantity(qty)
        Endif	
        !///////////////////////////////////////////
        ! 3. terms involving B_z/J_z/dTvardz/dPdz/
     
        !///////////////////////////////////////////
        ! 4. terms involving cylindrical radius. have fun.
       
        Deallocate(tmp2, tmp3)

    End Subroutine Compute_Axial_Field

End Module Diagnostics_Axial_Field
