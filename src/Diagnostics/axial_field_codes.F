    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Axial Outputs
    Integer, Parameter :: axial_off = turbke_offset+100 ! :OFFSET CODE: 

    Integer, Parameter :: v_z              = axial_off+1  ! :tex: $v_z$
    Integer, Parameter :: vm_z             = axial_off+2  ! :tex: $\overline{v_z}$
    Integer, Parameter :: vp_z             = axial_off+3  ! :tex: $v^\prime_z$

    Integer, Parameter :: dvzdz            = axial_off+4  ! :tex: $\frac{\partial v_z}{\partial z}$    
    Integer, Parameter :: dvzdz_m          = axial_off+5  ! :tex: $\overline{\frac{\partial v_z}{\partial z}}$
    Integer, Parameter :: dvzdz_p          = axial_off+6  ! :tex: $\frac{\partial v_z}{\partial z}^\prime$

    Integer, Parameter :: vort_z           = axial_off+7  ! :tex: $\omega_z$
    Integer, Parameter :: vortm_z          = axial_off+8  ! :tex: $\overline{\omega_z}$
    Integer, Parameter :: vortp_z          = axial_off+9  ! :tex: $\omega_z^\prime$

    Integer, Parameter :: kin_helicity_z    = axial_off+10 ! :tex: $v_z \omega_z$
    Integer, Parameter :: kin_helicity_z_mm = axial_off+11 ! :tex: $\overline{v_z} \overline{\omega_z} $
    Integer, Parameter :: kin_helicity_z_pp = axial_off+12 ! :tex: $v_z^{\prime} \omega_z^\prime $
    Integer, Parameter :: kin_helicity_z_mp = axial_off+13 ! :tex: $\overline{v_z} \omega_z^\prime $
    Integer, Parameter :: kin_helicity_z_pm = axial_off+14 ! :tex: $v_z^\prime \overline{\omega_z} $

    Integer, Parameter :: B_z              = axial_off+15  ! :tex: $B_z$
    Integer, Parameter :: Bm_z             = axial_off+16  ! :tex: $\overline{B_z}$
    Integer, Parameter :: Bp_z             = axial_off+17  ! :tex: $B_z^\prime$

    Integer, Parameter :: J_z              = axial_off+18  ! :tex: $\mathcal{J}_z$    
    Integer, Parameter :: Jm_z             = axial_off+19  ! :tex: $\overline{\mathcal{J}_z}$
    Integer, Parameter :: Jp_z             = axial_off+20  ! :tex: $\mathcal{J}_z^\prime$

    Integer, Parameter :: dTvardz          = axial_off+21  ! :tex: $\frac{\partial \Theta}{\partial z}$    
    Integer, Parameter :: dTvardz_m        = axial_off+22  ! :tex: $\overline{\frac{\partial \Theta}{\partial z}}$
    Integer, Parameter :: dTvardz_p        = axial_off+23  ! :tex: $\frac{\partial \Theta}{\partial z}^\prime$

    Integer, Parameter :: dPdz             = axial_off+24  ! :tex: $\frac{\partial P}{\partial z}$    
    Integer, Parameter :: dPdz_m           = axial_off+25  ! :tex: $\overline{\frac{\partial P}{\partial z}}$
    Integer, Parameter :: dPdz_p           = axial_off+26  ! :tex: $\frac{\partial P^\prime}{\partial z}$
