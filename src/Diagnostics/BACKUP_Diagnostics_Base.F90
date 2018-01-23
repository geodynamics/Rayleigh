#define DO_PSI Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;	Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define PSI k,r,t
#define PSI2 r,t
Module Diagnostics_Base
    !//////////////////////////////////////////////////////////
    ! This module holds common variables that may be accessed
    ! by any diagnostics routine.  These variables are primarily
    ! temporary, allocatable arrays and the output menu codes. 
    Use ProblemSize
    Use Spherical_IO
    Use Fields
    Use Math_Constants
    Use ReferenceState
    Use TransportCoefficients

    Implicit None

    !/////////////////////////////////////////////////////////
    !  Diagnostic Quantity Codes 
    !  Reynolds decompositions are often used in the outputs.
    !  As a result, some shorthand is used as follows:
    !   "m" and "< >" denote the azimuthal OR spherical mean.
    !   "p" and " ' " denote perturbations about that mean

    !////////////////////////////////////////////////////////
    ! Velocity field components and it's derivatives

    Integer, Parameter :: voffset = 0

    !------------ Field Components ----------!
    ! Full
    Integer, Parameter :: v_r      = voffset+1 ! :tex: $v_r$ 
    Integer, Parameter :: v_theta  = voffset+2 ! :tex: $v_\theta$    
    Integer, Parameter :: v_phi    = voffset+3 ! :tex: $v_\phi$         

    ! Fluctuating
    Integer, Parameter :: vp_r     = voffset+4 
    Integer, Parameter :: vp_theta = voffset+5     
    Integer, Parameter :: vp_phi   = voffset+6    

    ! Mean
    Integer, Parameter :: vm_r     = voffset+7 
    Integer, Parameter :: vm_theta = voffset+8     
    Integer, Parameter :: vm_phi   = voffset+9    

    !------------ Radial Derivatives -------------!
    ! Full
    Integer, Parameter :: dv_r_dr      = voffset+10 ! :tex: $\frac{\partial v_r}{\partial r}$
    Integer, Parameter :: dv_theta_dr  = voffset+11 ! :tex: $\frac{\partial v_\theta}{\partial r}$
    Integer, Parameter :: dv_phi_dr    = voffset+12 ! :tex: $\frac{\partial v_\phi}{\partial r}$ 

    Integer, Parameter :: dvp_r_dr     = voffset+13 ! Fluctuating 
    Integer, Parameter :: dvp_theta_dr = voffset+14  
    Integer, Parameter :: dvp_phi_dr   = voffset+15    
  
    Integer, Parameter :: dvm_r_dr     = voffset+16 ! Mean    
    Integer, Parameter :: dvm_theta_dr = voffset+17  
    Integer, Parameter :: dvm_phi_dr   = voffset+18    

    !------------ Theta Derivatives --------------!
    Integer, Parameter :: dv_r_dt      = voffset+19 ! Full 
    Integer, Parameter :: dv_theta_dt  = voffset+20  
    Integer, Parameter :: dv_phi_dt    = voffset+21   

    Integer, Parameter :: dvp_r_dt     = voffset+22 ! Fluctuating 
    Integer, Parameter :: dvp_theta_dt = voffset+23    
    Integer, Parameter :: dvp_phi_dt   = voffset+24    

    Integer, Parameter :: dvm_r_dt     = voffset+25 ! Mean 
    Integer, Parameter :: dvm_theta_dt = voffset+26   
    Integer, Parameter :: dvm_phi_dt   = voffset+27     

    !------------ Phi Derivatives ----------------!
    Integer, Parameter :: dv_r_dp      = voffset+28  ! Full
    Integer, Parameter :: dv_theta_dp  = voffset+29
    Integer, Parameter :: dv_phi_dp    = voffset+30

    Integer, Parameter :: dvp_r_dp     = voffset+31 ! Fluctuating
    Integer, Parameter :: dvp_theta_dp = voffset+32 ! (same as full if    
    Integer, Parameter :: dvp_phi_dp   = voffset+33 !  mean is azimuthal)

    Integer, Parameter :: dvm_r_dp     = voffset+34 ! We keep the phi derivatives of the mean    
    Integer, Parameter :: dvm_theta_dp = voffset+35 ! here, even though they are zero.
    Integer, Parameter :: dvm_phi_dp   = voffset+36 ! Somewhat useful placeholders...

    !------------ (1/r) * Theta Derivatives -------!
    Integer, Parameter :: dv_r_dtr      = voffset+37 ! Full 
    Integer, Parameter :: dv_theta_dtr  = voffset+38 
    Integer, Parameter :: dv_phi_dtr    = voffset+39   

    Integer, Parameter :: dvp_r_dtr     = voffset+40 ! Fluctuating
    Integer, Parameter :: dvp_theta_dtr = voffset+41
    Integer, Parameter :: dvp_phi_dtr   = voffset+42

    Integer, Parameter :: dvm_r_dtr     = voffset+43 ! Mean  
    Integer, Parameter :: dvm_theta_dtr = voffset+44
    Integer, Parameter :: dvm_phi_dtr   = voffset+45


    !------(1/{r sintheta})* Phi Derivatives ---!

    Integer, Parameter :: dv_r_dprs      = voffset+46 ! Full
    Integer, Parameter :: dv_theta_dprs  = voffset+47  
    Integer, Parameter :: dv_phi_dprs    = voffset+48

    Integer, Parameter :: dvp_r_dprs     = voffset+49 ! Fluctuating   
    Integer, Parameter :: dvp_theta_dprs = voffset+50
    Integer, Parameter :: dvp_phi_dprs   = voffset+51
 
    Integer, Parameter :: dvm_r_dprs     = voffset+52 ! Mean   
    Integer, Parameter :: dvm_theta_dprs = voffset+53
    Integer, Parameter :: dvm_phi_dprs   = voffset+54  


    !--------- Radial Second Derivatives --------
    Integer, Parameter :: dv_r_d2r      = voffset+55 ! Full 
    Integer, Parameter :: dv_theta_d2r  = voffset+56  
    Integer, Parameter :: dv_phi_d2r    = voffset+57   

    Integer, Parameter :: dvp_r_d2r     = voffset+58 ! Fluctuating 
    Integer, Parameter :: dvp_theta_d2r = voffset+59  
    Integer, Parameter :: dvp_phi_d2r   = voffset+60    
  
    Integer, Parameter :: dvm_r_d2r     = voffset+61 ! Mean    
    Integer, Parameter :: dvm_theta_d2r = voffset+62  
    Integer, Parameter :: dvm_phi_d2r   = voffset+63   

    !--------- Theta Second Derivatives --------
    Integer, Parameter :: dv_r_d2t      = voffset+64 ! Full 
    Integer, Parameter :: dv_theta_d2t  = voffset+65  
    Integer, Parameter :: dv_phi_d2t    = voffset+66   

    Integer, Parameter :: dvp_r_d2t     = voffset+67 ! Fluctuating 
    Integer, Parameter :: dvp_theta_d2t = voffset+68  
    Integer, Parameter :: dvp_phi_d2t   = voffset+69    
  
    Integer, Parameter :: dvm_r_d2t     = voffset+70 ! Mean    
    Integer, Parameter :: dvm_theta_d2t = voffset+71  
    Integer, Parameter :: dvm_phi_d2t   = voffset+72

     !--------- Phi Second Derivatives --------
    Integer, Parameter :: dv_r_d2p      = voffset+73 ! Full 
    Integer, Parameter :: dv_theta_d2p  = voffset+74  
    Integer, Parameter :: dv_phi_d2p    = voffset+75   

    Integer, Parameter :: dvp_r_d2p     = voffset+76 ! Fluctuating 
    Integer, Parameter :: dvp_theta_d2p = voffset+77  
    Integer, Parameter :: dvp_phi_d2p   = voffset+78    
  
    Integer, Parameter :: dvm_r_d2p     = voffset+79 ! Mean    
    Integer, Parameter :: dvm_theta_d2p = voffset+80  
    Integer, Parameter :: dvm_phi_d2p   = voffset+81  

    !--------- Radial-Theta Second Derivatives --------
    Integer, Parameter :: dv_r_d2rt      = voffset+82 ! Full 
    Integer, Parameter :: dv_theta_d2rt  = voffset+83  
    Integer, Parameter :: dv_phi_d2rt    = voffset+84   

    Integer, Parameter :: dvp_r_d2rt     = voffset+85 ! Fluctuating 
    Integer, Parameter :: dvp_theta_d2rt = voffset+86  
    Integer, Parameter :: dvp_phi_d2rt   = voffset+87    
  
    Integer, Parameter :: dvm_r_d2rt     = voffset+88 ! Mean    
    Integer, Parameter :: dvm_theta_d2rt = voffset+89  
    Integer, Parameter :: dvm_phi_d2rt   = voffset+90   

    !--------- Radial-phi Second Derivatives --------
    Integer, Parameter :: dv_r_d2rp      = voffset+91 ! Full 
    Integer, Parameter :: dv_theta_d2rp  = voffset+92  
    Integer, Parameter :: dv_phi_d2rp    = voffset+93   

    Integer, Parameter :: dvp_r_d2rp     = voffset+94 ! Fluctuating 
    Integer, Parameter :: dvp_theta_d2rp = voffset+95  
    Integer, Parameter :: dvp_phi_d2rp   = voffset+96    
  
    Integer, Parameter :: dvm_r_d2rp     = voffset+97 ! Mean    
    Integer, Parameter :: dvm_theta_d2rp = voffset+98  
    Integer, Parameter :: dvm_phi_d2rp   = voffset+99  

    !--------- theta-phi Second Derivatives --------
    Integer, Parameter :: dv_r_d2tp      = voffset+100 ! Full 
    Integer, Parameter :: dv_theta_d2tp  = voffset+101  
    Integer, Parameter :: dv_phi_d2tp    = voffset+102   

    Integer, Parameter :: dvp_r_d2tp     = voffset+103 ! Fluctuating 
    Integer, Parameter :: dvp_theta_d2tp = voffset+104  
    Integer, Parameter :: dvp_phi_d2tp   = voffset+105    
  
    Integer, Parameter :: dvm_r_d2tp     = voffset+106 ! Mean    
    Integer, Parameter :: dvm_theta_d2tp = voffset+107  
    Integer, Parameter :: dvm_phi_d2tp   = voffset+108 

    !------------ Mass Flux ---------------------!
    Integer, Parameter :: rhov_r      = voffset+109 ! Full 
    Integer, Parameter :: rhov_theta  = voffset+110
    Integer, Parameter :: rhov_phi    = voffset+111

    Integer, Parameter :: rhovp_r     = voffset+112 ! Fluctuating 
    Integer, Parameter :: rhovp_theta = voffset+113
    Integer, Parameter :: rhovp_phi   = voffset+114

    Integer, Parameter :: rhovm_r     = voffset+115 ! Mean
    Integer, Parameter :: rhovm_theta = voffset+116
    Integer, Parameter :: rhovm_phi   = voffset+117


    !//////////////////////////////////////////////////////////////////
    !  Pressure, Entropy or Temperature, and Their Derivatives
    !  Note:  In Boussinesq Mode, Temperature is Output Instead of Entropy
    Integer, Parameter :: pt_off = voffset+117 ! = 117

    !------------ Fields ---------------------!    
    Integer, Parameter :: entropy    = pt_off+1 ! Full
    Integer, Parameter :: pressure   = pt_off+2

    Integer, Parameter :: entropy_p  = pt_off+3 ! Fluctuating
    Integer, Parameter :: pressure_p = pt_off+4

    Integer, Parameter :: entropy_m  = pt_off+5 ! Mean
    Integer, Parameter :: pressure_m = pt_off+6

    !------------ Radial Derivatives --------------!
    Integer, Parameter :: entropy_dr     = pt_off+7  ! Full
    Integer, Parameter :: pressure_dr    = pt_off+8

    Integer, Parameter :: entropy_p_dr   = pt_off+9  ! Fluctuating
    Integer, Parameter :: pressure_p_dr  = pt_off+10

    Integer, Parameter :: entropy_m_dr   = pt_off+11 ! Mean
    Integer, Parameter :: pressure_m_dr = pt_off+12

    !------------ Theta Derivatives ---------------!
    Integer, Parameter :: entropy_dtheta     = pt_off+13 ! Full
    Integer, Parameter :: pressure_dtheta    = pt_off+14

    Integer, Parameter :: entropy_p_dtheta   = pt_off+15 ! Fluctuating
    Integer, Parameter :: pressure_p_dtheta  = pt_off+16

    Integer, Parameter :: entropy_m_dtheta   = pt_off+17 ! Mean
    Integer, Parameter :: pressure_m_dtheta  = pt_off+18

    !------------ Phi Derivatives -----------------!
    Integer, Parameter :: entropy_dphi     = pt_off+19 ! Full
    Integer, Parameter :: pressure_dphi    = pt_off+20

    Integer, Parameter :: entropy_p_dphi   = pt_off+21 ! Fluctuating
    Integer, Parameter :: pressure_p_dphi  = pt_off+22

    Integer, Parameter :: entropy_m_dphi   = pt_off+23 ! Mean
    Integer, Parameter :: pressure_m_dphi  = pt_off+24

    !------------ (1/r) * Theta Derivatives --------!
    Integer, Parameter :: entropy_dtr     = pt_off+25 ! Full
    Integer, Parameter :: pressure_dtr    = pt_off+26

    Integer, Parameter :: entropy_p_dtr   = pt_off+27 ! Fluctuating
    Integer, Parameter :: pressure_p_dtr  = pt_off+28

    Integer, Parameter :: entropy_m_dtr   = pt_off+29 ! Mean
    Integer, Parameter :: pressure_m_dtr  = pt_off+30

    !--- (1/{r sintheta}) * Phi Derivatives ---------!
    Integer, Parameter :: entropy_dprs     = pt_off+31 ! Full
    Integer, Parameter :: pressure_dprs    = pt_off+32

    Integer, Parameter :: entropy_p_dprs   = pt_off+33 ! Fluctuating
    Integer, Parameter :: pressure_p_dprs  = pt_off+34

    Integer, Parameter :: entropy_m_dprs   = pt_off+35 ! Mean
    Integer, Parameter :: pressure_m_dprs  = pt_off+36

    !------------ Radial Second Derivatives --------------!
    Integer, Parameter :: entropy_d2r     = pt_off+37  ! Full
    Integer, Parameter :: pressure_d2r    = pt_off+38

    Integer, Parameter :: entropy_p_d2r   = pt_off+39  ! Fluctuating
    Integer, Parameter :: pressure_p_d2r  = pt_off+40

    Integer, Parameter :: entropy_m_d2r   = pt_off+41 ! Mean
    Integer, Parameter :: pressure_m_d2r  = pt_off+42

    !------------ Theta Second Derivatives --------------!
    Integer, Parameter :: entropy_d2t     = pt_off+43  ! Full
    Integer, Parameter :: pressure_d2t    = pt_off+44

    Integer, Parameter :: entropy_p_d2t   = pt_off+45  ! Fluctuating
    Integer, Parameter :: pressure_p_d2t  = pt_off+46

    Integer, Parameter :: entropy_m_d2t   = pt_off+47 ! Mean
    Integer, Parameter :: pressure_m_d2t  = pt_off+48

    !------------ Phi Second Derivatives --------------!
    Integer, Parameter :: entropy_d2p     = pt_off+49  ! Full
    Integer, Parameter :: pressure_d2p    = pt_off+50

    Integer, Parameter :: entropy_p_d2p   = pt_off+51  ! Fluctuating
    Integer, Parameter :: pressure_p_d2p  = pt_off+52

    Integer, Parameter :: entropy_m_d2p   = pt_off+53 ! Mean
    Integer, Parameter :: pressure_m_d2p  = pt_off+54

    !------------ Radial-Theta Second Derivatives --------------!
    Integer, Parameter :: entropy_d2rt     = pt_off+55  ! Full
    Integer, Parameter :: pressure_d2rt    = pt_off+56

    Integer, Parameter :: entropy_p_d2rt   = pt_off+57  ! Fluctuating
    Integer, Parameter :: pressure_p_d2rt  = pt_off+58

    Integer, Parameter :: entropy_m_d2rt   = pt_off+59 ! Mean
    Integer, Parameter :: pressure_m_d2rt  = pt_off+60

    !------------ Radial-Phi Derivatives --------------!
    Integer, Parameter :: entropy_d2rp     = pt_off+61  ! Full
    Integer, Parameter :: pressure_d2rp    = pt_off+62

    Integer, Parameter :: entropy_p_d2rp   = pt_off+63  ! Fluctuating
    Integer, Parameter :: pressure_p_d2rp  = pt_off+64

    Integer, Parameter :: entropy_m_d2rp   = pt_off+65 ! Mean
    Integer, Parameter :: pressure_m_d2rp  = pt_off+66

    !------------ Theta-Phi Derivatives --------------!
    Integer, Parameter :: entropy_d2tp     = pt_off+67  ! Full
    Integer, Parameter :: pressure_d2tp    = pt_off+68

    Integer, Parameter :: entropy_p_d2tp   = pt_off+69  ! Fluctuating
    Integer, Parameter :: pressure_p_d2tp  = pt_off+70

    Integer, Parameter :: entropy_m_d2tp   = pt_off+71 ! Mean
    Integer, Parameter :: pressure_m_d2tp  = pt_off+72

    !--------- Also compute terms of the form rho_bar d_by_dr(P/rho_bar)
    Integer, Parameter :: rhopressure_dr  = pt_off+73
    Integer, Parameter :: rhopressurep_dr = pt_off+74
    Integer, Parameter :: rhopressurem_dr = pt_off+75

    !//////////////////////////////////////////////////////////////////////////
    !///////////////////////////////////////////////////
    !       Vorticity Outputs
    Integer, Parameter :: vort_off = pt_off+75 ! = 192

    Integer, Parameter :: vort_r      = vort_off+1  ! Full
    Integer, Parameter :: vort_theta  = vort_off+2
    Integer, Parameter :: vort_phi    = vort_off+3

    Integer, Parameter :: vortp_r     = vort_off+4  ! Fluctuating   
    Integer, Parameter :: vortp_theta = vort_off+5
    Integer, Parameter :: vortp_phi   = vort_off+6

    Integer, Parameter :: vortm_r     = vort_off+7  ! Mean 
    Integer, Parameter :: vortm_theta = vort_off+8
    Integer, Parameter :: vortm_phi   = vort_off+9

    Integer, Parameter :: enstrophy    = vort_off+10 ! Enstrophy
    Integer, Parameter :: enstrophy_pm = vort_off+11 ! (fluctuating-mean)
    Integer, Parameter :: enstrophy_mm = vort_off+12 ! (mean-mean)
    Integer, Parameter :: enstrophy_pp  = vort_off+13 ! (fluct-fluct)


    !//////////////////////////////////////////////////////////
    !               Radial Energy Fluxes
    Integer, Parameter :: eoffset = vort_off+ 13 ! =205
    Integer, Parameter :: ecrossb_r            = eoffset+1 ! [ExB]_r (un-normalized Poynting flux)
    Integer, Parameter :: ke_flux_radial       = eoffset+2 ! vr*KE
    Integer, Parameter :: thermalE_flux_radial = eoffset+3 ! vr*rho_bar*T_bar*S OR vr*T
    Integer, Parameter :: enth_flux_radial     = eoffset+4 ! vr*cp*rho_bar*T
    Integer, Parameter :: visc_flux_r          = eoffset+5 ! -[Div dot D]_r
    Integer, Parameter :: vol_heat_flux        = eoffset+6 ! "Flux" associated with volumetric heating
    Integer, Parameter :: cond_flux_r          = eoffset+7 ! Thermal conductive flux
    Integer, Parameter :: vol_heating          = eoffset+8 ! Volumetric Heating Function
    Integer, Parameter :: visc_heating         = eoffset+9 ! Viscous Heating


    !///////////////////////////////////////////////////////////
    !           Kinetic Energies
    Integer, Parameter :: keoffset = eoffset+9 ! = 214
    Integer, Parameter :: kinetic_energy = keoffset+1   ! 1/2 rho_bar v^2
    Integer, Parameter :: radial_ke      = keoffset+2   ! 1/2 rho_bar {v_r}^2
    Integer, Parameter :: theta_ke       = keoffset+3   ! 1/2 rho_bar {v_theta}^2
    Integer, Parameter :: phi_ke         = keoffset+4   ! 1/2 rho_bar {v_phi}^2

    Integer, Parameter :: mkinetic_energy = keoffset+5 ! 1/2 rho_bar <v>^2
    Integer, Parameter :: radial_mke      = keoffset+6 ! 1/2 rho_bar <v_r>^2
    Integer, Parameter :: theta_mke       = keoffset+7 ! 1/2 rho_bar <v_theta>^2
    Integer, Parameter :: phi_mke         = keoffset+8 ! 1/2 rho_bar <v_phi>^2

    Integer, Parameter :: pkinetic_energy = keoffset+9  ! 1/2 rho_bar {v'}^2
    Integer, Parameter :: radial_pke      = keoffset+10 ! 1/2 rho_bar {v_r'}^2
    Integer, Parameter :: theta_pke       = keoffset+11 ! 1/2 rho_bar {v_theta'}^2
    Integer, Parameter :: phi_pke         = keoffset+12 ! 1/2 rho_bar {v_phi'}^2

    !--- Since density varies with radius, it can be useful to output the 
    !--- squared fields, sans density, as well.
    Integer, Parameter :: vsq         = keoffset+13   ! v^2
    Integer, Parameter :: radial_vsq  = keoffset+14   ! {v_r}^2
    Integer, Parameter :: theta_vsq   = keoffset+15   ! {v_theta}^2
    Integer, Parameter :: phi_vsq     = keoffset+16   ! {v_phi}^2

    Integer, Parameter :: mvsq        = keoffset+17   ! <v>^2
    Integer, Parameter :: radial_mvsq = keoffset+18   ! <v_r>^2
    Integer, Parameter :: theta_mvsq  = keoffset+19   ! <v_theta>^2
    Integer, Parameter :: phi_mvsq    = keoffset+20   ! <v_phi>^2

    Integer, Parameter :: pvsq        = keoffset+21   ! {v'}^2
    Integer, Parameter :: radial_pvsq = keoffset+22   ! {v_r'}^2
    Integer, Parameter :: theta_pvsq  = keoffset+23   ! {v_theta'}^2
    Integer, Parameter :: phi_pvsq    = keoffset+24   ! {v_phi'}^2


    !--- Thermal Energies
    Integer, Parameter :: teoffset=keoffset+24 ! 238
    Integer, Parameter :: thermal_energy_full = teoffset+1  ! rho_bar T_bar S
    Integer, Parameter :: thermal_energy_p    = teoffset+2  ! rho_bar T_bar S'
    Integer, Parameter :: thermal_energy_m    = teoffset+3  ! rho_bar T_bar <S>


    !////////////////////////////////////////////////////////////////////////////////////////
    ! Momentum Equation ... 


    !////////////////////////  Advection Terms ////////////////////
    ! Reynolds decomposition about the azimuthal mean may also be output
    ! NOTE:  ADVECTION TERMS ARE SCALED BY DENSITY (so that they represent a force density)

    Integer, Parameter :: mom_eq_off = teoffset+3 ! = 241    ! Output offset for advection terms  
    Integer, Parameter :: v_grad_v_r       = mom_eq_off+1 ! radial component of v dot grad v
    Integer, Parameter :: v_grad_v_theta   = mom_eq_off+2 !  theta component of v dot grad v
    Integer, Parameter :: v_grad_v_phi     = mom_eq_off+3 !    phi component of v dot grad v

    Integer, Parameter :: vp_grad_vm_r     = mom_eq_off+4 ! radial component of v' dot grad <v>
    Integer, Parameter :: vp_grad_vm_theta = mom_eq_off+5 !  theta component of v' dot grad <v>
    Integer, Parameter :: vp_grad_vm_phi   = mom_eq_off+6 !    phi component of v' dot grad <v>

    Integer, Parameter :: vm_grad_vp_r     = mom_eq_off+7 ! radial component of <v> dot grad v'
    Integer, Parameter :: vm_grad_vp_theta = mom_eq_off+8 !  theta component of <v> dot grad v'
    Integer, Parameter :: vm_grad_vp_phi   = mom_eq_off+9 !    phi component of <v> dot grad v'

    Integer, Parameter :: vp_grad_vp_r     = mom_eq_off+10 ! radial component of v' dot grad v'
    Integer, Parameter :: vp_grad_vp_theta = mom_eq_off+11 !  theta component of v' dot grad v'
    Integer, Parameter :: vp_grad_vp_phi   = mom_eq_off+12 !    phi component of v' dot grad v'

    Integer, Parameter :: vm_grad_vm_r     = mom_eq_off+13 ! radial component of <v> dot grad <v>
    Integer, Parameter :: vm_grad_vm_theta = mom_eq_off+14 !  theta component of <v> dot grad <v>
    Integer, Parameter :: vm_grad_vm_phi   = mom_eq_off+15 !    phi component of <v> dot grad <v>

    !/////////////////////////////////////////////////////////////
    !  Linear Forces 
    !  Note:  the pressure gradient diagnostic codes are  above

    Integer, Parameter :: buoyancy_force  =  mom_eq_off+16
    Integer, Parameter :: buoyancy_pforce =  mom_eq_off+17
    Integer, Parameter :: buoyancy_mforce =  mom_eq_off+18

    Integer, Parameter :: Coriolis_Force_r      = mom_eq_off+19
    Integer, Parameter :: Coriolis_Force_theta  = mom_eq_off+20
    Integer, Parameter :: Coriolis_Force_phi    = mom_eq_off+21

    Integer, Parameter :: Coriolis_pForce_r     = mom_eq_off+22
    Integer, Parameter :: Coriolis_pForce_theta = mom_eq_off+23
    Integer, Parameter :: Coriolis_pForce_phi   = mom_eq_off+24

    Integer, Parameter :: Coriolis_mForce_r     = mom_eq_off+25
    Integer, Parameter :: Coriolis_mForce_theta = mom_eq_off+26
    Integer, Parameter :: Coriolis_mForce_phi   = mom_eq_off+27

    ! Viscous forces   
    Integer, Parameter :: viscous_Force_r       = mom_eq_off+28
    Integer, Parameter :: viscous_Force_theta   = mom_eq_off+29
    Integer, Parameter :: viscous_Force_phi     = mom_eq_off+30

    Integer, Parameter :: viscous_pForce_r      = mom_eq_off+31
    Integer, Parameter :: viscous_pForce_theta  = mom_eq_off+32
    Integer, Parameter :: viscous_pForce_phi    = mom_eq_off+33

    Integer, Parameter :: viscous_mForce_r      = mom_eq_off+34
    Integer, Parameter :: viscous_mForce_theta  = mom_eq_off+35
    Integer, Parameter :: viscous_mForce_phi    = mom_eq_off+36

    ! Pressure forces   
    Integer, Parameter :: pressure_Force_r       = mom_eq_off+37
    Integer, Parameter :: pressure_Force_theta   = mom_eq_off+38
    Integer, Parameter :: pressure_Force_phi     = mom_eq_off+39

    Integer, Parameter :: pressure_pForce_r      = mom_eq_off+40
    Integer, Parameter :: pressure_pForce_theta  = mom_eq_off+41
    Integer, Parameter :: pressure_pForce_phi    = mom_eq_off+42

    Integer, Parameter :: pressure_mForce_r      = mom_eq_off+43
    Integer, Parameter :: pressure_mForce_theta  = mom_eq_off+44
    Integer, Parameter :: pressure_mForce_phi    = mom_eq_off+45

    ! ell=0 pressure and buoyancy forces (r-direction only)
    ! These are substracted out from the radial terms above
    ! for obvious reasons.    
    Integer, Parameter :: buoyancy_force_ell0 = mom_eq_off+46
    Integer, Parameter :: pressure_force_ell0_r = mom_eq_off+47         


    !////////////////////////////////////////////////////////////////////////
    !  Thermal Energy Equation
    Integer, Parameter :: thrm_eq_off   = mom_eq_off+32 ! 169

    ! Advection terms  (all terms are scaled by rho_bar T_bar in anelastic mode)
    Integer, Parameter :: v_grad_s_r       = thrm_eq_off+1 ! radial component of v dot grad S
    Integer, Parameter :: v_grad_s_theta   = thrm_eq_off+2 !  theta component of v dot grad S
    Integer, Parameter :: v_grad_s_phi     = thrm_eq_off+3 !    phi component of v dot grad S

    Integer, Parameter :: vp_grad_sm_r     = thrm_eq_off+4 ! radial component of v' dot grad <s>
    Integer, Parameter :: vp_grad_sm_theta = thrm_eq_off+5 !  theta component of v' dot grad <s>
    Integer, Parameter :: vp_grad_sm_phi   = thrm_eq_off+6 !    phi component of v' dot grad <s>

    Integer, Parameter :: vm_grad_sp_r     = thrm_eq_off+7 ! radial component of <v> dot grad S'
    Integer, Parameter :: vm_grad_sp_theta = thrm_eq_off+8 !  theta component of <v> dot grad S'
    Integer, Parameter :: vm_grad_sp_phi   = thrm_eq_off+9 !    phi component of <v> dot grad S'

    Integer, Parameter :: vp_grad_sp_r     = thrm_eq_off+10 ! radial component of v' dot grad S'
    Integer, Parameter :: vp_grad_sp_theta = thrm_eq_off+11 !  theta component of v' dot grad S'
    Integer, Parameter :: vp_grad_sp_phi   = thrm_eq_off+12 !    phi component of v' dot grad S'

    Integer, Parameter :: vm_grad_sm_r     = thrm_eq_off+13 ! radial component of <v> dot grad <S>
    Integer, Parameter :: vm_grad_sm_theta = thrm_eq_off+14 !  theta component of <v> dot grad <S>
    Integer, Parameter :: vm_grad_sm_phi   = thrm_eq_off+15 !    phi component of <v> dot grad <S>




    !///////////////////////////////////////////////////////////////////////////////////
    !Angular Momentum Transport Diagnostics
    !  Reynolds decomposition of the azimuthally-averaged angular momentum fluxes.
    Integer, Parameter :: amoff = thrm_eq_off+ 15        ! = 238
    Integer, Parameter :: amom_fluct_r     = amoff+1 ! rho_bar * r * sintheta * {v_r'} * {v_phi'}
    Integer, Parameter :: amom_fluct_theta = amoff+2 ! rho_bar * r * sintheta * {v_theta'} * {v_phi'}
    Integer, Parameter :: amom_dr_r        = amoff+3 ! rho_bar * r * sintheta * <v_r> * <v_phi>
    Integer, Parameter :: amom_dr_theta    = amoff+4 ! rho_bar * r * sintheta * <v_theta> * <v_phi>
    Integer, Parameter :: amom_mean_r      = amoff+5 ! rho_bar * r^2 * sintheta^2 * <v_r> * Omega
    Integer, Parameter :: amom_mean_theta  = amoff+6 ! rho_bar * r^2 * sintheta^2 * <v_theta> * Omega

    ! Quantity codes for magnetic torques are defined with the Lorentz Forces



    !//////////////////////////////////////////////////////////
    !  Turbulent kinetic energy generation
    Integer, Parameter :: turbke_offset = amoff+6 ! 201


    Integer, Parameter :: production_buoyant_pKE   = turbke_offset + 1	! Buoyant Production of turbulent kinetic energy
    Integer, Parameter :: production_shear_pKE     = turbke_offset + 2	    ! Shear Production of turbulent kinetic energy
    Integer, Parameter :: dissipation_viscous_pKE  = turbke_offset + 3	! Viscous Dissipation of turbulent kinetic energy
    Integer, Parameter :: transport_pressure_pKE   = turbke_offset + 4	! Pressure Transport of turbulent kinetic energy
    Integer, Parameter :: transport_viscous_pKE    = turbke_offset + 5	    ! Viscous Transport of turbulent kinetic energy    
    Integer, Parameter :: transport_turbadvect_pKE = turbke_offset + 6	! Turbulent Advective Transport of turbulent kinetic energy
    Integer, Parameter :: transport_meanadvect_pKE = turbke_offset + 7	! Mean Advective Transport of turbulent kinetic energy
    Integer, Parameter :: rflux_pressure_pKE       = turbke_offset + 8		! Radial Pressure Flux of turbulent kinetic energy
    Integer, Parameter :: rflux_viscous_pKE        = turbke_offset + 9		    ! Radial Viscous Flux of turbulent kinetic energy    
    Integer, Parameter :: rflux_turbadvect_pKE     = turbke_offset + 10	    ! Radial Turbulent Advective Flux of turbulent kinetic energy
    Integer, Parameter :: rflux_meanadvect_pKE     = turbke_offset + 11	    ! Radial Mean Advective Flux of turbulent kinetic energy
    Integer, Parameter :: thetaflux_pressure_pKE   = turbke_offset + 12	! Colatitudinal Pressure Flux of turbulent kinetic energy
    Integer, Parameter :: thetaflux_viscous_pKE    = turbke_offset + 13	! Colatitudinal Viscous Flux of turbulent kinetic energy    
    Integer, Parameter :: thetaflux_turbadvect_pKE = turbke_offset + 14	! Colatitudinal Turbulent Advective Flux of turbulent kinetic energy
    Integer, Parameter :: thetaflux_meanadvect_pKE = turbke_offset + 15	! Colatitudinal Mean Advective Flux of turbulent kinetic energy

    ! We have some "known" outputs as well that allow us to verify that
    ! the spherical_io interface is functional
    Integer, Parameter :: dcheck_off = turbke_offset+ 15 ! 216
    Integer, Parameter :: diagnostic1 = dcheck_off+1
    Integer, Parameter :: diagnostic2 = dcheck_off+2
    Integer, Parameter :: test_y11 = dcheck_off+3
    Integer, Parameter :: test_y22 = dcheck_off+4
    Integer, Parameter :: test_y22_sq = dcheck_off+5

    Integer, Parameter :: test_dvrdrdr = dcheck_off+6
    Integer, Parameter :: test_dvrdtdt = dcheck_off+7
    Integer, Parameter :: test_dvrdpdp = dcheck_off+8

    Integer, Parameter :: test_dvtdrdr = dcheck_off+9
    Integer, Parameter :: test_dvpdpdp = dcheck_off+10
    !//////////////////////////////////////////////////////////
    !  Custom Hydo Outputs:  range from 301 through 400
    Integer, Parameter :: custom_hydro_offset = 300
    !Integer, Parameter :: v_grad_s = custom_hydro_offset + 1  ! {Entropy or T} advection


    !//////////////////////////////////////////////////////////
    !    Magnetic Outputs.  
    !    Start at 400 to leave ample room for additional hydro


    !//////////////////////////////////////////////////
    !               Magnetic field components.
    ! Fluctuations (denoted by "p") and azimuthal means
    ! (denoted by "m") may also be output.
    Integer, Parameter :: boffset = 400

    !------------ Field Components ----------!
    Integer, Parameter :: b_r      = boffset+1 ! Full 
    Integer, Parameter :: b_theta  = boffset+2     
    Integer, Parameter :: b_phi    = boffset+3          

    Integer, Parameter :: bp_r     = boffset+4 ! Fluctuating 
    Integer, Parameter :: bp_theta = boffset+5     
    Integer, Parameter :: bp_phi   = boffset+6    


    Integer, Parameter :: bm_r     = boffset+7 ! Mean
    Integer, Parameter :: bm_theta = boffset+8     
    Integer, Parameter :: bm_phi   = boffset+9    

    !------------ Radial Derivatives -------------!
    Integer, Parameter :: db_r_dr      = boffset+10 ! Full 
    Integer, Parameter :: db_theta_dr  = boffset+11  
    Integer, Parameter :: db_phi_dr    = boffset+12   

    Integer, Parameter :: dbp_r_dr     = boffset+13 ! Fluctuating 
    Integer, Parameter :: dbp_theta_dr = boffset+14  
    Integer, Parameter :: dbp_phi_dr   = boffset+15    
  
    Integer, Parameter :: dbm_r_dr     = boffset+16 ! Mean    
    Integer, Parameter :: dbm_theta_dr = boffset+17  
    Integer, Parameter :: dbm_phi_dr   = boffset+18    

    !------------ Theta Derivatives --------------!
    Integer, Parameter :: db_r_dt      = boffset+19 ! Full 
    Integer, Parameter :: db_theta_dt  = boffset+20  
    Integer, Parameter :: db_phi_dt    = boffset+21   

    Integer, Parameter :: dbp_r_dt     = boffset+22 ! Fluctuating 
    Integer, Parameter :: dbp_theta_dt = boffset+23    
    Integer, Parameter :: dbp_phi_dt   = boffset+24    

    Integer, Parameter :: dbm_r_dt     = boffset+25 ! Mean 
    Integer, Parameter :: dbm_theta_dt = boffset+26   
    Integer, Parameter :: dbm_phi_dt   = boffset+27     

    !------------ Phi Derivatives ----------------!
    Integer, Parameter :: db_r_dp      = boffset+28  ! Full
    Integer, Parameter :: db_theta_dp  = boffset+29
    Integer, Parameter :: db_phi_dp    = boffset+30

    Integer, Parameter :: dbp_r_dp     = boffset+31 ! Fluctuating
    Integer, Parameter :: dbp_theta_dp = boffset+32 ! (same as full if    
    Integer, Parameter :: dbp_phi_dp   = boffset+33 !  mean is azimuthal)

    Integer, Parameter :: dbm_r_dp     = boffset+34 ! Mean    
    Integer, Parameter :: dbm_theta_dp = boffset+35 ! (nonzero only when mean is 
    Integer, Parameter :: dbm_phi_dp   = boffset+36 !  spherical, not azimuthal)

    !------------ (1/r) * Theta Derivatives -------!
    Integer, Parameter :: db_r_dtr      = boffset+37 ! Full 
    Integer, Parameter :: db_theta_dtr  = boffset+38 
    Integer, Parameter :: db_phi_dtr    = boffset+39   

    Integer, Parameter :: dbp_r_dtr     = boffset+40 ! Fluctuating
    Integer, Parameter :: dbp_theta_dtr = boffset+41
    Integer, Parameter :: dbp_phi_dtr   = boffset+42

    Integer, Parameter :: dbm_r_dtr     = boffset+43 ! Mean  
    Integer, Parameter :: dbm_theta_dtr = boffset+44
    Integer, Parameter :: dbm_phi_dtr   = boffset+45


    !------(1/{r sintheta})* Phi Derivatives ---!

    Integer, Parameter :: db_r_dprs      = boffset+46 ! Full
    Integer, Parameter :: db_theta_dprs  = boffset+47  
    Integer, Parameter :: db_phi_dprs    = boffset+48

    Integer, Parameter :: dbp_r_dprs     = boffset+49 ! Fluctuating   
    Integer, Parameter :: dbp_theta_dprs = boffset+50
    Integer, Parameter :: dbp_phi_dprs   = boffset+51
 
    Integer, Parameter :: dbm_r_dprs     = boffset+52 ! Mean   
    Integer, Parameter :: dbm_theta_dprs = boffset+53
    Integer, Parameter :: dbm_phi_dprs   = boffset+54  


    !///////////////////////////////////////////////////
    !       Current Density Outputs (Including Ohmic Heating)
    !       This is Curl B -- rename accordingly
    Integer, Parameter :: joffset = boffset+54 ! = 454

    Integer, Parameter :: j_r  = joffset+1      ! Radial Current Density
    Integer, Parameter :: jp_r = joffset+2    
    Integer, Parameter :: jm_r = joffset+3 

    Integer, Parameter :: j_theta  = joffset+4  ! Theta Current Density
    Integer, Parameter :: jp_theta = joffset+5    
    Integer, Parameter :: jm_theta = joffset+6 

    Integer, Parameter :: j_phi  = joffset+7    ! Phi Current Density
    Integer, Parameter :: jp_phi = joffset+8    
    Integer, Parameter :: jm_phi = joffset+9 

    Integer, Parameter :: j_r_sq      = joffset+10 ! (j_r)^2
    Integer, Parameter :: jp_r_sq     = joffset+11 ! (jp_r)^2
    Integer, Parameter :: j_theta_sq  = joffset+12 ! (j_theta)^2
    Integer, Parameter :: jp_theta_sq = joffset+13 ! (jp_theta)^2
    Integer, Parameter :: j_phi_sq    = joffset+14 ! (j_theta)^2
    Integer, Parameter :: jp_phi_sq   = joffset+15 ! (jp_theta)^2
    Integer, Parameter :: j_sq        = joffset+16 ! j dot j
    Integer, Parameter :: jp_sq       = joffset+17 ! j' dot j'
    Integer, Parameter :: ohmic_heat    = joffset+18 ! eta{  j  dot  j}
    Integer, Parameter :: ohmic_heat_pp = joffset+19 ! eta{  j' dot  j'}
    Integer, Parameter :: ohmic_heat_mm = joffset+20 ! eta{ <j> dot <j>}

    !///////////////////////////////////////////////////////////
    !           Magnetic Energies
    Integer, Parameter :: meoffset = joffset+20 ! = 474

    Integer, Parameter :: magnetic_energy = meoffset+1 ! B^2
    Integer, Parameter :: radial_me       = meoffset+2 ! {B_r}^2
    Integer, Parameter :: theta_me        = meoffset+3 ! {B_theta}^2
    Integer, Parameter :: phi_me        = meoffset+4 ! {B_phi}^2

    Integer, Parameter :: mmagnetic_energy = meoffset+5 ! <B>^2
    Integer, Parameter :: radial_mme       = meoffset+6 ! <B_r>^2
    Integer, Parameter :: theta_mme        = meoffset+7 ! <B_theta>^2
    Integer, Parameter :: phi_mme        = meoffset+8 ! <B_phi>^2

    Integer, Parameter :: pmagnetic_energy = meoffset+9  ! {B'}^2
    Integer, Parameter :: radial_pme       = meoffset+10 ! {B_r'}^2
    Integer, Parameter :: theta_pme        = meoffset+11 ! {B_theta'}^2
    Integer, Parameter :: phi_pme        = meoffset+12 ! {B_phi'}^2



    !/////////////////////////// Lorentz Forces ///////////////////////////////
    !  ref%Lorentz_Coeff * (del x B) x B 
    !  ref%Lorentz_Coeff = 1/4pi when dimensional, Pr/(Pr_m E) when nondimesional
    !  j (below) is shorthand for ref%Lorentz_Coeff*delxB  (not quite the current density)
    Integer, Parameter :: loff = meoffset+12 ! =486
    Integer, Parameter :: j_cross_b_r       = loff+1  ! radial component of j x B
    Integer, Parameter :: j_cross_b_theta   = loff+2  !  theta component of j x B
    Integer, Parameter :: j_cross_b_phi     = loff+3  !    phi component of j x B

    Integer, Parameter :: jp_cross_bm_r     = loff+4  ! radial component of j' x <B>  
    Integer, Parameter :: jp_cross_bm_theta = loff+5  !  theta component of j' x <B>
    Integer, Parameter :: jp_cross_bm_phi   = loff+6  !    phi component of j' x <B>

    Integer, Parameter :: jm_cross_bp_r     = loff+7  ! radial component of <j> x B'
    Integer, Parameter :: jm_cross_bp_theta = loff+8  !  theta component of <j> x B'
    Integer, Parameter :: jm_cross_bp_phi   = loff+9  !    phi component of <j> x B'

    Integer, Parameter :: jm_cross_bm_r     = loff+10 ! radial component of <j> x <B>
    Integer, Parameter :: jm_cross_bm_theta = loff+11 !  theta component of <j> x <B>  
    Integer, Parameter :: jm_cross_bm_phi   = loff+12 !    phi component of <j> x <B> 

    Integer, Parameter :: jp_cross_bp_r     = loff+13 ! radial component of j' x B'  
    Integer, Parameter :: jp_cross_bp_theta = loff+14 !  theta component of j' x B'
    Integer, Parameter :: jp_cross_bp_phi   = loff+15 !    phi component of j' x B'


    Integer, Parameter :: maxwell_stress_r     = loff+16 ! -rsintheta {B_r'}{B_phi'}*Lorentz_Coeff
    Integer, Parameter :: maxwell_stress_theta = loff+17 ! -rsintheta {B_theta'}{B_phi'}*Lorentz_Coeff

    Integer, Parameter :: magnetic_torque_r     = loff+18 ! -rsintheta <B_r><B_phi>*Lorentz_Coeff
    Integer, Parameter :: magnetic_torque_theta = loff+19 ! -rsintheta <B_theta><B_phi>*Lorentz_Coeff

    !////////////////////////////// Induction Terms ///////////////////////////
    Integer, Parameter :: indoff = loff + 19 ! = 505

    !--------------- Terms involving v x B (full)
    Integer, Parameter :: induction_shear_r          = indoff+1  ! radial component of {B dot grad v}
    Integer, Parameter :: induction_comp_r           = indoff+2  ! radial component of -{div dot v}B
    Integer, Parameter :: induction_advec_r          = indoff+3  ! radial component of -{v dot grad B}
    Integer, Parameter :: induction_r                = indoff+4  ! radial component of del cros {v x B}
    Integer, Parameter :: induction_diff_r           = indoff+5  ! radial component of -del x (eta {del x B})

    Integer, Parameter :: induction_shear_theta      = indoff+6  ! theta component of {B dot grad v}
    Integer, Parameter :: induction_comp_theta       = indoff+7  ! theta component of -{div dot v}B
    Integer, Parameter :: induction_advec_theta      = indoff+8  ! theta component of -{v dot grad B}
    Integer, Parameter :: induction_theta            = indoff+9  ! theta component of del cros {v x B}
    Integer, Parameter :: induction_diff_theta       = indoff+10 ! theta component of -del x (eta {del x B})

    Integer, Parameter :: induction_shear_phi        = indoff+11 ! phi component of {B dot grad v}
    Integer, Parameter :: induction_comp_phi         = indoff+12 ! phi component of -{div dot v}B
    Integer, Parameter :: induction_advec_phi        = indoff+13 ! phi component of -{v dot grad B}
    Integer, Parameter :: induction_phi              = indoff+14 ! phi component of del cros {v x B}
    Integer, Parameter :: induction_diff_phi         = indoff+15 ! phi component of -del x (eta {del x B})

    !--------------- Terms involving <v> x <B> 
    Integer, Parameter :: induction_shear_vmbm_r     = indoff+16 ! radial component of {<B> dot grad <v>}
    Integer, Parameter :: induction_comp_vmbm_r      = indoff+17 ! radial component of -{div dot <v>}<B>
    Integer, Parameter :: induction_advec_vmbm_r     = indoff+18 ! radial component of -{<v> dot grad <B>}
    Integer, Parameter :: induction_vmbm_r           = indoff+19 ! radial component of del cros {<v> x <B>}
    Integer, Parameter :: induction_diff_bm_r        = indoff+20 ! radial component of -del x (eta {del x <B>})

    Integer, Parameter :: induction_shear_vmbm_theta = indoff+21 ! theta component of {<B> dot grad <v>}
    Integer, Parameter :: induction_comp_vmbm_theta  = indoff+22 ! theta component of -{div dot <v>}<B>
    Integer, Parameter :: induction_advec_vmbm_theta = indoff+23 ! theta component of -{<v> dot grad <B>}
    Integer, Parameter :: induction_vmbm_theta       = indoff+24 ! theta component of del cros {<v> x <B>}
    Integer, Parameter :: induction_diff_bm_theta    = indoff+25 ! theta component of -del x (eta {del x <B>})

    Integer, Parameter :: induction_shear_vmbm_phi   = indoff+26 ! phi component of {<B> dot grad <v>}
    Integer, Parameter :: induction_comp_vmbm_phi    = indoff+27 ! phi component of -{div dot <v>}<B>
    Integer, Parameter :: induction_advec_vmbm_phi   = indoff+28 ! phi component of -{<v> dot grad <B>}
    Integer, Parameter :: induction_vmbm_phi         = indoff+29 ! phi component of del cros {<v> x <B>}
    Integer, Parameter :: induction_diff_bm_phi      = indoff+30 ! phi component of -del x (eta {del x <B>})

    !--------------- Terms involving <v> x B' 
    Integer, Parameter :: induction_shear_vmbp_r     = indoff+31 ! radial component of {B' dot grad <v>}
    Integer, Parameter :: induction_comp_vmbp_r      = indoff+32 ! radial component of -{div dot <v>}B'
    Integer, Parameter :: induction_advec_vmbp_r     = indoff+33 ! radial component of -{<v> dot grad B'}
    Integer, Parameter :: induction_vmbp_r           = indoff+34 ! radial component of del cros {<v> x B'}
    Integer, Parameter :: induction_diff_bp_r        = indoff+35 ! radial component of -del x (eta {del x B'})

    Integer, Parameter :: induction_shear_vmbp_theta = indoff+36 ! theta component of {B' dot grad <v>}
    Integer, Parameter :: induction_comp_vmbp_theta  = indoff+37 ! theta component of -{div dot <v>}B'
    Integer, Parameter :: induction_advec_vmbp_theta = indoff+38 ! theta component of -{<v> dot grad B'}
    Integer, Parameter :: induction_vmbp_theta       = indoff+39 ! theta component of del cros {<v> x B'}
    Integer, Parameter :: induction_diff_bp_theta    = indoff+40 ! theta component of -del x (eta {del x B'})

    Integer, Parameter :: induction_shear_vmbp_phi   = indoff+41 ! phi component of {B' dot grad <v>}
    Integer, Parameter :: induction_comp_vmbp_phi    = indoff+42 ! phi component of -{div dot <v>}B'
    Integer, Parameter :: induction_advec_vmbp_phi   = indoff+43 ! phi component of -{<v> dot grad B'}
    Integer, Parameter :: induction_vmbp_phi         = indoff+44 ! phi component of del cros {<v> x B'}
    Integer, Parameter :: induction_diff_bp_phi      = indoff+45 ! phi component of -del x (eta {del x B'})

    !--------------- Terms involving v' x <B> 
    Integer, Parameter :: induction_shear_vpbm_r     = indoff+46 ! radial component of {<B> dot grad v'}
    Integer, Parameter :: induction_comp_vpbm_r      = indoff+47 ! radial component of -{div dot v'}<B>
    Integer, Parameter :: induction_advec_vpbm_r     = indoff+48 ! radial component of -{v' dot grad <B>}
    Integer, Parameter :: induction_vpbm_r           = indoff+49 ! radial component of del cros {v' x <B>}

    Integer, Parameter :: induction_shear_vpbm_theta = indoff+50 ! theta component of {<B> dot grad v'}
    Integer, Parameter :: induction_comp_vpbm_theta  = indoff+51 ! theta component of -{div dot v'}<B>
    Integer, Parameter :: induction_advec_vpbm_theta = indoff+52 ! theta component of -{v' dot grad <B>}
    Integer, Parameter :: induction_vpbm_theta       = indoff+53 ! theta component of del cros {v' x <B>}

    Integer, Parameter :: induction_shear_vpbm_phi   = indoff+54 ! phi component of {<B> dot grad v'}
    Integer, Parameter :: induction_comp_vpbm_phi    = indoff+55 ! phi component of -{div dot v'}<B>
    Integer, Parameter :: induction_advec_vpbm_phi   = indoff+56 ! phi component of -{v' dot grad <B>}
    Integer, Parameter :: induction_vpbm_phi         = indoff+57 ! phi component of del cros {v' x <B>}

    !--------------- Terms involving v' x B' 
    Integer, Parameter :: induction_shear_vpbp_r     = indoff+58 ! radial component of {B' dot grad v'}
    Integer, Parameter :: induction_comp_vpbp_r      = indoff+59 ! radial component of -{div dot v'}B'
    Integer, Parameter :: induction_advec_vpbp_r     = indoff+60 ! radial component of -{v' dot grad B'}
    Integer, Parameter :: induction_vpbp_r           = indoff+61 ! radial component of del cros {v' x B'}

    Integer, Parameter :: induction_shear_vpbp_theta = indoff+62 ! theta component of {B' dot grad v'}
    Integer, Parameter :: induction_comp_vpbp_theta  = indoff+63 ! theta component of -{div dot v'}B'
    Integer, Parameter :: induction_advec_vpbp_theta = indoff+64 ! theta component of -{v' dot grad B'}
    Integer, Parameter :: induction_vpbp_theta       = indoff+65 ! theta component of del cros {v' x B'}

    Integer, Parameter :: induction_shear_vpbp_phi   = indoff+66 ! phi component of {B' dot grad v'}
    Integer, Parameter :: induction_comp_vpbp_phi    = indoff+67 ! phi component of -{div dot v'}B'
    Integer, Parameter :: induction_advec_vpbp_phi   = indoff+68 ! phi component of -{v' dot grad B'}
    Integer, Parameter :: induction_vpbp_phi         = indoff+69 ! phi component of del cros {v' x B'}


    !///////////////////////////////////////////////////////////////////////////////////////////////
    ! Magnetic Diffusion Terms -- To Be Implemented

    !/////////////////////////////////////////////////////////////////////////////////////////////


    !//////////////////////////////////////////////////////////////////////////////////////////////
    ! User custom magnetic outputs:  Numbers range from 701 to 800

    Integer, Parameter :: custom_mag_offset   = 700
    Integer, Parameter :: cross_helicity      = custom_mag_offset + 1 ! v dot B
    Integer, Parameter :: turb_cross_helicity = custom_mag_offset+2
    !Integer, Parameter :: vB_angle       = ??? ! {v dot B}/{|v||B|} - cosine of angle between v and B



    !///////////////////////////////////
    Real*8, Allocatable :: qty(:,:,:)   ! This variable holds each quantity that we output
    Real*8, Allocatable :: tmp1(:,:,:)  ! A work array
    Real*8, Allocatable :: rweights(:), tweights(:), tmp1d(:)

    !//////////////////////////////////
    ! The ell0 and m0 _ values arrays contain, yes, the ell = 0 and m = 0 values of
    ! everything in buffer at output time.
    Real*8, Allocatable :: ell0_values(:,:), m0_values(:,:,:)

    ! This array will hold fluctuating quantities from the buffer { q - <q>}      
    Real*8, Allocatable :: fbuffer(:,:,:,:)	

    Logical :: azimuthal_mean = .true. ! when false, the m0_values are overwritten with the ell0_values

    !///////////////////////////////////////////////////////////////////////////
    ! A special buffer used for holding second derivatives at output time
    Type(SphericalBuffer) :: d2buffer
    ! ell0 and m0 values of those variables stored in d2buffer
    Real*8, Allocatable :: d2_ell0(:,:), d2_m0(:,:,:)

    ! This array will hold fluctuating quantities from the d2buffer { q - <q>}      
    Real*8, Allocatable :: d2_fbuffer(:,:,:,:)	


    ! Indices within the d2buffer
    Integer :: dvrdrdr, dvrdtdt, dvrdpdp, dvrdrdt, dvrdrdp, dvrdtdp
    Integer :: dvtdrdr, dvtdtdt, dvtdpdp, dvtdrdt, dvtdrdp, dvtdtdp
    Integer :: dvpdrdr, dvpdtdt, dvpdpdp, dvpdrdt, dvpdrdp, dvpdtdp

    Integer :: dtdrdr, dtdtdt, dtdpdp, dtdrdt, dtdrdp, dtdtdp
    Integer :: dpdrdr, dpdtdt, dpdpdp, dpdrdt, dpdrdp, dpdtdp

    Integer :: dbrdrdr, dbrdtdt, dbrdpdp, dbrdrdt, dbrdrdp, dbrdtdp
    Integer :: dbtdrdr, dbtdtdt, dbtdpdp, dbtdrdt, dbtdrdp, dbtdtdp
    Integer :: dbpdrdr, dbpdtdt, dbpdpdp, dbpdrdt, dbpdrdp, dbpdtdp

    Logical :: need_second_derivatives = .false.
Contains

    Subroutine Generate_Diagnostic_Labels()
        ! Define labels for our quantity codes
        Write(6,*)'A line of code.'
        !Call Load_Label(v_r,'V_r')
        !Call Load_Label(v_theta,'V_theta')
        !Call Load_Label(v_phi, 'V_phi')
    End Subroutine Generate_Diagnostic_Labels

    Subroutine Initialize_Diagnostics_Buffer()
        Logical :: dbtrans, dbconfig
        Logical :: test_reduce

        dbtrans = .false.
        dbconfig = .false.
        test_reduce = .false.

        Call cobuffer%init(field_count = cbfcount, config = 'p1a', &
            dynamic_transpose =dbtrans, dynamic_config = dbconfig, &
            hold_cargo = test_reduce, padding = pad_alltoall)        
    End Subroutine Initialize_Diagnostics_Buffer


    Subroutine Compute_Fluctuations(buffer)
        Implicit None
        Integer :: r,k, t, j,jmax
        Real*8, Intent(InOut) :: buffer(1:,my_r%min:,my_theta%min:,1:)
        jmax = size(buffer,4)
        Allocate(fbuffer(1:n_phi,my_r%min:my_r%max,my_theta%min:my_theta%max,1:jmax))
        

        Do j = 1, jmax
            DO_PSI
                fbuffer(PSI,j) = buffer(PSI,j) - m0_values(PSI2,j) 
            END_DO
        Enddo

    End Subroutine Compute_Fluctuations

    Subroutine DeAllocate_Fluctuations()
        Implicit None
        DeAllocate(fbuffer)
    End Subroutine DeAllocate_Fluctuations


End Module Diagnostics_Base
