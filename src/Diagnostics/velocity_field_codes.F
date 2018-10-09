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

    Integer, Parameter :: voffset = 0  ! :OFFSET CODE:

    !------------ Field Components ----------!
    ! Full
    Integer, Parameter :: v_r      = voffset+1 ! :tex: $v_r$
    Integer, Parameter :: v_theta  = voffset+2 ! :tex: $v_\theta$
    Integer, Parameter :: v_phi    = voffset+3 ! :tex: $v_\phi$

    ! Fluctuating
    Integer, Parameter :: vp_r     = voffset+4 ! :tex: $v_r'$
    Integer, Parameter :: vp_theta = voffset+5 ! :tex: $v_\theta'$
    Integer, Parameter :: vp_phi   = voffset+6  ! :tex: $v_\phi'$

    ! Mean
    Integer, Parameter :: vm_r     = voffset+7  ! :tex: $\overline{v_r}$
    Integer, Parameter :: vm_theta = voffset+8  ! :tex: $\overline{v_\theta}$
    Integer, Parameter :: vm_phi   = voffset+9  ! :tex: $\overline{v_\phi}$

    !------------ Radial Derivatives -------------!
    ! Full
    Integer, Parameter :: dv_r_dr      = voffset+10 ! :tex: $\frac{\partial v_r}{\partial r}$
    Integer, Parameter :: dv_theta_dr  = voffset+11 ! :tex: $\frac{\partial v_\theta}{\partial r}$
    Integer, Parameter :: dv_phi_dr    = voffset+12 ! :tex: $\frac{\partial v_\phi}{\partial r}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_dr     = voffset+13 ! :tex: $\frac{\partial v_r'}{\partial r}$
    Integer, Parameter :: dvp_theta_dr = voffset+14 ! :tex: $\frac{\partial v_\theta'}{\partial r}$
    Integer, Parameter :: dvp_phi_dr   = voffset+15 ! :tex: $\frac{\partial v_\phi'}{\partial r}$
    ! Mean
    Integer, Parameter :: dvm_r_dr     = voffset+16 !  :tex: $\frac{\partial \overline{v_r}}{\partial r}$
    Integer, Parameter :: dvm_theta_dr = voffset+17 !  :tex: $\frac{\partial \overline{v_\theta}}{\partial r}$
    Integer, Parameter :: dvm_phi_dr   = voffset+18 !  :tex: $\frac{\partial \overline{v_\phi}}{\partial r}$

    !------------ Theta Derivatives --------------!
    ! Full
    Integer, Parameter :: dv_r_dt      = voffset+19 ! :tex: $\frac{\partial v_r}{\partial \theta}$
    Integer, Parameter :: dv_theta_dt  = voffset+20 ! :tex: $\frac{\partial v_\theta}{\partial \theta}$
    Integer, Parameter :: dv_phi_dt    = voffset+21 ! :tex: $\frac{\partial v_\phi}{\partial \theta}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_dt     = voffset+22 ! :tex: $\frac{\partial v_r'}{\partial \theta}$
    Integer, Parameter :: dvp_theta_dt = voffset+23 ! :tex: $\frac{\partial v_\theta'}{\partial \theta}$
    Integer, Parameter :: dvp_phi_dt   = voffset+24 ! :tex: $\frac{\partial v_\phi'}{\partial \theta}$
    ! Mean
    Integer, Parameter :: dvm_r_dt     = voffset+25 !  :tex: $\frac{\partial \overline{v_r}}{\partial \theta}$
    Integer, Parameter :: dvm_theta_dt = voffset+26 !  :tex: $\frac{\partial \overline{v_\theta}}{\partial \theta}$
    Integer, Parameter :: dvm_phi_dt   = voffset+27 !  :tex: $\frac{\partial \overline{v_\phi}}{\partial \theta}$

    !------------ Phi Derivatives ----------------!
    ! Full
    Integer, Parameter :: dv_r_dp      = voffset+28 ! :tex: $\frac{\partial v_r}{\partial \phi}$
    Integer, Parameter :: dv_theta_dp  = voffset+29 ! :tex: $\frac{\partial v_\theta}{\partial \phi}$
    Integer, Parameter :: dv_phi_dp    = voffset+30 ! :tex: $\frac{\partial v_\phi}{\partial \phi}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_dp     = voffset+31 ! :tex: $\frac{\partial v_r'}{\partial \phi}$
    Integer, Parameter :: dvp_theta_dp = voffset+32 ! :tex: $\frac{\partial v_\theta'}{\partial \phi}$
    Integer, Parameter :: dvp_phi_dp   = voffset+33 ! :tex: $\frac{\partial v_\phi'}{\partial \phi}$
    ! Mean
    Integer, Parameter :: dvm_r_dp     = voffset+34 ! :tex: $\frac{\partial \overline{v_r}}{\partial \phi}$
    Integer, Parameter :: dvm_theta_dp = voffset+35 ! :tex: $\frac{\partial \overline{v_\theta}}{\partial \phi}$
    Integer, Parameter :: dvm_phi_dp   = voffset+36 ! :tex: $\frac{\partial \overline{v_\phi}}{\partial \phi}$

    !------------ (1/r) * Theta Derivatives -------!
    ! Full
    Integer, Parameter :: dv_r_dtr      = voffset+37 ! :tex: $\frac{1}{r}\frac{\partial v_r}{\partial \theta}$
    Integer, Parameter :: dv_theta_dtr  = voffset+38 ! :tex: $\frac{1}{r}\frac{\partial v_\theta}{\partial \theta}$
    Integer, Parameter :: dv_phi_dtr    = voffset+39 ! :tex: $\frac{1}{r}\frac{\partial v_\phi}{\partial \theta}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_dtr     = voffset+40 ! :tex: $\frac{1}{r}\frac{\partial v_r'}{\partial \theta}$
    Integer, Parameter :: dvp_theta_dtr = voffset+41 ! :tex: $\frac{1}{r}\frac{\partial v_r'}{\partial \theta}$
    Integer, Parameter :: dvp_phi_dtr   = voffset+42 ! :tex: $\frac{1}{r}\frac{\partial v_r'}{\partial \theta}$
    ! Mean
    Integer, Parameter :: dvm_r_dtr     = voffset+43 ! :tex: $\frac{1}{r}\frac{\partial \overline{v_r}}{\partial \theta}$
    Integer, Parameter :: dvm_theta_dtr = voffset+44 ! :tex: $\frac{1}{r}\frac{\partial \overline{v_\theta}}{\partial \theta}$
    Integer, Parameter :: dvm_phi_dtr   = voffset+45 ! :tex: $\frac{1}{r}\frac{\partial \overline{v_\phi}}{\partial \theta}$


    !------(1/{r sintheta})* Phi Derivatives ---!
    ! Full
    Integer, Parameter :: dv_r_dprs      = voffset+46 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial v_r}{\partial \phi}$
    Integer, Parameter :: dv_theta_dprs  = voffset+47 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial v_\theta}{\partial \phi}$
    Integer, Parameter :: dv_phi_dprs    = voffset+48 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial v_\phi}{\partial \phi}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_dprs     = voffset+49 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial v_r'}{\partial \phi}$
    Integer, Parameter :: dvp_theta_dprs = voffset+50 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial v_\theta'}{\partial \phi}$
    Integer, Parameter :: dvp_phi_dprs   = voffset+51 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial v_\phi'}{\partial \phi}$
    ! Mean
    Integer, Parameter :: dvm_r_dprs     = voffset+52 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial \overline{v_r}}{\partial \phi}$
    Integer, Parameter :: dvm_theta_dprs = voffset+53 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial \overline{v_\theta}}{\partial \phi}$
    Integer, Parameter :: dvm_phi_dprs   = voffset+54 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial \overline{v_\phi}}{\partial \phi}$


    !--------- Radial Second Derivatives --------
    ! Full
    Integer, Parameter :: dv_r_d2r      = voffset+55 ! :tex: $\frac{\partial^2 v_r}{\partial r^2}$
    Integer, Parameter :: dv_theta_d2r  = voffset+56 ! :tex: $\frac{\partial^2 v_\theta}{\partial r^2}$
    Integer, Parameter :: dv_phi_d2r    = voffset+57 ! :tex: $\frac{\partial^2 v_\phi}{\partial r^2}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_d2r     = voffset+58 ! :tex: $\frac{\partial^2 v_r'}{\partial r^2}$
    Integer, Parameter :: dvp_theta_d2r = voffset+59 ! :tex: $\frac{\partial^2 v_\theta'}{\partial r^2}$
    Integer, Parameter :: dvp_phi_d2r   = voffset+60 ! :tex: $\frac{\partial^2 v_\phi'}{\partial r^2}$
    ! Mean
    Integer, Parameter :: dvm_r_d2r     = voffset+61 ! :tex: $\frac{\partial^2 \overline{v_r}}{\partial r^2}$
    Integer, Parameter :: dvm_theta_d2r = voffset+62 ! :tex: $\frac{\partial^2 \overline{v_\theta}}{\partial r^2}$
    Integer, Parameter :: dvm_phi_d2r   = voffset+63 ! :tex: $\frac{\partial^2 \overline{v_\phi}}{\partial r^2}$

    !--------- Theta Second Derivatives --------
    ! Full
    Integer, Parameter :: dv_r_d2t      = voffset+64 ! :tex: $\frac{\partial^2 v_r}{\partial \theta^2}$
    Integer, Parameter :: dv_theta_d2t  = voffset+65 ! :tex: $\frac{\partial^2 v_\theta}{\partial \theta^2}$
    Integer, Parameter :: dv_phi_d2t    = voffset+66 ! :tex: $\frac{\partial^2 v_\phi}{\partial \theta^2}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_d2t     = voffset+67 ! :tex: $\frac{\partial^2 v_r'}{\partial \theta^2}$
    Integer, Parameter :: dvp_theta_d2t = voffset+68 ! :tex: $\frac{\partial^2 v_\theta'}{\partial \theta^2}$
    Integer, Parameter :: dvp_phi_d2t   = voffset+69 ! :tex: $\frac{\partial^2 v_\phi'}{\partial \theta^2}$
    ! Mean
    Integer, Parameter :: dvm_r_d2t     = voffset+70 ! :tex: $\frac{\partial^2 \overline{v_r}}{\partial \theta^2}$
    Integer, Parameter :: dvm_theta_d2t = voffset+71 ! :tex: $\frac{\partial^2 \overline{v_\theta}}{\partial \theta^2}$
    Integer, Parameter :: dvm_phi_d2t   = voffset+72 ! :tex: $\frac{\partial^2 \overline{v_\phi}}{\partial \theta^2}$

    !--------- Phi Second Derivatives --------
    ! Full
    Integer, Parameter :: dv_r_d2p      = voffset+73 ! :tex: $\frac{\partial^2 v_r}{\partial \phi^2}$
    Integer, Parameter :: dv_theta_d2p  = voffset+74 ! :tex: $\frac{\partial^2 v_\theta}{\partial \phi^2}$
    Integer, Parameter :: dv_phi_d2p    = voffset+75 ! :tex: $\frac{\partial^2 v_\phi}{\partial \phi^2}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_d2p     = voffset+76 ! :tex: $\frac{\partial^2 v_r'}{\partial \phi^2}$
    Integer, Parameter :: dvp_theta_d2p = voffset+77 ! :tex: $\frac{\partial^2 v_\theta'}{\partial \phi^2}$
    Integer, Parameter :: dvp_phi_d2p   = voffset+78 ! :tex: $\frac{\partial^2 v_\phi'}{\partial \phi^2}$
    ! Mean
    Integer, Parameter :: dvm_r_d2p     = voffset+79 ! :tex: $\frac{\partial^2 \overline{v_r}}{\partial \phi^2}$
    Integer, Parameter :: dvm_theta_d2p = voffset+80 ! :tex: $\frac{\partial^2 \overline{v_\theta}}{\partial \phi^2}$
    Integer, Parameter :: dvm_phi_d2p   = voffset+81 ! :tex: $\frac{\partial^2 \overline{v_\phi}}{\partial \phi^2}$

    !--------- Radial-Theta Second Derivatives --------
    ! Full
    Integer, Parameter :: dv_r_d2rt      = voffset+82 ! :tex: $\frac{\partial^2 v_r}{\partial r \partial \theta}$
    Integer, Parameter :: dv_theta_d2rt  = voffset+83 ! :tex: $\frac{\partial^2 v_\theta}{\partial r \partial \theta}$
    Integer, Parameter :: dv_phi_d2rt    = voffset+84 ! :tex: $\frac{\partial^2 v_\phi}{\partial r \partial    \theta}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_d2rt     = voffset+85 ! :tex: $\frac{\partial^2 v_r'}{\partial r \partial \theta}$
    Integer, Parameter :: dvp_theta_d2rt = voffset+86 ! :tex: $\frac{\partial^2 v_\theta'}{\partial r \partial \theta}$
    Integer, Parameter :: dvp_phi_d2rt   = voffset+87 ! :tex: $\frac{\partial^2 v_\phi'}{\partial r \partial \theta}$
    ! Mean
    Integer, Parameter :: dvm_r_d2rt     = voffset+88 ! :tex: $\frac{\partial^2 \overline{v_r}}{\partial r \partial \theta}$
    Integer, Parameter :: dvm_theta_d2rt = voffset+89 ! :tex: $\frac{\partial^2 \overline{v_\theta}}{\partial r \partial \theta}$
    Integer, Parameter :: dvm_phi_d2rt   = voffset+90 ! :tex: $\frac{\partial^2 \overline{v_\phi}}{\partial r \partial \theta}$

    !--------- Radial-phi Second Derivatives --------
    ! Full
    Integer, Parameter :: dv_r_d2rp      = voffset+91 ! :tex: $\frac{\partial^2 v_r}{\partial r \partial \phi}$
    Integer, Parameter :: dv_theta_d2rp  = voffset+92 ! :tex: $\frac{\partial^2 v_\theta}{\partial r \partial \phi}$
    Integer, Parameter :: dv_phi_d2rp    = voffset+93 ! :tex: $\frac{\partial^2 v_\phi}{\partial r \partial    \phi}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_d2rp     = voffset+94 ! :tex: $\frac{\partial^2 v_r'}{\partial r \partial \phi}$
    Integer, Parameter :: dvp_theta_d2rp = voffset+95 ! :tex: $\frac{\partial^2 v_\theta'}{\partial r \partial \phi}$
    Integer, Parameter :: dvp_phi_d2rp   = voffset+96 ! :tex: $\frac{\partial^2 v_\phi'}{\partial r \partial \phi}$
    ! Mean
    Integer, Parameter :: dvm_r_d2rp     = voffset+97 ! :tex: $\frac{\partial^2 \overline{v_r}}{\partial r \partial \phi}$
    Integer, Parameter :: dvm_theta_d2rp = voffset+98 ! :tex: $\frac{\partial^2 \overline{v_\theta}}{\partial r \partial \phi}$
    Integer, Parameter :: dvm_phi_d2rp   = voffset+99 ! :tex: $\frac{\partial^2 \overline{v_\phi}}{\partial r \partial \phi}$

    !--------- theta-phi Second Derivatives --------
    ! Full
    Integer, Parameter :: dv_r_d2tp      = voffset+100 ! :tex: $\frac{\partial^2 v_r}{\partial \theta \partial \phi}$
    Integer, Parameter :: dv_theta_d2tp  = voffset+101 ! :tex: $\frac{\partial^2 v_\theta}{\partial \theta \partial \phi}$
    Integer, Parameter :: dv_phi_d2tp    = voffset+102 ! :tex: $\frac{\partial^2 v_\phi}{\partial \theta \partial    \phi}$
    ! Fluctuating
    Integer, Parameter :: dvp_r_d2tp     = voffset+103 ! :tex: $\frac{\partial^2 v_r'}{\partial \theta \partial \phi}$
    Integer, Parameter :: dvp_theta_d2tp = voffset+104 ! :tex: $\frac{\partial^2 v_\theta'}{\partial \theta \partial \phi}$
    Integer, Parameter :: dvp_phi_d2tp   = voffset+105 ! :tex: $\frac{\partial^2 v_\phi'}{\partial \theta \partial \phi}$
    ! Mean
    Integer, Parameter :: dvm_r_d2tp     = voffset+106 ! :tex: $\frac{\partial^2 \overline{v_r}}{\partial \theta \partial \phi}$
    Integer, Parameter :: dvm_theta_d2tp = voffset+107 ! :tex: $\frac{\partial^2 \overline{v_\theta}}{\partial \theta \partial \phi}$
    Integer, Parameter :: dvm_phi_d2tp   = voffset+108 ! :tex: $\frac{\partial^2 \overline{v_\phi}}{\partial \theta \partial \phi}$



