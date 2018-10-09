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

    Integer, Parameter :: boffset = teoffset+100  ! :OFFSET CODE:

    !------------ Field Components ----------!
    ! Full
    Integer, Parameter :: b_r      = boffset+1 ! :tex: $B_r$
    Integer, Parameter :: b_theta  = boffset+2 ! :tex: $B_\theta$
    Integer, Parameter :: b_phi    = boffset+3 ! :tex: $B_\phi$

    ! Fluctuating
    Integer, Parameter :: bp_r     = boffset+4 ! :tex: $B_r'$
    Integer, Parameter :: bp_theta = boffset+5 ! :tex: $B_\theta'$
    Integer, Parameter :: bp_phi   = boffset+6  ! :tex: $B_\phi'$

    ! Mean
    Integer, Parameter :: bm_r     = boffset+7  ! :tex: $\overline{B_r}$
    Integer, Parameter :: bm_theta = boffset+8  ! :tex: $\overline{B_\theta}$
    Integer, Parameter :: bm_phi   = boffset+9  ! :tex: $\overline{B_\phi}$

    !------------ Radial Derivatives -------------!
    ! Full
    Integer, Parameter :: db_r_dr      = boffset+10 ! :tex: $\frac{\partial B_r}{\partial r}$
    Integer, Parameter :: db_theta_dr  = boffset+11 ! :tex: $\frac{\partial B_\theta}{\partial r}$
    Integer, Parameter :: db_phi_dr    = boffset+12 ! :tex: $\frac{\partial B_\phi}{\partial r}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_dr     = boffset+13 ! :tex: $\frac{\partial B_r'}{\partial r}$
    Integer, Parameter :: dbp_theta_dr = boffset+14 ! :tex: $\frac{\partial B_\theta'}{\partial r}$
    Integer, Parameter :: dbp_phi_dr   = boffset+15 ! :tex: $\frac{\partial B_\phi'}{\partial r}$
    ! Mean
    Integer, Parameter :: dbm_r_dr     = boffset+16 !  :tex: $\frac{\partial \overline{B_r}}{\partial r}$
    Integer, Parameter :: dbm_theta_dr = boffset+17 !  :tex: $\frac{\partial \overline{B_\theta}}{\partial r}$
    Integer, Parameter :: dbm_phi_dr   = boffset+18 !  :tex: $\frac{\partial \overline{B_\phi}}{\partial r}$

    !------------ Theta Derivatives --------------!
    ! Full
    Integer, Parameter :: db_r_dt      = boffset+19 ! :tex: $\frac{\partial B_r}{\partial \theta}$
    Integer, Parameter :: db_theta_dt  = boffset+20 ! :tex: $\frac{\partial B_\theta}{\partial \theta}$
    Integer, Parameter :: db_phi_dt    = boffset+21 ! :tex: $\frac{\partial B_\phi}{\partial \theta}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_dt     = boffset+22 ! :tex: $\frac{\partial B_r'}{\partial \theta}$
    Integer, Parameter :: dbp_theta_dt = boffset+23 ! :tex: $\frac{\partial B_\theta'}{\partial \theta}$
    Integer, Parameter :: dbp_phi_dt   = boffset+24 ! :tex: $\frac{\partial B_\phi'}{\partial \theta}$
    ! Mean
    Integer, Parameter :: dbm_r_dt     = boffset+25 !  :tex: $\frac{\partial \overline{B_r}}{\partial \theta}$
    Integer, Parameter :: dbm_theta_dt = boffset+26 !  :tex: $\frac{\partial \overline{B_\theta}}{\partial \theta}$
    Integer, Parameter :: dbm_phi_dt   = boffset+27 !  :tex: $\frac{\partial \overline{B_\phi}}{\partial \theta}$

    !------------ Phi Derivatives ----------------!
    ! Full
    Integer, Parameter :: db_r_dp      = boffset+28 ! :tex: $\frac{\partial B_r}{\partial \phi}$
    Integer, Parameter :: db_theta_dp  = boffset+29 ! :tex: $\frac{\partial B_\theta}{\partial \phi}$
    Integer, Parameter :: db_phi_dp    = boffset+30 ! :tex: $\frac{\partial B_\phi}{\partial \phi}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_dp     = boffset+31 ! :tex: $\frac{\partial B_r'}{\partial \phi}$
    Integer, Parameter :: dbp_theta_dp = boffset+32 ! :tex: $\frac{\partial B_\theta'}{\partial \phi}$
    Integer, Parameter :: dbp_phi_dp   = boffset+33 ! :tex: $\frac{\partial B_\phi'}{\partial \phi}$
    ! Mean
    Integer, Parameter :: dbm_r_dp     = boffset+34 ! :tex: $\frac{\partial \overline{B_r}}{\partial \phi}$
    Integer, Parameter :: dbm_theta_dp = boffset+35 ! :tex: $\frac{\partial \overline{B_\theta}}{\partial \phi}$
    Integer, Parameter :: dbm_phi_dp   = boffset+36 ! :tex: $\frac{\partial \overline{B_\phi}}{\partial \phi}$

    !------------ (1/r) * Theta Derivatives -------!
    ! Full
    Integer, Parameter :: db_r_dtr      = boffset+37 ! :tex: $\frac{1}{r}\frac{\partial B_r}{\partial \theta}$
    Integer, Parameter :: db_theta_dtr  = boffset+38 ! :tex: $\frac{1}{r}\frac{\partial B_\theta}{\partial \theta}$
    Integer, Parameter :: db_phi_dtr    = boffset+39 ! :tex: $\frac{1}{r}\frac{\partial B_\phi}{\partial \theta}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_dtr     = boffset+40 ! :tex: $\frac{1}{r}\frac{\partial B_r'}{\partial \theta}$
    Integer, Parameter :: dbp_theta_dtr = boffset+41 ! :tex: $\frac{1}{r}\frac{\partial B_r'}{\partial \theta}$
    Integer, Parameter :: dbp_phi_dtr   = boffset+42 ! :tex: $\frac{1}{r}\frac{\partial B_r'}{\partial \theta}$
    ! Mean
    Integer, Parameter :: dbm_r_dtr     = boffset+43 ! :tex: $\frac{1}{r}\frac{\partial \overline{B_r}}{\partial \theta}$
    Integer, Parameter :: dbm_theta_dtr = boffset+44 ! :tex: $\frac{1}{r}\frac{\partial \overline{B_\theta}}{\partial \theta}$
    Integer, Parameter :: dbm_phi_dtr   = boffset+45 ! :tex: $\frac{1}{r}\frac{\partial \overline{B_\phi}}{\partial \theta}$


    !------(1/{r sintheta})* Phi Derivatives ---!
    ! Full
    Integer, Parameter :: db_r_dprs      = boffset+46 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial B_r}{\partial \phi}$
    Integer, Parameter :: db_theta_dprs  = boffset+47 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial B_\theta}{\partial \phi}$
    Integer, Parameter :: db_phi_dprs    = boffset+48 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial B_\phi}{\partial \phi}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_dprs     = boffset+49 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial B_r'}{\partial \phi}$
    Integer, Parameter :: dbp_theta_dprs = boffset+50 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial B_\theta'}{\partial \phi}$
    Integer, Parameter :: dbp_phi_dprs   = boffset+51 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial B_\phi'}{\partial \phi}$
    ! Mean
    Integer, Parameter :: dbm_r_dprs     = boffset+52 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial \overline{B_r}}{\partial \phi}$
    Integer, Parameter :: dbm_theta_dprs = boffset+53 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial \overline{B_\theta}}{\partial \phi}$
    Integer, Parameter :: dbm_phi_dprs   = boffset+54 ! :tex: $\frac{1}{r\mathrm{sin}\theta}\frac{\partial \overline{B_\phi}}{\partial \phi}$


    !--------- Radial Second Derivatives --------
    ! Full
    Integer, Parameter :: db_r_d2r      = boffset+55 ! :tex: $\frac{\partial^2 B_r}{\partial r^2}$
    Integer, Parameter :: db_theta_d2r  = boffset+56 ! :tex: $\frac{\partial^2 B_\theta}{\partial r^2}$
    Integer, Parameter :: db_phi_d2r    = boffset+57 ! :tex: $\frac{\partial^2 B_\phi}{\partial r^2}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_d2r     = boffset+58 ! :tex: $\frac{\partial^2 B_r'}{\partial r^2}$
    Integer, Parameter :: dbp_theta_d2r = boffset+59 ! :tex: $\frac{\partial^2 B_\theta'}{\partial r^2}$
    Integer, Parameter :: dbp_phi_d2r   = boffset+60 ! :tex: $\frac{\partial^2 B_\phi'}{\partial r^2}$
    ! Mean
    Integer, Parameter :: dbm_r_d2r     = boffset+61 ! :tex: $\frac{\partial^2 \overline{B_r}}{\partial r^2}$
    Integer, Parameter :: dbm_theta_d2r = boffset+62 ! :tex: $\frac{\partial^2 \overline{B_\theta}}{\partial r^2}$
    Integer, Parameter :: dbm_phi_d2r   = boffset+63 ! :tex: $\frac{\partial^2 \overline{B_\phi}}{\partial r^2}$

    !--------- Theta Second Derivatives --------
    ! Full
    Integer, Parameter :: db_r_d2t      = boffset+64 ! :tex: $\frac{\partial^2 B_r}{\partial \theta^2}$
    Integer, Parameter :: db_theta_d2t  = boffset+65 ! :tex: $\frac{\partial^2 B_\theta}{\partial \theta^2}$
    Integer, Parameter :: db_phi_d2t    = boffset+66 ! :tex: $\frac{\partial^2 B_\phi}{\partial \theta^2}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_d2t     = boffset+67 ! :tex: $\frac{\partial^2 B_r'}{\partial \theta^2}$
    Integer, Parameter :: dbp_theta_d2t = boffset+68 ! :tex: $\frac{\partial^2 B_\theta'}{\partial \theta^2}$
    Integer, Parameter :: dbp_phi_d2t   = boffset+69 ! :tex: $\frac{\partial^2 B_\phi'}{\partial \theta^2}$
    ! Mean
    Integer, Parameter :: dbm_r_d2t     = boffset+70 ! :tex: $\frac{\partial^2 \overline{B_r}}{\partial \theta^2}$
    Integer, Parameter :: dbm_theta_d2t = boffset+71 ! :tex: $\frac{\partial^2 \overline{B_\theta}}{\partial \theta^2}$
    Integer, Parameter :: dbm_phi_d2t   = boffset+72 ! :tex: $\frac{\partial^2 \overline{B_\phi}}{\partial \theta^2}$

    !--------- Phi Second Derivatives --------
    ! Full
    Integer, Parameter :: db_r_d2p      = boffset+73 ! :tex: $\frac{\partial^2 B_r}{\partial \phi^2}$
    Integer, Parameter :: db_theta_d2p  = boffset+74 ! :tex: $\frac{\partial^2 B_\theta}{\partial \phi^2}$
    Integer, Parameter :: db_phi_d2p    = boffset+75 ! :tex: $\frac{\partial^2 B_\phi}{\partial \phi^2}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_d2p     = boffset+76 ! :tex: $\frac{\partial^2 B_r'}{\partial \phi^2}$
    Integer, Parameter :: dbp_theta_d2p = boffset+77 ! :tex: $\frac{\partial^2 B_\theta'}{\partial \phi^2}$
    Integer, Parameter :: dbp_phi_d2p   = boffset+78 ! :tex: $\frac{\partial^2 B_\phi'}{\partial \phi^2}$
    ! Mean
    Integer, Parameter :: dbm_r_d2p     = boffset+79 ! :tex: $\frac{\partial^2 \overline{B_r}}{\partial \phi^2}$
    Integer, Parameter :: dbm_theta_d2p = boffset+80 ! :tex: $\frac{\partial^2 \overline{B_\theta}}{\partial \phi^2}$
    Integer, Parameter :: dbm_phi_d2p   = boffset+81 ! :tex: $\frac{\partial^2 \overline{B_\phi}}{\partial \phi^2}$

    !--------- Radial-Theta Second Derivatives --------
    ! Full
    Integer, Parameter :: db_r_d2rt      = boffset+82 ! :tex: $\frac{\partial^2 B_r}{\partial r \partial \theta}$
    Integer, Parameter :: db_theta_d2rt  = boffset+83 ! :tex: $\frac{\partial^2 B_\theta}{\partial r \partial \theta}$
    Integer, Parameter :: db_phi_d2rt    = boffset+84 ! :tex: $\frac{\partial^2 B_\phi}{\partial r \partial    \theta}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_d2rt     = boffset+85 ! :tex: $\frac{\partial^2 B_r'}{\partial r \partial \theta}$
    Integer, Parameter :: dbp_theta_d2rt = boffset+86 ! :tex: $\frac{\partial^2 B_\theta'}{\partial r \partial \theta}$
    Integer, Parameter :: dbp_phi_d2rt   = boffset+87 ! :tex: $\frac{\partial^2 B_\phi'}{\partial r \partial \theta}$
    ! Mean
    Integer, Parameter :: dbm_r_d2rt     = boffset+88 ! :tex: $\frac{\partial^2 \overline{B_r}}{\partial r \partial \theta}$
    Integer, Parameter :: dbm_theta_d2rt = boffset+89 ! :tex: $\frac{\partial^2 \overline{B_\theta}}{\partial r \partial \theta}$
    Integer, Parameter :: dbm_phi_d2rt   = boffset+90 ! :tex: $\frac{\partial^2 \overline{B_\phi}}{\partial r \partial \theta}$

    !--------- Radial-phi Second Derivatives --------
    ! Full
    Integer, Parameter :: db_r_d2rp      = boffset+91 ! :tex: $\frac{\partial^2 B_r}{\partial r \partial \phi}$
    Integer, Parameter :: db_theta_d2rp  = boffset+92 ! :tex: $\frac{\partial^2 B_\theta}{\partial r \partial \phi}$
    Integer, Parameter :: db_phi_d2rp    = boffset+93 ! :tex: $\frac{\partial^2 B_\phi}{\partial r \partial    \phi}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_d2rp     = boffset+94 ! :tex: $\frac{\partial^2 B_r'}{\partial r \partial \phi}$
    Integer, Parameter :: dbp_theta_d2rp = boffset+95 ! :tex: $\frac{\partial^2 B_\theta'}{\partial r \partial \phi}$
    Integer, Parameter :: dbp_phi_d2rp   = boffset+96 ! :tex: $\frac{\partial^2 B_\phi'}{\partial r \partial \phi}$
    ! Mean
    Integer, Parameter :: dbm_r_d2rp     = boffset+97 ! :tex: $\frac{\partial^2 \overline{B_r}}{\partial r \partial \phi}$
    Integer, Parameter :: dbm_theta_d2rp = boffset+98 ! :tex: $\frac{\partial^2 \overline{B_\theta}}{\partial r \partial \phi}$
    Integer, Parameter :: dbm_phi_d2rp   = boffset+99 ! :tex: $\frac{\partial^2 \overline{B_\phi}}{\partial r \partial \phi}$

    !--------- theta-phi Second Derivatives --------
    ! Full
    Integer, Parameter :: db_r_d2tp      = boffset+100 ! :tex: $\frac{\partial^2 B_r}{\partial \theta \partial \phi}$
    Integer, Parameter :: db_theta_d2tp  = boffset+101 ! :tex: $\frac{\partial^2 B_\theta}{\partial \theta \partial \phi}$
    Integer, Parameter :: db_phi_d2tp    = boffset+102 ! :tex: $\frac{\partial^2 B_\phi}{\partial \theta \partial    \phi}$
    ! Fluctuating
    Integer, Parameter :: dbp_r_d2tp     = boffset+103 ! :tex: $\frac{\partial^2 B_r'}{\partial \theta \partial \phi}$
    Integer, Parameter :: dbp_theta_d2tp = boffset+104 ! :tex: $\frac{\partial^2 B_\theta'}{\partial \theta \partial \phi}$
    Integer, Parameter :: dbp_phi_d2tp   = boffset+105 ! :tex: $\frac{\partial^2 B_\phi'}{\partial \theta \partial \phi}$
    ! Mean
    Integer, Parameter :: dbm_r_d2tp     = boffset+106 ! :tex: $\frac{\partial^2 \overline{B_r}}{\partial \theta \partial \phi}$
    Integer, Parameter :: dbm_theta_d2tp = boffset+107 ! :tex: $\frac{\partial^2 \overline{B_\theta}}{\partial \theta \partial \phi}$
    Integer, Parameter :: dbm_phi_d2tp   = boffset+108 ! :tex: $\frac{\partial^2 \overline{B_\phi}}{\partial \theta \partial \phi}$



