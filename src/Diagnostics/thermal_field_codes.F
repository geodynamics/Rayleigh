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

    !//////////////////////////////////////////////////////////////////
    !  Pressure, Entropy or Temperature, and Their Derivatives
    !  Note:  In Boussinesq Mode, Temperature is Output Instead of Entropy
    Integer, Parameter :: pt_off = keoffset+100 ! :OFFSET CODE:

    !------------ Fields ---------------------!
    ! Full
    Integer, Parameter :: entropy    = pt_off+1 !  :tex:  $\Theta$
    Integer, Parameter :: pressure   = pt_off+2 !  :tex:  $P$
    ! Fluctuating
    Integer, Parameter :: entropy_p  = pt_off+3 !  :tex:  $\Theta'$
    Integer, Parameter :: pressure_p = pt_off+4 !  :tex:  $P'$
    ! Mean
    Integer, Parameter :: entropy_m  = pt_off+5 !  :tex:  $\overline{\Theta}$
    Integer, Parameter :: pressure_m = pt_off+6 !  :tex:  $\overline{P}$

    !------------ Radial Derivatives --------------!
    !  Full
    Integer, Parameter :: entropy_dr     = pt_off+7  ! :tex:  $\frac{\partial \Theta}{\partial r}$
    Integer, Parameter :: pressure_dr    = pt_off+8  ! :tex:  $\frac{\partial P}{\partial r}$
    ! Fluctuating
    Integer, Parameter :: entropy_p_dr   = pt_off+9  ! :tex:  $\frac{\partial \Theta'}{\partial r}$
    Integer, Parameter :: pressure_p_dr  = pt_off+10 ! :tex:  $\frac{\partial P'}{\partial r}$
    ! Mean
    Integer, Parameter :: entropy_m_dr   = pt_off+11 ! :tex:  $\frac{\partial \overline{\Theta}}{\partial r}$
    Integer, Parameter :: pressure_m_dr = pt_off+12  ! :tex:  $\frac{\partial \overline{P}}{\partial r}$

    !------------ Theta Derivatives ---------------!
    Integer, Parameter :: entropy_dtheta     = pt_off+13  ! :tex:  $\frac{\partial \Theta}{\partial \theta}$
    Integer, Parameter :: pressure_dtheta    = pt_off+14  ! :tex:  $\frac{\partial P}{\partial \theta}$

    Integer, Parameter :: entropy_p_dtheta   = pt_off+15  ! :tex:  $\frac{\partial \Theta'}{\partial \theta}$
    Integer, Parameter :: pressure_p_dtheta  = pt_off+16  ! :tex:  $\frac{\partial P'}{\partial \theta}$

    Integer, Parameter :: entropy_m_dtheta   = pt_off+17  ! :tex:  $\frac{\partial \overline{\Theta}}{\partial \theta}$
    Integer, Parameter :: pressure_m_dtheta  = pt_off+18  ! :tex:  $\frac{\partial \overline{P}}{\partial \theta}$

    !------------ Phi Derivatives -----------------!
    Integer, Parameter :: entropy_dphi     = pt_off+19  ! :tex:  $\frac{\partial \Theta}{\partial \phi}$
    Integer, Parameter :: pressure_dphi    = pt_off+20  ! :tex:  $\frac{\partial P}{\partial \phi}$

    Integer, Parameter :: entropy_p_dphi   = pt_off+21  ! :tex:  $\frac{\partial \Theta'}{\partial \phi}$
    Integer, Parameter :: pressure_p_dphi  = pt_off+22  ! :tex:  $\frac{\partial P'}{\partial \phi}$

    Integer, Parameter :: entropy_m_dphi   = pt_off+23  ! :tex:  $\frac{\partial \overline{\Theta}}{\partial \phi}$
    Integer, Parameter :: pressure_m_dphi  = pt_off+24  ! :tex:  $\frac{\partial \overline{P}}{\partial \phi}$

    !------------ (1/r) * Theta Derivatives --------!
    Integer, Parameter :: entropy_dtr     = pt_off+25 !  :tex:  $\frac{1}{r}\frac{\partial \Theta}{\partial \theta}$
    Integer, Parameter :: pressure_dtr    = pt_off+26 !  :tex:  $\frac{1}{r}\frac{\partial P}{\partial \theta}$

    Integer, Parameter :: entropy_p_dtr   = pt_off+27 !  :tex:  $\frac{1}{r}\frac{\partial \Theta'}{\partial \theta}$
    Integer, Parameter :: pressure_p_dtr  = pt_off+28 !  :tex:  $\frac{1}{r}\frac{\partial P'}{\partial \theta}$

    Integer, Parameter :: entropy_m_dtr   = pt_off+29 !  :tex:  $\frac{1}{r}\frac{\partial \overline{\Theta}}{\partial \theta}$
    Integer, Parameter :: pressure_m_dtr  = pt_off+30 !  :tex:  $\frac{1}{r}\frac{\partial \overline{P}}{\partial \theta}$

    !--- (1/{r sintheta}) * Phi Derivatives ---------!
    Integer, Parameter :: entropy_dprs     = pt_off+31 ! :tex: $\frac{1}{r\mathrm{sin}\theta} \frac{\partial \Theta}{\partial \phi}$
    Integer, Parameter :: pressure_dprs    = pt_off+32 ! :tex: $\frac{1}{r\mathrm{sin}\theta} \frac{\partial P}{\partial \phi}$

    Integer, Parameter :: entropy_p_dprs   = pt_off+33 ! :tex: $\frac{1}{r\mathrm{sin}\theta} \frac{\partial \Theta'}{\partial \phi}$
    Integer, Parameter :: pressure_p_dprs  = pt_off+34 ! :tex: $\frac{1}{r\mathrm{sin}\theta} \frac{\partial P'}{\partial \phi}$

    Integer, Parameter :: entropy_m_dprs   = pt_off+35 ! :tex: $\frac{1}{r\mathrm{sin}\theta} \frac{\partial \overline{\Theta}}{\partial \phi}$
    Integer, Parameter :: pressure_m_dprs  = pt_off+36 ! :tex: $\frac{1}{r\mathrm{sin}\theta} \frac{\partial \overline{P}}{\partial \phi}$

    !------------ Radial Second Derivatives --------------!
    Integer, Parameter :: entropy_d2r     = pt_off+37  ! :tex:  $\frac{\partial^2 \Theta}{\partial r^2}$
    Integer, Parameter :: pressure_d2r    = pt_off+38  ! :tex:  $\frac{\partial^2 P}{\partial r^2}$

    Integer, Parameter :: entropy_p_d2r   = pt_off+39  ! :tex:  $\frac{\partial^2 \Theta'}{\partial r^2}$
    Integer, Parameter :: pressure_p_d2r  = pt_off+40  ! :tex:  $\frac{\partial^2 P'}{\partial r^2}$

    Integer, Parameter :: entropy_m_d2r   = pt_off+41  !  :tex:  $\frac{\partial^2 \overline{\Theta}}{\partial r^2}$
    Integer, Parameter :: pressure_m_d2r  = pt_off+42  ! :tex:  $\frac{\partial^2 \overline{P}}{\partial r^2}$

    !------------ Theta Second Derivatives --------------!
    Integer, Parameter :: entropy_d2t     = pt_off+43  ! :tex:  $\frac{\partial^2 \Theta}{\partial \theta^2}$
    Integer, Parameter :: pressure_d2t    = pt_off+44  ! :tex:  $\frac{\partial^2 P}{\partial \theta^2}$

    Integer, Parameter :: entropy_p_d2t   = pt_off+45  ! :tex:  $\frac{\partial^2 \Theta'}{\partial \theta^2}$
    Integer, Parameter :: pressure_p_d2t  = pt_off+46  ! :tex:  $\frac{\partial^2 P'}{\partial \theta^2}$

    Integer, Parameter :: entropy_m_d2t   = pt_off+47  ! :tex:  $\frac{\partial^2 \overline{\Theta}}{\partial \theta^2}$
    Integer, Parameter :: pressure_m_d2t  = pt_off+48  ! :tex:  $\frac{\partial^2 \overline{P}}{\partial \theta^2}$

    !------------ Phi Second Derivatives --------------!
    Integer, Parameter :: entropy_d2p     = pt_off+49 ! :tex:  $\frac{\partial^2 \Theta}{\partial \phi^2}$
    Integer, Parameter :: pressure_d2p    = pt_off+50 ! :tex:  $\frac{\partial^2 P}{\partial \phi^2}$

    Integer, Parameter :: entropy_p_d2p   = pt_off+51  ! :tex:  $\frac{\partial^2 \Theta'}{\partial \phi^2}$
    Integer, Parameter :: pressure_p_d2p  = pt_off+52  ! :tex:  $\frac{\partial^2 P'}{\partial \phi^2}$

    Integer, Parameter :: entropy_m_d2p   = pt_off+53  ! :tex:  $\frac{\partial^2 \overline{\Theta}}{\partial \phi^2}$
    Integer, Parameter :: pressure_m_d2p  = pt_off+54  ! :tex:  $\frac{\partial^2 \overline{P}}{\partial \phi^2}$

    !------------ Radial-Theta Second Derivatives --------------!
    Integer, Parameter :: entropy_d2rt     = pt_off+55  ! :tex:  $\frac{\partial^2 \Theta}{\partial r \partial \theta}$
    Integer, Parameter :: pressure_d2rt    = pt_off+56  ! :tex:  $\frac{\partial^2 P}{\partial r \partial \theta}$

    Integer, Parameter :: entropy_p_d2rt   = pt_off+57  ! :tex:  $\frac{\partial^2 \Theta'}{\partial r \partial \theta}$
    Integer, Parameter :: pressure_p_d2rt  = pt_off+58  ! :tex:  $\frac{\partial^2 P'}{\partial r \partial \theta}$

    Integer, Parameter :: entropy_m_d2rt   = pt_off+59  ! :tex:  $\frac{\partial^2 \overline{\Theta}}{\partial r \partial \theta}$
    Integer, Parameter :: pressure_m_d2rt  = pt_off+60  ! :tex:  $\frac{\partial^2 \overline{P}}{\partial r \partial \theta}$

    !------------ Radial-Phi Derivatives --------------!
    Integer, Parameter :: entropy_d2rp     = pt_off+61  ! :tex:  $\frac{\partial^2 \Theta}{\partial r \partial \phi}$
    Integer, Parameter :: pressure_d2rp    = pt_off+62  ! :tex:  $\frac{\partial^2 P}{\partial r \partial \phi}$

    Integer, Parameter :: entropy_p_d2rp   = pt_off+63  ! :tex:  $\frac{\partial^2 \Theta'}{\partial r \partial \phi}$
    Integer, Parameter :: pressure_p_d2rp  = pt_off+64  ! :tex:  $\frac{\partial^2 P'}{\partial r \partial \phi}$

    Integer, Parameter :: entropy_m_d2rp   = pt_off+65  ! :tex:  $\frac{\partial^2 \overline{\Theta}}{\partial r \partial \phi}$
    Integer, Parameter :: pressure_m_d2rp  = pt_off+66  ! :tex:  $\frac{\partial^2 \overline{P}}{\partial r \partial \phi}$

    !------------ Theta-Phi Derivatives --------------!
    Integer, Parameter :: entropy_d2tp     = pt_off+67  ! :tex:  $\frac{\partial^2 \Theta}{\partial \theta \partial \phi}$
    Integer, Parameter :: pressure_d2tp    = pt_off+68  ! :tex:  $\frac{\partial^2 P}{\partial \theta \partial \phi}$

    Integer, Parameter :: entropy_p_d2tp   = pt_off+69  ! :tex:  $\frac{\partial^2 \Theta'}{\partial \theta \partial \phi}$
    Integer, Parameter :: pressure_p_d2tp  = pt_off+70  ! :tex:  $\frac{\partial^2 P'}{\partial \theta \partial \phi}$

    Integer, Parameter :: entropy_m_d2tp   = pt_off+71  ! :tex:  $\frac{\partial^2 \overline{\Theta}}{\partial \theta \partial \phi}$
    Integer, Parameter :: pressure_m_d2tp  = pt_off+72  ! :tex:  $\frac{\partial^2 \overline{P}}{\partial \theta \partial \phi}$

    !--------- Also compute terms of the form rho_bar d_by_dr(P/rho_bar)
    Integer, Parameter :: rhopressure_dr  = pt_off+73 ! :tex:  $\frac{\partial}{\partial r} \left( \frac{P}{\hat{\rho}}\right)$
    Integer, Parameter :: rhopressurep_dr = pt_off+74 ! :tex:  $\frac{\partial}{\partial r} \left( \frac{P'}{\hat{\rho}}\right)$
    Integer, Parameter :: rhopressurem_dr = pt_off+75 ! :tex:  $\frac{\partial}{\partial r} \left( \frac{\overline{P}}{\hat{\rho}}\right)$
