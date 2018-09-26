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

#define DO_IDX Do t = my_theta%min, my_theta%max;    Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_IDX2 Do t = my_theta%min, my_theta%max;    Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define END_DO enddo; enddo; enddo
#define IDX k,r,t
#define IDX2 r,t


#define DO_IDX3 Do mp = my_mp%min, my_mp%max; m = m_values(mp); Do imi = 1, 2; Do r = my_r%min, my_r%max
#define END_DO enddo; enddo; enddo
#define IDX3 m:l_max,r,imi
Module STABLE_Plugin

  Use ProblemSize
  Use Math_Constants
  Use ReferenceState
  Use TransportCoefficients
  Use General_MPI
  Use Legendre_Transforms, Only : Legendre_Transform
  Use Spectral_Derivatives
  Use Spherical_Buffer
  Use Fourier_Transform

  Implicit None

  !====================================================================
  ! STABLE-specific input parameters

  ! Flag for DR and MC
  ! < 0: No imposed mean flows
  ! = 1: After benchmark of Jouve et al (2007)
  ! = 2: After Dikpati & Gilman models
  ! = 3: After Bangalore models
  Integer, Public, save :: MeanFlows = -1

  ! Flag for turbulent eta
  ! < 0: No change - just use what is defined in Rayleigh
  ! 1: one-step profile as in the Jouve et al (2007) benchmark
  ! 2: two-step profile
  Integer, Public, save :: kindy_eta = -1

  ! Parameters for alpha effect and spotmaker
  ! < 0: Nothing: no poloidal source imposed
  ! = 1: nonlocal axisymmetric alpha-effect
  ! = 2: SpotMaker-D
  ! = 3: SpotMaker-F
  Integer, Public, save :: Poloidal_source = -1

  ! amplitude and profile of the meridional flow
  Real*8, Public, Save :: kindy_mc = 1.d0
  Real*8, Public, Save :: kindy_theta0 = 0.d0
  Real*8, Public, Save :: kindy_rb = 0.69d0

  ! amplitude and profile of the DR
  Real*8, Public, Save :: kindy_omega1 = -1.d0
  Real*8, Public, Save :: kindy_omega2 = -1.d0
  Real*8, Public, Save :: kindy_omega3 = -1.d0
  Real*8, Public, Save :: kindy_omega4 = -1.d0

  ! used for the profiles of various things
  Real*8, Public, Save :: kindy_r1 = -1.d0
  Real*8, Public, Save :: kindy_r2 = -1.d0
  Real*8, Public, Save :: kindy_d1 = -1.d0
  Real*8, Public, Save :: kindy_d2 = -1.d0

  ! used for the turbulent diffusivity eta
  Real*8, Public, Save :: kindy_etaCZ = 1.d12
  Real*8, Public, Save :: kindy_etaSZ = 1.d9
  Real*8, Public, Save :: kindy_eta_r1 = 0.7d0
  Real*8, Public, Save :: kindy_eta_d1 = 0.02d0
  Real*8, Public, Save :: kindy_eta_r2 = 0.95d0
  Real*8, Public, Save :: kindy_eta_d2 = 0.02d0

  ! used for magnetic pumping
  Real*8, Public, Save :: kindy_PumpingCZ = 0.d0
  Real*8, Public, Save :: kindy_PumpingS = 0.d0
  Real*8, Public, Save :: kindy_BotSurPump = 0.d0
  Real*8, Public, Save :: kindy_BotCZPump = 0.d0

  ! input parameters for SpotMaker
  Real*8, Public, Save :: kindy_alpha0 = 1.d0

  ! used for the poloidal source
  Real*8, Public, Allocatable, Save :: kindy_alpha(:)
  Real*8, Public, Allocatable, Save :: kindy_gtheta(:)
  Real*8, Public, Allocatable, Save :: Bmean(:)

  !Spherical buffer object for containing su and its r- and theta- derivatives
  Type(SphericalBuffer) :: grad_su

  Namelist /STABLE_Namelist/ STABLE_flag, MeanFlows, Poloidal_Source, kindy_eta, &
          kindy_mc, kindy_theta0, kindy_rb, kindy_omega1, kindy_omega2, kindy_omega3, &
          kindy_omega4, kindy_r1, kindy_r2, kindy_d1, kindy_d2, kindy_PumpingCZ, &
          kindy_PumpingS, kindy_BotSurPump, kindy_BotCZPump, kindy_alpha0

  !====================================================================
  ! Define public arrays

  ! Imposed DR and MC
  Real*8, allocatable, public, save :: Vmean_r(:,:)
  Real*8, allocatable, public, save :: Vmean_theta(:,:)
  Real*8, allocatable, public, save :: Vmean_phi(:,:)
  Real*8, allocatable, public, save :: dVmean_rdr(:,:)
  Real*8, allocatable, public, save :: dVmean_rdt(:,:)
  Real*8, allocatable, public, save :: dVmean_tdr(:,:)
  Real*8, allocatable, public, save :: dVmean_tdt(:,:)
  Real*8, allocatable, public, save :: dVmean_pdr(:,:)
  Real*8, allocatable, public, save :: dVmean_pdt(:,:)

Contains
  !====================================================================
  ! This routine computes the imposed DR and MC

  Subroutine Initialize_MeanFlows()
     Real*8, allocatable :: rr(:), x(:), chi(:), rstar(:), vpr(:)
     Real*8 :: rb, rc, omega_c, c2, d_tacho, omega_eq, u0, Comega
     Real*8 :: norm, fr1, fr2, fr, gtheta, arg, omega_k, dd, rfr1, xx
     Real*8 :: nu_eq, a2, a4, omega_s, twopi, mu2, mu4, lambda
     Real*8 :: rm2,dm2, wfun1,wfun2, omega_nsg, omc1, omA
     Real*8 :: f1,f2,f1_prime,f2_prime,h1,h2,h1_prime,h2_prime
     Real*8 :: dh1_dtheta,dh2_dtheta,dh_dtheta,hmc,fmc
     Real*8 :: hmc_prime, rth,zMC
     Real*8 :: psi_0, theta_0, r0, gamma_m, beta1, beta2
     Real*8 :: epsilon, midx, vscale, Rnorm, bvr, bvt, Lscale
     Real*8 :: fmc_prime, ee, ee1, ee2, th, pot, tt, rfac,drfac
     Real*8 :: gamCZ, gamS, gamr1,gamr2,gamd1,gamd2,gfr1,gfr2,gx
     Real*8 :: Solar_Radius, InnerRadius
     Integer :: r, t, k, theta

     ! Allocate mean fields
     Allocate(Vmean_r(my_r%min:my_r%max,my_theta%min:my_theta%max))
     Allocate(Vmean_theta(my_r%min:my_r%max,my_theta%min:my_theta%max))
     Allocate(Vmean_phi(my_r%min:my_r%max,my_theta%min:my_theta%max))

     Allocate(dVmean_rdr(my_r%min:my_r%max,my_theta%min:my_theta%max))
     Allocate(dVmean_rdt(my_r%min:my_r%max,my_theta%min:my_theta%max))
     Allocate(dVmean_tdr(my_r%min:my_r%max,my_theta%min:my_theta%max))
     Allocate(dVmean_tdt(my_r%min:my_r%max,my_theta%min:my_theta%max))
     Allocate(dVmean_pdr(my_r%min:my_r%max,my_theta%min:my_theta%max))
     Allocate(dVmean_pdt(my_r%min:my_r%max,my_theta%min:my_theta%max))

     Solar_Radius = 6.96d10
     InnerRadius = Radius(N_r)

     !===============================================================
     If (MeanFlows == 1) Then ! Benchmark setup
        ! Case SC' from Jouve et al (2008)

        ! Note - this is not in general divergenceless when combined
        ! with ref%density

        Allocate(rr(N_R))
        rr = Radius / Solar_Radius

        Comega = 1.4d5
        c2      = 0.2d0
        omega_c = 0.92d0

        rb = 0.65d0
        rc = 0.7d0
        d_tacho = 0.02d0

        omega_eq = Comega * eta(1) / (Solar_Radius**2)
        ! Note that eta(1)=10^11 cm^2/s in the Benchmark, so different eta will produce differnet Omega

        u0 = kindy_mc * 700.d0 * eta(1) / Solar_Radius

        !--------------------------------------------
        Allocate(x(N_R))

        dd = 1.d0 - rb

        !--------------------------------------------
        ! Vr
        norm = - u0 * 2.d0 / Pi

        x = (rr - rb)/(1.d0 - rb)

        DO_IDX2
           arg = Pi * x(r)
           fr = Sin(arg) * x(r)**2 * dd / rr(r)
           gtheta = 3.d0 * costheta(t)**2 - 1.d0
           Vmean_r(IDX2) = norm * fr * gtheta
        END_DO2

        ! Vtheta
        norm = u0 * 2.d0 / Pi
        DO_IDX2
           arg = Pi * x(r)
           fr1 = Sin(arg) * (3.d0*rr(r) - rb)
           fr2 = Cos(arg) * rr(r) * Pi * x(r)
           gtheta = costheta(t) * sintheta(t)
           Vmean_theta(IDX2) = norm * gtheta * x(r) * (fr1 + fr2) / rr(r)
        END_DO2

        ! Vphi
        x = (rr-rc)/d_tacho
        DO_IDX2
           fr = 0.5d0*(1.d0 + erf(x(r)))
           gtheta = 1.d0 - omega_c - c2*CosTheta(t)**2
           omega_k = omega_eq*(omega_c + fr * gtheta) - Angular_Velocity
           Vmean_phi(IDX2) = Radius(r) * Sintheta(t) * Omega_k
        END_DO2

        !--------------------------------------------

        ! Compute derivatives needed in the momentum advection due to these flows.
        ! Vr/dr
        x = (rr - rb)/(1.d0 - rb)
        DO_IDX2
           arg = Pi * x(r)
           xx = (rr(r) + rb) / ((rr(r) - rb) * rr(r))
           fr = Cos(arg) / Sin(arg)
           dVmean_rdr(IDX2) = Vmean_r(IDX2) * ( xx + Pi * fr / (1.d0-rb) ) / Solar_Radius
        END_DO2

        ! dVr/dtheta
        DO_IDX2
           gtheta = costheta(t)*sintheta(t)/(3.d0 * costheta(t)**2 - 1.d0)
           dVmean_rdt(IDX2) = -6.d0 * gtheta * Vmean_r(IDX2)
          !print*, rr(r), 90.d0 - 180.d0/Pi*acos(CosTheta(t)), dVmean_rdt(IDX2)
          ! looks fine
        END_DO2

        ! dVtheta/dr
        norm = u0 * 2.d0 / Pi
        DO_IDX2
           arg = Pi * x(r)
           rfr1 = rr(r) / (rr(r) - rb)
           gtheta = costheta(t) * sintheta(t)
           dVmean_tdr(IDX2) = ( norm * ( (5.d0*rr(r) - 2.d0*rb)/(1.d0 - rb) * Pi * Cos(arg) &
                                         + 3.d0 * Sin(arg) - &
                                         Pi**2 * rr(r) * x(r) * Sin(arg) / (1.d0-rb) &
                                       ) * gtheta * x(r)  / rr(r)  &
                                + Vmean_theta(IDX2) * rb / ( rr(r)*(rr(r)-rb) ) )/Solar_Radius
          ! looks fine
        END_DO2


        ! dVtheta/dt
        norm = u0 * 2.d0 / Pi
        DO_IDX2
           arg = Pi * x(r)
           fr1 = Sin(arg) * (3.d0*rr(r) - rb)
           fr2 = Cos(arg) * rr(r) * Pi * x(r)
           gtheta = costheta(t) * costheta(t) -  sintheta(t) * sintheta(t)
           dVmean_tdt(IDX2) = norm * gtheta * x(r) * (fr1 + fr2) / rr(r)
          ! looks fine
        END_DO2


        ! dVphi/dr
        x = (rr-rc)/d_tacho
        DO_IDX2
           fr = (1.d0/sqrt(pi)) * exp( - x(r)**2 )
           gtheta = SinTheta(t) * (1.d0 - omega_c - c2*CosTheta(t)**2)
           dVmean_pdr(IDX2) = (Vmean_phi(IDX2)/Radius(r)) + &
                              rr(r) * omega_eq * fr * gtheta / d_tacho
          ! looks fine
        END_DO2

        ! dVphi/dt
        DO_IDX2
           fr = 1.d0 + erf(x(r))
           gtheta = CosTheta(t)/SinTheta(t)
           dVmean_pdt(IDX2) = gtheta * ( Vmean_phi(IDX2) + &
                              c2 * omega_eq *Radius(r) * Sintheta(t)**3  * fr)
          ! looks fine
        END_DO2

        ! Clean up
        Deallocate(x,rr)

     !===============================================================
     Else If (MeanFlows == 2) Then ! Dikpati-Gilman setup

       twopi = 2.d0 * Pi

       If (kindy_omega1 < 0) Then
          ! default values
          omega_c = twopi * 432.8d-9
          nu_eq = 460.7d-9
          a2 = -62.9d-9
          a4 = 67.13d-9
          rc = 0.7d0
          d_tacho = 0.05d0

       Else

          omega_c = twopi * kindy_omega1
          nu_eq = kindy_omega2
          a2 = kindy_omega3
          a4 = kindy_omega4
          rc = kindy_r2
          d_tacho = kindy_d2

       EndIf

       !--------------------------------------
       ! Two different normalizations are used for radius

       ! hardwired in
       Lscale = 1.09d10
       Rnorm = Solar_Radius / Lscale

       Allocate(rr(N_R),rstar(N_R))

       rstar = Radius / Solar_Radius
       rr = Radius / Lscale

       !--------------------------------------
       ! radial profile for Omega
       Allocate(x(N_R),chi(N_R))

       x = 2.d0*(rstar-rc)/d_tacho
       chi = 0.5d0*(1.d0+erf(x))

       !--------------------------------------

       do theta = my_theta%min,my_theta%max
          mu2 = Costheta(theta)**2
          mu4 = Costheta(theta)**4
          omega_s = 2.d0*Pi*(nu_eq+a2*mu2+a4*mu4)

          do r = my_r%min, my_r%max

             omega_k = omega_c + chi(r)*(omega_s - omega_c) &
                   & - Angular_Velocity

             lambda = Radius(r) * Sintheta(theta)

             Vmean_phi(IDX2) = lambda * omega_k
         enddo
       enddo

       !--------------------------------------
       ! Now the MC

       If (kindy_rb < 0) Then

          ! default values
          psi_0 = 40.d0
          theta_0 = 0.d0
          rb = 0.65d0 * Rnorm

       Else
          ! defined such that positive kindy_mc gives CCW
          ! cells in the NH
          psi_0 = - kindy_mc

          ! assume theta_0 is input in degrees
          theta_0 = Pi * kindy_theta0 / 180.d0

          rb = kindy_rb * Rnorm

       EndIf

       !--------------------------------------
       ! hardwired in
       vscale = 1.09d2 / 1.1d0
       r0 = (Rnorm - rb)/5.d0
       gamma_m = 3.d0
       beta1 = 0.1d0
       beta2 = 0.3d0
       epsilon = 2.d0 + 1.d-8

       ! first compute rho and put it in chi

       ! kinematic profile
       !midx = 1.5d0
       !chi = ((1.d0/rstar) - 0.97d0)**midx

       chi = ref%density

       ! also define x
       x = (rr - rb)/(Rnorm-rb)

       ! and the normalization
       norm = vscale * psi_0

       pot = Pi / 2.d0

       !--------------------------------------
       ! now loop over indices to define vr and vtheta

       do r = my_r%min, my_r%max

          f1 = sin(Pi*x(r))
          f1_prime = Pi * cos(Pi*x(r)) / (Rnorm - rb)

          ee = - ((rr(r)-r0)/Gamma_m)**2
          f2 = exp(ee)
          f2_prime = -2.d0*(rr(r)-r0)*f2 / (Gamma_m**2)

          fmc = f1 * f2
          fmc_prime = f1_prime*f2 + f1*f2_prime

          do theta = my_theta%min, my_theta%max

             rth = acos(costheta(theta))

             If (rth <= pot) Then
                th = rth
             Else
                th = Pi - rth
             EndIf

             tt = th - theta_0

             ee1 = - beta1 * rr(r) * th**epsilon
             ee2 =   beta2 * rr(r) * (th-pot)

             h1 = 1.d0 - exp(ee1)
             h2 = 1.d0 - exp(ee2)
             hmc = h1*h2

             h1_prime = beta1 * (th**epsilon)*exp(ee1)
             h2_prime = - beta2 * (th-pot)*exp(ee2)
             hmc_prime = h1_prime*h2 + h1*h2_prime

             dh1_dtheta = beta1*rr(r)*epsilon*(th**(epsilon-1.d0)) &
                         & * exp(ee1)
             dh2_dtheta = - beta2*rr(r)*exp(ee2)
             dh_dtheta = dh1_dtheta*h2 + h1*dh2_dtheta

             bvr = hmc + tt*dh_dtheta
             bvt = fmc_prime*hmc + fmc*hmc_prime

             ! flip for southern hemisphere
             if (rth > pot) bvt = - bvt

             lambda = rr(r) * Sin(th)

             Vmean_r(IDX2) = norm*fmc*bvr / (chi(r)*rr(r)*lambda)

             Vmean_theta(IDX2) = - norm*tt*bvt / (chi(r)*lambda)
          enddo
       enddo

      !--------------------------------------
      ! clean up

      Deallocate(x,chi,rr,rstar)

     !===============================================================
     Else If (MeanFlows == 3) Then !Bangalore Meanfield set up

        twopi = 2.d0 * Pi

        If (kindy_omega1 < 0) Then
           ! default values
           omega_c = twopi * 432.8d-9
           nu_eq = 460.7d-9
           a2 = -62.9d-9
           a4 = 67.13d-9
           rc = 0.7d0
           d_tacho = 0.05d0

        Else

           omega_c = twopi * kindy_omega1
           nu_eq = kindy_omega2
           a2 = kindy_omega3
           a4 = kindy_omega4
           rc = kindy_r2
           d_tacho = kindy_d2

        EndIf

        !--------------------------------------
        ! Two different normalizations are used for radius

        ! hardwired in
        Lscale = 1.00d10
        Rnorm = Solar_Radius / Lscale

        Allocate(rr(N_R),rstar(N_R))

        rstar = Radius / Solar_Radius
        rr = Radius / Lscale

        Allocate(x(N_R),chi(N_R))

        ! Dikpati & Charbonneau (199) analytical profile
        !--------------------------------------
        ! radial profile for Omega

        x = 2.d0*(rstar-rc)/d_tacho
        chi = 0.5d0*(1.d0+erf(x))

        !--------------------------------------
        do theta = my_theta%min,my_theta%max
           mu2 = Costheta(theta)**2
           mu4 = Costheta(theta)**4
           omega_s = 2.d0*Pi*(nu_eq+a2*mu2+a4*mu4)

           do r = my_r%min, my_r%max

              omega_k = omega_c + chi(r)*(omega_s - omega_c) &
                        & - Angular_Velocity

              lambda = Radius(r) * Sintheta(theta)

              Vmean_phi(IDX2) = lambda * omega_k

           enddo
         enddo

     !--- Parameters to be used for pumping ---
     !  pumping falls rapidly to zero around 0.7R
         gamCZ = kindy_PumpingCZ
         gamS = kindy_PumpingS
         gamr1 = kindy_BotCZPump * Rnorm
         gamr2 = kindy_BotSurPump * Rnorm
         gamd1 = 0.01d0 * Rnorm
         gamd2 = 0.02d0 * Rnorm

     ! Magnetic pumping coded. At present, we have only radial pumping...
         Allocate(vpr(N_R))
         vpr = gamCZ * 0.5d0*(1.0d0 + erf((rr-gamr1)/gamd1)) &
            & + (gamS-gamCZ)*0.5d0*(1.0d0+erf((rr-gamr2)/gamd2))

     !-- Now the meridional circulation-----------------------------
         If (kindy_rb < 0) Then

            ! default values
            psi_0 = 30.d0
            theta_0 = 0.d0
            rb = 0.67d0 * Rnorm

           Else
            ! defined such that positive kindy_mc gives CCW
            ! cells in the NH
            psi_0 = - kindy_mc

            ! assume theta_0 is input in degrees
            theta_0 = Pi * kindy_theta0 / 180.d0

            rb = kindy_rb * Rnorm

         EndIf

         !--------------------------------------
         ! hardwired in
         vscale = 1.284d0 * 100.0d0
         r0 = 0.45 * Rnorm / 3.5d0
         gamma_m = 3.47d0
         beta1 = 1.0d0
         beta2 = 1.3d0
         epsilon = 2.d0 + 1.d-7

         ! first compute rho and put it in chi

         ! kinematic density profile
         ! midx = 1.5d0
         !chi = ((1.d0/rstar) - 0.95d0)**midx

       chi = ref%density

         ! also define x
         x = (rr - rb)/(Rnorm-rb)

         ! and the normalization
         norm = vscale * psi_0

         pot = Pi / 2.d0

         !--------------------------------------
         ! now loop over indices to define vr and vtheta

         do r = my_r%min, my_r%max

            f1 = sin(Pi*x(r))
            f1_prime = Pi * cos(Pi*x(r)) / (Rnorm - rb)

            ee = - ((rr(r)-r0)/Gamma_m)**2
            f2 = exp(ee)
            f2_prime = -2.d0*(rr(r)-r0)*f2 / (Gamma_m**2)

            fmc = f1 * f2
            fmc_prime = f1_prime*f2 + f1*f2_prime

            ! Makes sure that MC falls smoothly to zero at Rb
            rfac = (rr(r)-rb)
            !drfac = 1.0d0
            !rfac = (rr(r)-rb)**0.3 !If you want less smooth fall at Rb
            !drfac = 0.3d0*(rr(r)-rb)**(-0.7) !If you want less smooth fall at Rb

            do theta = my_theta%min, my_theta%max

               rth = acos(costheta(theta))

               If (rth <= pot) Then
                  th = rth
               Else
                  th = Pi - rth
               EndIf

               tt = th - theta_0

               !ee1 = - beta1 * rr(r) * th**epsilon
               !ee2 =   beta2 * rr(r) * (th-pot)
               ee1 = - beta1 * th**epsilon
               ee2 =   beta2 * (th-pot)

               h1 = 1.d0 - exp(ee1)
               h2 = 1.d0 - exp(ee2)
               hmc = h1*h2

               !h1_prime = beta1 * (th**epsilon)*exp(ee1)
               !h2_prime = - beta2 * (th-pot)*exp(ee2)
               !hmc_prime = h1_prime*h2 + h1*h2_prime

               !dh1_dtheta = beta1*rr(r)*epsilon*(th**(epsilon-1.d0)) &
               !             & * exp(ee1)
               !dh2_dtheta = - beta2*rr(r)*exp(ee2)
               dh1_dtheta = beta1 * epsilon * th**(epsilon-1.d0) * exp(ee1)
               dh2_dtheta = - beta2 * exp(ee2)
               dh_dtheta = dh1_dtheta * h2 + h1 * dh2_dtheta

               bvr = rfac*fmc*dh_dtheta
               bvt = (fmc + rfac*fmc_prime) * hmc

               ! flip for southern hemisphere
               if (rth > pot) bvt = - bvt

               lambda = rr(r) * Sin(th)

               Vmean_r(IDX2) = norm * bvr / (chi(r)*rr(r)*lambda) + vpr(r)
               Vmean_theta(IDX2) = - norm * bvt / (chi(r)*lambda)

           enddo
         enddo

     EndIf

  End Subroutine Initialize_MeanFlows

  !====================================================================
  Subroutine STABLE_eta()
     ! This just re-defines the turbulent magnetic diffusivity, eta to be consistent
     ! with other flux-transport dynamo models

     ! This routine over-writes the existing arrays eta and its logarithmic
     ! derivative, dlneta

     Real*8 :: eta_0, eta_c, rrc, dd, norm, Solar_Radius
     Real*8, allocatable :: x(:), arg(:), rr(:), chi1(:), chi2(:)
     Real*8 :: eta_mid, OuterRadius

     Solar_Radius = 6.96d10

     If (kindy_eta == 1) Then

       eta_0 = 1.d11
       eta_c = 1.d9 / eta_0
       rrc = 0.7d0
       dd = 0.02d0

       norm = 0.5d0*(1.d0 - eta_c)

       Allocate(x(N_R),arg(N_R))
       x = ((Radius/Solar_Radius) - rrc)/dd

       eta = eta_0 * (eta_c + norm * (1.d0 + erf(x)))

       arg = - x**2
       dlneta = eta_0 * norm * 2.d0 * exp(arg) / &
             & (Solar_Radius * dd * sqrt(Pi))
       dlneta = dlneta / eta

     ElseIf (kindy_eta == 2) Then

       eta_mid = kindy_etaCZ
       eta_c = kindy_etaSZ
       OuterRadius = Radius(1)

       Allocate(rr(N_R),x(N_R),chi1(N_R),chi2(N_R),arg(N_R))

       rr = Radius / OuterRadius

       !------------------------------------
       ! two step eta profile

       ! first step
       x = 2.d0 * (rr - kindy_eta_r1) / kindy_eta_d1
       chi1 = 0.5d0 * (1.d0 + erf(x))

       arg = - x**2
       dlneta = eta_mid * &
            2.d0 * exp(arg) / (sqrt(Pi)*kindy_eta_d1*OuterRadius)

       ! second step
       x = 2.d0 * (rr - kindy_eta_r2) / kindy_eta_d2
       chi2 = 0.5d0 * (1.d0 + erf(x))

       arg = - x**2
       dlneta = dlneta + eta_top * &
            2.d0 * exp(arg) / (sqrt(Pi)*kindy_eta_d2*OuterRadius)

       !------------------------------

       eta = eta_c + eta_mid*chi1 + eta_top*chi2

       dlneta = dlneta / eta

     EndIf


  End Subroutine STABLE_eta

  !====================================================================
  Subroutine Init_Poloidal_Source()
     Real*8 :: Solar_Radius, s0, r1, d1
     Real*8, allocatable :: x(:), arg(:)

     If ((Poloidal_Source) < 0) Return

     Allocate(kindy_alpha(N_R),kindy_gtheta(N_Theta))

     Solar_Radius = 6.96d10

     Select Case (Poloidal_Source)

        !-------------------------------------------------------
        Case (1)
           ! Flux transport dynamo benchmark
           ! axisymmetric alpha-effect

       !s0 = kindy_alpha0 * eta_top/Solar_Radius
       s0 = 35.d0 * eta_top/Solar_Radius

           r1 = kindy_r1
           d1 = kindy_d1

           allocate(x(N_R),arg(N_R))

           x = ((Radius/Solar_Radius) - r1)/d1
           arg = ((Radius/Solar_Radius) - 1.d0)/d1
           kindy_alpha = s0 * 0.5d0 * (1.d0+erf(x))*(1.d0-erf(arg))

           kindy_gtheta = SinTheta * CosTheta

           Deallocate(x,arg)

           ! This is the mean Bphi at r=rc, for use in computing
           ! the alpha-effect
           Allocate(Bmean(my_theta%min:my_theta%max))

        !-------------------------------------------------------
        Case (2)
           ! SpotMaker-D ; spot deposition

       ! Not yet implemented

        !-------------------------------------------------------
        Case (3)
           ! SpotMaker-F ; Lifting and twisting flow
       ! Not yet implemented

        !-------------------------------------------------------
        Case Default

       print*,'ERROR: Poloidal Source incorrectly specified'
       stop

        End Select

  End Subroutine Init_Poloidal_Source

  !====================================================================
  Subroutine Compute_SU()
     Real*8 :: Solar_Radius, s0, r1, d1
     Real*8, allocatable :: x(:), arg(:)

     Integer :: SU_B_COUNT(3,2)
    Integer :: k,r,t, mp, m, imi




    ! 1) initialize the buffer
    SU_B_COUNT(:,:) = 3
    Call grad_su%init(field_count = su_b_count, config = 'p3b')  ! Initiliaize in config p3b -- physical space
    Call grad_su%construct('p3b') ! Allocate space for p3b (not done via init)


    ! 2) Load (1/sin(theta)) ds_dtheta  == B/sin(theta)
    !$OMP PARALLEL DO PRIVATE(t,r,k)
    s0 = 2.d0 * Angular_Velocity * pressure_specific_heat
    DO_IDX
        !grad_su%p3b(IDX,2) = s0 * (Radius(r)*cottheta(t)*dVmean_pdr(r,t) - dVmean_pdt(r,t))/ref%gravity(r)  !ds_dtheta
        grad_su%p3b(IDX,2) = s0 * (Radius(r)*cottheta(t)*dVmean_pdr(IDX2) - dVmean_pdt(IDX2))/ref%gravity(r)  !ds_dtheta
        !print*,'88888888',grad_su%p3b(IDX,2)
    END_DO
    !$OMP END PARALLEL DO

    ! 3) FFT
    Call fft_to_spectral(grad_su%p3b, rsc = .true.)

    ! 4) Move to the p2b configuration: theta {r,m} configuration (theta in process)
    ! p3b was deallocated by reform, and config was updated to 'p2b' by reform
    Call grad_su%reform()

    ! 5) Legendre transform
    Call grad_su%construct('s2b')  ! allocate space for spectral configuration
    Call Legendre_Transform(grad_su%p2b,grad_su%s2b)    ! transform


    Call grad_su%deconstruct('p2b') ! free up physical configuration space
    grad_su%config='s2b'  ! update configuration  --- now in ell {r,m}


    ! 6) Take 1/sin(theta) d_by_dtheta ( sin^2 [B/sin(theta)])
    Call d_by_sdtheta(grad_su%s2b,2,3)  ! Take derivative of field index 2 and store in field index 3

    ! 7.)  Compute SU
    DO_IDX3

        grad_su%s2b(mp)%data(IDX3,3) = grad_su%s2b(mp)%data(IDX3,3)/over_l_l_plus1(m:l_max)
        If (m .eq. 0) THEN
            grad_su%s2b(mp)%data(0,r,imi,3) = 0.0d0
        Endif
    END_DO
    ! At this point, SU is contained in grad_su%s2b(:)%data(:,:,:,3)
    ! (1/sintheta) d_by_dtheta{SU} is contained in grad_su%s2b(:)%data(:,:,:,2)

    !8.)  Move to p1b  :  r {ell, m}  -- s2b was deallocated and config updated by reform
    Call grad_su%reform()

    !9.)  Chebyshev transform.  This transform is not in place, so we need a work array. We will use p1a...
    Call grad_su%construct('p1a')

    Call gridcp%To_Spectral(grad_su%p1b,grad_su%p1a)
    grad_su%p1b(:,:,:,:) = grad_su%p1a(:,:,:,:)

    Call gridcp%dealias_buffer(grad_su%p1b)

    ! 10.)  Take 1st radial derivative of thing in index 3,  Store in index 1  (last 1 is derivative order)
    Call gridcp%d_by_dr_cp(3,1,grad_su%p1b,1)

    ! 11.)  Inverse transform and free up p1b space
    Call gridcp%From_Spectral(grad_su%p1b,grad_su%p1a)
    Call grad_su%deconstruct('p1b')


    ! 12.)  Move to r (ell,m) space
    grad_su%config='p1a'
    Call grad_su%reform()

    ! 13. ) Inverse Legendre Transform (we are now in config s2a)
    Call grad_su%construct('p2a')
    Call Legendre_Transform(grad_su%s2a,grad_su%p2a)
    Call grad_su%deconstruct('s2a')
    grad_su%config='p2a'

    ! 14.)  Move back to physical space ('p3a')
    Call grad_su%reform()

    ! 15.) FFT and Fix the sin(theta) so that (1/r)*d_by_dtheta{SU} is stored in index 2
    Call fft_to_spectral(grad_su%p3a, rsc = .true.)
    DO_IDX
        grad_su%p3a(IDX,2) = grad_su%p3a(IDX,2)*sintheta(t)*One_Over_R(r)
    END_DO

    ! NOW
    ! d_by_dr{SU} is stored in grad_su%p3a(:,:,:,1)
    ! (1/r) d_by_dtheta{SU} is stored in grad_su%p3a(:,:,:2)
    ! SU is stored in grad_su%p3a(:,:,:,3)
  End Subroutine Compute_SU



End Module STABLE_Plugin
