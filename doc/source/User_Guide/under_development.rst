.. raw:: latex

   \clearpage

.. _under_development:

Under Development
=================

.. _scalar_fields:

Arbitrary Scalar Fields
-----------------------

Rayleigh can solve for additional active, :math:`\chi_{a_i}`, (coupled to the momentum equation through buoyancy) or
passive, :math:`\chi_{p_i}`, scalar fields (where :math:`i` can range up to 50 for each type of scalar).  Both types of field follow a simple advection-diffusion equation:

.. math::
   :label: scalar_evol

   \frac{\partial \chi_{a,p_i}}{\partial t}  + \boldsymbol{v}\cdot\boldsymbol{\nabla}\chi_{a,p_i}  = 0

The number of each type of field can be set using, e.g.:

::

   &physical_controls_namelist
    n_active_scalars = 2
    n_passive_scalars = 2
   /

Other model parameters follow the same convention as temperature but using the prefix `chi_a` or `chi_p` for active and passive
scalars respectively.

See `tests/chi_scalar` for example input files.


Pseudo-Incompressibility
-----------------------

Rayleigh can solve the fluid equations under the pseudo-incompressible approximation. The equation set is as follows:

.. math::
    \begin{aligned}
        \rho_*(r) \left[\frac{\partial\boldsymbol{v}}{\partial t} + \boldsymbol{v \cdot \nabla v}   % Advection
        + 2\Omega_0\hat{\boldsymbol{z}}\times\boldsymbol{v} \right]  =\; % Coriolis
        & \frac{\rho_*(r) g(r)}{c_P} \Theta\, \hat{\boldsymbol{r}} + \frac{\rho_*(r)}{c_P} \frac{d\hat{S}}{dr} P\, \hat{\boldsymbol{r}} % Buoyancy
        - \rho_*(r)\boldsymbol{\nabla}\left(\frac{P}{\hat{\rho}(r)}\right) \\ % Pressure Forces
        &+ \frac{\rho_*(r)}{4\pi\hat{\rho}(r)}\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
       + \frac{\rho_*(r)}{\hat{\rho}(r)}\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\\ % Viscous Forces
        %
        %
       \hat{\rho}(r)\,\hat{T}(r)\left[\frac{\partial \Theta}{\partial t} +\boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta + v_r\frac{d\hat{S}}{dr}\right] =\;
       &\boldsymbol{\nabla}\cdot\left[\hat{\rho}(r)\,\hat{T}(r)\,\kappa(r)\,\boldsymbol{\nabla}\Theta \right] % diffusion
       +Q(r)   % Internal heating
       \\ &+\Phi(r,\theta,\phi)
       +\frac{\eta(r)}{4\pi}\left[\boldsymbol{\nabla}\times\boldsymbol{B}\right]^2\\ % Ohmic Heating
       %
       %
        \boldsymbol{\nabla}\cdot\left[\rho_*(r)\boldsymbol{v}\right] = 0\end{aligned}

where :math:`\rho_*` is the pseudo-density of the background atmsophere. This pseudo-density is a thermodynamic state variable that depends on the mass density and the specific entropy density,

.. math::
   :label: definition_rho*
   
   \rho_* \equiv \hat{\rho} e^{s/c_p}

