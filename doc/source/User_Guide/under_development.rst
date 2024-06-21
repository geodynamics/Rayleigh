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

Rayleigh can solve the fluid equations under a simple form of the pseudo-incompressible approximation. The equation set is identical to the anelastic equations except for the momentum and continuity equations,

.. math::
    \begin{aligned}
        \hat{\rho}(r) \left[\frac{\partial\boldsymbol{v}}{\partial t} + \boldsymbol{v \cdot \nabla v}   % Advection
        + 2\Omega_0\hat{\boldsymbol{z}}\times\boldsymbol{v} \right]  =\; % Coriolis
        & \frac{\hat{\rho}(r)}{c_P} \left(g(r)\Theta + \frac{d\hat{S}}{dr} \frac{P}{\hat{\rho}(r)}\right) \, \hat{\boldsymbol{r}} % Buoyancy
        - \hat{\rho}(r)\boldsymbol{\nabla}\left(\frac{P}{\hat{\rho}(r)}\right) % Pressure Forces
    \\ 
        &+ \frac{1}{4\pi}\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
       + \boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}  % Viscous Forces
    \\
        \boldsymbol{\nabla}\cdot\left[\hat{\rho}_*(r)\boldsymbol{v}\right] =\; &0.  % Continuity
    \end{aligned}

The momentum equation solved in anelastic mode is augmented by a buoyancy term that is normally ignored under the LBR (Lantz-Braginsky-Roberts) formulation of the anelastic approximation.  This new term is proportional to the background entropy gradient.  Hence, in a convection zone, this term is quite small and its exclusion is well justified.  However, in a stable layer, this term is not small and should be considered.

The continuity equation is still a solenoidal constraint, but instead of the mass flux being divergenceless, a "pseudo-mass" flux is divergenceless.  In the solenoidal constraint above, the quantity :math:`\hat{\rho}_*` is the pseudo-density of the background atmosphere. This pseudo-density is a thermodynamic state variable that depends on the mass density and the specific entropy density,

.. math::
   :label: definition_rho*
   
   \hat{\rho}_*(r) \equiv \hat{\rho}(r) \, e^{\hat{S}(r)/c_P}

