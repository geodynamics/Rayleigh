.. raw:: latex

   \clearpage

.. _under_development:

Under Development
=================

The following features are large undertakings that are being incrementally improved.  They are provided here with no guarantee of functionality or accuracy.  After testing and further development they will become supported features and this documentation will be moved to the main section.  Currently this page is primarily intended as a reference for the Rayleigh developers.

.. _scalar_fields:

Arbitrary Scalar Fields
-----------------------

.. include:: arbitrary_scalar_fields.txt 


.. _coupled_bcs:

Coupled Boundary Conditions
---------------------------

.. include:: coupled_bcs.txt 



Pseudo-Incompressibility
-----------------------

Rayleigh can solve the fluid equations under a simple form of the pseudo-incompressible approximation. The equation set is identical to the anelastic equations except for the momentum and continuity equations,

.. math::
    \begin{aligned}
        \hat{\rho}(r) \left[\frac{\partial\boldsymbol{v}}{\partial t} + \boldsymbol{v \cdot \nabla v}   % Advection
        + 2\Omega_0\hat{\boldsymbol{z}}\times\boldsymbol{v} \right]  =\; % Coriolis
        & \frac{\hat{\rho}(r)}{c_P} \left[g(r)\Theta + \frac{d\hat{S}}{dr} \frac{P}{\hat{\rho}(r)}\right] \, \hat{\boldsymbol{r}} % Buoyancy
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



