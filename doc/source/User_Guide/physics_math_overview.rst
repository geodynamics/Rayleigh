.. raw:: latex

   \clearpage

.. _physics_math:

Underlying Physics 
=======================================

Rayleigh solves the MHD equations in spherical geometry under the
Boussinesq and anelastic approximations. This section will provide a 
basic overview of those equations as well as the mathematical approach 
Rayleigh uses to solve them. 

.. _notation:

Notation and Conventions
---------------------------

Vector and Tensor Notation
^^^^^^^^^^^^^^^^^^^^^^^^^^

All vector quantities are represented in bold italics. Components of a
vector are indicated in non-bold italics, along with a subscript
indicating the direction associated with that component. Unit vectors
are written in lower-case, bold math font and are indicated by the use
of a *hat* character. For example, a vector quantity
:math:`\boldsymbol{a}` would represented as

.. math::
   :label: vcomp

       \boldsymbol{a} = a_r\boldsymbol{\hat{r}}+a_\theta\boldsymbol{\hat{\theta}}+a_\phi\boldsymbol{\hat{\phi}}.

The symbols (:math:`\boldsymbol{\hat{r}}`,
:math:`\boldsymbol{\hat{\theta}}`, :math:`\boldsymbol{\hat{\phi}}`)
indicate the unit vectors in the
(:math:`r`,\ :math:`\theta`,\ :math:`\phi`) directions, and
(:math:`a_r`, :math:`a_\theta`, :math:`a_\phi`) indicate the components
of :math:`\boldsymbol{a}` along those directions respectively.

Vectors may be written in lower case, as with the velocity field
:math:`\boldsymbol{v}`, or in upper case as with the magnetic field
:math:`\boldsymbol{B}`. Tensors are indicated by bold, upper-case,
script font, as with the viscous stress tensor
:math:`\boldsymbol{\mathcal{D}}`. Tensor components are indicated in
non-bold, and with directional subscripts (i.e.,
:math:`\mathcal{D}_{r\theta}`).

Reference-State Values
^^^^^^^^^^^^^^^^^^^^^^

The *hat* notation is also used to indicate reference-state quantities.
These quantities are scalar, and they are not written in bold font. They
vary only in radius and have no :math:`\theta`-dependence or
:math:`\phi`-dependence. The reference-state density is indicated by
:math:`\hat{\rho}` and the reference-state temperature by
:math:`\hat{T}`, for instance.

Averaged and Fluctuating Values
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most of the output variables have been decomposed into a
zonally-averaged value, and a fluctuation about that average. The
average is indicated by an overbar, such that

.. math::
   :label: avging

       \overline{a}\equiv \frac{1}{2\pi}\int_{0}^{2\pi} a(r,\theta,\phi)\, \mathrm{d}\phi.

Fluctations about that average are indicated by a *prime* superscript,
such that

.. math::
   :label: prime

       a'(r,\theta,\phi)\equiv a(r,\theta,\phi)-\overline{a}(r,\theta)

Finally, some quantities are averaged over the full sphere. These are
indicated by a double-zero subscript (i.e. :math:`\ell=0,\,m=0`), such
that

.. math::
   :label: fullsph

   a_{00}\equiv \frac{1}{4\pi}\int_{0}^{2\pi}\int_{0}^{\pi} a(r,\theta,\phi)\, r\mathrm{sin}\,\theta\mathrm{d}\theta\mathrm{d}\phi.

.. _equations_solved:

The System of Equations Solved in Rayleigh
------------------------------------

Rayleigh solves the Boussinesq or anelastic MHD equations in spherical
geometry. Both the equations that Rayleigh solves and its diagnostics
can be formulated either dimensionally or nondimensionally. A
nondimensional Boussinesq formulation, as well as dimensional and
nondimensional anelastic formulations (based on a polytropic reference
state) are provided as part of Rayleigh. The user may employ alternative
formulations via the custom Reference-state interface. To do so, they
must specify the functions :math:`\mathrm{f}_i` and the constants
:math:`c_i` in Equations :eq:`momentum`-:eq:`induction` at
input time (*in development*).

The general form of the momentum equation solved by Rayleigh is given by

.. math::
   :label: momentum

       \mathrm{f}_1(r)\left[\frac{\partial \boldsymbol{v}}{\partial t}  + \boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}  %advection
        + c_1\boldsymbol{\hat{z}}\times\boldsymbol{v} \right]  =\ % Coriolis
       &c_2\,\mathrm{f}_2(r)\Theta\,\boldsymbol{\hat{r}} % buoyancy
        - c_3\,\mathrm{f}_1(r)\boldsymbol{\nabla}\left(\frac{P}{\mathrm{f}_1(r)}\right) % pressure
        \\
        &+ c_4\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
        + c_5\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}},

where the stress tensor :math:`\mathcal{D}` is given by

.. math::
   :label: stress_tensor

       \mathcal{D}_{ij} = 2\mathrm{f}_1(r)\,\mathrm{f}_3(r)\left[e_{ij} - \frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)\delta_{ij}\right].
Here :math:`e_{ij}` and :math:`\delta_{ij}` refer to the standard rate-of-strain tensor and Kronecker delta, respectively. 

The velocity field is denoted by :math:`\boldsymbol{v}`, the thermal
anomoly by :math:`\Theta`, the pressure by :math:`P`, and the magnetic
field by :math:`\boldsymbol{B}`. All four of these quantities (eight, if you count the three components each for :math:`\boldsymbol{v}` and :math:`\boldsymbol{B}`) 
are 3-dimensional functions of position, in contrast to the 1-dimensional
functions of radius :math:`\mathrm{f}_i(r)`. The velocity and magnetic
fields are subject to the constraints

.. math::
   :label: v_constrain

       \boldsymbol{\nabla}\cdot\left[\mathrm{f}_1(r)\,\boldsymbol{v}\right] = 0

and

.. math::
   :label: divB

       \boldsymbol{\nabla}\cdot\boldsymbol{B}=0,

respectively. The evolution of :math:`\Theta` is described by

.. math::
   :label: theta_evol

   \mathrm{f}_1(r)\,\mathrm{f}_4(r)\left[\frac{\partial \Theta}{\partial t}  + \boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta + c_{11}\,\mathrm{f}_{14}(r)v_r\right] =\
       c_6\,\boldsymbol{\nabla}\cdot\left[\mathrm{f}_1(r)\,\mathrm{f}_4(r)\,\mathrm{f}_5(r)\,\boldsymbol{\nabla}\Theta \right] \\
        +\ c_{10}\,\mathrm{f}_6(r)
        + c_8\,\Phi(r,\theta,\phi)
        + c_9\,\mathrm{f}_7(r)|\boldsymbol{\nabla}\times\boldsymbol{B}|^2,

where the viscous heating :math:`\Phi` is given by

.. math::
   :label: vischeat

       \Phi(r,\theta,\phi) = c_5\,\mathcal{D}_{ij}e_{ij} &= 2\,c_5\,\mathrm{f}_1(r)\mathrm{f}_3(r)\left[e_{ij}e_{ij} - \frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)^2\right] \\
       &= 2\,c_5\,\mathrm{f}_1(r)\mathrm{f}_3(r)\left[e_{ij} - \frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)\delta_{ij}\right]^2.

Finally, the evolution of :math:`\boldsymbol{B}` is described by the
induction equation

.. math::
   :label: induction

       \frac{\partial \boldsymbol{B}}{\partial t} = \boldsymbol{\nabla}\times\left[\boldsymbol{v}\times\boldsymbol{B} - c_7\,\mathrm{f}_7(r)\boldsymbol{\nabla}\times\boldsymbol{B}\,\right].

Note that when Rayleigh actually solves the equations, the following additional derivative functions are used:

.. math::
    \mathrm{f}_8(r) &= \frac{d\ln{\mathrm{f}_1}}{dr}\\
    \mathrm{f}_9(r) &= \frac{d^2\ln{\mathrm{f}_1}}{dr^2}\\
    \mathrm{f}_{10}(r) &= \frac{d\ln{\mathrm{f}_4}}{dr}\\
    \mathrm{f}_{11}(r) &= \frac{d\ln{\mathrm{f}_3}}{dr}\\
    \mathrm{f}_{12}(r) &= \frac{d\ln{\mathrm{f}_5}}{dr}\\
    \mathrm{f}_{13}(r) &= \frac{d\ln{\mathrm{f}_7}}{dr}.

When supplying a custom reference state, the user may specify the six derivative functions "by hand." If the user fails to do so, Rayleigh will compute the required derivatives (only if the user supplies the function whose derivative is to be taken) from the function's Chebyshev coefficients. 

Note that equations :eq:`momentum`-:eq:`induction` could have been formulated in other ways. For instance, we could combine
:math:`\mathrm{f}_1` and :math:`\mathrm{f}_3` into a single function in
Equation :eq:`vischeat`. The form of the equations
presented here has been chosen to reflect that actually used in the
code, which was originally written dimensionally. 

We now describe the nondimensional Boussinesq, and dimensional/nondimensional anelastic formulations used in the code.

.. _boussinesq:

Nondimensional Boussinesq Formulation of the MHD Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rayleigh can be run using a nondimensional, Boussinesq formulation
of the MHD equations (``reference_type=1``). The nondimensionalization
employed is as follows:

.. math::

   \begin{aligned}
       \mathrm{Length} &\rightarrow L &\;\;\;\; \mathrm{(Shell\; Depth)} \\
       \mathrm{Time} &\rightarrow   \frac{L^2}{\nu_o} &\;\;\;\; \mathrm{(Viscous\; Timescale)}\\
       \mathrm{Temperature} &\rightarrow \Delta T&\;\;\;\; \mathrm{(Temperature\; Contrast\; Across\; Shell)} \\
       \mathrm{Magnetic\; Field} &\rightarrow \sqrt{\hat{\rho}\mu\eta_o\Omega_0} \\
       \mathrm{Reduced\; Pressure} &\rightarrow \nu_o\Omega_0&\;\;\; (\mathrm{[Thermodynamic\; Pressure]}/\hat{\rho}),\end{aligned}

where :math:`\Omega_0` is the rotation rate of the frame, :math:`\hat{\rho}`
is the (constant) density of the fluid, :math:`\eta_o` is the magnetic diffusivity at the top of the domain (i.e., at :math:`r=r_o`), :math:`\nu_o`
is the kinematic viscosity at the top of the domain, and :math:`\mu` is the magnetic permeability. Note that in Gaussian units for vacuum, :math:`\mu=4\pi`. After nondimensionalizing, the following
nondimensional numbers appear in our equations:

.. math::

   \begin{aligned}
       Pr &=\frac{\nu_o}{\kappa_o}                          &\;\;\;\;\;\; \mathrm{Prandtl\; Number} \\
       Pm &=\frac{\nu_o}{\eta_o}                            &\;\;\;\;\;\; \mathrm{Magnetic\; Prandtl\; Number} \\
       E  &=\frac{\nu_o}{\Omega_0\,L^2}                   &\;\;\;\;\;\; \mathrm{Ekman\; Number} \\
       Ra &=\frac{\alpha g_o \Delta T\,L^3}{\nu_o\kappa_o}  &\;\;\;\;\;\; \mathrm{Rayleigh\; Number},\end{aligned}

where :math:`\alpha` is coefficient of thermal expansion, :math:`g_o`
is the gravitational acceleration at the top of the domain, and :math:`\kappa` is the thermal
diffusivity. Adopting this nondimensionalization is equivalent to
assigning the following to the functions :math:`\mathrm{f}_i(r)` and the constants :math:`c_i`:

.. math::

   \begin{aligned}
   \mathrm{f}_1(r) &\rightarrow 1\; &c_1 &\rightarrow \frac{2}{E} \\
   \mathrm{f}_2(r) &\rightarrow \left(\frac{r}{r_o}\right)^n \; &c_2 &\rightarrow \frac{Ra}{Pr} \\
   \mathrm{f}_3(r) &\rightarrow \tilde{\nu}(r)\; &c_3 &\rightarrow \frac{1}{E}\\
   \mathrm{f}_4(r) &\rightarrow 1\; &c_4 &\rightarrow \frac{1}{E\,Pm} \\
   \mathrm{f}_5(r) &\rightarrow \tilde{\kappa}(r)\; &c_5 &\rightarrow 1 \\
   \mathrm{f}_6(r) &\rightarrow 0\; &c_6 &\rightarrow \frac{1}{Pr}  \\
   \mathrm{f}_7(r) &\rightarrow \tilde{\eta}(r)\; &c_7 &\rightarrow \frac{1}{Pm} \\
    &\vdots &c_8&\rightarrow 0\\ 
    &\vdots &c_9&\rightarrow 0 \\
    &\vdots &c_{10}&\rightarrow 0 \\
    \mathrm{f}_{14}(r)&\rightarrow 0\; &c_{11}&\rightarrow 0.\end{aligned}

Here the tildes denote nondimensional radial profiles, e.g., :math:`\tilde{\nu}(r) = \nu(r)/\nu_o`. 

Our choice of :math:`\mathrm{f}_{14}(r)\rightarrow 0` sets the default atmosphere in non-dimensional Boussinesq to be neutrally stable. For other choices (i.e., convectively stable or unstable), one must use the custom-reference-state framework. 

Our choice of :math:`\mathrm{f}_2(r)` allows gravity to vary
with radius based on the value of the exponent :math:`n`, which has a
default value of :math:`0` in the code. Note also that our definition of
:math:`Ra` assumes fixed-temperature boundary conditions. We might specify fixed-flux boundary conditions and/or an internal heating
through a suitable choice :math:`c_{10}\mathrm{f}_6(r)`, in which case the
meaning of :math:`Ra` in our equation set changes, with :math:`Ra`
denoting a flux Rayleigh number instead. In addition, ohmic and viscous
heating, which do not appear in the Boussinesq formulation, are turned
off when this nondimensionalization is specified at runtime. When these
substitutions are made, Equations :eq:`momentum`-:eq:`induction`
transform as follows.

.. math::

   \begin{aligned}
       \left[\frac{\partial \boldsymbol{v}}{\partial t}  + \boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}  %advection
        + \frac{2}{E}\boldsymbol{\hat{z}}\times\boldsymbol{v} \right]  &= % Coriolis
       \frac{Ra}{Pr}\left(\frac{r}{r_o}\right)^n\Theta\,\boldsymbol{\hat{r}} % buoyancy
        - \frac{1}{E}\boldsymbol{\nabla}P % pressure
        + \frac{1}{E\,Pm}\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
        + \boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}& \\ 
       %
       %
       \left[\frac{\partial \Theta}{\partial t}  + \boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta \right] &=
       \frac{1}{Pr}\boldsymbol{\nabla}\cdot\left[\tilde{\kappa}(r)\boldsymbol{\nabla}\Theta\right] \\ % Diffusion
       %
       %
       \frac{\partial \boldsymbol{B}}{\partial t} &= \boldsymbol{\nabla}\times\left[\boldsymbol{v}\times\boldsymbol{B} - \frac{1}{Pm}\tilde{\eta}(r)\boldsymbol{\nabla}\times\boldsymbol{B}\right]\\
       \mathcal{D}_{ij} &= 2\tilde{\nu}(r)e_{ij} \\
       %
       %
       \boldsymbol{\nabla}\cdot\boldsymbol{v}&=0\\
       \boldsymbol{\nabla}\cdot\boldsymbol{B}&=0 \end{aligned}

Here :math:`\Theta` refers to the temperature (perturbation from the background) and :math:`P` to the reduced pressure (ratio of the thermodynamic pressure to the constant density). 

.. _dim_anelastic:

Dimensional Anelastic Formulation of the MHD Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

When run in dimensional, anelastic mode (cgs units; ``reference_type=2``),
the following values are assigned to the functions
:math:`\mathrm{f}_i` and the constants :math:`c_i`:

.. math::

   \begin{aligned}
       \mathrm{f}_1(r) &\rightarrow \hat{\rho}(r)\; &c_1 &\rightarrow 2\Omega_0 \\
       \mathrm{f}_2(r) &\rightarrow \frac{\hat{\rho}(r)}{c_P}g(r)\; &c_2 &\rightarrow 1 \\
       \mathrm{f}_3(r) &\rightarrow \nu(r)\; &c_3 &\rightarrow 1\\
       \mathrm{f}_4(r) &\rightarrow \hat{T}(r)\; &c_4 &\rightarrow \frac{1}{4\pi} \\
       \mathrm{f}_5(r) &\rightarrow \kappa(r)\; &c_5 &\rightarrow 1 \\
       \mathrm{f}_6(r) &\rightarrow \frac{Q(r)}{L_*}\; &c_6 &\rightarrow 1  \\
       \mathrm{f}_7(r) &\rightarrow \eta(r)\; &c_7 &\rightarrow 1 \\
       &\vdots &c_8&\rightarrow 1\\ 
       &\vdots &c_9&\rightarrow \frac{1}{4\pi} \\
       &\vdots &c_{10}&\rightarrow L_* \\
       \mathrm{f}_{14}(r)&\rightarrow \frac{d\hat{S}}{dr }&c_{11}&\rightarrow 1.\end{aligned}

Here :math:`\hat{\rho}(r)`, :math:`\hat{T}(r)`, and :math:`d\hat{S}/dr` are the spherically symmetric, time-independent reference-state
density, temperature, and entropy gradient, respectively. The thermal variables satisfy the
linearized equation of state

.. math:: \frac{P}{\hat{P}}= \frac{T}{\hat{T}} + \frac{\rho}{\hat{\rho}}

:math:`g(r)` is the gravitational
acceleration, :math:`c_P` is the specific heat at constant pressure, and
:math:`\Omega_0` is the frame rotation rate. The viscous, thermal, and
magnetic diffusivities (also assumed to be spherically symmetric and time-independent) are given by :math:`\nu(r)`, :math:`\kappa(r)`, and
:math:`\eta(r)`, respectively. Note that the entropy gradient term :math:`f_{14}(r)v_r` is only used in Equation :eq:`theta_evol` if ``advect_reference_state=.true.``. Finally, :math:`Q(r)` is an internal heating
function; it might represent radiative heating or heating due to nuclear
fusion, for instance. In our convention, the volume integral of :math:`\mathrm{f}_6(r)` equals unity, and :math:`c_{10}` equals the ``luminosity`` or ``heating_integral`` :math:`L_*` specified in the main_input file. When using a custom reference state, this allows easy adjustment of the luminosity using the **override_constants** formalism, e.g.,

``override_constants(10) = T``

``ra_constants(10) = 3.846d33``

specified in the in the ``reference_namelist``.

Note that in the anelastic formulation, the
thermal variable :math:`\Theta` is interpreted as the entropy perturbation,
rather than the temperature perturbation. When these substitutions are made,
Equations :eq:`momentum`-:eq:`induction` transform as follows.

.. math::

   \begin{aligned}
       \hat{\rho}(r)\left[\frac{\partial \boldsymbol{v}}{\partial t} +\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}  %advection
       +2\Omega_0\boldsymbol{\hat{z}}\times\boldsymbol{v} \right]  =\; % Coriolis
       &\frac{\hat{\rho}(r)}{c_P}g(r)\Theta\,\boldsymbol{\hat{r}} % buoyancy
       +\hat{\rho}(r)\boldsymbol{\nabla}\left(\frac{P}{\hat{\rho}(r)}\right) % pressure
       \\
       &+\frac{1}{4\pi}\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
       +\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\\
       %
       %
       \hat{\rho}(r)\,\hat{T}(r)\left[\frac{\partial \Theta}{\partial t} +\boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta + v_r\frac{d\hat{S}}{dr}\right] =\;
       &\boldsymbol{\nabla}\cdot\left[\hat{\rho}(r)\,\hat{T}(r)\,\kappa(r)\,\boldsymbol{\nabla}\Theta \right] % diffusion
       +Q(r)   % Internal heating
       \\ &+\Phi(r,\theta,\phi)
       +\frac{\eta(r)}{4\pi}\left[\boldsymbol{\nabla}\times\boldsymbol{B}\right]^2\\ % Ohmic Heating
       %
       %
       \frac{\partial \boldsymbol{B}}{\partial t} =\; &\boldsymbol{\nabla}\times\left[\,\boldsymbol{v}\times\boldsymbol{B}-\eta(r)\boldsymbol{\nabla}\times\boldsymbol{B}\,\right] \\
       %
       %
       \mathcal{D}_{ij} =\; &2\hat{\rho}(r)\,\nu(r)\left[e_{ij}-\frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)\delta_{ij}\right] \\
       %
       %
       \Phi(r,\theta,\phi) =\; &2\,\hat{\rho}(r)\nu(r)\left[e_{ij}e_{ij}-\frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)^2\right] \\
       %
       %
       \boldsymbol{\nabla}\cdot\left[\hat{\rho}(r)\,\boldsymbol{v}\right] =\; &0 \\
       \boldsymbol{\nabla}\cdot\boldsymbol{B} =\; &0. \end{aligned}

.. _nondim_anelastic:

Nondimensional Anelastic MHD Equations
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To run in nondimensional anelastic mode, you must set
``reference_type=3`` in the Reference_Namelist. The reference state is
assumed to be polytropic with a :math:`\frac{1}{r^2}` profile for
gravity. When this mode is
active, the following nondimensionalization is used
(following `Heimpel et al., 2016, Nat. Geo., 9, 19 <https://www.nature.com/articles/ngeo2601/>`_ ):

.. math::

   \begin{aligned}
       \mathrm{Length} &\rightarrow L \equiv r_o - r_i &\;\;\;\; \mathrm{(Shell\; Depth)} \\
       \mathrm{Time} &\rightarrow   \frac{1}{\Omega_0} &\;\;\;\; \mathrm{(Rotational\; Timescale)}\\
       \mathrm{Temperature} &\rightarrow T_o\equiv\hat{T}(r_o)&\;\;\;\; \mathrm{(Reference\; Temperature\; at\; Upper\; Boundary)} \\
       \mathrm{Density} &\rightarrow \rho_o\equiv\hat{\rho}(r_o)&\;\;\;\; \mathrm{(Reference\; Density\; at\; Upper\; Boundary)} \\
       \mathrm{Entropy} &\rightarrow \Delta{s}&\;\;\;\; \mathrm{(Entropy\; Constrast\; Across\; Shell)} \\
       \mathrm{Magnetic~Field} &\rightarrow \sqrt{\hat{\rho}_o\mu\eta_o\Omega_0} \\
       \mathrm{Pressure} &\rightarrow \rho_oL^2\Omega_0^2.\end{aligned}

When run in this mode, Rayleigh employs a polytropic background state,
with an assumed :math:`\frac{1}{r^2}` variation in gravity. These
choices result in the functions :math:`\mathrm{f}_i` and the constants
:math:`c_i` (tildes indicate nondimensional reference-state variables):

.. math::

   \begin{aligned}
       \mathrm{f}_1(r) &\rightarrow \tilde{\rho}(r)\; &c_1 &\rightarrow 2 \\
       \mathrm{f}_2(r) &\rightarrow \tilde{\rho}(r)\frac{r_\mathrm{max}^2}{r^2}\; &c_2 &\rightarrow \mathrm{Ra}^* \\
       \mathrm{f}_3(r) &\rightarrow \tilde{\nu}(r)\; &c_3 &\rightarrow 1\\
       \mathrm{f}_4(r) &\rightarrow \tilde{T}(r)\; &c_4 &\rightarrow \frac{\mathrm{E}}{\mathrm{Pm}} \\
       \mathrm{f}_5(r) &\rightarrow \tilde{\kappa}(r)\; &c_5 &\rightarrow \mathrm{E} \\
       \mathrm{f}_6(r) &\rightarrow \frac{\tilde{Q}(r)}{L_*}; &c_6 &\rightarrow \frac{\mathrm{E}}{\mathrm{Pr}}  \\
       \mathrm{f}_7(r) &\rightarrow \tilde{\eta}(r) \; &c_7 &\rightarrow \frac{\mathrm{E}}{\mathrm{Pm}} \\
       &\vdots &c_8&\rightarrow \frac{\mathrm{E}\,\mathrm{Di}}{\mathrm{Ra}^*}\\ 
       &\vdots &c_9&\rightarrow  \frac{\mathrm{E}^2\,\mathrm{Di}}{\mathrm{Pm}^2\mathrm{Ra}^*}\\
       &\vdots &c_{10}&\rightarrow L_* \\
       \mathrm{f}_{14}(r)&\rightarrow 0&c_{11}&\rightarrow 0.\end{aligned}

As in the Boussinesq case, the nondimensional diffusivities are defined according to, e.g., :math:`\tilde{\nu}(r) \equiv \nu(r)/\nu_o`. The nondimensional heating :math:`\tilde{Q}(r)` is defined such that its volume integral equals the nondimensional ``luminosity`` or ``heating_integral`` set in the *main_input* file. As in the dimensional anelastic case, the volume integral of :math:`\mathrm{f}_6(r)` equals unity, and :math:`\mathrm{c}_{10} = L_*`. The unit for luminosity in this nondimensionalization (to get a dimensional luminosity from the nondimensional :math:`L_*`) is :math:`\rho_oL^3T_o\Delta s\Omega_0`. 

Two new nondimensional numbers appear in our equations, in addition to those defined for the Boussinesq case. :math:`\mathrm{Di}`, the
dissipation number, is defined by

.. math::
   :label: Di

       \mathrm{Di}= \frac{g_o\,\mathrm{L}}{c_\mathrm{P}\,T_o},

where :math:`g_o` and :math:`T_o` are the gravitational acceleration
and temperature at the outer boundary respectively. Once more, the
thermal anomaly :math:`\Theta` should be interpreted as (nondimensional) entropy. The symbol :math:`\mathrm{Ra}^*` is the modified Rayleigh number,
given by

.. math::
   :label: Ra

   \mathrm{Ra}^*=\frac{g_o}{c_\mathrm{P}\Omega_0^2}\frac{\Delta s}{L}

We thus arrive at the following nondimensionalized equations:

.. math::

   \begin{aligned}
       \tilde{\rho}(r)\left[\frac{\partial \boldsymbol{v}}{\partial t}  + \boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}  %advection
        + 2\boldsymbol{\hat{z}}\times\boldsymbol{v}\right]  =\; % Coriolis
       &\mathrm{Ra}^*\tilde{\rho}(r)\left(\frac{r_\mathrm{max}^2}{r^2}\right)\Theta\,\boldsymbol{\hat{r}} % buoyancy
        + \tilde{\rho}(r)\boldsymbol{\nabla}\left(\frac{P}{\tilde{\rho}(r)}\right) % pressure
       \\
       &+ \frac{\mathrm{E}}{\mathrm{Pm}}\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
        + \mathrm{E}\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}}\\
       %
       %
       \tilde{\rho}(r)\,\tilde{T}(r)\left[\frac{\partial \Theta}{\partial t} + \boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta\right] =\;
       &\frac{\mathrm{E}}{\mathrm{Pr}}\boldsymbol{\nabla}\cdot\left[\tilde{\kappa}(r)\tilde{\rho}(r)\,\tilde{T}(r)\,\boldsymbol{\nabla}\Theta \right] % diffusion
       + \tilde{Q}(r)   % Internal heating
       \\
       &+ \frac{\mathrm{E}\,\mathrm{Di}}{\mathrm{Ra}^*}\Phi(r,\theta,\phi)
       + \frac{\mathrm{Di\,E^2}}{\mathrm{Pm}^2\mathrm{Ra}^*}\tilde{\eta}(r)|\boldsymbol{\nabla}\times\boldsymbol{B}|^2 \\ % Ohmic Heating
       %
       %
       \frac{\partial \boldsymbol{B}}{\partial t} =\;  &\boldsymbol{\nabla}\times\left[\,\boldsymbol{v}\times\boldsymbol{B}-\frac{\mathrm{E}}{\mathrm{Pm}}\tilde{\eta}(r)\boldsymbol{\nabla}\times\boldsymbol{B}\,\right] \\
       %
       %
       \mathcal{D}_{ij} =\; &2\tilde{\rho}(r)\tilde{\nu}(r)\left[e_{ij} - \frac{1}{3}\boldsymbol{\nabla}\cdot\boldsymbol{v}\right] \\
       %
       %
       \Phi(r,\theta,\phi) =\; &2\tilde{\rho}(r)\tilde{\nu}(r)\left[e_{ij}e_{ij} - \frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)^2\right] \\
       %
       %
       \boldsymbol{\nabla}\cdot\left[\tilde{\rho}(r)\boldsymbol{v}\right]=\; &0 \\
       \boldsymbol{\nabla}\cdot\boldsymbol{B} =\; &0. \end{aligned}


.. _streamfunctions:

The Streamfunction Formulation
------------------------------

The velocity field in Rayleigh is evolved subject to the solenoidal constraint

.. math::
   :label: solenoidal

       \boldsymbol{\nabla}\cdot\left[\mathrm{f}_1(r)\,\boldsymbol{v}\right] = 0.

This is accomplished by casting :math:`\mathrm{f}_1\boldsymbol{v}` in terms of streamfunctions such that

.. math::
   :label: streamdecomp

       \mathrm{f_1}\,\boldsymbol{v} = \boldsymbol{\nabla}\times\boldsymbol{\nabla}\times\left( W\,\boldsymbol{\hat{r}}\right)+\boldsymbol{\nabla}\times\left( Z\,\boldsymbol{\hat{r}}\right),  
       
where *W* and *Z* are referred to as the poloidal and toroidal stream functions respectively.  Rather than evolving the three components of :math:`\boldsymbol{v}` directly, the momentum equations are cast in terms of these variables before advancing the timestep.  The velocity components are related to the streamfunctions via the relations:

.. math::
   :label: vrstream
   
       \mathrm{f_1}v_r = - \frac{1}{r^2\mathrm{sin}\theta}\frac{\partial}{\partial\theta}\left(\mathrm{sin}\theta\frac{\partial W}{\partial\theta} \right)-\frac{1}{r^2\mathrm{sin}^2\theta}\frac{\partial^2 W}{\partial\phi^2},   

.. math::
   :label: vtstream
   
        \mathrm{f_1}\,v_{\theta} = \frac{1}{r}\frac{\partial^2 W}{\partial r\partial\theta}+ \frac{1}{r\mathrm{sin}\theta}\frac{\partial Z}{\partial\phi},
    
and

.. math::
   :label: vpstream

       \mathrm{f_1}v_{\phi} = \frac{1}{r\mathrm{sin}\theta}\frac{\partial^2 W}{\partial r\partial\phi} - \frac{1}{r}\frac{\partial Z}{\partial\theta}.
       
When the velocity field and streamfunctions are projected onto a spherical harmonic basis, two additional useful relations are given by

.. math::
   :label: vrWstream
   
       \left[\mathrm{f_1}\,v_r\right]_\ell^m = \frac{\ell(\ell+1)}{r^2}W_\ell^m
       
and

.. math::
   :label: curlstream
   
       \left[ \left\{\boldsymbol{\nabla}\times\left(\mathrm{f_1}\,\boldsymbol{v}\right)\right\}_r \right]_\ell^m = \frac{\ell(\ell+1)}{r^2}Z_\ell^m.
       

The equations that are solved are then equations for the radial component of the momentum equation (21): 

.. math::
   :label: Radial Component of the Momentum Equation

    \tiny\begin{aligned}    
        \frac{\partial}{\partial t}\left(\overline{\rho}v_{r}\right)_{l}^{m}=\frac{\ell\left(\ell+1\right)}{r^2}\frac{\partial  W_{l}^{m}}{\partial t}=-\rho\frac{\partial{P_{l}^{m}}}{\partial r}-\overline{g}\left(\frac{\partial\overline{\rho}}{\partial \Theta}\right)_{p,\xi}
        \\
        +\frac{2\Omega}{r}\left[im\frac{\partial W_{l}^{m}}{\partial r}+\left(\ell+2\right)d_{l}^{m}Z_{l+1}^{m}-\left(\ell-1\right)d_{l}^{m}Z_{\ell-1}^{m}\right]
        \\
        +\frac{\overline{\nu}\ell\left(\ell+
        1\right)}{r^2}\left[\frac{{\partial^2 W_{l}^{m}}}{{\partial r^2}}+\left(2 h_{\nu}-\frac{h_{\rho}}{3}\right) \frac{{\partial W_{l}^{m}}}{{\partial r}}\right.
        \\
        \left.-\left(\frac{4}{3}\left(\left(\frac{h_{\rho}}{r}+\frac{dh_{\rho}}{dr}\right)+h_{\nu}\left(\frac{3}{r}+h_{\rho}\right)\right)+\frac{\ell\left(\ell+1\right)}{r^2}\right)W_{l}^{m}\right]
        \\
        +\frac{FLMW1_l^m}{r^2}
    \end{aligned}


    

the radial component of the curl of the momentum equation (22)

.. math:: 
    :label: Radial Component of the Curl of the Momentum Equation

     \tiny\begin{aligned}
     \frac{\partial \left(\nabla\times \overline{\rho}\bf{v}\right)_{r,l}^{m}}{\partial t}=\frac{\ell\left(\ell+1\right)}{r^2}\frac{\partial Z_{l}^{m}}{\partial t}=\frac{2\Omega}{r^2}\left[im Z_{l}^{m} + \right .
    \\
    \left . \ell\left(\ell+1\right)d_{l+1}^{m}\left(\frac{\partial W_{l+1}^{m}}{\partial r}+
    \frac{\left(l+1\right)}{r^2}W_{l+1}^{m}\right)+\left(\ell+1\right)\left(\ell-1\right)d_{l}^m \left(\frac{\partial W_{l-1}^m}{\partial r}-\frac{\ell}{r}W_{l-1}^m\right)\right]
    \\
    +\frac{\nu\ell\left(\ell+1\right)}{r^2}\left[\frac{\partial^2 Z_{l}^m}{\partial r^2}+\left(h_{\nu}-h_{\rho}\right)\frac{\partial Z_{l}^m}{\partial r} 
    -\left(\frac{2h_{\rho}}{r}+\frac{dh_{\rho}}{dr}+h_{\nu}\left(\frac{2}{r}+h_{\rho}\right)+\frac{\ell\left(\ell+1\right)}{r^2}\right)Z_{l}^m\right]
    \\
    \left(\ell+1\right)C_{l}^m FLMW3_{l-1}^m-\ell C_{l+1}^m FLMW3_{l+1}^m-im FLMW2_{l}^m
    \end{aligned}

and the Horizontal Divergence of the Momentum Equation (23)

.. math:: 
    :label: Horizontal Divergence of the Momentum Equation

    \tiny\begin{aligned}
    \frac{\partial \left(\nabla\cdot\overline{\rho} \bf{v}\right)_{l}^m}{\partial t}=-\frac{\ell\left(\ell+1\right)}{r^2}\frac{\partial}{\partial t}\left(\frac{\partial W_l^m}{\partial r}\right)=\frac{\ell\left(\ell+1\right)}{r^2}\overline{\rho}P_l^m+ 
    \\
    \frac{2\Omega}{r^2}\left[\ell\left(\ell+2\right)d_{l+1}^mZ_{l+1}^m +\left(\ell+1\right)\left(\ell-1\right)d_l^m Z_{l-1}^m-im\left(\frac{\partial W_l^m}{\partial r}+\frac{\ell\left(\ell+1\right)}{r}W_l^m\right)\right]
    \\
    +\frac{\nu\ell\left(\ell+1\right)}{r^2}\left[-\frac{\partial^3{W_l^m}}{\partial r^3} -\left(h_{\nu}-h_{\rho}\right)\frac{\partial^2 W_l^m}{r^2}\right.
    \\
    \left. +\left(\frac{2h_\rho}{r}+\frac{\partial h_{\rho}}{\partial r}+h_{\nu}\left(\frac{2}{r}+h_{\rho}\right)+\frac{\ell\left(\ell+1\right)}{r^2}\right)\frac{\partial W_l^m}{\partial r}\right .
    \\
    \left . -\frac{\ell\left(\ell+1\right)}{r^2}\left(h_{\nu}+\frac{2}{3}h_{\rho}+\frac{2}{r}\right)W_l^m\right]
    \\
    +\left[\left(\ell+1\right)C_l^mFLMW2_{l-1}^m-\ell C_{l+1}^m FLMW2_{l+1}^m+im FLMW3_l^m\right]
    \end{aligned}


A similar decomposition is performed on the magnetic field to ensure it remains divergence free.  In that case, the magnetic field is projected onto flux functions such that

.. math::
   :label: fluxdecomp

       \boldsymbol{B} = \boldsymbol{\nabla}\times\boldsymbol{\nabla}\times\left( C\,\boldsymbol{\hat{r}}\right)+\boldsymbol{\nabla}\times\left( A\,\boldsymbol{\hat{r}}\right),   
       
where *C* and *A* are the poloidal and toroidal flux functions respectively.  Similar to the velocity field, the components of :math:`\boldsymbol{B}` satisfy

.. math::
   :label: Brstream
   
       B_r = - \frac{1}{r^2\mathrm{sin}\theta}\frac{\partial}{\partial\theta}\left(\mathrm{sin}\theta\frac{\partial C}{\partial\theta} \right)-\frac{1}{r^2\mathrm{sin}^2\theta}\frac{\partial^2 C}{\partial\phi^2},   

.. math::
   :label: Btstream
   
        B_{\theta} = \frac{1}{r}\frac{\partial^2 C}{\partial r\partial\theta}+ \frac{1}{r\mathrm{sin}\theta}\frac{\partial A}{\partial\phi},
    
.. math::
   :label: Bpstream

       B_{\phi} = \frac{1}{r\mathrm{sin}\theta}\frac{\partial^2 C}{\partial r\partial\phi} - \frac{1}{r}\frac{\partial A}{\partial\theta},

.. math::
   :label: BrCflux
   
       \left[B_r\right]_\ell^m = \frac{\ell(\ell+1)}{r^2}C_\ell^m,
       
and

.. math::
   :label: curlflux
   
       \left[ \left\{\boldsymbol{\nabla}\times\boldsymbol{B}\right\}_r \right]_\ell^m = \frac{\ell(\ell+1)}{r^2}A_\ell^m.


.. _pseudospectral:

The equations for C and A, which are solved by Rayleigh are then the Radial Component of the Magnetic Induction Equation (30): 

.. math:: 
    :label: Radial Component of the Magnetic Induction Equation

    \tiny\begin{aligned}    
    \frac{\partial B_{r,l}^m}{\partial t}=\frac{\ell\left(\ell+1\right)}{r^2}\frac{\partial C_l^m}{\partial t}  =\overline{\eta}\frac{\ell\left(\ell+1\right)}{r^2}\left(\frac{\partial^2 C_l^m}{\partial r^2}-\frac{\ell\left(\ell+1\right)}{r^2}C_l^m\right)
    \\
    +\left[\left(\ell+1\right)d_{l}^mFLMB3_{l-1}^m-\ell d_{l+1}^mFLMB3_{l+1}^m-imFLMB2_l^m\right]   
    \end{aligned}

and the radial component of the curl of the magnetic induction equation ():

.. math:: 
    :label: Radial Component of the Curl of the Magnetic Induction Equation

    \tiny\begin{aligned}
    \frac{\partial\left(\nabla\times B\right)_{r,l}^m}{\partial t}=\frac{\ell\left(\ell+1\right)}{r^2}\frac{\partial A_l^m}{\partial t}=\overline{\eta}\frac{\ell\left(\ell+1\right)}{r^2}\left(\frac{\partial^2 A_l^m}{\partial r^2}+h_{\eta}\frac{\partial A_l^m}{\partial r}-\frac{\ell\left(\ell+1\right)}{r^2} A_l^m \right)+
    \\
    \frac{1}{r^2}\left[\frac{\ell\left(\ell+1\right)}{r^2}FLMB1_l^m+\frac{\partial}{\partial r}\left(r^2\left(\left(\ell+1\right)d_l^mFLMB2_{l-1}^m-\ell d_{l+1}^mFLMB2_{l+1}^m+imFLMB3_l^m\right)\right)\right]
    \end{aligned}


Where the "FLM*" terms refer to nonlinear terms, defined as: 

.. math::
    :label: FLMW1

    \scriptsize FLMW1=r^2\left[-\left(\nabla\cdot\overline{\rho}\bf{v}\bf{v}\right)_r+\frac{1}{\mu}\left(\left(\nabla\times\bf{B}\right)\times\bf{B}\right)_r +\Omega^2\rho r\sin^2\theta\right]_l^m

.. math::
    :label: FLMW2

    \scriptsize FLMW2=\left[\frac{-\nabla\cdot\left(\overline{\rho}\bf{v}\bf{v}\right)_{\phi}}{r\sin\theta}+\frac{1}{\mu}\frac{\left(\left(\nabla\times\bf{B}\right)\times\bf{B}\right)_{\theta}}{r\sin\theta}+\Omega^2\rho\cos\theta\right]_l^m

.. math:: 
    :label: FLMW3

    \scriptsize FLMW3=\left[\frac{-\left(\nabla\cdot\overline{\rho}\bf{v}\bf{v}\right)_{\phi}}{r\sin\theta}+\frac{1}{\mu}\frac{\left(\left(\nabla\times\bf{B}\right)\times\bf{B}\right)_{\phi}}{r\sin\theta}\right]_l^m

.. math:: 
    :label: FLMB1

    \scriptsize FLMB1=\left[r^2\left(\bf{v}\times\bf{B}\right)_r\right]_l^m

.. math:: 
    :label: FLMB2

    \scriptsize FLMB2=\left[\frac{\left(\bf{v}\times\bf{B}\right)_{\theta}}{r\sin\theta}\right]_l^m    

.. math:: 
    :label: FLMB3 

    \scriptsize FLMB3=\left[\frac{\left(\bf{v}\times\bf{B}\right)_{\phi}}{r\sin\theta}\right]_l^m
       

The Pseudospectral Approach
---------------------------

Section needed.

.. _performance:

Parallelization and Performance
-------------------------------

Section needed.
