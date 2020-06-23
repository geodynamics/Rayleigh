.. raw:: latex

   \clearpage

.. _physics:

Physics Controls
================

Rayleigh solves the MHD equations in spherical geometry under the
Boussinesq and anelastic approximations. Both the equations that
Rayleigh solves and its diagnostics can be formulated either
dimensionally or nondimensionally. A nondimensional Boussinesq
formulation, as well as dimensional and non-dimensional anelastic
formulations (based on a polytropic reference state) are provided as
part of Rayleigh.

In this section, we present the equation sets solved when running in
each of these three modes, and discuss the relevant control parameters
for each mode. We also discuss the boundary conditions available in
Rayleigh and those namelist variables that can be used to modify the
code’s behavior in any of these three modes.

Anelastic Mode (dimensional)
----------------------------

**Example Input: Rayleigh/input_examples/main_input_sun**

When run in dimensional, anelastic mode, **reference_type=2** must be
specified in the Reference_Namelist. In that case, Rayleigh solves the
following form of the MHD equations:

.. math::

   \begin{aligned}
   \hat{\rho}(r)\left[\frac{\partial \boldsymbol{v}}{\partial t} +\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}  %advection
                                                            +2\Omega_0\boldsymbol{\hat{z}}\times\boldsymbol{v} \right]  &= % Coriolis
                                                            \frac{\hat{\rho}(r)}{c_P}g(r)\Theta\,\boldsymbol{\hat{r}} % buoyancy
                                                            +\hat{\rho}(r)\boldsymbol{\nabla}\left(\frac{P}{\hat{\rho}(r)}\right) % pressure
                                                            +\frac{1}{4\pi}\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
                                                            +\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}} \;\;\; &\mathrm{Momentum}\\
   %
   %
   \hat{\rho}(r)\,\hat{T}(r)\left[\frac{\partial \Theta}{\partial t} +\boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta \right] &=
                                                \boldsymbol{\nabla}\cdot\left[\hat{\rho}(r)\,\hat{T}(r)\,\kappa(r)\,\boldsymbol{\nabla}\Theta \right] % diffusion
                                                +Q(r)   % Internal heating
                                                +\Phi(r,\theta,\phi)
                                                +\frac{\eta(r)}{4\pi}\left[\boldsymbol{\nabla}\times\boldsymbol{B}\right]^2 &\mathrm{Thermal\; Energy}\\ % Ohmic Heating
   %
   %
   \frac{\partial \boldsymbol{B}}{\partial t} &= \boldsymbol{\nabla}\times\left(\,\boldsymbol{v}\times\boldsymbol{B}-\eta(r)\boldsymbol{\nabla}\times\boldsymbol{B}\,\right) &\mathrm{Induction} \\
   %
   %
   \mathcal{D}_{ij} &= 2\hat{\rho}(r)\,\nu(r)\left[e_{ij}-\frac{1}{3}\boldsymbol{\nabla}\cdot\boldsymbol{v}\right] &\mathrm{Viscous\; Stress\; Tensor}\\
   %
   %
   \Phi(r,\theta,\phi) &= 2\,\hat{\rho}(r)\nu(r)\left[e_{ij}e_{ij}-\frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)^2\right] &\mathrm{Viscous\; Heating} \\
   %
   %
   \boldsymbol{\nabla}\cdot\left(\hat{\rho}(r)\,\boldsymbol{v}\right)&=0 &\mathrm{Solenoidal\; Mass\; Flux}\\
   \boldsymbol{\nabla}\cdot\boldsymbol{B}&=0 &\mathrm{Solenoidal\; Magnetic\; Field}\end{aligned}

Here, :math:`\hat{\rho}` and :math:`\hat{T}` are the reference-state
density and temperature respectively. :math:`g` is the gravitational
acceleration, :math:`c_P` is the specific heat at constant pressure, and
:math:`\Omega_0` is the frame rotation rate. The velocity field vector
is denoted by :math:`\boldsymbol{v}`, the magnetic field vector by
:math:`\boldsymbol{B}`, and the pressure by :math:`P`. The thermal
anomoly is denoted by :math:`\Theta` and should be interpreted is as
entropy :math:`s` in this formulation. The thermal variables satisfy the
linearized equation of state

.. math:: \frac{P}{\hat{P}}= \frac{T}{\hat{T}} + \frac{\rho}{\hat{\rho}}

The kinematic viscosity, thermal diffusivity, and magnetic diffusivity
are given by :math:`\nu`, :math:`\kappa`, and :math:`\eta` respectively.
Finally, :math:`Q(r)` is an internal heating function; it might
represent radiative heating or heating due to nuclear fusion, for
instance.

When running in anelastic mode, the **reference_type** variable in the
Reference_Namelist must be set to 2.

Moreover, certain variables in the **Reference_Namelist** and the
**Transport_Namelist** must be specified. The Reference_Namelist
variables are described in Table table_anelastic_ and the Transport_Namelist
variables are described in Table table_anelastic_trans_. Default values
indicated in brackets.

.. _table_anelastic:

.. centered:: **Table. Anelastic.**

Variables in the Reference_Namelist that
must be specified when running in dimensional anelastic mode. In
addition, reference_type=2 must also be specified.

   +-----------------------------------+-----------------------------------+
   | Variable                          | Description                       |
   +===================================+===================================+
   | poly_n [0]                        | polytropic index                  |
   |                                   | (:math:`P\propto\rho^n`)          |
   +-----------------------------------+-----------------------------------+
   | poly_Nrho [0]                     | number of density scaleheights    |
   |                                   | spanning the domain               |
   +-----------------------------------+-----------------------------------+
   | poly_mass [0]                     | mass interior to :math:`rmin`     |
   +-----------------------------------+-----------------------------------+
   | poly_rho_i [0]                    | density at rmin,                  |
   |                                   | :math:`\rho(r=rmin)`              |
   +-----------------------------------+-----------------------------------+
   | pressure_specific_heat [0]        | specific heat at constant         |
   |                                   | pressure                          |
   +-----------------------------------+-----------------------------------+
   | angular_velocity [1.0]            | cyclic frequency of the rotating  |
   |                                   | frame                             |
   +-----------------------------------+-----------------------------------+

   .. _table_anelastic_trans:


.. centered:: **Table. Anelastic Transport.**

Variables in the Transport_Namelist
that must be specified when running in dimensional anelastic mode. In
addition, reference_type=2 must also be specified in the
Reference_Namelist.

   +-----------------------------------+-----------------------------------+
   | Variable                          | Description                       |
   +===================================+===================================+
   | nu_top [1.0]                      | kinematic viscosity at rmax,      |
   |                                   | :math:`\nu(rmax)`                 |
   +-----------------------------------+-----------------------------------+
   | nu_type [1]                       | determines whether :math:`\nu` is |
   |                                   | constant with radius (1) or       |
   |                                   | varies with density (2)           |
   +-----------------------------------+-----------------------------------+
   | nu_power [0.0]                    | exponent in :                     |
   |                                   | :math:`\nu(r) = \left( \frac{\rho |
   |                                   | (r)}{\rho(r=rmax)} \right)^       |
   |                                   | {nu\_power}`;                     |
   |                                   | use with nu_type=2                |
   +-----------------------------------+-----------------------------------+
   | kappa_top [1.0]                   | thermal diffusivity at rmax,      |
   |                                   | :math:`\kappa(rmax)`              |
   +-----------------------------------+-----------------------------------+
   | kappa_type [1]                    | determines whether :math:`\kappa` |
   |                                   | is constant with radius (1) or    |
   |                                   | varies with density (2)           |
   +-----------------------------------+-----------------------------------+
   | kappa_power [0.0]                 | exponent in :                     |
   |                                   | :math:`\kappa(r) = \left( \frac{\ |
   |                                   | rho(r)}{\rho(r=rmax)} \right)^    |
   |                                   | {kappa\_power}`;                  |
   |                                   | use with kappa_type=2             |
   +-----------------------------------+-----------------------------------+
   | eta_top [1.0]                     | magnetic diffusivity at rmax,     |
   |                                   | :math:`\eta(rmax)`                |
   +-----------------------------------+-----------------------------------+
   | eta_type [1]                      | determines whether :math:`\eta`   |
   |                                   | is constant with radius (1) or    |
   |                                   | varies with density (2)           |
   +-----------------------------------+-----------------------------------+
   | eta_power [0.0]                   | exponent in :                     |
   |                                   | :math:`\eta(r) = \left( \frac{    |
   |                                   | \rho(r)}{\rho(r=rmax)} \right)^   |
   |                                   | {eta\_power}`;                    |
   |                                   | use with eta_type=2               |
   +-----------------------------------+-----------------------------------+

The polytropic reference state is the same as that used in the
benchmarks and is described in detail in Jones et al. (2011).

See the example input file **main_input_sun** for a an example of how to
run a solar-like model using Rayleigh’s dimensional, anelastic
formulation.

.. raw:: latex

   \clearpage

.. _physics_boussinesq_nondimensional:

Boussinesq Mode (nondimensional)
--------------------------------

**Example Input: Rayleigh/input_examples/c2001_case1_input**

When run in nondimensional Boussinesq mode, **reference_type=1** must be
specified in the Reference_Namelist. In that case, Rayleigh employs the
nondimensionalization

.. math::

   \begin{aligned}
   \mathrm{Length} &\rightarrow L &\;\;\;\; \mathrm{(Shell\; Depth)} \\
   \mathrm{Time} &\rightarrow   \frac{L^2}{\nu} &\;\;\;\; \mathrm{(Viscous\; Timescale)}\\
   \mathrm{Temperature} &\rightarrow \Delta T&\;\;\;\; \mathrm{(Temperature\; Contrast\; Across\; Shell)} \\
   \mathrm{Magnetic~Field} &\rightarrow \sqrt{\rho\mu\eta\Omega_0},\end{aligned}

where :math:`\Omega_0` is the rotation rate of the frame, :math:`\rho`
is the (constant) density of the fluid, :math:`\mu` is the magnetic
permeability, :math:`\eta` is the magnetic diffusivity, and :math:`\nu`
is the kinematic viscosity. After nondimensionalizing, the following
nondimensional numbers appear in our equations

.. math::

   \begin{aligned}
   Pr &=\frac{\nu}{\kappa}                          &\;\;\;\;\;\; \mathrm{Prandtl\; Number} \\
   Pm &=\frac{\nu}{\eta}                            &\;\;\;\;\;\; \mathrm{Magnetic\; Prandtl\; Number} \\
   E  &=\frac{\nu}{\Omega_0\,L^2}                   &\;\;\;\;\;\; \mathrm{Ekman\; Number} \\
   Ra &=\frac{\alpha g_0 \Delta T\,L^3}{\nu\kappa}  &\;\;\;\;\;\; \mathrm{Rayleigh\; Number}, \\\end{aligned}

where :math:`\alpha` is coefficient of thermal expansion, :math:`g_0`
is the gravitational acceleration at the top of the domain, and
:math:`\kappa` is the thermal diffusivity.

In addition, ohmic and viscous heating, which do not appear in the
Boussinesq formulation, are turned off when this nondimensionalization
is specified at runtime. Rayleigh solves the following equations when
running in nondimensional Boussinesq mode:

.. math::

   \begin{aligned}
   \left[\frac{\partial \boldsymbol{v}}{\partial t} +\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}  %advection
                                                            +\frac{2}{E}\boldsymbol{\hat{z}}\times\boldsymbol{v} \right]  &= % Coriolis
                                                            \frac{Ra}{Pr}\left(\frac{r}{r_o}\right)^n\Theta\,\boldsymbol{\hat{r}} % buoyancy
                                                            -\frac{1}{E}\boldsymbol{\nabla}P % pressure
                                                            +\frac{1}{E\,Pm}\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
                                                            +\boldsymbol{\nabla}^2\boldsymbol{v} \;\;\; &\mathrm{Momentum}\\
   %
   %
   \left[\frac{\partial \Theta}{\partial t} +\boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta \right] &=
                                                \frac{1}{Pr}\boldsymbol{\nabla}^2\Theta  &\mathrm{Thermal\; Energy}\\ % Diffusion
   %
   %
   \frac{\partial \boldsymbol{B}}{\partial t} &= \boldsymbol{\nabla}\times\left(\,\boldsymbol{v}\times\boldsymbol{B}\right)+\frac{1}{Pm}\boldsymbol{\nabla}^2\boldsymbol{B} &\mathrm{Induction} \\
   %
   %
   %
   %
   %
   %
   \boldsymbol{\nabla}\cdot\boldsymbol{v}&=0 &\mathrm{Solenoidal\; Velocity\; Field}\\
   \boldsymbol{\nabla}\cdot\boldsymbol{B}&=0 &\mathrm{Solenoidal\; Magnetic\; Field},\end{aligned}

where :math:`r_0 \equiv rmax`. In this formulation, :math:`\Theta`
should be interpreted as the temperature perturbation :math:`T`. Those
Reference_Namelist variables that must be set for this model are
indicated in Table table_boussinesq_.

Note that our choice for the temperature scale assumes fixed-temperature
boundary conditions. We might choose to specify fixed-flux boundary
conditions and/or an internal heating, in which case the meaning of
:math:`\Delta T` in our equation set changes, with
:math:`\Delta T \equiv L\frac{\partial T}{\partial r}` instead, for some
fiducial value of :math:`\frac{\partial T}{\partial r}`. Which regard to
the temperature scaling, it is up to the user to select boundary
conditions appropriate for their desired values of :math:`\Delta T`. If
:math:`\Delta T` denotes the temperature contrast across the domain,
then their boundary condition variables should look like:

::

   &boundary\_conditions\_namelist
   T_Top    = 0.0d0
   T_Bottom = 1.0d0
   fix_tvar_top = .true.
   fix_tvar_bottom = .true.
   /

Alternatively, if the temperature scale is determined by a gradient at
one boundary, the user should ensure that the amplitude of the
temperature gradient at that boundary is 1. For example:

::

   &boundary\_conditions\_namelist
   dTdr_bottom = -1.0d0
   fix_dtdr_bottom = .true.
   /

Boundary conditions and internal heating are discussed in
§\ :ref:`boundary_conditions`. Finally, in Boussinesq mode, the
namelist variables **nu_type**, **kappa_type**, and **eta_type** should
be set to 1. Their values will be determined by Pr and Pm, instead of
nu_top, kappa_top, or eta_top.

   .. _table_boussinesq:

.. centered:: **Table. Boussinesq.**

Variables in the Reference_Namelist that
must be specified when running in nondimensional Boussinesq mode. In
addition, reference_type=1 must also be specified.

   +-----------------------------------+-----------------------------------+
   | Variable                          | Description                       |
   +===================================+===================================+
   | Ekman_Number                      | The Ekman Number :math:`E`        |
   +-----------------------------------+-----------------------------------+
   | Rayleigh_Number                   | The Rayleigh Number :math:`Ra`    |
   +-----------------------------------+-----------------------------------+
   | Prandtl_Number                    | The Prandtl Number :math:`Pr`     |
   +-----------------------------------+-----------------------------------+
   | Magnetic_Prandtl_Number           | The Magnetic Prandtl Number       |
   |                                   | :math:`Pm`                        |
   +-----------------------------------+-----------------------------------+
   | Gravity_Power                     | Buoyancy coefficient =            |
   |                                   | :math:`\frac{\mathrm{Ra}}{\mathrm |
   |                                   | {Pr}}\left(\frac{r}{rmax} \right) |
   |                                   | ^\mathrm{gravity\_power}`         |
   +-----------------------------------+-----------------------------------+

.. raw:: latex

   \clearpage

Anelastic Mode (nondimensional)
-------------------------------

**Example Input: Rayleigh/input_examples/main_input_jupiter**

When running in nondimensional anelastic mode, you must set
**reference_type=3** in the Reference_Namelist. When this parameter is
set, the following nondimensionalization is used (following :cite:`Heimpel_etal2016`):

.. math::

   \begin{aligned}
   \mathrm{Length} &\rightarrow L &\;\;\;\; \mathrm{(Shell\; Depth)} \\
   \mathrm{Time} &\rightarrow   \frac{1}{\Omega_0} &\;\;\;\; \mathrm{(Rotational\; Timescale)}\\
   \mathrm{Temperature} &\rightarrow T_o\equiv\hat{T}(r_\mathrm{max})&\;\;\;\; \mathrm{(Reference-State\; Temperature\; at\; Upper\; Boundary)} \\
   \mathrm{Density} &\rightarrow \rho_o\equiv\hat{\rho}(r_\mathrm{max})&\;\;\;\; \mathrm{(Reference-State\; Density\; at\; Upper\; Boundary)} \\
   \mathrm{Entropy} &\rightarrow \Delta{s}&\;\;\;\; \mathrm{(Entropy\; Constrast\; Across\; Shell)} \\
   \mathrm{Magnetic~Field} &\rightarrow \sqrt{\rho_o\mu\eta\Omega_0}.\end{aligned}

We assume a polytropic background state (similar to dimensional
anelastic mode), with gravity varying as :math:`\frac{1}{r^2}`. We
further assume that the transport coefficients :math:`\nu`,
:math:`\kappa`, and :math:`\eta` do not vary with radius. The results in
the nondimensionalized equations (tildes used to indicated
nondimensional reference-state values):

.. math::

   \begin{aligned}
   \frac{\partial \boldsymbol{v}}{\partial t} +\boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}  %advection
                                                            +2\boldsymbol{\hat{z}}\times\boldsymbol{v}  &= % Coriolis
                                                            \mathrm{Ra}^*\frac{r_\mathrm{max}^2}{r^2}\Theta\,\boldsymbol{\hat{r}} % buoyancy
                                                            +\boldsymbol{\nabla}\left(\frac{P}{\tilde{\rho}(r)}\right) % pressure
                                                            +\frac{\mathrm{E}}{\mathrm{Pm}\,\tilde{\rho}}\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
                                                            +\frac{\mathrm{E}}{\tilde{\rho(r)}}\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}} \;\;\; &\mathrm{Momentum}\\
   %
   %
   \tilde{\rho}(r)\,\tilde{T}(r)\left[\frac{\partial \Theta}{\partial t} +\boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta \right] &=
                                                \frac{\mathrm{E}}{\mathrm{Pr}}\boldsymbol{\nabla}\cdot\left[\tilde{\rho}(r)\,\tilde{T}(r)\,\boldsymbol{\nabla}\Theta \right] % diffusion
                                                +Q(r)   % Internal heating
                                                +\frac{\mathrm{E}\,\mathrm{Di}}{\mathrm{Ra}^*}\Phi(r,\theta,\phi)
                                                +\frac{\mathrm{Di\,E^2}}{\mathrm{Pm}^2\mathrm{R}^*}\left[\boldsymbol{\nabla}\times\boldsymbol{B}\right]^2 &\mathrm{Thermal\; Energy}\\ % Ohmic Heating
   %
   %
   \frac{\partial \boldsymbol{B}}{\partial t} &= \boldsymbol{\nabla}\times\left(\,\boldsymbol{v}\times\boldsymbol{B}-\frac{\mathrm{E}}{\mathrm{Pm}}\boldsymbol{\nabla}\times\boldsymbol{B}\,\right) &\mathrm{Induction} \\
   %
   %
   \mathcal{D}_{ij} &= 2\tilde{\rho}(r)\left[e_{ij}-\frac{1}{3}\boldsymbol{\nabla}\cdot\boldsymbol{v}\right] &\mathrm{Viscous\; Stress\; Tensor}\\
   %
   %
   \Phi(r,\theta,\phi) &= 2\,\tilde{\rho}(r)\left[e_{ij}e_{ij}-\frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)^2\right] &\mathrm{Viscous\; Heating} \\
   %
   %
   \boldsymbol{\nabla}\cdot\left(\tilde{\rho}(r)\,\boldsymbol{v}\right)&=0 &\mathrm{Solenoidal\; Mass\; Flux}\\
   \boldsymbol{\nabla}\cdot\boldsymbol{B}&=0. &\mathrm{Solenoidal\; Magnetic\; Field}\end{aligned}

In the equations above, Di is the dissipation number, defined by

.. math:: \mathrm{Di}= \frac{g_o\,\mathrm{L}}{c_\mathrm{P}\,T_o},

where :math:`g_o` and :math:`T_o` are the gravitational acceleration
and temperature at the outer boundary respectively. Once more, the
thermal anomoly :math:`\Theta` should be interpreted as entropy
:math:`s`. The symbol Ra\ :math:`^*` is the modified Rayleigh number,
given by

.. math:: \mathrm{Ra}^*=\frac{g_o}{c_\mathrm{P}\Omega_0^2}\frac{\Delta s}{L}

Those Reference_Namelist variables that must be set for this model are
indicated in Table table_anelastic_nd_. As
with :math:`\Delta T` in the nondimensional Boussinesq mode, the user
must choose boundary conditions suitable for their definition of
:math:`\Delta s`. As with the dimensional anelastic formulation, the
background state is polytropic and is described through a polytropic
index and number of density scale heights.

**Note:** As with the Boussinesq mode, please set the variables
**nu_type**, **kappa_type**, **eta_type** in the Transport_Namelist.

   .. _table_anelastic_nd:

.. centered:: **Table. Anelastic_nd.**

Variables in the Reference_Namelist that
must be specified when running in nondimensional anelastic mode. In
addition, reference_type=3 must also be specified.

   +-----------------------------------+-----------------------------------+
   | Variable                          | Description                       |
   +===================================+===================================+
   | Ekman_Number                      | The Ekman Number E                |
   +-----------------------------------+-----------------------------------+
   | Modified_Rayleigh_Number          | The Modified Rayleigh Number      |
   |                                   | Ra\ :math:`^*`                    |
   +-----------------------------------+-----------------------------------+
   | Prandtl_Number                    | The Prandtl Number Pr             |
   +-----------------------------------+-----------------------------------+
   | Magnetic_Prandtl_Number           | The Magnetic Prandtl Number Pm    |
   +-----------------------------------+-----------------------------------+
   | poly_n [0]                        | polytropic index                  |
   |                                   | (:math:`P\propto\rho^n`)          |
   +-----------------------------------+-----------------------------------+
   | poly_Nrho [0]                     | number of density scaleheights    |
   |                                   | spanning the domain               |
   +-----------------------------------+-----------------------------------+

.. _boundary_conditions:

Boundary Conditions & Internal Heating
--------------------------------------

Boundary conditions are controlled through the
**Boundary_Conditions_Namelist**. All Rayleigh simulations are run with
impenetrable boundaries. These boundaries may be either no-slip or
stress-free (default). If you want to employ no-slip conditions at both
boundaries, set **no_slip_boundaries = .true.**. If you wish to set
no-slip conditions at only one boundary, set **no_slip_top=.true.** or
**no_slip_bottom=.true.** in the Boundary_Conditions_Namelist.

By default, magnetic boundary conditions are set to match to a potential field at
each boundary.

By default, the thermal anomoly :math:`\Theta` is set to a fixed value
at each boundary. The upper and lower boundary-values are specified by
setting **T_top** and **T_bottom** respectively in the
Boundary_Conditions_Namelist. Their defaults values are 1 and 0
respectively.

Alternatively, you may specify a constant value of
:math:`\partial\Theta/\partial r` at each boundary. This is accomplished
by setting the variables **fix_dTdr_top** and **fix_dTdr_bottom**.
Values of the gradient may be enforced by setting the namelist variables
**dTdr_top** and **dTdr_bottom**. Both default to a value of zero.

Generic Boundary Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Boundary conditions for temperature, :math:`T`, and the magnetic poloidal potential, :math:`C`,
may also be set using generic spectral input.  As with initial conditions (see :ref:`sec:generic_ic`)
this is done by generating a generic input file using the **rayleigh_spectral_input.py** script.
The file output from this script can then be applied using:

-  **fix_Tvar_top** and **T_top_file** (overrides **T_top** for a constant boundary condition)
   to set a fixed upper :math:`T` boundary condition
-  **fix_dTdr_top** and **dTdr_top_file** (overrides **dTdr_top**) to set a fixed upper :math:`T` gradient boundary condition
-  **fix_Tvar_bottom** and **T_bottom_file** (overrides **T_bottom**) to set a fixed lower :math:`T` boundary condition
-  **fix_dTdr_bottom** and **dTdr_bottom_file** (overrides **dTdr_bottom**) to set a fixed lower :math:`T` gradient boundary
   condition
-  **fix_poloidal_top** and **C_top_file** (overrides **impose_dipole_field**) to set a fixed upper :math:`C` boundary condition
-  **fix_poloidal_bottom** and **C_bottom_file** (overrides **impose_dipole_field**) to set a fixed lower :math:`C` boundary condition

For example, to set a :math:`C` boundary condition on both boundaries with modes (l,m) = (1,0) and (1,1) set to pre-calculated
values run:

::
   
   rayleigh_spectral_input.py -m 1 0 2.973662220170157 -m 1 1 0.5243368809294343+0.j -o ctop_init_bc
   rayleigh_spectral_input.py -m 1 0 8.496177771914736 -m 1 1 1.4981053740840984+0.j -o cbottom_init_bc

which will generate generic spectral input files `ctop_init_bc` and `cbottom_init_bc`.  Set these to be used as the boundary
conditions in `main_input` using:

::

   &Boundary_Conditions_Namelist
   fix_poloidalfield_top = .true.
   fix_poloidalfield_bottom = .true.
   C_top_file = 'ctop_init_bc'
   C_bottom_file = 'cbottom_init_bc'
   /

This can be seen being applied in `tests/generic_input`.

Internal Heating
~~~~~~~~~~~~~~~~

The internal heating function :math:`Q(r)` is activated and described by
two variables in the **Reference_Namelist**. These are **Luminosity**
and **heating_type**. Note that these values are part of the
**Reference_Namelist** and not the **Boundary_Conditions** namelist.
Three heating types (0,1, and 4) are fully supported at this time.
Heating type zero corresponds to no heating. This is the default.

**Heating_type=1:** This heating type is given by :

.. math::

   \label{eq:heating}
   %\frac{\partial \Theta}{\partial t}=\gamma\left( 1 -\frac{\hat{\rho}(r_\mathrm{max})\,\hat{T}(r_\mathrm{max})  }{\hat{\rho}(r)\, \hat{T}(r)} \right),
   Q(r)= \gamma\,\hat{\rho}(r)\, \hat{T}(r)

where :math:`\gamma` is a normalization constant defined such that

.. _eq_lum:

.. math::


   %4\pi r_o^2 \hat{\rho}\hat{T}\kappa(r)\frac{\partial \Theta}{\partial r}=\mathrm{Luminosity}
   4\pi\int_{r=r_\mathrm{min}}^{r=r_\mathrm{max}} Q(r)\,  r^2 dr = \mathrm{Luminosity}.

This heating profile is particularly useful for emulating radiative
heating in a stellar convection zone.

**Heating_type=4:** This heating type corresponds a heating that is
variable in radius, but constant in *energy density*. Namely

.. math:: \hat{\rho}\hat{T}\frac{\partial \Theta}{\partial t}=\gamma.

The constant :math:`\gamma` in this case is also set by enforcing
Equation eq_lum_.

General Physics Controls
------------------------

A number of logical variables can be used to turn certain physics on
(value = .true.) or off ( value = .false.). These variables are
described in Table table_logicals_, with default
values indicated in brackets.

  .. _table_logicals:

.. centered:: **Table. Logicals.**

Variables in the Physical_Controls_Namelist
that may be specified to control run behavior (defaults indicated in
brackets)

   +-----------------------------------+-----------------------------------+
   | Variable                          | Description                       |
   +===================================+===================================+
   | magnetism [.false.]               | Turn magnetism on or off          |
   +-----------------------------------+-----------------------------------+
   | rotation [.false.]                | Turn rotation on or off (pressure |
   |                                   | is not scaled by E when off)      |
   +-----------------------------------+-----------------------------------+
   | lorentz_forces [.true.]           | Turn Lorentz forces on or off     |
   |                                   | (magnetism must be .true.)        |
   +-----------------------------------+-----------------------------------+
   | viscous_heating [.true.]          | Turn viscous heating on or off    |
   |                                   | (inactive in Boussinesq mode)     |
   +-----------------------------------+-----------------------------------+
   | ohmic_heating [.true.]            | Turn ohmic heating off or on      |
   |                                   | (inactive in Boussinesq mode)     |
   +-----------------------------------+-----------------------------------+

Initializing a Model
--------------------

A Rayleigh simulation may be initialized with a random thermal and/or
magnetic field, or it may be restarted from an existing checkpoint file
(see §\ :ref:`checkpointing` for a detailed
discussion of checkpointing). This behavior is controlled through the
**initial_conditions_namelist** and the **init_type** and
**magnetic_init_type** variables. The init_type variable controls the
behavior of the velocity and thermal fields at initialization time.
Available options are:

-  init_type=-1 ; read velocity and thermal fields from a checkpoint
   file

-  init_type=1 ; Christensen et al. (2001) case 0 benchmark init (
   {:math:`\ell=4,m=4`} temperature mode)

-  init_type=6 ; Jones et al. (2011) steady anelastic benchmark (
   {:math:`\ell=19,m=19`} entropy mode)

-  init_type=7 ; random temperature or entropy perturbation

-  init_type=8 ; user generated temperature or entropy perturbation
   (see Generic Initial Conditions below)

When initializing a random thermal field, all spherical harmonic modes
are independently initialized with a random amplitude whose maximum
possible value is determined by the namelist variable **temp_amp**. The
mathematical form of of this random initialization is given by

.. _eq_init:

.. math::

   T(r,\theta,\phi) = \sum_\ell \sum_m  c_\ell^m f(r)g(\ell)\mathrm{Y}_\ell^m(\theta,\phi),

where the :math:`c_\ell^m`\ ’s are (complex) random amplitudes,
distributed normally within the range [-temp_amp, temp_amp]. The radial
amplitude :math:`f(r)` is designed to taper off to zero at the
boundaries and is given by

.. math:: f(r) = \frac{1}{2}\left[1-\mathrm{cos}\left( 2\pi\frac{r-rmin}{rmax-rmin} \right)   \right].

The amplitude function :math:`g(\ell)` concentrates power in the
central band of spherical harmonic modes used in the simulation. It is
given by

.. math:: g(\ell) = \mathrm{exp}\left[  - 9\left( \frac{ 2\,\ell-\ell_\mathrm{max} }{ \ell_\mathrm{max} }  \right)^2 \right],

which is itself, admittedly, a bit random.

When initializing using a random thermal perturbation, it is important
to consider whether it makes sense to separately initialize the
spherically-symmetric component of the thermal field with a profile that
is in conductive balance. This is almost certainly the case when running
with fixed temperature conditions. The logical namelist variable
**conductive_profile** can be used for this purpose. It’s default value
is .false. (off), and its value is ignored completely when restarting
from a checkpoint. To initialize a simulation with a random temperature
field superimposed on a spherically-symmetric, conductive background
state, something similar to the following should appear in your
main_input file:

::

   &initial_conditions_namelist
   init_type=7
   temp_amp = 1.0d-4
   conductive_profile=.true.
   /

Magnetic-field initialization follows a similar pattern. Available
values for magnetic_input type are:

-  magnetic_init_type = -1 ; read magnetic field from a checkpoint file

-  magnetic_init_type = 1 ; Christensen et al. (2001) case 0 benchmark
   init

-  magnetic_init_type = 7 ; randomized vector potential

-  magnetic_init_type=8 ; user generated magnetic potential fields
   (see Generic Initial Conditions below)

For the randomized magnetic field, both the poloidal and toroidal
vector-potential functions are given a random power distribution
described by Equation eq_init_. Each mode’s random
amplitude is then determined by namelist variable **mag_amp**. This
variable should be interpreted as an approximate magnetic field strength
(it’s value is rescaled appropriately for the poloidal and toroidal
vector potentials, which are differentiated to yield the magnetic
field).

When initializing all fields from scratch, a main_input file should
contain something similar to:

::

   &initial_conditions_namelist
   init_type=7
   temp_amp = 1.0d-4
   conductive_profile=.true.  ! Not always necessary (problem dependent) ...
   magnetic_init_type=7
   mag_amp = 1.0d-1
   /


.. _sec:generic_ic:

Generic Initial Conditions
~~~~~~~~~~~~~~~~~~~~~~~~~~

The user can input any initial conditions from data files generated by
a python routine "rayleigh_spectral_input.py", which can be called as
a script or imported as a python class.

The available generic initial conditions options are

::

   &initial_conditions_namelist
   init_type=8
   T_init_file = '<filename>'  !! Temperature
   W_init_file = '<filename>'  !! Poloidal velocity potential
   Z_init_file = '<filename>'  !! Toroidal velocity potential
   P_init_file = '<filename>'  !! `Pressure` potential
   magneic_init_type=8
   C_init_file = '<filename>'  !! Poloidal magnetic potential
   A_init_file = '<filename>'  !! Toroidal magnetic potential

   /

where `T_init_file` is a user generated initial temperature field and
<filename> is the name of the file generated by the python script.  If
`T_init_file` is not specified the initial field will be zero by
default.  The same for the other fields.  Fields T, W, Z, and P are
only initialized from the file if `init_type=8`.  Fields C and A are
only initialized from file if `magnetic_init_type=8`.

To generate a generic initial condition input file, for example, if a user wanted to specify a single mode in that input file then they could just run the script:

::

   rayleigh_spectral_input.py -m 0 0 0 1.+0.j -o example


to specify (n,l,m) = (0,0,0) to have a coefficient 1.+0.j and output it to the file example.

This could also be done using the python as a module. In a python
shell this would look like:

::

   from rayleigh_spectral_input import *
   si = SpectralInput()
   si.add_mode(1., n=0, l=0, m=0)
   si.write('example')


For a more complicated example, e.g. the hydrodynamic benchmark from
Christensen et al. 2001, the user can specify functions of theta, phi
and radius that the python will convert to spectral:

::

   rayleigh_spectral_input.py -ar 0.35 -sd 1.0 -nt 96 -nr 64 -o example \
    -e 'import numpy as np; x = 2*radius - rmin - rmax;
    rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) -
    x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

in "script" mode.

Alternatively, in "module" mode in a python shell:

::

   from rayleigh_spectral_input import *
   si = SpectralInput(n_theta=96, n_r=64)
   rmin, rmax = radial_extents(aspect_ratio=0.35, shell_depth=1.0)
   def func(theta, phi, radius):
      x = 2*radius - rmin - rmax
      return rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) - x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)
   si.transform_from_rtp_function(func, aspect_ratio=0.35, shell_depth=1.0)
   si.write('example')


The above commands will generate a file called `example` which can be
called by

::

   &initial_conditions_namelist
   init_type=8
   T_init_file = 'example'

Note that these two examples will have produced different data formats - the first one sparse (listing only the mode specified) and the second one dense (listing all modes).

For more examples including magnetic potentials see `tests/generic_input`.
