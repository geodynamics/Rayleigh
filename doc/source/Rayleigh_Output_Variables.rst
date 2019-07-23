.. raw:: latex

   \maketitle

.. raw:: latex

   \tableofcontents

.. raw:: latex

   \clearpage

Overview of Diagnostic Outputs in Rayleigh
==========================================

The purpose of this document is to describe Rayleigh’s internal menu
system used for specifying diagnostic outputs. Rayleigh’s design
includes an onboard diagnostics package that allows a user to output a
variety of system quantities as the run evolves. These include system
state variables, such as velocity and entropy, as well as derived
quantities, such as the vector components of the Lorentz force and the
kinetic energy density. Each diagnostic quantity is requested by adding
its associated menu number to the *main_input* file. Radial velocity,
for instance, has menu code 1, :math:`\theta`-component of velocity has
menu code 2, etc.

A few points to keep in mind are

-  **This document is intended to describe the diagnostics output menu
   only.** A complete description of Rayleigh’s diagnostic package is
   provided in Rayleigh/doc/Diagnostic_Plotting.pdf. A more in-depth
   description of the anelastic and Boussinesq modes available in
   Rayleigh is provided in Rayleigh/doc/user_guide.pdf.

-  A number of *output methods* may be used to output any system
   diagnostic. No diagnostic is linked to a particular *output method*.
   The same diagnostic might be output in volume-averaged,
   azimuthally-averaged, and fully 3-D form, for instance.

-  You may notice a good deal of redundancy in the available outputs.
   For instance, the azimuthal velocity, :math:`v_\phi`, and its zonal
   average, :math:`\overline{v_\phi}`, are both available as outputs.
   Were the user to output both of these in an azimuthally-averaged
   format, the result would be the same. 3-D output, however, would not
   yield the same result. This redundancy has been added to help with
   post-processing calculations in which it can be useful to have all
   data products in a similar format.

-  Given the degree of redundancy found in the list below, you may be
   surprised to notice that several values are not available for output
   at all. Some of these are best added as custom-user diagnostics and
   may be included in a future release. Many, however, may be obtained
   by considering either the sum, or difference, of those outputs
   already available.

Definitions and Conventions
===========================

Vector and Tensor Notation
--------------------------

All vector quantities are represented in bold italics. Components of a
vector are indicated in non-bold italics, along with a subscript
indicating the direction associated with that component. Unit vectors
are written in lower-case, bold math font and are indicated by the use
of a *hat* character. For example, a vector quantity
:math:`\boldsymbol{a}` would represented as

.. math::
   :label: vcomp

       \boldsymbol{a} = a_r\boldsymbol{\hat{a}}+a_\theta\boldsymbol{\hat{\theta}}+a_\phi\boldsymbol{\hat{\phi}}.

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
----------------------

The *hat* notation is also used to indicate reference-state quantities.
These quantities are scalar, and they are not written in bold font. They
vary only in radius and have no :math:`\theta`-dependence or
:math:`\phi`-dependence. The reference-state density is indicated by
:math:`\hat{\rho}` and the reference-state temperature by
:math:`\hat{T}`, for instance.

Averaged and Fluctuating Values
-------------------------------

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

The Equation Sets Solved by Rayleigh
====================================

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
        + c_1\boldsymbol{\hat{z}}\times\boldsymbol{v} \right]  = % Coriolis
       c_2\,\mathrm{f}_2(r)\Theta\,\boldsymbol{\hat{r}} % buoyancy
        - c_3\,\mathrm{f}_1(r)\boldsymbol{\nabla}\left(\frac{P}{\mathrm{f}_1(r)}\right) % pressure
        + c_4\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
        + c_5\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}},

where the stress tensor :math:`\mathcal{D}` is given by

.. math:: 
   :label: stress_tensor

       \mathcal{D}_{ij} = 2\mathrm{f}_1(r)\,\mathrm{f}_3(r)\left[e_{ij} - \frac{1}{3}\boldsymbol{\nabla}\cdot\boldsymbol{v}\right].

The velocity field is denoted by :math:`\boldsymbol{v}`, the thermal
anomoly by :math:`\Theta`, the pressure by :math:`P`, and the magnetic
field by :math:`\boldsymbol{B}`. All four of these quantities are
3-dimensional functions of position, in contrast to the 1-dimensional
coefficient functions :math:`\mathrm{f}_i`. The velocity and magnetic
fields are subject to the constraints

.. math::
   :label: v_constrain
 
       \boldsymbol{\nabla}\cdot\left(\mathrm{f}_1(r)\,\boldsymbol{v}\right) = 0

and

.. math::
   :label: divB
 
       \boldsymbol{\nabla}\cdot\boldsymbol{B}=0

respectively. The evolution of :math:`\Theta` is described

.. math::
   :label: theta_evol

   \mathrm{f}_1(r)\,\mathrm{f}_4(r)\left[\frac{\partial \Theta}{\partial t}  + \boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta \right] =
       c_6\,\boldsymbol{\nabla}\cdot\left[\mathrm{f}_1(r)\,\mathrm{f}_4(r)\,\mathrm{f}_5(r)\,\boldsymbol{\nabla}\Theta \right] % diffusion
        + c_{10}\mathrm{f}_6(r)   % Internal heating
        + c_8\Phi(r,\theta,\phi)
        + c_9\mathrm{f}_7(r)\left[\boldsymbol{\nabla}\times\boldsymbol{B}\right]^2,  % Ohmic Heating

where the viscous heating :math:`\Phi` is given by

.. math::
   :label: vischeat

       \Phi(r,\theta,\phi) = 2\,\mathrm{f}_1(r)\mathrm{f}_3(r)\left[e_{ij}e_{ij} - \frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)^2\right].

Finally, the evolution of :math:`\boldsymbol{B}` is described by the
induction equation

.. math::
   :label: induction

       \frac{\partial \boldsymbol{B}}{\partial t} = \boldsymbol{\nabla}\times\left(\,\boldsymbol{v}\times\boldsymbol{B} - c_7\,\mathrm{f}_7(r)\boldsymbol{\nabla}\times\boldsymbol{B}\,\right).

Equations :eq:`momentum`-:eq:`induction` could have been formulated in other ways. For instance, we could combine
:math:`\mathrm{f}_1` and :math:`\mathrm{f}_3` into a single function in
Equation :eq:`vischeat`. The form of the equations
presented here has been chosen to reflect that actually used in the
code, which was originally written dimensionally. We now describe the
dimensional anelastic and nondimensional Boussinesq formulations used in
the code.

Dimensional Anelastic Formulation of the MHD Equations
------------------------------------------------------

When run in dimensional, anelastic mode (cgs units; **reference_type=2**
), the following values are assigned to the functions
:math:`\mathrm{f}_i` and the constants :math:`c_i`:

.. math::

   \begin{aligned}
       \mathrm{f}_1(r) &\rightarrow \hat{\rho}(r)\; &c_1 &\rightarrow 2\Omega_0 \\
       \mathrm{f}_2(r) &\rightarrow \frac{\hat{\rho(r)}}{c_P}g(r)\; &c_2 &\rightarrow 1 \\
       \mathrm{f}_3(r) &\rightarrow \nu(r)\; &c_3 &\rightarrow 1\\
       \mathrm{f}_4(r) &\rightarrow \hat{T}(r)\; &c_4 &\rightarrow \frac{1}{4\pi} \\
       \mathrm{f}_5(r) &\rightarrow \kappa(r)\; &c_5 &\rightarrow 1 \\
       \mathrm{f}_6(r) &\rightarrow Q(r)\; &c_6 &\rightarrow 1  \\
       \mathrm{f}_7(r) &\rightarrow \eta(r)\; &c_7 &\rightarrow 1 \\
       c_8&\rightarrow 1 &c_9 &\rightarrow \frac{1}{4\pi} \\
       c_{10}&\rightarrow 1.\end{aligned}

Here, :math:`\hat{\rho}` and :math:`\hat{T}` are the reference-state
density and temperature respectively. :math:`g` is the gravitational
acceleration, :math:`c_P` is the specific heat at constant pressure, and
:math:`\Omega_0` is the frame rotation rate. The viscous, thermal, and
magnetic diffusivities are given by :math:`\nu`, :math:`\kappa`, and
:math:`\eta` respectively. Finally, :math:`Q(r)` is an internal heating
function; it might represent radiative heating or heating due to nuclear
fusion, for instance. Note that in the anelastic formulation, the
thermal variable :math:`\Theta` is interpreted is as entropy :math:`s`,
rather than temperature :math:`T`. When these substitutions are made,
Equations :eq:`momentum`-:eq:`induction` transform as follows.

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

Nondimensional Boussinesq Formulation of the MHD Equations
----------------------------------------------------------

Rayleigh can also be run using a nondimensional, Boussinesq formulation
of the MHD equations (**reference_type=1**). The nondimensionalization
employed is as follows:

.. math::

   \begin{aligned}
       \mathrm{Length} &\rightarrow L &\;\;\;\; \mathrm{(Shell\; Depth)} \\
       \mathrm{Time} &\rightarrow   \frac{L^2}{\nu} &\;\;\;\; \mathrm{(Viscous\; Timescale)}\\
       \mathrm{Temperature} &\rightarrow \Delta T&\;\;\;\; \mathrm{(Temperature\; Contrast\; Across\; Shell)} \\
       \mathrm{Magnetic\; Field} &\rightarrow \sqrt{\rho\mu\eta\Omega_0},\end{aligned}

where :math:`\Omega_0` is the rotation rate of the frame, :math:`\rho`
is the (constant) density of the fluid, :math:`\mu` is the magnetic
permeability, :math:`\eta` is the magnetic diffusivity, and :math:`\nu`
is the kinematic viscosity. After nondimensionalizing, the following
nondimensional numbers appear in our equations:

.. math::

   \begin{aligned}
       Pr &=\frac{\nu}{\kappa}                          &\;\;\;\;\;\; \mathrm{Prandtl\; Number} \\
       Pm &=\frac{\nu}{\eta}                            &\;\;\;\;\;\; \mathrm{Magnetic\; Prandtl\; Number} \\
       E  &=\frac{\nu}{\Omega_0\,L^2}                   &\;\;\;\;\;\; \mathrm{Ekman\; Number} \\
       Ra &=\frac{\alpha g_0 \Delta T\,L^3}{\nu\kappa}  &\;\;\;\;\;\; \mathrm{Rayleigh\; Number},\end{aligned}

where :math:`\alpha` is coefficient of thermal expansion, :math:`g_0`
is the gravitational acceleration, and :math:`\kappa` is the thermal
diffusivity. Adopting this nondimensionalization is equivalent to
assigning values to :math:`\mathrm{f}_i` and the constants :math:`c_i`:

.. math::

   \begin{aligned}
   \mathrm{f}_1(r) &\rightarrow 1\; &c_1 &\rightarrow \frac{2}{E} \\
   \mathrm{f}_2(r) &\rightarrow \left(\frac{r}{r_o}\right)^n \; &c_2 &\rightarrow \frac{Ra}{E\,Pr} \\
   \mathrm{f}_3(r) &\rightarrow 1\; &c_3 &\rightarrow \frac{1}{E}\\
   \mathrm{f}_4(r) &\rightarrow 1\; &c_4 &\rightarrow \frac{1}{E\,Pm} \\
   \mathrm{f}_5(r) &\rightarrow 1\; &c_5 &\rightarrow 0 \\
   \mathrm{f}_6(r) &\rightarrow 0\; &c_6 &\rightarrow \frac{1}{Pr}  \\
   \mathrm{f}_7(r) &\rightarrow 1\; &c_7 &\rightarrow \frac{1}{Pm} \\
   c_8&\rightarrow 0 &c_9 &\rightarrow 0 \\ 
   c_{10}&\rightarrow 0.\end{aligned}

Note that our choice of :math:`\mathrm{f}_2(r)` allows gravity to vary
with radius based on the value of the exponent :math:`n`, which has a
default value of 0 in the code. Note also that our definition of
:math:`Ra` assumes fixed-temperature boundary conditions. We might
choose specify fixed-flux boundary conditions and/or an internal heating
through a suitable choice :math:`\mathrm{f}_6(r)`, in which case the
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
        + \boldsymbol{\nabla}^2\boldsymbol{v} \;\;\; &\mathrm{Momentum}\\
       %
       %
       \left[\frac{\partial \Theta}{\partial t}  + \boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta \right] &=
       \frac{1}{Pr}\boldsymbol{\nabla}^2\Theta  &\mathrm{Thermal\; Energy}\\ % Diffusion
       %
       %
       \frac{\partial \boldsymbol{B}}{\partial t} &= \boldsymbol{\nabla}\times\left(\,\boldsymbol{v}\times\boldsymbol{B}\right)+\frac{1}{Pm}\boldsymbol{\nabla}^2\boldsymbol{B} &\mathrm{Induction} \\
       %
       %
       \boldsymbol{\nabla}\cdot\boldsymbol{v}&=0 &\mathrm{Solenoidal\; Velocity\; Field}\\
       \boldsymbol{\nabla}\cdot\boldsymbol{B}&=0 &\mathrm{Solenoidal\; Magnetic\; Field}\end{aligned}

Nondimensional Anelastic MHD Equations
--------------------------------------

To run in nondimensional anelastic mode, you must set
**reference_type=3** in the Reference_Namelist. The reference state is
assumed to be polytropic with a :math:`\frac{1}{r^2}` profile for
gravity. Transport coefficients :math:`\nu`, :math:`\kappa`,
:math:`\eta` are assumed to be constant in radius. When this mode is
active, the following nondimensionalization is used 
(following `Heimpel et al., 2016, Nat. Geo., 9, 19 <https://www.nature.com/articles/ngeo2601/>`_ ):

.. math::

   \begin{aligned}
       \mathrm{Length} &\rightarrow L &\;\;\;\; \mathrm{(Shell\; Depth)} \\
       \mathrm{Time} &\rightarrow   \frac{1}{\Omega_0} &\;\;\;\; \mathrm{(Rotational\; Timescale)}\\
       \mathrm{Temperature} &\rightarrow T_o\equiv\hat{T}(r_\mathrm{max})&\;\;\;\; \mathrm{(Reference-State\; Temperature\; at\; Upper\; Boundary)} \\
       \mathrm{Density} &\rightarrow \rho_o\equiv\hat{\rho}(r_\mathrm{max})&\;\;\;\; \mathrm{(Reference-State\; Density\; at\; Upper\; Boundary)} \\
       \mathrm{Entropy} &\rightarrow \Delta{s}&\;\;\;\; \mathrm{(Entropy\; Constrast\; Across\; Shell)} \\
       \mathrm{Magnetic~Field} &\rightarrow \sqrt{\tilde{\rho}(r_\mathrm{max})\mu\eta\Omega_0}.\end{aligned}

When run in this mode, Rayleigh employs a polytropic background state,
with an assumed :math:`\frac{1}{r^2}` variation in gravity. These
choices result in the functions :math:`\mathrm{f}_i` and the constants
:math:`c_i` (tildes indicate nondimensional reference-state variables):

.. math::

   \begin{aligned}
       \mathrm{f}_1(r) &\rightarrow \tilde{\rho}(r)\; &c_1 &\rightarrow 2 \\
       \mathrm{f}_2(r) &\rightarrow \tilde{\rho(r)}\frac{r_\mathrm{max}^2}{r^2}\; &c_2 &\rightarrow \mathrm{Ra}^* \\
       \mathrm{f}_3(r) &\rightarrow 1\; &c_3 &\rightarrow 1\\
       \mathrm{f}_4(r) &\rightarrow \tilde{T}(r)\; &c_4 &\rightarrow \frac{\mathrm{E}}{\mathrm{Pm}} \\
       \mathrm{f}_5(r) &\rightarrow 1\; &c_5 &\rightarrow \mathrm{E} \\
       \mathrm{f}_6(r) &\rightarrow Q(r)\; &c_6 &\rightarrow \frac{\mathrm{E}}{\mathrm{Pr}}  \\
       \mathrm{f}_7(r) &\rightarrow 1 \; &c_7 &\rightarrow \frac{\mathrm{E}}{\mathrm{Pm}} \\
       c_8&\rightarrow \frac{\mathrm{E}\,\mathrm{Di}}{\mathrm{Ra}^*} &c_9 &\rightarrow \frac{\mathrm{E}^2\,\mathrm{Di}}{\mathrm{Pm}^2\mathrm{Ra}^*} \\
       c_{10}&\rightarrow 1.\end{aligned}

Two new nondimensional numbers appear in our equations. Di, the
dissipation number, is defined by

.. math::
   :label: Di
 
       \mathrm{Di}= \frac{g_o\,\mathrm{L}}{c_\mathrm{P}\,T_o},

where :math:`g_o` and :math:`T_o` are the gravitational acceleration
and temperature at the outer boundary respectively. Once more, the
thermal anomaly :math:`\Theta` should be interpreted as entropy
:math:`s`. The symbol Ra\ :math:`^*` is the modified Rayleigh number,
given by

.. math:: 
   :label: Ra

   \mathrm{Ra}^*=\frac{g_o}{c_\mathrm{P}\Omega_0^2}\frac{\Delta s}{L}   %\frac{\partial \Theta}{\partial r}|_{r=rmin}

We arrive at the following nondimensionalized equations:

.. math::

   \begin{aligned}
       \frac{\partial \boldsymbol{v}}{\partial t}  + \boldsymbol{v}\cdot\boldsymbol{\nabla}\boldsymbol{v}  %advection 
        + 2\boldsymbol{\hat{z}}\times\boldsymbol{v}  &= % Coriolis
       \mathrm{Ra}^*\frac{r_\mathrm{max}^2}{r^2}\Theta\,\boldsymbol{\hat{r}} % buoyancy
        + \boldsymbol{\nabla}\left(\frac{P}{\tilde{\rho}(r)}\right) % pressure
        + \frac{\mathrm{E}}{\mathrm{Pm}\,\tilde{\rho}}\left(\boldsymbol{\nabla}\times\boldsymbol{B}\right)\times\boldsymbol{B} % Lorentz Force
        + \frac{\mathrm{E}}{\tilde{\rho(r)}}\boldsymbol{\nabla}\cdot\boldsymbol{\mathcal{D}} \;\;\; &\mathrm{Momentum}\\
       %
       %
       \tilde{\rho}(r)\,\tilde{T}(r)\left[\frac{\partial \Theta}{\partial t} + \boldsymbol{v}\cdot\boldsymbol{\nabla}\Theta \right] &=
       \frac{\mathrm{E}}{\mathrm{Pr}}\boldsymbol{\nabla}\cdot\left[\tilde{\rho}(r)\,\tilde{T}(r)\,\boldsymbol{\nabla}\Theta \right] % diffusion
       + Q(r)   % Internal heating
       + \frac{\mathrm{E}\,\mathrm{Di}}{\mathrm{Ra}^*}\Phi(r,\theta,\phi)
       + \frac{\mathrm{Di\,E^2}}{\mathrm{Pm}^2\mathrm{Ra}^*}\left[\boldsymbol{\nabla}\times\boldsymbol{B}\right]^2 &\mathrm{Thermal\; Energy}\\ % Ohmic Heating
       %
       %
       \frac{\partial \boldsymbol{B}}{\partial t} &= \boldsymbol{\nabla}\times\left(\,\boldsymbol{v}\times\boldsymbol{B}-\frac{\mathrm{E}}{\mathrm{Pm}}\boldsymbol{\nabla}\times\boldsymbol{B}\,\right) &\mathrm{Induction} \\
       %
       %
       \mathcal{D}_{ij} &= 2\tilde{\rho}(r)\left[e_{ij} - \frac{1}{3}\boldsymbol{\nabla}\cdot\boldsymbol{v}\right] &\mathrm{Viscous\; Stress\; Tensor}\\
       %
       %
       \Phi(r,\theta,\phi) &= 2\,\tilde{\rho}(r)\left[e_{ij}e_{ij} - \frac{1}{3}\left(\boldsymbol{\nabla}\cdot\boldsymbol{v}\right)^2\right] &\mathrm{Viscous\; Heating} \\
       %
       %
       \boldsymbol{\nabla}\cdot\left(\tilde{\rho}(r)\,\boldsymbol{v}\right)&=0 &\mathrm{Solenoidal\; Mass\; Flux}\\
       \boldsymbol{\nabla}\cdot\boldsymbol{B}&=0. &\mathrm{Solenoidal\; Magnetic\; Field}\end{aligned}
