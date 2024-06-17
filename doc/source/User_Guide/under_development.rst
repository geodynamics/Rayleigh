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




