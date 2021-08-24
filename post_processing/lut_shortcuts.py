"""
Module to hold user defined shortcuts for accessing the quantity codes

The standard way to get the quantity code for the full kinetic energy
is to call:

    ke_index = lookup('kinetic_energy')

This module allows you to define your own shortcuts such that the above
code simplifies to:

    ke_index = lookup('ke')

In order to define a shortcut, you must know the quantity name that you
want to map. Then simply define a dictionary entry with the format:

    shortcuts[key] = value

The key can be anything, but must be a lowercase string that is unique
(prevents overwriting previous definitions). The value must be an
available quantity name.

Some common mappings are included in the source code below, they include:

    ur,ut,up ---> full velocity components

    vr,vt,vp ---> full velocity components (some people prefer 'v' over 'u')

    ke ---> full kinetic energy

    me ---> full magnetic energy

"""
from collections import OrderedDict

shortcuts = OrderedDict() # this must be named "shortcuts"

# full velocity
shortcuts['ur'] = 'v_r'
shortcuts['ut'] = 'v_theta'
shortcuts['up'] = 'v_phi'
shortcuts['vr'] = 'v_r'
shortcuts['vt'] = 'v_theta'
shortcuts['vp'] = 'v_phi'

# full magnetic field
shortcuts['br'] = 'b_r'
shortcuts['bt'] = 'b_theta'
shortcuts['bp'] = 'b_phi'

# full kinetic/magnetic energies
shortcuts['ke'] = 'kinetic_energy'
shortcuts['me'] = 'magnetic_energy'

