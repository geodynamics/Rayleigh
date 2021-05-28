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

The key can be anything, but must be unique (prevents overwriting previous
definitions). The value should be an available quantity name.

Best practice:
    1) use your favorite source to find the quantity code: qcode
    2) use lookup function to convert qcode into a usable quantity name.
       This cannot be done from within this file in the source code below,
       because python will complain about a circular dependency. The best
       way to perform this function is with the steps:
           a) change directories to where ever this file lives:
                  cd /path/to/Rayleigh/post_processing
           b) start python: python or python3 (will be system dependent)
           c) from within python, import the lut module:
                  import lut
           d) convert your chosen quantity code to a valid name:
                  qname = lut.lookup(qcode)
    3) define your shortcut in the source code below:
           shortcuts['my_first_shortcut_name'] = 'contents_of_variable_qname'

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

