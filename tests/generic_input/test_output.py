#!/usr/bin/env python

from rayleigh_diagnostics import Spherical_3D_multi
import numpy as np
import sys

base = Spherical_3D_multi('00000001_0501', path='base/Spherical_3D/')
script = Spherical_3D_multi('00000001_0501', path='script/Spherical_3D/')
maxabsdiff = 0.0
for k in base.vals.keys():
  maxabsdiffk = np.abs(base.vals[k] - script.vals[k]).max() 
  print("Maximum difference ({}) = {}".format(k, maxabsdiffk,))
  maxabsdiff = max(maxabsdiff, maxabsdiffk)
if maxabsdiff > 1.e-10:
  print("ERROR: [magnetic_]init_type 1 and [magnetic_]init_type 8 produced different initial conditions (within a tolerance of 1.e-10)!")
  sys.exit(1)
sys.exit(0)



