#!/usr/bin/env python

from rayleigh_diagnostics import Spherical_3D
import numpy as np
import sys

base = Spherical_3D('00000001_0501', path='base/Spherical_3D/')
script = Spherical_3D('00000001_0501', path='script/Spherical_3D/')
maxabsdiff = np.abs(base.vals - script.vals).max() 
print("Maximum difference = {}".format(maxabsdiff,))
if maxabsdiff > 1.e-10:
  print("ERROR: init_type 1 and init_type 8 produced different initial conditions (within a tolerance of 1.e-10)!")
  sys.exit(1)
sys.exit(0)



