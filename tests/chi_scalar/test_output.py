#!/usr/bin/env python

from rayleigh_diagnostics import G_Avgs, build_file_list
from rayleigh_spectral_input import SpectralInput, radial_extents
import numpy as np
import sys
import vtk

error = False

#benchfiles = build_file_list(0,1000000,path='bench/G_Avgs')
#bench = G_Avgs(benchfiles[-1], path='')
#kebench = bench.vals[-1, bench.lut[401]]

Tfiles = build_file_list(0,1000000,path='T/G_Avgs')
T = G_Avgs(Tfiles[-1], path='')
keT = T.vals[-1, T.lut[401]]

chifiles = build_file_list(0,1000000,path='chi/G_Avgs')
chi = G_Avgs(chifiles[-1], path='')
kechi = chi.vals[-1, chi.lut[401]]

print("times: ", T.time[-1], chi.time[-1])
print("KEs: ", keT, kechi)

absdiff = abs(keT-kechi)

if absdiff > 1.e-10:
  print("ERROR: chi scalar produced a different answer to the temperature case (within a tolerance of 1.e-10)!")
  error = True

if error:
  sys.exit(1)
sys.exit(0)

