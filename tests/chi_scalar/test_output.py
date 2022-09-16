#!/usr/bin/env python3

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
TT = T.vals[-1, T.lut[501]]
chia1T = T.vals[-1, T.lut[10001]]
chia2T = T.vals[-1, T.lut[10201]]
chip1T = T.vals[-1, T.lut[20001]]
chip2T = T.vals[-1, T.lut[20201]]

chifiles = build_file_list(0,1000000,path='chi/G_Avgs')
chi = G_Avgs(chifiles[-1], path='')
kechi = chi.vals[-1, chi.lut[401]]
Tchi = chi.vals[-1, chi.lut[501]]
chia1chi = chi.vals[-1, chi.lut[10001]]
chia2chi = chi.vals[-1, chi.lut[10201]]
chip1chi = chi.vals[-1, chi.lut[20001]]
chip2chi = chi.vals[-1, chi.lut[20201]]

print("times:  ", T.time[-1], chi.time[-1])
print("KEs:    ", keT, kechi)
print("Ts:     ", TT, Tchi)
print("chia1s: ", chia1T, chia1chi)
print("chia2s: ", chia2T, chia2chi)
print("chip1s: ", chip1T, chip1chi)
print("chip2s: ", chip2T, chip2chi)

keabsdiff  = abs(keT-kechi)
Tabsdiff   = abs(TT-Tchi)
chia1absdiff = abs(chia1T-chia1chi)
chia2absdiff = abs(chia2T-chia2chi)
chip1absdiff = abs(chip1T-chip1chi)
chip2absdiff = abs(chip2T-chip2chi)

if keabsdiff > 1.e-10:
  print("ERROR: chi scalar produced a different KE to the temperature case (within a tolerance of 1.e-10)!")
  error = True

if Tabsdiff > 1.e-10:
  print("ERROR: chi scalar produced a different T to the temperature case (within a tolerance of 1.e-10)!")
  error = True

if chia1absdiff > 1.e-10:
  print("ERROR: chi scalar produced a different chia1 to the temperature case (within a tolerance of 1.e-10)!")
  error = True

if chia2absdiff > 1.e-10:
  print("ERROR: chi scalar produced a different chia2 to the temperature case (within a tolerance of 1.e-10)!")
  error = True

if chip1absdiff > 1.e-10:
  print("ERROR: chi scalar produced a different chip1 to the temperature case (within a tolerance of 1.e-10)!")
  error = True

if chip2absdiff > 1.e-10:
  print("ERROR: chi scalar produced a different chip2 to the temperature case (within a tolerance of 1.e-10)!")
  error = True

if error:
  sys.exit(1)
sys.exit(0)

