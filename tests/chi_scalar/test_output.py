#!/usr/bin/env python

from rayleigh_diagnostics import G_Avgs, build_file_list
from rayleigh_spectral_input import SpectralInput, radial_extents
import numpy as np
import sys
import vtk

error = False

benchfiles = build_file_list(0,1000000,path='bench/G_Avgs')
bench = G_Avgs(benchfiles[-1], path='')
kebench = bench.vals[-1, bench.lut[401]]

Tfiles = build_file_list(0,1000000,path='T/G_Avgs')
T = G_Avgs(Tfiles[-1], path='')
keT = T.vals[-1, T.lut[401]]

chifiles = build_file_list(0,1000000,path='chi/G_Avgs')
chi = G_Avgs(chifiles[-1], path='')
kechi = chi.vals[-1, chi.lut[401]]

print(bench.time[-1], T.time[-1], chi.time[-1])
print(kebench, keT, kechi, keT-kechi)

