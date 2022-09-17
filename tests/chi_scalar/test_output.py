#!/usr/bin/env python3

from rayleigh_diagnostics import G_Avgs, build_file_list
from rayleigh_spectral_input import SpectralInput, radial_extents
import numpy as np
import os
import sys

def check_results(dirs, tol=1.e-10):
  error = False
  results = {
             "times"     : [],
             "kesgav"    : [],
             "Tsgav"     : [],
             "chia1sgav" : [],
             "chia2sgav" : [],
             "chip1sgav" : [],
             "chip2sgav" : [],
            }
  for d in dirs:
    files = build_file_list(0,1000000,path=os.path.join(d,'G_Avgs'))
    a = G_Avgs(files[-1], path='')
    results["times"].append(a.time[-1])
    results["kesgav"].append(a.vals[-1, a.lut[401]])
    results["Tsgav"].append(a.vals[-1, a.lut[501]])
    results["chia1sgav"].append(a.vals[-1, a.lut[10001]])
    results["chia2sgav"].append(a.vals[-1, a.lut[10201]])
    results["chip1sgav"].append(a.vals[-1, a.lut[20001]])
    results["chip2sgav"].append(a.vals[-1, a.lut[20201]])

  for k,v in results.items():
    print(k+":\t", v)
    if np.any(np.abs(v-v[0]) > tol):
      print("ERROR: different "+k+" produced (within a tolerance of "+repr(tol)+")!")
      error = True

  return error

error = check_results(["T", "chi"]) or check_results(["T.check", "chi.check"])

if error: sys.exit(1)
sys.exit(0)

