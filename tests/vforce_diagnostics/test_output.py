#!/usr/bin/env python3

from rayleigh_diagnostics import G_Avgs, Point_Probes, build_file_list
import numpy as np
import os
import sys
import json

regen = "-r" in sys.argv or "--regenerate" in sys.argv

vals = {
        "viscous_force_r"      : 1228,
        "viscous_force_theta"  : 1229,
        "viscous_force_phi"    : 1230,
        "viscous_pforce_r"     : 1231,
        "viscous_pforce_theta" : 1232,
        "viscous_pforce_phi"   : 1233,
        "viscous_mforce_r"     : 1234,
        "viscous_mforce_theta" : 1235,
        "viscous_mforce_phi"   : 1236,
       }

def get_results(type, dir):
  results = {
              "times"     : [],
            }
  for k in vals.keys(): results[k] = []
  files = build_file_list(0,1000000,path=dir)
  a = type(files[-1], path='')
  results["times"].append(a.time[-1])
  for k,v in vals.items():
      results[k].append(a.vals[-1, a.lut[v]].tolist())
  return results

gavgs = get_results(G_Avgs, 'G_Avgs')

if regen:
  # write results
  with open("gavgs.json", "w", encoding='utf-8') as f:
    json.dump(gavgs, f)

ptps = get_results(Point_Probes, 'Point_Probes')

if regen:
  # write results
  with open("ptps.json", "w", encoding='utf-8') as f:
    json.dump(ptps, f)

if regen: sys.exit(0)

def check_results(results, filename, tol=1.e-8):
  error = False
  
  # open old results
  with open(filename, "r") as f:
    results_old = json.load(f)
  
  # compare new and old results
  for k, v in results.items():
    err = np.abs(np.asarray(v) - np.asarray(results_old[k])).max()
    if err > tol:
      print("ERROR: different "+k+" produced compared to "+filename+" ("+repr(err)+" > "+repr(tol)+")!")
      error = True
  
  return error

error = check_results(gavgs, "gavgs.json")
error = check_results(ptps, "ptps.json") or error

if error: sys.exit(1)
sys.exit(0)
