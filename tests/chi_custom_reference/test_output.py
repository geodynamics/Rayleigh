#!/usr/bin/env python3

from rayleigh_diagnostics import G_Avgs, build_file_list
import numpy as np
import os
import sys

#aspect ratio
ar = 0.35e0
# shell depth 
d=1.e0 

#outer radial boundary
rmax=d/(1-ar)

#inner radial boundary
rmin=ar*rmax

vals = {
        "kesgav"       : (401, 0.0),
        "Tsgav"        : (501, (rmax**3*rmin/6 - rmax*rmin**3/2 + rmin**4/3)*3/(rmax**3-rmin**3)),
        "Tpsgav"       : (503, 0.0),
        "Tmsgav"       : (505, (rmax**3*rmin/6 - rmax*rmin**3/2 + rmin**4/3)*3/(rmax**3-rmin**3)),
        "Tdrsgav"      : (507, -rmax*rmin*(rmax-rmin)*3/(rmax**3-rmin**3)),
        "Tpdrsgav"     : (509, 0.0),
        "Tmdrsgav"     : (511, -rmax*rmin*(rmax-rmin)*3/(rmax**3-rmin**3)),
        "Tdtsgav"      : (513, 0.0),
        "Tpdtsgav"     : (515, 0.0),
        "Tmdtsgav"     : (517, 0.0),
        "Tdpsgav"      : (519, 0.0),
        "Tpdpsgav"     : (521, 0.0),
        "Tmdpsgav"     : (523, 0.0),
        "chia1sgav"    : (10001, (rmax**3*rmin/6 - rmax*rmin**3/2 + rmin**4/3)*3/(rmax**3-rmin**3)),
        "chia1psgav"   : (10002, 0.0),
        "chia1msgav"   : (10003, (rmax**3*rmin/6 - rmax*rmin**3/2 + rmin**4/3)*3/(rmax**3-rmin**3)),
        "chia1drsgav"  : (10004, -rmax*rmin*(rmax-rmin)*3/(rmax**3-rmin**3)),
        "chia1pdrsgav" : (10005, 0.0),
        "chia1mdrsgav" : (10006, -rmax*rmin*(rmax-rmin)*3/(rmax**3-rmin**3)),
        "chia1dtsgav"  : (10007, 0.0),
        "chia1pdtsgav" : (10008, 0.0),
        "chia1mdtsgav" : (10009, 0.0),
        "chia1dpsgav"  : (10010, 0.0),
        "chia1pdpsgav" : (10011, 0.0),
        "chia1mdpsgav" : (10012, 0.0),
       }

def check_results(dirs, tol=1.e-10):
  error = False
  results = {
             "times"     : [],
            }
  for k in vals.keys(): results[k] = []
  for d in dirs:
    files = build_file_list(0,1000000,path=os.path.join(d,'G_Avgs'))
    a = G_Avgs(files[-1], path='')
    results["times"].append(a.time[-1])
    for k,(v,off) in vals.items():
      try:
        val = a.vals[-1, a.lut[v]]
        if d.endswith('full'): val -= off
        results[k].append(val)
      except IndexError:
        import ipdb; ipdb.set_trace()

  for k,v in results.items():
    print(k+":\t", v)
    if np.any(np.abs(v-v[0]) > tol):
      print("ERROR: different "+k+" produced between runs (within a tolerance of "+repr(tol)+")!")
      error = True
    if k.startswith("chi"):
      Tk = "T"+ k[5:]
      if np.any(np.abs(np.asarray(v)-np.asarray(results[Tk])) > tol):
        print("ERROR: different "+k+" and "+Tk+" produced (with a tolerance of "+repr(tol)+")!")
        error = True

  return error

error = check_results(["T.full", "chi.full"]) or check_results(["T.split", "chi.split"])

if error: sys.exit(1)
sys.exit(0)

