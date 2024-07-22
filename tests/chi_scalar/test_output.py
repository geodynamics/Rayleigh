#!/usr/bin/env python3

from rayleigh_diagnostics import G_Avgs, build_file_list
import numpy as np
import os
import sys

vals = {
        "kesgav"       : 401,
        "Tsgav"        : 501,
        "Tpsgav"       : 503,
        "Tmsgav"       : 505,
        "Tdrsgav"      : 507,
        "Tpdrsgav"     : 509,
        "Tmdrsgav"     : 511,
        "Tdtsgav"      : 513,
        "Tpdtsgav"     : 515,
        "Tmdtsgav"     : 517,
        "Tdpsgav"      : 519,
        "Tpdpsgav"     : 521,
        "Tmdpsgav"     : 523,
        "chia1sgav"    : 10001,
        "chia1psgav"   : 10002,
        "chia1msgav"   : 10003,
        "chia1drsgav"  : 10004,
        "chia1pdrsgav" : 10005,
        "chia1mdrsgav" : 10006,
        "chia1dtsgav"  : 10007,
        "chia1pdtsgav" : 10008,
        "chia1mdtsgav" : 10009,
        "chia1dpsgav"  : 10010,
        "chia1pdpsgav" : 10011,
        "chia1mdpsgav" : 10012,
        "chia2sgav"    : 10201,
        "chia2psgav"   : 10202,
        "chia2msgav"   : 10203,
        "chia2drsgav"  : 10204,
        "chia2pdrsgav" : 10205,
        "chia2mdrsgav" : 10206,
        "chia2dtsgav"  : 10207,
        "chia2pdtsgav" : 10208,
        "chia2mdtsgav" : 10209,
        "chia2dpsgav"  : 10210,
        "chia2pdpsgav" : 10211,
        "chia2mdpsgav" : 10212,
        "chip1sgav"    : 20001,
        "chip1psgav"   : 20002,
        "chip1msgav"   : 20003,
        "chip1drsgav"  : 20004,
        "chip1pdrsgav" : 20005,
        "chip1mdrsgav" : 20006,
        "chip1dtsgav"  : 20007,
        "chip1pdtsgav" : 20008,
        "chip1mdtsgav" : 20009,
        "chip1dpsgav"  : 20010,
        "chip1pdpsgav" : 20011,
        "chip1mdpsgav" : 20012,
        "chip2sgav"    : 20201,
        "chip2psgav"   : 20202,
        "chip2msgav"   : 20203,
        "chip2drsgav"  : 20204,
        "chip2pdrsgav" : 20205,
        "chip2mdrsgav" : 20206,
        "chip2dtsgav"  : 20207,
        "chip2pdtsgav" : 20208,
        "chip2mdtsgav" : 20209,
        "chip2dpsgav"  : 20210,
        "chip2pdpsgav" : 20211,
        "chip2mdpsgav" : 20212,
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
    for k,v in vals.items():
      try:
        results[k].append(a.vals[-1, a.lut[v]])
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

error = check_results(["T", "chi"]) or check_results(["T.check", "chi.check"])

if error: sys.exit(1)
sys.exit(0)

