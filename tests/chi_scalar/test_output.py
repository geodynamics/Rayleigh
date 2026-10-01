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

# The keys in `vals` that make sense for a run with no active/passive
# scalars declared at all (T.restart_no_scalars).
hydro_keys = [k for k in vals.keys() if not k.startswith("chi")]

def check_same_trajectory(dir_a, dir_b, iteration, tol=1.e-10, keys=None):
  # Unlike check_results (which cross-checks a "T-active" run against a
  # "chi-active" run with swapped Rayleigh numbers), this does a direct,
  # unswapped per-key comparison of dir_a and dir_b's G_Avgs at the same
  # iteration -- for two directories expected to produce the identical
  # trajectory. `keys` restricts the comparison to a subset of `vals` (e.g.
  # hydro_keys).
  if keys is None: keys = vals.keys()
  error = False
  fname = "{:08d}".format(iteration)
  results = {"times": []}
  for k in keys: results[k] = []
  for d in [dir_a, dir_b]:
    a = G_Avgs(os.path.join(d, 'G_Avgs', fname), path='')
    results["times"].append(a.time[-1])
    for k in keys:
      results[k].append(a.vals[-1, a.lut[vals[k]]])

  for k,v in results.items():
    print(k+":\t", v)
    if np.any(np.abs(np.asarray(v)-v[0]) > tol):
      print("ERROR: "+dir_b+" does not match "+dir_a+"'s "+k+
            " at iteration "+str(iteration)+" (tolerance "+repr(tol)+")!")
      error = True

  return error

error_hydro = check_results(["T", "chi"])
error_hydro_restart = check_results(["T.check", "chi.check"])

# chi.magnetic/chi.magnetic.check repeat the chi/chi.check setup with
# magnetism = .true. (magnetic field held at zero) to test the
# magnetism + active/passive scalar branch of Read_Checkpoint, which
# chi.check alone does not (magnetism=.false. there). 
# chi.magnetic runs straight through to
# iteration 18 (forcing a checkpoint at iteration 17);
# chi.magnetic.check restarts from that checkpoint and takes one more step to
# iteration 18. The two are compared directly.
error_magnetic_restart = check_same_trajectory("chi.magnetic", "chi.magnetic.check", 18)

# With B held at zero, chi.magnetic should reproduce chi's trajectory
# exactly (bit-for-bit): magnetism=.true. only adds the induction
# equation's own equations/coefficients, none of which should touch the
# hydro/scalar fields when B == 0. This catches bugs that corrupt every
# active-scalar run under magnetism regardless of restarting. 
# chi.magnetic.check vs chi.magnetic alone would NOT catch
# this, since both sides suffer the identical corruption and still agree
# with each other.
error_magnetic_vs_hydro = check_same_trajectory("chi", "chi.magnetic", 17)

# chi.restart_magnetism_on restarts from chi's checkpoint (magnetism=.false.
# there) with magnetism turned on at the restart (B held at zero, since
# chi's checkpoint has no C/A data to restart from anyway). Since B == 0 on
# both sides, this should be indistinguishable from chi.magnetic (magnetism
# on for the whole run) at the same iteration -- testing Read_Checkpoint
# reading a hydro+scalar-only checkpoint into a magnetism=.true. run
# (read_magnetism=0, so the checkpoint's absent C/A files are never
# attempted). Compared at iteration 19, not 18: iteration 18 just re-reports
# the checkpoint's own restarted-from state, so a genuinely new step (where a
# bug would actually show up) only appears at iteration 19.
error_restart_add_magnetism = check_same_trajectory("chi.magnetic", "chi.restart_magnetism_on", 19)

# chi.magnetic.restart_magnetism_off is the reverse: restarts from
# chi.magnetic's checkpoint (magnetism=.true., B == 0, there) with magnetism
# turned off. magnetism=.false. takes the unconditional read_var(:) = 1
# path, but the checkpoint_suffix built for *this* run has no C/A entries at
# all, so its extra C/A files are simply never looked at. Also compared
# against chi.magnetic's iteration-19 result, for the same reason as above.
error_restart_drop_magnetism = check_same_trajectory("chi.magnetic", "chi.magnetic.restart_magnetism_off", 19)

# T.restart_no_scalars restarts from T's checkpoint (n_active_scalars=2,
# n_passive_scalars=2, both Ra=0/passive) into a run with no scalars
# at all. Compared against T.check's iteration-18 result, hydro
# keys only.
error_restart_drop_scalars = check_same_trajectory("T.check", "T.restart_no_scalars", 18, keys=hydro_keys)

error = (error_hydro or error_hydro_restart or error_magnetic_restart or
          error_magnetic_vs_hydro or error_restart_add_magnetism or
          error_restart_drop_magnetism or error_restart_drop_scalars)

if error: sys.exit(1)
sys.exit(0)

