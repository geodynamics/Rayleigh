#!/usr/bin/env python

from rayleigh_diagnostics import Spherical_3D_multi, Shell_Slices
from rayleigh_spectral_input import SpectralInput, radial_extents
import numpy as np
import sys
import vtk

error = False

# test using the raw data imported from the output files
base = Spherical_3D_multi('00000001_0501', path='base/Spherical_3D/')
script = Spherical_3D_multi('00000001_0501', path='script/Spherical_3D/')
maxabsdiff = 0.0
for k in base.vals.keys():
  maxabsdiffk = np.abs(base.vals[k] - script.vals[k]).max() 
  print("Maximum difference ({}) = {}".format(k, maxabsdiffk,))
  maxabsdiff = max(maxabsdiff, maxabsdiffk)

if maxabsdiff > 1.e-10:
  print("ERROR: [magnetic_]init_type 1 and [magnetic_]init_type 8 produced different initial conditions (within a tolerance of 1.e-10)!")
  error = True

# test using the data written to vtk by the convert script
gridreader = vtk.vtkXMLUnstructuredGridReader()

gridreader.SetFileName("base/Spherical_3D/full3d_00000001.vtu")
gridreader.Update()
ugradbase = gridreader.GetOutput()

gridreader.SetFileName("script/Spherical_3D/full3d_00000001.vtu")
gridreader.Update()
ugridscript = gridreader.GetOutput()

pdatabase = ugradbase.GetPointData()
fnames = [pdatabase.GetArrayName(i) for i in range(pdatabase.GetNumberOfArrays())]

pdatascript = ugridscript.GetPointData()

maxabsdiff = 0.0
for k in fnames:
  sdatabase = pdatabase.GetScalars(k)
  database = np.asarray([sdatabase.GetTuple1(i) for i in range(sdatabase.GetNumberOfTuples())])
  sdatascript = pdatascript.GetScalars(k)
  datascript = np.asarray([sdatascript.GetTuple1(i) for i in range(sdatascript.GetNumberOfTuples())])
  maxabsdiffk = np.abs(database - datascript).max()
  maxabsdiff = max(maxabsdiff, maxabsdiffk)
  print("Maximum vtk difference ({}) = {}".format(k, maxabsdiffk,))

if maxabsdiff > 1.e-10:
  print("ERROR: [magnetic_]init_type 1 and [magnetic_]init_type 8 produced different initial vtks (within a tolerance of 1.e-10)!")
  error = True

# test bcs using generic input and hard coded values
bcs_base = Shell_Slices('00000002', path='bcs_base/Shell_Slices/')
bcs_script = Shell_Slices('00000002', path='bcs_script/Shell_Slices/')

maxabsdiff = 0.0
for k in bcs_base.qv:
  i_b = bcs_base.lut[k]
  i_s = bcs_script.lut[k]
  for ri in range(2):
    maxabsdiffk = np.abs(bcs_base.vals[:,:,ri,i_b,0] - bcs_script.vals[:,:,ri,i_s,0]).max()
    print("Maximum bc {} difference from base ({}) = {}".format(ri, k, maxabsdiffk,))
    maxabsdiff = max(maxabsdiff, maxabsdiffk)
    if k == 501:
      if ri == 0:
        maxabsdiffk = np.abs(bcs_script.vals[:,:,ri,i_s,0] + 5.).max()
      else:
        maxabsdiffk = np.abs(bcs_script.vals[:,:,ri,i_s,0] - 5.).max()
      print("Maximum T bc {} difference from input ({}) = {}".format(ri, k, maxabsdiffk,))

if maxabsdiff > 1.e-10:
  print("ERROR: generic input bc produced an unexpected result (within a tolerance of 1.e-10)!")
  error = True

# test using the raw data imported from the output files for radial case
base_r = Spherical_3D_multi('00000001_0501', path='radial_base/Spherical_3D/')
sparse_r = Spherical_3D_multi('00000001_0501', path='radial_sparse/Spherical_3D/')
maxabsdiff = 0.0
for k in base_r.vals.keys():
  maxabsdiffk = np.abs(base_r.vals[k] - sparse_r.vals[k]).max() 
  print("Maximum difference ({}) = {}".format(k, maxabsdiffk,))
  maxabsdiff = max(maxabsdiff, maxabsdiffk)

if maxabsdiff > 1.e-10:
  print("ERROR: init_type 7 and init_type 8 (sparse) produced different initial conditions (within a tolerance of 1.e-10)!")
  error = True

# test using the raw data imported from the output files for radial case
dense_r = Spherical_3D_multi('00000001_0501', path='radial_dense/Spherical_3D/')
maxabsdiff = 0.0
for k in base_r.vals.keys():
  maxabsdiffk = np.abs(base_r.vals[k] - dense_r.vals[k]).max() 
  print("Maximum difference ({}) = {}".format(k, maxabsdiffk,))
  maxabsdiff = max(maxabsdiff, maxabsdiffk)

if maxabsdiff > 1.e-10:
  print("ERROR: init_type 7 and init_type 8 (dense) produced different initial conditions (within a tolerance of 1.e-10)!")
  error = True

# test spectral input reading and inverse transforming
maxabsdiff = 0.0

si_w = SpectralInput()
def five():
  return 5.0
si_w.transform_from_rtp_function(five)
si_r = SpectralInput()
si_r.read('bcs_script/five_init_bc')
maxabsdiffk = np.abs(si_r.coeffs-si_w.coeffs).max()
print("Maximum difference on sparse coefficient read in = {}".format(maxabsdiffk,))
maxabsdiff = max(maxabsdiff, maxabsdiffk)
maxabsdiffk = np.abs(si_r.inverse_transform() - 5.0).max()
print("Maximum difference on sparse coefficient read in and inverse transform = {}".format(maxabsdiffk,))
maxabsdiff = max(maxabsdiff, maxabsdiffk)

si_w = SpectralInput(n_theta=64, n_r=48)
ar = 0.35
sd = 1.0
rmin, rmax = radial_extents(aspect_ratio=ar, shell_depth=sd)
def t_bench(radius, theta, phi):
  x = 2*radius - rmin - rmax
  return rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) - x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)
si_w.transform_from_rtp_function(t_bench, aspect_ratio=ar, shell_depth=sd)
si_r = SpectralInput(n_theta=64, n_r=48)
si_r.read('script/bench_t_init')
maxabsdiffk = np.abs(si_w.coeffs-si_r.coeffs).max()
print("Maximum difference on dense coefficient read in = {}".format(maxabsdiffk,))
maxabsdiff = max(maxabsdiff, maxabsdiffk)

if maxabsdiff > 1.e-10:
  print("ERROR: generic input read in produced an unexpected result (within a tolerance of 1.e-10)!")
  error = True

if error:
  sys.exit(1)
sys.exit(0)

