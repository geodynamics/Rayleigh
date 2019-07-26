#!/usr/bin/env python

from rayleigh_diagnostics import Spherical_3D_multi
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

if error:
  sys.exit(1)
sys.exit(0)



