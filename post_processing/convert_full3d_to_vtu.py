#!/usr/bin/env python

###############################################
#
#  Converts Spherical_3D output to vtus for visualization in paraview

from rayleigh_diagnostics import Spherical_3D
import vtk
import numpy as np
import glob
import argparse
import os

parser = argparse.ArgumentParser( \
                       description="""Convert binary spherical 3d output to vtus.""")
parser.add_argument('-i', '--initial', action='store', type=float, dest='initial', required=False, default=None, 
                    help='initial index to plot from')
parser.add_argument('-f', '--final', action='store', type=float, dest='final', required=False, default=None, 
                    help='final index to plot to')
args = parser.parse_args()

files = glob.glob1("Spherical_3D", "*_grid")
files = [f for f in files if (args.initial is None or int(f.split("_")[0]) >= args.initial) \
                          and (args.final is None or int(f.split("_")[0]) <= args.final)]
files = sorted(files, key = lambda f: int(f.split("_")[0]))

f3d = Spherical_3D(files[0])
ntheta = f3d.ntheta
nphi = f3d.nphi
nr = f3d.nr

rs = f3d.rs
thetas = f3d.thetas
phis = np.linspace(0, 2*np.pi, nphi+1)

VTK_HEXAHEDRON = vtk.vtkHexahedron().GetCellType()

ugrid = vtk.vtkUnstructuredGrid()

# Add the points
points = vtk.vtkPoints()
points.SetDataTypeToDouble()
index = 0

xyzToNode = [[[] for j in xrange(ntheta)] for i in xrange(nr)]
for i, r in enumerate(rs):
  for j, theta in enumerate(thetas):
    for k, phi in enumerate(phis[:-1]):
      points.InsertNextPoint(r*np.cos(phi)*np.sin(theta), \
                             r*np.sin(phi)*np.sin(theta), \
                             r*np.cos(theta))
      xyzToNode[i][j].append(index)
      index += 1
for i in xrange(nr):
  for j in xrange(ntheta):
    xyzToNode[i][j].append(xyzToNode[i][j][0])

ugrid.SetPoints(points)

# Add the volume elements
for i in xrange(nr-1):
  for j in xrange(ntheta-1):
    for k in xrange(nphi):
      idList = vtk.vtkIdList()
      idList.InsertNextId(xyzToNode[i][j][k])
      idList.InsertNextId(xyzToNode[i][j+1][k])
      idList.InsertNextId(xyzToNode[i][j+1][k+1])
      idList.InsertNextId(xyzToNode[i][j][k+1])
      idList.InsertNextId(xyzToNode[i+1][j][k])
      idList.InsertNextId(xyzToNode[i+1][j+1][k])
      idList.InsertNextId(xyzToNode[i+1][j+1][k+1])
      idList.InsertNextId(xyzToNode[i+1][j][k+1])
      ugrid.InsertNextCell(VTK_HEXAHEDRON, idList)

# Add the fields
pointdata = ugrid.GetPointData()

# NOTE: we cache the grid to speed up the conversion - not valid with variable meshes

for f in files:

  f3d = Spherical_3D(f)

  for index, vals in f3d.vals.iteritems():
    data = vtk.vtkDoubleArray()
    data.SetNumberOfValues(len(vals))
    data.SetName('f'+index)
    for i, v in enumerate(vals):
      data.SetValue(i, v)
    pointdata.AddArray(data)
    pointdata.SetActiveScalars('f'+index)

  gridwriter = vtk.vtkXMLUnstructuredGridWriter()
  gridwriter.SetFileName(os.path.join("Spherical_3D", "full3d_"+f3d.basefilename+".vtu"))
  if vtk.vtkVersion.GetVTKMajorVersion() <= 5:
    gridwriter.SetInput(ugrid)
  else:
    gridwriter.SetInputData(ugrid)

  gridwriter.Write()


