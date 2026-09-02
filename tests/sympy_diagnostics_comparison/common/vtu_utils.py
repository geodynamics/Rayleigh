"""
Convert Rayleigh's Spherical_3D output and analytic expression (and their difference) to a VTU for debugging.
"""

import os
import sys
import numpy as np
import vtk
from vtk.util import numpy_support

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..', '..', 'post_processing'))
from rayleigh_diagnostics import Spherical_3D_multi  # noqa: E402


def write_comparison_vtu(quantity_codes, numeric, full3d_path='Spherical_3D/', out_path=None, index='00000001'):
    """quantity_codes: {code: name}. numeric: {name: callable(r, theta, phi)}.

    Writes a .vtu (default: <full3d_path>/comparison_<index>.vtu) with, for
    each name, three point-data arrays: rayleigh_<name>, analytic_<name>, and
    diff_<name> = rayleigh_<name> - analytic_<name>. Returns the path written.
    """
    f3d = Spherical_3D_multi(index, path=full3d_path)
    rs, thetas = f3d.rs, f3d.thetas
    nr, ntheta, nphi = len(rs), len(thetas), f3d.nphi
    phis = np.linspace(0, 2 * np.pi, nphi + 1)

    VTK_HEXAHEDRON = vtk.vtkHexahedron().GetCellType()
    ugrid = vtk.vtkUnstructuredGrid()

    # -------------------------------------------------------------
    # Grid (points + hexahedral cells), identical to convert_full3d_to_vtu.py
    # -------------------------------------------------------------
    points = vtk.vtkPoints()
    points.SetDataTypeToDouble()
    idx = 0
    xyzToNode = [[[] for _ in range(ntheta)] for _ in range(nr)]
    for i, r in enumerate(rs):
        for j, theta in enumerate(thetas):
            for k, phi in enumerate(phis[:-1]):
                points.InsertNextPoint(r * np.cos(phi) * np.sin(theta),
                                        r * np.sin(phi) * np.sin(theta),
                                        r * np.cos(theta))
                xyzToNode[i][j].append(idx)
                idx += 1
    for i in range(nr):
        for j in range(ntheta):
            xyzToNode[i][j].append(xyzToNode[i][j][0])

    ugrid.SetPoints(points)

    for i in range(nr - 1):
        for j in range(ntheta - 1):
            for k in range(nphi):
                idList = vtk.vtkIdList()
                idList.InsertNextId(xyzToNode[i][j][k])
                idList.InsertNextId(xyzToNode[i][j + 1][k])
                idList.InsertNextId(xyzToNode[i][j + 1][k + 1])
                idList.InsertNextId(xyzToNode[i][j][k + 1])
                idList.InsertNextId(xyzToNode[i + 1][j][k])
                idList.InsertNextId(xyzToNode[i + 1][j + 1][k])
                idList.InsertNextId(xyzToNode[i + 1][j + 1][k + 1])
                idList.InsertNextId(xyzToNode[i + 1][j][k + 1])
                ugrid.InsertNextCell(VTK_HEXAHEDRON, idList)

    # -------------------------------------------------------------
    # Fields: rayleigh_/analytic_/diff_ per quantity. f3d.vals[code] has
    # shape (nphi, ntheta, nr) (see rayleigh_diagnostics.py / compare_utils.py's
    # meshgrid convention); flatten(order='F') makes phi vary fastest, then
    # theta, then r -- matching the point-insertion order above exactly, as
    # in convert_full3d_to_vtu.py.
    # -------------------------------------------------------------
    pointdata = ugrid.GetPointData()
    PHI, THETA, RADIUS = np.meshgrid(phis[:-1], thetas, rs, indexing='ij')

    for code, name in sorted(quantity_codes.items()):
        rayleigh_vals = np.asarray(f3d.vals[f'{code:05d}'], dtype='float64')
        analytic_vals = np.broadcast_to(numeric[name](RADIUS, THETA, PHI), rayleigh_vals.shape).astype('float64')
        diff_vals = rayleigh_vals - analytic_vals

        for prefix, vals in (('rayleigh_', rayleigh_vals), ('analytic_', analytic_vals), ('diff_', diff_vals)):
            flat = np.ascontiguousarray(vals.flatten(order='F'))
            data = numpy_support.numpy_to_vtk(flat, deep=True)
            data.SetName(prefix + name)
            pointdata.AddArray(data)

    if out_path is None:
        out_path = os.path.join(full3d_path, f'comparison_{f3d.basefilename}.vtu')

    gridwriter = vtk.vtkXMLUnstructuredGridWriter()
    gridwriter.SetFileName(out_path)
    gridwriter.SetInputData(ugrid)
    gridwriter.Write()
    print(f"wrote {out_path}")
    return out_path
