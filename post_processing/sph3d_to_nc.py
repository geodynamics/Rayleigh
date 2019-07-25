#!/usr/bin/env python
import sys
import os
import numpy as np
import netCDF4 as nc

import rayleigh_diagnostics as rd

if len(sys.argv) < 2:
    exe = os.path.basename(sys.argv[0])
    print("{} [Spherical_3D files to convert]".format(exe))
    print("converts Spherical_3D output to netCDF4 files with CF convections for spherical coordinates")
    print("example call:")
    print("    {} Spherical_3D/*_0*".format(sys.argv[0]))

def indexonly(f):
    fsplit = os.path.basename(f).split('_')
    if len(fsplit) > 1:
        fsplit = fsplit[:-1]
    return os.path.join(os.path.dirname(f), ''.join(fsplit))

files = [indexonly(x) for x in sys.argv[1:]]
files = list(set(files))
files.sort()

dims = ('rad', 'lat', 'lon')

for f in files:
    print('processing ' + f)
    # The prefix is needed for automatic file grouping in ParaView.
    ncfile = os.path.join(os.path.dirname(f), 'spherical.' + os.path.basename(f) + '.nc')
    with nc.Dataset(ncfile, mode='w', format='NETCDF4') as d:
        sp = rd.Spherical_3D_multi(os.path.basename(f), os.path.dirname(f) + '/')
        d.createDimension('rad', sp.nr)
        d.createDimension('lat', sp.ntheta)
        d.createDimension('lon', sp.nphi)
        d.createDimension('nv', 2)

        rvar = d.createVariable('rad', np.float64, dimensions=('rad',))
        rvar[:] = sp.rs[::-1]
        rvar.units = 'cm'
        rvar.standard_name = 'height'

        latvar = d.createVariable('lat', np.float64, dimensions=('lat',))
        latvar.units = 'degrees_north'
        latvar.standard_name = 'latitude'
        latvar.bounds = 'lat_bnds'
        latvar[:] = sp.thetas[::-1] * 180.0 / np.pi - 90.0

        latbndsvar = d.createVariable('lat_bnds', np.float64, dimensions=('lat','nv'))
        latbndsvar.units = 'degrees_north'
        latbnds = 0.5 * (latvar[1:] + latvar[:-1])
        latbndsvar[1:, 0] = latbnds
        latbndsvar[:-1, 1] = latbnds
        latbndsvar[0, 0] = -90.0
        latbndsvar[-1, 1] = 90.0

        lonbndsvar = d.createVariable('lon_bnds', np.float64, dimensions=('lon','nv'))
        lonbndsvar.units = 'degrees_north'
        lonbnds = np.linspace(-180., 180., sp.nphi+1)
        lonbndsvar[:,0] = lonbnds[:-1] # left edge
        lonbndsvar[:,1] = lonbnds[1:] # right edge

        lonvar = d.createVariable('lon', np.float64, dimensions=('lon',))
        lonvar.units = 'degrees_east'
        lonvar.standard_name = 'longitude'
        lonvar.bounds = 'lon_bnds'
        lonvar[:] = 0.5 * (lonbnds[1:] + lonbnds[:-1])

        for v in sp.vals:
            var = d.createVariable(v, np.float32, dimensions=dims)
            var[:] = np.swapaxes(sp.vals[v][:, ::-1, ::-1], 0, 2)
