#!/usr/bin/env python

#
#  Copyright (C) 2019 by the authors of the RAYLEIGH code.
#
#  This file is part of RAYLEIGH.
#
#  RAYLEIGH is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3, or (at your option)
#  any later version.
#
#  RAYLEIGH is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with RAYLEIGH; see the file LICENSE.  If not see
#  <http://www.gnu.org/licenses/>.

###############################################
#
#  Converts Spherical_3D output to netCDF files for visualization in ParaView
#  and other software.

import sys
import os
import argparse
import glob
import numpy as np
import netCDF4 as nc

import rayleigh_diagnostics as rd

parser = argparse.ArgumentParser(description="Converts Spherical_3D output to netCDF4 files according to CF convections for spherical coordinates.")
parser.add_argument('-a', '--all', action='store_true', help='run on all files in Spherical_3D')
parser.add_argument('-u', '--update', action='store_true', help='forces update on already existing files')
parser.add_argument('-v', '--verbose', action='store_true', help='print processed files')
parser.add_argument('files', metavar='filename', nargs='*', help='files to convert')
args = parser.parse_args()

if args.all:
    files = glob.glob('Spherical_3D/*_grid') + args.files
else:
    files = []
files += args.files
files = list(set(files)) # filter out duplicates
files.sort()

if len(files) == 0:
    parser.print_help()
    parser.exit()

dims = ('rad', 'lat', 'lon')

for f in files:
    # The prefix is needed for automatic file grouping in ParaView.
    ncfile = os.path.join(os.path.dirname(f), 'spherical.' + os.path.basename(f) + '.nc')

    if not args.update and os.path.exists(ncfile):
        if args.verbose:
            print('skipping already existing file ' + f)
        continue

    if args.verbose:
        print('processing ' + f)

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
            var = d.createVariable('f' + v, np.float32, dimensions=dims)
            var[:] = np.swapaxes(sp.vals[v][:, ::-1, ::-1], 0, 2)
