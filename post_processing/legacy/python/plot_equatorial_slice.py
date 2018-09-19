#
#  Copyright (C) 2018 by the authors of the RAYLEIGH code.
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
#

#Equatorial Slice Plotting Example.
#Displays a single variable from the equatorial_slice data structure
#
# Equatorial Slice Data stucture format:
#    """Rayleigh Equatorial Slice Structure
#    ----------------------------------
#    self.niter                                    : number of time steps
#    self.nq                                       : number of diagnostic quantities output
#    self.nr                                       : number of radial points
#    self.nphi                                     : number of phi points
#    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
#    self.radius[0:nr-1]                           : radial grid
#    self.vals[0:phi-1,0:nr-1,0:nq-1,0:niter-1]    : The equatorial_slices
#    self.phi[0:nphi-1]                            : phi values (in radians)
#    self.iters[0:niter-1]                         : The time step numbers stored in this output file
#    self.time[0:niter-1]                          : The simulation time corresponding to each time step
#    self.version                                  : The version code for this particular output (internal use)
#    self.lut                                      : Lookup table for the different diagnostics output
#    """

from diagnostic_reading import Equatorial_Slice
import numpy as np
import matplotlib.pyplot as plt

timestep = '00005000'
quantity_code = 64  # read in temperature
remove_mean = True  # remove the m=0 mean
tindex = 0          # Display the first timestep from the file

a = Equatorial_Slice(timestep)

#Set up the grid
nr = a.nr
nphi = a.nphi
r = a.radius/np.max(a.radius)
phi = np.zeros(nphi+1,dtype='float64')
phi[0:nphi] = a.phi
phi[nphi] = np.pi*2  # For display purposes, it is best to have a redunant data point at 0,2pi

radius_matrix, phi_matrix = np.meshgrid(r,phi)
X = radius_matrix * np.cos(phi_matrix)
Y = radius_matrix * np.sin(phi_matrix)

qindex = a.lut[64]
field = np.zeros((nphi+1,nr),dtype='float64')
field[0:nphi,:] =a.vals[:,:,qindex,tindex]
field[nphi,:] = field[0,:]  #replicate phi=0 values at phi=2pi

#remove the mean if desired
if (remove_mean):
    for i in range(nr):
        the_mean = np.mean(field[:,i])
        field[:,i] = field[:,i]-the_mean

#Plot
plt.figure(1)
img = plt.pcolormesh(X,Y,field,cmap='jet')
plt.show()
