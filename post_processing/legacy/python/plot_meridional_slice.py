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

# Meridional Slice Plotting Example.
# Displays a single variable from the meridional_slice data structure
# (in the r-theta plane)
#
#    """Rayleigh Meridional Slice Structure
#    ----------------------------------
#    self.niter                                    : number of time steps
#    self.nq                                       : number of diagnostic quantities output
#    self.nr                                       : number of radial points
#    self.ntheta                                   : number of theta points
#    self.nphi                                     : number of phi points sampled
#    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
#    self.radius[0:nr-1]                           : radial grid
#    self.costheta[0:ntheta-1]                     : cos(theta grid)
#    self.sintheta[0:ntheta-1]                     : sin(theta grid)
#    self.phi[0:nphi-1]                            : phi values (radians)
#    self.phi_indices[0:nphi-1]                    : phi indices (from 1 to nphi)
#    self.vals[0:nphi-1,0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] : The meridional slices 
#    self.iters[0:niter-1]                         : The time step numbers stored in this output file
#    self.time[0:niter-1]                          : The simulation time corresponding to each time step
#    self.version                                  : The version code for this particular output (internal use)
#    self.lut                                      : Lookup table for the different diagnostics output
#    """

from diagnostic_reading import Meridional_Slice
import numpy as np
import matplotlib.pyplot as plt

timestep = '00005000'
quantity_code = 64  # read in temperature
remove_mean = True  # remove the ell=0 mean
tindex = 0          # Display the first timestep from the file
pindex = 4          # Display the 5th phi-value output

a = Meridional_Slice(timestep)


#Set up the grid
nr = a.nr
ntheta = a.ntheta
r = a.radius/np.max(a.radius)
# We apply a shift in theta so that theta=0 corresponds to equator and the poles appear where we expect them to...
theta = np.arccos(a.costheta)-np.pi/2  

radius_matrix, theta_matrix = np.meshgrid(r,theta)

X = radius_matrix * np.cos(theta_matrix)
Y = radius_matrix * np.sin(theta_matrix)


qindex = a.lut[quantity_code]
field = np.zeros((ntheta,nr),dtype='float64')
field[:,:] =a.vals[pindex,:,:,qindex,tindex]


#remove the mean if desired
if (remove_mean):
    for i in range(nr):
        the_mean = np.mean(field[:,i])
        field[:,i] = field[:,i]-the_mean
radtodeg = 180.0/np.pi
print 'Displaying meridional slice at phi (radians, degrees) = ', a.phi[pindex], a.phi[pindex]*radtodeg

#Plot
plt.figure(1)
img = plt.pcolormesh(X,Y,field,cmap='jet')
plt.show()
