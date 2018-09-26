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

#  Reads in a Rayleigh shell slice file and renders data from that file
#  in 3-D using the Mayavi volume-rendering package.  You will need to have Mayavi 2
#  installed for this example to work.   In this example, one quantity is used
#  to generate a 3-D elevation map, and another quantity is used to color that map.
#  See vis_option = 1 and 2 (in vis_option 2, the same quantity is used for both)
#  You may find Mayavi 2 here:
#  http://docs.enthought.com/mayavi/mayavi/

from mayavi import mlab
import matplotlib.pylab as plt
import numpy as np
from diagnostic_reading import ShellSlice



def scale_data(indata, scale_factor = 1.5, no_mean = False):
    # Rescales indata so that:
    #   1. indata has a zero spherical mean
    #   2. indata is rescaled to have min at -1, max at +1
    #   3. indata "saturates" at absolute values > scale_factor * stdev(indata) 
    # indata should be dimensioned [theta,phi]
    ntheta = indata.shape[0]
    nphi   = indata.shape[1]
    #Subtract the spherical mean
    if (no_mean):
        indata[:,:] = indata[:,:] - np.mean(indata)
    sigma = np.std(indata)
    data_lim = scale_factor*sigma
    for i in range(ntheta):
        for j in range(nphi):
            if (indata[i,j] > data_lim):
                indata[i,j] = data_lim
            if (indata[i,j] < -data_lim):
                indata[i,j] = -data_lim
            indata[i,j] = indata[i,j]/data_lim

save_rendering = False  # Set to True to save an image (you lose the 3-D interactive window)
figure_file    = 'shell_rendering.jpg'  # The file we save to
shellfile      = 'r4_2_shell'  # A Rayleigh Shell_Slice file (assumed to lie within current working directory)
vis_option     = 2

if (vis_option == 1):
    # Radial component of velocity gives the elevation
    # Phi component velocity gives the color 
    hillsize = 0.002             # hillsize sets how much the surface of the globe may vary
    elevation_spec = [0,1,0]     # grab time index 0, quantity code 1, radial index 0
    elevation_scale_factor = 1.5 # scale_factor (used in tandem with stddev -- see scale_data above)

    color_spec = [0,3,0] #grab time index 0, quantity code 3, radial index 0
    color_scale_factor = 2.5
if (vis_option == 2):
    # Radial component of velocity gives the elevation and the color
    hillsize = 0.003  
    elevation_spec = [0,1,0]
    elevation_scale_factor = 1.5

    color_spec = [0,1,0]
    color_scale_factor = 1.5


shell1 = ShellSlice(shellfile, path='./',slice_spec = elevation_spec)
elevation_data = shell1.vals[:,:].reshape(shell1.nphi,shell1.ntheta)
elevation_data = np.transpose(elevation_data)

shell2 = ShellSlice(shellfile, path='./',slice_spec = color_spec)
color_data = shell2.vals[:,:].reshape(shell1.nphi,shell1.ntheta)
color_data = np.transpose(color_data)

nphi = shell1.nphi
ntheta = shell2.ntheta


# Scale both the elevation and the color data as desired
# If rendering something with a spherical mean, set no_mean = True 
# in the call to scale_data (as with temperature/entropy)
scale_data(elevation_data, scale_factor = elevation_scale_factor)
scale_data(color_data, scale_factor = color_scale_factor)



# Create a spherical grid
rmean = 0.3     # Mean radius of the sphere in Mayavi's units (leave this alone at first)
pi = np.pi
cos = np.cos
sin = np.sin
phi, theta = np.mgrid[0:pi:ntheta*1j, 0:2 * pi:nphi*1j]


# Add elevation to our globe (perturbations about rmean)
r = hillsize*elevation_data+rmean

# Recast in terms of x,y,z
x = r * sin(phi) * cos(theta)
y = r * sin(phi) * sin(theta)
z = r * cos(phi)

# Set up the window
winsize = 512
winx = winsize
winy = winsize
mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(winx, winy))
mlab.clf()

# Render topology using elevation_data (x,y,z) and color using color_data
mlab.mesh(x , y , z, scalars=-color_data, colormap='RdYlBu')


# Set the scene parameters and render
azimuth_angle = 0
elevation_angle = 70
view_distance = 1.3
mlab.view(azimuth = azimuth_angle, elevation = elevation_angle, distance = view_distance)
if (save_rendering):
    mlab.savefig(figure_file) 
else:
    mlab.show() 



