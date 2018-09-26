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
#  in 2-D orthographic projection.  You will need to have Basemap installed.
#  This example was tested using BaseMap 1.0.7.
#  You may find BaseMap here:
#  http://matplotlib.org/basemap/users/download.html



from mpl_toolkits.basemap import Basemap, addcyclic
import numpy as np
import matplotlib.pyplot as plt
from diagnostic_reading import ShellSlice


#Specify which data to plot
shellfile = 'r4_2_shell'
rec_spec = [0,3,1]  # grab time index 0, quantity code 3, and radial index 1  


save_figure = False # Set to true to save to a figure_file (below)
figure_file = 'shell_slice_basemap.png'

#The resolution of the png file or the view window in pixels
xpixels = 1024  
ypixels = 1024

# Initialize the projection.
# lon_0, lat_0 are the center point of the projection.
# resolution = 'l' means use low resolution coastlines.
m = Basemap(projection='ortho',lon_0=-20,lat_0=35,resolution='l')


#Read in the data
a = ShellSlice(shellfile, path='./',slice_spec = rec_spec)
data = a.vals[:,:,0,0,0].reshape(a.nphi,a.ntheta)
data = np.transpose(data)
maxabs = 3*np.std(data)  # saturate data at +/- 3 sigma
dlon = 360.0/(a.nphi)
dlat = 180.0/a.ntheta

#Scale the data
for i in range(a.ntheta):
    for j in range(a.nphi):
        if (data[i,j] > maxabs):
            data[i,j] = maxabs
        if (data[i,j] < -maxabs):
            data[i,j] = -maxabs
        data[i,j] = data[i,j]/maxabs 


#Generate 1-D grids of latitude and longitude
lons = np.zeros(a.nphi)
for i in range(a.nphi):
    lons[i] = dlon*i-180.0

lats = np.zeros(a.ntheta)
for i in range(a.ntheta):
    lats[i] = 90.0-np.arccos(a.costheta[i])*180.0/np.pi
    lats[i] = i*dlat-90


# Convert to 2-D grids
llons, llats = np.meshgrid(lons, lats)

# Get x-y projection points on the plane
x, y = m(llons, llats)

# Do some interpolation 
# This is necessary for moderately-sized and large shell slices (ell_max >= 255)
nx = int((m.xmax-m.xmin)/2000.)+1; ny = int((m.ymax-m.ymin)/2000.)+1
nx = 1024
ny = 1024
print 'interpolating'
topodat,x,y =\
m.transform_scalar(data,lons,lats,nx,ny,returnxy=True,masked=True,order=1)
print 'done'

#Initialize the window
plt.figure(figsize=(xpixels/100.0, ypixels/100.0))

#View the data
my_cmap = plt.cm.RdYlBu_r
m.pcolormesh(x,y,topodat,cmap=my_cmap)

# draw parallels and meridians.
m.drawparallels(np.arange(-90.,120.,30.))
m.drawmeridians(np.arange(0.,420.,60.))
m.drawmapboundary(fill_color='aqua')

if (save_figure):
    plt.savefig(figure_file)
else:
    plt.show()

