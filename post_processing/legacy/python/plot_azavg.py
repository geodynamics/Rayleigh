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

####################################################################################################
#
#  Azimuthal-Averages (AZ_Avgs) plotting example
#  Time-averages data from time steps 3 million 
#  through 3.3 million.
#  Plots Entropy, Differential Rotation, and Mass Flux
#
#  This example routine makes use of the AZAverage
#  data structure associated with the AZ_Avg output.
#  Upon initializing an AZAverage object, the 
#  object will contain the following attributes:
#    ----------------------------------
#    self.niter                                    : number of time steps
#    self.nq                                       : number of diagnostic quantities output
#    self.nr                                       : number of radial points
#    self.ntheta                                   : number of theta points
#    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
#    self.radius[0:nr-1]                           : radial grid
#    self.costheta[0:ntheta-1]                     : cos(theta grid)
#    self.sintheta[0:ntheta-1]                     : sin(theta grid)
#    self.vals[0:ntheta-1,0:nr-1,0:nq-1,0:niter-1] : The phi-averaged diagnostics 
#    self.iters[0:niter-1]                         : The time step numbers stored in this output file
#    self.time[0:niter-1]                          : The simulation time corresponding to each time step
#    self.version                                  : The version code for this particular output (internal use)
#    self.lut                                      : Lookup table for the different diagnostics output
#    ----------------------------------

from diagnostic_reading import AzAverage, build_file_list, TimeAvg_AZAverages
import matplotlib.pyplot as plt
import pylab as p 
import numpy as np
from azavg_util import *
#################################################################
# First define a few logical control flags
# Set saveplot to True to save to a .png file. 
# Set to False to view plots interactively on-screen.
saveplot = False
savefile = 'az_average.png'

# Once the data has been time-averaged, the result may be read in
#  without re-averaging by setting timeavg to False
timeavg = True


########################################################################
# Read in the data

# If time-averaging needs to be performed, we read in
# the indicated files from AZ_Avgs and a time-averaget to avg_file
avg_file = 'time_averaged_azavg'
if (timeavg):
    files = build_file_list(3000000,3300000,path='AZ_Avgs')
    TimeAvg_AZAverages(files,avg_file)

az = AzAverage(filename=avg_file,path='')


# time index to grab (this file only has one time since it is an average)
ind = 0 

#Find the indices associated with the quantities we want to plot
vphi_index    = az.lut[3]
entropy_index = az.lut[64]
rhovr_index   = az.lut[61]
rhovt_index   = az.lut[62]

#Grab quantities of interest from az.vals

n_r = az.nr
n_t = az.ntheta

vphi    = az.vals[:,:,az.lut[3],ind].reshape(n_t,n_r)
entropy = az.vals[:,:,az.lut[64],ind].reshape(n_t,n_r)

rhovr = az.vals[:,:,az.lut[61],ind].reshape(n_t,n_r)
rhovt = az.vals[:,:,az.lut[62],ind].reshape(n_t,n_r)

sintheta = az.sintheta
costheta = np.cos(np.arcsin(sintheta))
radius = az.radius

costheta[0:n_t/2] = -costheta[0:n_t/2] #sintheta is symmetric, so arcsin gives redundant values across equator
#costheta grid in Rayleigh runs from -1 to 1 (south pole to north pole..)

#Compute omega and subtract mean from entropy
Omega=np.zeros((n_t,n_r))

#Subtrace the ell=0 component from entropy at each radius
for i in range(n_r):
    entropy[:,i]=entropy[:,i] - np.mean(entropy[:,i])

#Convert v_phi to an Angular velocity
for i in range(n_r):
    Omega[:,i]=vphi[:,i]/(radius[i]*sintheta[:])
Omega = Omega*1.e9/(2.*np.pi) # s^-1 -> nHz

#Generate a streamfunction from rhov_r and rhov_theta
psi = streamfunction(rhovr,rhovt,radius,costheta,order=0)

#contours of mass flux are overplotted on the streamfunction PSI
rhovm = np.sqrt(rhovr**2+rhovt**2)*np.sign(psi)

##############################################
#   Here we set up the actual figure.
#   We do a single row of 3 images 
#   Spacing is default spacing set up by subplot

if (saveplot):
    f1 = p.figure(figsize=(5.5*3, 3*3), dpi=300)
    plt.rcParams.update({'font.size': 12})
else:
    f1 = p.figure(figsize=(5.5*3, 3*3), dpi=80)
    plt.rcParams.update({'font.size': 14})


#Entropy
tsize = 20
lsize = 15
ax1 = f1.add_subplot(1,3,1)
units = r'erg g$^{-1}$ K$^{-1}$'
plot_azav(f1,ax1,entropy,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 1.5, boundstype='rms', units=units, fontsize = lsize)
plt.title('Specific Entropy',fontsize=tsize)

#Differential Rotation
ax1 = f1.add_subplot(1,3,2)
units = 'nHz'
plot_azav(f1,ax1,Omega,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 0.05, units=units, fontsize = lsize)
plt.title(r'$\Omega-\Omega_0$',fontsize=tsize)

#Mass Flux
ax1 = f1.add_subplot(1,3,3)
units = r'g cm$^{-2}$ s$^{-1}$'
plot_azav(f1,ax1,psi,radius,costheta,sintheta,mycmap='RdYlBu_r',boundsfactor = 1.5, boundstype='rms', units=units, fontsize = lsize, underlay = rhovm)
plt.title('Mass Flux',fontsize = tsize)

if (saveplot):
    p.savefig(savefile)  
else:
    plt.show()

