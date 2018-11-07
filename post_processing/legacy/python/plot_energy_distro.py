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

###############################################
#
#  Shell-Averages (Shell_Avgs) plotting example
#  Reads in time steps 3 million through 3.3 million
#  Plots average KE and v' vs. radius
#
#  This example routine makes use of the ShellAverage
#  data structure associated with the Shell_Avg output.
#  Upon initializing a ShellAverage object, the 
#  object will contain the following attributes:
#
#    ----------------------------------
#    self.niter                         : number of time steps
#    self.nq                            : number of diagnostic quantities output
#    self.nr                            : number of radial points
#    self.qv[0:nq-1]                    : quantity codes for the diagnostics output
#    self.radius[0:nr-1]                : radial grid
#
#    self.vals[0:nr-1,0:3,0:nq-1,0:niter-1] : The spherically averaged diagnostics
#                                             0-3 refers to moments (index 0 is mean, index 3 is kurtosis)    
#    self.iters[0:niter-1]              : The time step numbers stored in this output file
#    self.time[0:niter-1]               : The simulation time corresponding to each time step
#    self.version                       : The version code for this particular output (internal use)
#    self.lut                           : Lookup table for the different diagnostics output
#   -------------------------------------

from diagnostic_reading import ShellAverage, build_file_list
import matplotlib.pyplot as plt
import numpy
# Set saveplot to True to save to a file. 
# Set to False to view plots interactively on-screen.
saveplot = False 
savefile = 'energy_distro.pdf'  #If pdf gives issues, try .png instead.

#Time-average the data generated between timesteps 3 million and 3.5 million
# Note that this averaging scheme assumes taht the radial grid didn't
# change during the course of the run. 
files = build_file_list(3000000,3300000,path='Shell_Avgs')
icount = 0.0
for i,f in enumerate(files):
    a = ShellAverage(f,path='')
    if (i == 0):
        data = numpy.zeros((a.nr,4,a.nq),dtype='float64')
    for j in range(a.niter):
        data[:,:,:] = data[:,:,:]+a.vals[:,:,:,j]
        icount = icount+1.0
data = data/icount  # The average is now stored in data
radius = a.radius
rnorm = radius/radius[0] # normalize radius so that upper boundary is r=1



#For plot 1, we plot the three components of the kinetic energy as a function of radius
rke_index = a.lut[126]  # KE associated with radial motion
tke_index = a.lut[127]  # KE associated with theta motion
pke_index = a.lut[128]  # KE associated with azimuthal motion

rke = data[:,0,rke_index] # Grab appropriate slices of the data array
tke = data[:,0,tke_index]
pke = data[:,0,pke_index]

#For plot 2, we plot the square root of the variance of each velocity component.
#In this case, the mean used in the variance is computed on spherical shells.
# Note that since we are grabbing the 2nd moment of the distribution, we
# grab index 1, not 0.  Skewness is stored in index 2, and Kurtosis in 3.

vr_index = a.lut[1]  # V_r
vt_index = a.lut[2]  # V_theta
vp_index = a.lut[3]  # V_phi

vr = numpy.sqrt(data[:,1,vr_index])*0.01 # Grab appropriate slices of the data array
vt = numpy.sqrt(data[:,1,vt_index])*0.01 # convert to m/s
vp = numpy.sqrt(data[:,1,vp_index])*0.01

lw = 1.5

if (saveplot):
    plt.figure(1,figsize=(7.5, 4.0), dpi=300)
    plt.rcParams.update({'font.size': 12})
else:
    plt.figure(1,figsize=(15,5),dpi=100)
    plt.rcParams.update({'font.size': 14})

plt.subplot(121)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.plot(rnorm,rke,label = 'KE'+r'$_{rad}$',linewidth=lw)
plt.plot(rnorm,tke,label = 'KE'+r'$_\theta$',linewidth=lw)
plt.plot(rnorm,pke,label = 'KE'+r'$_\phi$',linewidth=lw)
plt.xlim([0.75,1])
plt.ylabel('Energy Density '+r'( erg cm$^{-3}$)')
plt.xlabel('Radius '+r'(r/r$_{o}$)')
legend = plt.legend(loc='upper left', shadow=True, ncol = 1) 


plt.subplot(122)
plt.plot(rnorm,vr,label = "v'"+r'$_{r}$',linewidth=lw)
plt.plot(rnorm,vt,label = "v'"+r'$_\theta$',linewidth=lw)
plt.plot(rnorm,vp,label = "v'"+r'$_\phi$',linewidth=lw)
plt.xlim([0.75,1])
plt.ylabel('Speed '+r'( m s$^{-1}$)')
plt.xlabel('Radius '+r'(r/r$_{o}$)')
legend = plt.legend(loc='upper left', shadow=True, ncol = 1) 
plt.tight_layout()

if (saveplot):
    plt.savefig(savefile)
else:
    plt.show()


