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
#  Plots average Energy flux (normalized by 4pi r^2)
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
savefile = 'energy_flux.pdf'  #If pdf gives issues, try .png instead.

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
fpr = 4.0*numpy.pi*radius*radius  # four pi r^2
rnorm = radius/radius[0] # normalize radius so that upper boundary is r=1


#Next, we identify and grab the various fluxes we would like to plot
cf_index  = a.lut[122]  # Conductive energy flux
vhf_index = a.lut[121]  # The effective flux associated with volumetric heating
enf_index = a.lut[119]  # Enthalpy flux
vsf_index = a.lut[120]  # Viscous energy flux
kef_index = a.lut[117]  # Kinetic energy flux

cflux = data[:,0,cf_index] # Grab appropriate slices of the data array
hflux = data[:,0,vhf_index]
eflux = data[:,0,enf_index]
vflux = data[:,0,vsf_index]
kflux = data[:,0,kef_index]

tflux = cflux+eflux+vflux+kflux+hflux # total flux


#Create the plot

lw = 1.5

if (saveplot):
    plt.figure(1,figsize=(7.0, 5.0), dpi=300)
    plt.rcParams.update({'font.size': 12})
else:
    plt.figure(1,figsize=(7,5),dpi=100)
    plt.rcParams.update({'font.size': 14})

plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.plot(rnorm,tflux*fpr,label = 'F'+r'$_{total}$',linewidth=lw,color='black')
plt.plot(rnorm,eflux*fpr,label = 'F'+r'$_{enth}$',linewidth=lw)
plt.plot(rnorm,cflux*fpr,label = 'F'+r'$_{cond}$',linewidth=lw)
plt.plot(rnorm,vflux*fpr,label = 'F'+r'$_{visc}$',linewidth=lw)
plt.plot(rnorm,hflux*fpr,label = 'F'+r'$_{heat}$',linewidth=lw)
plt.plot(rnorm,kflux*fpr,label = 'F'+r'$_{KE}$',linewidth=lw)

plt.xlim([0.75,1])
plt.ylim([-2e33,5e33])
plt.ylabel(r'4$\pi r^2$'+'x Energy Flux '+r'( erg s$^{-1}$)')
plt.xlabel('Radius '+r'(r/r$_{o}$)')
legend = plt.legend(loc='lower left', shadow=True, ncol = 3) 
plt.tight_layout()


if (saveplot):
    plt.savefig(savefile)
else:
    plt.show()


