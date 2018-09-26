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

########################################################################3
#
#   Plotting Example:  Shell_Spectra2
#
#   We plot the velocity power spectrum from one
#   shell spectrum output using the ShellSpectra
#   class from diagnostic_reading.py.  
#   
#   Ideally, one should loop over several spectra
#   and take a time-average, similar to the approach
#   taken in plot_energy_flux.py.  This routine plots
#   a SINGLE SNAPSHOT of the power.
#
#   One attribute of the ShellSpectra Object is lpower - that's what we use here
#   Lpower is the m-integrated power at each ell-value for all data in the ShellSpectra 
#   structure
#  
#    """Rayleigh Shell Spectrum Structure
#    ----------------------------------
#    self.niter                                    : number of time steps
#    self.nq                                       : number of diagnostic quantities output
#    self.nr                                       : number of shell slices output
#    self.nell                                     : number of ell values
#    self.nm                                       : number of m values
#    self.lmax                                     : maximum spherical harmonic degree l
#    self.mmax                                     : maximum spherical harmonic degree m
#    self.qv[0:nq-1]                               : quantity codes for the diagnostics output
#    self.radius[0:nr-1]                           : radii of the shell slices output
#    self.inds[0:nr-1]                             : radial indices of the shell slices output
#    self.vals[0:lmax,0:mmax,0:nr-1,0:nq-1,0:niter-1] 
#                                                  : The complex spectra of the shells output 
#    self.lpower[0:lmax,0:nr-1,0:nq-1,0:niter-1,3]    : The power as a function of ell, integrated over m
#                                                     :  index indicates (0:total,1:m=0, 2:total-m=0 power)
#    self.iters[0:niter-1]                         : The time step numbers stored in this output file
#    self.time[0:niter-1]                          : The simulation time corresponding to each time step
#    self.version                                  : The version code for this particular output (internal use)
#    self.lut                                      : Lookup table for the different diagnostics output
#    """
#    -------------------------------------
##################################

from diagnostic_reading import ShellSpectra
import matplotlib.pyplot as plt
import numpy as np

#Set savefig to True to save to savefile.  Interactive plot otherwise.
savefig = False
savefile = 'power_spectrum.pdf'
if (savefig):
    plt.figure(1,figsize=(7.5, 4.0), dpi=300)
    plt.rcParams.update({'font.size': 12})
else:
    pass
    #plt.figure(1,figsize=(15,5),dpi=100)
    #plt.rcParams.update({'font.size': 14})


spec = ShellSpectra('00002000')
vr_index = spec.lut[1]
vt_index = spec.lut[2]
vp_index = spec.lut[3]
test_index = spec.lut[194]
power = spec.lpower[:,:,vr_index,:,:]+ spec.lpower[:,:,vt_index,:,:] + spec.lpower[:,:,vp_index,:,:]

time_index = 0  #plot the first record
rad_index = 3 # plot the mid-depth radius (4 radius output in this case)


plt.figure(1)
lw = 1.5


plt.plot(power[:,rad_index,time_index,0], label ='Total Power',linewidth=lw)
plt.plot(power[:,rad_index,time_index,1], label ='m=0 Power',linewidth=lw)
plt.plot(power[:,rad_index,time_index,2], label ='Convective Power ( Total - {m=0} )',linewidth=lw)
plt.yscale('log')
plt.xscale('log')

plt.xlim(xmin=1,xmax=spec.lmax-1)
plt.xlabel('Spherical Harmonic Degree '+r'$\ell$')
plt.ylabel('Velocity Power ')
legend = plt.legend(loc='lower left', shadow=True, ncol = 2) 
plt.tight_layout()
if (savefig):
    plt.savefig(savefile)
else:
    plt.show()

