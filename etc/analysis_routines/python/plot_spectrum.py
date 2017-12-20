########################################################################3
#
#   Plotting Example:  Shell_Spectra
#
#   We plot the velocity power spectrum from one
#   shell spectrum output using the PowerSpectrum
#   class from diagnostic_reading.py.  
#   
#   Ideally, one should loop over several spectra
#   and take a time-average, similar to the approach
#   taken in plot_energy_flux.py.  This routine plots
#   a SINGLE SNAPSHOT of the power.
#
#   When a PowerSpectrum object is initialized, 
#   only power is computed from the Shell_Spectra
#   file and saved.
#  
#   If we want to look at individual m-values,
#   we can use the ShellSpectra class instead.
#   See diagnostic_reading.py for a brief 
#   description of that data structure.
#
#   The power spectrum structure is simpler to
#   work with.  It contains the following attributes:
#    ----------------------------------
#    self.niter                                    : number of time steps
#    self.nr                                       : number of radii at which power spectra are available
#    self.lmax                                     : maximum spherical harmonic degree l
#    self.radius[0:nr-1]                           : radii of the shell slices output
#    self.inds[0:nr-1]                             : radial indices of the shell slices output
#    self.power[0:lmax,0:nr-1,0:niter-1,0:2]       : the velocity power spectrum.  The third
#                                                  : index indicates (0:total,1:m=0, 2:total-m=0 power)
#    self.mpower[0:lmax,0:nr-1,0:niter-1,0:2]      : the magnetic power spectrum
#    self.iters[0:niter-1]                         : The time step numbers stored in this output file
#    self.time[0:niter-1]                          : The simulation time corresponding to each time step
#    self.magnetic                                 : True if mpower exists
#    -------------------------------------
##################################

from diagnostic_reading import PowerSpectrum
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

# spec.power will contain the power spectrum
# this object can also be initialized with magnetic=True
# to generate a magnetic power spectrum stored 
# in spec.mpower
spec = PowerSpectrum('Shell_Spectra/03290000')

time_index = 0  #plot the first record
rad_index = 3 # plot the mid-depth radius (4 radius output in this case)


plt.figure(1)
lw = 1.5
plt.plot(spec.power[:,rad_index,time_index,0], label ='Total Power',linewidth=lw)
plt.plot(spec.power[:,rad_index,time_index,1], label ='m=0 Power',linewidth=lw)
plt.plot(spec.power[:,rad_index,time_index,2], label ='Convective Power ( Total - {m=0} )',linewidth=lw)
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

