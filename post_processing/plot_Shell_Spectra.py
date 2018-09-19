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

# VIII.  Spherical Harmonic Spectra
# ===================
# 
# **Summary:**    Spherical Harmonic Spectra sampled at discrete radii. 
# 
# **Subdirectory:**  Shell_Spectra
# 
# **main_input prefix:** shellspectra
# 
# **Python Classes:** 
# 
# * Shell_Spectra :  Complete data structure associated with Shell_Spectra outputs.
# * PowerSpectrum :  Reduced data structure -- contains power spectrum of velocity and/or magnetic fields only.
# 
# **Additional Namelist Variables:**  
# 
# * shellspectra_levels (indicial) : indices along radial grid at which to output spectra.
# 
# * shellspectra_levels_nrm (normalized) : normalized radial grid coordinates at which to output spectra.
# 
# 
# The shell-spectra output type allows us to examine the spherical harmonic decomposition of output variables at discrete radii.
# 
# Examining the *main_input* file, we see that the following output values have been denoted for the Shell Spectra (see *rayleigh_output_variables.pdf* for mathematical formulae):
# 
# 
# | Menu Code  | Description |
# |------------|-------------|
# | 1          | Radial Velocity |
# | 2          | Theta Velocity |
# | 3          | Phi Velocity  |
# 
# 
# 
# Spherical harmonic spectra can be read into Python using either the **Shell_Spectra** or **PowerSpectrum** classes.  
# 
# The **Shell_Spectra** class provides the full complex spectra, as a function of degree ell and azimuthal order m, for each specified output variable.   It possesses an attribute named *lpower* that contains the associated power for each variable, along with its m=0 contributions separated and removed.
# 
# The **Power_Spectrum** class can be used to read a Shell_Spectra file and quickly generate a velocity or magnetic-field power spectrum.   For this class to work correctly, your file must contain all three components of either the velocity or magnetic field.   Other variables are ignored (use Shell_Spectrum's lpower for those).
# 
# We illustrate how to use these two classes below.  As usual, we call the help() function to display the docstrings that describe the different data structures embodied by each class.

# In[23]:

import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy
from rayleigh_diagnostics import Shell_Spectra, Power_Spectrum
istring = '00040000'


tind = 0
rind = 0

vpower = Power_Spectrum(istring)

power = vpower.power

fig, ax = plt.subplots(nrows=3, figsize=(6,6))
ax[0].plot(power[:,rind,tind,0])
ax[0].set_xlabel(r'Degree $\ell$')
ax[0].set_title('Velocity Power (total)')


ax[1].plot(power[:,rind,tind,1])
ax[1].set_xlabel(r'Degree $\ell$')
ax[1].set_title('Velocity Power (m=0)')

ax[2].plot(power[:,rind,tind,2])
ax[2].set_xlabel(r'Degree $\ell$')
ax[2].set_title('Velocity Power ( total - {m=0} )')

plt.tight_layout()
savefile = 'Power_1D.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)

fig, ax = plt.subplots()
ss = Shell_Spectra(istring)

mmax = ss.mmax
lmax = ss.lmax
power_spectrum = numpy.zeros((lmax+1,mmax+1),dtype='float64')

for i in range(1,4):   # i takes on values 1,2,3
    qind=ss.lut[i]
    complex_spectrum = ss.vals[:,:,rind,qind,tind]
    power_spectrum = power_spectrum+numpy.real(complex_spectrum)**2 + numpy.imag(complex_spectrum)**2

power_spectrum = numpy.transpose(power_spectrum)

tiny = 1e-6
img=ax.imshow(numpy.log10(power_spectrum+tiny), origin='lower')
ax.set_ylabel('Azimuthal Wavenumber m')
ax.set_xlabel(r'Degree $\ell$')
ax.set_title('Velocity Power Spectrum')

#colorbar ...
cbar = plt.colorbar(img) # ,shrink=0.5, aspect = 15)
cbar.set_label('Log Power')
        
tick_locator = ticker.MaxNLocator(nbins=5)
cbar.locator = tick_locator
cbar.update_ticks()
cbar.ax.tick_params()   #font size for the ticks


savefile = 'Power_2D.pdf'
print('Saving figure to: ', savefile)
plt.savefig(savefile)
