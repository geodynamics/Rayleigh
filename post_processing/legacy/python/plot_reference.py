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
#    Plot the reference state
#    This routine reads in a ReferenceState file
#    It plots the density, temperature, and density scaleheight.
#
#    The Attributes of the ReferenceState Structure are:
#    ----------------------------------
#    self.n_r         : number of radial points
#    self.radius      : radial coordinates
#    self.density     : density
#    self.dlnrho      : logarithmic derivative of density
#    self.d2lnrho     : d_by_dr of dlnrho
#    self.pressure    : pressure
#    self.temperature : temperature
#    self.dlnt        : logarithmic derivative of temperature
#    self.dsdr        : entropy gradient (radial)
#    self.entropy     : entropy
#    self.gravity     : gravity
#    --------------------------------------
from diagnostic_reading import ReferenceState
import matplotlib.pyplot as plt
ref = ReferenceState() #looks for ./reference by default
saveplot = False
savefile = 'reference.pdf'
plt.figure(1)

plt.subplot(221)
plt.plot(ref.radius/ref.radius[1],ref.density)
plt.xlabel('Radius (r/r_outer)')
plt.ylabel('Density'+r' g cm$^{-3}$')

plt.subplot(222)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.plot(ref.radius/ref.radius[1],ref.temperature)
plt.xlabel('Radius (r/r_outer)')
plt.ylabel('Temperature (K)')

plt.subplot(223)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.plot(ref.radius/ref.radius[1],ref.gravity)
plt.xlabel('Radius (r/r_outer)')
plt.ylabel('Gravity '+r'(cm s$^{-2}$)')

plt.subplot(224)
plt.plot(ref.radius/ref.radius[1],-1.0/ref.dlnrho)
plt.xlabel('Radius (r/r_outer)')
plt.ylabel('Density Scale Height (cm)')
plt.tight_layout()
if (saveplot):
    plt.savefig(savefile)
else:
    plt.show()
