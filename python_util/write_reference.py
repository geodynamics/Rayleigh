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

# The routine write_reference in this module takes a filename and
# writes radial profiles of the reference state in binary form
# to that filename, such that Rayleigh can interpret that binary file
# specified through 
#
# reference_type=4 
# custom_reference_file = [filename]
#
# ...in the Reference_Namelist of main_input.

import numpy as np

def write_reference(filename, radius, rho, dlnrho, d2lnrho, pressure,\
        temperature, dlnT, dsdr, entropy, gravity):

    # first write the
    # two int32 binary numbers (unformatted): 314 and n_r. Then write
    # 10*n_r float64 numbers, representing (in order):
    # radius, density, dlnrho, d2lnrho, pressure, 
    # temperature dlnT, dsdr, entropy, gravity
    # n_r here should be >> n_r from Rayleigh so that
    # Rayleigh can interpolate your possibly fine-structured 
    # reference state onto its grid
    # the radius array may be linearly spaced; Rayleigh will interpolate

    f = open(filename, "wb")
    # may need to specify the data type for a successful read on
    # Rayleigh's end
    sigpi = np.array(314, dtype=np.int32)
    nr = np.array(len(radius), dtype=np.int32)
    f.write(sigpi.tobytes())
    f.write(nr.tobytes())
    f.write(radius.tobytes())
    f.write(rho.tobytes())
    f.write(dlnrho.tobytes())
    f.write(d2lnrho.tobytes())
    f.write(pressure.tobytes())
    f.write(temperature.tobytes())
    f.write(dlnT.tobytes())
    f.write(dsdr.tobytes())
    f.write(entropy.tobytes())
    f.write(gravity.tobytes())
    f.close()
