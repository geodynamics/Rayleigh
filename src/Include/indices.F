!
!  Copyright (C) 2018 by the authors of the RAYLEIGH code.
!
!  This file is part of RAYLEIGH.
!
!  RAYLEIGH is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 3, or (at your option)
!  any later version.
!
!  RAYLEIGH is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with RAYLEIGH; see the file LICENSE.  If not see
!  <http://www.gnu.org/licenses/>.
!

#define RHSP wsp%p3b
#define DO_IDX Do t = my_theta%min, my_theta%max;    Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define END_DO enddo; enddo; enddo
#define IDX k,r,t
#define FIELDSP wsp%p3a
#define IDXM r,t
#define DO_IDX2 Do mp = my_mp%min, my_mp%max; m = m_values(mp); Do imi = 1, 2; Do r = my_r%min, my_r%max
#define IDX2 m:l_max,r,imi
#define SBUFFA wsp%s2a(mp)%data
#define SBUFFB wsp%s2b(mp)%data
#define ASBUFFA cobuffer%s2a(mp)%data
#define DO_PSI Do t = my_theta%min, my_theta%max;    Do r = my_r%min, my_r%max ;Do k = 1, n_phi
#define DO_PSI2 Do t = my_theta%min, my_theta%max;    Do r = my_r%min, my_r%max
#define END_DO2 enddo; enddo
#define PSI k,r,t
#define PSI2 r,t
#define DDBUFF d2buffer%p3a

