
import numpy as np
import scipy.special as scisp
import inspect
import math

def compute_plms(lm_max, costheta):
  """
  Compute the normalized spherical harmonics up to degree and order `lm_max` and colocation points `costheta`.
  """
  # FIXME: this will blow up rapidly because of its naive use of math.factorial
  plms = np.zeros((lm_max+1, lm_max+1, len(costheta)))
  for l in range(lm_max+1):
    for m in range(l+1):
      norm = np.sqrt(((2.*l+1.)/(4.*np.pi))*(float(math.factorial(l-m))/float(math.factorial(l+m))))
      plms[m,l,:] = norm*scisp.lpmv(m,l,costheta)
  return plms

def compute_gamma(n_r):
  """
  Compute the Chebyshev colocation points (zero-crossing of degree n_r).
  """
  return np.asarray([np.pi*(2*k + 1)/2./n_r for k in range(n_r)])

def compute_tns(n_max, gamma):
  """
  Compute Chebyshev polynomials up to degree `n_max` at colocation points `gamma`.
  """
  return np.asarray([np.cos(j*gamma) for j in range(n_max+1)])

def dealias_m2g(m_max):
  """
  Return the dealiased number of grid points.
  """
  g = int(3*(m_max + 1)/2)
  return g

def dealias_g2m(g):
  """
  Return the dealiased maximum order/degree.
  """
  m_max = int((2*g - 1)/3)
  if m_max < 0: m_max = 0
  return m_max

class SpectralInput(object):
  """
  Rayleigh class describing and writing generic spectral input for boundary/initial conditions.
  """

  version = 1
  lm_max  = None
  n_max   = None
  indices = None
  coeffs  = None

  def __init__(self, n_theta=None, n_r=None, \
               lm_max=None, n_max=None):
    """
    Optional arguments:
      - `n_theta` : Number of co-latitudinal grid points.
                    Ignored if lm_max is supplied.
      - `n_r`     : Number of radial grid point.
                    Ignored if n_max is supplied.
      - `lm_max`  : Maximum Legendre order and degree.
                    If supplied a multidimensional coefficient array is allocated, otherwise a sparse list structure will be used.

      - `n_max`   : Maximum Chebyshev polynomial degree.
                    Ignored if `lm_max` is not supplied.  
                    Assumed to be 0 if not supplied.
    """
    self.lm_max = lm_max
    self.n_max = n_max
    if self.lm_max is None and n_theta is None:
      # no dimensions provided so just set up a 1-D sparse storage to add modes to
      self.indices = []
      self.coeffs = np.zeros((0,), dtype='complex')
    else:
      # dimensions provided, set up a (n,l,m) array of coefficients
      if self.lm_max is None: self.lm_max = dealias_g2m(n_theta)
      if self.n_max is None:
        if n_r is None: 
          self.n_max = 0
        else:
          self.n_max = dealias_g2m(n_r)
      self.coeffs = np.zeros((self.n_max+1, self.lm_max+1, self.lm_max+1), dtype='complex')

  def add_mode(self, coeff, n=None, l=None, m=None, mode='add'):
    """
    Add a mode coefficient to the class.

    The coefficient (`coeff`) may be scalar, a list or 1-D array, or a multi-dimensional array (indexed by [l, m] in 2-D or [n, l,
    m] in 3-D).  
    If `coeff` is scalar or a list/1-D array, corresponding lists of indices for `l` and `m` must be provided. If `n` is not provided it is 
    assumed to be 0.  
    If `coeff` is a multi-dimensional array, `l`, `m` and `n` are ignored and the entries are assumed to be indexed by [l, m]
    (assiuming n=0) or [n, l, m] in increasing order from 0.

    Optional arguments:
      - `n`    : A list of Chebyshev n modes corresponding to the values in `coeff`.
                 Assumed to be 0 if not supplied. Ignored if `coeff` is 2- or 3-D.
      - `l`    : A list of Legendre function degrees corresponding to the values in `coeff`.
                 Ignored if `coeff` is 2- or 3-D.
      - `m`    : A list of Legendre function orders corresponding to the values in `coeff`.
                 Ignored if `coeff` is 2- or 3-D.
      - `mode` : How coefficients are added, either "replace", which overwrites, or "add", which adds, to existing entries.
    """

    if np.ndim(coeff) > 3:
      raise Exception("Don't know how to deal with coefficients of specified dimensions (({})).".format(",".join(map(str,coeff.shape))))

    def check_args(coeff, l, m, n):
      """Check arguments are consistent with scalar or list addition."""
      # Check the necessary arguments have been supplied
      if l is None or m is None: raise Exception("Specify l and m (n will be assumed to be 0 if not supplied).")
      # Elevate scalars to lists/arrays
      if np.isscalar(l): l = [l]
      if np.isscalar(m): m = [m]
      if n is None:
        n = [0]*len(l)
      else:
        if np.isscalar(n): n = [n]
      if np.isscalar(coeff): coeff = [coeff]
      coeff = np.asarray(coeff, dtype='complex')
      # Check lengths and indices
      if not len(l) == len(m) == len(n) == len(coeff): raise Exception("Length of l, m, n (if supplied) and coeff lists must match.")
      for i in range(len(l)):
        if n[i] < 0 or l[i] < 0 or m[i] < 0: raise Exception("All indicies must be non-negative.")
        if m[i] > l[i]: raise Exception("At i = {0}, m[i] ({1}) > l[i] ({2}), which is not allowed".format(i, m[i], l[i]))
      return coeff, l, m, n

    def check_dims(coeff):
      """Check dimensions are consistent with multi-dimensional addition."""
      # we assume we've been given a full array of coefficients
      # and that if ndim is 2, n dependence has been dropped
      # i.e. l, m and n are ignored
      if not (l is None and m is None and n is None): 
        raise Warning("l, m and n arguments will be ignored for rank > 1 arrays (assumed to be ordered).")
      coeff = np.asarray(coeff, dtype='complex')
      # Get dimensions of coeff
      l_maxp1, m_maxp1 = coeff.shape[-2:]
      # Fix 2-D arrays to be 3-D
      n_maxp1 = 1
      if np.ndim(coeff) == 3: 
        n_maxp1 = coeff.shape[0]
      else:
        coeff = coeff.reshape((1,)+coeff.shape)
      return coeff, l_maxp1, m_maxp1, n_maxp1
      
    # first, the case where lm_max hasn't been supplied to the class so arrays are 1-D
    if np.ndim(self.coeffs) == 1:
      # lists to store new entries
      newi = []
      newc = []
      # multi-dimensional coeffs
      if np.ndim(coeff) > 1:
        # we assume we've been given a full array of coefficients
        # and that if ndim is 2, n dependence has been dropped
        # i.e. l, m and n are ignored
        coeff, l_maxp1, m_maxp1, n_maxp1 = check_dims(coeff)
        for mj in range(m_maxp1):
          for lj in range(mj, l_maxp1):
            for nj in range(n_maxp1):
              nlm = (nj, lj, mj)
              i = None
              try:
                i = self.indices.index(nlm)
              except ValueError:
                pass
              if i is None:
                newi.append(nlm)
                newc.append(coeff[nlm])
              else:
                if mode == 'replace':
                  self.coeffs[i] = coeff[nlm]
                elif mode == 'add':
                  self.coeffs[i] += coeff[nlm]
                else:
                  raise Exception("Unknown addition mode ({}).".format(mode,))
      # scalar and 1-D lists/arrays
      else:
        # use l, m and n to locate coeffs in array
        coeff, l, m, n = check_args(coeff, l, m, n)
        for j in range(len(l)):
          nlm = (n[j], l[j], m[j])
          i = None
          try:
            i = self.indices.index(nlm)
          except ValueError:
            pass
          if i is None:
            newi.append(nlm)
            newc.append(coeff[j])
          else:
            if mode == 'replace':
              self.coeffs[i] = coeff[j]
            elif mode == 'add':
              self.coeffs[i] += coeff[j]
            else:
              raise Exception("Unknown addition mode ({}).".format(mode,))
      # if there are new entries to add, add them now (doesn't matter what the mode is)
      if len(newc) > 0:
        self.indices += newi
        self.coeffs = np.append(self.coeffs, np.asarray(newc, dtype='complex'))
    else:
      # multi-dimensional coeffs
      if np.ndim(coeff) > 1:
        coeff, l_maxp1, m_maxp1, n_maxp1 = check_dims(coeff)
        if n_maxp1 > self.coeffs.shape[0] \
           or l_maxp1 > self.coeffs.shape[1] \
           or m_maxp1 > self.coeffs.shape[2]:
          raise Exception("Number of coefficients exceed dimensions set by lm_max and n_max.")
        if mode == 'replace':
          self.coeffs[:n_maxp1, :l_maxp1, :m_maxp1] = coeff
        elif mode == 'add':
          self.coeffs[:n_maxp1, :l_maxp1, :m_maxp1] += coeff
      # scalar and 1-D lists/arrays
      else:
        # use l, m and n to locate coeffs in array
        coeff, l, m, n = check_args(coeff, l, m, n)
        for j in range(len(l)):
          if l[j] > self.lm_max or m[j] > self.lm_max or n[j] > self.n_max:
            raise Exception("Specified index exceeds lm_max or n_max, (n,l,m)=({},{},{})".format(str(n[j]), str(l[j]), str(m[j])))
          if mode == 'replace':
            self.coeffs[n[j], l[j], m[j]] = coeff[j]
          elif mode == 'add':
            self.coeffs[n[j], l[j], m[j]] += coeff[j]
          else:
            raise Exception("Unknown addition mode ({}).".format(mode,))

  def transform_from_rtp_function(self, func, func_kwargs={}, \
                                  n_theta=None, n_phi=None, n_r=None, \
                                  rmin=None, rmax=None, shell_depth=None, aspect_ratio=None, \
                                  mode='replace'):
    """
    Transform the provided function of spherical coordinates (radius, theta, phi) into Chebyshev/spherical harmonics (n, l, m).
    
    The provided function (`func(theta, phi, radius, ...)`) should accept any combination (including none) of the input parameters:
      - `theta`  : co-latitude
      - `phi`    : longitude
      - `radius` : radius
    These arguments must be named following this scheme.

    Optional arguments:
      - `func_kwargs`  : Supply a dictionary of additional keyword arguments to `func`.
                         Cannot contain keys for `theta`, `phi`, or `radius`.
      - `n_theta`      : Specify the number of co-latitudinal grid points.  
                         Set based on `lm_max` if not provided and `func` is a function of `theta`.
                         Ignored if `func` is not a function of `theta`.
      - `n_phi`        : Specify the number of longitudinal grid points.
                         Set based on `n_theta` if not provided and `func` is a function of `phi`.
                         Ignored if `func` is not a function of `phi`.
      - `n_r`          : Specify the number of radial grid points.
                         Set based on `n_max` if not provided and `func` is a function of `radius`.
                         Ignored if `func` is not a function of `radius`.
      - `rmin`         : Specify the minimum radius.
                         Ignored if `func` is not a function of `radius`.
      - `rmax`         : Specify the maximum radius.
                         Ignored if `func` is not a function of `radius`.
      - `shell_depth`  : Specify the shell depth. 
                         Ignored if `rmax` and either `rmin` or `aspect_ratio` are supplied.
      - `aspect_ratio` : Specify the aspect ratio.  
                         Ignored if  `rmax` and `rmin` are supplied.
      - `mode`         : How coefficients are added, either "replace", which overwrites, 
                         or "add", which adds, to existing entries
                         of the same order and degree.
    """
    func_param = inspect.signature(func).parameters
    # check what our function is a function of
    func_theta = False
    if 'theta' in func_param: func_theta = True
    func_phi = False
    if 'phi' in func_param: func_phi = True
    func_radius = False
    if 'radius' in func_param: func_radius = True

    # sanity check that nothing we want to set has been set by the user
    for k in ('theta', 'phi', 'radius'):
      if k in func_kwargs:
        raise Exception("'{}' supplied in func_kwargs.  You probably didn't mean to do this.".format((k,)))

    # some sanity checks that we have enough info for the radial domain
    if func_radius:
      if rmax is None:
        if shell_depth is None or aspect_ratio is None:
          raise Exception("Must supply shell_depth and aspect_ratio if rmax is not set.")
        rmax = shell_depth/(1.-aspect_ratio)
      if rmin is None:
        if aspect_ratio is None:
          raise Exception("Must supply aspect_ratio if rmin is not set.")
        rmin = rmax*aspect_ratio
    else:
      rmax = None
      rmin = None

    # work out the number of points...
    lm_max = self.lm_max
    if func_theta or func_phi:
      # theta (co-latitude)
      if n_theta is None:
        if lm_max is None:
          raise Exception("Cannot evaluate n_theta (or n_phi) unless lm_max is set.  Supply n_theta or set lm_max.")
        n_theta = dealias_m2g(lm_max)
      # phi (longitude)
      if func_phi and n_phi is None: n_phi = n_theta*2
      if lm_max is None: lm_max = dealias_g2m(n_theta)
      if not func_theta: n_theta = 1
    else:
      n_theta = 1
      n_phi = 1
      lm_max = 0

    n_max = self.n_max
    if func_radius:
      # radius (if we care about it)
      if n_r is None:
        if n_max is None: 
          raise Exception("Cannot evaluate n_r unless n_max is set.  Supply n_r or set n_max.")
        n_r = dealias_m2g(n_max)
      if n_max is None: n_max = dealias_g2m(n_r)
    else:
      n_r = 1
      n_max = 0
    
    # get Legendre Gauss integration points and weights
    # (we do this regardless of whether func_theta is true)
    costheta, legweights = np.polynomial.legendre.leggauss(n_theta)
    # set up Chebyshev integration points
    # (we do this regardless of whether func_radius is true)
    gamma = compute_gamma(n_r)

    theta, phi, rx, radius = [None]*4
    if func_theta:
      # set up co-latitudinal colocation points
      theta = np.arccos(costheta)
    if func_phi:
      # set up longitudinal colocation points
      phi = np.linspace(0, 2*np.pi, n_phi+1)[:-1]
    if func_radius:
      # set up radial colocation points
      rx = np.cos(gamma)
      # rescale to physical space
      radius = (rx-rx[-1])*(rmax-rmin)/(rx[0]-rx[-1]) + rmin

    # evaluate the function
    try:
      # set up a (2/3D) grid of points to evaluate the function at
      # hopefully most efficient if function is vectorized
      coord = [None]*3
      if func_theta:  coord[0] = theta
      if func_phi:    coord[1] = phi
      if func_radius: coord[2] = radius
      thetag, phig, radiusg = np.meshgrid(*coord)
      coordg = {}
      if func_theta:  coordg['theta']  = thetag
      if func_phi:    coordg['phi']    = phig
      if func_radius: coordg['radius'] = radiusg
      data_rtp = func(**{**coordg, **func_kwargs}).transpose()
    except:
      # more annoying but maybe more forgiving...
      data_rtp = np.zeros((n_r, n_theta, n_phi))
      coordg = {}
      for t in range(n_theta):
        if func_theta: coordg['theta'] = theta[t]
        for p in range(n_phi):
          if func_phi: coordg['phi'] = phi[p]
          for r in range(n_r):
            if func_radius: coordg['radius'] = radius[r]
            data_rtp[r,t,p] = func(**{**coordg, **func_kwargs})

    self.transform_from_rtp_data(data_rtp, \
                                 gamma=gamma, costheta=costheta, weights=legweights, \
                                 mode=mode)

  def transform_from_rtp_data(self, data_rtp, \
                              gamma=None, costheta=None, weights=None, \
                              mode='replace'):
    """
    Transform the provided array in spherical coordinates [radius, theta, phi] into Chebyshev/spherical harmonics (n, l, m).
    
    The provided array (`data_rtp`) must be indexed in [radius, theta, phi] in that order (radius, co-latitude, longitude).  
    It is assumed that the data is
    distributed in a structured grid in a spherical shell following an appropriate distribution of points: Chebyshev zero points in
    radius, Legendre-Gauss quadrature points in theta, and evenly spaced in phi.  Optional arguments may be supplied to describe the
    grid and integration weights.

    Optional arguments:
      - `gamma`    : Chebyshev colocation points in [-1,1].  Dimension must match data_rtp.shape[0].
      - `costheta` : Colocation points in [-1,1] (cosine of co-latitude).  Dimension must match data_rtp.shape[1].
      - `weights`  : Integration weights.  Dimension must match data_rtp.shape[1].
      - `mode`     : How coefficients are added, either "replace", which overwrites, or "add", which adds, to existing entries.
    """
    
    n_r = data_rtp.shape[0]
    n_theta = data_rtp.shape[1]
    n_phi = data_rtp.shape[2]

    # work out lm_max & n_max
    lm_max = self.lm_max
    if lm_max is None: lm_max = dealias_g2m(n_theta)
    n_max = self.n_max
    if n_max is None: n_max = dealias_g2m(n_r)

    # get Legendre Gauss integration points and weights
    if costheta is None or weights is None:
      if costheta is not None:
        raise Exception("Supplied weights but not costheta.")
      if weights is not None:
        raise Exception("Supplied costheta but not weights.")
      costheta, weights = np.polynomial.legendre.leggauss(n_theta)
    else:
      if len(costheta) != n_theta or len(weights) != n_theta:
        raise Exception("If supplied, length of costheta and weights must match input n_theta (data_rtp.shape[1]).")
     
    # set up Chebyshev integration points
    if gamma is None:
      gamma = compute_gamma(n_r)
    else:
      if len(gamma) != n_r:
        raise Exception("If supplied, length of gamma must match input n_r (data_rtp.shape[0]).")

    # Fourier transform: phi -> m
    data_rtm = np.fft.rfft(data_rtp)/n_phi
    data_rtm[:,:,1:] = 2.*data_rtm[:,:,1:]

    # compute the spherical harmonics
    plms = compute_plms(lm_max, costheta)
    data_rlm = np.zeros((n_r, lm_max+1, lm_max+1), dtype='complex')
    # Legendre transform: theta -> l
    for l in range(lm_max+1):
      for m in range(l+1):
        for r in range(n_r):
          data_rlm[r, l, m] = sum(2.*np.pi*data_rtm[r, :, m]*plms[m, l, :]*weights)

    # compute the Chebyshev polynomials
    tns = compute_tns(n_max, gamma)
    data_nlm = np.zeros((n_max+1, lm_max+1, lm_max+1), dtype='complex')
    # Chebyshev transform: radius -> n
    for l in range(lm_max+1):
      for m in range(l+1):
        for n in range(n_max+1):
          data_nlm[n, l, m] = (2./n_r)*sum(data_rlm[:, l, m]*tns[n, :])
    
    self.add_mode(data_nlm, mode=mode)


  def inverse_transform(self, n_theta=None, n_phi=None, n_r=None):
    """
    Inverse transform back to (r, theta, phi) space from the stored [n, l, m] coefficients.
    """
    if np.ndim(self.coeffs)==3:
      n_max = self.n_max
      lm_max = self.lm_max
      data_nlm = self.coeffs
    else:
      n,l,m = zip(*self.indices)
      n_max = max(n)
      lm_max = max(l)
      m_max = max(m)
      if m_max > lm_max: raise Exception("Found m_max > lm_max.  This should not happen.")
      data_nlm = np.zeros((n_max+1, lm_max+1, lm_max+1), dtype='complex')
      # FIXME: not the most efficient way of handling this but makes code below easier
      # (with lots of potential multiplying by zeros below)
      for i, nlm in enumerate(self.indices):
        data_nlm[nlm] = self.coeffs[i]

    # work out the number of points...
    # theta (co-latitude)
    if n_theta is None: n_theta = dealias_m2g(lm_max)
    # phi (longitude)
    if n_phi is None: n_phi = n_theta*2
    # radius
    if n_r is None: n_r = dealias_m2g(n_max)
    
    # get Legendre Gauss integration points and weights and compute the harmonics
    costheta, legweights = np.polynomial.legendre.leggauss(n_theta)
    plms = compute_plms(lm_max, costheta)
    phis = np.linspace(0, 2*np.pi, n_phi+1)[:-1]
    # set up Chebyshev integration points and compute the polynomials
    gamma = compute_gamma(n_r)
    tns = compute_tns(n_max, gamma)
    
    # Perform the inverse transform
    data_rtp = np.zeros((n_r, n_theta, n_phi))
    for r in range(n_r):
      for l in range(lm_max+1):
        for m in range(l+1):
          data_rlm = sum(tns[:,r]*data_nlm[:,l,m]) - 0.5*data_nlm[0,l,m]
          for p, phi in enumerate(phis):
            data_rtp[r,:,p] += (data_rlm.real*np.cos(m*phi) - data_rlm.imag*np.sin(m*phi))*plms[m,l,:]

    return data_rtp

  def sort(self):
    """
    Order the coefficients in column-major ordering, i.e. given indices (n,l,m), n varies fastest then l then m.
    
    Only affects sparsely stored coefficients as densely stored values are automatically ordered.
    """
    if np.ndim(self.coeffs)==1:
      order = [j[0] for j in sorted(enumerate(self.indices), key=lambda x: x[1][::-1])]
      self.coeffs = self.coeffs[order]
      self.indices = [self.indices[i] for i in order]

  def write(self, filename, byteswap=False):
    """
    Write spectral coefficients to file `filename`.
    """
    # FIXME: byteswap currently does nothing
    fd = open(filename, 'wb')

    if np.ndim(self.coeffs)==1:
      # sparsely stored coefficients (requires larger header)
      self.sort()
      header = np.ndarray((4+3*len(self.indices)),dtype='int32')
      header[0] = 314          # endian tag
      header[1] = self.version # file version (hard-coded above)
      header[2] = 0            # mode
      header[3] = len(self.indices)
      header[4:] = [i for nlm in zip(*self.indices) for i in nlm]
      header.tofile(fd)
      self.coeffs.real.tofile(fd)
      self.coeffs.imag.tofile(fd)
    else:
      flatcoeffs = np.asarray([self.coeffs[n,l,m] for m in range(self.lm_max+1) \
                                                  for l in range(m,self.lm_max+1) \
                                                  for n in range(self.n_max)])
      header = np.ndarray((5), dtype='int32')
      header[0] = 314          # endian tag
      header[1] = self.version # file version (hard-coded above)
      header[2] = 1            # mode
      header[3] = self.n_max
      header[4] = self.lm_max
      header.tofile(fd)
      flatcoeffs.real.tofile(fd)
      flatcoeffs.imag.tofile(fd)

    fd.close()
    



