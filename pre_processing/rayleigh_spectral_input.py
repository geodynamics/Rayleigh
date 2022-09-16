#!/usr/bin/env python3

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

import numpy as np
import scipy.special as scisp
import funcsigs
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

def radial_extents(rmin=None, rmax=None, aspect_ratio=None, shell_depth=None):
  """
  Return rmin and rmax if they aren't set using aspect_ratio and shell_depth.
  """
  if rmax is None:
    if shell_depth is None or aspect_ratio is None:
      raise Exception("Must supply shell_depth and aspect_ratio if rmax is not set.")
    rmax = shell_depth/(1.-aspect_ratio)
  if rmin is None:
    if aspect_ratio is None:
      raise Exception("Must supply aspect_ratio if rmin is not set.")
    rmin = rmax*aspect_ratio
  return rmin, rmax

def swapwrite(vals,fd,byteswap=False):
  """
  Write to file handle fd with the option of byteswapping.
  """
  # this is a tidied up copy of post_processing/rayleigh_diagnostics.py
  # reproduced here to avoid dependency on post_processing and to tidy it up without conflict
  # FIXME: unify python routines to avoid this duplication
  valsout = vals
  if byteswap: valsout = vals.newbyteorder()
  valsout.tofile(fd)

def swapread(fd,dtype='float64',count=1,byteswap=False):
  """
  Read from file handle fd with the option of byteswapping.
  """
  # this is a tidied up copy of post_processing/rayleigh_diagnostics.py
  # reproduced here to avoid dependency on post_processing and to tidy it up without conflict
  # FIXME: unify python routines to avoid this duplication
  vals = np.fromfile(fd, dtype=dtype, count=count)
  if byteswap: vals.byteswap()
  return vals

def check_byteswap(fd):
  """
  Check endianess of file fd by reading an integer 314 and testing it matches.  
  Returns True if byteswapping is necessary.
  """
  # this is a tidied up copy of post_processing/rayleigh_diagnostics.py
  # reproduced here to avoid dependency on post_processing and to tidy it up without conflict
  # FIXME: unify python routines to avoid this duplication
  chk = np.fromfile(fd, dtype='int32', count=1)
  if (chk == 314): return False
  return True

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
                    Ignored if `lm_max` is supplied.
                    If this or `lm_max` is supplied a multidimensional coefficient array is allocated, otherwise a sparse list structure will be used.
      - `n_r`     : Number of radial grid point.
                    Ignored if `n_max` is supplied.
      - `lm_max`  : Maximum Legendre order and degree.
                    If this or `n_theta` is supplied a multidimensional coefficient array is allocated, otherwise a sparse list structure will be used.
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
    func_param = funcsigs.signature(func).parameters
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
      rmin, rmax = radial_extents(rmin, rmax, aspect_ratio, shell_depth)
    else:
      rmin = None
      rmax = None

    # work out the number of points...
    if func_theta or func_phi:
      # theta (co-latitude)
      if n_theta is None:
        if self.lm_max is None:
          raise Exception("Cannot evaluate n_theta (or n_phi) unless lm_max is set.  Supply n_theta or set lm_max.")
        n_theta = dealias_m2g(self.lm_max)
      # phi (longitude)
      if n_phi is None:
        if func_phi: 
          n_phi = n_theta*2
        else:
          n_phi = 1
      if not func_theta: n_theta = 1
    else:
      n_theta = 1
      n_phi = 1

    if func_radius:
      # radius (if we care about it)
      if n_r is None:
        if self.n_max is None: 
          raise Exception("Cannot evaluate n_r unless n_max is set.  Supply n_r or set n_max.")
        n_r = dealias_m2g(self.n_max)
    else:
      n_r = 1
    
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
      coordg.update(func_kwargs)
      if func_theta:  coordg['theta']  = thetag
      if func_phi:    coordg['phi']    = phig
      if func_radius: coordg['radius'] = radiusg
      data_rtp = func(**coordg).transpose()
    except:
      # more annoying but maybe more forgiving...
      data_rtp = np.zeros((n_r, n_theta, n_phi))
      coordg = {}
      coordg.update(func_kwargs)
      for t in range(n_theta):
        if func_theta: coordg['theta'] = theta[t]
        for p in range(n_phi):
          if func_phi: coordg['phi'] = phi[p]
          for r in range(n_r):
            if func_radius: coordg['radius'] = radius[r]
            data_rtp[r,t,p] = func(**coordg)

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

    # work out lm_max, l_max, m_max & n_max
    lm_max = self.lm_max
    if lm_max is None: lm_max = dealias_g2m(n_theta)
    l_max = min(lm_max, dealias_g2m(n_theta))
    m_max = min(lm_max, dealias_g2m(n_phi))
    n_max = self.n_max
    if n_max is None: n_max = dealias_g2m(n_r)
    n_max = min(n_max, dealias_g2m(n_r))

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
    data_rlm = np.zeros((n_r, l_max+1, m_max+1), dtype='complex')
    # Legendre transform: theta -> l
    for l in range(l_max+1):
      for m in range(min(l+1, m_max+1)):
        for r in range(n_r):
          data_rlm[r, l, m] = sum(2.*np.pi*data_rtm[r, :, m]*plms[m, l, :]*weights)

    # compute the Chebyshev polynomials
    tns = compute_tns(n_max, gamma)
    data_nlm = np.zeros((n_max+1, l_max+1, m_max+1), dtype='complex')
    # Chebyshev transform: radius -> n
    for l in range(l_max+1):
      for m in range(min(l+1, m_max+1)):
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

    Optional arguments:
      - byteswap: byte swap the data being written to switch endianness (default: False)
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
      swapwrite(header, fd, byteswap)
      swapwrite(self.coeffs.real, fd, byteswap)
      swapwrite(self.coeffs.imag, fd, byteswap)
    else:
      flatcoeffs = np.asarray([self.coeffs[n,l,m] for m in range(self.lm_max+1) \
                                                  for l in range(m,self.lm_max+1) \
                                                  for n in range(self.n_max+1)])
      header = np.ndarray((5), dtype='int32')
      header[0] = 314          # endian tag
      header[1] = self.version # file version (hard-coded above)
      header[2] = 1            # mode
      header[3] = self.n_max
      header[4] = self.lm_max
      swapwrite(header, fd, byteswap)
      swapwrite(flatcoeffs.real, fd, byteswap)
      swapwrite(flatcoeffs.imag, fd, byteswap)

    fd.close()

  def read(self, filename, mode='add'):
    """
    Read spectral coefficients from file `filename`.  Assumed to be in spectral input format.

    Optional arguments:
      - mode: either 'add' to or 'replace' existing coefficients (default: 'add')
    """

    fd = open(filename, 'rb')
    bs = check_byteswap(fd)
    version = swapread(fd,dtype='int32',count=1,byteswap=bs)[0]
    fmode   = swapread(fd,dtype='int32',count=1,byteswap=bs)[0]
    if fmode == 0:
      f_n_lmn = swapread(fd,dtype='int32',count=1,byteswap=bs)[0]
      f_ns = swapread(fd,dtype='int32',count=f_n_lmn,byteswap=bs)
      f_ls = swapread(fd,dtype='int32',count=f_n_lmn,byteswap=bs)
      f_ms = swapread(fd,dtype='int32',count=f_n_lmn,byteswap=bs)
    elif fmode == 1:
      f_n_max = swapread(fd,dtype='int32',count=1,byteswap=bs)[0]
      f_lm_max = swapread(fd,dtype='int32',count=1,byteswap=bs)[0]
      f_n_lmn = int((f_n_max+1)*(f_lm_max+1)*(f_lm_max+2)/2)
      f_flatindices = [(n,l,m) for m in range(f_lm_max+1) \
                               for l in range(m,f_lm_max+1) \
                               for n in range(f_n_max+1)]
      f_ns, f_ls, f_ms = [i for i in zip(*f_flatindices)]
    else:
      raise Exception("Unknown file mode in read: {:d}".format(fmode,))
    f_flatcoeffs = swapread(fd,dtype='float64',count=f_n_lmn,byteswap=bs).astype('complex')
    f_flatcoeffs.imag = swapread(fd,dtype='float64',count=f_n_lmn,byteswap=bs)
    self.add_mode(f_flatcoeffs, n=f_ns, l=f_ls, m=f_ms, mode=mode)

def main(fformat=None, n_theta=None, lm_max=None, n_r=None, n_max=None, n_phi=None, \
                       rmin=None, rmax=None, aspect_ratio=None, shell_depth=None, \
                       modes=None, expressions=None, filename=None):
  """
  Main function to process and write spectral input to file.
  """
  if fformat == "dense":
    l_lm_max = lm_max
    if lm_max is None and n_theta is None: l_lm_max = 0
    si = SpectralInput(n_theta=n_theta, lm_max=l_lm_max, n_r=n_r, n_max=n_max)
  elif fformat == "sparse":
    si = SpectralInput()
  else:
    raise Exception("Unknown file format: {:s}".format(fformat,))

  # deal with individual modes
  if modes is not None:
    for index, coeff in modes:
      si.add_mode(coeff, n=index[0], l=index[1], m=index[2])

  # deal with any expressions provided
  if expressions is not None:
    import ast, os

    # from https://stackoverflow.com/questions/39379331/python-exec-a-code-block-and-eval-the-last-line
    def exec_then_eval(theta=None, phi=None, radius=None, \
                       rmin=None, rmax=None, aspect_ratio=None, shell_depth=None, \
                       code=None):
      block = ast.parse(code, mode='exec')

      # assumes last node is an expression
      last = ast.Expression(block.body.pop().value)

      _globals = {'theta':theta, 'phi':phi, 'radius':radius,\
                  'rmin':rmin, 'rmax':rmax, \
                  'aspect_ratio':aspect_ratio, 'shell_depth':shell_depth}
      _locals = {}
      exec(compile(block, '<string>', mode='exec'), _globals, _locals)
      return eval(compile(last, '<string>', mode='eval'), _globals, _locals)

    func_kwargs = {'rmin':rmin, 'rmax':rmax, 'aspect_ratio':aspect_ratio, 'shell_depth':shell_depth}
    for expr in expressions:
      func_kwargs['code'] = expr

      func_args = []
      if expr.find('theta') >= 0: func_args.append("theta")
      if expr.find('phi') >= 0: func_args.append("phi")
      if expr.find('radius') >= 0:
        func_args.append("radius")
        if func_kwargs['rmin'] is None or func_kwargs['rmax'] is None:
          func_kwargs['rmin'], func_kwargs['rmax'] = radial_extents(rmin, rmax, aspect_ratio, shell_depth)
      func_argstr, exec_argstr = "", ""
      if len(func_args) > 0: 
        func_argstr = ", ".join([fa+"=None" for fa in func_args])+", "
        exec_argstr = ", ".join([fa+"="+fa for fa in func_args])+", "
      # from https://stackoverflow.com/questions/1409295/set-function-signature-in-python
      func_str = "def func({:s}**kwargs):".format(func_argstr)+os.linesep+\
                 "  return exec_then_eval({:s}**kwargs)".format(exec_argstr)+os.linesep
      _globals, _locals = {'exec_then_eval': exec_then_eval}, {}
      exec(compile(func_str, '<string>', mode='exec'), _globals, _locals)
      func = _locals["func"]

      si.transform_from_rtp_function(func, n_theta=n_theta, n_r=n_r, \
                                     rmin=rmin, rmax=rmax, aspect_ratio=aspect_ratio, shell_depth=shell_depth,
                                     func_kwargs=func_kwargs)

  si.write(filename)

###################################
# main script below               #
###################################
if __name__ == "__main__":
  import argparse
  import re, os

  def parse_mode(modein):
    modestr = " ".join(modein)
    # regular expression to split up values string into a list
    # the list may be comma(,), semicolon(;), space ( ) or newline (\n) delimited
    # if the list is of bracketed items (e.g. tuples) then any delimiters within 
    # the brackets are preserved
    # brackets may be (), [] or {}
    #NODE                     EXPLANATION
    #--------------------------------------------------------------------------------
    #  (?:                      group, but do not capture (1 or more times
    #                           (matching the most amount possible)):
    #--------------------------------------------------------------------------------
    #    [^,; \n([{]              any character except: ',', ';', ' ',
    #                             '\n' (newline), '(', '[', '{'
    #--------------------------------------------------------------------------------
    #   |                        OR
    #--------------------------------------------------------------------------------
    #    \(                       '('
    #--------------------------------------------------------------------------------
    #    [^)]*                    any character except: ')' (0 or more
    #                             times (matching the most amount
    #                             possible))
    #--------------------------------------------------------------------------------
    #    \)                       ')'
    #--------------------------------------------------------------------------------
    #   |                        OR
    #--------------------------------------------------------------------------------
    #    \[                       '['
    #--------------------------------------------------------------------------------
    #    [^]]*                    any character except: ']' (0 or more
    #                             times (matching the most amount
    #                             possible))
    #--------------------------------------------------------------------------------
    #    \]                       ']'
    #--------------------------------------------------------------------------------
    #   |                        OR
    #--------------------------------------------------------------------------------
    #    \{                       '{'
    #--------------------------------------------------------------------------------
    #    [^}]*                    any character except: '}' (0 or more
    #                             times (matching the most amount
    #                             possible))
    #--------------------------------------------------------------------------------
    #    \}                       '}'
    #--------------------------------------------------------------------------------
    #  )+                       end of grouping
    # - from http://rick.measham.id.au/paste/explain.pl
    r = re.compile(r'(?:[^,; '+os.linesep+'([{]|\([^)]*\)|\[[^]]*\]|\{[^}]*\})+')
    modelist = r.findall(modestr)
    for index in modelist[:min(3,len(modelist)-1)]:
      if not index.isdigit():
        raise argparse.ArgumentTypeError("Expected non-negative integer as index, not {:s}.".format(index,))
    modecoeff = complex(modelist[-1])
    modeindex = [0]*3
    if len(modelist) == 2:
      # interpret this as a Chebyshev mode, with only radial dependence
      # the first entry is then interpretted as a mode index (int) and the second the coefficient value (complex)
      modeindex[0] = int(modelist[0])
    elif len(modelist) == 3:
      # interpret this as a spherical mode, with only (theta, phi) dependence
      # the first two entries are l and m and the third the coefficient value (complex)
      modeindex[1:] = [int(index) for index in modelist[:2]]
    elif len(modelist) == 4:
      # interpret this as a mode with both radial and spherical dependence
      # the first three entries are n, l and m and the fourth the coefficient value (complex)
      modeindex[:] = [int(index) for index in modelist[:3]]
    else:
      raise argparse.ArgumentTypeError("Expected 2 (radial), 3 (spherical) or 4 (radial & spherical) arguments describing a mode, "
                                       "got {:d} ({:s}).".format(len(modelist), modestr))
    return (tuple(modeindex), modecoeff)

  epilog = '''
EXAMPLES
To write a single constant mode, (n,l,m)=(0,0,0), with coefficient 1.+0.j, 
to the file `example` run:

 > %(prog)s -m 0 0 0 1.+0.j -o example

or:

 > %(prog)s -m 0 0 1.+0.j -o example

where n is assumed to be 0 when not supplied, or:

 > %(prog)s -m 0 1.+0.j -o example

where (l,m) is assumed to be (0,0) when not supplied.

To write spectral input matching the Christensen et al. 2001 hydrodynamic 
benchmark initial condition to the file `example`, run:
              
 > %(prog)s -ar 0.35 -sd 1.0 -nt 96 -nr 64 -o example \\
    -e 'import numpy as np; x = 2*radius - rmin - rmax; 
rmax*rmin/radius - rmin + 210*0.1*(1 - 3*x*x + 3*(x**4) - x**6)*(np.sin(theta)**4)*np.cos(4*phi)/np.sqrt(17920*np.pi)'

           '''
  parser = argparse.ArgumentParser( \
                         description="""Generate generic spectral input for Rayleigh.""", epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument("-m", "--mode", nargs='+', type=str, dest='modes', metavar='mode', default=[], action='append', required=False,
                      help='''Add a mode to the spectral input.  This should be formatted as "-m index coefficient"
                              where index is either 1, 2 or 3 integers depending on the desired mode (see below).
                              The coefficient should be parseable as a complex number.  
                              For a pure Chebyshev mode supply 1 integer as the n index followed by the coefficient value, 
                              e.g. "-m 0 1.+1.j".  
                              For a purely spherical harmonic mode supply 2 integers representing the l and m degree and order
                              respectively, followed by the coefficient value, e.g. "-m 0 0 1.+1.j".  
                              For a mode with both radial and spherical dependence supply 3 integers representing n, l and m Chebyshev
                              index and spherical harmonic degree and order respectively, followed by the coefficient value,
                              e.g. "-m 0 0 0 1.+1.j".  
                              All three examples here add the coefficient 1.+1.j to the mode (n,l,m) = (0,0,0).  Multiple modes can
                              be added by using "-m index coefficient" repeatedly.''')
  parser.add_argument("-e", "--expr", type=str, dest='expressions', metavar='expr', default=None, action='append', required=False,
                      help='''Transform the given expression into Chebyshev-spectral space.  The expression can depend on any
                              combination (or none) of the variables `radius`, `theta` (co-latitude), `phi` (longitude).  
                              In addition it may use `rmin`, `rmax`, `aspect_ratio` and
                              `shell_depth`.  The expression should return the field value at the given radius, theta and phi.  
                              It may be vectorized to process multiple radii, thetas and phis at once.  If multiple expressions are
                              supplied their modes will be added.  Similarly any modes supplied (-m) will be added.''')
  parser.add_argument("-rn", "--rmin", type=float, dest='rmin', default=None, action='store', required=False,
                      help='''Supply the minimum radius of the domain.  Required if transforming from an expression that depends on
                              radius and `aspect_ratio` is not supplied along with either `rmax` or `shell_depth`.  Ignored if no
                              expression supplied (-e) or the expression does not depend on radius.''')
  parser.add_argument("-rx", "--rmax", type=float, dest='rmax', default=None, action='store', required=False,
                      help='''Supply the maximum radius of the domain.  Required if transforming from an expression that depends on
                              radius and `aspect_ratio` and `shell_depth` are not supplied.  Ignored if no
                              expression supplied (-e) or the expression does not depend on radius.''')
  parser.add_argument("-sd", "--shell_depth", type=float, dest='shell_depth', default=None, action='store', required=False,
                      help='''Supply the shell depth of the domain.  Required if transforming from an expression that depends on
                              radius and `rmax` is not supplied.  Ignored if no
                              expression supplied (-e) or the expression does not depend on radius.''')
  parser.add_argument("-ar", "--aspect_ratio", type=float, dest='aspect_ratio', default=None, action='store', required=False,
                      help='''Supply the shell depth of the domain.  Required if transforming from an expression that depends on
                              radius and `rmax` and `rmin` are not supplied.  Ignored if no
                              expression supplied (-e) or the expression does not depend on radius.''')
  parser.add_argument("-nt", "--n_theta", type=int, dest='n_theta', default=None, action='store', required=False,
                      help='''Specify the number of co-latitudinal grid points.  Required if `lm_max` is not supplied and either
                              a dense format is requested or an expression that depends on theta is supplied.''')
  parser.add_argument("-np", "--n_phi", type=int, dest='n_phi', default=None, action='store', required=False,
                      help='''Specify the number of longitudinal grid points.  Not required.  Set from `n_theta` if not
                              supplied.''')
  parser.add_argument("-nr", "--n_r", type=int, dest='n_r', default=None, action='store', required=False,
                      help='''Specify the number of radial grid points.  Required if an expression that depends on radius is
                              supplied and n_max is not specified.''')
  parser.add_argument("-lm", "--lm_max", type=int, dest='lm_max', default=None, action='store', required=False,
                      help='''Specify the maximum Legendre order and degree.  Required if `n_theta` is not supplied and either a
                              dense format is requested or an expression that depends on theta is supplied.''')
  parser.add_argument("-nm", "--n_max", type=int, dest='n_r', default=None, action='store', required=False,
                      help='''Specify the maximum Chebyshev polynomial degree.  Required if an expression that depends on radius is
                              supplied.''')
  parser.add_argument("-f", "--format", type=str, choices=["dense", "sparse"], dest='fformat', default="sparse", action='store', required=False,
                      help='''Storage format, either `dense` or `sparse`.  Defaults to `%(default)s`.''')
  parser.add_argument("-o", "--output", type=str, dest='filename', required=True, action='store',
                      help='''Specify the filename of the output file.''')
  args = parser.parse_args()

  dargs = vars(args)
  dargs['modes'] = [parse_mode(mode) for mode in args.modes]

  main(**dargs)

