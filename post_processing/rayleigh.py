#
#  Copyright (C) 2021 by the authors of the RAYLEIGH code.
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

import sys
import os
import mmap
import collections.abc
import abc

import numpy as np
import matplotlib.pyplot as plt

import lut

if sys.maxsize < 2**63 - 1:
    # We don't want mmap on 32-bit systems where virtual memory is limited.
    use_mmap = False
else:
    use_mmap = True

try:
    import tqdm
    try:
        in_notebook = get_ipython().__class__.__name__ == 'ZMQInteractiveShell'
    except NameError:
        in_notebook = False

    if in_notebook:
        progress = tqdm.tqdm_notebook
    else:
        progress = tqdm.tqdm
except ImportError:
    progress = lambda x: x

def get_bounds(a, start, end):
    a = 0.5 * (a[:-1] + a[1:])
    return np.concatenate([[start], a, [end]], axis=0)

def format_time(t: float, unit=None, digits=2, return_scaled_t=False):
    factors = {'s': 1., 'min': 60., 'hour': 3600., 'day': 24.*3600., 'year': 31556926.}
    if unit is None:
        fac = 1.
        name = ''
    elif isinstance(unit, str):
        fac = factors[unit]
        name = ' ' + unit
    else:
        raise ValueError('unknown value for unit "{}"'.format(unit))

    if return_scaled_t:
        return t/fac
    else:
        return ('{0:.%df}{1}'%digits).format(t/fac, name)

class BaseFile(object):
    @staticmethod
    def get_endian(fd, sig: int, sigtype) -> str:
        """returns > if the file is big endian and < if the file is little endian"""
        dtype = np.dtype(sigtype)
        buf = fd.read(dtype.itemsize)
        if np.frombuffer(buf, dtype="<"+sigtype, count=1)[0] == sig:
            return "<"
        elif np.frombuffer(buf, dtype=">"+sigtype, count=1)[0] == sig:
            return ">"
        else:
            raise IOError("could not determine endianness")

    def __init__(self, filename: str, endian=None, memmap=use_mmap):

        self.fh = open(filename, "rb")
        self._memmap = memmap

        if self._memmap:
            buf = mmap.mmap(self.fh.fileno(), 0, mmap.MAP_SHARED, mmap.ACCESS_READ)
            self.fh.close()
            self.fh = buf

        if endian is None:
            self.endian = self.get_endian(self.fh, 314, 'i4')
        else:
            self.endian = endian

    def get_value(self, dtype: str, shape=None):
        dtype = np.dtype(dtype).newbyteorder(self.endian)
        if shape is None:
            scalar = True
            shape = [1]
        else:
            scalar = False
        size = np.product(shape)
        if self._memmap:
            out = np.ndarray(shape, dtype=dtype, buffer=self.fh,
                             offset=self.fh.tell(), order='F')
            self.fh.seek(dtype.itemsize * size, os.SEEK_CUR)
        else:
            out = np.fromfile(self.fh, dtype=dtype, count=size)
            out.shape = shape
        if scalar:
            return out[0]
        else:
            return out

    @property
    def qvmap(self):
        return {v: i for i, v in enumerate(self.qv)}


class TimeSeries(object):
    pass

class Spherical_3D_grid(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.nr, self.ntheta, self.nphi = self.get_value('i4', shape=[3])
        assert(self.nphi == 2 * self.ntheta)

        self.radius = self.get_value('f8', shape=[self.nr])
        self.thetas = self.get_value('f8', shape=[self.ntheta])

class Spherical_3D_value(np.ndarray):
    def __new__(cls, filename, nr, ntheta, nphi, endian, **kwargs):
        f = BaseFile(filename, endian=endian, **kwargs)

        return f.get_value('f8', shape=[nphi, ntheta, nr]).view(type=cls)

class Spherical_3D_TimeSeries(TimeSeries):
    def __init__(self, directory, qcode, snaps):
        super().__init__()

        self.directory = directory
        self.qcode = qcode
        self.snaps = snaps

    def __getitem__(self, ind):
        def getone(i):
            f = os.path.join(self.directory, "{:08d}_grid".format(i))
            grid = Spherical_3D_grid(f)
            f = os.path.join(self.directory, "{:08d}_{:04d}".format(i, self.qcode))
            return Spherical_3D_value(f, grid.nr, grid.ntheta, grid.nphi,
                                      endian=grid.endian)
        if np.isscalar(ind):
            return getone(self.snaps[ind])
        else:
            return [getone(i) for i in self.snaps[ind]]

class Spherical_3D_Snapshot(object):
    def __init__(self, directory, snap):
        self.directory = directory
        self.snap = snap

        f = os.path.join(self.directory, "{:08d}_grid".format(snap))
        grid = Spherical_3D_grid(f)

        self.radius = grid.radius
        self.thetas = grid.thetas
        self.phi_edge = np.linspace(0., 2 * np.pi, grid.nphi)
        self.phis = 0.5 * (self.phi_edge[1:] + self.phi_edge[:-1])

        self.endian = grid.endian

    def q(self, q):
        f = os.path.join(self.directory, "{:08d}_{:04d}".format(self.snap, q))
        return Spherical_3D_value(f, len(self.radius), len(self.thetas),
                                  len(self.phis), endian=self.endian)

    def __getattr__(self, q):
        qcode = lut.parse_quantity(q)[0]
        if qcode is None:
            raise AttributeError("unknown quantity ({})".format(q))
        return self.q(qcode)


class Spherical_3D(object):
    def __init__(self, directory='Spherical_3D'):
        super().__init__()

        self.directory = directory
        files = os.listdir(directory)
        self.snaps = set()
        self.quants = set()

        for f in progress(files):
            snap, quant = f.split('_')
            if quant == 'grid':
                continue
            s = int(snap)
            q = int(quant)
            assert("{:08d}".format(s) == snap)
            assert("{:04d}".format(q) == quant)
            self.snaps.add(s)
            self.quants.add(q)

        self.snaps = list(self.snaps)
        self.snaps.sort()
        self.quants = list(self.quants)
        self.quants.sort()

    def q(self, q):
        return Spherical_3D_TimeSeries(self.directory, q, self.snaps)

    def __getitem__(self, ind):
        return Spherical_3D_Snapshot(self.directory, ind)

    def __getattr__(self, q):
        qcode = lut.parse_quantity(q)[0]
        if qcode is None:
            raise AttributeError("unknown quantity ({})".format(q))
        return self.q(qcode)

class Rayleigh_TimeSeries(collections.abc.Sequence):
    def __init__(self, base, qcode):
        self.base = base
        self.qcode = qcode

    def __getitem__(self, i):
        return self.base.get_q(i, self.qcode)

    def __len__(self):
        return len(self.base)


class Rayleigh_TimeStep(object):
    def __init__(self, base, i):
        self.base = base
        self.i = i

    def __getattr__(self, name):
        qcode = lut.parse_quantity(name)[0]
        if qcode is None:
            raise AttributeError("unknown quantity '{}'".format(name))
        return self.base.get_q(self.i, qcode)


class Rayleigh_Output(collections.abc.Sequence):
    @abc.abstractmethod
    def get_q(self, i, qcode):
        pass

    def __len__(self):
        return len(self.val)

    def __getattr__(self, name):
        qcode = lut.parse_quantity(name)[0]
        if qcode is None:
            raise AttributeError("unknown quantity '{}'".format(name))
        return Rayleigh_TimeSeries(self, qcode)

    def __getitem__(self, i):
        return Rayleigh_TimeStep(self, i)

    def __init__(self, filecls, directory, subrange=None):
        super().__init__()

        self.directory = directory
        files = os.listdir(directory)
        files.sort()

        if subrange is not None:
            if isinstance(subrange, int):
                subrange = range(0, len(files), subrange)
            files = files[subrange]

        self.val = []
        self.time = []
        self.iter = []
        self.gridpointer = []

        for a in self.attrs:
            setattr(self, a, [])

        for i, f in enumerate(progress(files)):
            m = filecls(os.path.join(directory, f))
            self.val += m.val
            self.time += m.time
            self.iter += m.iter

            self.gridpointer += len(m.val) * [i]

            for a in self.attrs:
                getattr(self, a).append(getattr(m, a))

class Plot2D(abc.ABC):
    @abc.abstractmethod
    def get_coords(self, i):
        pass

    @abc.abstractmethod
    def get_coord_labels(self):
        pass

    def pcolor(self, i, q, Clear=True, iv=0, Colorbar=True, projection=None,
               tunit=None, **kwargs):
        qcode = lut.parse_quantity(q)[0]
        if qcode is None:
            raise AttributeError("unknown quantity ({})".format(q))
        fig = plt.gcf()
        if Clear:
            fig.clear()

        ax = fig.add_subplot(111, projection=projection)

        X, Y = self.get_coords(i)
        Xl, Yl = self.get_coord_labels()
        im = ax.pcolormesh(X, Y, self.get_q(i, qcode)[iv, :, :], **kwargs)
        if Colorbar:
            cbar = plt.colorbar(im, ax=ax)
            cbar.set_label(f"{lut.latex_formula(q)}")

        ax.set_xlabel(Xl)
        ax.set_ylabel(Yl)

        ax.set_title(f"$t=${format_time(self.time[i], unit=tunit)}")
        if projection is None:
            ax.set_aspect('equal')


class Meridional_Slices_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')

        self.nr = self.get_value('i4')
        self.ntheta = self.get_value('i4')
        self.nphi = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.radius = self.get_value('f8', shape=[self.nr])
        self.costheta = self.get_value('f8', shape=[self.ntheta])
        self.sintheta = np.sqrt(1.0 - self.costheta**2)
        self.phi_inds = self.get_value('i4', shape=[self.nphi]) - 1
        if self.nphi == 1:
            self.phi_inds = np.array([self.phi_inds])
        self.phi = np.zeros(self.nphi,dtype='float64')

        dphi = (2*np.pi)/(self.ntheta*2)
        for i in range(self.nphi):
            self.phi[i] = self.phi_inds[i]*dphi

        self.val = []
        self.time = []
        self.iter = []
        for i in range(self.nrec):
            self.val.append(self.get_value('f8', shape=[self.nphi, self.ntheta, self.nr, self.nq]))
            self.time.append(self.get_value('f8'))
            self.iter.append(self.get_value('i4'))


class Meridional_Slices(Rayleigh_Output, Plot2D):
    attrs = ("radius", "costheta", "sintheta", "qvmap")

    def __init__(self, directory='Meridional_Slices'):
        super().__init__(Meridional_Slices_file, directory)

        self.theta = [np.arccos(x) for x in self.costheta]
        self.costheta_bounds = [np.cos(get_bounds(t, np.pi, 0.)) for t in self.theta]
        self.sintheta_bounds = [np.sqrt(1.0 - ct**2) for ct in self.costheta_bounds]
        self.radius_bounds = [get_bounds(r, r[0] + 0.5 * (r[0] - r[1]),
                                         r[-1] - 0.5 * (r[-2] - r[-1]))
                              for r in self.radius]

    def get_coords(self, i):
        igrid = self.gridpointer[i]
        X = self.sintheta_bounds[igrid][:, None] * self.radius_bounds[igrid][None, :]
        Y = self.costheta_bounds[igrid][:, None] * self.radius_bounds[igrid][None, :]
        return X, Y

    def get_coord_labels(self):
        return "$x$", "$y$"

    def get_q(self, i, qcode):
        igrid = self.gridpointer[i]
        return self.val[i][:, :, :, self.qvmap[igrid][qcode]]


class Equatorial_Slices_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')

        self.nphi = self.get_value('i4')
        self.nr = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.radius = self.get_value('f8', shape=[self.nr])

        dphi = 2 * np.pi / self.nphi
        self.phi = np.arange(self.nphi) * dphi

        self.val = []
        self.time = []
        self.iter = []
        for i in range(self.nrec):
            self.val.append(self.get_value('f8', shape=[self.nphi, self.nr, self.nq]))
            self.time.append(self.get_value('f8'))
            self.iter.append(self.get_value('i4'))


class Equatorial_Slices(Rayleigh_Output, Plot2D):
    attrs = ("radius", "phi", "qvmap")

    def __init__(self, directory='Equatorial_Slices'):
        super().__init__(Equatorial_Slices_file, directory)

        self.phi_bounds = [get_bounds(p, 0., 2. * np.pi) for p in self.phi]
        self.radius_bounds = [get_bounds(r, r[0] + 0.5 * (r[0] - r[1]),
                                         r[-1] - 0.5 * (r[-2] - r[-1]))
                              for r in self.radius]

    def get_coords(self, i):
        igrid = self.gridpointer[i]
        X = np.cos(self.phi_bounds[igrid][:, None]) * self.radius_bounds[igrid][None, :]
        Y = np.sin(self.phi_bounds[igrid][:, None]) * self.radius_bounds[igrid][None, :]
        return X, Y

    def get_coord_labels(self):
        return "$x$", "$y$"

    def get_q(self, i, qcode):
        igrid = self.gridpointer[i]
        return self.val[i][None, :, :, self.qvmap[igrid][qcode]]


class Point_Probes_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')

        self.nr = self.get_value('i4')
        self.ntheta = self.get_value('i4')
        self.nphi = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.radius = self.get_value('f8', shape=[self.nr])
        self.rad_inds = self.get_value('i4', shape=[self.nr]) - 1
        self.costheta = self.get_value('f8', shape=[self.ntheta])
        self.sintheta = np.sqrt(1.0 - self.costheta**2)
        self.theta_inds = self.get_value('i4', shape=[self.ntheta]) - 1
        self.phi = self.get_value('f8', shape=[self.nphi])
        self.phi_inds = self.get_value('i4', shape=[self.nphi]) - 1

        self.val = []
        self.time = []
        self.iter = []
        for i in range(self.nrec):
            self.val.append(self.get_value('f8', shape=[self.nphi, self.ntheta, self.nr, self.nq]))
            self.time.append(self.get_value('f8'))
            self.iter.append(self.get_value('i4'))

class Point_Probes(Rayleigh_Output):
    attrs = ("radius", "costheta", "sintheta", "phi", "qvmap")

    def __init__(self, directory='Point_Probes'):
        super().__init__(Point_Probes_file, directory)

        self.theta = [np.arccos(x) for x in self.costheta]

    def get_q(self, i, qcode):
        igrid = self.gridpointer[i]
        return self.val[i][:, :, :, self.qvmap[igrid][qcode]]


class AZ_Avgs_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')

        self.nr = self.get_value('i4')
        self.ntheta = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.radius = self.get_value('f8', shape=[self.nr])
        self.costheta = self.get_value('f8', shape=[self.ntheta])
        self.sintheta = np.sqrt(1.0 - self.costheta**2)

        self.val = []
        self.time = []
        self.iter = []
        for i in range(self.nrec):
            self.val.append(self.get_value('f8', shape=[self.ntheta, self.nr, self.nq]))
            self.time.append(self.get_value('f8'))
            self.iter.append(self.get_value('i4'))


class AZ_Avgs(Rayleigh_Output, Plot2D):
    attrs = ("radius", "costheta", "sintheta", "qvmap")

    def __init__(self, directory='AZ_Avgs'):
        super().__init__(AZ_Avgs_file, directory)

        self.theta = [np.arccos(x) for x in self.costheta]
        self.costheta_bounds = [np.cos(get_bounds(t, np.pi, 0.)) for t in self.theta]
        self.sintheta_bounds = [np.sqrt(1.0 - ct**2) for ct in self.costheta_bounds]
        self.radius_bounds = [get_bounds(r, r[0] + 0.5 * (r[0] - r[1]),
                                         r[-1] - 0.5 * (r[-2] - r[-1]))
                              for r in self.radius]

    def get_coords(self, i):
        igrid = self.gridpointer[i]
        X = self.sintheta_bounds[igrid][:, None] * self.radius_bounds[igrid][None, :]
        Y = self.costheta_bounds[igrid][:, None] * self.radius_bounds[igrid][None, :]
        return X, Y

    def get_coord_labels(self):
        return "$x$", "$y$"

    def get_q(self, i, qcode):
        igrid = self.gridpointer[i]
        return self.val[i][None, :, :, self.qvmap[igrid][qcode]]


class Shell_Slices_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')

        self.ntheta = self.get_value('i4')
        self.nphi = 2 * self.ntheta
        self.nr = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.radius = self.get_value('f8', shape=[self.nr])
        self.inds = self.get_value('i4', shape=[self.nr]) - 1
        self.costheta = self.get_value('f8', shape=[self.ntheta])
        self.sintheta = np.sqrt(1.0 - self.costheta**2)

        dphi = 2 * np.pi / self.nphi
        self.phi = np.arange(self.nphi) * dphi

        self.val = []
        self.time = []
        self.iter = []
        for i in range(self.nrec):
            self.val.append(self.get_value('f8', shape=[self.nphi, self.ntheta, self.nr, self.nq]))
            self.time.append(self.get_value('f8'))
            self.iter.append(self.get_value('i4'))


class Shell_Slices(Rayleigh_Output, Plot2D):
    attrs = ("radius", "costheta", "sintheta", "phi", "qvmap")

    def __init__(self, directory='Shell_Slices'):
        super().__init__(Shell_Slices_file, directory)

        self.theta = [np.arccos(x) for x in self.costheta]
        self.theta_bounds = [get_bounds(t, np.pi, 0.) for t in self.theta]
        self.costheta_bounds = [np.cos(t) for t in self.theta_bounds]
        self.sintheta_bounds = [np.sqrt(1.0 - ct**2) for ct in self.costheta_bounds]
        self.phi_bounds = [get_bounds(p, 0., 2. * np.pi) for p in self.phi]

    def get_coords(self, i):
        igrid = self.gridpointer[i]
        X = np.repeat(self.phi_bounds[igrid][:, None], len(self.theta_bounds[igrid]), axis=1)
        Y = np.repeat(self.theta_bounds[igrid][None, :], len(self.phi_bounds[igrid]), axis=0)
        # Adjust to bounds used for projection.
        X = X - np.pi
        Y = 0.5 * np.pi - Y
        return X, Y

    def get_coord_labels(self):
        return r"$\phi$", r"$\theta$"

    def get_q(self, i, qcode):
        igrid = self.gridpointer[i]
        x = self.val[i][:, :, :, self.qvmap[igrid][qcode]]
        # Move radial coordinate to first position so it can be selected in pcolor.
        return np.moveaxis(x, 2, 0)


class SPH_Modes_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')

        self.nell = self.get_value('i4')
        self.nr = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.radius = self.get_value('f8', shape=[self.nr])
        self.inds = self.get_value('i4', shape=[self.nr]) - 1
        self.lvals = self.get_value('i4', shape=[self.nell])

        nm = np.max(self.lvals) + 1

        self.val = []
        self.time = []
        self.iter = []
        for i in range(self.nrec):
            self.val.append(np.zeros([nm, self.nell, self.nr, self.nq], dtype='c16'))
            for qv in range(self.nq):
                for p in range(2):
                    for rr in range(self.nr):
                        for lv in range(self.nell):
                            nm = self.lvals[lv] + 1
                            tmp = self.get_value('f8', shape=[nm])
                            if p == 0:
                                self.val[-1][:nm, lv, rr, qv].real = tmp
                            else:
                                self.val[-1][:nm, lv, rr, qv].imag = tmp
            if (self.version < 4):
                # The m>0 --power-- is too high by a factor of 2
                # We divide the --complex amplitude-- by sqrt(2)
                self.vals[-1][1:,:,:,:] /= np.sqrt(2.0)
            self.time.append(self.get_value('f8'))
            self.iter.append(self.get_value('i4'))


class SPH_Modes(Rayleigh_Output):
    attrs = ("radius", "lvals", "qvmap")

    def __init__(self, directory='SPH_Modes'):
        super().__init__(SPH_Modes_file, directory)

    def get_q(self, i, qcode):
        igrid = self.gridpointer[i]
        return self.val[i][:, :, :, self.qvmap[igrid][qcode]]


class Shell_Avgs_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')

        self.nr = self.get_value('i4')
        self.nq = self.get_value('i4')

        if self.version >= 6:
            npcol = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.radius = self.get_value('f8', shape=[self.nr])

        self.val = []
        self.time = []
        self.iter = []

        for i in range(self.nrec):
            if self.version == 1:
                self.val.append(self.get_value('f8', shape=[self.nr, self.nq]))
            elif self.version < 6:
                self.val.append(self.get_value('f8', shape=[self.nr, 4, self.nq]))
            else:
                self.val.append(np.zeros([self.nr, 4, self.nq], dtype='f8'))
                rind = 0
                nr_base = self.nr // npcol
                nr_mod = self.nr % npcol
                for j in range(npcol):
                    nrout= nr_base
                    if (j < nr_mod) :
                        nrout=nrout+1
                    self.val[-1][rind:rind+nrout,:,:] = self.get_value('f8', shape=[nrout, 4, self.nq])
                    rind=rind+nrout
            self.time.append(self.get_value('f8'))
            self.iter.append(self.get_value('i4'))


class Shell_Avgs(Rayleigh_Output):
    attrs = ("radius", "qvmap")

    def __init__(self, directory='Shell_Avgs'):
        super().__init__(Shell_Avgs_file, directory)

    def get_q(self, i, qcode):
        igrid = self.gridpointer[i]
        return self.val[i][:, :, self.qvmap[igrid][qcode]]


class G_Avgs_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.val = []
        self.time = []
        self.iter = []
        for i in range(self.nrec):
            self.val.append(self.get_value('f8', shape=[self.nq]))
            self.time.append(self.get_value('f8'))
            self.iter.append(self.get_value('i4'))


class G_Avgs(Rayleigh_Output):
    attrs = ("qvmap",)

    def __init__(self, directory='G_Avgs'):
        super().__init__(G_Avgs_file, directory)

    def get_q(self, i, qcode):
        igrid = self.gridpointer[i]
        return self.val[i][self.qvmap[igrid][qcode]]


class Shell_Spectra_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')

        self.lmax = self.get_value('i4')
        self.nell = self.lmax + 1
        self.nm = self.nell
        self.mmax = self.nm - 1
        self.nr = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.radius = self.get_value('f8', shape=[self.nr])
        self.inds = self.get_value('i4', shape=[self.nr]) - 1

        self.val = []
        self.time = []
        self.iter = []
        for i in range(self.nrec):
            self.val.append(np.zeros([self.nell, self.nm, self.nr, self.nq], dtype='c16'))
            self.val[-1].real = self.get_value('f8', shape=[self.nell, self.nm, self.nr, self.nq])
            self.val[-1].imag = self.get_value('f8', shape=[self.nell, self.nm, self.nr, self.nq])
            if (self.version != 4):
                # The m>0 --power-- is too high by a factor of 2
                # We divide the --complex amplitude-- by sqrt(2)
                self.val[-1][:,1:,:,:] /= np.sqrt(2.0)
            self.time.append(self.get_value('f8'))
            self.iter.append(self.get_value('i4'))


class Shell_Spectra(Rayleigh_Output):
    attrs = ("radius", "lmax", "mmax", "qvmap")

    def __init__(self, directory='Shell_Spectra'):
        super().__init__(Shell_Spectra_file, directory)

    def lpower(self, i, qcode):
        val = self.get_q(i, qcode)

        # m = 0 power
        x = np.abs(val[:, 0, :, None])**2

        # m != 0 (convective) power
        y = (np.abs(val[:, 1:, :, None])**2).sum(axis=1)

        return np.concatenate((x + y, x, y), axis=2)

    def get_q(self, i, qcode):
        igrid = self.gridpointer[i]
        return self.val[i][:, :, :, self.qvmap[igrid][qcode]]


class PDE_Coefficients(BaseFile):
    nconst = 10
    nfunc = 14

    def __init__(self, filename='equation_coefficients', **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.cset = self.get_value('i4', shape=[self.nconst])
        self.fset = self.get_value('i4', shape=[self.nfunc])

        self.constants = self.get_value('f8', shape=[self.nconst])
        self.nr = self.get_value('i4')
        self.radius = self.get_value('f8', shape=[self.nr])
        self.functions = self.get_value('f8', shape=[self.nr, self.nfunc])

        # aliases
        self.density = self.rho = self.functions[:,1-1]
        self.dlnrho  = self.functions[:,8-1]
        self.d2lnrho = self.functions[:,9-1]

        self.temperature = self.T = self.functions[:,4-1]
        self.dlnT        = self.functions[:,10-1]

        self.dsdr = self.functions[:,14-1]

        self.heating = self.functions[:,6-1]*self.constants[10-1]/self.rho/self.T

        self.nu   = self.functions[:,3-1]
        self.dlnu = self.functions[:,11-1]
        self.kappa    = self.functions[:,5-1]
        self.dlnkappa = self.functions[:,12-1]
        self.eta    = self.functions[:,7-1]
        self.dlneta = self.functions[:,13-1]
