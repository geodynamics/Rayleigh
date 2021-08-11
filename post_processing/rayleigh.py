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
import glob

import numpy as np
import matplotlib.pyplot as plt

import lut

if sys.maxsize < 2**63 - 1:
    # We don't want mmap on 32-bit systems where virtual memory is limited.
    use_mmap = False
else:
    use_mmap = True

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

    def get_value(self, dtype: str, shape=[1]):
        dtype = np.dtype(dtype).newbyteorder(self.endian)
        size = np.product(shape)
        if self._memmap:
            out = np.ndarray(shape, dtype=dtype, buffer=self.fh,
                             offset=self.fh.tell(), order='F')
            self.fh.seek(dtype.itemsize * size, os.SEEK_CUR)
        else:
            out = np.fromfile(self.fh, dtype=dtype, count=size)
            out.shape = shape
        if size == 1:
            return out[0]
        else:
            return out

class TimeSeries(object):
    pass

class Spherical_3D_grid(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.nr, self.ntheta, self.nphi = self.get_value('i4', shape=[3])
        assert(self.nphi == 2 * self.ntheta)

        self.rs = self.get_value('f8', shape=[self.nr])
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

        self.rs = grid.rs
        self.thetas = grid.thetas
        self.phi_edge = np.linspace(0., 2 * np.pi, grid.nphi)
        self.phis = 0.5 * (self.phi_edge[1:] + self.phi_edge[:-1])

        self.endian = grid.endian

    def q(self, q):
        f = os.path.join(self.directory, "{:08d}_{:04d}".format(self.snap, q))
        return Spherical_3D_value(f, len(self.rs), len(self.thetas),
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

        for f in files:
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

class Shell_Slice_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')

        self.ntheta = self.get_value('i4')
        self.nphi = 2 * self.ntheta
        self.nr = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.rs = self.get_value('f8', shape=[self.nr])
        self.inds = self.get_value('i4', shape=[self.nr]) - 1
        self.costheta = self.get_value('f8', shape=[self.ntheta])
        self.sintheta = np.sqrt(1.0 - self.costheta**2)

        self.val = self.get_value('f8', shape=[self.nphi, self.ntheta, self.nr, self.nq, self.nrec])

class G_Avgs_file(BaseFile):
    def __init__(self, filename, **kwargs):
        super().__init__(filename, **kwargs)

        self.version = self.get_value('i4')
        self.nrec = self.get_value('i4')
        self.nq = self.get_value('i4')

        self.qv = self.get_value('i4', shape=[self.nq])

        self.vals = np.empty((self.nrec, self.nq), dtype='f8')
        self.time = np.empty((self.nrec,), dtype='f8')
        self.iters = np.empty((self.nrec,), dtype='i4')

        buf = self.get_value(np.dtype([('val', 'f8', (self.nq,)),
                                       ('time', 'f8'),
                                       ('iters', 'i4')]),
                              shape=[self.nrec])

        self.vals = buf['val']
        self.time = buf['time']
        self.iters = buf['iters']

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

        if self.version == 1:
            dtype = np.dtype([('val', 'f8', (self.nr, self.nq)),
                              ('time', 'f8'),
                              ('iters', 'i4')])
        elif self.version < 6:
            dtype = np.dtype([('val', 'f8', (self.nr, 4, self.nq)),
                              ('time', 'f8'),
                              ('iters', 'i4')])
        else:
            # fixme
            dtype = np.dtype([('val', 'f8', (self.nr, 4, self.nq)),
                              ('time', 'f8'),
                              ('iters', 'i4')])

        buf = self.get_value(dtype, shape=[self.nrec])
        self.val = buf['val']
        self.time = buf['time']
        self.iters = buf['iters']
