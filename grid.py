#!/usr/bin/env python
import numpy as np


class griddata(object):
    def __init__(self, origin, pmin, pmax, npts):
        self.Nx = npts[0]
        self.Ny = npts[1]
        self.Nz = npts[2]
        self.origin = np.array([pmin[0], pmin[1], pmin[2]], dtype=np.float64)
        self.xx = (pmax[0] - pmin[0]) / npts[0]
        self.yy = (pmax[1] - pmin[1]) / npts[1]
        self.zz = (pmax[2] - pmin[2]) / npts[2]
        self.xy = 0.0
        self.xz = 0.0
        self.yx = 0.0
        self.yz = 0.0
        self.zx = 0.0
        self.zy = 0.0
        self.xrange = np.linspace(pmin[0], pmax[0], npts[0])
        self.yrange = np.linspace(pmin[1], pmax[1], npts[1])
        self.zrange = np.linspace(pmin[2], pmax[2], npts[2])
        self.xyzgrid = np.zeros((Nx*Ny*Nz,3))


    def make_xyz_grid(self):
        counter = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    self.xyzgrid[counter, 0] = self.xrange[i]
                    self.xyzgrid[counter, 1] = self.xrange[j]
                    self.xyzgrid[counter, 2] = self.xrange[k]
                    counter += 1
