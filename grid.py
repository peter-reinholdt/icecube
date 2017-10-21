#!/usr/bin/env python
from __future__ import division
import numpy as np
import horton


class griddata(object):
    def __init__(self, qmfilename, npts=np.array([120, 120, 120]), bufsize=10.0):
        #load qm data
        IO                  = horton.IOData.from_file(qmfilename)
        self.dm             = IO.dm
        self.icharges       = IO.numbers
        self.fcharges       = IO.numbers.astype(np.float64)
        self.obasis         = IO.obasis
        self.natoms         = len(self.icharges)
        self.coordinates    = IO.coordinates

        #get origin pmin, pmax
        pmin = np.min(self.coordinates, axis=0) - bufsize
        pmax = np.max(self.coordinates, axis=0) + bufsize

        origin = pmin

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
        self.xrange  = np.linspace(pmin[0], pmax[0], npts[0])
        self.yrange  = np.linspace(pmin[1], pmax[1], npts[1])
        self.zrange  = np.linspace(pmin[2], pmax[2], npts[2])
        self.xyzgrid = np.zeros((Nx*Ny*Nz,3), dtype=np.float64)
        self.make_xyz_grid()
        self.data    = np.zeros((Nx*Ny*Nz),   dtype=np.float64)



    def make_xyz_grid(self):
        counter = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    self.xyzgrid[counter, 0] = self.xrange[i]
                    self.xyzgrid[counter, 1] = self.xrange[j]
                    self.xyzgrid[counter, 2] = self.xrange[k]
                    counter += 1
    

    def compute_density(self):
        self.data = self.compute_grid_esp_dm(self.dm, self.xyzgrid, self.data)


    def compute_potential(self):
        self.data = self.compute_grid_esp_dm(self.dm, self.xyzgrid, self.fcharges, self.data)
