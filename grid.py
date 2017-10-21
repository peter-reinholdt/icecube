#!/usr/bin/env python
from __future__ import division
import numpy as np
import horton


class griddata(object):
    def __init__(self, qmfilename, npts=np.array([120, 120, 120]), bufsize=5.0):
        #load qm data
        IO                  = horton.IOData.from_file(qmfilename)
        self.dm             = IO.get_dm_full()
        self.icharges       = IO.numbers
        self.fcharges       = IO.numbers.astype(np.float64)
        self.obasis         = IO.obasis
        self.natoms         = len(self.icharges)
        self.coordinates    = IO.coordinates

        self.center_of_charge = np.sum(self.coordinates*self.fcharges[:,np.newaxis], axis=0) / np.sum(self.fcharges)
        
        pmin = self.center_of_charge - bufsize
        pmax = self.center_of_charge + bufsize
        print(pmin, pmax)
        self.origin = pmin

        self.Nx = npts[0]
        self.Ny = npts[1]
        self.Nz = npts[2]
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
        self.xyzgrid = np.zeros((self.Nx*self.Ny*self.Nz,3), dtype=np.float64)
        self.make_xyz_grid()
        self.data    = np.zeros((self.Nx*self.Ny*self.Nz),   dtype=np.float64)



    def make_xyz_grid(self):
        counter = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    self.xyzgrid[counter, 0] = self.xrange[i]
                    self.xyzgrid[counter, 1] = self.yrange[j]
                    self.xyzgrid[counter, 2] = self.zrange[k]
                    counter += 1
    

    def compute_density(self):
        self.data = self.obasis.compute_grid_density_dm(self.dm, self.xyzgrid, self.data)


    def compute_potential(self):
        self.data = self.obasis.compute_grid_esp_dm(self.dm, self.xyzgrid, self.fcharges, self.data)
