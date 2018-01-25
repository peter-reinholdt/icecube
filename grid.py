#!/usr/bin/env python
from __future__ import division
import numpy as np
import horton
import sys
from mpi4py import MPI


class griddata(object):
    def __init__(self, qmfilename, grid_density=6.0, bufsize=5.0):
        #load qm data
        IO                  = horton.IOData.from_file(qmfilename)
        self.dm             = IO.get_dm_full()
        self.icharges       = IO.numbers
        self.fcharges       = IO.numbers.astype(np.float64)
        self.obasis         = IO.obasis
        self.natoms         = len(self.icharges)
        self.coordinates    = IO.coordinates
        self.name           = qmfilename

        pmin = np.min(self.coordinates, axis=0) - bufsize
        pmax = np.max(self.coordinates, axis=0) + bufsize
        npts = ((pmax - pmin)*grid_density).astype(np.int)
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
    

    def compute_density(self, nprocs=1):
        print(nprocs, type(nprocs))
        if nprocs == 1:
            self.data = self.obasis.compute_grid_density_dm(self.dm, self.xyzgrid)
        else:
            comm = MPI.COMM_SELF.Spawn(sys.executable, 
                    args="{}/density_worker.py".format(sys.path[0]), 
                    maxprocs=nprocs)
            for iproc in range(nprocs):
                comm.send(self.name, dest=iproc, tag=11)
            chunksize = self.xyzgrid.shape[0] // nprocs
            for iproc in range(nprocs):
                start = iproc * chunksize
                end   = (iproc+1) * chunksize
                if iproc == (nprocs - 1):
                    end = self.xyzgrid.shape[0] 
                comm.Recv(self.data[start:end], source=iproc, tag=13)
            comm.Disconnect()


    def compute_potential(self, nprocs=1):
        if nprocs == 1:
            self.data = self.obasis.compute_grid_esp_dm(self.dm, self.coordinates, self.fcharges, self.xyzgrid)
        else:
            comm = MPI.COMM_SELF.Spawn(sys.executable, 
                    args="{}/potential_worker.py".format(sys.path[0]), 
                    maxprocs=nprocs)
            for iproc in range(nprocs):
                comm.send(self.name, dest=iproc, tag=11)
            chunksize = self.xyzgrid.shape[0] // nprocs
            for iproc in range(nprocs):
                start = iproc * chunksize
                end   = (iproc+1) * chunksize
                if iproc == (nprocs - 1):
                    end = self.xyzgrid.shape[0] 
                comm.Recv(self.data[start:end], source=iproc, tag=13)
            comm.Disconnect()
