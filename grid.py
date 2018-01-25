#!/usr/bin/env python
from __future__ import division
import numpy as np
import horton
from mpi4py import MPI


class griddata(object):
    def __init__(self, qmfilename, npts=np.array([120, 120, 120]), bufsize=2.5):
        #load qm data
        IO                  = horton.IOData.from_file(qmfilename)
        self.dm             = IO.get_dm_full()
        self.icharges       = IO.numbers
        self.fcharges       = IO.numbers.astype(np.float64)
        self.obasis         = IO.obasis
        self.natoms         = len(self.icharges)
        self.coordinates    = IO.coordinates

        #self.center_of_charge = np.sum(self.coordinates*self.fcharges[:,np.newaxis], axis=0) / np.sum(self.fcharges)
        pmin = np.min(self.coordinates, axis=0) - bufsize
        pmax = np.max(self.coordinates, axis=0) + bufsize
        #pmin = self.center_of_charge - bufsize
        #pmax = self.center_of_charge + bufsize
        #print(pmin, pmax)
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
        if MPI.COMM_WORLD.Get_size() == 1:
            self.data = self.obasis.compute_grid_density_dm(self.dm, self.xyzgrid)
        else:
            self.data = _compute_density(self.obasis, self.dm, self.xyzgrid)


    def compute_potential(self):
        if MPI.COMM_WORLD.Get_size() == 1:
            self.data = self.obasis.compute_grid_esp_dm(self.dm, self.coordinates, self.fcharges, self.xyzgrid)
        else:
            self.data = _compute_potential(self.obasis, self.dm, self.coordinates, self.fcharges, self.xyzgrid)

def _compute_density(obasis, dm, xyzgrid):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    #each rank is responsible for a subset of the grid
    chunksize = xyzgrid.shape[0] // size
    start = rank * chunksize
    end   = (rank+1) * chunksize
    if rank == (size - 1):
        end = xyzgrid.shape[0] 
    density_chunk = obasis.compute_grid_density_dm(dm, xyzgrid[start:end,:])
    
    #receive data from slaves
    if rank == 0:
        data = np.zeros(xyzgrid.shape[0])
        loc = 0
        data[loc:loc+density_chunk.shape[0]] = density_chunk
        loc += density_chunk.shape[0]
        for i in range(1,size):
            print("I am rank 0 and want to received something from rank {}".format(i))
            comm.Recv(density_chunk, source=i, tag=13)
            data[loc:loc+density_chunk.shape[0]] = density_chunk
            loc += density_chunk.shape[0]
    else:
        comm.Send(density_chunk, dest=0, tag=13)
    if rank == 0:
        return data
    else:
        return


def _compute_potential(obasis, dm, coordinates, charges, xyzgrid):
    comm = MPI.COMM_WORLD
    size = comm.Get_size()
    rank = comm.Get_rank()
    #each rank is responsible for a subset of the grid
    chunksize = xyzgrid.shape[0] // size
    start = rank * chunksize
    end   = (rank+1) * chunksize
    if rank == (size - 1):
        end = xyzgrid.shape[0] 

    potential_chunk = obasis.compute_grid_esp_dm(dm, coordinates, charges, xyzgrid[start:end,:])
    
    #receive data from slaves
    if rank == 0:
        data = np.zeros(xyzgrid.shape[0])
        loc = 0
        data[loc:loc+potential_chunk.shape[0]] = potential_chunk
        loc += potential_chunk.shape[0]
        for i in range(1,size):
            comm.Recv(potential_chunk, source=i, tag=13)
            data[loc:loc+potential_chunk.shape[0]] = potential_chunk
            loc += potential_chunk.shape[0]
    else:
        comm.Send(potential_chunk, dest=0, tag=13)
    if rank == 0:
        return data
    else:
        return
