#!/usr/bin/env python
from __future__ import division
import numpy as np
import horton
import sys
from mpi4py import MPI
from surfaces import compute_vdW_surface
from utilities import container


class griddata(object):
    def __init__(self, qm_file_name):
        IO = horton.IOData.from_file(qm_file_name)
        self.dm = IO.get_dm_full()
        self.icharges = IO.numbers
        self.fcharges = IO.numbers.astype(np.float64)
        self.obasis = IO.obasis
        self.natoms = len(self.icharges)
        self.coordinates = IO.coordinates
        self.qm_file_name = qm_file_name
        self.xyzgrid = np.zeros((0,3), dtype=np.float64)
        self.data = np.zeros((0,),  dtype=np.float64)
        self.grid_type = None


    def get_cubic_grid(self, grid_buffer, grid_density):
        pmin = np.min(self.coordinates, axis=0) - grid_buffer
        pmax = np.max(self.coordinates, axis=0) + grid_buffer
        npts = ((pmax - pmin)*grid_density).astype(np.int)
        self.origin = pmin
        self.Nx = npts[0]
        self.Ny = npts[1]
        self.Nz = npts[2]
        self.xy = 0.0
        self.xz = 0.0
        self.yx = 0.0
        self.yz = 0.0
        self.zx = 0.0
        self.zy = 0.0
        self.xrange  = np.linspace(pmin[0], pmax[0], npts[0])
        self.yrange  = np.linspace(pmin[1], pmax[1], npts[1])
        self.zrange  = np.linspace(pmin[2], pmax[2], npts[2])
        self.xx = self.xrange[1] - self.xrange[0]
        self.yy = self.yrange[1] - self.yrange[0]
        self.zz = self.zrange[1] - self.zrange[0]
        self.xyzgrid = np.zeros((self.Nx*self.Ny*self.Nz,3), dtype=np.float64)
        #self.make_xyz_grid()
        self.data    = np.zeros((self.Nx*self.Ny*self.Nz),   dtype=np.float64)
        counter = 0
        for i in range(self.Nx):
            for j in range(self.Ny):
                for k in range(self.Nz):
                    self.xyzgrid[counter, 0] = self.xrange[i]
                    self.xyzgrid[counter, 1] = self.yrange[j]
                    self.xyzgrid[counter, 2] = self.zrange[k]
                    counter += 1
        self.grid_type = 'cubic'
        self.grid_buffer = grid_buffer
        self.grid_density = grid_density



    def get_vdW_surface(self, surface_vdW_scale, surface_point_density):
        self.xyzgrid = compute_vdW_surface(self.icharges, self.coordinates, surface_point_density=surface_point_density, surface_vdW_scale=surface_vdW_scale)
        self.data = np.zeros((self.xyzgrid.shape[0]), dtype=np.float64) 
        self.grid_type = 'surface'
        self.surface_vdW_scale = surface_vdW_scale
        self.surface_point_density = surface_point_density


    def compute_density(self, nprocs=1):
        if nprocs == 1:
            self.data = self.obasis.compute_grid_density_dm(self.dm, self.xyzgrid)
        else:
            comm = MPI.COMM_SELF.Spawn(sys.executable, 
                    args="{}/density_worker.py".format(sys.path[0]), 
                    maxprocs=nprocs)
            for iproc in range(nprocs):
                if self.grid_type == 'cubic':
                    work_info = container(qm_file_name=self.qm_file_name, grid_type=self.grid_type, grid_buffer=self.grid_buffer, grid_density=self.grid_density)
                elif self.grid_type == 'surface':
                    work_info = container(qm_file_name=self.qm_file_name, grid_type=self.grid_type, surface_vdW_scale=self.surface_vdW_scale, surface_point_density=self.surface_point_density)
                comm.send(work_info, dest=iproc, tag=11)
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
                if self.grid_type == 'cubic':
                    work_info = container(qm_file_name=self.qm_file_name, grid_type=self.grid_type, grid_buffer=self.grid_buffer, grid_density=self.grid_density)
                elif self.grid_type == 'surface':
                    work_info = container(qm_file_name=self.qm_file_name, grid_type=self.grid_type, surface_vdW_scale=self.surface_vdW_scale, surface_point_density=self.surface_point_density)
                comm.send(work_info, dest=iproc, tag=11)
            chunksize = self.xyzgrid.shape[0] // nprocs
            for iproc in range(nprocs):
                start = iproc * chunksize
                end   = (iproc+1) * chunksize
                if iproc == (nprocs - 1):
                    end = self.xyzgrid.shape[0] 
                comm.Recv(self.data[start:end], source=iproc, tag=13)
            comm.Disconnect()
