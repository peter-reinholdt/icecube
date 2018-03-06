#!/usr/bin/env python

from mpi4py import MPI
import numpy as np
import horton
import sys
import grid
from utilities import container


comm = MPI.Comm.Get_parent()
size = comm.Get_size()
rank = comm.Get_rank()


work_info = comm.recv(source=0, tag=11)
gr = grid.griddata(work_info.qm_file_name)
if work_info.grid_type == 'cubic':
    gr.get_cubic_grid(work_info.grid_buffer, work_info.grid_density)
elif work_info.grid_type == 'surface':
    gr.get_vdW_surface(work_info.surface_vdW_scale, work_info.surface_point_density)


chunksize = gr.xyzgrid.shape[0] // size
start = rank * chunksize
end   = (rank+1) * chunksize
if rank == (size - 1):
    end = gr.xyzgrid.shape[0] 
chunk = gr.obasis.compute_grid_esp_dm(gr.dm, gr.coordinates, gr.fcharges, gr.xyzgrid[start:end,:])

comm.Send(chunk, dest=0, tag=13)
comm.Disconnect()
