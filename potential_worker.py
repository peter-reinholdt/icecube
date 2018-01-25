#!/usr/bin/env python

from mpi4py import MPI
import numpy as np
import horton
import sys
import grid

#get these obasis, dm, xyzgrid

comm = MPI.Comm.Get_parent()
size = comm.Get_size()
rank = comm.Get_rank()

print(rank, size)
qmfilename = comm.recv(source=0, tag=11)
gr = grid.griddata(qmfilename)


chunksize = gr.xyzgrid.shape[0] // size
start = rank * chunksize
end   = (rank+1) * chunksize
if rank == (size - 1):
    end = gr.xyzgrid.shape[0] 
chunk = gr.obasis.compute_grid_esp_dm(gr.dm, gr.coordinates, gr.fcharges, gr.xyzgrid[start:end,:])
comm.Send(chunk, dest=0, tag=13)

comm.Disconnect()
