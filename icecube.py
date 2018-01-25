#!/usr/bin/env python
import sys
import grid
import cube
import numpy as np
import time
from mpi4py import MPI

if __name__ == "__main__":
    if len(sys.argv) < 2:
        if MPI.COMM_WORLD.Get_rank() == 0:
            print("Usage: icecube [qmfile]")
        exit()
    qmfilename = sys.argv[1] 

    #write density
    gr = grid.griddata(qmfilename)
    gr.compute_density()
    if MPI.COMM_WORLD.Get_rank() == 0:
        cube.write_cube("{}_rho.cube".format(qmfilename), gr)
