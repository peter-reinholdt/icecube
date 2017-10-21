#!/usr/bin/env python
import sys
import grid
import cube
import numpy as np

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: icecube [qmfile]")
        exit()
    qmfilename = sys.argv[1] 

    #write density
    gr = grid.griddata(qmfilename)
    gr.compute_density()
    cube.write_cube("{}_rho.cube".format(qmfilename), gr)

    #write log(density) -- note the base
    gr.data = np.log(gr.data)
    cube.write_cube("{}_logrho.cube".format(qmfilename), gr)

