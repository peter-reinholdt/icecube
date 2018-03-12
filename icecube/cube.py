#!/usr/bin/env python


import numpy as np


def get_xyzgrid_cube(origin, Nx, Ny, Nz, xx, xy, xz, yx, yy, yz, zx, zy, zz):
    xyzgrid = np.zeros((Nx*Ny*Nz,3), dtype=np.float64)
    ki = np.array([xx,xy,xz])
    kj = np.array([yx,yy,yz])
    kk = np.array([zx,zy,zz])
    print(ki, kj, kk)
    counter = 0
    for i in range(Nx):
        for j in range(Ny):
            for k in range(Nz):
                point = origin + i*ki + j*kj + k*kk 
                xyzgrid[counter, 0:3] = point
                counter += 1
    return xyzgrid
