#!/usr/bin/env python


from utilities import split_to_float
import numpy as np
import grid


def read_cube(filename, griddata):
    """Loads cube file into a griddata object"""
    with open(filename, "r") as f:
        #real header
        f.readline()
        f.readline()
        #read atoms, origin
        data = split_to_float(f.readline())
        natoms = int(data[0])
        origin = data[1:4]
        #read bodata info
        data = split_to_float(f.readline())
        Nx = int(data[0])
        x,xy,xz = data[1:4]

        data = split_to_float(f.readline())
        Ny = int(data[0])
        yx,yy,yz = data[1:4]

        data = split_to_float(f.readline())
        Nz = int(data[0])
        zx,zy,zz = data[1:4]
        
        icharges = np.zeros(natoms, dtype=np.int64)
        fcharges = np.zeros(natoms, dtype=np.float64)
        coordinates = np.zeros((natoms,3), dtype=np.float64)

        #read atom data
        for i in range(natoms):
            data = split_to_float(f.readline())
            icharges[i] = data[0]
            fcharges[i] = data[1]
            coordinates[i] = data[2:5]
        
        #read in main data
        cube_data = np.array(split_to_float(f.read())).reshape(Nx,Ny,Nz)
    if type(griddata) == grid.griddata:
        griddata.natoms = natoms
        griddata.origin = origin
        griddata.Nx = Nx
        griddata.Ny = Ny
        griddata.Nz = Nz
        griddata.xx = xx
        griddata.xy = xy
        griddata.xz = xz
        griddata.yx = yx
        griddata.yy = yy
        griddata.yz = yz
        griddata.zx = xz
        griddata.zy = zy
        griddata.zz = zz
        griddata.icharges = icharges
        griddata.fcharges = fcharges
        griddata.coordinates = coordinates
        self.grid_type = 'cubic'
        griddata.data = cube_data
    else:
        pass
    return griddata
