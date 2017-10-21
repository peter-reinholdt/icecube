#!/usr/bin/env python


def write_cube(filename, griddata)
    with open(filename, "w") as f:
        #header, comments
        f.write("Cube file generated by icecube\n")
        f.write("Title card\n")
        #header 
        f.write("{} {} {} {} {}\n".format(griddata.natoms, gridata.origin[0], griddata.origin[1], griddata.origin[2]))
        #header, box info. N is number of points, ij is increments. Data for X/Y/Z axis
        f.write("{} {} {} {}\n".format(griddata.Nx, griddata.xx, griddata.xy, griddata.xz))
        f.write("{} {} {} {}\n".format(griddata.Ny, griddata.yy, griddata.yy, griddata.yz))
        f.write("{} {} {} {}\n".format(griddata.Nz, griddata.zy, griddata.zy, griddata.zz))
        #header atoms
        for i in range(griddata.natoms):
            f.write("{} {} {} {} {}\n".format(griddata.icharge[i], griddata.fcharge[i],\
                griddata.coordinates[i,0], griddata.coordinates[i,1], griddata.coordinates[i,2])
 
        #main data
        counter = 0 
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    f.write("{}".format(griddata.data[counter]))
                    counter += 1
                    if counter % 6 == 0:
                        f.write("\n")
                f.write("\n")
