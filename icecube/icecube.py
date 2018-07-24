#!/usr/bin/env python
import sys
import grid
import cube
import numpy as np
import time
import argparse
import classic
import surfaces

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make cube files from commom QM restart formats')
    parser.add_argument('-qm', '--qm-file',         dest='qm_file_name', help='Name of the QM restart file to process (.molden/.fchk/...)')
    parser.add_argument('-cl', '--classic-file',    dest='classicfile', default=None, help='Name of file containing classical parameters (coordinates, charges, dipoles)')
    parser.add_argument('-N',  '--nprocs',          dest='nprocs', default=1, type=int, help='Number of MPI processes to use')
    parser.add_argument('--density',                dest='do_density', action='store_true', help='Request density cube file (true/false).')
    parser.add_argument('--potential',              dest='do_potential', action='store_true', help='Request potential cube file (true/false).')
    parser.add_argument('--cube-density',           dest='cube_density', default=4.0, type=float, help='Points/bohr to output to cube file.')
    parser.add_argument('--cube-buffer',            dest='cube_buffer', default=5.0, type=float, help='Buffer (in bohr) to add around min/max value of coordinates')
    parser.add_argument('--surface-potential-qm',   dest='do_surface_potential_qm', action='store_true', help='Request calculation of QM ESP at molecular vdW surface')
    parser.add_argument('--surface-potential-cl',   dest='do_surface_potential_classic', action='store_true', help='Request calculation of classic ESP at molecular vdW surface')
    parser.add_argument('--surface-vdW-scale',      dest='surface_vdW_scale', default=2.0, type=float, help='Set the vdw radius scale parameter.') 
    parser.add_argument('--surface-point-density',  dest='surface_point_density', default=20.0, type=float, help='Set the vdw surface point density') 
    parser.add_argument('--read-grid',              dest='read_grid', default="", type=str, help='Read a manually generated grid instead of using the internal grid generation')

    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    
    if args.qm_file_name is not None:
        gr = grid.griddata(args.qm_file_name)
    else:
        raise ValueError("Missing QM input file")
    if args.classicfile:
        data = np.loadtxt(args.classicfile, dtype=np.float64)
        coordinates = data[:,0:3]
        charges     = data[:,3]
        dipoles     = data[:,4:7]
    

    if args.do_density:
        print("Computing the density...")
        gr.get_cubic_grid(args.cube_buffer, args.cube_density)
        gr.compute_density(nprocs=args.nprocs)
        print("...Finished!")
        print("Writing data to cube...")
        cube.write_cube("{}_rho.cube".format(args.qm_file_name), gr)
        print("...Finished!")
   

    if args.do_potential:
        print("Computing the QM ESP...")
        gr.get_cubic_grid(args.cube_buffer, args.cube_density)
        gr.compute_potential(nprocs=args.nprocs)
        print("...Finished!")
        print("Writing data to cube...")
        cube.write_cube("{}_QMESP.cube".format(args.qm_file_name), gr)
        print("...Finished!")
        if args.classicfile:
            print("Computing contribution due to classical point charges...")
            esp = classic.charge_potential(charges, coordinates, gr.xyzgrid)
            gr.data -= esp
            print("...Finished!")
            print("Writing data to cube...")
            cube.write_cube("{}_QMESP_q.cube".format(args.qm_file_name), gr)

            print("Computing contribution due to classical point dipoles...")
            esp = classic.dipole_potential(dipoles, coordinates, gr.xyzgrid)
            gr.data -= esp
            print("...Finished!")
            print("Writing data to cube...")
            cube.write_cube("{}_QMESP_qmu.cube".format(args.qm_file_name), gr)


    if args.do_surface_potential_qm or args.do_surface_potential_classic:
        if not args.read_grid:
            print("Computing molecular surface of {}*vdW with a {} point density".format(args.surface_vdW_scale, args.surface_point_density))
            gr.get_vdW_surface(args.surface_vdW_scale, args.surface_point_density)
        else:
            gr.xyzgrid = np.loadtxt(args.read_grid, skiprows=1)
        if args.do_surface_potential_qm:
            print("Computing ESP due to QM density and nuclei at the gridpoints")
            gr.compute_potential(nprocs=args.nprocs)
            with open("{}_{}_{}_qm.dat".format(args.qm_file_name, args.surface_vdW_scale, args.surface_point_density), "w") as f:
                f.write("#Rx,Ry,Rz,QM_ESP(R)\n")
                for i in range(gr.xyzgrid.shape[0]):
                    f.write("{} {} {} {}\n".format(gr.xyzgrid[i,0], gr.xyzgrid[i,1], gr.xyzgrid[i,2], gr.data[i]))
        if args.do_surface_potential_classic and args.classicfile:
            print("Computing ESP due to classic charges at the gridpoints")
            charge_esp = classic.charge_potential(charges, coordinates, gr.xyzgrid)
            with open("{}_{}_{}_charge.dat".format(args.qm_file_name, args.surface_vdW_scale, args.surface_point_density), "w") as f:
                f.write("#Rx,Ry,Rz,charge_ESP(R)\n")
                for i in range(gr.xyzgrid.shape[0]):
                    f.write("{} {} {} {}\n".format(gr.xyzgrid[i,0], gr.xyzgrid[i,1], gr.xyzgrid[i,2], charge_esp[i]))
            print("Computing ESP due to classic dipoles at the same gridpoints")
            dipole_esp = classic.dipole_potential(dipoles, coordinates, gr.xyzgrid)
            with open("{}_{}_{}_dipole.dat".format(args.qm_file_name, args.surface_vdW_scale, args.surface_point_density), "w") as f:
                f.write("#Rx,Ry,Rz,dipole_ESP(R)\n")
                for i in range(gr.xyzgrid.shape[0]):
                    f.write("{} {} {} {}\n".format(gr.xyzgrid[i,0], gr.xyzgrid[i,1], gr.xyzgrid[i,2], dipole_esp[i]))
