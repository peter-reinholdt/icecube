#!mpirun python
import sys
import grid
import cube
import numpy as np
import time
import argparse
import classic

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Make cube files from commom QM restart formats')
    parser.add_argument('-qm', '--qm-file', dest='qmfile', help='Name of the QM restart file to process (.molden/.fchk/...)')
    parser.add_argument('-cl', '--classic-file', dest='classicfile', default=None, help='Name of file containing classical parameters (coordinates, charges, dipoles)')
    parser.add_argument('-N',  '--nprocs', dest='nprocs', default=1, type=int, help='Number of MPI processes to use')
    parser.add_argument('--density', dest='do_density', action='store_true', help='Request density cube file (true/false).')
    parser.add_argument('--potential', dest='do_potential', action='store_true', help='Request potential cube file (true/false).')
    parser.add_argument('--cube-density', dest='cube-density', default=1.0, help='Points/bohr to output to cube file.')

    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        exit()
    
    if args.qmfile is not None:
        gr = grid.griddata(args.qmfile)
    else:
        print("Missing QM file")
        exit()
    if args.classicfile:
        data = np.loadtxt(args.classicfile, dtype=np.float64)
        coordinates = data[:,0:3]
        charges     = data[:,3]
        dipoles     = data[:,4:7]
    
    if args.do_density:
        print("Computing the density...")
        gr.compute_density(nprocs=args.nprocs)
        print("...Finished!")
        print("Writing data to cube...")
        cube.write_cube("{}_rho.cube".format(args.qmfile), gr)
        print("...Finished!")
   
    if args.do_potential:
        print("Computing the QM ESP...")
        gr.compute_potential(nprocs=args.nprocs)
        print("...Finished!")
        print("Writing data to cube...")
        cube.write_cube("{}_QMESP.cube".format(args.qmfile), gr)
        print("...Finished!")
        if args.classicfile:
            print("Computing contribution due to classical point charges...")
            esp = charge_potential(charges, coordinates, gr.xyzgrid)
            gr.data -= esp
            print("...Finished!")
            print("Writing data to cube...")
            cube.write_cube("{}_QMESP_q.cube".format(args.qmfile), gr)

            print("Computing contribution due to classical point dipoles...")
            esp = dipole_potential(dipoles, coordinates, gr.xyzgrid)
            gr.data -= esp
            print("...Finished!")
            print("Writing data to cube...")
            cube.write_cube("{}_QMESP_qmu.cube".format(args.qmfile), gr)
