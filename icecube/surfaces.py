#!/usr/bin/env python

import numpy as np
from numba import jit, boolean
from utilities import vdw_radii

@jit('float64(float64[:], float64[:])', nopython=True, cache=True)
def dist(x,y):
    """
    Compute distance between vectors x,y
    Assumes that x,y have the same length
    """
    result = 0.0
    for i in range(x.shape[0]):
        result += (x[i] - y[i])**2
    return np.sqrt(result)


@jit(nopython=True, cache=True)
def compute_vdW_surface(atomic_charges, coordinates, surface_point_density=5.0, surface_vdW_scale=1.4):
    """
    Generates apparent uniformly spaced points on a vdw_radii
    surface of a molecule.
    
    vdw_radii   = van der Waals radius of atoms
    points      = number of points on a sphere around each atom
    grid        = output points in x, y, z
    idx         = used to keep track of index in grid, when generating 
                  initial points
    density     = points per area on a surface
    chkrm       = (checkremove) used to keep track in index when 
                  removing points
    """
    
    natoms = coordinates.shape[0]
    points = np.zeros(natoms, dtype=np.int64)
    for i in range(natoms):
        #        area of sphere is               4*pi    *        r**2
        points[i] = np.int(surface_point_density*4*np.pi*(surface_vdW_scale*vdw_radii[atomic_charges[i]])**2)
    # grid = [x, y, z]
    grid = np.zeros((np.sum(points), 3), dtype=np.float64)
    idx = 0
    for i in range(natoms):
        N = points[i]
        #Saff & Kuijlaars algorithm
        for k in range(N):
            h = -1.0 +2.0*k/(N-1)
            theta = np.arccos(h)
            if k == 0 or k == (N-1):
                phi = 0.0
            else:
                #phi_k  phi_{k-1}
                phi = ((phi + 3.6/np.sqrt(N*(1-h**2)))) % (2*np.pi)
            x = surface_vdW_scale*vdw_radii[atomic_charges[i]]*np.cos(phi)*np.sin(theta)
            y = surface_vdW_scale*vdw_radii[atomic_charges[i]]*np.sin(phi)*np.sin(theta)
            z = surface_vdW_scale*vdw_radii[atomic_charges[i]]*np.cos(theta)
            grid[idx, 0] = x + coordinates[i,0]
            grid[idx, 1] = y + coordinates[i,1]
            grid[idx, 2] = z + coordinates[i,2]
            idx += 1
            
   
    
    #This is the distance points have to be apart
    #since they are from the same atom
    grid_spacing = dist(grid[1,:], grid[2,:])
    
    #Remove overlap all points to close to any atom
    not_near_atom = np.ones(grid.shape[0], dtype=boolean)
    for i in range(natoms):
        for j in range(grid.shape[0]):
            r = dist(grid[j,:], coordinates[i,:])
            if r < surface_vdW_scale*0.99*vdw_radii[atomic_charges[i]]:
                not_near_atom[j] = False
    grid = grid[not_near_atom]

    # Double loop over grid to remove close lying points
    not_overlapping = np.ones(grid.shape[0], dtype=boolean)
    for i in range(grid.shape[0]):
        for j in range(i+1, grid.shape[0]):
            if (not not_overlapping[j]): continue #already marked for removal
            r = dist(grid[i,:], grid[j,:])
            if 0.80 * grid_spacing > r:
                not_overlapping[j] = False
    grid = grid[not_overlapping]
    return grid
