#!/usr/bin/env python

from __future__ import division
from numba import jit
import numpy as np


@jit(nopython=True)
def charge_potential(charges, coordinates, gridpoints):
    """
    Evaluate potential of charges on specified gridpoints
    charges     : 1D array of charges              shape=(ncharges,)
    coordinates : positions of the charges,        shape=(ncharges,3)
    gridpoints  : where to evaluate the potential, shape=(ngridpoints,3)
    """
    ncharges    = coordinates.shape[0]
    ngridpoints = gridpoints.shape[0]
    grid        = np.zeros(ngridpoints, dtype=np.float64)
    for i in range(ncharges):
        ri = coordinates[i,:]
        for j in range(ngridpoints):
            rj = gridpoints[j,:]
            rij = ri - rj
            grid[j] += charges[i] / np.sqrt(np.dot(rij,rij))
    return grid


@jit(nopython=True)
def dipole_potential(dipoles, coordinates, gridpoints):
    """
    Evaluate potential of dipoles on specified gridpoints
    dipoles     : 2D array of dipoles              shape=(ndipoles,3)
    coordinates : positions of the dipoles,        shape=(ndipoles,3)
    gridpoints  : where to evaluate the potential, shape=(ngridpoints,3)
    """
    ndipoles    = coordinates.shape[0]
    ngridpoints = gridpoints.shape[0]
    grid        = np.zeros(ngridpoints, dtype=np.float64)
    for i in range(ndipoles):
        ri = coordinates[i,:]
        for j in range(ngridpoints):
            rj = gridpoints[j,:]
            rij = ri - rj
            rij_inv_cubed = (1.0/np.sqrt(np.dot(rij,rij)))**3 #1/|r|**3
            for k in range(3):
                grid[j] -= dipoles[i,k] * rij[k] * rij_inv_cubed #-mu_x * x / |r|**3
    return grid
