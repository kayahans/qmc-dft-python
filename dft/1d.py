#!/usr/bin/env pythonw

from numpy import array, abs, pi, sqrt, power, log, arctan, zeros, diag, ones
from numpy.linalg import eigh
from util.basic import obj
from util.box import Box, Grid
from functionals import Xc
from copy import copy
import pdb
import numpy as np

class Kinetic(obj):
    # 1-D only
    def __init__(self):
        self.T_gg = None
    #end def
    
    def set_matrix(self, grid):
        ng = grid.shape
        T_gg = zeros((ng[0], ng[0])) #1-D
        pdb.set_trace()
        for i in range(ng[0]): #1-D
            T_gg[i,i] = -2.0
            if i > 0:
                T_gg[i, i-1] = 1.0
                T_gg[i-1, i] = 1.0
            #end if
        #end if
        T_gg *= 0.5/grid.dx**2
        self.T = T_gg
    #end def
#end def

class Density(Grid):
    def __init__(self, Nn, grid):
        ng = grid.shape[0] #1-D
        dx = grid.dx
        self.shape = ng
        self.dx    = dx
        self.data  = 2.0 * Nn / (ng * dx) * ones(ng)
    #end def

#end def




def Loop():

    b = Box()
    b.from_len(12., dim=1)
    g = Grid(200, b)
    dx = g.dx
    Nn = 4
    d = Density(Nn, g)
    n_g = d.data
    Ng  = d.shape
    x_g = g.points
    
    vext_g     = .5 * g.points**2
    vhartree_g = zeros(g.shape)
    vx_g       = zeros(g.shape)

    T_gg = np.zeros((Ng, Ng))  # Kinetic operator
    for i in range(Ng):
        T_gg[i, i] = -2.0
        if i > 0:
            T_gg[i, i - 1] = 1.0
            T_gg[i - 1, i] = 1.0
    T_gg *= -0.5 / dx**2

    # Initialize density as even:
    n_g = 2.0 * Nn / (Ng * dx) * np.ones(Ng)
    print('Initial charge', n_g.sum() * dx)
    
    # Nn states, each one doubly occupied.
    # Initialize as constant density:
    vhartree_g = np.zeros(Ng)
    vx_g = np.zeros(Ng)


    def soft_poisson_solve(n_g):
        vhartree_g = np.zeros(Ng)
        for i in range(Ng):
            for j in range(Ng):
                vhartree_g[i] += n_g[j] / np.sqrt(1.0 + (x_g[i] - x_g[j])**2)
        vhartree_g *= dx
        return vhartree_g


    density_change_integral = 1.0
    while density_change_integral > 1e-6:
        # Calculate Hamiltonian
        veff_g = vext_g + vhartree_g + vx_g
        H_gg = T_gg + np.diag(veff_g)  # Hamiltonian
        pdb.set_trace()
        # Solve KS equations
        eps_n, psi_gn = np.linalg.eigh(H_gg)
        print('Energies', ' '.join('{:4f}'.format(eps) for eps in eps_n[:Nn]))

        # Normalize states (states are normalized, but not in our dx metric)
        psi_gn /= np.sqrt(dx)

        # Update density
        nold_g = n_g
        n_g = 2.0 * (psi_gn[:, :Nn]**2).sum(axis=1)
        density_change_integral = np.abs(nold_g - n_g).sum() * dx
        
        charge = n_g.sum() * dx
        print('Number of electrons', charge)
        print('Convergence err', density_change_integral)
        assert abs(charge - 2.0 * Nn) < 1e-14

        # Calculate Hartree potential
        vhartree_g = soft_poisson_solve(n_g)
        Ehartree = 0.5 * (vhartree_g * n_g).sum() * dx
        print('Electrostatic energy', Ehartree)

        # Calculate exchange potential (we won't bother with correlation!)
        vx_g = -(3. / np.pi * n_g)**(1. / 3.)
        Ex = -3. / 4. * (3. / np.pi)**(1. / 3.) * (n_g**(4. / 3.)).sum() * dx
        print('Exchange energy', Ex)

        Eks = 2.0 * eps_n[:Nn].sum()  # Band structure energy
        Ekin = Eks - (veff_g * n_g).sum() * dx
        print('Ekin', Ekin)
        Epot = Ehartree + Ex + (vext_g * n_g).sum() * dx
        print('Epot', Epot)
        Etot = Ekin + Epot
        print('Energy', Etot)
    
if __name__ == '__main__':
    Loop()
