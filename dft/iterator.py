#!/usr/bin/env pythonw

from numpy import array, meshgrid, linspace, ravel, ones, ceil
from scipy.sparse.linalg import eigsh

from util.basic import obj
from util.geom  import Molecule
from psi        import Wavefunction
from copy import copy

import sys
sys.dont_write_bytecode = True

import pdb

class Iterator(obj):
    def __init__(self, ham, etol=10**-6):
        self.etol   = etol
        self.ham    = ham
    #end def

    def solve(self):
        nelect = self.ham.system.nelect
        
        density_update = False
        if hasattr(self.ham, 'hartree') or hasattr(self.ham, 'xc'):
            density_update = True
            print 'e, e_T,  e_ext, e_hartree, e_xc, ediff'
        #end if

        T    = self.ham.k
        Vext = self.ham.vext
        V_nn = self.ham.vext.nn
        Vtot = T.v + Vext.v
        
        i = 0
        Eprev = 0
        ediff = 10**6
        psi = None
        while ediff > self.etol:
            E, psi = eigsh(Vtot, k=int(ceil(nelect*1.0/2)), which='SA')

            if density_update:
                i += 1
                psi = Wavefunction(psi, self.ham)
                rho = psi.get_density()

                hartree = self.ham.hartree
                xc      = self.ham.xc

                hartree.update(rho)
                xc.update(rho)

                e_T       = T.get_ek(psi)
                e_ext     = Vext.get_e_ne(rho)
                e_hartree = hartree.get_eh(rho)
                e_xc      = xc.get_e_vxc(rho).e_x + xc.get_e_vxc(rho).e_c
                E         = e_T + e_ext + e_hartree + e_xc
                ediff     = abs(Eprev-E)
                Eprev     = E
            
                print 'Iteration #'+str(i), E, e_T,  e_ext, e_hartree, e_xc, ediff
                Vtot = hartree.v + xc.v + T.v + Vext.v
               
            else:
                psi = Wavefunction(psi, self.ham)
                E     = E[0]
                ediff = 0
            #end if
        self.psi = psi
        self.d   = psi.get_density()
            
        #end while

        
        print '\nFinal energy: ', E + V_nn, '\n'
