#!/usr/bin/env python

from util.basic import obj

from interactions import Kinetic, Coulomb, Hartree, XC

import sys
sys.dont_write_bytecode = True

class Hamiltonian(obj):
    def __init__(self, system, grid, xc=None):
        self.system = system
        self.grid   = grid
        self.k      = Kinetic(system, grid)
        self.vext   = Coulomb(system, grid)
        if xc is not None and self.system.nelect > 1:
            self.hartree = Hartree(system, grid)
            self.xc      = XC(system, grid,name=xc)
        #end def
        
if __name__ == '__main__':
    from interactions import Coulomb, Kinetic
    from util.geom    import h, he
    from grid         import Grid
    from iterator     import Iterator

    import pdb
    
    system = h
    grid   = Grid(h, 30, pad=3.)
    
    ham = Hamiltonian(system, grid)

    it = Iterator(ham)
    it.solve()

    system = he
    grid   = Grid(h, 30, pad=3.)
    
    ham = Hamiltonian(system, grid, xc='ldax')

    it = Iterator(ham)
    it.solve()
