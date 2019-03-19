#!/usr/bin/env python


from util.basic import obj
from grid import Grid
from hamiltonian import Hamiltonian
from iterator import Iterator

import sys
sys.dont_write_bytecode = True

class System(obj):
    def __init__(self, structure, **kwargs):
        kw_grid = dict()
        kw_ham  = dict()
        self.structure = structure
        self.grid      = Grid(structure, 30, **kw_grid)
        self.ham       = Hamiltonian(self.structure, self.grid, **kw_ham)
        self.it        = Iterator(self.ham)
        self.solved    = False
    #end def

    def solve(self):
        self.it.solve()
        self.solved = True
    #end def

    def get_density(self):
        from numpy import reshape
        if not self.solved:
            self.solve()
        density = self.it.d.d
        density = reshape(density, self.grid.shape)
        return density
    #end def

    def plot_density(self, isosurface=0.1):
        #incomplete
        from mpl_toolkits import mplot3d
        import matplotlib.pyplot as plt
        fig = plt.figure()
        ax = plt.axes(projection='3d')
        x,y,z = self.grid.get_grid()
        d = self.get_density()
        pdb.set_trace()
        surf = ax.plot_surface(x[0],y[0],d[0])
        plt.show()
    #end def
    
if __name__ == '__main__':
    from util.geom    import h, he
    import pdb
    
    aa = System(h)
    aa.solve()
    dd = aa.get_density()
    aa.plot_density()
    pdb.set_trace()
