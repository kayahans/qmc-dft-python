#!/usr/bin/env pythonw

from util.basic   import obj
from util.geom    import Molecule
from grid         import Grid
from psi          import Density, Wavefunction
from functionals  import XCFun

from numpy import array, abs, pi, sqrt, power, log, arctan
from scipy.sparse.linalg import cgs
from scipy.sparse import spdiags, csr_matrix,csc_matrix
from itertools    import combinations

import sys
sys.dont_write_bytecode = True

class Interaction(obj):
    def __init__(self, system, grid):
        assert isinstance(system, Molecule)
        assert isinstance(grid, Grid)
        self.system = system
        self.grid   = grid
        self.v      = []
        self.e      = 0.
    #end def
#end class

class Coulomb(Interaction):
    def __init__(self, system, grid):
        super(Coulomb, self).__init__(system, grid)

        if isinstance(system, Molecule):
            self.nn  = self.get_e_nn()
            self.v   = self.get_v_ne()
        else:
            self.error('System must be type Molecule!')
        #end if
    #end def

    def get_e_nn(self):
        for pi in self.system.pos:
            nn     = 0.0
            npos   = self.system.pos
            nchg   = self.system.chg
            natoms = len(npos)
            if natoms > 1:
                for i,j in combinations(range(natoms), 2):
                    d = abs(norm(npos[i]-npos[j]))
                    nn += nchg[i]*nchg[j]/d**2
                #end for
            #end if
            return nn
    #end def calc_nn

    def get_v_ne(self):
        grid   = self.grid
        npos   = self.system.pos
        nchg   = self.system.chg
        natoms = len(npos) 
        ngrids = grid.shape[0]*grid.shape[1]*grid.shape[2]
        
        xn, yn, zn = grid.get_gridv()

        v_ne = 0
        for i in range(natoms):
            xi, yi, zi = npos[i]
            z          = nchg[i]
            v_ne      += -z/sqrt((xn-xi)**2+(yn-yi)**2+(zn-zi)**2)
        #end for
        v_ne = spdiags(v_ne, 0, ngrids, ngrids)
        self.v = v_ne
        return v_ne
    #end def

    def get_e_ne(self, density):
        assert isinstance(density, Density)
        density = density.get_density()
        grid = self.grid
        v_ne = self.v
        dv   = grid.dv
            
        e_ne = sum(density*v_ne)*dv
        return e_ne
    #end def

    def update(self,density):
        self.en  = self.get_e_ne(density)
    #end def
    
#end class

class Hartree(Interaction):
    def __init__(self, system, grid, comp_chg=0.):
        super(Hartree, self).__init__(system, grid)
        self.e        = 0
        self.lapn     = grid.get_lapn()
    #end def

    def get_vh(self, density):
        assert isinstance(density, Density)
        density = density.d
        grid = self.grid
        ngrids = grid.shape[0]*grid.shape[1]*grid.shape[2]
        #pdb.set_trace()
        vh     = cgs(self.lapn, -4*pi*(density))[0]
        vh     = spdiags(vh, 0, ngrids, ngrids)
        return vh
    #end def 
            
    def get_eh(self, density):
        assert isinstance(density, Density)
        density = density.d
        grid = self.grid
        dv   = grid.dv
        
        eh = 0.5*sum(density*self.v)*dv
        return eh
    #end def

    def update(self,density):
        self.v = self.get_vh(density)
        self.e = self.get_eh(density)
    #end def
    
#end class

class Kinetic(Interaction):
    def __init__(self, system, grid):
        super(Kinetic, self).__init__(system, grid)
        self.v    = self.get_vk()
        self.e    = 0.
    #end def

    def get_vk(self):
        lapn   = self.grid.get_lapn() #Laplacian
        return -0.5*lapn
    #end def

    def get_ek(self, psi):
        assert isinstance(psi, Wavefunction)
        v      = self.v
        dv     = self.grid.dv
        nelect = self.system.nelect
        psi    = psi.get_psi()

        # T = psi'*Lap3*psi
        T1   = csr_matrix.dot(v,psi)
        T2   = csc_matrix.dot(csc_matrix(psi), T1)
        T    = T2*nelect*dv
        T    = T[0]
        return T
    #end def

    def update(self,psi):
        self.e = get_ek(psi)
    #end def
    
#end class

class XC(Interaction):
    def __init__(self,system,grid,name='ldax'):
        super(XC, self).__init__(system, grid)
        self.name    = name
        self.xcfun   = self.get_xcfun()
    #end def
    
    def get_xcfun(self):
        xcfun = XCFun(name = self.name)
        return xcfun
    #end def

    def get_e_vxc(self, density):
        density = density.d
        grid    = self.grid
        dv      = grid.dv

        xc = self.xcfun.get_xc(density)
        
        ngrids = grid.shape[0]*grid.shape[1]*grid.shape[2]
        
        x = xc.x
        c = xc.c

        o = obj()
        o.v_x = spdiags(x.f, 0, ngrids, ngrids)
        o.v_c = spdiags(c.f, 0, ngrids, ngrids)

        o.e_x = x.F*dv
        o.e_c = c.F*dv
        
        return o
    #end def

    def update(self, density):
        e_v_xc = self.get_e_vxc(density)
        self.e = e_v_xc.e_x + e_v_xc.e_c
        self.v = e_v_xc.v_x + e_v_xc.v_c
    #end def
    
#end class

if __name__ == '__main__':
    from util.geom import h
    system = h
    grid   = Grid(h, 40, 3)
    cc     = Coulomb(system, grid)
    hh     = Hartree(system, grid)
    kk     = Kinetic(system, grid)
    xc     = XC(system, grid)
