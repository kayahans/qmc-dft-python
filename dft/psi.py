#!/usr/bin/env python

from numpy import ravel, sqrt

from util.basic   import obj
from util.geom    import Molecule
from grid         import Grid

import sys
sys.dont_write_bytecode = True

import pdb

class Density(obj):
    def __init__(self,d):
        self.d      = d
    #end def

    def set_density(self,rho):
        self.d = rho
    #end def

    def get_density(self):
        return self.d
    #end def
        
    
#end class

class Wavefunction(obj):
    def __init__(self, psi, ham):
        self.ham    = ham
        self.psi    = self.set_psi(psi)
        self.nelect = ham.system.nelect
        self.d      = Density(self.nelect*self.psi**2)
    #end def

    def get_density(self):
        return self.d    
    #end def

    def set_density(self, rho):
        d = self.d
        d.set_density(rho)
    #end def

    def update_density(self):
        d = self.d
        d.set_density(self.nelect*self.psi**2)
    #end def
    
    def get_psi(self):
        return self.psi
    #end def

    def set_psi(self, psi):
        dv = self.ham.grid.dv
        self.psi = ravel(psi)/sqrt(dv)
        return self.psi
    #end def
    
    def update(self, psi):
        self.psi = psi
        self.d   = self.update_density()
    #end def

#end class

    
        
