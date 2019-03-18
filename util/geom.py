#!/usr/bin/env pythonw

from numpy import array, place, array2string, amin, amax
from itertools import combinations
from numpy.linalg import norm
from numpy.random import randn
from copy import copy, deepcopy
from collections import OrderedDict
from basic import obj
from elements import radius
import pdb
import sys

sys.dont_write_bytecode = True

class Particle(obj):
    def __init__(self,chg=0,pos=[0,0,0], w=0.):
        self.chg = int(chg)
        self.pos = tuple(pos)
        self.w   = float(w)   
    #end def __init__

    def error(self, str):
        print str + '\n'
        exit()
    #end def error
        
#end class Particle

class Element(Particle):
    def __init__(self, chg=1.0, w=1.0, pos=[0,0,0]):
        super(Element, self).__init__(chg=chg, w=w, pos=pos)
        self.r = radius[chg]
        self.units = 'B'
        self.chg = self.chg
        self.pos = self.pos
        self.w   = self.w
    #end def

    def bbox(self, pad = 5.):
        boxmax = array([pad, pad, pad])
        boxmin = -boxmax
        o = obj()
        o.min = boxmin
        o.max = boxmax       
        return o 
#end class

class Molecule(Element):
    def __init__(self, atoms=[],**kwargs):
        self.natoms = len(atoms)
        self.pos    = []
        self.chg    = []
        self.w      = []
        self.tot_chg= kwargs.get('charge', 0)
        self.mult   = kwargs.get('multiplicity')
        self.name   = kwargs.get('name', 'a molecule')
        self.units  = 'B' 
        nelect = -self.tot_chg
        
        for at in atoms:
            self.pos.append(at.pos)
            self.chg.append(at.chg)
            self.w.append(at.w)
            nelect += at.chg
        #end for

        self.nelect = nelect
        
        if self.mult == None:
            self.set_mult()
        #end if

        self.nn = self.get_nn()
        #self.box = self.bbox()
        
    def set_mult(self):
        if self.nelect % 2 == 0:
            self.mult = 1
        else:
            self.mult = 2
        #end if
    #end def

    def get_nn(self):
        ''' Get nuclear repulsion'''
        return None
    #end def

    def bbox(self, pad = 5.):
        minxyz = amin(self.pos, axis=0)
        maxxyz = amax(self.pos, axis=0)

        boxmin = minxyz - pad
        boxmax = maxxyz + pad
        o = obj()
        o.min = boxmin
        o.max = boxmax
        return o
    #end def

#end class
h   = Molecule(atoms=[Element(chg=1, w=1, pos=[0,0,0])])
h2  = Molecule(atoms=[Element(chg=1, w=1, pos=[-0.34,0,0]), Element(chg=1, w=1, pos=[0.34,0,0])])
h2p1= Molecule(atoms=[Element(chg=1, w=1, pos=[-0.34,0,0]), Element(chg=1, w=1, pos=[0.34,0,0])], charge=1)
he  = Molecule(atoms=[Element(chg=2, w=4, pos=[0,0,0])])         

        

if __name__ == '__main__':
    print h.__repr__()
    pdb.set_trace()
    
