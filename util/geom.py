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
        self.pos = array(pos)
        self.w   = float(w)   
    #end def __init__

    def error(self, str):
        print str + '\n'
        exit()
    #end def error
        
#end class Particle

class Atom(Particle):
    def __init__(self, chg=1.0, w=1.0, pos=[0,0,0]):
        super(Atom, self).__init__(chg=chg, w=w, pos=pos)
        self.r = radius[chg]
        self.units = 'B'
    #end def
#end class

class Molecule(obj):
    def __init__(self, atoms=[],**kwargs):
        self.atoms  = []
        self.charge = kwargs.get('charge', 0)
        self.mult   = kwargs.get('multiplicity')
        self.name   = kwargs.get('name', 'a molecule')
        self.units  = 'B'

        ne = -self.charge
        
        for at in atoms:
            self.atoms.append(Atom(chg=at.chg,pos=at.pos,w=at.w))
            ne += at.chg
        #end for

        self.ne = ne
        
        if self.mult == None:
            self.set_mult()
        #end if

        self.nn = self.get_nn()
        self.box = self.bbox()
        
    def set_mult(self):
        if self.ne % 2 == 0:
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
        pos = []
        for at in self.atoms:
            pos.append(at.pos)
        #end for
        minxyz = amin(pos, axis=0)
        maxxyz = amax(pos, axis=0)

        boxmin = minxyz - pad
        boxmax = maxxyz + pad
        o = obj()
        o.min = boxmin
        o.max = boxmax
        return o
    #end def

#end class
h   = Atom(chg=1, w=1, pos=[0,0,0])
h2  = Molecule(atoms=[Atom(chg=1, w=1, pos=[-0.34,0,0]), Atom(chg=1, w=1, pos=[0.34,0,0])])
h2p1= Molecule(atoms=[Atom(chg=1, w=1, pos=[-0.34,0,0]), Atom(chg=1, w=1, pos=[0.34,0,0])], charge=1)
             

        

if __name__ == '__main__':

    pdb.set_trace()
    
