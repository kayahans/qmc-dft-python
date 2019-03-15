#!/usr/bin/env python -i

from numpy import max, min, array, place
from itertools import combinations
from numpy.linalg import norm
from numpy.random import randn
import pdb

class Particle():
    def __init__(self,chg=None,pos=None):
        self.chg = chg
        self.pos = array(pos)
    #end def __init__

    def __repr__(self):
        c =  'Charge = ' +  str(self.p) + '\n'
        c += 'Position      = ' +  str(self.pos)
        return c
    #end def __repr__

    def error(self, str):
        print str + '\n'
        exit()
    #end def error
    
#end class Atom


class Atom(Particle):
    def __init__(self, p=None, e=None, pos=None):
        self.p   = Particle(chg=p, pos=pos)
        self.e   = []
        for i in range(-e):
            self.e.append(Particle(chg=-1, pos=randn(3))

        if p is not None or p < 0.0:
            self.error('Atoms should have positive charge')
        #end if
    #end def __init__

    def __repr__(self):
        c =  'Atomic number       = ' +  str(self.p) + '\n'
        c += 'Number of electrons = ' +  str(self.e) + '\n'
        c += 'Position            = ' +  str(self.pos)
        return c
    #end def __repr__
    
#end class Atom
