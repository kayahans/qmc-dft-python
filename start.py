#!/usr/bin/env python -i

from numpy import max, min, array, place
from itertools import combinations
from numpy.linalg import norm
import pdb

class Atom():
    def __init__(self, p=None, w=None, pos=None):
        self.p   = p
        self.w   = w
        self.pos = array(pos)
    #end def __init__

    def __repr__(self):
        c =  'Atomic number = ' +  str(self.p) + '\n'
        c += 'Weight        = ' +  str(self.w) + '\n'
        c += 'Position      = ' +  str(self.pos)
        return c
    #end def __repr__
#end class Atom

class Box():
    def __init__(self):
        self.origin = array([-10,-10,-10])
        self.corner = array([10,10,10])
        self.pad    = 5
    #end def __init__

    def __repr__(self):
        c =  'Origin = ' +  str(self.origin) + '\n'
        c += 'Corner = ' +  str(self.corner) + '\n'
        c += 'Pad    = ' +  str(self.pad)
        return c
    #end def __repr__
    
    def from_atoms(self, atoms):
        origin = array([1e6, 1e6, 1e6])
        corner = array([-1e6, -1e6, -1e6])
        for a in atoms:
            if isinstance(a, Atom):
                if a.pos[0] <= origin[0]:
                    origin[0] = a.pos[0]
                #end if
                if a.pos[1] <= origin[1]:
                    origin[1] = a.pos[1]
                #end if
                if a.pos[2] <= origin[2]:
                    origin[2] = a.pos[2]
                #end if
                
                if a.pos[0] >= corner[0]:
                    corner[0] = a.pos[0]
                #end if
                if a.pos[1] >= corner[1]:
                    corner[1] = a.pos[1]
                #end if
                if a.pos[2] >= corner[2]:
                    corner[2] = a.pos[2]
                #end if
            #end if
        #end for
        origin -= self.pad
        corner += self.pad
        self.origin = origin
        self.corner = corner
   #def from_atoms
#end class Box

class H():
    def __init__(self):
        self.system = None
        self.h      = None
        self.phi    = None
    #end def
    def __repr__(self):
        c =  'System = ' +  str(self.system) + '\n'
        c += 'Ham = ' +  str(self.h) + '\n'
        c += 'Phi    = ' +  str(self.phi)
        return c
    #end def __repr__
    def add_system(self, system):
        #h = -1/2*Lap(phi)
        v = 0

        if not isinstance(system, list):
            self.error('System must be a list')
        #end if
        
        if isinstance(system[0], Atom):
            atoms = system
            if len(atoms) > 1:
                for a1, a2 in combinations(atoms, 2):
                    v+=(a1.p*a2.p)/(norm(a1.pos-a2.pos))
                #end for
            else:
                v = 0
            #end if
        else:
            self.error('Ion-ion potential can\'t be added')
        #end if
            
        self.h=v

    #end def add_system

    def error(str):
        print str
        exit()
#end class H
                        
                
h = Atom(1, 1, [0,0,0])
h2 = Atom(1,1,[0,0,1])
s = [h, h2]
b = Box()
b.from_atoms(s)
ham=H()
ham.add_system(s)
