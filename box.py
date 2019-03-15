#!/usr/bin/env python

from numpy import max, min, array, place
from itertools import combinations
from numpy.linalg import norm
import pdb

class Box():
    def __init__(self):
        self.origin   = array([-10,-10,-10])
        self.corner   = array([10,10,10])
        self.pad      = 5
        self.periodicity = [False, False, False]
    #end def __init__

    def __repr__(self):
        c =  'Origin = ' +  str(self.origin) + '\n'
        c += 'Corner = ' +  str(self.corner) + '\n'
        c += 'Pad    = ' +  str(self.pad) + '\n'
        c += 'Periodicity ' + str(self.periodicity) + '\n'
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
    #end def from_atoms

    def from_len(self, len):
       self.origin = -len*1.0/2
       self.corner = len*1.0/2
    #end def from_len
    
#end class Box
