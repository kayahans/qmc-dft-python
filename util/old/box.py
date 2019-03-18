#!/usr/bin/env python

from numpy import max, min, array, place, zeros, linspace
from itertools import combinations
from numpy.linalg import norm
from util.basic import obj
import pdb

class Box(obj):
    def __init__(self):
        self.origin   = array([0,0,0])
        self.corner   = array([0,0,0])
        self.periodicity = [False, False, False]
    #end def __init__
    
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

    def from_len(self, length, dim=None):
       if dim == 3:
           self.origin = array([0., 0., 0.])
           self.corner = -length*array([1., 1., 1.])/2
       elif dim == 1:
           self.origin = array([0.])
           self.corner = -length*array([1.])/2
        #end if
        
    #end def from_len
    
#end class Box


class Grid(obj):
    def __init__(self, shape, box):
        b = box
        diag  = b.origin-2*b.corner
        pdb.set_trace()
        
        self.shape  = [int(shape)]
        self.data   = zeros(self.shape)
        self.points = linspace(b.corner,  b.origin-b.corner, self.shape[0])
        self.dx     = self.points[1]-self.points[0]
    #end def

