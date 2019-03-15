#!/usr/bin/env python

from numpy import max, min, array, place, array2string
from itertools import combinations
from numpy.linalg import norm
from numpy.random import randn
from copy import copy, deepcopy
from collections import OrderedDict
from qmc_basic import obj
import pdb


class Particle(obj):
    def __init__(self,chg=0,pos=[0,0,0], w=0):
        self.chg = chg
        self.pos = array(pos)
        self.w   = w   
    #end def __init__

    def error(self, str):
        print str + '\n'
        exit()
    #end def error
        
#end class Atom


class Structure(obj):
    def __init__(self, center=None, pos=None, chg=None, name = ''):
        self.center = center
        self.pos    = pos
        if chg is None and pos is not None:
            self.chg = array([1.0]*len(pos))
        else:
            self.chg = chg
        #end if
        
        if name == 'h':
            self.h()
        elif name == 'h2+':
            self.h2()
        #end def
    #end def

    def h(self):
        self.center = array([0,0,0])
        self.pos    = array([[0,0,0]])
        self.ch     = array([1.0])
    #end def

    def h2(self):
        self.center = array([0,0,0])
        self.pos    = array([[-0.34,0,0], [0.34,0,0]])
        self.chg    = array([1.0, 1.0])
    #end def

#end class

        

if __name__ == '__main__':
    aa = Structure(name='h')
    pdb.set_trace()
    
