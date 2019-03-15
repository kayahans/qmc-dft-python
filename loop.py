#!/usr/bin/env python

from numpy import max, min, array, place, random
from numpy.random import randn
from itertools import combinations
from numpy.linalg import norm
import pdb

from atoms import Atom
from box import Box
from stats import Stats

class vmc_loop():
    def __init__(self, nstep = 20, nblock = 10):
        self.eq        = False
        self.nstep     = nstep
        self.block     = nblock
        self.particles = []
        self.d         = 0.2
        self.box       = None
        self.wave_func = None
        self.stats     = Stats()        
    #end def
    
    def error(self, str):
        print str + '\n'
        exit()
    #end def error

    def warning(self, str):
        print 'Warning: ' + str + '\n'
    #end def warning
    
    def move(self, p):
        pnew = p.copy()
        d = self.d
        newpos = p.pos + d*(randn(3)-0.5)
        if (newpos < self.box.origin).any() or (newpos > self.box.corner).any():
            move(p, d)
        else:
            pnew.pos = newpos
            return pnew
        #end if
    #end def move
    
    def evaluate_energy(self, eq=False):
        for nb in range(self.block):
            for ns in range(self.nstep):
                for p in self.particles:
                    if p.chg < 0.0:
                        pnew = move(p, 
