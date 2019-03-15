#!/usr/bin/env python

from numpy import max, min, array, place, random
from numpy.random import randn
from itertools import combinations
from numpy.linalg import norm
import pdb

from atoms import Atom
from box import Box
from stats import Stats

class Particles():
    def __init__(self, n = 1):
        self.n = n
        self.box = Box()
        self.density = Density()
        self.p = generate_particles()
    #end def
    def set_particles(self, n):
        self.n = n
    #end def
    def generate_particles():
        n = self.n
        box = self.box
        return random(n)
    
class vmc_loop():
    def __init__(self, nstep = 20, nblock = 10):
        self.eq        = False
        self.nstep     = nstep
        self.block     = nblock
        self.particles = Particles()
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
    
    def evaluate_energy(self):
        for nb in range(self.block):
            for ns in range(self.nstep):
                for p in self.particles:
                    if p.chg < 0.0:
                        pnew = move(p)
                        

for i in range(n):
    x[i] = random.rand(1)
