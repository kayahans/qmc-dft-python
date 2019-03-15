#!/usr/bin/env python
from itertools import combinations
from numpy.linalg import norm
from particles import Particle

class Coulomb():
    def __init__(self, particles=[]):
        if particles is not []:
            for particle in particles:
                if not isinstance(particle, Particle):
                    self.error('{0} must be an instance of Particle()'.format(particle))
                #end if
            #end for
        #end if
        self.ee        = 0.0
        self.en        = 0.0
        self.nn        = 0.0
        self.e         = 0.0
    #end def __init__
    
    def v(self,particles):
        for i,j in itertools.combinations(particles, 2):
            d = abs(norm(i.pos-j-pos))
            energy = -i.chg*j*chg/d**2
            self.e += energy
            if i.chg < 0.0:
                if j.chg < 0.0:
                    self.ee += energy
                else:
                    self.en += energy
                #end if
            else:
                if j.chg < 0.0:
                    self.en += energy
                else:
                    self.nn += energy
                #end if
            #end if
    #end def v

    def error(self, str):
        print str + '\n'
        exit()
    #end def error
    
#end class Coulomb
    
