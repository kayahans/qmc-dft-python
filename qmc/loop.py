#!/usr/bin/env python

from numpy import max, min, array, place, random
from numpy.random import rand, randn, seed, uniform, normal
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
    
class vmc():
    def __init__(self, nstep = 20, nblock = 10, nwalkers=1000):
        self.eq        = False
        self.nstep     = nstep
        self.nwalkers  = nwalkers
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

    def loop(self):
        #Start loop
        eave   = 0
        esqu   = 0
        accept = 0
        nstep    = self.nstep
        nwalkers = self.nwalkers
        p        = self.particles
        nelect   = p.nelect
        X        = p.pos
        delta    = self.d
        
        for istep in range(nstep):
            for k in range(nwalkers):
                #Call Psi
                for i in range(nelect):
                    Y = zeros(3)
                    for j in range(ndim-1):
                        Y[j] = X[i, j, k] + delta*(uniform()-0.5) #Moves are random here!!!
                    #end for
                    #Call newPsi
                    #Acceptance prob.
                    q = min(PsiY**2/Psi**X, 1)
                    
                    if q >= rand():
                        #Call copy X(k), Y
                        accept = accept+1/nelect
                    #end if
                #end for
                #Ex = elocal(X(k))
                eave = eave + Ex
                esqu = esqu + Ex**2
            #end for
        #end for
        aratio = accept/nwalkers/nstep
        emean  = eave/nwalkers/nstep
        esigma = sqrt(esqu/nwalker/nstep-emean**2)/(sqrt(nwalker*nstep-1))
    #end def loop
    
    def elocal(self,Xk):
        return 0
    #end def

    def callPsi(self,Xk):
        PsiX = 0
        return PsiX
    #end def

    def newPsi(self,Y,k):
        PsiY = 0
        return PsiY
    #end def

    def copy(self, X, k, Y):
        X[k] = Y
    #end def
    
        
   

for i in range(n):
    x[i] = random.rand(1)
