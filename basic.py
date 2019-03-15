#!/usr/bin/env python

from numpy import max, min, array, place, random, sqrt
from numpy.random import rand, randn, seed
from itertools import combinations
from numpy.linalg import norm
import pdb

from atoms import Atom
from box import Box
from stats import Stats

def metropolis(n, nstep):
    seed(120398127)
    x = []
    for i in range(n):
        x.append(rand())
    #end for
    fave = 0
    fsqu = 0
    for i in range(nstep):
        for j in range(n):
            # Move particles
            y = x[i] + delta*(rand()-0.5))

            # Acceptance prob.
            # pstar:: 
            q = min(pstar(y)/pstar(x[i]), 1.)

            if q >= rand():
                x[i] = y
            #end if
            
            # Calculate 
            fave = fave + f(x[i])
            fsqu = fsqu + f(x[i])**2
        #end for
    #end for

    fmean = fave/nstep/n
    sigma = sqrt(fsqu/nstep/n - fmean**2)/sqrt(nstep*n - 1)

#end def metropolis

def gmetropolis(n, nstep):
    seed(120398127)
    x = []
    for i in range(n):
        x.append(rand())
    #end for
    fave = 0
    fsqu = 0
    for i in range(nstep):
        for j in range(n):
            # Move particles
            y = x[i] + tran

            #Acceptance prob.
            # G:: transition function (green's function?)
            q = G(x[i], y, t)*pstar(y)/G(y, x[i], t)/pstar(x[i])
            b = min(q, 1.)

            if 
    
                
            
