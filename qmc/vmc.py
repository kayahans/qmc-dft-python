#!/usr/bin/env python

## Algorithm 2.1

from numpy import max, min, array, place, random, zeros
from numpy.random import rand, randn, seed, uniform, normal
from itertools import combinations
from numpy.linalg import norm
import pdb

def vmc():
    nwalkers = 1
    nelect   = 1
    min_max  =(-5., 5.) #Cube box
    delta = 0.1
    ndim = 3
    
    X = zeros((nelect, 3, nwalkers)) #particles initialize
    for k in range(nwalkers):
        for i in range(nelect):
            for j in range(ndim-1):
                X[i, j, k] = uniform(min_max[0],min_max[1])
            #end for
        #end for
    #end for

    #Start loop
    eave   = 0
    esqu   = 0
    accept = 0
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
            #Ex = Elocal(X(k))
            eave = eave + Ex
            esqu = esqu + Ex**2
        #end for
    #end for

    aratio = accept/nwalkers/nstep
    emean  = eave/nwalkers/nstep
    esigma = sqrt(esqu/nwalker/nstep-emean**2)/(sqrt(nwalker*nstep-1))
#end def vmc

    
