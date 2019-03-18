#!/usr/bin/env python

from numpy import max, min, array, place, random, sqrt, zeros, mean, pi, exp, std
from numpy.random import rand, randn, seed, uniform, normal
from itertools import combinations
from numpy.linalg import norm
import pdb
from collections import OrderedDict
import scipy.integrate as integrate
from scipy import stats
import sys

sys.dont_write_bytecode = True

class obj(object):
    "General all-purpose object"
    def __repr__(self):
        attr = self.__dict__
        attr = OrderedDict(sorted(attr.items()))
        c = ''
        for i in attr.keys():
            s = str(type(attr[i]))
            s = s.rsplit('\'', 2)[1]
            c+= '\n'+ i+' \t\t' + s 
        #end for
        return c[1:] #remove heading whitespace
    #end def __str__

    def error(self, text):
        print 'ERROR: ', str(text)
        exit()
    #end def

    def warning(self, text):
        print 'WARNING: ', str(text)
    #end def
    
#end class obj

def quad(x):
    x = array(x)
    y = x**2
    return y
#end def

def gauss(x, sigma=100, mu=0, check=False):
    x = array(x)
    y = 1./(2*pi*sigma**2)*exp(-(x-mu)**2)/(2*sigma**2)
    if check:
        int = integrate.quad(lambda x:  1./(2*pi*sigma**2)*exp(-(x-mu)**2)/(2*sigma**2), -5, 5)
        print 'Gaussian integral from scipy.quad: '+'\t'+str(int[0])+' +/- ' + str(int[1])
    #end if
    return y
#end def

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
            y = x[i] + delta*(rand()-0.5)

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

    return fmean, sigma

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

            if b >= rand():
                x[i] = y
            #end if
            fave = fave + f(x[i])
            fsqu = fsqu + f(x[i])**2
        #end for
    #end for

    fmean = fave/nstep/n
    sigma = sqrt(fsqu/nstep/n - fmean**2)/sqrt(nstep*n - 1)

    return fmean, sigma
#end def gmetropolis

def uniform_sampling(f, min_max=None, nsample=1000):
    b = min_max[1]
    a = min_max[0]
    X = uniform(a, b, size=nsample)
    fmean = (b-a)/nsample*sum(f(X))
    
    #Compare to importance sampling
    ##h = lambda x : 1.0 #Use later for expectation values
    #p = lambda x : stats.uniform.pdf(x, loc=a, scale=b-a) 
    #g = lambda x : f(x)/p(x)
    #print 1./nsample*sum(h(X)*g(X))
    
    return fmean
#end def uniform_sampling

def importance_sampling(f, min_max=None, scale = None, nsample=1000):
    b = min_max[1]
    a = min_max[0]
    mean = (a+b)/2

    X = normal(mean,scale=scale,size=nsample)
    h = lambda x : 1.0  #Use later for expectation values
    p = lambda x : stats.norm.pdf(x, loc=mean, scale=scale)
    g = lambda x : f(x)/p(x)
    fmean = 1./nsample*sum(h(X)*g(X))

    return fmean       
#end def importance_sampling


if __name__ == '__main__':
    print '='*80
    print 'Uniform_sampling vs. Importance_sampling'
    print '='*80
    mu = []; mi=[]
    for i in range(100):
        mu.append(uniform_sampling(gauss, min_max=(-5., 5.), nsample=100))
        mi.append(importance_sampling(gauss, min_max=(-5., 5.), scale = 1, nsample=100))
    #end for

    gauss(0, check=True)
    print 'Uniform sampling: ', '\t\t\t', mean(mu), '+/-', std(mu)
    print 'Importance samplling: ', '\t\t\t', mean(mi), '+/-', std(mi)
