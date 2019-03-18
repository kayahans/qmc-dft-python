#!/usr/bin/env pythonw

from numpy import array, abs, pi, sqrt,linspace, zeros, exp
from itertools import combinations
from numpy.linalg import norm
import pdb
from math import factorial
from scipy.misc import factorial2
from basic import obj
import matplotlib.pyplot as plt

class Orbital(obj):
    def __init__():
        self.type = None
    #end def __init__
    
    def error(self,str):
        print str + '\n'
        exit()
    #end def error

    def set_type(self, t):
        self.type = t
    #end def set_type
    
#end class Orbital

class STO(Orbital):
    # STO is of the type:
    # R(r) = N*r^{n-1}*e^-{alpha*r} in radial coordinates
    # N is the normalizing factor
    # Eta is related to screening, 1 is optimal 
    # n is the principle quantum number
    # Defaults are optimal for Hydrogen atom
    def __init__(self, alpha=1., n=1, eta=1, center=array([0.0, 0.0, 0.0])):
        self.n      = n
        self.eta    = eta
        self.center = center
        self.normf  = alpha**n/factorial(n-1)
        self.set_type('STO')
    #end def __init__

    # Return value
    def v(self, pos):
        rij   = array(norm(pos-self.center, axis=1))
        normf = self.normf
        n     = self.n
        eta   = self.eta
        
        return normf*rij**(n-1)*exp(-eta*rij)
    #end def v
    
    #Return derivative
    def d(self, pos):
        rij   = array(norm(pos-self.center, axis=1))
        normf = self.normf
        n     = self.n
        eta   = self.eta
        
        return ((n-1)/rij - eta)*self.v(pos)
    #end def d

    #Return Laplacian
    def lap(self, pos):
        rij   = array(norm(pos-self.center, axis=1))
        normf = self.normf
        n     = self.n
        eta   = self.eta
        
        return (n*(n-1)/rij**2 - 2*n*eta/rij+eta**2)*self.v(pos) 
    #end def lap

    def plot_1D(self):
        x = []
        for i in linspace(-3,3,101):
            x.append([i,0.,0.])
        #end for
        x = array(x)
        v = self.v(x)
        d = self.d(x)
        lap = self.lap(x)
        rij = x[:,0]
        plt.plot(rij,v, label='f')
        plt.plot(rij,d, label='f\'')
        plt.plot(rij,lap, label='f\'\'')
        plt.ylim(min(d)-1, max(v)+1)
        plt.axhline(y=0, xmin=-3, xmax=3, linestyle='--', color='k')
        plt.axhline(y=1, xmin=-3, xmax=3, linestyle='--', color='k')
        plt.axvline(x=0, ymin=min(d)-1, ymax=max(v)+1, linestyle='--', color='k')
        plt.legend()
        plt.show()
#end class STO

class GTO(Orbital):
    # Incomplete
    
    # Cartesian GTO is of the type:
    # Phi(x,y,z) = N*x^{a}*y^{b}*z^{c}*e^-{alpha*r^2} where L = a+b+c is the angular momentum
    # l and m follows from nlm notation
    # Used http://www.chem.unifr.ch/cd/lectures/files/module5.pdf
    def __init__(self, l=0, m=0, eta=1, center=array([0.0, 0.0, 0.0])):
        self.l      = l
        self.m      = m
        self.eta    = eta
        self.center = center
        self.normf  = ((2/pi)**(3.0/4.0))*((2.0**l)*eta**((2.0*l+3.0)/4.0)) #denominator is always 1 for s and p orbitals
        self.set_type('GTO')
    #end def __init__

    # Return value
    def v(self, pos):
        rij   = array(norm(pos-self.center, axis=1))
        vec   = pos-self.center
        l     = self.l
        m     = self.m
        eta   = self.eta
        normf = self.normf
        return normf*(vec[m+1]**l)*exp(-eta*r**2)
    #end def v
    
    #Return derivative
    def d(self, pos):
        rij   = array(norm(pos-self.center, axis=1))
        normf = self.normf
        n     = self.n
        eta   = self.eta
        
        return -2*v(r)
    #end def d

    #Return Laplacian
    def lap(self, pos):
        rij   = array(norm(pos-self.center, axis=1))
        normf = self.normf
        n     = self.n
        eta   = self.eta
        
        return (n*(n-1)/r**2 - 2*n*eta/r+eta**2)*v(r) 
    #end def lap
#end class STO

if __name__  == '__main__':
    ss = STO(alpha=1.)
    ss.plot_1D()
    pdb.set_trace()                                       
