#!/usr/bin/env pythonw

from numpy import array, abs, pi, sqrt, power, log, arctan
from util.basic import obj
from copy import copy
import pdb
import sys

sys.dont_write_bytecode = True

def zero_low_rho(rho, tol=1e-10):
    rho[rho<tol]=0.
    return rho
#end def

class Functional(obj):
    def __init__(self, rhoa, rhob):
        self.rhoa = rhoa
        self.rhob = rhob
        self.e   = 0.
        self.v   = 0.
    #end def
#end class

class XCFun(Functional):
    def __init__(self, name = 'ldax', pol= False):
        self.name = name
        self.pol  = pol
    #end def

    def get_xc(self, density):
        o = obj()
        name = self.name
        if name == 'ldax':
            o.x = self.calc_x(density)
            o.c = self.empty(density)
        elif name == 'lda':
            o.x = self.calc_xs(density)
            o.c = self.calc_vwn5(density)
        #end if

        return o
    #end def
    
            
    def get_xaf(self):
        return self.xa.f
    #end def

    def get_xaF(self):
        return self.xa.F
    #end def

    def get_xbf(self):
        return self.xb.f
    #end def

    def get_xbF(self):
        return self.xb.F
    #end def
    
    def get_caf(self):
        return self.ca.f
    #end def

    def get_caf(self):
        return self.ca.F
    #end def

    def empty(self,rho):
        o = obj()
        o.F = 0
        o.f = rho*0.
        return o
    #end def
        
    def calc_xs(self, rho):
        '''
        2/3 alpha exchange
        '''
        alpha = 2./3
        fac = -2.25*alpha*power(.75/pi, 1./3)
        rho_3 = power(rho, 1./3)
        F = fac*rho*rho_3
        f = 4./3*fac*rho_3
        o = obj()
        o.F = F
        o.f = f
        return o
    #end def

    def calc_x(self, rho):
        '''
        3/4 slater exchange
        '''
        fac = power(3./pi, 1./3)
        f   = lambda fac, rho : -fac*power(rho,1./3)
        F   = lambda fac, rho : -sum(3./4*fac*power(rho, 4./3))
        F=F(fac,rho)
        f=f(fac,rho)
    
        o = obj()
        o.F = F
        o.f = f
        return o
    #end def

    def calc_vwn5(self, rhoa, rhob):
        '''
        https://www.molpro.net/info/2010.1/doc/manual/node772.html
        '''
        rhoa = zero_low_rho(rhoa)
        rhob = zero_low_rho(rhob)
        rho  = rhoa+rhob
        zeta = (rhoa-rhob)/rho

        p1 = [0.0310907,-0.10498,3.72744,12.9352]
        p2 = [0.01554535,-0.32500,7.06042,18.0578]
        p3 = [-1./6*pi**-2, -0.0047584, 1.13107, 13.0045]
        
        x = lambda rho: pow(3./4./pi/rho,1/6.)
        y = lambda zeta: 9./8*(1+zeta)**(4./3) + 9./8*(1-zeta)**(4./3)-9./4

        x = x(rho)
        y = y(zeta)
        
        X = lambda i,c,d: i**2 + c*i +d
        Q = lambda c,d: sqrt(4*d-c**2)
        q = lambda A, p, c, d: A*(log(x**2/X(x,c,d)) +
                                      2*c*arctan(Q(c,d)/(2*x+c))/Q(c,d) -
                                      c*p*(log((x-p)**2/X(x,c,d))+ 2*(c+2*p)*arctan(Q(c,d)/(2*x+c))/Q(c,d))/X(p,c,d))

        Lam = q(p1)
        lam = q(p2)
        alpha = q(p3)

        h = lambda lam, Lam, alpha: 4./9*(lam-Lam)/((2**(1./3)-1)*alpha) - 1
        e = lambda Lam, alpha, y, h, zeta: Lam + alpha*y*(1+h*zeta**4)
        
        
        h = h(lam, Lam, alpha)
        e = e(Lam, alpha, y, h, zeta)
    #end def
#end class

    

