#!/usr/bin/env pythonw

from numpy import array, abs, pi, sqrt,linspace, zeros, exp, reshape, isclose
from util.basic import obj
from scipy.special import factorial2, factorial
from grid import Grid
from collections import OrderedDict

import sys, pdb
sys.dont_write_bytecode = True

class Orbital(obj):
    def __init__(self, type):
        self.type = None
    #end def __init__

    def __repr__(self):
        attr = self.__dict__
        attr = OrderedDict(sorted(attr.items()))
        c = ''
        for k,v in attr.items():
            c += '\n' + str(k) + '\t' + str(getattr(self,k))
        #end for
        return c[1:] #remove heading whitespace
    #end def repr
    
    def error(self,str):
        print str + '\n'
        exit()
    #end def error

    def set_type(self, t):
        self.type = t
    #end def set_type

    def plot_1D(self):
        import matplotlib.pyplot as plt
        x = []
        for i in linspace(0,3,101):
            x.append([i,0.,0.])
        #end for
        x = array(x)
        v = self.value(x)
        rij = x[:,0]
        plt.plot(rij,v, label='pow= '+str(self.pow))
        plt.axhline(y=0, linestyle='--', color='k')
        plt.axhline(y=1, linestyle='--', color='k')
        plt.axvline(x=0, linestyle='--', color='k')
        plt.axvline(x=3, linestyle='--', color='k')
        plt.xlabel('Distance (r)')
        plt.ylabel('Amplitude')
        plt.legend()
        return plt
    #end def
    
#end class BasisSet

class pgto(Orbital):
    '''
    Primitive Gaussian Type Orbital in Cartesian Basis 
    https://journals.jps.jp/doi/pdf/10.1143/JPSJ.21.2313

    g(x,y,z) = A*(x^l)*(y^m)*(z^n)*exp{-a*(r-r_0)^2}
    '''    
    def __init__(self, origin=(0,0,0), pow=(0,0,0), exp=1.0):
        assert len(origin) == 3
        assert len(pow) == 3
        
        self.type   = 'Primitive GTO'
        self.origin = origin
        self.pow    = pow
        self.exp    = float(exp)
        self._normalize()
    #end def
    
    def _normalize(self):
        l,m,n = self.pow
        spow  = l+m+n

        self.norm = ((2**(2*spow+3./2)*self.exp**(spow+3./2)) /
                     (factorial2(2*l-1)*factorial2(2*m-1)*factorial2(2*n-1)*pi**(3./2)))**(1./2)
    #end def
    
    def _eval(self, func, grid):
        '''
        Evaluate any basisset function over the grid
        '''
        x0, y0, z0 = self.origin
        if isinstance(grid, Grid):
            xn, yn, zn = grid.get_gridv()
            shape      = grid.shape
        else:
            grid=array(grid)
            if grid.ndim == 1:
                xn = grid[0]
                yn = grid[1]
                zn = grid[2]
            elif grid.ndim == 2:
                xn = grid[:,0]
                yn = grid[:,1]
                zn = grid[:,2]
            else:
                self.error('Not Grid or list of points')
            #end if
        #end if
        
        dx, dy, dz = xn-x0, yn-y0, zn-z0
        r_sq =dx**2+dy**2+dz**2
    
        g = func(dx, dy, dz, r_sq)

        if isinstance(grid, Grid):
            g = reshape(g, shape)

        return g
    #end def

    def value(self, grid):

        # g(x,y,z)
        l, m, n = self.pow
        A       = self.norm
        a       = self.exp
        
        g = lambda x, y, z, r_sq : A*(x**l)*(y**m)*(z**n)*exp(-a*r_sq)
        # Evaluate function g over grid
        return self._eval(g, grid)
    #end def

    def S(self, orb):
        '''
        Overlap integral
        Page 2315
        https://journals.jps.jp/doi/pdf/10.1143/JPSJ.21.2313
        '''
        assert isinstance(orb, pgto)
        l1, m1, n1 = self.pow
        l2, m2, n2 = orb.pow
        a1 = self.exp
        a2 = orb.exp

        A     = self.origin
        B     = orb.origin
        gamma = a1+a2
        
        P     = (a1*A+a2*B)/(a1+a2)
        AB_sq = dot(A-B, A-B)
        PA    = sqrt(dot(P-A, P-A))
        PB    = sqrt(dot(P-B, P-B))
        
        ## fj(l,m,a,b)
        f = lambda l,m,a,b,j,i : comb(l,i)*comb(m,j-i)*a**i*b**(j-i) 

        # Sum over 0 <= i <= j for all possible j <=l1+l2
        fj = []
        for j in range(l1+l2+1):
            fj_temp = []
            for i in range(j+1):
                fj_temp.append(f(l1, l2, PA, PB))
            #end for
            fj.append(sum(fj_temp))
        #end for
        
        ## fj complete

        # start integral
        
        pre = (pi/gamma)**3./2*exp(-a1*a2*AB_sq/gamma)

        int_out = 0
        for q in [l1+l2, m1+m2, n1+n2]:
            int_in = 1
            for i in range(int(ceil(0.5*(q)))):
                int_in *= f[2*i]*factorial2(2*i-1)/(2*gamma)**i
            #end for
            int_out += int_in
        #end for
        
        result = pre*int_out
        
        return result        
    #end def
    
    def T(self, orb):
        '''
        Kinetic energy
        '''
        assert isinstance(orb, pgto)
        a2 = orb.exp
        l1, m1, n1 = self.pow
        l2, m2, n2 = orb.pow
        a1 = self.exp
        a2 = orb.exp

        result = a2*(2*(l2+m2+n2)+3)*self.S(orb) \
            -2*a2**2(self.S(orb.get_new_pow((l2+2, m2, n2))) + self.S(orb.get_new_pow((l2, m2+2, n2))) + self.S(orb.get_new_pow((l2, m2, n2+2)))) \
            - 1./2*(l2*(l2-1)*self.S(orb.get_new_pow((l2-2, m2, n2))) + m2*(m2-1)*self.S(orb.get_new_pow((l2, m2-2, n2))) + n2*(n2-1)*self.S(orb.get_new_pow((l2, m2, n2-2))))

        return result
        
    #end def
    
    def V(self, orb, atom):
        '''
        Nuclear Attraction

        Incomplete
        '''
        assert isinstance(orb, pgto)
        a2 = orb.exp
        l1, m1, n1 = self.pow
        l2, m2, n2 = orb.pow
        a1 = self.exp
        a2 = orb.exp
        
        A     = self.origin
        B     = orb.origin
        C      = atom.pos
        gamma = a1+a2
        
        P     = (a1*A+a2*B)/(a1+a2)
        AB_sq = dot(A-B, A-B)
        PA    = sqrt(dot(P-A, P-A))
        PB    = sqrt(dot(P-B, P-B))
        PC    = sqrt(dot(P-C, P-C))
        
        ## fj(l,m,a,b)
        f = lambda l,m,a,b,j,i : comb(l,i)*comb(m,j-i)*a**i*b**(j-i)

        pre = 1./(2*pi**2)*exp(-a1*a2*AB_sq/gamma)
        
        def A_term(i,r,u,l1,l2,PAx,PBx,CPx,gamma):

            return pow(-1,i)*f(i,l1,l2,PAx,PBx)*\
                pow(-1,u)*factorial(i)*pow(CPx,i-2*r-2*u)*\
                pow(0.25/gamma,r+u)/factorial(r)/factorial(u)/factorial(i-2*r-2*u)
        #end def

        def A_array(l1,l2,PA,PB,CP,g):

            Imax = l1+l2+1
            A = [0]*Imax
            for i in range(Imax):
                for r in range(int(floor(i/2)+1)):
                    for u in range(int(floor((i-2*r)/2)+1)):
                        I = i-2*r-u
                        A[I] = A[I] + A_term(i,r,u,l1,l2,PA,PB,CP,g)
            
            return A
        #end def

        def Fgamma(m,x):
            eps=1e-12
            x = max(x,eps)
            return 0.5*pow(x,-m-0.5)*gamm_inc(m+0.5,x)
        
        Ax = A_array(l1,l2,PA[0],PB[0],PC[0],gamma)
        Ay = A_array(m1,m2,PA[1],PB[1],PC[1],gamma)
        Az = A_array(n1,n2,PA[2],PB[2],PC[2],gamma)
        
        total = 0.
        for I in range(l1+l2+1):
            for J in range(m1+m2+1):
                for K in range(n1+n2+1):
                    total += Ax[I]*Ay[J]*Az[K]*Fgamma(I+J+K,rcp2*gamma)
                    
        result= pre*total
        return result
    #end def

    def get_new_pow(self, pow):
        a = pgto(origin=self.origin, pow=pow, exp=self.exp)
        return a
    #end def
    
#end class

class cgto(obj):
    '''
    Contracted Gaussian Type Orbital in Cartesian Basis 
    
    Not implemented yet
    
    '''
    
    def __init__(self, origin=(0,0,0), pow=(0,0,0), exps=[], coeffs=[]):
        assert len(origin) == 3
        assert len(pow) == 3

        self.origin = origin
        self.pow    = pow
        self.pgtos  = [] # Made of pgtos
        self.coeffs = []
        for expn, coeffn in zip(exps, coeffs):
            self.add_pgto(expn, coef)

        self.normalize()
    #end def
        
    def add_pgto(self, expn, coefn):
        self.pgtos.append(pgto(origin=self.origin, pow=self.pow, exp=expn))
        self.coeffs.append(coefn)
    #end def

    def normalize(self):
        pass
    #end def
        
class BasisSet(obj):
    def __init__(self,molecule, name='sto3g'):
        self.bfs = []
        self.shells = []
            
if __name__  == '__main__':
    print '1s orbital '
    s1 = pgto()
    assert(isclose(s1.value([0,0,0]), 0.7127054))
    print s1
    plt  = s1.plot_1D()
    print
    print 
    print 'px orbital '
    px = pgto(pow=(1,0,0))
    assert(isclose(px.value([0,0,0]), 0.0))
    print px
    px.plot_1D()

    plt.show()
