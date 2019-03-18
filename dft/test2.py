#!/usr/bin/env pythonw

from numpy import linspace, meshgrid, sqrt, ones, ravel, pi, power, dot, real, exp
from scipy.sparse import spdiags, eye, kron, csr_matrix,csc_matrix
from scipy.sparse.linalg import eigs,eigsh, cgs, lsqr
from scipy.linalg import eigh
from scipy.special import erf
import pdb

#Grid
g  = 40
g3 = g**3
p = linspace(-3, 3, g)
h = p[1]-p[0]
x, y, z = meshgrid(p,p,p)
xn=ravel(x)
yn=ravel(y)
zn=ravel(z)

#external potential
r = sqrt(xn**2 + yn**2 + zn**2)
vext = -1./r
Vext3 = spdiags(vext, 0, g3, g3)

#Kinetic energy
e = ones(g)
Lap = spdiags([e, -2*e, e], [-1,0,1], g, g)/h**2 # Kinetic in 1-D
I = eye(g)
Lap3 = kron(kron(Lap,I),I) + kron(kron(I,Lap),I) + kron(kron(I,I),Lap)

#Solve
E, psi = eigs(-0.5*Lap3+Vext3, k=1, which='SR')

print 'Hydrogen', E


#Helium
#external potential
r = sqrt(xn**2 + yn**2 + zn**2)
vext = -2./r

Eprev = 0.
ediff = 10**6
vtot = vext

comp_chg = True

if comp_chg:
    ncomp = exp(-r**2/2)
    ncomp = -2.*ncomp/sum(ncomp)/h**3
    vcomp = -2./r*erf(r/sqrt(2))
else:
    ncomp = 0.
    vcomp = 0.
#end if

while ediff > 10**-6:
    Vtot3 = spdiags(vtot, 0, g3, g3)

    E, psi = eigsh(-0.5*Lap3+Vtot3, k=1, which='SA') #Real eigenval and vectors
    #pdb.set_trace()
    psi = ravel(psi)
    psi = psi/h**(3./2)
    n = 2*psi**2
    
    #exchange
    vex = -(3./pi)**(1./3)*n**(1./3)
    
    #hartree
    #Lap3 Vh = -4*pi*n
    vh = cgs(Lap3, -4*pi*(n+ncomp))[0] - vcomp#Poisson solver
    
    vtot = vex+vh+vext
    
    T1   = csr_matrix.dot(-0.5*Lap3,psi)
    T2   = csc_matrix.dot(csc_matrix(psi), T1)
    T    = T2*2*h**3
    T    = T[0]

    Eext = sum(n*vext)*h**3
    Eh   = 0.5 * sum(n*vh)*h**3
    Ex   = sum(-(3./4)*(3./pi)**(1./3)*n**(4./3))*h**3

    Etot = T + Eext + Eh + Ex

    ediff = abs(Eprev-Etot)
    Eprev = Etot
    E     = E[0]
    #pdb.set_trace()
    print 'Helium', Etot, T, Eext, Eh, Ex, E, ediff
    
