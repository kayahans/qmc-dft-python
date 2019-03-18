#!/usr/bin/env pythonw

from numpy import array, abs, pi, sqrt, power, log, arctan
from util.basic import obj
from functionals import Xc
from copy import copy
import pdb


def get_xc(grid, D, **kwargs):
    rho = grid.getdens(D)
    
    xc = Xc(rho, name='ldat', pol=False)
    w = grid.points[:,3]

    fx   = xc.get_xaF()
    fc   = xc.get_caF()
    
    dfxa = xc.get_xaf()
    dfca = xc.get_caf()

    Vxc = np.einsum('g,g,gI,gJ->IJ',w,dfxa+dfca,grid.bfamps,grid.bfamps)

    Exc = np.dot(w, 2*fx+fc)

    o = obj()
    o.e = Exc
    o.v = Vxc

    return o
#end def

