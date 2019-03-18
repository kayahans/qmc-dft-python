#!/usr/bin/env pythonw

from numpy import array, meshgrid, linspace, ravel, ones
from scipy.sparse import spdiags, eye, kron, csr_matrix,csc_matrix
from util.basic import obj
from util.geom import Molecule
from copy import copy


class Grid(obj):
    def __init__(self, system, shape, pad = 5.):
        assert isinstance(system, Molecule)
        pad = float(pad)
        self.bbox   = system.bbox(pad=pad)

        if isinstance(shape, int):
            shape = [shape]*3
        #end if

        self.shape   = array(shape)
    #end def

    def get_grid(self):
        min_corner  = self.bbox.min
        max_corner  = self.bbox.max
        shape = self.shape
        px = linspace(min_corner[0], max_corner[0], shape[0])
        py = linspace(min_corner[1], max_corner[1], shape[1])
        pz = linspace(min_corner[2], max_corner[2], shape[2])

        self.dv = (px[1]-px[0])*(py[1]-py[0])*(pz[1]-pz[0])
        self.dx = (px[1]-px[0], py[1]-py[0], pz[1]-pz[0])
        return meshgrid(px,py,pz)
    #end def

    def get_gridv(self):
        x, y, z = self.get_grid()
        
        xn = ravel(x)
        yn = ravel(y)
        zn = ravel(z)
        
        return xn, yn, zn
    #end def

    def get_lapn(self):
        gx, gy, gz  = self.shape
        self.get_grid()
        dx, dy, dz = self.dx
        e = ones(gx)
        Lap = spdiags([e, -2*e, e], [-1,0,1], gx, gx)/dx**2
        iy = eye(gy)
        iz = eye(gz)
        Lap3 =  kron(kron(Lap,iy),iz) + kron(kron(iz,Lap),iy) + kron(kron(iy,iz),Lap)

        return Lap3
    #end def
        
if __name__ == '__main__':
    from util.geom import h, h2
    g  = Grid(h2, 30)
    g2 = Grid(h, 30)
