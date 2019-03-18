#!/usr/bin/env pythonw

from numpy import array, abs, pi, sqrt, power, log, arctan
from util.basic import obj
from grid import Grid 

class DFT():
    def __init__(self, geo, basisset, xcname = 'ldat'):
        self.grid = Grid(geo)
        self.basisset = basisset


if __name__ == '__main__':
    from utils.geom import h
    from utils.orbitals import STO
    basisset = STO()
    solver = DFT(h, basisset, xcname = 'ldat')
