#!/usr/bin/env pythonw
from itertools import combinations
from numpy import max, min, array, place, array2string, zeros, mean
from numpy.random import rand, randn, seed, uniform, normal
from numpy.linalg import norm
from particles import Particle, Structure
import pdb
from qmc_basic import obj

class Coulomb(obj):
    def __init__(self,s):
        if isinstance(s, Structure):
            self.npos = s.pos
            self.chg  = s.chg
        self.e    = 0.
        self.en   = 0.
        self.nn   = 0.
        self.ee   = 0.
        #self.box  = box_bc.box_nopbc
        self.mat  = None
        self._calc_nn()
    #end def __init__

    def _calc_nn(self):
        nn = 0.0
        natoms = len(self.npos)
        for i,j in combinations(range(natoms), 2):
            d = abs(norm(self.npos[i]-self.npos[j]))
            nn += self.chg[i]*self.chg[j]/d**2
        #end for
        self.e +=nn
        self.nn = nn
    #end def get_nn

    def _calc_en(self, epos=None):
        epos = array(epos)
        en = 0.0
        nelect = len(epos)
        if epos is not None:
            for j in range(nelect):
                for i in self.npos:
                    d = abs(norm(epos[j]-self.npos[j]))
                    en += -1.0*self.chg[j]/d**2
                #end for
                self.mat[j,j] += en
            #end for
        #end if
        self.e +=en
        self.en = en
    #end def

    def _calc_ee(self, epos=None):
        epos = array(epos)
        ee = 0.0
        nelect = len(epos)
        if nelect > 1:
            if epos is not None:
                for j in range(nelect):
                    for i in range(nelect):
                        d = abs(norm(epos[j]-epos[j]))
                        ee += 1.0/d**2
                        self.mat[i,j] += ee
                    #end for
                #end for
            #end if
        self.e += ee
        self.ee = ee
    #end def

    def v(self, epos=None):
        nelect = len(epos)
        self.mat = zeros((nelect, nelect))
        self.e = 0
        self._calc_ee(epos)
        self._calc_en(epos)
        self._calc_nn()
    #end def v

    def error(self, str):
        print str + '\n'
        exit()
    #end def error
    
#end class Coulomb

if __name__ == '__main__':
    #pp = Particle(chg=1.0,  pos=[0,0,0], w=1.0)
    #pe = Particle(chg=-1.0, pos=[0.5,0,0], w=1.0)
    #pp = [pp, pe]
    #cc = Coulomb()
    #cc.v(pp)
    ss = Structure(name='h2+1')
    cc = Coulomb(ss)
    energy = []
    scale = 0.1
    for i in range(100):
        cc.v(epos=[[normal(scale=scale),normal(scale=scale),normal(scale=scale)]])
        energy.append(cc.e)
    #end for
    print mean(energy)
    pdb.set_trace()
