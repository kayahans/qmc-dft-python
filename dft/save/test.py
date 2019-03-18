#!/usr/bin/env python

from utils.orbitals import STO
from hamiltonians import DFT

class test_dft():
    def test_h_lda_sto(self):
        from utils.geom import h
        basisset = STO()
        solver   = DFT(h, basisset, 'xs')
        ens      = solver.converge()
        self.assertAlmostEqual(solver.energy,-0.5)
      #end def
#end class


