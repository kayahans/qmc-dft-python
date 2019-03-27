#!/usr/bin/env python

from util.geom    import h, he
from system import System
import pdb
    
print 'LDAx Hydrogen no basis-set'
s = System(h)
s.solve()
print 'LDAx Helium no basis-set'
s = System(he)
s.solve()
pdb.set_trace()
    
