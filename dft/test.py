#!/usr/bin/env python

from system    import System
from util.geom import h, he

hydrogen = System(h)
hydrogen.solve()
h_density = hydrogen.get_density()
h_density.plot_density()

helium = System(he)
helium.solve()
h_density = helium.get_density()
h_density.plot_density()
