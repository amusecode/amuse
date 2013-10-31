import sys
import os
import numpy.random

from amuse.test import amusetest
from amuse.units import units, nbody_system
from amuse.ext.boss_bodenheimer import bb79_cloud

numpy.random.seed(1234567)

class BossBodenheimerTests(amusetest.TestCase):
    def test1(self):
        mc=bb79_cloud(targetN=1000).result
        self.assertEqual(len(mc),1000)
        ek=mc.kinetic_energy()
        ep=mc.potential_energy(G=nbody_system.G)
        eth=mc.thermal_energy()
        self.assertAlmostRelativeEqual(eth/ep, -0.25,2)
        self.assertAlmostRelativeEqual(ek/ep, -0.2,1)

    def test2(self):
        convert=nbody_system.nbody_to_si(1. | units.MSun,3.2e16| units.cm)
        mc=bb79_cloud(targetN=1000,convert_nbody=convert).result
        self.assertEqual(len(mc),1000)
        ek=mc.kinetic_energy()
        ep=mc.potential_energy()
        eth=mc.thermal_energy()
        self.assertAlmostRelativeEqual(eth/ep, -0.25,2)
        self.assertAlmostRelativeEqual(ek/ep, -0.2,2)

