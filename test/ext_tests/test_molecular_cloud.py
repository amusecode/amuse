import sys
import os
import numpy.random

from amuse.test import amusetest
from amuse.units import units, nbody_system
from amuse.ext.molecular_cloud import molecular_cloud, ism_cube

numpy.random.seed(1234567)

class MolecularCloudTests(amusetest.TestCase):
    def test1(self):
        mc=molecular_cloud(targetN=1000).result
        self.assertEqual(len(mc),1000)
        ek=mc.kinetic_energy()
        ep=mc.potential_energy(G=nbody_system.G)
        eth=mc.thermal_energy()
        self.assertAlmostRelativeEqual(eth/ep, -0.01,2)
        self.assertAlmostRelativeEqual(ek/ep, -1.,2)

    def test2(self):
        mc=molecular_cloud(targetN=1000,ethep_ratio=0.05).result
        self.assertEqual(len(mc),1000)
        ek=mc.kinetic_energy()
        ep=mc.potential_energy(G=nbody_system.G)
        eth=mc.thermal_energy()
        self.assertAlmostRelativeEqual(eth/ep, -0.05,2)
        self.assertAlmostRelativeEqual(ek/ep, -1.,2)

    def test3(self):
        mc=molecular_cloud(targetN=1000,ekep_ratio=2.).result
        self.assertEqual(len(mc),1000)
        ek=mc.kinetic_energy()
        ep=mc.potential_energy(G=nbody_system.G)
        eth=mc.thermal_energy()
        self.assertAlmostRelativeEqual(eth/ep, -0.01,2)
        self.assertAlmostRelativeEqual(ek/ep, -2.,2)
