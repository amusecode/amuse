from legacy_support import TestWithMPI


from amuse.legacy.seba.interface import SebaInterface, Seba
from amuse.support.units import units
from amuse.support.data.core import Particle

class TestMPIInterface(TestWithMPI):

    def test1(self):
        instance = SebaInterface()
        endtime, mass, radius, error = instance.evolve_star(5, 130, 0.02)
        self.assertEquals(error, 0)
        self.assertTrue( endtime <= 130.0)
        self.assertAlmostRelativeEqual(mass, 0.9906, 4)
        
class TestOOInterface(TestWithMPI):

    def test1(self):
        instance = Seba()
        endtime, mass, radius = instance.evolve(5 | units.MSun, 130 | units.Myr, 0.02 | units.none)
        
        self.assertTrue( endtime <= 130 | units.Myr)
        self.assertAlmostRelativeEqual(mass, 0.9906 | units.MSun, 4)
        self.assertAlmostRelativeEqual(radius, 0.0079867 | units.RSun, 4)
        
    def test2(self):
        instance = Seba()
        p = Particle()
        p.mass = 5 | units.MSun
        p.metallicity = 0.02 | units.none
        
        p = instance.particles.add_particle(p)
        instance.evolve_model(130 | units.Myr)

        self.assertAlmostRelativeEqual(p.mass, 0.9906 | units.MSun, 4)
        
