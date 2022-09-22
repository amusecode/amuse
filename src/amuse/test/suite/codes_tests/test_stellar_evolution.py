from amuse.test.amusetest import TestWithMPI

from amuse.community.sse.interface import SSE
from amuse.community.mesa.interface import MESA
from amuse.community.evtwin.interface import EVtwin
from amuse.community.seba.interface import SeBa
from amuse.community.bse.interface import BSE

import numpy
import time

from amuse.units import units
from amuse import datamodel


class _TestStellarEvolutionCodes(TestWithMPI):
    
    def code_factory(self):
        self.skip("abstract test")
    
    def test1(self):
        instance=self.code_factory()

    def test2(self):
        for m in ([0.2,1.,5.,25.]|units.MSun):
          instance=self.code_factory()()
          p=datamodel.Particle(mass=m)
          p2=instance.particles.add_particle(p)
          self.assertAlmostEqual(p2.mass,p.mass)
          self.assertTrue( hasattr(p2, "radius") )
          self.assertTrue( hasattr(p2, "luminosity") )
          self.assertTrue( hasattr(p2, "age") )
          self.assertTrue( hasattr(p2, "stellar_type") )

    def test3(self):
        instance=self.code_factory()()
        p=datamodel.Particle(mass=1. | units.MSun)
        p2=instance.particles.add_particle(p)
        self.assertAlmostEqual(p2.mass,p.mass)
        self.assertTrue( hasattr(p2, "radius") )
        self.assertTrue( hasattr(p2, "luminosity") )
        self.assertTrue( hasattr(p2, "age") )
        self.assertTrue( hasattr(p2, "stellar_type") )
        self.assertEqual( str(p2.stellar_type),"Main Sequence star")

    def test4(self):
        for m in ([0.2,1.,5.]|units.MSun):
          instance=self.code_factory()()
          p=datamodel.Particle(mass=m)
          p2=instance.particles.add_particle(p)
          self.assertAlmostEqual(p2.age, 0. | units.Myr)
          instance.evolve_model(1.| units.Myr)
          self.assertGreaterEqual(p2.age,  (1. | units.Myr) )
          self.assertGreaterEqual(instance.model_time,  (1. | units.Myr) )
          

    def test5(self):
        instance=self.code_factory()()
        p1=datamodel.Particle(mass=1. | units.MSun)
        p2=datamodel.Particle(mass=1. | units.MSun)
        p1=instance.particles.add_particle(p1)
        instance.evolve_model( 1.|units.Gyr)
        p2=instance.particles.add_particle(p2)
        instance.evolve_model( 2.|units.Gyr)
        self.assertGreaterEqual(p1.age, (1. | units.Gyr) )
        self.assertGreaterEqual(p2.age,(1. | units.Gyr) )
        self.assertLess(p2.age, (2. | units.Gyr) )

    def test6(self):
        instance1=self.code_factory()()
        p1=datamodel.Particle(mass=1. | units.MSun)
        p2=datamodel.Particle(mass=1. | units.MSun)
        p3=datamodel.Particle(mass=1. | units.MSun)
        p1=instance1.particles.add_particle(p1)
        instance1.evolve_model( 10.|units.Myr)
        p2=instance1.particles.add_particle(p2)
        instance1.evolve_model( 20.|units.Myr)
        
        instance2=self.code_factory()()
        p3=instance2.particles.add_particle(p3)
        instance2.evolve_model(10.|units.Myr)
        
        self.assertAlmostEqual(p2.mass,p3.mass)
        self.assertAlmostEqual(p2.radius,p3.radius)
        self.assertAlmostEqual(p2.luminosity,p3.luminosity)

    def test7(self):
        instance1=self.code_factory()()
        instance2=self.code_factory()()
        p1=datamodel.Particle(mass=1. | units.MSun)
        p1=instance1.particles.add_particle(p1)
        p1.evolve_one_step()
        p2=datamodel.Particle(mass=1. | units.MSun)
        p2=instance2.particles.add_particle(p2)
        instance2.evolve_model( p1.age)
        self.assertAlmostEqual(p1.mass,p2.mass)
        self.assertAlmostEqual(p1.radius,p2.radius)
        self.assertAlmostEqual(p1.luminosity,p2.luminosity)

    def test8(self):
        instance1=self.code_factory()()
        instance2=self.code_factory()()
        p1=datamodel.Particle(mass=1. | units.MSun)
        p1=instance1.particles.add_particle(p1)
        p1.evolve_for(10. | units.Myr)
        p2=datamodel.Particle(mass=1. | units.MSun)
        p2=instance2.particles.add_particle(p2)
        instance2.evolve_model( 10.|units.Myr)
        self.assertAlmostEqual(p1.mass,p2.mass)
        self.assertAlmostEqual(p1.radius,p2.radius)
        self.assertAlmostEqual(p1.luminosity,p2.luminosity)

        


"""
the following two tests will not work - this is to be fixed.

    def test6(self):
        instance=self.code_factory()()
        p1=datamodel.Particle(mass=1. | units.MSun)
        p2=datamodel.Particle(mass=1. | units.MSun,age=1.| units.Gyr)
        p1=instance.particles.add_particle(p1)
        p2=instance.particles.add_particle(p2)
        self.assertAlmostEqual(p1.age, (0. | units.Gyr) )
        self.assertAlmostEqual(p2.age, (1. | units.Gyr) )
        
    def test7(self):
        p1=datamodel.Particle(mass=8. | units.MSun)
        p2=datamodel.Particle(mass=8. | units.MSun,age=42. | units.Myr)

        instance1=self.code_factory()()
        p1=instance1.particles.add_particle(p1)
        instance2=self.code_factory()()
        p2=instance2.particles.add_particle(p2)
        
        instance1.evolve_model(44. | units.Myr)
        instance1.evolve_model(2. | units.Myr)

        self.assertalmostEqual( p1.mass,p2.mass)
"""

class TestSSECode(_TestStellarEvolutionCodes):
    def code_factory(self):
        return SSE

#class xTestBSECode(_TestStellarEvolutionCodes):
#    def code_factory(self):
#        return BSE

class TestMESACode(_TestStellarEvolutionCodes):
    def code_factory(self):
        try:
            MESA()
        except Exception as message:
            self.skip("Tried to instantiate a new object of the optional code with type '{0}', but this code is not available".format(MESA))
        return MESA

class TestSeBaCode(_TestStellarEvolutionCodes):
    def code_factory(self):
        return SeBa

class TestEVtwinCode(_TestStellarEvolutionCodes):
    def code_factory(self):
        return EVtwin
    

