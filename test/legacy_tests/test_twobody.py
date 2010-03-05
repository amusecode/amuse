

from amuse.support.data import core
from amuse.support.units import units
from amuse.support.units import nbody_system

import numpy

from amuse.legacy.twobody import twobody

from legacy_support import TestWithMPI


class TwoBodyCodeTests(TestWithMPI):
    
    def test_stumpff(self):
        self.assertAlmostEqual(twobody.stumpff_C(0),twobody.stumpff_C(0.0001),5)
        self.assertAlmostEqual(twobody.stumpff_C(0),twobody.stumpff_C(-0.0001),5)
        self.assertAlmostEqual(twobody.stumpff_S(0),twobody.stumpff_S(0.0001),5)
        self.assertAlmostEqual(twobody.stumpff_S(0),twobody.stumpff_S(-0.0001),5)

class TwoBodyInterfaceTests(TestWithMPI):
        
    def test1(self):
        
        instance = twobody.TwoBodyInterface()
        
        res1 = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        res2 = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        self.assertEquals(0, res1['index_of_the_particle'])
        self.assertEquals(1, res2['index_of_the_particle'])

        retrieved_state1 = instance.get_state(0)
        retrieved_state2 = instance.get_state(1)

        self.assertEquals(11.0,  retrieved_state1['mass'])
        self.assertEquals(21.0,  retrieved_state2['mass'])
        self.assertEquals(0.0,  retrieved_state1['x'])
        self.assertEquals(10.0,  retrieved_state2['x'])

        del instance


        
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(5.9742e24 | units.kg, 1e6| units.m)
        instance = twobody.TwoBody(convert_nbody)
        
        p = core.Particle()
        p.mass = 5.9742e24 | units.kg
        p.radius = 6.371e6 | units.m
        p.position = [0.,7.e6,-1.2124e7] | units.m
        p.velocity = [0.,2.6679e3,4.6210e3] | units.m/units.s
        
        instance.particles.add_particle(p)
        
        instance.evolve_model(3600.0 | units.s)
        
        position = instance.particles[0].position
        velocity = instance.particles[0].velocity
        
        self.assertAlmostEqual(position.x.value_in(units.m),0.,7)
        self.assertAlmostEqual(position.y.value_in(units.m)/(-3.30647600568e6),1.,7)
        self.assertAlmostEqual(position.z.value_in(units.m)/7.40831575351e6,1.,7)
        self.assertAlmostEqual(velocity.x.value_in(units.m / units.s),0.,7)
        self.assertAlmostEqual(velocity.y.value_in(units.m / units.s)/(-8.29821376206e3),1.,7)
        self.assertAlmostEqual(velocity.z.value_in(units.m / units.s)/(-0.972888312209e3),1.,7)

    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(5.9742e24 | units.kg, 1e6| units.m)
        instance = twobody.TwoBody(convert_nbody)
        
        p = core.Particle()
        p.mass = 5.9742e24 | units.kg
        p.radius = 7.1e6 | units.m
        p.position = [0.,7.e6,-1.2124e7] | units.m
        p.velocity = [0.,2.6679e3,4.6210e3] | units.m/units.s
        
        instance.particles.add_particle(p)
        
        instance.evolve_model(3600.0 | units.s)
        
        dt = convert_nbody.to_si(instance.model_time)
        self.assertAlmostEqual(dt.value_in(units.s)/2583.44780926,1.,7)
        
        position = instance.particles[0].position
        self.assertAlmostEqual(((position.x**2+position.y**2+position.z**2)/(7.1e6)**2).value_in(units.m**2),1.,7)
