import numpy

from amuse.test import amusetest
from amuse.support.exceptions import AmuseException
from amuse.units import nbody_system
from amuse.units import units
from amuse.ic.kingmodel import new_king_model

class TestKingModel(amusetest.TestCase):
    def test1(self):
        print "First test: making a King model."
        number_of_particles = 10
        particles = new_king_model(number_of_particles, 6.0)
        self.assertAlmostEqual(particles.mass.sum(), 1 | nbody_system.mass)
        self.assertAlmostEqual(particles[0].mass , 0.1 | nbody_system.mass)
        print particles
        self.assertFalse(particles[0].mass.is_vector())
    
    def test2(self):
        print "Testing kinetic and potential energy of a King model realisation."
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        number_of_particles = 500
        particles = new_king_model(number_of_particles, 6.0, convert_nbody, do_scale = True)
        self.assertEquals(len(particles), number_of_particles)
        self.assertAlmostEqual(particles[0].mass, (1.0|units.MSun)/number_of_particles, 3, in_units=units.MSun)
        self.assertAlmostEqual(convert_nbody.to_nbody(particles.kinetic_energy()), 0.25 | nbody_system.energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(particles.potential_energy()), -0.5 | nbody_system.energy)
    
    def slowtest3(self):
        print "King models with varying King dimensionless depth W0."
        number_of_particles = 10
        for w_0 in [1.0, 6.0, 11.0, 16.0]:
            particles = new_king_model(number_of_particles, W0 = w_0)
            self.assertEquals(len(particles), number_of_particles)
    
    def test4(self):
        print "Testing maximal/minimal value of King dimensionless depth W0."
        number_of_particles = 10
        self.assertRaises(AmuseException, new_king_model, number_of_particles, W0 = 16.5, 
            expected_message = "makeking: must specify w0 < 16")
        
        self.assertRaises(ZeroDivisionError, new_king_model, number_of_particles, W0 = 0.0)
    
    def test5(self):
        print "Testing a specific King model realisation."
        numpy.random.seed(345672)
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        particles = new_king_model(500, 6.0, convert_nbody)
        self.assertEquals(len(particles), 500)
        self.assertAlmostEqual(particles.total_mass(), 1.0 | units.MSun)
        self.assertAlmostEqual(particles.mass, 1.0 / 500 | units.MSun)
        self.assertAlmostEqual(particles.center_of_mass(),          [0,0,0] | units.AU)
        self.assertAlmostEqual(particles.center_of_mass_velocity(), [0,0,0] | units.km / units.s)
        self.assertAlmostEqual(particles[:3].position, [[-0.23147381,-0.19421449,-0.01165137],
            [-0.09283025,-0.06444658,-0.07922396], [-0.44189946,0.23786357,0.39115629]] | units.AU)
    
