from amuse.test import amusetest

from amuse.ext.kingmodel import MakeKingModel
from amuse.legacy.hermite0.interface import Hermite
from amuse.support.units import nbody_system
from amuse.support.units import units

class TestKingModel(amusetest.TestCase):
    def test1(self):
        print "First test: making a King model."
        number_of_particles = 10
        instance = MakeKingModel(number_of_particles, None, W0=6)
        stars = instance.result
        self.assertEquals(len(stars), number_of_particles)
        print stars
    
    def test2(self):
        print "Testing kinetic and potential energy of a King model realisation."
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        number_of_particles = 500
        instance = MakeKingModel(number_of_particles, convert_nbody, W0=6)
        stars = instance.result
        self.assertEquals(len(stars), number_of_particles)
        instance.move_particles_to_center_of_mass(stars)
        stars.radius = 0.0 | units.RSun
        gravity = Hermite(convert_nbody)
        gravity.particles.add_particles(stars)
        self.assertAlmostEqual(gravity.particles[0].mass, (1.0|units.MSun)/number_of_particles, 3, in_units=units.MSun)
        print convert_nbody.to_nbody(gravity.kinetic_energy), convert_nbody.to_nbody(gravity.potential_energy)
        self.assertAlmostEqual(convert_nbody.to_nbody(gravity.kinetic_energy), 0.25 | nbody_system.energy, 1)
        self.assertAlmostEqual(convert_nbody.to_nbody(gravity.potential_energy), -0.5 | nbody_system.energy, 1)
        gravity.stop()
    
    def xtest3(self):
        # Turned off by default for speed.
        print "King models with varying King dimensionless depth W0."
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        number_of_particles = 10
        for w_0 in range(1,17):
            stars = MakeKingModel(number_of_particles, convert_nbody, W0=w_0).result
            self.assertEquals(len(stars), number_of_particles)
    
    def test4(self):
        print "Testing maximal/minimal value of King dimensionless depth W0."
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        number_of_particles = 10
        try:
            stars = MakeKingModel(number_of_particles, convert_nbody, W0=17).result
            self.fail("Should not be able to create models with King dimensionless depth W0 >= 16.")
        except Exception as ex:
            self.assertEquals("makeking: must specify w0 < 16", str(ex))
        try:
            stars = MakeKingModel(number_of_particles, convert_nbody, W0=0).result
            self.fail("Should not be able to create models with King dimensionless depth W0 = 0.")
        except Exception as ex:
            self.assertEquals("float division", str(ex))
