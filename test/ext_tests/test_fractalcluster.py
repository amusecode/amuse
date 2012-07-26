from amuse.test import amusetest

from amuse.ic.fractalcluster import new_fractal_cluster_model

from amuse.units import nbody_system
from amuse.units import units

class TestEvrardModel(amusetest.TestCase):
    def test1(self):
        print "First test: making a fractal cluster."
        target_number_of_particles = 1000
        parts = new_fractal_cluster_model(N=target_number_of_particles)
        self.assertEquals(len(parts), 1000)

    def test2(self):
        print "test 2: test energy."
        target_number_of_particles = 1000
        parts = new_fractal_cluster_model(N=target_number_of_particles)
        ek=parts.kinetic_energy()
        ep=parts.potential_energy(G=nbody_system.G)
        self.assertAlmostEqual(ek/abs(ep),0.5,12)

    def test3(self):
        print "test 2: test energy physical units."
        target_number_of_particles = 1000
        convert_nbody = nbody_system.nbody_to_si(1000 | units.MSun, 1 | units.parsec) 
        parts = new_fractal_cluster_model(N=target_number_of_particles,convert_nbody=convert_nbody)
        ek=parts.kinetic_energy()
        ep=parts.potential_energy()
        self.assertAlmostEqual(ek/abs(ep),0.5,12)




