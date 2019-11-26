from amuse.test import amusetest
from amuse.units import units, nbody_system, constants

from amuse.ic.fractalcluster import new_fractal_cluster_model


class TestFractalCluster(amusetest.TestCase):
    
    def test1(self):
        print("First test: making a fractal cluster.")
        target_number_of_particles = 100
        parts = new_fractal_cluster_model(N=target_number_of_particles)
        self.assertEqual(len(parts), 100)

    def test2(self):
        print("test 2: test energy.")
        target_number_of_particles = 100
        parts = new_fractal_cluster_model(N=target_number_of_particles)
        ek=parts.kinetic_energy()
        ep=parts.potential_energy(G=nbody_system.G)
        self.assertAlmostEqual(ek/abs(ep),0.5,12)
        self.assertAlmostRelativeEqual(ep, -0.5 | nbody_system.energy, 12)
    
    def test3(self):
        print("test 3: test energy physical units.")
        target_number_of_particles = 100
        convert_nbody = nbody_system.nbody_to_si(1000 | units.MSun, 1 | units.parsec) 
        parts = new_fractal_cluster_model(N=target_number_of_particles,convert_nbody=convert_nbody)
        ek=parts.kinetic_energy()
        ep=parts.potential_energy()
        self.assertAlmostEqual(ek/abs(ep),0.5,12)
        self.assertAlmostRelativeEqual(ep, -0.5 * constants.G * (1000|units.MSun)**2 / (1.0|units.parsec), 12)
    
    def test4(self):
        print("Test with masses")
        target_number_of_particles = 100
        masses = (list(range(1,11)) | units.MSun) * 1.0
        convert_nbody = nbody_system.nbody_to_si(1000 | units.MSun, 1 | units.parsec) 
        particles = new_fractal_cluster_model(masses=masses, convert_nbody=convert_nbody, do_scale=True)
        
        ek = particles.kinetic_energy()
        ep = particles.potential_energy()
        self.assertEqual(len(particles), 10)
        self.assertAlmostEqual(particles.total_mass(), 1000 | units.MSun) # Note: total_mass == converter's mass unit!
        self.assertAlmostEqual(masses.sum(), 55 | units.MSun) # Note: total_mass != masses.sum()
        self.assertAlmostEqual(particles.center_of_mass(), [0, 0, 0] | units.parsec)
        self.assertAlmostEqual(particles.center_of_mass_velocity(), [0, 0, 0] | units.km / units.s)
        self.assertAlmostEqual(ek/ep, -0.5, 12)
        self.assertAlmostRelativeEqual(ek, (0.25 * constants.G * (1000|units.MSun)**2 / (1.0|units.parsec)).as_quantity_in(ek.unit), 12)
    
    def test5(self):
        print("Test with masses, with correct mass unit in converter")
        target_number_of_particles = 100
        masses = (list(range(1,11)) | units.MSun) * 1.0
        convert_nbody = nbody_system.nbody_to_si(masses.sum(), 1 | units.parsec) 
        particles = new_fractal_cluster_model(masses=masses, convert_nbody=convert_nbody, do_scale=True)
        
        ek = particles.kinetic_energy()
        ep = particles.potential_energy()
        self.assertEqual(len(particles), 10)
        self.assertAlmostEqual(particles.total_mass(), 55 | units.MSun) # Note: total_mass == converter's mass unit!
        self.assertAlmostEqual(particles.center_of_mass(), [0, 0, 0] | units.parsec)
        self.assertAlmostEqual(particles.center_of_mass_velocity(), [0, 0, 0] | units.km / units.s)
        self.assertAlmostEqual(ek/ep, -0.5, 12)
        self.assertAlmostRelativeEqual(ek, (0.25 * constants.G * (55|units.MSun)**2 / (1.0|units.parsec)).as_quantity_in(ek.unit), 12)
    
    def test6(self):
        print("Test fractal dimension.")
        number_of_particles = 1000
        for target_fractal_dimension in [1.6, 2.0, 2.5, 3.0]:
            particles = new_fractal_cluster_model(
                N=number_of_particles, 
                fractal_dimension=target_fractal_dimension, 
                do_scale=False, random_seed=1234321)
            self.assertAlmostRelativeEquals(particles.box_counting_dimension(), 
                target_fractal_dimension, 1)
            self.assertAlmostRelativeEquals(particles.correlation_dimension(), 
                target_fractal_dimension, 1)
    

    
