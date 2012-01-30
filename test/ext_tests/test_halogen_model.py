from amuse.test.amusetest import TestWithMPI

from amuse.support.exceptions import AmuseException
from amuse.ext.halogen_model import new_halogen_model
from amuse.units import nbody_system, units

class NewHalogenModelTests(TestWithMPI):
    
    def test1(self):
        number_of_particles = 100
        particles = new_halogen_model(number_of_particles, alpha = 2.0, beta = 5.0, gamma = 0.0, random_seed = 1.0)
        
        self.assertEquals(len(particles), number_of_particles)
        self.assertAlmostEquals(particles.total_mass(), 1.0 | nbody_system.mass)
        self.assertAlmostEquals(particles.kinetic_energy(), 
            0.17345836639 | nbody_system.energy) # for number_of_particles = 100
        
        self.assertRaises(AmuseException, new_halogen_model, number_of_particles, expected_message = 
            "Error when calling 'commit_parameters' of a 'Halogen', errorcode is -2, error is "
            "'Missing or bad parameter for halo (see amuse/community/halogen/src/doc for details on required parameters).'")
    
    def test2(self):
        number_of_particles = 1000
        black_hole_mass = 1.0e6 | units.MSun
        stellar_mass = number_of_particles | units.MSun
        scale_radius = 0.1 | units.parsec
        
        
        converter = nbody_system.nbody_to_si(black_hole_mass, scale_radius)
        
        particles = new_halogen_model(number_of_particles, convert_nbody=converter, 
            alpha=1.0, beta=3.0, gamma=1.0, # NFW
            random_seed=1.0, redirection='none',
            black_hole_mass=black_hole_mass, total_mass=stellar_mass,
            cutoff_radius=10.0*scale_radius, scale_radius=scale_radius)
        
        self.assertEquals(len(particles), number_of_particles + 1)
        self.assertAlmostEquals(particles[-1].mass, black_hole_mass)
        self.assertAlmostEquals(particles.total_mass(), black_hole_mass + stellar_mass)
        self.assertAlmostRelativeEquals(particles.kinetic_energy(), 
            2.27538127277e+43 | units.J, 10)
    

