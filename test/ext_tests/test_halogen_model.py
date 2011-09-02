from amuse.test.amusetest import TestWithMPI

from amuse.support.exceptions import AmuseException
from amuse.ext.halogen_model import new_halogen_model
from amuse.units import nbody_system
class NewHalogenModelTests(TestWithMPI):
    
    def test1(self):
        number_of_particles = 100
        particles = new_halogen_model(number_of_particles, alpha = 2.0, beta = 5.0, gamma = 0.0, random_seed = 1.0)
        
        self.assertEquals(len(particles), number_of_particles)
        self.assertAlmostEquals(particles.total_mass(), 1.0 | nbody_system.mass)
        self.assertAlmostEquals(particles.kinetic_energy(), 
            0.164458665999 | nbody_system.energy) # for number_of_particles = 100
        
        self.assertRaises(AmuseException, new_halogen_model, number_of_particles, expected_message = 
            "Error when calling 'commit_parameters' of a 'Halogen', errorcode is -2, error is "
            "'Missing or bad parameter for halo (see amuse/community/halogen/src/doc for details on required parameters).'")
    

