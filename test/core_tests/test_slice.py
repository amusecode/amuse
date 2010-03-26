from amuse.test import amusetest

from amuse.support.data import core

class TestSliceParticles(amusetest.TestCase):
    
    def test1(self):
        print "Test: slice a particle set."
        number_of_particles = 10
        original_set = core.Particles(number_of_particles)
        self.assertTrue(isinstance(original_set, core.Particles))
        self.assertEqual(len(original_set),number_of_particles)
        print "Defining all kind of slices of the particle set..."
        subset1 = original_set[:2]  # contains first two particles
        subset2 = original_set[2:]  # contains the rest of the particles
        odd = original_set[1::2]    # contains all particles with odd indices
        even = original_set[::2]    # contains all particles with even indices
        reverse = original_set[::-1]# contains all particles in reverse order
        all = original_set[:]       # contains all particles
        one = original_set[3]       # contains one particle (Particle)
        another = original_set[5:6] # contains one particle (ParticlesSubset)
        empty = original_set[2:2]   # contains no particle
        print "Checking their type (slicing returns a subset, indexing returns a particle)... ",
        self.assertTrue(isinstance(subset1, core.ParticlesSubset))
        self.assertTrue(isinstance(subset2, core.ParticlesSubset))
        self.assertTrue(isinstance(odd,     core.ParticlesSubset))
        self.assertTrue(isinstance(even,    core.ParticlesSubset))
        self.assertTrue(isinstance(reverse, core.ParticlesSubset))
        self.assertTrue(isinstance(all,     core.ParticlesSubset))
        self.assertTrue(isinstance(one,     core.Particle))
        self.assertTrue(isinstance(another, core.ParticlesSubset))
        self.assertTrue(isinstance(empty,   core.ParticlesSubset))
        print "ok!"
        print "Checking their length... ",
        self.assertEqual(len(subset1), 2)
        self.assertEqual(len(subset2), number_of_particles-2)
        self.assertEqual(len(odd),     int((number_of_particles/2.0)))
        self.assertEqual(len(even),    int((0.5+number_of_particles/2.0)))
        self.assertEqual(len(reverse), number_of_particles)
        self.assertEqual(len(all),     number_of_particles)
        self.assertEqual(len([one]),   1)
        self.assertEqual(len(another), 1)
        self.assertEqual(len(empty),   0)
        print "ok!"
    
