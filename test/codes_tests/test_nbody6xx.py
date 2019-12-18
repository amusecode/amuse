from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.nbody6xx.interface import Nbody6xxInterface
from amuse.community.nbody6xx.interface import Nbody6xx

import math

class Nbody6xxInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = Nbody6xxInterface()
        instance.initialize_code()
        res1 = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        res2 = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        self.assertEqual(1, res1['index_of_the_particle'])
        self.assertEqual(2, res2['index_of_the_particle'])
    
        retrieved_state1 = instance.get_state(1)
        retrieved_state2 = instance.get_state(2)
    
        self.assertEqual(11.0,  retrieved_state1['mass'])
        self.assertEqual(21.0,  retrieved_state2['mass'])
        self.assertEqual(0.0,  retrieved_state1['x'])
        self.assertEqual(10.0,  retrieved_state2['x'])
        #self.assertEquals(2.0,  retrieved_state1['radius'])
        #self.assertEquals(5.0,  retrieved_state2['radius'])
    
        instance.cleanup_code()
        instance.stop()

    
    def test2(self):
        instance = Nbody6xxInterface()
        instance.initialize_code()
        
        instance.new_particle([10,20],[0,0],[0,0], [0,0], [0,0], [0,0], [0,0],[1,1])
        retrieved_state = instance.get_state(1)
        
        self.assertEqual(10.0,  retrieved_state['mass'])
        #self.assertEquals(1, retrieved_state['radius'])
    
        retrieved_state = instance.get_state([1,2])
        
        self.assertEqual(20.0,  retrieved_state['mass'][1])
        self.assertEqual(instance.get_number_of_particles()['number_of_particles'], 2)
        instance.cleanup_code() 
        instance.stop()
        
        
    def xtest3(self):
        instance = Nbody6xxInterface(redirection="none")#,debugger="ddd")
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        
        # Set up an equal-mass binary on a circular orbit:
        self.assertEqual([1, 0], list(instance.new_particle(0.5,  0.5, 0, 0,  0, 0.5, 0, 0.01).values()))
        self.assertEqual([2, 0], list(instance.new_particle(0.5,  -0.5, 0, 0,  0,-0.5, 0, 0.01).values()))
        self.assertEqual([3, 0], list(instance.new_particle(100.5,  -0.5, 0, 0,  0,-0.5, 0, 0.01).values()))
        self.assertEqual(0, instance.commit_particles())
        
        self.assertEqual(0, instance.evolve_model(math.pi))
        for result, expected in zip(instance.get_position(1).values(), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)
        for result, expected in zip(instance.get_position(2).values(), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)
        
        self.assertEqual(0, instance.evolve_model(2 * math.pi))
        for result, expected in zip(instance.get_position(1).values(), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)
        for result, expected in zip(instance.get_position(2).values(), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 3)
        
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

