from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.nbody6xx.interface import Nbody6xxInterface
from amuse.community.nbody6xx.interface import Nbody6xx

class Nbody6xxInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = Nbody6xxInterface(redirection="none")
        instance.initialize_code()
        res1 = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        res2 = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        self.assertEquals(1, res1['index_of_the_particle'])
        self.assertEquals(2, res2['index_of_the_particle'])
    
        retrieved_state1 = instance.get_state(1)
        retrieved_state2 = instance.get_state(2)
    
        self.assertEquals(11.0,  retrieved_state1['mass'])
        self.assertEquals(21.0,  retrieved_state2['mass'])
        self.assertEquals(0.0,  retrieved_state1['x'])
        self.assertEquals(10.0,  retrieved_state2['x'])
        #self.assertEquals(2.0,  retrieved_state1['radius'])
        #self.assertEquals(5.0,  retrieved_state2['radius'])
    
        instance.cleanup_code()
        instance.stop()

    
