import numpy

from amuse.test.amusetest import TestWithMPI
from amuse.community.higpus.interface import HiGPUsInterface, HiGPUs



class HiGPUsInterfaceTests(TestWithMPI):
    
    def test0(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        instance.initialize_code()
        instance.cleanup_code()
        instance.stop()
    
    def test1(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        instance.initialize_code()
        instance.set_number_of_Threads(32)
        instance.set_number_of_GPU(1)
        error = instance.commit_parameters()
        index, error = instance.new_particle(1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        self.assertEquals(error, 0)
        self.assertEquals(index, 0)
        index, error = instance.new_particle(0.000003003, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0)
        self.assertEquals(error, 0)
        self.assertEquals(index, 1)
        error = instance.commit_particles()
        self.assertEquals(error, 0)
        retrieved_state = instance.get_state(index)
        self.assertEquals(retrieved_state['__result'], 0)
        self.assertEquals(0.000003003,  retrieved_state['mass'])
        self.assertEquals(0.0, retrieved_state['radius'])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 2)
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        instance.initialize_code()
        for x in [0.101, 4.0]:
            error = instance.set_eps(x)
            self.assertEquals(error, 0)            
            value, error = instance.get_eps()
            self.assertEquals(error, 0)
            self.assertEquals(x, value)
        instance.cleanup_code()
        instance.stop()

    def test3(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        instance.initialize_code()
        instance.set_number_of_Threads(32)
        instance.set_number_of_GPU(1)
        error = instance.commit_parameters()
        
        instance.new_particle([11.0,12.0,13.0,14.0],
            [2.0,3.0,4.0,5.0],
            [2.1,3.1,4.1,5.1],
            [2.2,3.2,4.2,5.2],
            [2.3,3.3,4.3,5.3],
            [2.4,3.4,4.4,5.4],
            [2.5,3.5,4.5,5.5],
            [2.6,3.6,4.6,5.6])
        error = instance.commit_particles()
        retrieved_state = instance.get_state(0)
        self.assertEquals(11.0,  retrieved_state['mass'])
        retrieved_state = instance.get_state([2,3,4])
        self.assertEquals(14.0,  retrieved_state['mass'][1])
        self.assertEquals(instance.get_number_of_particles()['number_of_particles'], 4)
        instance.cleanup_code()
        instance.stop()
    
    def test4(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        instance.initialize_code()
        instance.set_number_of_Threads(128)
        instance.set_number_of_GPU(1)
        error = instance.commit_parameters()
        
        n = 4000
        ids = [i for i in range(1,n)]
        values = [1.0 * i for i in range(1,n)]
        instance.new_particle(values, values, values, values, values, values, values, values)
        error = instance.commit_particles()
        retrieved_state = instance.get_state(0)
        self.assertEquals(1.0,  retrieved_state['mass'])
        retrieved_state = instance.get_state(3998)
        self.assertEquals(3999.0,  retrieved_state['mass'])
        instance.cleanup_code()
        instance.stop()
        
    def test5(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        instance.initialize_code()
        instance.set_number_of_Threads(32)
        instance.set_number_of_GPU(1)
        instance.set_eps(0.0)
        instance.set_gpu_name("GeForce GTX 480")
        error = instance.commit_parameters()
        
        instance.new_particle( 
            [1.0,1.0,1.0],
            [0.0,0.0,0.0],
            [1.0,0.0,-1.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0],
            [0.0,1.0,0.0],
            [0.0,0.0,0.0],
            [0.0,0.0,0.0] )
        instance.commit_particles()
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
        self.assertEqual( Ek, 0.5)
        self.assertEqual( Ep, -2.5)    
        instance.delete_particle(1)
        instance.recommit_particles()
        n=instance.get_number_of_particles()['number_of_particles']
        Ep=instance.get_potential_energy()['potential_energy']
        Ek=instance.get_kinetic_energy()['kinetic_energy']
        
        self.assertEqual( n, 2)
        self.assertEqual( Ek, 0.)
        self.assertEqual( Ep, -0.5)
        instance.cleanup_code()
        instance.stop()
    
    def test6(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        instance.initialize_code()
        instance.set_number_of_Threads(32)
        instance.set_number_of_GPU(1)
        instance.set_eps(0.1)
        instance.set_gpu_name("GeForce GTX 480")
        error = instance.commit_parameters()
        
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        id2,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -10.0 / (2.0**2 + 0.1**2)**0.5)
        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 10.0]) / 2.0)
        instance.cleanup_code()
        instance.stop()
        
    
    def test7(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        instance.initialize_code()
        instance.set_number_of_Threads(32)
        instance.set_number_of_GPU(1)
        instance.set_eps(0.0)
        instance.set_gpu_name("GeForce GTX480")
        error = instance.commit_parameters()
        
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        id2,errorcode = instance.new_particle(mass = 1.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -1.0 / numpy.sqrt(2.0**2), 8)
        
        potential, errorcode = instance.get_potential(id2)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -10.0 / numpy.sqrt(2.0**2), 8)
        
        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * [10.0, 1.0]) / 2.0)
        instance.cleanup_code()
        instance.stop()
    


