import numpy

from amuse.test.amusetest import TestWithMPI
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.community.higpus.interface import HiGPUsInterface, HiGPUs
from amuse.ic.plummer import new_plummer_model



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
            error = instance.set_eta4(x)
            self.assertEquals(error, 0)            
            value, error = instance.get_eta4()
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
            [2.6,3.6,4.6,5.6],
            [0.001,0.001,0.001,0.001])
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
        instance.set_gpu_name("GeForce GTX 480")
        error = instance.commit_parameters()
        
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0, soft = 0.0)
        id2,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0, soft = 0.0)
        
        instance.commit_particles()
        potential, errorcode = instance.get_potential(id1)
        self.assertEquals(errorcode, 0)
        self.assertAlmostRelativeEquals(potential,  -10.0 / (2.0**2 + 0.1**2)**0.5, 2)
        total_potential, errorcode = instance.get_potential_energy()
        potentials, errorcode = instance.get_potential([id1, id2])
        
        self.assertAlmostRelativeEquals(total_potential, numpy.sum(potentials * 10.0) / 2.0)
        instance.cleanup_code()
        instance.stop()
        
    
    def test7(self):
        instance = self.new_instance_of_an_optional_code(HiGPUsInterface)
        instance.initialize_code()
        instance.set_number_of_Threads(32)
        instance.set_number_of_GPU(1)
        instance.set_gpu_name("GeForce GTX 480")
        error = instance.commit_parameters()
        
        id1,errorcode = instance.new_particle(mass = 10.0, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0, soft = 0.0)
        id2,errorcode = instance.new_particle(mass = 1.0, radius = 1.0, x = 2.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0, soft = 0.0)
        
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
    


class TestHiGPUs(TestWithMPI):
    def new_system_of_sun_and_earth(self):
        stars = datamodel.Stars(2)
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)
        sun.soft =  units.m (0.0)

        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371.) 
        earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(numpy.array((0.0,29800.,0.0)))
        earth.soft = units.m (0.0)
        return stars


    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        instance = self.new_instance_of_an_optional_code(HiGPUs, convert_nbody)
        instance.initialize_code()
    
        instance.parameters.eta_6 = 0.01
        instance.commit_parameters()
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
    
        instance.particles.add_particles(stars)
        instance.commit_particles()
        instance.evolve_model(365 | units.day)
    
        instance.particles.copy_values_of_all_attributes_to(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        instance.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
                
        instance.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        instance.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        instance.cleanup_code()
        
        instance.stop()

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = self.new_instance_of_an_optional_code(HiGPUs, convert_nbody)

        instance.initialize_code()
        instance.commit_parameters()
    
        particles = datamodel.Particles(2)
        self.assertEquals(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEquals(len(instance.particles), 2)
         
        instance.particles.mass =  [17.0, 33.0] | units.kg
        
        instance.commit_particles()

        self.assertEquals(instance.get_mass(0), 17.0| units.kg) 
        self.assertEquals(instance.get_mass(1), 33.0| units.kg)  
        
        instance.cleanup_code()
        instance.stop()
        
    def test3(self):
        instance = self.new_instance_of_an_optional_code(HiGPUs)
        instance.initialize_code()
        instance.commit_parameters()  
        particles = datamodel.Particles(6)
        particles.mass = nbody_system.mass.new_quantity(range(1,7))
        particles.radius =   0.00001 | nbody_system.length
        particles.position = [[-1.0,0.0,0.0],[1.0,0.0,0.0],[0.0,-1.0,0.0],[0.0,1.0,0.0],[0.0,0.0,-1.0],[0.0,0.0,1.0]] | nbody_system.length
        particles.velocity = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]] | nbody_system.speed
        particles.soft = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0] | nbody_system.length
        instance.particles.add_particles(particles)
        
        instance.commit_particles()
        copyof = instance.particles.copy()
        
        instance.cleanup_code()
        instance.stop()
        
        self.assertEquals(2 | nbody_system.mass, copyof[1].mass)  
        
        
    def test4(self):
        self.skip("Need to implement get_gravity_at_point for HiGPUs")
        instance = self.new_instance_of_an_optional_code(HiGPUs)
        instance.initialize_code()
        instance.commit_parameters()  

        particles = datamodel.Particles(2)
        particles.mass = [1.0, 1.0] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        particles.soft = [0.0, 0.0] | nbody_system.length
        instance.particles.add_particles(particles)
        
        instance.commit_particles()

        
        zero = 0.0 | nbody_system.length
        fx, fy, fz = instance.get_gravity_at_point(0.0000001 | nbody_system.length, 2.0 | nbody_system.length, zero, zero)
        self.assertAlmostEqual(fx, -0.25 | nbody_system.acceleration)
        self.assertAlmostEqual(fy, 0.0 | nbody_system.acceleration)
        self.assertAlmostEqual(fz, 0.0 | nbody_system.acceleration)
        for x in (0.25, 0.5, 0.75):
            x0 = x | nbody_system.length
            x1 = (2.0 - x) | nbody_system.length
            potential0 = instance.get_potential_at_point(0.0000001 | nbody_system.length, x0, zero, zero)
            potential1 = instance.get_potential_at_point(0.0000001 | nbody_system.length, x1, zero, zero)
            fx0, fy0, fz0 = instance.get_gravity_at_point(0.0000001 | nbody_system.length, x0, zero, zero)
            fx1, fy1, fz1 = instance.get_gravity_at_point(0.0000001 | nbody_system.length, x1, zero, zero)
            
            self.assertAlmostEqual(fy0, 0.0 | nbody_system.acceleration)
            self.assertAlmostEqual(fz0, 0.0 | nbody_system.acceleration)
            self.assertAlmostEqual(fy1, 0.0 | nbody_system.acceleration)
            self.assertAlmostEqual(fz1, 0.0 | nbody_system.acceleration)
            
            self.assertAlmostEqual(fx0, -1.0 * fx1, 5)
            fx = (-1.0 / (x0**2) + 1.0 / (x1**2)) * (1.0 | nbody_system.length ** 3 / nbody_system.time ** 2)
            self.assertAlmostEqual(fx, fx0, 2)
            self.assertAlmostEqual(potential0, potential1, 5)

        instance.cleanup_code()
        instance.stop()
        
    def test5(self):
        instance = self.new_instance_of_an_optional_code(HiGPUs)
        instance.initialize_code()
    
        instance.parameters.eta_6 = 0.01
        instance.commit_parameters()

        stars = new_plummer_model(100)
        
        instance.particles.add_particles(stars)
        instance.commit_particles()
        
        instance.evolve_model(0.001 | nbody_system.time)
    
        e0 = instance.get_kinetic_energy() + instance.get_potential_energy()
        
        stars.mass *= 0.9
        
        instance.synchronize_model()
        
        e1 = instance.get_kinetic_energy() + instance.get_potential_energy()
        
        instance.cleanup_code()
        instance.stop()
        
        delta_e = e1 - e0
        
        self.assertTrue(e1 != e0)
    
    
    def test6(self):
        particles = datamodel.Particles(2)
        particles.mass = [1.0, 0.000003003] | nbody_system.mass
        particles.radius =  [0.0001, 0.0001] | nbody_system.length
        particles.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 1.0, 0.0]] | nbody_system.speed
        particles.soft = [0.0, 0.0] | nbody_system.length

        particle = datamodel.Particles(1)
        particle.mass = 0.000003003 | nbody_system.mass
        particle.radius =  0.0001 | nbody_system.length
        particle.position = [0.0, 1.0, 0.0] | nbody_system.length
        particle.velocity = [-1.0, 0.0, 0.0] | nbody_system.speed
        particles.soft = 0.0 | nbody_system.length
        
        instance = self.new_instance_of_an_optional_code(HiGPUs)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.particles.add_particles(particles)
        instance.commit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.evolve_model(1 | nbody_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.particles.remove_particle(particles[1])
        instance.particles.add_particle(particle)
        self.assertEquals(instance.get_name_of_current_state(), 'UPDATE')
        instance.recommit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.evolve_model(1 | nbody_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.synchronize_model()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()


    def test7(self):
       
        instance = self.new_instance_of_an_optional_code(HiGPUs)
        instance.parameters.eta_6 = 0.5
        
        instance.commit_parameters()        
        
        self.assertEquals( instance.parameters.eta_6 ,  0.5 )
        self.assertEquals( instance.parameters.eta_4 ,  0.01)
        self.assertEquals( instance.parameters.begin_time ,  0.0 | nbody_system.time) 
        self.assertEquals( instance.parameters.r_core_plummer ,  0.0 | nbody_system.length)
        self.assertEquals( instance.parameters.mass_plummer ,  0.0 | nbody_system.mass)
        self.assertEquals( instance.parameters.Threads ,  128)
        self.assertEquals( instance.parameters.n_Print ,  1000000)
        self.assertEquals( instance.parameters.dt_Print ,  1000000.0 | nbody_system.time)
        self.assertEquals( instance.parameters.max_step ,  pow(2.,-3.0) | nbody_system.time)
        self.assertEquals( instance.parameters.min_step ,  pow(2.,-30.0) | nbody_system.time)
        self.assertEquals( instance.parameters.gpu_name ,  "")
        self.assertEquals( instance.parameters.output_path_name ,  "./data/")
        self.assertEquals( instance.parameters.n_gpu ,  1)
        instance.cleanup_code()
        instance.stop()
