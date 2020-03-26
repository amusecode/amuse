import os
import os.path
import math

from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from amuse.units import units, nbody_system, constants
from amuse.datamodel import Particles

from amuse.community.brutus.interface import BrutusInterface, Brutus

import random

try:
    import mpmath
    HAS_MPMATH=True
except ImportError:
    HAS_MPMATH=False


class TestBrutusInterface(TestWithMPI):
    
    def test1(self):
        print("Test BrutusInterface initialization")
        instance =  self.new_instance_of_an_optional_code(BrutusInterface)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.set_brutus_output_directory(instance.output_directory))
        self.assertEqual(0, instance.commit_parameters())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test2(self):
        print("Test BrutusInterface new_particle / get_state")
        instance =  self.new_instance_of_an_optional_code(BrutusInterface)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.set_brutus_output_directory(instance.output_directory))
        self.assertEqual(0, instance.commit_parameters())
        
        id, error = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        self.assertEqual(0, error)
        self.assertEqual(0, id)
        id, error = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        self.assertEqual(0, error)
        self.assertEqual(1, id)
        self.assertEqual(0, instance.commit_particles())
        
        retrieved_state1 = instance.get_state(0)
        retrieved_state2 = instance.get_state(1)
        self.assertEqual(0,  retrieved_state1['__result'])
        self.assertEqual(0,  retrieved_state2['__result'])
        self.assertEqual(11.0,  retrieved_state1['mass'])
        self.assertEqual(21.0,  retrieved_state2['mass'])
        self.assertEqual( 0.0,  retrieved_state1['x'])
        self.assertEqual(10.0,  retrieved_state2['x'])
        
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test4(self):
        print("Test BrutusInterface particle property getters/setters")
        instance =  self.new_instance_of_an_optional_code(BrutusInterface)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.set_brutus_output_directory(instance.output_directory))
        self.assertEqual(0, instance.commit_parameters())

        self.assertEqual([0, 0], list(instance.new_particle(0.01,  1, 0, 0,  0, 1, 0, 0.1).values()))
        self.assertEqual([1, 0], list(instance.new_particle(0.02, -1, 0, 0,  0,-1, 0, 0.1).values()))
        self.assertEqual(0, instance.commit_particles())
        
        # getters
        mass, result = instance.get_mass(0)
        self.assertAlmostEqual(0.01, mass)
        self.assertEqual(0,result)
        radius, result = instance.get_radius(1)
        self.assertAlmostEqual(0.1, radius)
        self.assertEqual(0,result)
        #self.assertEquals(-3, instance.get_mass(2)['__result']) # Particle not found
        self.assertEqual([ 1, 0, 0,  0], list(instance.get_position(0).values()))
        self.assertEqual([-1, 0, 0,  0], list(instance.get_position(1).values()))
        self.assertEqual([ 0, 1, 0,  0], list(instance.get_velocity(0).values()))
        self.assertEqual([ 0,-1, 0,  0], list(instance.get_velocity(1).values()))
        
        # setters
        self.assertEqual(0, instance.set_state(0, 0.01, 1,2,3, 4,5,6, 0.1))
        self.assertEqual([0.01, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.1, 0], list(instance.get_state(0).values()))
        self.assertEqual(0, instance.set_mass(0, 0.02))
        self.assertEqual([0.02, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.1, 0], list(instance.get_state(0).values()))
        self.assertEqual(0, instance.set_radius(0, 0.2))
        self.assertEqual([0.02, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.2, 0], list(instance.get_state(0).values()))
        self.assertEqual(0, instance.set_position(0, 10,20,30))
        self.assertEqual([0.02, 10.0,20.0,30.0, 4.0,5.0,6.0, 0.2, 0], list(instance.get_state(0).values()))
        self.assertEqual(0, instance.set_velocity(0, 40,50,60))
        self.assertEqual([0.02, 10.0,20.0,30.0, 40.0,50.0,60.0, 0.2, 0], list(instance.get_state(0).values()))

        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test5(self):
        print("Test BrutusInterface parameters")
        instance =  self.new_instance_of_an_optional_code(BrutusInterface)
        self.assertEqual(0, instance.initialize_code())
        
        # word length
        self.assertEqual([64, 0], list(instance.get_word_length().values()))
        self.assertEqual(0, instance.set_word_length(80))
        self.assertEqual([80, 0], list(instance.get_word_length().values()))

        # bs tolerance, default (double) implementation
        self.assertEqual([1.0e-6, 0], list(instance.get_bs_tolerance().values()))
        self.assertEqual(0, instance.set_bs_tolerance(1.0e-8))
        self.assertEqual([1.0e-8, 0], list(instance.get_bs_tolerance().values()))

        # bs tolerance, string implementation for values requiring higher precision (note: actual accuracy depends on word_length)
        #self.assertEquals(1e-8, eval(instance.get_bs_tolerance_string()[""]))
        #self.assertEquals(0, instance.set_bs_tolerance_string("1e-10"))
        #self.assertEquals(["1e-10", 0], instance.get_bs_tolerance_string().values())

        # eta, float64
        self.assertEqual([0.24, 0], list(instance.get_eta().values()))
        self.assertEqual(0, instance.set_eta(0.10))
        self.assertEqual([0.10, 0], list(instance.get_eta().values()))

        # eta, string
        #self.assertEquals(["0.10", 0], instance.get_eta_string().values())
        self.assertEqual(0, instance.set_eta_string("123"))
        self.assertEqual(["123", 0], list(instance.get_eta_string().values()))

        # output dir
        #self.assertEquals(["./", 0], instance.get_brutus_output_directory().values())
        self.assertEqual(0, instance.set_brutus_output_directory("./out"))
        self.assertEqual(["./out/", 0], list(instance.get_brutus_output_directory().values()))
        self.assertEqual(0, instance.set_brutus_output_directory(instance.output_directory))
        self.assertEqual([instance.output_directory+"/", 0], list(instance.get_brutus_output_directory().values()))
        
        self.assertEqual(0, instance.commit_parameters())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test6(self):
        print("Test BrutusInterface evolve_model, equal-mass binary")

        instance =  self.new_instance_of_an_optional_code(BrutusInterface)
        self.assertEqual(0, instance.initialize_code())

        self.assertEqual(0, instance.set_bs_tolerance(1.0e-10))
        self.assertEqual(0, instance.set_word_length(72))

        self.assertEqual(0, instance.commit_parameters())
        
        self.assertEqual([0, 0], list(instance.new_particle(0.5,  0.5, 0, 0,  0, 0.5, 0).values()))
        self.assertEqual([1, 0], list(instance.new_particle(0.5, -0.5, 0, 0,  0,-0.5, 0).values()))
        self.assertEqual(0, instance.commit_particles())

        self.assertEqual(0, instance.evolve_model(math.pi)) # half an orbit
        for result, expected in zip(instance.get_position(0).values(), [-0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 5)

        self.assertEqual(0, instance.evolve_model(2 * math.pi)) # full orbit
        for result, expected in zip(instance.get_position(0).values(), [0.5, 0.0, 0.0, 0]):
            self.assertAlmostEqual(result, expected, 5)

        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test7(self):
        print("Test BrutusInterface evolve_model, pythagorean problem")

        instance =  self.new_instance_of_an_optional_code(BrutusInterface)
        self.assertEqual(0, instance.initialize_code())

        self.assertEqual(0, instance.set_bs_tolerance(1.0e-6))
        self.assertEqual(0, instance.set_word_length(56))

        self.assertEqual(0, instance.commit_parameters())

        self.assertEqual([0, 0], list(instance.new_particle("3",  "1",  "3", "0", "0", "0", "0").values()))
        self.assertEqual([1, 0], list(instance.new_particle("4", "-2", "-1", "0", "0", "0", "0").values()))
        self.assertEqual([2, 0], list(instance.new_particle("5",  "1", "-1", "0", "0", "0", "0").values()))

        self.assertEqual(0, instance.commit_particles())

        self.assertEqual(0, instance.evolve_model(10))
        
        ## add a check for assertequal final coordinates
        for result, expected in zip(instance.get_position(0).values(), [0.778480410138085492274810667212415, 0.141392300290086165745727207379442, 0, 0]):
            self.assertAlmostEqual(result, expected, 3)

        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test8(self):
        print("Test BrutusInterface string parameters")
        instance =  self.new_instance_of_an_optional_code(BrutusInterface)
        instance.initialize_code()
        
        instance.set_word_length(128)
        
        for i in range(100):
          x=random.random()
          x=str(x)
          instance.set_eta_string(x)
          x_,err=instance.get_eta_string()
          instance.set_eta_string(x_)
          x__,err=instance.get_eta_string()
          #~ assert x==x_ 
          self.assertEqual(x_,x__) 
        
        instance.stop()
        

class TestBrutus(TestWithMPI):
    
    def new_sun_earth_system(self):
        particles = Particles(2)
        particles.mass = [1.0, 3.0037e-6] | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * particles.total_mass() / (1.0 | units.AU)).sqrt()
        return particles

    def test1(self):
        print("Testing Brutus initialization")
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance =  self.new_instance_of_an_optional_code(Brutus, convert_nbody)
        instance.initialize_code()
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()

    def test2(self):
        print("Testing Brutus parameters")

        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)

        instance = self.new_instance_of_an_optional_code(Brutus,convert_nbody)
        instance.initialize_code()
        
#        print instance.parameters
        self.assertEqual(instance.parameters.bs_tolerance, 1.0e-6)
        instance.parameters.bs_tolerance = 1.0e-9
        self.assertEqual(instance.parameters.bs_tolerance, 1.0e-9)
        
        self.assertEqual(instance.parameters.word_length, 64)
        instance.parameters.word_length = 128
        self.assertEqual(instance.parameters.word_length, 128)
        
        self.assertEqual(instance.parameters.dt_param, 0.24)
        instance.parameters.dt_param = 0.10
        self.assertEqual(instance.parameters.dt_param, 0.10)        

        self.assertEqual(instance.parameters.brutus_output_directory, instance.output_directory + os.sep)
        instance.parameters.brutus_output_directory = "./out"
        self.assertEqual(instance.parameters.brutus_output_directory, "./out/")
        instance.parameters.brutus_output_directory = instance.output_directory
        self.assertEqual(instance.parameters.brutus_output_directory, instance.output_directory + os.sep)
                
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):
        print("Testing Brutus particles")

        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)

        instance = self.new_instance_of_an_optional_code(Brutus,convert_nbody)
        instance.initialize_code()

        instance.commit_parameters()

        instance.particles.add_particles(self.new_sun_earth_system())
        instance.commit_particles()
        
        self.assertAlmostEqual(instance.particles.mass, [1.0, 3.0037e-6] | units.MSun)
        self.assertAlmostEqual(instance.particles.radius, 1.0 | units.RSun)
        self.assertAlmostEqual(instance.particles.position, 
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU)
        self.assertAlmostEqual(instance.particles.velocity, 
            [[0.0, 0.0, 0.0], [0.0, 29.7885, 0.0]] | units.km / units.s, 3)
        
        instance.cleanup_code()
        instance.stop()
 
    def test4(self):
        print("Testing Brutus evolve_model, 2 particles")

        particles = Particles(2)
        particles.mass = 0.5 | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * (1.0 | units.MSun) / (1.0 | units.AU)).sqrt()
        particles.move_to_center()
        
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)

        instance = self.new_instance_of_an_optional_code(Brutus, convert_nbody)
        instance.initialize_code()

        instance.parameters.bs_tolerance = 1e-6
        instance.parameters.word_length = 56

        instance.commit_parameters()

        instance.particles.add_particles(particles)

        instance.commit_particles()

        primary = instance.particles[0]
        
        P = 2 * math.pi * primary.x / primary.vy
        
        position_at_start = primary.position.x
        instance.evolve_model(P / 4.0)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.y, 6)
        
        instance.evolve_model(P / 2.0)
        self.assertAlmostRelativeEqual(position_at_start, -primary.position.x, 6)
        
        instance.evolve_model(P)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.x, 6)
        
        instance.cleanup_code()
        instance.stop()
        
    def sun_and_planets(self):
        particles = Particles(9)
        sun = particles[0] 
        mercury = particles[1]
        venus = particles[2]
        earth = particles[3]
        mars = particles[4]
        jupiter = particles[5]
        saturn = particles[6]
        uranus = particles[7]
        neptune = particles[8]
        
        sun.mass = 1047.517| units.MJupiter  
        sun.radius = 1.0 | units.RSun      
        sun.position = ( 0.005717 , -0.00538 , -2.130e-5 ) | units.AU
        sun.velocity = ( 0.007893 , 0.01189 , 0.0002064 )  | units.kms 
            
        mercury.mass = 0.000174 | units.MJupiter  
        mercury.radius =  0  | units.RSun    
        mercury.position = ( -0.31419 , 0.14376 , 0.035135 ) | units.AU 
        mercury.velocity = ( -30.729 , -41.93 , -2.659 )  | units.kms  
            
        venus.mass = 0.002564 | units.MJupiter     
        venus.radius =   0    | units.RSun 
        venus.position = ( -0.3767 , 0.60159 , 0.0393 ) | units.AU 
        venus.velocity = ( -29.7725 , -18.849 , 0.795 )  | units.kms 
            
        earth.mass = 0.003185 | units.MJupiter     
        earth.radius =   0    | units.RSun 
        earth.position = ( -0.98561 , 0.0762 , -7.847e-5 ) | units.AU 
        earth.velocity = ( -2.927 , -29.803 , -0.000533 )  | units.kms 
            
        mars.mass = 0.000338 | units.MJupiter     
        mars.radius =     0  | units.RSun 
        mars.position = ( -1.2895 , -0.9199 , -0.048494 ) | units.AU 
        mars.velocity = ( 14.9 , -17.721 , 0.2979 )  | units.kms 
            
        jupiter.mass = 1 | units.MJupiter     
        jupiter.radius =  0    | units.RSun  
        jupiter.position = ( -4.9829 , 2.062 , -0.10990 ) | units.AU 
        jupiter.velocity = ( -5.158 , -11.454 , -0.13558 )  | units.kms 
        
        saturn.mass = 0.29947 | units.MJupiter     
        saturn.radius =   0    | units.RSun 
        saturn.position = ( -2.075 , 8.7812 , 0.3273 ) | units.AU 
        saturn.velocity = ( -9.9109 , -2.236 , -0.2398 )  | units.kms 
        
        uranus.mass = 0.045737 | units.MJupiter     
        uranus.radius =    0   | units.RSun 
        uranus.position = ( -12.0872 , -14.1917 , 0.184214 ) | units.AU 
        uranus.velocity = ( 5.1377 , -4.7387 , -0.06108 )  | units.kms 
            
        neptune.mass = 0.053962 | units.MJupiter 
        neptune.radius =   0    | units.RSun 
        neptune.position = ( 3.1652 , 29.54882 , 0.476391 ) | units.AU 
        neptune.velocity = ( -5.443317 , 0.61054 , -0.144172 )  | units.kms 
        
        particles.move_to_center()
        return particles

    def test5(self):
        if not HAS_MPMATH:
            self.skip("mpmath not available")
        print("MPmath available -> Doing tests")
        bodies = self.sun_and_planets()
        convert_nbody = nbody_system.nbody_to_si(bodies.mass.sum(),bodies[1].position.length())
        gravity = Brutus(convert_nbody,number_of_workers=1)
        gravity.parameters.bs_tolerance = 1e-30
        gravity.parameters.word_length = 180
        gravity.parameters.dt_param = 0.0000000000010
        gravity.particles.add_particles(bodies)
        Etot_init = gravity.kinetic_energy + gravity.potential_energy
        Ein = gravity.get_total_energy_p_si()
        gravity.evolve_model(gravity.model_time + (30| units.day))
        Eout = gravity.get_total_energy_p_si()
        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        Loss_double = ((Etot_init-Etot)/gravity.get_time())
        Loss_mp = (Ein - Eout)/gravity.get_time_p_si()
        print("Loss with \"normal\" double =",Loss_double.number," (W)")
        print("Loss with multiprecision =",Loss_mp," (W)")
        gravity.stop()
        self.assertTrue((Loss_mp <= 0.0000007) and (Loss_mp > 0.0000006))
    


