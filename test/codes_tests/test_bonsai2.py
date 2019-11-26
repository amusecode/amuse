from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.bonsai2.interface import BonsaiInterface2, Bonsai2

import os
import sys 
import numpy

from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi import channel
from amuse.ic.plummer import new_plummer_model
from amuse.support.exceptions import AmuseException

default_options = dict()
#default_options = dict(redirection="none")
#default_options = dict(redirection="none", debugger="gdb")


class TestBonsaiInterface(TestWithMPI):
    
    
    def test1(self):
        plummer_size = 500
        plummer =  new_plummer_model(plummer_size)
        mass=plummer.mass.number
        radius=plummer.radius.number
        x=plummer.x.number
        y=plummer.y.number
        z=plummer.z.number
        vx=plummer.vx.number
        vy=plummer.vy.number
        vz=plummer.vz.number


        instance = self.new_instance_of_an_optional_code(BonsaiInterface2, **default_options)
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual([0, 0], list(instance.get_number_of_particles().values()))
        ids, errors = instance.new_particle(mass,x,y,z,vx,vy,vz,radius)
        self.assertEqual(0, errors)
        self.assertEqual(list(range(plummer_size)), ids)
        self.assertEqual(0, instance.commit_particles())
        
        self.assertEqual([500, 0], list(instance.get_number_of_particles().values()))
        masses, errors = instance.get_mass(list(range(500)))
        self.assertEqual(0, errors)
        self.assertAlmostEqual(0.002, masses)
        masses,xs,ys,zs,vxs,vys,vzs,radii, errors = instance.get_state(list(range(500)))
        self.assertEqual(0, errors)
        self.assertAlmostRelativeEquals(xs, x, 6)
        
        self.assertEqual(0, instance.evolve_model(0.00001))
        energy_total_init = instance.get_potential_energy()["potential_energy"] + instance.get_kinetic_energy()["kinetic_energy"]
        self.assertEqual(0, instance.evolve_model(1))
        energy_total_final = instance.get_potential_energy()["potential_energy"] + instance.get_kinetic_energy()["kinetic_energy"]
        self.assertAlmostRelativeEqual(energy_total_init, energy_total_final, 2)
        instance.stop()


class TestBonsai(TestWithMPI):
    
    def test1(self):
        print("Testing Bonsai initialization")
        instance = self.new_instance_of_an_optional_code(Bonsai2, **default_options)
        instance.initialize_code()
        print("cleaning the code")
        instance.cleanup_code()
        print("done")
        instance.stop()
    
    def test2(self):
        print("Testing Bonsai parameters")
        instance = self.new_instance_of_an_optional_code(Bonsai2,  **default_options)
        instance.initialize_code()
        self.assertAlmostEqual(instance.parameters.epsilon_squared, 0.0025 | nbody_system.length**2)
        self.assertAlmostEqual(instance.parameters.timestep, 1.0 / 64 | nbody_system.time)
        instance.parameters.epsilon_squared = 0.01 | nbody_system.length**2
        instance.parameters.timestep = 0.001 | nbody_system.time
        instance.commit_parameters()
        print((instance.parameters.epsilon_squared))
        print((instance.parameters.timestep))
        self.assertAlmostEqual(instance.parameters.epsilon_squared, 0.01 | nbody_system.length**2)
        self.assertAlmostEqual(instance.parameters.timestep, 0.001 | nbody_system.time)
        instance.stop()
    
    def test3(self):
        print("Testing Bonsai, N-body units")
        instance = self.new_instance_of_an_optional_code(Bonsai2, **default_options)
        instance.initialize_code()
        plummer = new_plummer_model(500)
        instance.particles.add_particles(plummer)
        self.assertAlmostEqual(instance.particles.mass, 0.002 | nbody_system.mass)
        instance.evolve_model(1.0 | nbody_system.time)
        self.assertAlmostEqual(instance.model_time, 1.0 | nbody_system.time)
        instance.stop()
    
    def test4(self):
        print("Testing Bonsai, SI units")
        convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1.0 | units.parsec)
        instance = self.new_instance_of_an_optional_code(Bonsai2, convert_nbody, **default_options)
        instance.initialize_code()
        plummer = new_plummer_model(500, convert_nbody = convert_nbody)
        instance.particles.add_particles(plummer)
        instance.commit_particles()
        instance.evolve_model(1 | nbody_system.time)
        instance.stop()
    
    def test5(self):
        print("Testing Bonsai remove_particle")
        convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1.0 | units.parsec)
        instance = self.new_instance_of_an_optional_code(Bonsai2, convert_nbody, **default_options)
        instance.initialize_code()
        plummer = new_plummer_model(500, convert_nbody = convert_nbody)
        instance.particles.add_particles(plummer)
        instance.commit_particles()
        instance.particles.remove_particle(plummer[0])
#    self.assertRaises(AmuseException, instance.particles.remove_particle, plummer[0], 
#            expected_message = "Error when calling 'delete_particle' of a 'Bonsai', errorcode is -2")
        print(instance.particles[0].mass)
        instance.stop()
    
    def test6(self):
        print("Testing Bonsai states")
        plummer = new_plummer_model(500)
        
        print("First do everything manually:")
        instance = self.new_instance_of_an_optional_code(Bonsai2, **default_options)
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.particles.add_particles(plummer)
        instance.commit_particles()
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        mass = instance.particles[0].mass
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertEqual(instance.get_name_of_current_state(), 'EVOLVED')
        instance.cleanup_code()
        self.assertEqual(instance.get_name_of_current_state(), 'END')
        instance.stop()

        print(("initialize_code(), commit_parameters(), commit_particles(), and cleanup_code() should be called " +
            "automatically before editing parameters, new_particle(), get_xx(), and stop():"))
        instance = self.new_instance_of_an_optional_code(Bonsai2, **default_options)
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        timestep = instance.parameters.timestep
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particles(plummer)
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        mass = instance.particles[0].mass
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertEqual(instance.get_name_of_current_state(), 'EVOLVED')
        instance.stop()
        self.assertEqual(instance.get_name_of_current_state(), 'STOPPED')
    
    def test7(self):
        print("Testing Bonsai properties")
        numpy.random.seed(12345)
        plummer = new_plummer_model(500)
        instance = self.new_instance_of_an_optional_code(Bonsai2, **default_options)
        instance.particles.add_particles(plummer)
        
        self.assertEqual(instance.model_time, 0.0 | nbody_system.time)
        instance.evolve_model(0.0001 | nbody_system.time)
        print(instance.model_time)
        self.assertTrue(instance.model_time >= 0.0001 | nbody_system.time)
        self.assertEqual(instance.model_time, 1.0 * instance.parameters.timestep)
        
        self.assertAlmostEqual(instance.potential_energy, -0.50625962019 | nbody_system.energy)
        self.assertAlmostEqual(instance.kinetic_energy,    0.244611829519 | nbody_system.energy)
        self.assertAlmostEqual(instance.total_mass, 1.0 | nbody_system.mass)
        
        E0 = instance.kinetic_energy + instance.potential_energy

        self.assertRaises(AmuseException, getattr, instance, "total_radius", expected_message = 
            "Error when calling 'get_total_radius' of a 'Bonsai2', errorcode is -2, "
            "error is 'Called function is not implemented.'")
        self.assertRaises(AmuseException, getattr, instance, "center_of_mass_position", expected_message = 
            "Error when calling 'get_center_of_mass_position' of a 'Bonsai2', "
            "errorcode is -2, error is 'Called function is not implemented.'")
        self.assertRaises(AmuseException, getattr, instance, "center_of_mass_velocity", expected_message = 
            "Error when calling 'get_center_of_mass_velocity' of a 'Bonsai2', "
            "errorcode is -2, error is 'Called function is not implemented.'")

        instance.evolve_model(1.0 | nbody_system.time)
        self.assertAlmostEqual(instance.model_time, 1.0 | nbody_system.time)
        self.assertAlmostEqual(instance.potential_energy, -0.5115 | nbody_system.energy, 4)
        self.assertAlmostEqual(instance.kinetic_energy,    0.24986 | nbody_system.energy, 4)
        self.assertAlmostEqual(
            instance.kinetic_energy + instance.potential_energy, 
            E0, 4)
        instance.particles.remove_particle(plummer[2])
        instance.evolve_model(2.0 | nbody_system.time)
        instance.particles.remove_particle(plummer[20])
        instance.particles.remove_particle(plummer[30])
        instance.evolve_model(3.0 | nbody_system.time)
        instance.particles.remove_particle(plummer[35])
        instance.stop()


"""
    def test8(self):
        print "Testing Bonsai2 collision_detection"
        particles = datamodel.Particles(7)
        particles.mass = 0.001 | nbody_system.mass
        particles.radius = 0.01 | nbody_system.length
        particles.x = [-101.0, -100.0, -0.5, 0.5, 100.0, 101.0, 104.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [[2, 0, 0], [-2, 0, 0]]*3 + [[-4, 0, 0]] | nbody_system.speed
        
        instance = self.new_instance_of_an_optional_code(Bonsai2)
        instance.initialize_code()
        instance.parameters.set_defaults()
        instance.particles.add_particles(particles)
        collisions = instance.stopping_conditions.collision_detection
        collisions.enable()
        instance.evolve_model(1.0 | nbody_system.time)
        
        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 0.5 | nbody_system.time)
        self.assertEquals(len(collisions.particles(0)), 3)
        self.assertEquals(len(collisions.particles(1)), 3)
        self.assertEquals(len(particles - collisions.particles(0) - collisions.particles(1)), 1)
        self.assertEquals(abs(collisions.particles(0).x - collisions.particles(1).x) < 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True, True, True])
        
        sticky_merged = datamodel.Particles(len(collisions.particles(0)))
        sticky_merged.mass = collisions.particles(0).mass + collisions.particles(1).mass
        sticky_merged.radius = collisions.particles(0).radius
        for p1, p2, merged in zip(collisions.particles(0), collisions.particles(1), sticky_merged):
            merged.position = (p1 + p2).center_of_mass()
            merged.velocity = (p1 + p2).center_of_mass_velocity()
        
        print instance.model_time
        print instance.particles
        instance.particles.remove_particles(collisions.particles(0) + collisions.particles(1))
        instance.particles.add_particles(sticky_merged)
        
        instance.evolve_model(1.0 | nbody_system.time)
        print
        print instance.model_time
        print instance.particles
        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 1.0 | nbody_system.time)
        self.assertEquals(len(collisions.particles(0)), 1)
        self.assertEquals(len(collisions.particles(1)), 1)
        self.assertEquals(len(instance.particles - collisions.particles(0) - collisions.particles(1)), 2)
        self.assertEquals(abs(collisions.particles(0).x - collisions.particles(1).x) < 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True])
        instance.stop()
    
    def test9(self):
        print "Testing Bonsai2 tree build exception"
        plummer = new_plummer_model(50)
        #instance = self.new_instance_of_an_optional_code(Bonsai2, **default_options)
        instance = self.new_instance_of_an_optional_code(Bonsai2,debugger="gdb", **default_options)
        instance.particles.add_particles(plummer)
        
        instance.particles[0].position -= [1e9, 0, 0] | nbody_system.length
        self.assertRaises(AmuseException, instance.evolve_model, 1.0 | nbody_system.time, expected_message = 
            "Error when calling 'evolve_model' of a 'Bonsai2', errorcode is -4, error is "
            "'The tree has become too deep, consider the removal of far away particles to prevent a too large box.'")
        
        instance.particles.remove_particle(instance.particles[0])
        instance.evolve_model(0.1 | nbody_system.time)
        instance.stop()
"""    

