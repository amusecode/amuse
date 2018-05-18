from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.fdps.interface import FDPSInterface, FDPS

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


class TestFDPSInterface(TestWithMPI):
    
    
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


        instance = self.new_instance_of_an_optional_code(FDPSInterface, **default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals([0, 0], instance.get_number_of_particles().values())
        ids, errors = instance.new_particle(mass,x,y,z,vx,vy,vz,radius)
        self.assertEquals(0, errors)
        self.assertEquals(range(plummer_size), ids)
        self.assertEquals(0, instance.commit_particles())
        
        self.assertEquals([500, 0], instance.get_number_of_particles().values())
        masses, errors = instance.get_mass(range(500))
        self.assertEquals(0, errors)
        self.assertAlmostEquals(0.002, masses)
        masses,xs,ys,zs,vxs,vys,vzs,radii, errors = instance.get_state(range(500))
        self.assertEquals(0, errors)
        self.assertAlmostRelativeEquals(xs, x, 6)
        
        self.assertEquals(0, instance.evolve_model(0.00001))
        energy_total_init = instance.get_potential_energy()["potential_energy"] + instance.get_kinetic_energy()["kinetic_energy"]
        self.assertEquals(0, instance.evolve_model(1))
        energy_total_final = instance.get_potential_energy()["potential_energy"] + instance.get_kinetic_energy()["kinetic_energy"]
        self.assertAlmostRelativeEqual(energy_total_init, energy_total_final, 2)
        instance.stop()


class TestFDPS(TestWithMPI):
    
    def test1(self):
        print "Testing FDPS initialization"
        instance = self.new_instance_of_an_optional_code(FDPS, **default_options)
        instance.initialize_code()
        print "cleaning the code"
        instance.cleanup_code()
        print "done"
        instance.stop()
    
    def test2(self):
        print "Testing FDPS parameters"
        instance = self.new_instance_of_an_optional_code(FDPS, **default_options)
        instance.initialize_code()
        self.assertAlmostEquals(instance.parameters.epsilon_squared, 0.0025 | nbody_system.length**2)
        self.assertAlmostEquals(instance.parameters.timestep, 1.0 / 64 | nbody_system.time)
        instance.parameters.epsilon_squared = 0.01 | nbody_system.length**2
        instance.parameters.timestep = 0.001 | nbody_system.time
        instance.commit_parameters()
        self.assertAlmostEquals(instance.parameters.epsilon_squared, 0.01 | nbody_system.length**2)
        self.assertAlmostEquals(instance.parameters.timestep, 0.001 | nbody_system.time)
        instance.stop()
    
    def test3(self):
        print "Testing FDPS, N-body units"
        instance = self.new_instance_of_an_optional_code(FDPS, **default_options)
        instance.initialize_code()
        plummer = new_plummer_model(500)
        instance.particles.add_particles(plummer)
        self.assertAlmostEquals(instance.particles.mass, 0.002 | nbody_system.mass)
        instance.evolve_model(1.0 | nbody_system.time)
        self.assertAlmostEquals(instance.model_time, 1.0 | nbody_system.time)
        instance.stop()
    
    def test4(self):
        print "Testing FDPS, SI units"
        convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1.0 | units.parsec)
        instance = self.new_instance_of_an_optional_code(FDPS, convert_nbody, **default_options)
        instance.initialize_code()
        plummer = new_plummer_model(500, convert_nbody = convert_nbody)
        instance.particles.add_particles(plummer)
        instance.commit_particles()
        instance.evolve_model(1 | nbody_system.time)
        instance.stop()
    
    def test5(self):
        print "Testing FDPS remove_particle"
        convert_nbody = nbody_system.nbody_to_si(100.0 | units.MSun, 1.0 | units.parsec)
        instance = self.new_instance_of_an_optional_code(FDPS, convert_nbody, **default_options)
        instance.initialize_code()
        plummer = new_plummer_model(500, convert_nbody = convert_nbody)
        instance.particles.add_particles(plummer)
        instance.commit_particles()
        self.assertEquals(instance.get_number_of_particles(), 500)
        self.assertEquals(instance.particles[0].key, plummer[0].key)
        self.assertAlmostEquals(instance.particles[0].x, plummer[0].x)
        instance.particles.remove_particle(plummer[0])
        self.assertEquals(instance.get_number_of_particles(), 499)
        self.assertAlmostEquals(instance.particles[0].x, plummer[1].x)
        instance.stop()
    
    def test6(self):
        print "Testing FDPS states"
        plummer = new_plummer_model(500)
        
        print "First do everything manually:"
        instance = self.new_instance_of_an_optional_code(FDPS, **default_options)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.particles.add_particles(plummer)
        instance.commit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        mass = instance.particles[0].mass
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()

        print ("initialize_code(), commit_parameters(), commit_particles(), and cleanup_code() should be called " +
            "automatically before editing parameters, new_particle(), get_xx(), and stop():")
        instance = self.new_instance_of_an_optional_code(FDPS, **default_options)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        timestep = instance.parameters.timestep
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particles(plummer)
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        mass = instance.particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'STOPPED')
    
