
from amuse.test.amusetest import TestWithMPI

from amuse.legacy.bhtree.interface import BHTree
from amuse.legacy.hermite0.interface import Hermite
from amuse.legacy.phiGRAPE.interface import PhiGRAPE
from amuse.legacy.fi.interface import Fi

from amuse.support.data import core
from amuse.support.units import units
from amuse.support.units import nbody_system

from amuse.ext.plummer import new_plummer_sphere

import numpy

class TestGravityCodes(TestWithMPI):
    
    def gravity_code_factory(self):
        self.skip("abstract test")
        
    def test1(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        particles = new_plummer_sphere(100)
        particles.radius = 0 | nbody_system.length
        particles.move_to_center()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        self.assertEquals(len(instance.particles), 100)
        outer = particles.select(lambda r: r.length() > (0.5 | nbody_system.length), ["position"])
        print len(outer)
        self.assertTrue(len(outer) > 0)
        self.assertTrue(len(outer) < 100)
        instance.synchronize_model()
        instance.particles.remove_particles(outer)
        self.assertEquals(len(instance.particles), 100-len(outer))
        number_of_particles_in_module = instance.legacy_interface.get_number_of_particles()['number_of_particles']
        self.assertEquals(number_of_particles_in_module, 100-len(outer))
        instance.stop()
            
    def getset_attribute(self, attributename, attributevalue1, attributevalue2):
        
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        try:
            particles = new_plummer_sphere(100)
            particles.move_to_center()
            particles.radius = 0.0 | nbody_system.length
            instance.particles.add_particles(particles)
            instance.commit_particles()
            setattr(instance.particles, attributename, attributevalue1)
            self.assertTrue(numpy.all(getattr(instance.particles, attributename) == attributevalue1), "attribute '{0}' not equal to expected value '{1}'".format(attributename, attributevalue1))
            setattr(instance.particles, attributename, attributevalue2)
            self.assertTrue(numpy.all(getattr(instance.particles, attributename) == attributevalue2), "attribute '{0}' not equal to expected value '{1}'".format(attributename, attributevalue2))
        finally:
            instance.stop()
            
    def test2(self):
        self.getset_attribute("radius", 1.0 | nbody_system.length, 2.0 | nbody_system.length)
        
    def test3(self):
        self.getset_attribute("position", [1.0, 2.0, 3.0] | nbody_system.length, [5.0, 6.0, 7.0] | nbody_system.length)
    
    def test4(self):
        self.getset_attribute("velocity", [1.0, 2.0, 3.0] | nbody_system.speed, [5.0, 6.0, 7.0] | nbody_system.speed)
    
    def test5(self):
        self.getset_attribute("mass", 1.0 | nbody_system.mass, 2.0 | nbody_system.mass)
        
    def test6(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        try:
            particles = new_plummer_sphere(100)
            more_particles = new_plummer_sphere(50)
            particles.radius = 0 | nbody_system.length
            more_particles.radius = 1.0 | nbody_system.length
            particles.move_to_center()
            more_particles.move_to_center()
            
            instance.particles.add_particles(particles)
            instance.commit_particles()
            self.assertEquals(len(instance.particles), 100)
            #instance.synchronize_model()
            instance.particles.add_particles(more_particles)
            number_of_particles_in_module = instance.legacy_interface.get_number_of_particles()['number_of_particles']
            self.assertEquals(len(instance.particles), 150)
            self.assertEquals(number_of_particles_in_module, 150)
            instance.synchronize_model()
            instance.particles.remove_particles(particles)
            number_of_particles_in_module = instance.legacy_interface.get_number_of_particles()['number_of_particles']
            self.assertEquals(len(instance.particles), 50)
            self.assertEquals(number_of_particles_in_module, 50) 
        finally:
            instance.stop()
    


class TestBHTreeGravityCode(TestGravityCodes):
    pass
    
    def gravity_code_factory(self):
        return BHTree


class TestHermiteGravityCode(TestGravityCodes):
    pass
    
    def gravity_code_factory(self):
        return Hermite
        
class TestPhiGRAPEGravityCode(TestGravityCodes):
    pass
    
    def gravity_code_factory(self):
        return PhiGRAPE

class TestFiGravityCode(TestGravityCodes):
    pass
    
    def gravity_code_factory(self):
        return Fi
