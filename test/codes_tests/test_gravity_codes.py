from amuse.test.amusetest import TestWithMPI

from amuse.community.bhtree.interface import BHTree
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi

from amuse.ext.plummer import new_plummer_sphere

import numpy
import time
from amuse.units import units
from amuse.units import nbody_system
from amuse import datamodel
class _TestGravityCodes(TestWithMPI):
    length_unit = nbody_system.length
    speed_unit = nbody_system.speed
    mass_unit = nbody_system.mass
    time_unit = nbody_system.time
    
    @property
    def nbody_converter(self):
        return None
        
    def gravity_code_factory(self):
        self.skip("abstract test")
        
    def test1(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        particles = new_plummer_sphere(100, convert_nbody = self.nbody_converter)
        particles.radius = 0 | self.length_unit
        particles.move_to_center()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        self.assertEquals(len(instance.particles), 100)
        outer = particles.select(lambda r: r.length() > (0.5 | self.length_unit), ["position"])
        print len(outer)
        self.assertTrue(len(outer) > 0)
        self.assertTrue(len(outer) < 100)
        instance.synchronize_model()
        instance.particles.remove_particles(outer)
        instance.recommit_particles()
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
            actual = getattr(instance.particles, attributename)
            self.assertAlmostRelativeEqual(actual, attributevalue1)
            setattr(instance.particles, attributename, attributevalue2)
            actual = getattr(instance.particles, attributename)
            print actual.as_quantity_in(attributevalue2.unit)
            self.assertAlmostRelativeEqual(actual, attributevalue2)
        finally:
            instance.stop()
            
    def test2(self):
        self.getset_attribute("radius", 1.0 | self.length_unit, 2.0 | self.length_unit)
        
    def test3(self):
        self.getset_attribute("position", [1.0, 2.0, 3.0] | self.length_unit, [5.0, 6.0, 7.0] | self.length_unit)
    
    def test4(self):
        self.getset_attribute("velocity", [1.0, 2.0, 3.0] | self.speed_unit, [5.0, 6.0, 7.0] | self.speed_unit)
    
    def test5(self):
        self.getset_attribute("mass", 1.0 | self.mass_unit, 2.0 | self.mass_unit)
        
    def test6(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory)
        try:
            particles = new_plummer_sphere(100, convert_nbody = self.nbody_converter)
            more_particles = new_plummer_sphere(50, convert_nbody = self.nbody_converter)
            particles.radius = 0 | self.length_unit
            more_particles.radius = 1.0 | self.length_unit
            particles.move_to_center()
            more_particles.move_to_center()
            
            instance.particles.add_particles(particles)
            instance.commit_particles()
            self.assertEquals(len(instance.particles), 100)
            instance.synchronize_model()
            instance.particles.add_particles(more_particles)
            instance.recommit_particles()
            number_of_particles_in_module = instance.legacy_interface.get_number_of_particles()['number_of_particles']
            self.assertEquals(len(instance.particles), 150)
            self.assertEquals(number_of_particles_in_module, 150)
            instance.synchronize_model()
            instance.particles.remove_particles(particles)
            self.assertEquals(len(instance.particles), 50)
            instance.recommit_particles()
            instance.synchronize_model()
            number_of_particles_in_module = instance.legacy_interface.get_number_of_particles()['number_of_particles']
            self.assertEquals(len(instance.particles), 50)
            self.assertEquals(number_of_particles_in_module, 50)
        finally:
            instance.stop()
    



    def new_gravity_code(self):
        self.gravity_code_factory()
    
    
class TestBHTreeGravityCode(_TestGravityCodes):
    
    def gravity_code_factory(self):
        return BHTree


class TestHermiteGravityCode(_TestGravityCodes):
    
    def gravity_code_factory(self):
        return Hermite
        
class TestPhiGRAPEGravityCode(_TestGravityCodes):
    
    def gravity_code_factory(self):
        return PhiGRAPE
        
    def test7(self):
        factory = self.gravity_code_factory()
        instance = self.new_instance_of_an_optional_code(factory, mode='gpu')
        try:
            sc = instance.stopping_conditions.collision_detection
            
            n = 1000
            instance.parameters.initialize_gpu_once=1
            sc.enable()
            
            particles = new_plummer_sphere(n, convert_nbody = self.nbody_converter)
            particles.radius = 1.0 / n | self.length_unit
            particles.move_to_center()
            
            instance.particles.add_particles(particles)
            instance.commit_particles()
            self.assertEquals(len(instance.particles), n)
            instance.synchronize_model()
            for t in range(1, 101):
                instance.evolve_model((t * 0.01) | self.time_unit)
                print sc.is_set()
                if sc.is_set():
                    particle1 = sc.particles(0)[0]
                    particle2 = sc.particles(1)[0] 
                    newparticle = datamodel.Particles(1)
                    newparticle.mass = particle1.mass + particle2.mass
                    newparticle.radius = particle2.radius
                    newparticle.position = (particle1.position + particle2.position) /2
                    newparticle.velocity = (particle1.velocity + particle2.velocity) /2
                    instance.particles.remove_particle(particle1)
                    instance.particles.remove_particle(particle2)
                    merged = instance.particles.add_particles(newparticle)
                    print 'Remnant:\n',merged
            
        finally:
            instance.stop()

class TestFiGravityCode(_TestGravityCodes):
    
    def gravity_code_factory(self):
        return Fi

class TestGadget2GravityCode(_TestGravityCodes):
    length_unit = units.parsec
    speed_unit = units.parsec / units.Myr
    mass_unit = units.MSun

    @property
    def nbody_converter(self):
        return nbody_system.nbody_to_si(1.0 | units.parsec, 100.0 | units.MSun)
        
    def gravity_code_factory(self):
        return Gadget2
    
    

    def test2(self):
        self.skip("no support for setting of radius")
    
    
