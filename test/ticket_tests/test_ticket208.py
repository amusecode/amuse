import numpy
from amuse.test import amusetest

import os.path

from amuse.community.phigrape.interface import PhiGRAPEInterface, PhiGRAPE
from amuse.community.hermite.interface import Hermite
from amuse.community.bhtree.interface import BHTree
from amuse.community.ph4.interface import ph4
from amuse.community.fi.interface import Fi
from amuse.community.gadget2.interface import Gadget2
from amuse import io
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
class TestsForTicket208(amusetest.TestCase):
    
    def _run_addition_removal_test(
            self,
            instance, 
            length_unit = nbody_system.length,
            speed_unit = nbody_system.speed,
            mass_unit = nbody_system.mass
        ):
        
        instance.initialize_code()
    
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(1,11) | mass_unit
        particles.radius = numpy.arange(1,11) | length_unit
        particles.x = numpy.arange(1,11) | length_unit
        particles.y = numpy.arange(1,11) | length_unit
        particles.z = numpy.arange(1,11) | length_unit
        particles.vx = numpy.arange(1,11) | speed_unit
        particles.vy = numpy.arange(1,11) | speed_unit
        particles.vz = numpy.arange(1,11) | speed_unit
        
        instance.particles.add_particles(particles)
        instance.commit_particles()
        
        self.assertEqual(len(instance.particles), 10)
        self.assertAlmostRelativeEquals(instance.particles.mass.as_quantity_in(mass_unit), list(numpy.arange(1,11)) | mass_unit)
                
        instance.particles.remove_particle(particles[2])
        instance.particles.remove_particle(particles[5])
        
        self.assertEqual(len(instance.particles), 8)
        self.assertAlmostRelativeEquals(instance.particles.mass.as_quantity_in(mass_unit), [1,2,4,5,7,8,9,10] | mass_unit)
        
        particles_new = datamodel.Particles(1)
        particles_new.mass = 20 | mass_unit
        particles_new.radius = 21 | length_unit
        particles_new.x = 22 | length_unit
        particles_new.y = 23 | length_unit
        particles_new.z = 24 | length_unit
        particles_new.vx = 25 | speed_unit
        particles_new.vy = 26 | speed_unit
        particles_new.vz = 27 | speed_unit
        
        instance.particles.add_particles(particles_new)
        self.assertEqual(len(instance.particles), 9)
        
        self.assertAlmostRelativeEquals(instance.particles.mass.as_quantity_in(mass_unit), [1,2,4,5,7,8,9,10,20] | mass_unit)
        self.assertAlmostRelativeEquals(instance.particles.x.as_quantity_in(length_unit), [1,2,4,5,7,8,9,10,22] | length_unit)
        
        instance.cleanup_code()
        instance.stop()
        
    def test1(self):
        
        instance = Hermite()
        self._run_addition_removal_test(instance)
        
    def test2(self):
        
        instance = PhiGRAPE()
        self._run_addition_removal_test(instance)
        
    def test3(self):
        
        instance = BHTree()
        self._run_addition_removal_test(instance)
        
    def test4(self):
        
        instance = ph4()
        self._run_addition_removal_test(instance)
        
    def test5(self):
        instance = Fi()
        self._run_addition_removal_test(instance)
        
    def test6(self):
        length_unit = units.parsec
        speed_unit = units.parsec / units.Myr
        mass_unit = units.MSun
        instance = Gadget2()
        self._run_addition_removal_test(instance, length_unit, speed_unit, mass_unit)
