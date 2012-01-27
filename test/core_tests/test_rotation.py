from amuse.test import amusetest
import os.path
import numpy
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.units import quantities
from amuse.ic.plummer import new_plummer_sphere
from amuse.datamodel import rotation
from amuse.datamodel import particles

class TestRotations(amusetest.TestCase):

    def test01(self):
        x = rotation.new_rotation_matrix(0,0,0)
        self.assertAlmostRelativeEquals(x, ((1,0,0),(0,1,0),(0,0,1)))

        x = rotation.new_rotation_matrix(numpy.pi,0,0)
        print x
        self.assertAlmostRelativeEquals(x, ((1,0,0),(0,-1,0),(0,0,-1)))
        x = rotation.new_rotation_matrix(0,numpy.pi,0)
        self.assertAlmostRelativeEquals(x, ((-1,0,0),(0,1,0),(0,0,-1)))
        x = rotation.new_rotation_matrix(0, 0, numpy.pi)
        self.assertAlmostRelativeEquals(x, ((-1,0,0),(0,-1,0),(0,0,1)))
    def test02(self):
        x = rotation.new_rotation_matrix(numpy.pi / 2.0,0,0)
        print x
        self.assertAlmostRelativeEquals(x, ((1,0,0),(0,0, -1),(0,1,0)))
        x = rotation.new_rotation_matrix(0,numpy.pi / 2.0,0)
        print x
        self.assertAlmostRelativeEquals(x, ((0,0,1),(0,1,0),(-1,0,0)))
        x = rotation.new_rotation_matrix(0, 0, numpy.pi / 2.0)
        print x
        self.assertAlmostRelativeEquals(x, ((0,-1,0),(1,0,0),(0,0,1)))
        x = rotation.new_rotation_matrix(numpy.pi / 2.0,numpy.pi / 2.0,0)
        print x
        self.assertAlmostRelativeEquals(x, ((0,1,0),(0,0,-1),(-1,0,0)))

    def test03(self):
        positions = [ [1.0, 2.0, 3.0 ] ] | units.m
        rotated = rotation.rotated(positions, 0.0, 0.0, 0.0)
    
        print rotated
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, 2.0, 3.0 ] ] | units.m)
        rotated = rotation.rotated(positions, numpy.pi, 0.0, 0.0)
    
        print rotated
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, -2.0, -3.0 ] ] | units.m)

    def test04(self):
        positions = [ [1.0, 2.0, 3.0 ] , [4.0, 5.0, 6.0] ] | units.m
        rotated = rotation.rotated(positions, 0.0, 0.0, 0.0)
    
        print rotated
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, 2.0, 3.0 ], [4.0, 5.0, 6.0] ] | units.m)
    
        rotated = rotation.rotated(positions, numpy.pi, 0.0, 0.0)
        print rotated
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, -2.0, -3.0 ], [4.0, -5.0, -6.0] ] | units.m)

    def test05(self):
        positions = [ [1.0, 2.0, 3.0 ] , [4.0, 5.0, 6.0] ] | units.m
        rotated = rotation.rotated(positions,  numpy.pi/2, 0.0, 0.0)
    
        print rotated
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, -3.0, 2.0 ], [4.0, -6.0, 5.0] ] | units.m)

    def test06(self):
        particles = new_plummer_sphere(100)
        kinetic_energy0 = particles.kinetic_energy()
        potential_energy0 = particles.potential_energy()
        particles.position = rotation.rotated(particles.position,  numpy.pi/2, 0.0, 0.0)
        particles.velocity = rotation.rotated(particles.velocity,  numpy.pi/2, 0.0, 0.0)
    
    
        kinetic_energy1 = particles.kinetic_energy()
        potential_energy1 = particles.potential_energy()
        self.assertAlmostRelativeEquals(kinetic_energy1, kinetic_energy0)
        self.assertAlmostRelativeEquals(potential_energy1, potential_energy0)

    def test07(self):
        particles = new_plummer_sphere(100)
        kinetic_energy0 = particles.kinetic_energy()
        potential_energy0 = particles.potential_energy()
        particles.position = rotation.rotated(particles.position,  numpy.pi/3, numpy.pi/2, 0.0)
        particles.velocity = rotation.rotated(particles.velocity,  numpy.pi/3, numpy.pi/2, 0.0)
    
    
        kinetic_energy1 = particles.kinetic_energy()
        potential_energy1 = particles.potential_energy()
        self.assertAlmostRelativeEquals(kinetic_energy1, kinetic_energy0)
        self.assertAlmostRelativeEquals(potential_energy1, potential_energy0)
    
    def test08(self):
        particles = new_plummer_sphere(100)
        kinetic_energy0 = particles.kinetic_energy()
        potential_energy0 = particles.potential_energy()
        
        particles.move_to_center()
        particles.position += [3, 0, 2] | nbody_system.length
        particles.rotate(numpy.pi/4, numpy.pi/2, 0.0)
        self.assertAlmostRelativeEquals(particles.center_of_mass(), 
            [numpy.sqrt(2), -numpy.sqrt(2), -3] | nbody_system.length)
        
        kinetic_energy1 = particles.kinetic_energy()
        potential_energy1 = particles.potential_energy()
        self.assertAlmostRelativeEquals(kinetic_energy1, kinetic_energy0)
        self.assertAlmostRelativeEquals(potential_energy1, potential_energy0)

