from amuse.test import amusetest
import os.path
import numpy
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.units import quantities
from amuse.ic.plummer import new_plummer_model
from amuse.datamodel import rotation
from amuse.datamodel import particles

class TestRotations(amusetest.TestCase):

    def test01(self):
        x = rotation.new_rotation_matrix(0,0,0)
        self.assertAlmostRelativeEquals(x, ((1,0,0),(0,1,0),(0,0,1)))

        x = rotation.new_rotation_matrix(numpy.pi,0,0)
        print(x)
        self.assertAlmostRelativeEquals(x, ((1,0,0),(0,-1,0),(0,0,-1)))
        x = rotation.new_rotation_matrix(0,numpy.pi,0)
        self.assertAlmostRelativeEquals(x, ((-1,0,0),(0,1,0),(0,0,-1)))
        x = rotation.new_rotation_matrix(0, 0, numpy.pi)
        self.assertAlmostRelativeEquals(x, ((-1,0,0),(0,-1,0),(0,0,1)))
    def test02(self):
        x = rotation.new_rotation_matrix(numpy.pi / 2.0,0,0)
        print(x)
        self.assertAlmostRelativeEquals(x, ((1,0,0),(0,0, -1),(0,1,0)))
        x = rotation.new_rotation_matrix(0,numpy.pi / 2.0,0)
        print(x)
        self.assertAlmostRelativeEquals(x, ((0,0,1),(0,1,0),(-1,0,0)))
        x = rotation.new_rotation_matrix(0, 0, numpy.pi / 2.0)
        print(x)
        self.assertAlmostRelativeEquals(x, ((0,-1,0),(1,0,0),(0,0,1)))
        x = rotation.new_rotation_matrix(numpy.pi / 2.0,numpy.pi / 2.0,0)
        print(x)
        self.assertAlmostRelativeEquals(x, ((0,1,0),(0,0,-1),(-1,0,0)))

    def test03(self):
        positions = [ [1.0, 2.0, 3.0 ] ] | units.m
        rotated = rotation.rotated(positions, 0.0, 0.0, 0.0)
    
        print(rotated)
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, 2.0, 3.0 ] ] | units.m)
        rotated = rotation.rotated(positions, numpy.pi, 0.0, 0.0)
    
        print(rotated)
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, -2.0, -3.0 ] ] | units.m)

    def test04(self):
        positions = [ [1.0, 2.0, 3.0 ] , [4.0, 5.0, 6.0] ] | units.m
        rotated = rotation.rotated(positions, 0.0, 0.0, 0.0)
    
        print(rotated)
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, 2.0, 3.0 ], [4.0, 5.0, 6.0] ] | units.m)
    
        rotated = rotation.rotated(positions, numpy.pi, 0.0, 0.0)
        print(rotated)
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, -2.0, -3.0 ], [4.0, -5.0, -6.0] ] | units.m)

    def test05(self):
        positions = [ [1.0, 2.0, 3.0 ] , [4.0, 5.0, 6.0] ] | units.m
        rotated = rotation.rotated(positions,  numpy.pi/2, 0.0, 0.0)
    
        print(rotated)
    
        self.assertAlmostRelativeEquals(rotated, [ [1.0, -3.0, 2.0 ], [4.0, -6.0, 5.0] ] | units.m)

    def test06(self):
        particles = new_plummer_model(100)
        kinetic_energy0 = particles.kinetic_energy()
        potential_energy0 = particles.potential_energy(G=nbody_system.G)
        particles.position = rotation.rotated(particles.position,  numpy.pi/2, 0.0, 0.0)
        particles.velocity = rotation.rotated(particles.velocity,  numpy.pi/2, 0.0, 0.0)
    
    
        kinetic_energy1 = particles.kinetic_energy()
        potential_energy1 = particles.potential_energy(G=nbody_system.G)
        self.assertAlmostRelativeEquals(kinetic_energy1, kinetic_energy0)
        self.assertAlmostRelativeEquals(potential_energy1, potential_energy0)

    def test07(self):
        particles = new_plummer_model(100)
        kinetic_energy0 = particles.kinetic_energy()
        potential_energy0 = particles.potential_energy(G=nbody_system.G)
        particles.position = rotation.rotated(particles.position,  numpy.pi/3, numpy.pi/2, 0.0)
        particles.velocity = rotation.rotated(particles.velocity,  numpy.pi/3, numpy.pi/2, 0.0)
    
    
        kinetic_energy1 = particles.kinetic_energy()
        potential_energy1 = particles.potential_energy(G=nbody_system.G)
        self.assertAlmostRelativeEquals(kinetic_energy1, kinetic_energy0)
        self.assertAlmostRelativeEquals(potential_energy1, potential_energy0)
    
    def test08(self):
        particles = new_plummer_model(100)
        kinetic_energy0 = particles.kinetic_energy()
        potential_energy0 = particles.potential_energy(G=nbody_system.G)
        
        particles.move_to_center()
        particles.position += [3, 0, 2] | nbody_system.length
        particles.rotate(numpy.pi/4, numpy.pi/2, 0.0)
        self.assertAlmostRelativeEquals(particles.center_of_mass(), 
            [numpy.sqrt(2), -numpy.sqrt(2), -3] | nbody_system.length, 7)
        
        kinetic_energy1 = particles.kinetic_energy()
        potential_energy1 = particles.potential_energy(G=nbody_system.G)
        self.assertAlmostRelativeEquals(kinetic_energy1, kinetic_energy0)
        self.assertAlmostRelativeEquals(potential_energy1, potential_energy0)
    
    def test09(self):
        print("Test add_spin particle attribute, to add rigid body rotation")
        numpy.random.seed(123456)
        particles = new_plummer_model(1000)
        kinetic_energy0 = particles.kinetic_energy()
        potential_energy0 = particles.potential_energy(G=nbody_system.G)
        
        particles.position += [3, 0, 2] | nbody_system.length
        particles.velocity += [1, 10, 100] | nbody_system.speed
        particles.add_spin(3.0 | nbody_system.time**-1)
        
        self.assertAlmostRelativeEquals(particles.center_of_mass(), [3, 0, 2] | nbody_system.length, 12)
        self.assertAlmostRelativeEquals(particles.potential_energy(G=nbody_system.G), potential_energy0, 12)
        self.assertAlmostRelativeEquals(particles.center_of_mass_velocity(), [1, 10, 100] | nbody_system.speed, 12)
        
        r = particles.position - particles.center_of_mass()
        v = particles.velocity - particles.center_of_mass_velocity()
        spin_direction = (r).cross(v).mean(axis=0)
        spin_direction /= spin_direction.length()
        R = r - r*spin_direction
        omega = ((R).cross(v) / R.lengths_squared().reshape((-1,1))).mean(axis=0).length()
        self.assertAlmostEqual(spin_direction, [0, 0, 1.0], 1)
        self.assertAlmostEqual(omega, 3.0 | nbody_system.time**-1, 1)
        
        particles.add_spin([1.0, 0, -3.0] | nbody_system.time**-1)
        v = particles.velocity - particles.center_of_mass_velocity()
        spin_direction = (r).cross(v).mean(axis=0)
        spin_direction /= spin_direction.length()
        R = r - r*spin_direction
        omega = ((R).cross(v) / R.lengths_squared().reshape((-1,1))).mean(axis=0).length()
        self.assertAlmostEqual(omega, 1.0 | nbody_system.time**-1, 1)
 
    def test10(self):
        print("test conservation of dot, transformation of cross")
        p=particles.Particles(1)
        p.position=[1.,2.,3.]
        p.velocity=[-4,5,6.]
        
        dot1=p[0].position.dot(p[0].velocity)
        cross1=numpy.cross(p[0].position,p[0].velocity)
        
        rm=rotation.new_rotation_matrix(0.1,0.5,3.5)
        p.rotate(0.1,0.5,3.5)
        
        dot2=p[0].position.dot(p[0].velocity)
        cross2=numpy.cross(p[0].position,p[0].velocity)
        
        self.assertAlmostRelativeEquals(dot1,dot2)
        self.assertAlmostRelativeEquals(cross2,cross1.dot(numpy.linalg.inv(rm)))

    def test11(self):
        print("test conservation of dot, transformation of cross with units")
        p=particles.Particles(5)
        p.position=[1.,2.,3.] | units.km
        p.velocity=[-4,5,6.] | units.kms
        
        p[1:].position*=0
        p[1:].velocity*=0
                
        dot1=p[0].position.dot(p[0].velocity)
        cross1=p[0].position.cross(p[0].velocity)
        
        rm=rotation.new_rotation_matrix(0.1,0.5,3.5)
        p.rotate(0.1,0.5,3.5)
        
        dot2=p[0].position.dot(p[0].velocity)
        cross2=p[0].position.cross(p[0].velocity)
                
        self.assertAlmostRelativeEquals(dot1.value_in(units.km**2/units.s),dot2.value_in(units.km**2/units.s))
        self.assertAlmostRelativeEquals(cross2.value_in(units.km**2/units.s),cross1.dot(numpy.linalg.inv(rm)).value_in(units.km**2/units.s))
        
        
        
        
      
      

