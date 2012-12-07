from amuse.test import amusetest

from amuse.support.exceptions import AmuseException
from amuse import datamodel
from amuse.units import units
from amuse.units import nbody_system

import numpy

class TestParticleLinkToParticle(amusetest.TestCase):
    """
    Tests One-to-One relation between particles
    """
    def test1(self):
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        
        parent = particles[0]
        child1 = particles[1]
        child2 = particles[2]
        
        parent.child1 = child1
        parent.child2 = child2
        child1.sibling = child2
        
        self.assertEquals(len(particles.child1), 3)
        self.assertEquals(len(particles.child1.as_set().compressed()), 1)
        self.assertEquals(len(particles.child2), 3)
        self.assertEquals(len(particles.child2.as_set().compressed()), 1)
        
        self.assertAlmostRelativeEqual(parent.child1.mass, 3 | units.kg)

    def test2(self):
        
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        
        parent = particles[0]
        child1 = particles[1]
        child2 = particles[2]
        
        parent.child1 = child1
        parent.child2 = child2
        child1.sibling = child2
        
        
        convert_nbody = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        particles = datamodel.ParticlesWithUnitsConverted(
            particles, 
            convert_nbody.as_converter_from_generic_to_si()
        )
        self.assertEquals(len(particles.child1), 3)
        self.assertEquals(len(particles.child1.as_set().compressed()), 1)
        self.assertEquals(len(particles.child2), 3)
        self.assertEquals(len(particles.child2.as_set().compressed()), 1)
        
        self.assertAlmostRelativeEqual(particles[0].child1.mass, 0.3 | nbody_system.mass)
        self.assertAlmostRelativeEqual(particles[0].child1.mass, 0.3 | nbody_system.mass)
    
    def test3(self):
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        
        parent = particles[0]
        child1 = particles[1]
        child2 = particles[2]
        
        parent.child1 = child1
        parent.child2 = child2
        child1.sibling = child2
        
        other = datamodel.Particles()
        particles.synchronize_to(other)
        root = None
        for x in other:
            if not x.child1 is None:
               root = x
               break
        
        self.assertFalse(root is None)
        
        self.assertAlmostRelativeEqual(root.child1.mass, 3 | units.kg)
        self.assertTrue(root.child1.as_set()._original_set() is other)
        self.assertTrue(root.child1.sibling.as_set()._original_set() is other)
        
    def test4(self):
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        other = particles.copy()
        
        parent = particles[0]
        child1 = particles[1]
        child2 = particles[2]
        
        parent.child1 = child1
        parent.child2 = child2
        child1.sibling = child2
        
        channel = particles.new_channel_to(other)
        channel.copy()
        
        self.assertTrue(other[0].child1.as_set()._original_set() is other)
        self.assertTrue(other[1].sibling.as_set()._original_set() is other)
        self.assertAlmostRelativeEqual(other[0].child1.mass, 3 | units.kg)
    
    def test5(self):
        
        particles = datamodel.Particles(3)
        parent = particles[0]
        child1 = particles[1]
        child2 = particles[2]
        
        particles.mass = [2,3,4] | units.kg
        
        parent.child1 = child1
        parent.child2 = child2
        child1.sibling = child2
        
        other = particles.copy()
        
        self.assertAlmostRelativeEqual(other[0].child1.mass, 3 | units.kg)
        self.assertTrue(other[0].child1.as_set()._original_set() is other)
        self.assertTrue(other[1].sibling.as_set()._original_set() is other)  
          
    def test6(self):
        
        binaries = datamodel.Particles(2)
        binaries.mass = [4,6] | units.kg
        
        stars = datamodel.Particles(4)
        stars.mass = [1,3,4,2] | units.kg
        
        binaries[0].child1 = stars[0]
        binaries[0].child2 = stars[1]
        binaries[1].child1 = stars[3]
        binaries[1].child2 = stars[2]
        
        other = binaries.copy()
        stars.mass = 2 * ([1,3,4,2] | units.kg)
        self.assertAlmostRelativeEqual(binaries[0].child1.mass, 2 | units.kg)
        self.assertAlmostRelativeEqual(other[0].child1.mass, 1 | units.kg)
        self.assertAlmostRelativeEqual(other[0].child2.mass, 3 | units.kg)
        self.assertAlmostRelativeEqual(binaries[1].child1.mass, 4 | units.kg)
        self.assertAlmostRelativeEqual(other[1].child1.mass, 2 | units.kg)
        self.assertAlmostRelativeEqual(other[1].child2.mass, 4 | units.kg)
    
          
    def test7(self):
        
        binaries = datamodel.Particles(2)
        binaries.mass = [4,6] | units.kg
        
        stars = datamodel.Particles(4)
        stars.mass = [1,3,4,2] | units.kg
        
        binaries[0].child1 = stars[0]
        binaries[0].child2 = stars[1]
        binaries[1].child1 = stars[3]
        binaries[1].child2 = stars[2]
        stars[0].binary = binaries[0]
        stars[1].binary = binaries[0]
        stars[2].binary = binaries[1]
        stars[3].binary = binaries[1]
        
        other = binaries.copy()
        stars.mass *= 2
        binaries.mass  *= 2
        
        self.assertAlmostRelativeEqual(other[0].child1.mass, 1 | units.kg)
        self.assertAlmostRelativeEqual(binaries[0].child1.mass, 2 | units.kg)
        self.assertAlmostRelativeEqual(other[0].child2.mass, 3 | units.kg)
        self.assertAlmostRelativeEqual(other[1].child1.mass, 2 | units.kg)
        self.assertAlmostRelativeEqual(other[1].child2.mass, 4 | units.kg)
        self.assertAlmostRelativeEqual(binaries[0].child1.binary.mass, 8 | units.kg)
        self.assertAlmostRelativeEqual(other[0].child1.binary.mass, 4 | units.kg)
        self.assertAlmostRelativeEqual(other[1].child1.binary.mass, 6 | units.kg)
        
    def test8(self):
        
        particles = datamodel.Particles(2)
        particles.mass = [2,3] | units.kg
        
        binary = datamodel.Particles(1)
        binary[0].mass = 5 | units.kg
        binary[0].child1 = particles[0]
        binary[0].child2 = particles[1]
        
        binaries = datamodel.Particles()
        binaries.add_particles(binary)
        
        self.assertAlmostRelativeEqual(binaries[0].child1.mass, 2 | units.kg)
        self.assertAlmostRelativeEqual(binaries[0].child2.mass, 3 | units.kg)
        
    
    
    def test9(self):
        particles1 = datamodel.Particles(keys=[9,10,11,12])
        particles1.mass = [1,2,3,4] | units.kg
        
        particles2 = datamodel.Particles(keys=[13,14,15])
        particles2.mass = [20,21,22] | units.kg
        particles1[0].particle2 = particles2[0]
        particles1[1].particle2 = particles2[0]
        
        particles1_copy = particles1.copy()
        
        self.assertEquals(particles1_copy[0].particle2.key,  13)
        self.assertEquals(particles1_copy[1].particle2.key,  13)
        
        particles1_copy[0].particle2.mass = 30 | units.kg
        
        self.assertAlmostRelativeEquals(particles1_copy[0].particle2.mass, 30 | units.kg)
        self.assertAlmostRelativeEquals(particles1_copy[1].particle2.mass, 30 | units.kg)
        

class TestParticleLinkToParticles(amusetest.TestCase):
    """
    Tests One-to-Many relation between a particle and particle-sets
    """
    def test0(self):
        particles1 = datamodel.Particles(keys=[9,10,11,12])
        particles1.mass = [1,2,3,4] | units.kg
        x = numpy.empty((10,), dtype=numpy.object)
        index = numpy.asarray([1])
        x[index] = particles1
        print x
        self.assertEquals(len(x[1]), 4)
        
    def test1(self):
        particles1 = datamodel.Particles(keys=[9,10,11,12])
        particles1.mass = [1,2,3,4] | units.kg
        particles2 = datamodel.Particles(keys=[13,14])
        particles2.mass = [20,21] | units.kg
        particles1[1].particles2 = particles2
        self.assertEquals(particles1[1].particles2.key,  [13,14])
        self.assertAlmostRelativeEquals(particles1[1].particles2.mass, [20,21] | units.kg)
        self.assertEquals(particles1[0].particles2, None)
        
    def test2(self):
        particles1 = datamodel.Particles(keys=[9,10,11,12])
        particles1.mass = [1,2,3,4] | units.kg
        particles2 = datamodel.Particles(keys=[13,14])
        particles2.mass = [20,21] | units.kg
        particles1[1].particles2 = particles2
        self.assertEquals(len(particles1.particles2), 4)
        self.assertEquals(particles1.particles2[1].key,  [13,14])
        self.assertEquals(particles1.particles2[0],  None)
        
    def test3(self):
        particles1 = datamodel.Particles(keys=[9,10,11,12])
        particles1.mass = [1,2,3,4] | units.kg
        particles2 = datamodel.Particles(keys=[13,14])
        particles2.mass = [20,21] | units.kg
        particles1[1].particles2 = particles2
        
        particles1_copy = particles1.copy()
        self.assertEquals(particles1_copy[1].particles2.key,  [13,14])
        self.assertAlmostRelativeEquals(particles1_copy[1].particles2.mass, [20,21] | units.kg)
        self.assertEquals(particles1_copy[0].particles2, None)
        self.assertNotEquals(id(particles1_copy[1].particles2), id(particles1[1].particles2))
        
    
    def test4(self):
        particles1 = datamodel.Particles(keys=[9,10,11,12])
        particles1.mass = [1,2,3,4] | units.kg
        particles1_copy = particles1.copy()
        
        particles2 = datamodel.Particles(keys=[13,14])
        particles2.mass = [20,21] | units.kg
        particles1[1].particles2 = particles2
        
        channel = particles1.new_channel_to(particles1_copy)
        channel.copy()
        
        self.assertEquals(particles1_copy[1].particles2.key,  [13,14])
        self.assertAlmostRelativeEquals(particles1_copy[1].particles2.mass, [20,21] | units.kg)
        self.assertEquals(particles1_copy[0].particles2, None)
        self.assertEquals(id(particles1_copy[1].particles2), id(particles1[1].particles2))
        
    def test5(self):
        particles1 = datamodel.Particles(keys=[9,10,11,12])
        particles1.mass = [1,2,3,4] | units.kg
        
        particles2 = datamodel.Particles(keys=[13,14])
        particles2.mass = [20,21] | units.kg
        particles1.particles2 = [None, particles2, None, particles2[0]]
        
        self.assertEquals(particles1[1].particles2.key,  [13,14])
        self.assertEquals(particles1[1].particles2.mass, [20,21] | units.kg)
        self.assertAlmostRelativeEquals(particles1[3].particles2.key,  13)
        self.assertAlmostRelativeEquals(particles1[3].particles2.mass, 20 | units.kg)
        self.assertEquals(particles1[0].particles2, None)
        self.assertEquals(particles1[2].particles2, None)
        
    
    def test6(self):
        particles1 = datamodel.Particles(keys=[9,10,11,12])
        particles1.mass = [1,2,3,4] | units.kg
        
        particles2 = datamodel.Particles(keys=[13,14,15])
        particles2.mass = [20,21,22] | units.kg
        particles1[0].particles2 = particles2[0:2]
        particles1[1].particles2 = particles2[1:]
        
        particles1_copy = particles1.copy(keep_structure=False)
        
        self.assertEquals(particles1_copy[0].particles2.key,  [13,14])
        self.assertEquals(particles1_copy[1].particles2.key,  [14,15])
        
        particles1_copy[0].particles2[1].mass = 30 | units.kg
        
        self.assertAlmostRelativeEquals(particles1_copy[0].particles2[1].mass, 30 | units.kg)
        self.assertAlmostRelativeEquals(particles1_copy[1].particles2[0].mass, 21 | units.kg)
        
        particles1_copy = particles1.copy(keep_structure=True)
        
        self.assertEquals(particles1_copy[0].particles2.key,  [13,14])
        self.assertEquals(particles1_copy[1].particles2.key,  [14,15])
        
        particles1_copy[0].particles2[1].mass = 30 | units.kg
        
        self.assertAlmostRelativeEquals(particles1_copy[0].particles2[1].mass, 30 | units.kg)
        self.assertAlmostRelativeEquals(particles1_copy[1].particles2[0].mass, 30 | units.kg)
        

class TestParticleLinkToGridPoint(amusetest.TestCase):
    """
    Tests One-to-One relation between particles and gridpoints
    """
    def test1(self):
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        particles[0].point11 = grid[1][1]
        
        self.assertAlmostRelativeEquals(grid[1][1].rho, 6 |  units.kg / units.m**3)
        self.assertAlmostRelativeEquals(particles[0].point11.rho, 6 |  units.kg / units.m**3)
        
    def test2(self):
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        particles[0].point11 = grid[1][1]
        
        particles_copy = particles.copy()
        
        self.assertAlmostRelativeEquals(grid[1][1].rho, 6 |  units.kg / units.m**3)
        self.assertAlmostRelativeEquals(particles_copy[0].point11.rho, 6 |  units.kg / units.m**3)
        
    def test3(self):
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        particles_copy = particles.copy()
        
        particles[0].point11 = grid[1][1]
        
        channel = particles.new_channel_to(particles_copy)
        channel.copy()
        
        self.assertAlmostRelativeEquals(grid[1][1].rho, 6 |  units.kg / units.m**3)
        self.assertAlmostRelativeEquals(particles_copy[0].point11.rho, 6 |  units.kg / units.m**3)


class TestParticleLinkToGrids(amusetest.TestCase):
    """
    Tests One-to-Many relation between particles and grids
    """
    
    def test1(self):
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        particles[0].grid = grid
        
        self.assertAlmostRelativeEquals(grid[1][1].rho, 6 |  units.kg / units.m**3)
        self.assertAlmostRelativeEquals(particles[0].grid[1][1].rho, 6 |  units.kg / units.m**3)
        
    def test2(self):
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        particles[0].grid = grid
        
        self.assertAlmostRelativeEquals(grid[1][1].rho, 6 |  units.kg / units.m**3)
        self.assertAlmostRelativeEquals(particles[0].grid[1][1].rho, 6 |  units.kg / units.m**3)
        
        particles_copy = particles.copy()
        
        grid[1][1].rho = 10 |  units.kg / units.m**3
        
        self.assertAlmostRelativeEquals(grid[1][1].rho, 10 |  units.kg / units.m**3)
        self.assertAlmostRelativeEquals(particles[0].grid[1][1].rho, 10 |  units.kg / units.m**3)
        self.assertAlmostRelativeEquals(particles_copy[0].grid[1][1].rho, 6 |  units.kg / units.m**3)
    
    
        
    def test3(self):
        
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        particles_copy = particles.copy()
        
        particles[0].grid = grid
        
        self.assertAlmostRelativeEquals(grid[1][1].rho, 6 |  units.kg / units.m**3)
        self.assertAlmostRelativeEquals(particles[0].grid[1][1].rho, 6 |  units.kg / units.m**3)
        
        channel = particles.new_channel_to(particles_copy)
        channel.copy()
        
        self.assertAlmostRelativeEquals(particles_copy[0].grid[1][1].rho, 6 |  units.kg / units.m**3)

class TestGridPointLinkToParticle(amusetest.TestCase):
    """
    Tests One-to-One relation between gridpoints and particles
    """
    def test1(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        
        grid[0][0].particle = particles[1]
        
        self.assertAlmostRelativeEquals(grid[0][0].particle.mass, 3 | units.kg)
        self.assertEquals(grid[0][0].particle, particles[1])
        self.assertEquals(grid[1][1].particle, None)

    def test2(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        
        grid[0][0].particle = particles[1]
        
        grid_copy = grid.copy()
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].particle.mass, 3 | units.kg)
        self.assertEquals(grid_copy[0][0].particle, particles[1])
        self.assertEquals(grid_copy[1][1].particle, None)

        grid[0][0].particle.mass = 5 | units.kg
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].particle.mass, 3 | units.kg)

    def test3(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        
    
        grid_copy = grid.copy()
    
        grid[0][0].particle = particles[1]
        
        channel = grid.new_channel_to(grid_copy)
        channel.copy()
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].particle.mass, 3 | units.kg)
        self.assertEquals(grid_copy[0][0].particle, particles[1])
        self.assertEquals(grid_copy[1][1].particle, None)

class TestGridPointLinkToParticles(amusetest.TestCase):
    """
    Tests One-to-Many relation between gridpoints and particles
    """
    def test1(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        
        grid[0][0].particles = particles
        
        self.assertAlmostRelativeEquals(grid[0][0].particles[1].mass, 3 | units.kg)
        self.assertEquals(grid[0][0].particles[1], particles[1])
        self.assertEquals(grid[1][1].particles, None)
        
    
    def test2(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        
        grid[0][0].particles = particles
        
        grid_copy = grid.copy()
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].particles[1].mass, 3 | units.kg)
        
        grid[0][0].particles[1].mass = 10 | units.kg
        
        
    
    def test3(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        particles = datamodel.Particles(3)
        particles.mass = [2,3,4] | units.kg
        
        
        grid_copy = grid.copy()
        grid[0][0].particles = particles
        
        channel = grid.new_channel_to(grid_copy)
        channel.copy()
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].particles[1].mass, 3 | units.kg)
        
        grid[0][0].particles[1].mass = 10 | units.kg
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].particles[1].mass, 10 | units.kg)
        

class TestGridPointLinkToGridPoint(amusetest.TestCase):
    """
    Tests One-to-One relation between gridpoints
    """
    def test1(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        grid[0][0].neighbour = grid[0][1]
        
        self.assertAlmostRelativeEquals(grid[0][0].neighbour.rho, 3  | units.kg / units.m**3)
        self.assertEquals(grid[0][0].neighbour, grid[0][1])
        self.assertEquals(grid[1][1].neighbour, None)

    def test2(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        grid[0][0].neighbour = grid[0][1]
        
        grid_copy = grid.copy()
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].neighbour.rho, 3  | units.kg / units.m**3)
        self.assertEquals(grid_copy[0][0].neighbour, grid_copy[0][1])
        self.assertEquals(grid_copy[1][1].neighbour, None)

        grid[0][0].neighbour.rho = 5  | units.kg / units.m**3
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].neighbour.rho, 3  | units.kg / units.m**3)
        
        grid_copy[0][1].rho = 6 | units.kg / units.m**3
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].neighbour.rho, 6  | units.kg / units.m**3)

    

    def test3(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        
        grid_copy = grid.copy()
        
        grid[0][0].neighbour = grid[0][1]
        
        channel = grid.new_channel_to(grid_copy)
        channel.copy()
        
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].neighbour.rho, 3  | units.kg / units.m**3)
        self.assertEquals(grid_copy[0][0].neighbour, grid_copy[0][1])
        self.assertEquals(grid_copy[1][1].neighbour, None)

        grid[0][0].neighbour.rho = 5  | units.kg / units.m**3
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].neighbour.rho, 3  | units.kg / units.m**3)
        
        grid_copy[0][1].rho = 6 | units.kg / units.m**3
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].neighbour.rho, 6  | units.kg / units.m**3)
        

class TestGridPointLinkToGrid(amusetest.TestCase):
    """
    Tests One-to-Many relation between gridpoints
    """
    def test1(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        grid[0][0].container = grid
        
        self.assertAlmostRelativeEquals(grid[0][0].container[0][1].rho, 3  | units.kg / units.m**3)
        self.assertEquals(grid[0][0].container[0][1], grid[0][1])
        self.assertEquals(grid[1][1].container, None)

    def test2(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        grid[0][0].container = grid
        
        grid_copy = grid.copy()
        
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].container[0][1].rho, 3  | units.kg / units.m**3)
        self.assertEquals(grid_copy[0][0].container[0][1], grid_copy[0][1])
        self.assertEquals(grid_copy[1][1].container, None)


        grid[0][0].container[0][1].rho = 5  | units.kg / units.m**3
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].container[0][1].rho, 3  | units.kg / units.m**3)
        
        grid_copy[0][1].rho = 6 | units.kg / units.m**3
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].container[0][1].rho, 6  | units.kg / units.m**3)

    

    def test3(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        
        grid_copy = grid.copy()
        
        grid[0][0].container = grid
        
        channel = grid.new_channel_to(grid_copy)
        channel.copy()
        
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].container[0][1].rho, 3  | units.kg / units.m**3)
        self.assertEquals(grid_copy[0][0].container[0][1], grid_copy[0][1])
        self.assertEquals(grid_copy[1][1].container, None)

        grid[0][0].container[0][1].rho = 5  | units.kg / units.m**3
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].container[0][1].rho, 3  | units.kg / units.m**3)
        
        grid_copy[0][1].rho = 6 | units.kg / units.m**3
        
        self.assertAlmostRelativeEquals(grid_copy[0][0].container[0][1].rho, 6  | units.kg / units.m**3)
        
    
    def test4(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        grid[...,0].container = grid

        self.assertEquals(id(grid[0][0].container), id(grid))
        self.assertEquals(id(grid[1][0].container), id(grid)) 
        self.assertEquals(grid[0][1].container, None)
        self.assertEquals(grid[1][1].container, None)
        
    
    def test5(self):
        
        grid = datamodel.Grid(2,3)
        grid.rho = [[2,3,4],[5,6,7]] | units.kg / units.m**3
        
        grid.container = grid

        for index in numpy.ndindex(*grid.shape):
            self.assertEquals(id(grid[index].container), id(grid))
            self.assertEquals(grid[index].rho, (index[0] * 3 + index[1] + 2) | units.kg / units.m**3) 