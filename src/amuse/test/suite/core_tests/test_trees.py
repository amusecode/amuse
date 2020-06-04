from amuse.test import amusetest

from amuse.support.exceptions import AmuseException
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.support.interface import InCodeComponentImplementation
import numpy
import time
import sys

from amuse.datamodel import trees
from amuse.datamodel import Particles
class TestBinaryTree(amusetest.TestCase):
    
    def test1(self):
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        
        self.assertEqual(particles[0].mass, 0 | units.kg)
        self.assertEqual(particles[0].child1.mass, 1 | units.kg)
        self.assertEqual(particles[1].child1, None)
    
        children1 = particles.child1.as_set().compressed()
        children2 = particles.child2.as_set().compressed()
        children = children1 + children2
        roots = particles - children
    
        print(len(roots))
        self.assertEqual(len(roots), 8)
        self.assertEqual(len(children), 2)
        
    def test2(self):
        n = 100000
        particles = Particles(n)
        particles.mass = list(range(n)) | units.kg
        particles[n-1].child1 = particles[0]
        particles[n-1].child2 = particles[1]
        
        self.assertEqual(particles[0].mass, 0 | units.kg)
        self.assertEqual(particles[n-1].child1.mass, 0 | units.kg)
        self.assertEqual(particles[n-1].child2.mass, 1 | units.kg)
    
        children1 = particles.child1.as_set().compressed()
        children2 = particles.child2.as_set().compressed()
        children = children1 + children2
        roots = particles - children
    
        self.assertEqual(len(roots), n - 2)
        self.assertEqual(len(children), 2)
        binaries = particles.select_array(lambda x : x != [None], ['child1',])
        self.assertEqual(len(binaries), 1)
    def test3(self):
        n = 10
        particles = Particles(n)
        particles.mass = list(range(n)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        particles[1].child1 = particles[3]
        particles[1].child2 = particles[4]
        
        self.assertEqual(particles[0].child1.child1.mass, 3 | units.kg)
    
        binaries = particles.select_array(lambda x : x != [None], ["child1",])
    
        print(len(binaries))
        self.assertEqual(len(binaries), 2)
    
        binaries_children1 = binaries.child1.as_set().compressed().select_array(lambda x : x != [None], ["child1",])
        binaries_children2 = binaries.child2.as_set().compressed().select_array(lambda x : x != [None], ["child1",])
        binaries_roots = binaries - (binaries_children1 + binaries_children2)
    
        self.assertEqual(len(binaries_roots), 1)

        self.assertEqual(binaries_roots[0].mass, 0 | units.kg)


    def test4(self):
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        particles[1].child1 = particles[3]
        particles[1].child2 = particles[4]
        
    
        x = trees.BinaryTreesOnAParticleSet(particles, "child1", "child2")
        roots = list(x.iter_roots())
    
        self.assertEqual(len(roots), 1)
        x = [x.mass.value_in(units.kg) for x in roots[0].iter_descendants()]
        self.assertEqual(x, [1,2,3,4])
        self.assertEqual(roots[0].get_descendants_subset().mass, [1,2,3,4] | units.kg)

    def test5(self):
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        particles[1].child1 = particles[3]
        particles[1].child2 = particles[4]
    
        x = trees.BinaryTreesOnAParticleSet(particles, "child1", "child2")
        roots = list(x.iter_roots())
    
        self.assertEqual(len(roots), 1)
        y = [(event, x.mass.value_in(units.kg)) for event, x in roots[0].iter_events()]
        self.assertEqual(y, 
            [
                ('start', 0.0), 
                    ('start', 1.0), 
                        ('start', 3.0), 
                        ('end', 3.0),
                        ('start', 4.0),
                        ('end', 4.0),
                    ('end', 1.0), 
                    ('start', 2.0), 
                    ('end', 2.0), 
                ('end', 0.0)
            ]
        )

    def test6(self):
        
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        particles[1].child1 = particles[3]
        particles[1].child2 = particles[4]
    
        x = trees.BinaryTreesOnAParticleSet(particles, "child1", "child2")
        roots = list(x.iter_roots())
    
        self.assertEqual(len(roots), 1)
        y = [(event, x.mass.value_in(units.kg)) for event, x in roots[0].iter_levels()]
        self.assertEqual(y, 
            [
                (0, 0.0), 
                    (1, 1.0), 
                        (2, 3.0), 
                        (2, 4.0),
                    (1, 2.0),
            ]
        )
        
    def test7(self):
        
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        particles[1].child1 = particles[3]
        particles[1].child2 = particles[4]
    
        x = trees.BinaryTreesOnAParticleSet(particles, "child1", "child2")
        roots = list(x.iter_roots())
    
        self.assertEqual(len(roots), 1)
        binary = roots[0]
        output = ''
        for level, particle in binary.iter_levels():
            output += '..' * level
            output += str(particle.mass.value_in(units.kg))
            output += '\n'
        
        print(output)
        self.assertEqual(output, """0.0
..1.0
....3.0
....4.0
..2.0
""")

    
    def test8(self):
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        particles[1].child1 = particles[3]
        particles[1].child2 = particles[4]
        
    
        x = trees.BinaryTreesOnAParticleSet(particles, "child1", "child2")
        self.assertEqual(len(list(x.iter_roots())), 1)
        self.assertEqual(len(x.particles_not_in_a_multiple()), 5)
            

class TestChildTree(amusetest.TestCase):
    
    def test1(self):
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        
        self.assertEqual(particles[0].mass, 0 | units.kg)
        self.assertEqual(particles[0].child1.mass, 1 | units.kg)
        self.assertEqual(particles[1].child1, None)
    
        tree = particles.as_binary_tree()
        self.assertFalse(tree.is_leaf())
        self.assertEqual(len(list(tree.iter_children())), 8)
        self.assertEqual(len(list(tree.iter_branches())), 1)
        self.assertEqual(len(list(tree.iter_leafs())), 7)
        branches = list(tree.iter_branches())
        self.assertEqual(len(list(branches[0].iter_children())), 2)
        self.assertEqual(len(list(branches[0].iter_branches())), 0)
        self.assertEqual(len(list(branches[0].iter_leafs())), 2)
        
        
    
    def test2(self):
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        
        self.assertEqual(particles[0].mass, 0 | units.kg)
        self.assertEqual(particles[0].child1.mass, 1 | units.kg)
        self.assertEqual(particles[1].child1, None)
    
        tree = particles.as_binary_tree()
        self.assertFalse(tree.is_leaf())
        self.assertEqual(len(list(tree.iter_descendant_leafs())), 9)
        self.assertEqual(len(list(tree.iter_descendant_branches())), 1)
        branches = list(tree.iter_branches())
        self.assertEqual(len(list(branches[0].iter_descendant_leafs())), 2)
        self.assertEqual(len(list(branches[0].iter_descendant_branches())), 0)


    def test3(self):
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        particles[1].child1 = particles[3]
        particles[1].child2 = particles[4]
        particles[2].child1 = particles[5]
        particles[2].child2 = particles[6]
        
    
        tree = particles.as_binary_tree()
        self.assertFalse(tree.is_leaf())
        self.assertEqual(len(list(tree.iter_children())), 4)
        self.assertEqual(len(list(tree.iter_branches())), 1)
        self.assertEqual(len(list(tree.iter_leafs())), 3)
        self.assertEqual(len(list(tree.iter_descendant_leafs())), 7)
        self.assertEqual(len(list(tree.iter_descendant_branches())), 3)
        branches = list(tree.iter_branches())
        self.assertEqual(len(list(branches[0].iter_children())), 2)
        self.assertEqual(len(list(branches[0].iter_branches())), 2)
        self.assertEqual(len(list(branches[0].iter_leafs())), 0)
        self.assertEqual(len(list(branches[0].iter_descendant_leafs())), 4)
        self.assertEqual(len(list(branches[0].iter_descendant_branches())), 2)
        
    def test4(self):
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[2]
        particles[0].child2 = particles[1]
        particles[1].child1 = particles[4]
        particles[1].child2 = particles[3]
        particles[2].child1 = particles[6]
        particles[2].child2 = particles[5]
        
        tree = particles.as_binary_tree()
        masses = [] | units.kg
        stack = list(reversed(tree.get_children()))
        while stack:
            x = stack.pop()
            masses.append(x.particle.mass)
            stack.extend(reversed(x.get_children()))
            
        self.assertEqual(masses, [0,2,6,5,1,4,3,7,8,9] | units.kg)
         
    def test5(self):
        
        particles = Particles(10)
        particles.mass = list(range(10)) | units.kg
        particles[0].child1 = particles[1]
        particles[0].child2 = particles[2]
        particles[1].child1 = particles[3]
        particles[1].child2 = particles[4]
    
        x = particles.as_binary_tree()
        roots = list(x.iter_branches())
    
        self.assertEqual(len(roots), 1)
        binary = roots[0]
        output = ''
        for level, node in binary.iter_levels():
            output += '..' * level
            output += str(node.particle.mass.value_in(units.kg))
            output += '\n'
        
        print(output)
        self.assertEqual(output, """0.0
..1.0
....3.0
....4.0
..2.0
""")

