from amuse.test import amusetest

from amuse.support.units import units
from amuse.support.data import core
import numpy

class TestAddParticles(amusetest.TestCase):
    
    def test1(self):
        print
        print "Test1: create a particle subset by adding a particle to a set."
        original_set = core.Particles(4)
        original_set.x = [1.0, 2.0, -789.0, 3.0] | units.m
        set = original_set[:2]
        particle = original_set[3]
        self.assertTrue(isinstance(set, core.ParticlesSubset))
        self.assertTrue(isinstance(particle, core.Particle))
        print "You can add a particle to a set by using 'add':"
        print "new_set = set.add(particle)"
        new_set = set.add(particle)
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),len(set)+1)
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0]|units.m))
        print "You can also add a particle to a set by using '+':"
        print "new_set = set + particle"
        new_set = set + particle
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),len(set)+1)
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0]|units.m))
        print "You can also add a particle to a set by using '+=':"
        print "set += particle"
        set += particle
        self.assertTrue(isinstance(set, core.ParticlesSubset))
        self.assertEqual(len(set),3)
        self.assertEqual(set.x, ([1.0, 2.0, 3.0]|units.m))
    
    def test2(self):
        print
        print "Test2: create a particle subset by adding a set to a set."
        original_set = core.Particles(5)
        original_set.x = [1.0, 2.0, -789.0, 3.0, 4.0] | units.m
        set1 = original_set[:2]
        set2 = original_set[3:]
        self.assertTrue(isinstance(set1, core.ParticlesSubset))
        self.assertTrue(isinstance(set2, core.ParticlesSubset))
        print "You can add a set to a set by using 'add':"
        print "new_set = set1.add(set2)"
        new_set = set1.add(set2)
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),len(set1)+len(set2))
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
        print "You can also add a set to a set by using '+':"
        print "new_set = set1 + set2"
        new_set = set1 + set2
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),len(set1)+len(set2))
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
        print "You can also add a set to a set by using '+=':"
        print "set1 += set2"
        set1 += set2
        self.assertTrue(isinstance(set1, core.ParticlesSubset))
        self.assertEqual(len(set1),4)
        self.assertEqual(set1.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
    
    def test3(self):
        print
        print "Test3: create a particle superset by adding a particle to a set."
        set = core.Particles(2)
        set.x = [1.0, 2.0] | units.m
        particle = core.Particle()
        particle.x = 3.0 | units.m
        print "You can add a particle to a set by using 'add':"
        print "superset = set.add(particle, creat_super=True)"
        superset = set.add(particle, creat_super=True)
        self.assertTrue(isinstance(superset, core.ParticlesSuperset))
        self.assertEqual(len(superset),len(set)+1)
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0]|units.m))
        print "You can add a set to a set by using 'add':"
        print "superset = set.add(set2, creat_super=True)"
        set2 = core.Particles(2)
        set2.x = [3.0, 4.0] | units.m
        superset = set.add(set2, creat_super=True)
        self.assertTrue(isinstance(superset, core.ParticlesSuperset))
        self.assertEqual(len(superset),len(set)+len(set2))
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
    
    def test4(self):
        print
        print "Test4: check if the particle is already part of the set."
        set = core.Particles(2)
        particle = core.Particle()
        set = set.add(particle, creat_super=True)
        print "Should not be able to add the same particle twice... ",
        try:
            incorrect_set = set.add(particle, creat_super=True)
            print "oops!"
            self.fail("Should not be able to add the same particle twice.")
        except Exception as ex:
            self.assertEquals("Unable to add a particle, because it was "
                "already part of this set.", str(ex))
            print "ok!"
        self.assertEqual(len(set),3)
        other_set = core.Particles(2)
        other_set = other_set.add(particle, creat_super=True)
        print "The particle is now a member of both sets, and thus the sets "
        print "can't be combined anymore... ",
        try:
            incorrect_set = set.add(other_set, creat_super=True)
            print "oops!"
            self.fail("Should not be able to add the same particle twice.")
        except Exception as ex:
            self.assertEquals("Unable to add a particle, because it was "
                "already part of this set.", str(ex))
            print "ok!"
    
    def test5(self):
        print
        print "Test5: recursive addition, create a new superset from supersets."
        particle = core.Particle()
        set1 = core.Particles(2)
        set2 = core.Particles(2)
        set3 = core.Particles(2)
        set4 = core.Particles(2)
        superset1 = set1.add(set2, creat_super=True)
        superset2 = set3.add(set4, creat_super=True)
        for x in [particle, set3, superset2]:
            supersuperset = superset1.add(x, creat_super=True)
            self.assertTrue(isinstance(supersuperset, core.ParticlesSuperset))
            self.assertEqual(len(supersuperset),len(superset1)+numpy.size(x.key))
    
    def test6(self):
        print
        print "Test6: check if the particle belongs to the same particle set as self."
        set1 = core.Particles(2)
        set2 = core.Particles(2)
        particle = set2[0]
        print "Should not be able to create a subset from particles " \
            "belonging to separate particle sets.. ",
        try:
            incorrect_set = set1.add(set2, creat_super=False)
            print "oops!"
            self.fail("Should not be able to create this subset.")
        except Exception as ex:
            self.assertEquals("Can't create new subset from particles belonging to "
                "separate particle sets. Try creating a superset instead.", str(ex))
            print "ok!"
        try:
            incorrect_set = set1.add(particle, creat_super=False)
            self.fail("Should not be able to create this subset.")
        except Exception as ex:
            self.assertEquals("Can't create new subset from particles belonging to "
                "separate particle sets. Try creating a superset instead.", str(ex))
    
    def test7(self):
        print
        print "Test7: add a particle (set) to a particle."
        original_set = core.Particles(4)
        particle1 = original_set[0]
        particle2 = original_set[1]
        set = original_set[2:]
        self.assertTrue(isinstance(particle1, core.Particle))
        self.assertTrue(isinstance(particle2, core.Particle))
        self.assertTrue(isinstance(set, core.ParticlesSubset))
        new_set = particle1.add(particle2)
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),2)
        new_set = particle1.add(set)
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),3)
        new_set = particle1 + particle2
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),2)
        new_set = particle1 + set
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),3)
    
class TestSubtractParticles(amusetest.TestCase):
    
    def test1(self):
        print
        print "Test1: create a particle subset by removing a particle from a set."
        set = core.Particles(4)
        set.x = [1.0, 2.0, -789.0, 3.0] | units.m
        particle = set[2]
        self.assertTrue(isinstance(particle, core.Particle))
        print "You can remove a particle from a set by using 'sub':"
        print "new_set = set.sub(particle)"
        new_set = set.sub(particle)
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),len(set)-1)
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0]|units.m))
        print "You can also remove a particle from a set by using '-':"
        print "new_set = set - particle"
        new_set = set - particle
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),len(set)-1)
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0]|units.m))
        print "You can also remove a particle from a set by using '-=':"
        print "set -= particle"
        set -= particle
        self.assertTrue(isinstance(set, core.ParticlesSubset))
        self.assertEqual(len(set),3)
        self.assertEqual(set.x, ([1.0, 2.0, 3.0]|units.m))
    
    def test2(self):
        print
        print "Test2: create a particle subset by removing a set from a set."
        set1 = core.Particles(5)
        set1.x = [1.0, 2.0, -789.0, 3.0, 4.0, -456.0] | units.m
        set2 = set1[2::3]
        self.assertTrue(isinstance(set1, core.Particles))
        self.assertTrue(isinstance(set2, core.ParticlesSubset))
        print "You can remove a set to a set by using 'sub':"
        print "new_set = set1.sub(set2)"
        new_set = set1.sub(set2)
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),len(set1)-len(set2))
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
        print "You can also remove a set to a set by using '-':"
        print "new_set = set1 - set2"
        new_set = set1 - set2
        self.assertTrue(isinstance(new_set, core.ParticlesSubset))
        self.assertEqual(len(new_set),len(set1)-len(set2))
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
        print "You can also remove a set to a set by using '-=':"
        print "set1 -= set2"
        set1 -= set2
        self.assertTrue(isinstance(set1, core.ParticlesSubset))
        self.assertEqual(len(set1),4)
        self.assertEqual(set1.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
    
    def test3(self):
        print
        print "Test3: check if the particle is actually part of the set."
        set = core.Particles(2)
        particle = core.Particle()
        print "Should not be able to subtract a particle that is not in the set... ",
        try:
            incorrect_set = set - particle
            print "oops!"
            self.fail("Should not be able to subtract a particle that is not in the set.")
        except Exception as ex:
            self.assertEquals("Unable to subtract a particle, because "
                "it is not part of this set.", str(ex))
            print "ok!"
    
    def test4(self):
        print
        print "Test4: recursive subtraction, remove particles until the set is empty."
        set = core.Particles(10)
        self.assertEqual(len(set), 10)
        while len(set):
            set -= set[0]
        self.assertEqual(len(set), 0)
    
    def test5(self):
        print
        print "Test5: check if it's possible to subtract particle(s) from a particle."
        particle = core.Particle()
        print "Should not be able to subtract particle(s) from a particle... ",
        try:
            incorrect_set = particle - particle
            print "oops!"
            self.fail("Should not be able to subtract particle(s) from a particle.")
        except Exception as ex:
            self.assertEquals("Cannot subtract particle(s) from a particle.", str(ex))
            print "ok!"
        particle2 = core.Particle()
        try:
            incorrect_set = particle - particle2
            self.fail("Should not be able to subtract particle(s) from a particle.")
        except Exception as ex:
            self.assertEquals("Cannot subtract particle(s) from a particle.", str(ex))
    
