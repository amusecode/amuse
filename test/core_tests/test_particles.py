from amuse.test import amusetest

from amuse.support.exceptions import AmuseException



from amuse.support.interface import InCodeComponentImplementation

import numpy
import time
import sys
import pickle

from amuse.units import units
from amuse.units import quantities
from amuse.units import constants
from amuse.units import nbody_system
from amuse import datamodel

class TestParticles(amusetest.TestCase):
    
    def test1(self):
        
        particles = datamodel.Particles(2)
        for i in range(3):
            particles.mass = (i * 1.0) | units.kg
            particles.savepoint((i + 1) * 1.0 | units.s)
        
        masses = particles[0].get_timeline_of_attribute("mass")
        self.assertEquals(len(masses), 3)
        self.assertEquals(masses[0][0], 1.0 | units.s)
        self.assertEquals(masses[1][0], 2.0 | units.s)
        self.assertEquals(masses[2][0], 3.0 | units.s)

    def test2(self):
        particles = datamodel.Particles(2)
        particles.mass = [1,1]|units.kg
        heavy = particles.select(lambda m: m>1|units.kg, ["mass"])
        self.assertTrue(heavy.is_empty())
        
    def test3(self):
        particles = datamodel.Particles(2)
        for i in range(3):
            particles.mass = (i * 2.0) | units.kg
            particles.savepoint((i + 1) * 1.0 | units.s)
        
        state0 = particles.get_state_at_timestamp(2.0 | units.s)
        self.assertEquals(state0[0].mass, 2.0 | units.kg)
        state0 = particles.get_state_at_timestamp(2.4 | units.s)
        self.assertEquals(state0[0].mass, 2.0 | units.kg)
        state0 = particles.get_state_at_timestamp(2.6 | units.s)
        self.assertEquals(state0[0].mass, 4.0 | units.kg)
        
    
    def test4(self):
        
        particles = datamodel.Particles(2)
        particles.mass = 1.0 | units.kg
        particles.vy = 1.0 | units.m / units.s
        particles.vx = 0.0 | units.m / units.s
        particles.vz = 0.0 | units.m / units.s
        particles.x = [0.0, 1.0] | units.m
        particles.y = 0.0 | units.m
        particles.z = 0.0 | units.m
        self.assertEquals(particles.kinetic_energy(), 1.0 | units.J)
        self.assertEquals(particles.potential_energy(), -1.0 * constants.G * (1.0 | units.kg ** 2 / units.m))
        self.assertEquals(particles.center_of_mass().x, 0.5 | units.m)
        self.assertEquals(particles.center_of_mass().y, 0.0 | units.m)
        self.assertEquals(particles.center_of_mass().z, 0.0 | units.m)
        self.assertEquals(particles.center_of_mass_velocity().x, 0.0 | units.m / units.s)
        self.assertEquals(particles.center_of_mass_velocity().y, 1.0 | units.m / units.s)
        self.assertEquals(particles.center_of_mass_velocity().z, 0.0 | units.m / units.s)
        
    def test5(self):
        
        particles = datamodel.Particles(4)
        particles.mass = 1.0 | units.kg
        particles[2].mass = 2.0 | units.kg
        subset = particles.select_array(lambda mass: mass > 1.0 | units.kg, ["mass"])
        self.assertEquals(len(subset), 1)
        copyof_subset = subset.copy()
        self.assertEquals(len(copyof_subset), 1)
        

    def test6(self):
        
        particles = datamodel.Particles(2)
        particles.mass = 1.0 | units.kg
        particles.vy = 1.0 | units.m / units.s
        particles.vx = 0.0 | units.m / units.s
        particles.vz = 0.0 | units.m / units.s
        particles.x = [0.0, 1.0] | units.m
        particles.y = 0.0 | units.m
        particles.z = 0.0 | units.m
        
        Ek=0. | units.J
        Ep=0. | units.J
        for x in particles:
            Ek+=x.mass*x.specific_kinetic_energy()
            Ep+=x.mass*x.potential()/2
        self.assertEquals(particles.kinetic_energy(), Ek)
        self.assertEquals(particles.potential_energy(), Ep)
        
        self.assertEquals((particles.mass*particles.specific_kinetic_energy()).sum(), Ek)
        self.assertEquals(0.5*(particles.mass*particles.potential()).sum(), Ep)

    
    def test7(self):
        
        particles = datamodel.Particles(10)
        for i in range(10):
            particles[i].mass = (i * 1.0) | units.kg

    def test8(self):
        particles = datamodel.Particles(4)
        particles.mass = [1.,2.,3.,4.] | units.kg
        subset = particles[0:2]
        self.assertEquals(len(subset), 2)
        self.assertTrue(str(subset).find('kg') > 0)

    def test9(self):
        particles = datamodel.Particles(4)
        particles.mass = [1.,2.,3.,4.] | units.kg
        particles.vy = [1.,2.,3.,4.] | units.m / units.s
        particles.vx = [0.,1.,2.,3.] | units.m / units.s
        particles.vz = [1.,2.,3.,4.] | units.m / units.s
        particles.x = [0., 1.,2.,3.] | units.m
        particles.y = [4., 3.,2.,1.] | units.m
        particles.z = [4., 3.,2.,1.] | units.m
        
        Ek=0. | units.J
        Ep=0. | units.J
        for x in particles:
            Ek+=x.mass*x.specific_kinetic_energy()
            Ep+=x.mass*x.potential()/2
        self.assertEquals(particles.kinetic_energy(), Ek)
        self.assertEquals(particles.potential_energy(), Ep)
        
        self.assertEquals((particles.mass*particles.specific_kinetic_energy()).sum(), Ek)
        self.assertEquals(0.5*(particles.mass*particles.potential()).sum(), Ep)

    def test10(self):
        particles = datamodel.Particles(4)
        particles.mass = [1.,2.,3.,4.] | units.kg
        subset = particles[0:2]
        self.assertEquals(len(subset), 2)
        masses, = subset.get_values_in_store(None, ['mass'])
        self.assertEquals(len(masses), 2)
        subset.mass = 5.0 | units.kg
        self.assertAlmostRelativeEquals( particles.mass , [5.,5.,3.,4.] | units.kg)
        subset.set_values_in_store(None, ['mass'], [6.0 | units.kg,])
        self.assertAlmostRelativeEquals( particles.mass , [6., 6.,3.,4.] | units.kg)
        
        
    def test11(self):
        particles = datamodel.Particles(3)
        print range(3)| units.kg
        particles.mass = range(3)| units.kg
        particles.nounit = range(3)
        print particles.mass
        self.assertAlmostRelativeEquals( particles.mass ,range(3)| units.kg)
        self.assertAlmostRelativeEquals( particles.nounit , [0, 1, 2])
        
        
        
    def xtest12(self):
        particles = datamodel.Particles(3)
        particles.mass = range(3)| units.kg
        particles[0].child = particles[1]
        particles[1].child = particles[2]
        print particles[0].child
        print particles[1].child
        print particles[2].child
        self.assertEquals(particles[0].child, particles[1])
        self.assertEquals(particles[1].child, particles[2])
        self.assertEquals(particles[2].child, None)

        
    def test14(self):
        
        particles1 = datamodel.Particles(2)
        particles1.mass = 1| units.kg
    
        particles2 = particles1.copy()
        particles2.stellar_mass = range(2) * (1.0 | units.kg)
        particles2.mass = 0 | units.kg
        
        self.assertFalse(hasattr(particles1, 'stellar_mass'))
        
        channel = particles2.new_channel_to(particles1)
        channel.copy_overlapping_attributes()
        
        self.assertFalse(hasattr(particles1, 'stellar_mass'))
        self.assertEquals(particles1.mass, [0,0] | units.kg)
        
        channel.copy_all_attributes()
        
        self.assertTrue(hasattr(particles1, 'stellar_mass'))
        self.assertEquals(particles1.mass, [0,0] | units.kg)
        
    def test15(self):
        particles = datamodel.Particles(10)
        
        # List of scalar quantities:
        masses = [i | nbody_system.mass for i in range(10)]
        particles.mass = masses
        print particles.mass
        self.assertEqual(particles.mass, range(10) | nbody_system.mass)
        
        # List of vector quantities:
        positions = [(i, 2*i, 3*i) | units.m for i in range(10)]
        particles.position = positions
        self.assertEqual(particles.position, [(i, 2*i, 3*i) for i in range(10)] | units.m)
        
        # Even lists of tuples of quantities work:
        positions = [(i | units.m, i | units.cm, i | units.km) for i in range(10)]
        particles.position = positions
        self.assertEqual(particles.position, [(i, 0.01*i, 1000*i) for i in range(10)] | units.m)
    
    def test16(self):
        particles = datamodel.Particles(2)
        particles.add_vector_attribute('unitless', ['u1', 'u2', 'u3'])
        particles.unitless = [[1,2,3],[4,5,6]]
        self.assertEquals(particles[0].u1, 1)
        self.assertEquals(particles[0].u3, 3)
        self.assertEquals(particles[1].u2, 5)
        self.assertEquals(particles[0].unitless, [1,2,3])
        self.assertEquals(particles.unitless, [[1,2,3],[4,5,6]])
        particles[0].unitless = [7,8,9]
        self.assertEquals(particles[0].unitless, [7,8,9])
        self.assertEquals(particles[0].u1, 7)
        self.assertEquals(particles.unitless, [ [7,8,9],[4,5,6]])
        
        
    def test17(self):
        particles = datamodel.Particles(2)
        particles.a = [1.0, 2.0]
        self.assertEquals(particles[0].a, 1.0)
        particles.b = 3.0
        self.assertEquals(particles[0].b, 3.0)
        self.assertEquals(particles[1].b, 3.0)
        self.assertAlmostRelativeEquals(particles.b, [3.0, 3.0])
        particles[0].b = 4.0
        self.assertAlmostRelativeEquals(particles.b, [4.0, 3.0])
        particles[1].a = 5.0
        self.assertAlmostRelativeEquals(particles.a, [1.0, 5.0])
        
    def test18(self):
        particles = datamodel.Particles(3)
        particles.mass = [1.0, 2.0, 3.0] | units.kg
        self.assertTrue(particles[1] in particles)
        self.assertFalse(datamodel.Particle in particles)
        
    def test19(self):
        particles = datamodel.Particles(mass = [1.0,2.0,3.0] | units.kg, radius = 1 | units.m)
        self.assertEquals(len(particles), 3)
        self.assertAlmostRelativeEquals(particles.mass,  [1.0,2.0,3.0] | units.kg)
        self.assertAlmostRelativeEquals(particles.radius,  [1.0,1.0,1.0] | units.m)
        particles = datamodel.Particles(b = 1, a = [1.0,2.0,3.0] )
        self.assertEquals(len(particles), 3)
        self.assertAlmostRelativeEquals(particles.a,  [1.0,2.0,3.0])
        self.assertAlmostRelativeEquals(particles.b,  [1.0,1.0,1.0])
        particles = datamodel.Particles(size = 3,b = 1, a = [1.0,2.0,3.0] )
        self.assertEquals(len(particles), 3)
        self.assertAlmostRelativeEquals(particles.a,  [1.0,2.0,3.0])
        self.assertAlmostRelativeEquals(particles.b,  [1.0,1.0,1.0])
        particles = datamodel.Particles(size = 3,b = 1, a = 2 )
        self.assertEquals(len(particles), 3)
        self.assertAlmostRelativeEquals(particles.a,  [2,2,2])
        self.assertAlmostRelativeEquals(particles.b,  [1.0,1.0,1.0])
    
    def test20(self):
        
        particles = datamodel.Particles(3)
        particles.a = [1.0, 2.0, 3,0]
        self.assertEquals(particles[0].a, 1.0)
        # should be able to set the attribute with a unitless array
        # this code will be obsoleted when units.none is completely gone
        particles.a = [4.0, 5.0, 6,0] | units.none 
        self.assertEquals(particles[0].a, 4.0)
        particles.b = [1, 2, 3] | units.none 
        self.assertEquals(particles[0].b, 1 | units.none)
        particles.b = [4, 5, 6]
        self.assertEquals(particles[0].b, 4 | units.none)
    
    def test21(self):
        particles = datamodel.Particles(3)
        particles.z = quantities.zero
        print particles.z, particles[0].z
        self.assertEquals(particles[0].z, quantities.zero)
        particles.z += 1 | units.kg
        self.assertEquals(particles[0].z, 1 | units.kg)
        print particles[0].z
        self.assertEquals(particles.z.unit, units.kg)
        
        
    def test22(self):
        particles = datamodel.Particles(3)
        particles.z = quantities.zero
        particles.z = 1 | units.kg
        self.assertEquals(particles[0].z, 1 | units.kg)
        self.assertEquals(particles.z.unit, units.kg)
        
class TestStars(amusetest.TestCase):

    def test1(self):
        stars = datamodel.Particles(2)
        stars[0].mass = 10 | units.g
        stars[0].position = units.m(numpy.array([1.0,2.0,1.0])) 
        stars[1].mass = 10 | units.g
        stars[1].position = units.m(numpy.array([0.0,0.0,0.0]))
        self.assertEquals(0.5 | units.m, stars.center_of_mass().x)

    def test2(self):
        stars = datamodel.Particles(2)
        stars[0].mass = 10 | units.g
        stars[0].velocity = (units.m / units.s)(numpy.array([1.0,2.0,1.0])) 
        stars[1].mass = 10 | units.g
        stars[1].velocity = (units.m / units.s)(numpy.array([0.0,0.0,0.0]))
        self.assertEquals(0.5 | units.m / units.s, stars.center_of_mass_velocity().x)
        self.assertEquals(1.0 | units.m / units.s, stars.center_of_mass_velocity().y)
        
    
    def test3(self):
        stars = datamodel.Particles(2)
        stars[0].mass = 10 | units.g
        stars[0].velocity = (units.m / units.s)(numpy.array([1.0,2.0,1.0])) 
        stars[0].position = units.m(numpy.array([1.0,2.0,1.0])) 
        stars[1].mass = 10 | units.g
        stars[1].velocity = (units.m / units.s)(numpy.array([0.0,0.0,0.0]))
        stars[1].position = units.m(numpy.array([0.0,0.0,0.0]))
        
        
        self.assertEquals(stars.mass[0], 10|units.g)
        
        self.assertEquals(stars.position[0], [1.0, 2.0, 1.0] | units.m)
        self.assertEquals(stars.velocity[0], [1.0, 2.0, 1.0] | units.m / units.s)


    
    def test4(self):
        stars = datamodel.Particles(2)
        stars[0].x = 1.0  | units.km
        stars[0].y = 2000.0 | units.m
        stars[0].z = 3500.0 | units.m
        
        self.assertEquals(stars.position[0], [1000.0, 2000.0, 3500.0] | units.m)

    def test5(self):
        stars = datamodel.Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[1,2,3],[4,5,6]] | units.km
        
        self.assertEquals(stars[0].md, [1,2,3] | units.km)
        self.assertEquals(stars[1].md, [4,5,6] | units.km)
        
        self.assertEquals(stars.md[0], [1,2,3] | units.km)
        self.assertEquals(stars.md[1], [4,5,6] | units.km)
        stars.md_nounit = [[1,2,3],[4,5,6]] 
        
        self.assertEquals(stars[0].md_nounit, [1,2,3])
        self.assertEquals(stars[1].md_nounit, [4,5,6])
        self.assertEquals(stars.md_nounit[0], [1,2,3])
        self.assertEquals(stars.md_nounit[1], [4,5,6])
        
        stars.md2 = [[[1,3],[2,4],[3,6]],[[4,2],[5,3],[6,1]]] | units.km
        
        self.assertEquals(stars[0].md2, [[1,3],[2,4],[3,6]] | units.km)
        self.assertEquals(stars[1].md2, [[4,2],[5,3],[6,1]] | units.km)
        self.assertEquals(stars.md2[0], [[1,3],[2,4],[3,6]] | units.km)
        self.assertEquals(stars.md2[1], [[4,2],[5,3],[6,1]] | units.km)
        
    def test6(self):
        stars = datamodel.Particles(2)
        stars.x = 1.0  | units.km
        stars.md = [[1,2,3],[4,5,6]] | units.km
        
        self.assertEquals(stars[0].md, [1,2,3] | units.km)
        self.assertEquals(stars[1].md, [4,5,6] | units.km)
        
        copy = stars.copy_to_memory()
        
        self.assertEquals(copy[0].md, [1,2,3] | units.km)
        self.assertEquals(copy[1].md, [4,5,6] | units.km)
        copy[0].md = [7,8,9] | units.km
        self.assertEquals(stars[0].md, [1,2,3] | units.km)
        self.assertEquals(copy[0].md, [7,8,9] | units.km)
        
class TestParticlesWithBinding(amusetest.TestCase):
    class TestLegacyCode(object):
            
        def __init__(self):
            self.masses = {}
            self.links = {}
            self.grids = {}
            
        def get_mass(self, id):
            masses = []
            errors = []
            for x in id:
                masses.append(self.masses[x])
                errors.append(0)
            return ( masses, errors, )
        
        def set_mass(self, id, mass):
            for i,m in zip(id,mass):
                self.masses[i] = m
                
            return ( [0] * len(id),)
            
            
        def get_link(self, id):
            result = []
            errors = []
            for x in id:
                result.append(self.links[x])
                errors.append(0)
            return ( result, errors, )
        
        def set_link(self, id, link):
            for i,l in zip(id,link):
                self.links[i] = l
                
            return ( [0] * len(id),)
            
            
            
        def get_grid(self, id, index1, index2):
            result = []
            errors = []
            for x,i1,i2 in zip(id, index1, index2):
                result.append(self.grids[x][i1][i2])
                errors.append(0)
            return ( result, errors, )
            
        def set_grid(self, id, index1, index2, value):
            errors = []
            for x,i1,i2, v in zip(id, index1, index2, value):
                self.grids[x][i1][i2] = v
                errors.append(0)
            return ( 0, )
        
        def get_grid_range(self):
            return ( 0, 3, 0, 2 )
            
        def new_particle(self, mass):
            ids = []
            errors = []
            
            for x in mass:
                id = len(self.masses)
                self.masses[id]  = x
                self.links[id]  = -1
                self.grids[id] = numpy.arange(4*3).reshape(4,3)
                ids.append(id)
                errors.append(0)
                
            return (ids, errors)
        
        def delete_particle(self, ids):
            errors = []
            for x in ids:
                del self.masses[x]
                del self.links[x]
                errors.append(0)
            return errors
            
        def get_number_of_particles(self):
            return (len(self.masses), 0)
            
        set_state = set_mass
        get_state = get_mass
        
    class TestInterface(InCodeComponentImplementation):
        
        def __init__(self):
            InCodeComponentImplementation.__init__(self, TestParticlesWithBinding.TestLegacyCode())
        
        def define_methods(self, handler):
            handler.add_method('get_mass',(handler.NO_UNIT,), (units.g, handler.ERROR_CODE))
            handler.add_method('set_mass',(handler.NO_UNIT, units.g,), (handler.ERROR_CODE,))
            handler.add_method('get_link',(handler.NO_UNIT,), (handler.LINK('particles'), handler.ERROR_CODE))
            handler.add_method('set_link',(handler.NO_UNIT, handler.LINK('particles'),), (handler.ERROR_CODE,))
            handler.add_method('get_grid',(handler.NO_UNIT,handler.NO_UNIT,handler.NO_UNIT), (units.g, handler.ERROR_CODE))
            handler.add_method('set_grid',(handler.NO_UNIT,handler.NO_UNIT,handler.NO_UNIT, units.g), ( handler.ERROR_CODE))
        
            handler.add_method('new_particle',(units.g,), (handler.INDEX, handler.ERROR_CODE))
            handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
            handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
            
            
        def define_particle_sets(self, handler):
            handler.define_set('particles', 'id')
            handler.set_new('particles', 'new_particle')
            handler.set_delete('particles', 'delete_particle')
            handler.add_setter('particles', 'set_mass')
            handler.add_getter('particles', 'get_mass', names = ('mass',))
            handler.add_setter('particles', 'set_link')
            handler.add_getter('particles', 'get_link', names = ('link',))
            handler.add_gridded_getter('particles', 'get_grid','get_grid_range', names = ('grid',))
            handler.add_gridded_setter('particles', 'set_grid','get_grid_range', names = ('grid',))
        
    
    def test1(self):
        interface = self.TestInterface()
        interface.particles.add_particles_to_store(
            [1,2],
            ["mass"],
            [[3.0, 4.0] | units.kg]
        )
        
        remote_particles = interface.particles
        local_particles = remote_particles.copy()
        
        channel = remote_particles.new_channel_to(local_particles)
        
        self.assertEquals(interface.masses[0], 3000)
        
        interface.masses[0] = 3500
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3.5)
        self.assertEquals(local_particles[0].mass.value_in(units.kg), 3.0)
        
        channel.copy_attributes(["mass"])
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3.5)
        self.assertEquals(local_particles[0].mass.value_in(units.kg), 3.5)
        
        
    
    
    def test2(self):
        interface = self.TestInterface()
        interface.particles.add_particles_to_store(
            [1, 2],
            ["mass"],
            [[3.0, 4.0] | units.kg]
        )
        
        remote_particles = interface.particles
        
        self.assertEquals(len(remote_particles), 2)
        
        interface.particles.add_particles_to_store(
            [3, 4],
            ["mass"],
            [[5.0, 6.0] | units.kg]
        )
        
        self.assertEquals(len(remote_particles), 4)
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3)
        self.assertEquals(remote_particles[2].mass.value_in(units.kg), 5)
        
        interface.particles.remove_particles_from_store((1,3))
        
        self.assertEquals(len(remote_particles), 2)
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 4)
        self.assertEquals(remote_particles[1].mass.value_in(units.kg), 6)
        
    
    def test3(self):
        interface = self.TestInterface()
        
        local_particles = datamodel.Particles(2)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEquals(len(remote_particles), 2)
        
        local_particles2 = datamodel.Particles(2)
        local_particles2.mass = units.kg.new_quantity([5.0, 6.0])
       
        remote_particles.add_particles(local_particles2)
        
        self.assertEquals(len(remote_particles), 4)
        
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 3)
        self.assertEquals(remote_particles[2].mass.value_in(units.kg), 5)
        
        keys = remote_particles.get_all_keys_in_store()
        interface.particles.remove_particles_from_store((keys[0], keys[2]))
        
        self.assertEquals(len(remote_particles), 2)
        self.assertEquals(remote_particles[0].mass.value_in(units.kg), 4)
        self.assertEquals(remote_particles[1].mass.value_in(units.kg), 6)
        
        
    
    def test4(self):
        interface = self.TestInterface()
        
        local_particles = datamodel.Particles(2)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        
        self.assertEquals(len(remote_particles), 2)
        
        particle = datamodel.Particle()
        particle.mass = units.g.new_quantity(10.0)
        local_particles.add_particle(particle)
        
        self.assertEquals(len(local_particles), 3)
        
    
    
    def test5(self):
        interface = self.TestInterface()
        
        local_particles1 = datamodel.Particles(2)
        local_particles1.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles1)
        
        local_particles2 = datamodel.Particles(3)
        local_particles2.mass = units.kg.new_quantity([5.0, 6.0, 7.0])
        
        local_particles1.add_particles(local_particles2)
        
        self.assertEquals(len(local_particles1), 5)
        self.assertEquals(len(local_particles2), 3)
        
        local_particles1.synchronize_to(remote_particles)
        
        local_particles1.remove_particles_from_store([local_particles1.get_all_keys_in_store()[0]])
        local_particles1.synchronize_to(remote_particles)
        self.assertEquals(len(remote_particles), 4)
        
    
    def test6(self):
        interface = self.TestInterface()
        
        local_particles1 = datamodel.Particles(2)
        local_particles1.mass = units.kg.new_quantity([3.0, 4.0])
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles1)
        
        local_particle = datamodel.Particle()
        local_particle.mass = units.kg.new_quantity(5.0)
        
        local_particles1.add_particle(local_particle)
        
        self.assertEquals(len(local_particles1), 3)
        self.assertEquals(len(remote_particles), 2)
        
        local_particles1.synchronize_to(remote_particles)
        
        self.assertEquals(len(remote_particles), 3)
        
        remote_particles.remove_particle(local_particle)
        self.assertEquals(len(remote_particles), 2)
        
        
    def test7(self):
        interface = self.TestInterface()
        
        local_particles = datamodel.Particles(2)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0])
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        
        local_particles.unknown_attribute = [3.0, 4.0] | units.m
        local_particles.mass = [1,3] | units.kg
        
        channel = local_particles.new_channel_to(remote_particles)
        self.assertRaises(Exception, channel.copy_all_attributes)
        
        channel.copy()
        
        self.assertEquals(remote_particles.mass, local_particles.mass)

            
    def test8(self):
        interface = self.TestInterface()
        
        local_particles1 = datamodel.Particles(2)
        local_particles1.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles1)
        
        interface.set_link([0], remote_particles[1].as_set())
        self.assertEquals(interface.legacy_interface.links[0], 1)
        p = interface.get_link([0])
        self.assertEquals(p[0], remote_particles[1])
        self.assertEquals(p[0], local_particles1[1])
        
        
    def test9(self):
        interface = self.TestInterface()
        
        local_particles1 = datamodel.Particles(2)
        local_particles1.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles1)
        self.assertEquals(remote_particles.link[0], None)
        self.assertEquals(remote_particles.link[1], None)
        remote_particles[0].link = remote_particles[1]
        self.assertEquals(remote_particles.link[0], remote_particles[1])
        self.assertEquals(remote_particles.link[0], local_particles1[1])
        self.assertEquals(remote_particles.link[1], None)
        channel = remote_particles.new_channel_to(local_particles1)
        channel.copy_all_attributes()
        self.assertEquals(local_particles1.link[0], remote_particles[1])
        self.assertEquals(local_particles1.link[0], local_particles1[1])
        self.assertEquals(local_particles1.link[1], None)
        
        
    def test10(self):
        interface = self.TestInterface()
        
        local_particles = datamodel.Particles(2)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        self.assertEquals(remote_particles[0].grid[1][2], 5 | units.g)
        self.assertEquals(remote_particles.grid[0][1][2], 5 | units.g)
        
    def test11(self):
        interface = self.TestInterface()
        
        local_particles = datamodel.Particles(2)
        local_particles.mass = units.kg.new_quantity([3.0, 4.0])
        
        
        remote_particles = interface.particles
        remote_particles.add_particles(local_particles)
        
        channel = remote_particles.new_channel_to(local_particles)
        channel.copy()
        
        self.assertEquals(local_particles[0].grid[1][2], 5 | units.g)
        self.assertEquals(local_particles.grid[0][1][2], 5 | units.g)
        
        channel = local_particles.new_channel_to(remote_particles)
        channel.copy()
        
class TestParticlesWithUnitsConverted(amusetest.TestCase):
    
    def test1(self):
        stars = datamodel.Particles(2)
        stars[0].mass = 10 | units.g
        stars[1].mass = 20 | units.g
        
        class LengthMassConverter(object):
            "source == length, m"
            "target == mass, g"
            source_unit = units.m
            target_unit = units.g
            
            def from_source_to_target(self, quantity):
                value = quantity.value_in(self.source_unit)
                return self.target_unit.new_quantity(value) 
                
            def from_target_to_source(self, quantity):
                value = quantity.value_in(self.target_unit)
                return self.source_unit.new_quantity(value) 
        
        converted_stars = datamodel.ParticlesWithUnitsConverted(stars, LengthMassConverter())
        
        self.assertEquals(stars[0].mass, 10 | units.g)
        self.assertEquals(converted_stars[0].mass, 10 | units.m)
        
        
        converted_stars[0].mass = 30 | units.m
        
        self.assertEquals(stars[0].mass, 30 | units.g)
        
    
    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        stars = datamodel.Particles(1)
        stars[0].mass = 10 | nbody_system.mass
        stars[0].x = 10.0 | nbody_system.length
        stars[0].y = 20.0 | nbody_system.length
        stars[0].z = 30.0 | nbody_system.length
        
        
        converted_stars = datamodel.ParticlesWithUnitsConverted(
            stars, 
            convert_nbody.as_converter_from_si_to_generic())
        
        self.assertEquals(stars[0].mass, 10 | nbody_system.mass)
        self.assertAlmostEqual(converted_stars[0].mass, 100.0 | units.kg, 5)
        
        converted_star = converted_stars[0]
        
        expected = units.m.new_quantity([50.0, 100.0, 150.0])
        for i in range(3):
            self.assertAlmostEqual(converted_star.position[i], expected[i], 6)
            
        converted_star.position = [100.0, 200.0, 300.0] | units.m
        star = stars[0]
        
        expected = nbody_system.length([20.0, 40.0, 60.0])
        for i in range(3):
            self.assertAlmostEqual(star.position[i], expected[i], 6)
            
        
        

class TestParticlesWithReferences(amusetest.TestCase):
    
    def test1(self):
        
        set = datamodel.Particles(3)
        parent = set[0]
        child1 = set[1]
        child2 = set[2]
        
        set.mass = [2,3,4] | units.kg
        parent.child1 = child1
        parent.child2 = child2
        child1.sibling = child2
        
        print set.child1
        self.assertEquals(len(set.child1), 3)
        self.assertEquals(len(set.child1.compress()), 1)
        self.assertEquals(len(set.child2), 3)
        self.assertEquals(len(set.child2.compress()), 1)
        
        self.assertAlmostRelativeEqual(parent.child1.mass, 3 | units.kg)

    def test2(self):
        
        
        set = datamodel.Particles(3)
        parent = set[0]
        child1 = set[1]
        child2 = set[2]
        
        set.mass = [2,3,4] | units.kg
        
        
        parent.child1 = child1
        parent.child2 = child2
        child1.sibling = child2
        
        
        convert_nbody = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        set = datamodel.ParticlesWithUnitsConverted(
            set, 
            convert_nbody.as_converter_from_generic_to_si()
        )
        self.assertEquals(len(set.child1), 3)
        self.assertEquals(len(set.child1.compress()), 1)
        self.assertEquals(len(set.child2), 3)
        self.assertEquals(len(set.child2.compress()), 1)
        
        self.assertAlmostRelativeEqual(set[0].child1.mass, 0.3 | nbody_system.mass)
        self.assertAlmostRelativeEqual(set[0].child1.mass, 0.3 | nbody_system.mass)
        
class TestParticlesWithChildren(amusetest.TestCase):
    
    def test1(self):
        
        all = datamodel.Particles(3)
        
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        
        parent.add_child(child1)
        parent.add_child(child2)
        
        
        self.assertEquals(child1.parent.key, parent.key)
        self.assertEquals(child2.parent.key, parent.key)
        
        children = parent.children()
        
        self.assertEquals(len(children), 2)
        
    def test2(self):
        
        code1 = TestParticlesWithBinding.TestInterface()
        code2 = TestParticlesWithBinding.TestInterface()
        
        
        all = datamodel.Particles(3)
        all.mass = [4.0, 3.0, 1.0] | units.kg
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        
        parent.add_child(child1)
        parent.add_child(child2)
        outputstr = str(all)
        print outputstr
        self.assertTrue("  4.000e+00           --" in outputstr)
        code1.particles.add_particles(parent.as_set())
        code2.particles.add_particles(parent.children())
        
        self.assertEquals(len(parent.children()), 2)
        self.assertEquals(len(code1.particles), 1)
        self.assertEquals(len(code2.particles), 2)
        
        code1.legacy_interface.set_mass([0], [10000.0])
        code2.legacy_interface.set_mass([0], [9000.0])
        
        self.assertEquals(parent.mass, 4.0 | units.kg)
        self.assertEquals(child1.mass, 3.0 | units.kg)
        
        code1.particles.copy_values_of_all_attributes_to(all)
        
        self.assertEquals(parent.mass, 10.0 | units.kg)
        self.assertEquals(child1.mass, 3.0 | units.kg)
        
        
        code2.particles.copy_values_of_all_attributes_to(all)
        
        self.assertEquals(parent.mass, 10.0 | units.kg)
        self.assertEquals(child1.mass, 9.0 | units.kg)
        
        
    
    def test3(self):
        
        code1 = TestParticlesWithBinding.TestInterface()
        code2 = TestParticlesWithBinding.TestInterface()
        
        
        all = datamodel.Particles(5)
        all.mass = [4.0, 3.0, 1.0, 6.0, 5.0] | units.kg
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        
        parent.add_child(child1)
        parent.add_child(child2)
        
        self.assertEquals(parent.parent, None)
        self.assertEquals(child1.parent, parent)
        self.assertEquals(child2.parent, parent)
        
        
        all_except_children = all.difference(parent.children())
        code1.particles.add_particles(all_except_children)
        code2.particles.add_particles(parent.children())
        
        self.assertEquals(len(code1.particles), 3)
        self.assertEquals(len(code2.particles), 2)
        
    def test4(self):
        all = datamodel.Particles(5)
        all.mass = [1.0, 2.0, 3.0, 4.0, 5.0] | units.kg
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        child3 = all[3]
        child4 = all[4]
        
        parent.add_child(child1)
        child1.add_child(child2)
        child2.add_child(child3)
        child3.add_child(child4)
        
        self.assertEquals(len(parent.children()), 1)
        self.assertEquals(len(parent.descendents()), 4)
        self.assertEquals(len(child1.descendents()), 3)
        self.assertEquals(len(child2.descendents()), 2)
        self.assertEquals(len(child3.descendents()), 1)
        
    def test5(self):
        all = datamodel.Particles(5)
        all.mass = [1.0, 2.0, 3.0, 4.0, 5.0] | units.kg
        parent = all[0]
        child1 = all[1]
        child2 = all[2]
        child3 = all[3]
        child4 = all[4]
        
        parent.add_child(child1)
        child1.add_child(child2)
        child2.add_child(child3)
        child3.add_child(child4)
        
        copy = all.copy()
        
        self.assertEquals(copy[0].parent, None)
        self.assertEquals(copy[1].parent, copy[0])
        self.assertEquals(copy[2].parent, copy[1])
        self.assertEquals(copy[3].parent, copy[2])
        self.assertEquals(copy[4].parent, copy[3])
        self.assertEquals(copy[1].parent, parent)
        self.assertEquals(len(copy[0].descendents()), 4)
        
class TestParticlesSuperset(amusetest.TestCase):
    
    def test1(self):
        print "Test1: getting attributes of a particle superset."
        superset = datamodel.Particles(2)
        superset.x = [1.0, 2.0] | units.m
        set2 = datamodel.Particles(2)
        set2.x = [4.0, 5.0] | units.m
        particle1 = datamodel.Particle()
        particle1.x = 3.0 | units.m
        particle2 = datamodel.Particle()
        particle2.x = 6.0 | units.m
        for x in [particle1, set2, particle2]:
            superset = datamodel.ParticlesSuperset([superset, x.as_set()])
        self.assertTrue(isinstance(superset, datamodel.ParticlesSuperset))
        self.assertEqual(len(superset),6)
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]|units.m))
        self.assertEqual(superset[1].x,  2.0 |units.m)
        self.assertEqual(superset[0:2][1].x,  2.0 |units.m)
        # Check whether it returns the right value for the right key when the key order is 'random'
        dictionary = dict(zip(superset.key, superset.x))
        sorted_keys = sorted(superset.get_all_keys_in_store())
        sorted_values = superset.get_values_in_store(sorted_keys,['x'])[0]
        for key, value in zip(sorted_keys, sorted_values):
            self.assertEqual(dictionary[key],value)
    
    def test2(self):
        print "Test2: setting attributes of a particle superset."
        superset = datamodel.Particles(2)
        set2 = datamodel.Particles(2)
        particle1 = datamodel.Particle()
        particle2 = datamodel.Particle()
        for x in [particle1, set2, particle2]:
            superset = datamodel.ParticlesSuperset([superset, x.as_set()])
        self.assertTrue(isinstance(superset, datamodel.ParticlesSuperset))
        self.assertEqual(len(superset),6)
        superset.x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0] | units.m
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0, 4.0, 5.0, 6.0]|units.m))
        superset.y = 9.0 | units.m
        self.assertEqual(superset.y, ([9.0, 9.0, 9.0, 9.0, 9.0, 9.0]|units.m))
        superset.z = [-1.0, 1.0] | units.m
        self.assertEqual(superset.z, ([-1.0, 1.0, -1.0, 1.0, -1.0, 1.0]|units.m))
        # Check whether it sets the value of the right particle when the key order is 'random'
        sorted_keys = sorted(superset.get_all_keys_in_store())
        sorted_values = superset.get_values_in_store(sorted_keys,['x'])[0]
        superset.set_values_in_store(sorted_keys,['zz'],[sorted_values])
        dictionary = dict(zip(superset.key, superset.zz))
        for key, value in zip(sorted_keys, sorted_values):
            print dictionary[key],value
            self.assertEqual(dictionary[key],value)
            
    def test3(self):
        set1 = datamodel.Particles(2)
        set2 = datamodel.Particles(3)
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [5.0 , 6.0, 7.0] | units.m 
        set2.y = [1.0 , 2.0, 3.0] | units.m
        set1.add_calculated_attribute("xtimesy", lambda x, y : x * y)
        set2.add_calculated_attribute("xtimesy", lambda x, y : x * y)
        superset = datamodel.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 5)
        self.assertEquals(superset.x, ([1.0, 2.0, 5.0, 6.0, 7.0] | units.m))
        self.assertEquals(superset.xtimesy, ([3.0, 8.0, 5.0, 12.0, 21.0] | units.m ** 2))
        
    def test4(self):
        set1 = datamodel.Particles(2)
        set2 = datamodel.Particles(3)
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [5.0 , 6.0, 7.0] | units.m 
        set2.y = [1.0 , 2.0, 3.0] | units.m
        set1.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y - o, lambda ap, p, o : p.x * p.y - o)
        set2.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y + o, lambda ap, p, o : p.x * p.y - o)
        superset = datamodel.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 5)
        self.assertEquals(superset.x, ([1.0, 2.0, 5.0, 6.0, 7.0] | units.m))
        self.assertEquals(superset.xtimesypluso(0.0 | units.m ** 2), ([3.0, 8.0, 5.0, 12.0, 21.0] | units.m ** 2))
        self.assertEquals(superset.xtimesypluso(2.0 | units.m ** 2), ([1.0, 6.0, 7.0, 14.0, 23.0] | units.m ** 2))
        self.assertEquals(superset[0].xtimesypluso(0.0 | units.m ** 2), 3.0 | units.m ** 2)
        self.assertEquals(superset[3].xtimesypluso(0.0 | units.m ** 2), 12.0 | units.m ** 2)
        
    def test5(self):
        set1 = datamodel.Particles(2)
        set2 = datamodel.Particles()
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [] | units.m 
        set2.y = [] | units.m 
        set1.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y - o)
        set2.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y + o)
        superset = datamodel.ParticlesSuperset([set2, set1])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesypluso(0.0 | units.m ** 2), ([3.0, 8.0] | units.m ** 2))
        self.assertEquals(superset.xtimesypluso(2.0 | units.m ** 2), ([1.0, 6.0] | units.m ** 2))
        superset = datamodel.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesypluso(0.0 | units.m ** 2), ([3.0, 8.0] | units.m ** 2))
        self.assertEquals(superset.xtimesypluso(2.0 | units.m ** 2), ([1.0, 6.0] | units.m ** 2))
        
    def test6(self):
        set1 = datamodel.Particles(2)
        set2 = datamodel.Particles()
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [] | units.m 
        set2.y = [] | units.m 
        set1.add_calculated_attribute("xtimesy", lambda x, y : x * y)
        set2.add_calculated_attribute("xtimesy", lambda x, y : x * y)
        superset = datamodel.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesy, ([3.0, 8.0] | units.m ** 2))
        superset = datamodel.ParticlesSuperset([set2, set1])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesy, ([3.0, 8.0] | units.m ** 2))
    
    def test7(self):
        def xtimesy1(x,y):
            raise Exception("error querying function on empty set")
            
        set1 = datamodel.Particles(2)
        set2 = datamodel.Particles()
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [] | units.m 
        set2.y = [] | units.m 
        set1.add_function_attribute("xtimesy", lambda p : p.x * p.y)
        set2.add_function_attribute("xtimesy", xtimesy1)
        superset = datamodel.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesy(), ([3.0, 8.0] | units.m ** 2))
        superset = datamodel.ParticlesSuperset([set2, set1])
        self.assertEquals(len(superset), 2)
        self.assertEquals(superset.x, ([1.0, 2.0] | units.m))
        self.assertEquals(superset.xtimesy(), ([3.0, 8.0] | units.m ** 2))
        
    def test8(self):
        set1 = datamodel.Particles(2)
        set2 = datamodel.Particles(3)
        set1.x = [1.0 , 2.0] | units.m 
        set1.y = [3.0 , 4.0] | units.m
        set2.x = [5.0 , 6.0, 7.0] | units.m 
        set2.y = [1.0 , 2.0, 3.0] | units.m
        set1.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y - o, lambda all, p, o : p.x * p.x - o)
        set2.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y + o, lambda all, p, o : p.x * p.x - o)
        superset = datamodel.ParticlesSuperset([set1, set2])
        self.assertEquals(len(superset), 5)
        self.assertEquals(superset.x, ([1.0, 2.0, 5.0, 6.0, 7.0] | units.m))
        self.assertEquals(superset[1].xtimesypluso(0.0 | units.m ** 2), (4.0 | units.m ** 2))
        self.assertEquals(list(superset)[1].xtimesypluso(0.0 | units.m ** 2), (4.0 | units.m ** 2))
        
    
    def test9(self):
        
        particles1 = datamodel.Particles(2)
        particles1.mass = 10 | units.kg
        particles2 = datamodel.Particles(2)
        particles2.mass = 20 | units.kg
        superset = datamodel.ParticlesSuperset([particles1, particles2])
            
        self.assertTrue(hasattr(superset, 'mass'))
        self.assertFalse(hasattr(superset, 'radius'))
        particles1.radius = 10 | units.m
        self.assertFalse(hasattr(superset, 'radius'))
        particles2.radius = 20 | units.m
        self.assertTrue(hasattr(superset, 'radius'))
        
    def test10(self):
        
        particles1 = datamodel.Particles(2)
        particles1.u1 = 10 
        particles2 = datamodel.Particles(2)
        particles2.u1 = 20 
        superset = datamodel.ParticlesSuperset([particles1, particles2])
            
        self.assertTrue(hasattr(superset, 'u1'))
        self.assertEquals(superset.u1 , [10,10,20,20])
        self.assertFalse(hasattr(superset, 'u2'))
        particles1.u2 = 30
        self.assertFalse(hasattr(superset, 'u2'))
        particles2.u2 = 20
        self.assertTrue(hasattr(superset, 'u2'))
        self.assertEquals(superset.u2 , [30,30,20,20])
        superset.u2 = 15
        self.assertEquals(superset.u2 , [15,15,15,15])
    
    def test11(self):
        print "ParticlesSuperset from one set, of one particle"
        particles1 = datamodel.Particles(1)
        particles1.u = 10
        superset = datamodel.ParticlesSuperset([particles1])
        
        self.assertTrue(hasattr(superset, 'u'))
        self.assertEquals(superset.u, [10])
        self.assertEquals(len(superset.u), 1)
        self.assertEquals(len(superset), 1)
        
        particles1.y = 3 | units.m
        self.assertTrue(hasattr(superset, 'y'))
        self.assertEquals(superset.y , [3] | units.m)
        self.assertEquals(len(superset.y), 1)
        
        particles1.x = [1.0] | units.m
        particles1.y = [3.0] | units.m
        particles1.add_function_attribute("xtimesypluso", lambda p, o : p.x * p.y + o, lambda ap, p, o : p.x * p.y + o)
        superset = datamodel.ParticlesSuperset([particles1])
        self.assertEquals(len(superset), 1)
        self.assertEquals(superset.x, [1.0] | units.m)
        self.assertEquals(superset.xtimesypluso(0.0 | units.m**2), [3.0] | units.m**2)
        self.assertEquals(superset.xtimesypluso(2.0 | units.m**2), [5.0] | units.m**2)
        self.assertEquals(superset[0].xtimesypluso(0.0 | units.m**2), 3.0 | units.m**2)
        self.assertEquals(superset[0].xtimesypluso(2.0 | units.m**2), 5.0 | units.m**2)
    
    def test12(self):
        print "ParticlesSuperset - query"
        set1 = datamodel.Particles(4)
        set1.x = [-1.0, 1.0, 2.0, 3.0] | units.m
        set1.add_function_attribute("greater_than", lambda p, o : p[p.x > o], lambda ap, p, o : p if p.x > o else None)
        superset = datamodel.ParticlesSuperset([set1])
        
        self.assertEqual(len(set1.greater_than(-2.0 | units.m)), 4)
        self.assertTrue(isinstance(set1.greater_than(-2.0 | units.m), datamodel.ParticlesSubset))
        self.assertEqual(len(set1.greater_than(0.0 | units.m)), 3)
        self.assertTrue(isinstance(set1.greater_than(0.0 | units.m), datamodel.ParticlesSubset))
        self.assertTrue(numpy.all(set1.greater_than(0.0 | units.m).x > 0.0 | units.m))
        self.assertEqual(len(set1.greater_than(3.0 | units.m)), 0)
        self.assertTrue(isinstance(set1.greater_than(3.0 | units.m), datamodel.ParticlesSubset))
        
        self.assertEqual(len(superset.greater_than(-2.0 | units.m)), 4)
        self.assertTrue(isinstance(superset.greater_than(-2.0 | units.m), datamodel.ParticlesSuperset))
        self.assertEqual(len(superset.greater_than(0.0 | units.m)), 3)
        self.assertTrue(isinstance(superset.greater_than(0.0 | units.m), datamodel.ParticlesSuperset))
        self.assertTrue(numpy.all(superset.greater_than(0.0 | units.m).x > 0.0 | units.m))
        self.assertEqual(len(superset.greater_than(3.0 | units.m)), 0)
        self.assertTrue(isinstance(superset.greater_than(3.0 | units.m), datamodel.ParticlesSuperset))
        
        self.assertEqual(superset[0].greater_than(-2.0 | units.m), superset[0])
        self.assertTrue(isinstance(superset[0].greater_than(-2.0 | units.m), datamodel.Particle))
        self.assertEqual(superset[0].greater_than(0.0 | units.m), None)
        
        set2 = datamodel.Particles(1)
        set2.x = [4.0] | units.m
        set2.add_function_attribute("greater_than", lambda p, o : p[p.x > o], lambda ap, p, o : p if p.x > o else None)
        superset = datamodel.ParticlesSuperset([set1, set2])
        self.assertEqual(len(superset.greater_than(-2.0 | units.m)), 5)
        self.assertTrue(isinstance(superset.greater_than(-2.0 | units.m), datamodel.ParticlesSuperset))
        self.assertEqual(len(superset.greater_than(3.0 | units.m)), 1)
        self.assertTrue(isinstance(superset.greater_than(3.0 | units.m), datamodel.ParticlesSuperset))
        self.assertEqual(superset.greater_than(2.0 | units.m).x, [3.0, 4.0] | units.m)
        
        superset.add_function_attribute("greater_than", lambda p, o : p[p.x > o], lambda ap, p, o : p if p.x > o else None)
        self.assertEqual(len(superset.greater_than(-2.0 | units.m)), 5)
        self.assertTrue(isinstance(superset.greater_than(-2.0 | units.m), datamodel.ParticlesSubset))
        self.assertEqual(len(superset.greater_than(3.0 | units.m)), 1)
        self.assertTrue(isinstance(superset.greater_than(3.0 | units.m), datamodel.ParticlesSubset))
        self.assertEqual(superset.greater_than(2.0 | units.m).x, [3.0, 4.0] | units.m)
        
    
    
class TestSliceParticles(amusetest.TestCase):
    
    def test1(self):
        print "Test: slice a particle set."
        number_of_particles = 10
        original_set = datamodel.Particles(number_of_particles)
        self.assertTrue(isinstance(original_set, datamodel.Particles))
        self.assertEqual(len(original_set),number_of_particles)
        print "Defining all kind of slices of the particle set..."
        subset1 = original_set[:2]  # contains first two particles
        subset2 = original_set[2:]  # contains the rest of the particles
        odd = original_set[1::2]    # contains all particles with odd indices
        even = original_set[::2]    # contains all particles with even indices
        reverse = original_set[::-1]# contains all particles in reverse order
        all = original_set[:]       # contains all particles
        one = original_set[3]       # contains one particle (Particle)
        another = original_set[5:6] # contains one particle (ParticlesSubset)
        empty = original_set[2:2]   # contains no particle
        print "Checking their type (slicing returns a subset, indexing returns a particle)... ",
        self.assertTrue(isinstance(subset1, datamodel.ParticlesSubset))
        self.assertTrue(isinstance(subset2, datamodel.ParticlesSubset))
        self.assertTrue(isinstance(odd,     datamodel.ParticlesSubset))
        self.assertTrue(isinstance(even,    datamodel.ParticlesSubset))
        self.assertTrue(isinstance(reverse, datamodel.ParticlesSubset))
        self.assertTrue(isinstance(all,     datamodel.ParticlesSubset))
        self.assertTrue(isinstance(one,     datamodel.Particle))
        self.assertTrue(isinstance(another, datamodel.ParticlesSubset))
        self.assertTrue(isinstance(empty,   datamodel.ParticlesSubset))
        print "ok!"
        print "Checking their length... ",
        self.assertEqual(len(subset1), 2)
        self.assertEqual(len(subset2), number_of_particles-2)
        self.assertEqual(len(odd),     int((number_of_particles/2.0)))
        self.assertEqual(len(even),    int((0.5+number_of_particles/2.0)))
        self.assertEqual(len(reverse), number_of_particles)
        self.assertEqual(len(all),     number_of_particles)
        self.assertEqual(len([one]),   1)
        self.assertEqual(len(another), 1)
        self.assertEqual(len(empty),   0)
        print "ok!"
    
class TestAddParticles(amusetest.TestCase):
    
    def test1(self):
        print "Test1: create a particle subset by adding a particle to a set."
        original_set = datamodel.Particles(4)
        original_set.x = [1.0, 2.0, -789.0, 3.0] | units.m
        set = original_set[:2]
        particle = original_set[3]
        self.assertTrue(isinstance(set, datamodel.ParticlesSubset))
        self.assertTrue(isinstance(particle, datamodel.Particle))
        
        new_set = set + particle
        self.assertTrue(isinstance(new_set, datamodel.ParticlesSubset))
        self.assertEqual(len(new_set),len(set)+1)
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0]|units.m))
        
        set += particle
        self.assertTrue(isinstance(set, datamodel.ParticlesSubset))
        self.assertEqual(len(set),3)
        self.assertEqual(set.x, ([1.0, 2.0, 3.0]|units.m))
    
    def test2(self):
        print "Test2: create a particle subset by adding a set to a set."
        original_set = datamodel.Particles(5)
        original_set.x = [1.0, 2.0, -789.0, 3.0, 4.0] | units.m
        set1 = original_set[:2]
        set2 = original_set[3:]
        self.assertTrue(isinstance(set1, datamodel.ParticlesSubset))
        self.assertTrue(isinstance(set2, datamodel.ParticlesSubset))
        
        new_set = set1 + set2
        self.assertTrue(isinstance(new_set, datamodel.ParticlesSubset))
        self.assertEqual(len(new_set),len(set1)+len(set2))
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
        
        set1 += set2
        self.assertTrue(isinstance(set1, datamodel.ParticlesSubset))
        self.assertEqual(len(set1),4)
        self.assertEqual(set1.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
    
    def test3(self):
        print "Test3: create a particle superset by adding a particle to a set."
        set = datamodel.Particles(2)
        set.x = [1.0, 2.0] | units.m
        particle = datamodel.Particle()
        particle.x = 3.0 | units.m
        
        superset = datamodel.ParticlesSuperset([set, particle.as_set()])
        self.assertTrue(isinstance(superset, datamodel.ParticlesSuperset))
        self.assertEqual(len(superset),len(set)+1)
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0]|units.m))
        
        set2 = datamodel.Particles(2)
        set2.x = [3.0, 4.0] | units.m
        superset = datamodel.ParticlesSuperset([set, set2])
        self.assertTrue(isinstance(superset, datamodel.ParticlesSuperset))
        self.assertEqual(len(superset),len(set)+len(set2))
        self.assertEqual(superset.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
    
    def test4(self):
        print "Test4: check if the particle is already part of the set."
        set = datamodel.Particles(2)
        particle = datamodel.Particle()
        set = datamodel.ParticlesSuperset([set, particle.as_set()])
        self.assertRaises(AmuseException, datamodel.ParticlesSuperset, [set, particle.as_set()], 
            expected_message = "Unable to add a particle, because it was already part of this set.")
        self.assertEqual(len(set),3)
        other_set = datamodel.Particles(2)
        other_set = datamodel.ParticlesSuperset([other_set, particle.as_set()])
        self.assertRaises(AmuseException, datamodel.ParticlesSuperset, [set, other_set], 
            expected_message = "Unable to add a particle, because it was already part of this set.")
    
    def test5(self):
        print "Test5: recursive addition, create a new superset from supersets."
        particle = datamodel.Particle()
        set1 = datamodel.Particles(2)
        set2 = datamodel.Particles(2)
        set3 = datamodel.Particles(2)
        set4 = datamodel.Particles(2)
        superset1 = datamodel.ParticlesSuperset([set1, set2])
        superset2 = datamodel.ParticlesSuperset([set3, set4])
        for x in [particle, set3, superset2]:
            supersuperset = datamodel.ParticlesSuperset([superset1, x.as_set()])
            self.assertTrue(isinstance(supersuperset, datamodel.ParticlesSuperset))
            self.assertEqual(len(supersuperset),len(superset1)+len(x.as_set()))
            supersuperset.mass = 1 | units.kg
            self.assertEqual(x.mass, 1 | units.kg)
    
    def test6(self):
        print "Test6: check if the particle belongs to the same particle set as self."
        set1 = datamodel.Particles(2)
        set2 = datamodel.Particles(2)
        particle = set2[0]
        self.assertRaises(AmuseException, lambda: set1 + set2, 
            expected_message = "Can't create new subset from particles belonging to "
            "separate particle sets. Try creating a superset instead.")
        self.assertRaises(AmuseException, lambda: set1 + particle, 
            expected_message = "Can't create new subset from particles belonging to "
            "separate particle sets. Try creating a superset instead.")
    
    def test7(self):
        print "Test7: add a particle (set) to a particle."
        original_set = datamodel.Particles(4)
        particle1 = original_set[0]
        particle2 = original_set[1]
        set = original_set[2:]
        self.assertTrue(isinstance(particle1, datamodel.Particle))
        self.assertTrue(isinstance(particle2, datamodel.Particle))
        self.assertTrue(isinstance(set, datamodel.ParticlesSubset))
        new_set = particle1 + particle2
        self.assertTrue(isinstance(new_set, datamodel.ParticlesSubset))
        self.assertEqual(len(new_set),2)
        new_set = particle1 + set
        self.assertTrue(isinstance(new_set, datamodel.ParticlesSubset))
        self.assertEqual(len(new_set),3)
    
class TestSubtractParticles(amusetest.TestCase):
    
    def test1(self):
        print "Test1: create a particle subset by removing a particle from a set."
        set = datamodel.Particles(4)
        set.x = [1.0, 2.0, -789.0, 3.0] | units.m
        particle = set[2]
        self.assertTrue(isinstance(particle, datamodel.Particle))
        
        new_set = set - particle
        self.assertTrue(isinstance(new_set, datamodel.ParticlesSubset))
        self.assertEqual(len(new_set),len(set)-1)
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0]|units.m))
        
        set -= particle
        self.assertTrue(isinstance(set, datamodel.ParticlesSubset))
        self.assertEqual(len(set),3)
        self.assertEqual(set.x, ([1.0, 2.0, 3.0]|units.m))
    
    def test2(self):
        print "Test2: create a particle subset by removing a set from a set."
        set1 = datamodel.Particles(5)
        set1.x = [1.0, 2.0, -789.0, 3.0, 4.0, -456.0] | units.m
        set2 = set1[2::3]
        self.assertTrue(isinstance(set1, datamodel.Particles))
        self.assertTrue(isinstance(set2, datamodel.ParticlesSubset))
        
        new_set = set1 - set2
        self.assertTrue(isinstance(new_set, datamodel.ParticlesSubset))
        self.assertEqual(len(new_set),len(set1)-len(set2))
        self.assertEqual(new_set.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
        
        set1 -= set2
        self.assertTrue(isinstance(set1, datamodel.ParticlesSubset))
        self.assertEqual(len(set1),4)
        self.assertEqual(set1.x, ([1.0, 2.0, 3.0, 4.0]|units.m))
    
    def test3(self):
        print "Test3: check if the particle is actually part of the set."
        set = datamodel.Particles(2)
        particle = datamodel.Particle()
        self.assertRaises(AmuseException, lambda: set - particle, 
            expected_message = "Unable to subtract a particle, because it is not part of this set.")
    
    def test4(self):
        print "Test4: recursive subtraction, remove particles until the set is empty."
        set = datamodel.Particles(10)
        self.assertEqual(len(set), 10)
        while len(set):
            set -= set[0]
        self.assertEqual(len(set), 0)
    
    def test5(self):
        print "Test5: check if it's possible to subtract particle(s) from a particle."
        particle = datamodel.Particle()
        self.assertRaises(AmuseException, lambda: particle - particle, 
            expected_message = "Cannot subtract particle(s) from a particle.")
        particle2 = datamodel.Particle()
        self.assertRaises(AmuseException, lambda: particle - particle2, 
            expected_message = "Cannot subtract particle(s) from a particle.")
    

class TestIterateOverParticles(amusetest.TestCase):
    
    def iterate_over_array(self, particles):
        for x in particles:
            x.radius
    
    def iterate_over_particles(self, particles):
        for x in particles:
            x.key


    def test1(self):
        self.total_number_of_points = 10000
    
        class Test(object):
            def __init__(self):
                self.radius = 1.0
    
        particles = [Test() for x in range(self.total_number_of_points)]
        particles = datamodel.Particles(self.total_number_of_points)
        particles.radius = 2.0 | nbody_system.length
        t0 = time.time()
        self.iterate_over_particles1(particles)
        t1 = time.time()
        dt0 = t1 - t0
        
        particles = datamodel.Particles(self.total_number_of_points)
        particles.radius = 2.0 | nbody_system.length
        t0 = time.time()
        self.iterate_over_particles2(particles)
        t1 = time.time()
        dt1 = t1 - t0
        
        print dt0, dt1,  dt1 / dt0
    
        self.assertTrue((dt1 / dt0) < 20)

    def iterate_over_particles1(self, particles):
        for x in particles:
            x.key
    
    def iterate_over_particles2(self, particles):
        for x in particles:
            x.radius
    
    def test2(self):
        self.total_number_of_points = 10000
    
        class Test(object):
            def __init__(self):
                self.radius = 1.0
    
        particles = [Test() for x in range(self.total_number_of_points)]
        
        particles = datamodel.Particles(self.total_number_of_points)
        particles.radius = 2.0 | nbody_system.length
        t0 = time.time()
        #self.iterate_over_array(particles)
        for key in particles.get_all_keys_in_store():
            key
        t1 = time.time()
        dt0 = t1 - t0
        
        particles = datamodel.Particles(self.total_number_of_points)
        particles.radius = 2.0 | nbody_system.length
        t0 = time.time()
        self.iterate_over_particles2(particles)
        t1 = time.time()
        dt1 = t1 - t0
        
        print  dt0, dt1, dt1 / dt0
    
        self.assertTrue((dt1 / dt0) < 400)   
             
class TestParticlesIndexing(amusetest.TestCase):

    def test1(self):
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        self.assertAlmostRelativeEquals(particles[5].mass, 5 | units.kg)
        self.assertAlmostRelativeEquals(particles[[1,3,2,6]].mass, [1,3,2,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,2,4,6,8] | units.kg)

    def test2(self):
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        
        self.assertAlmostRelativeEquals(particles[5:].mass,  [5,6,7,8,9] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3].mass,  [1,2] | units.kg)
        self.assertAlmostRelativeEquals(particles[:5].mass,  [0,1,2,3,4] | units.kg)

    def test3(self):
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        
        self.assertAlmostRelativeEquals(particles[5:][:2].mass,  [5,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3][:1].mass,  1 | units.kg)
        self.assertAlmostRelativeEquals(particles[:5][2:].mass,  [2,3,4] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,2] | units.kg)
        
    def test4(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5].mass, 0.5 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[[1,3,2,6]].mass, [0.1,0.3,0.2,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,.2,.4,.6,.8] | nbody_system.mass)
   
    def test5(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1|  nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5:].mass,  [0.5,0.6,0.7,0.8,0.9] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3].mass,  [0.1,0.2] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5].mass,  [0,0.1,0.2,0.3,0.4] | nbody_system.mass)
        
    def test6(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        self.assertAlmostRelativeEquals(converted[5:][:2].mass, [0.5,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3][:1].mass,  0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5][2:].mass,  [0.2,0.3,0.4] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,0.2] | nbody_system.mass)
    
    def test7(self):
        
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        self.assertAlmostRelativeEquals(particles[5].mass, 5 | units.kg)
        self.assertAlmostRelativeEquals(particles[[1,3,2,6]].mass, [1,3,2,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,2,4,6,8] | units.kg)


    def test8(self):
        
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])

        self.assertAlmostRelativeEquals(particles[5:].mass,  [5,6,7,8,9] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3].mass,  [1,2] | units.kg)
        self.assertAlmostRelativeEquals(particles[:5].mass,  [0,1,2,3,4] | units.kg)

    def test9(self):
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])

        
        self.assertAlmostRelativeEquals(particles[5:][0].mass,  5 | units.kg)
        self.assertAlmostRelativeEquals(particles[5:][:2].mass,  [5,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3][:1].mass,  1 | units.kg)
        self.assertAlmostRelativeEquals(particles[:5][2:].mass,  [2,3,4] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,2] | units.kg)
           
    def test10(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5].mass, 0.5 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[[1,2,3]].mass, [0.1,0.2,0.3] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,.2,.4,.6,.8] | nbody_system.mass)
 
    def test11(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        self.assertAlmostRelativeEquals(converted[5:][:2].mass, [0.5,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3][:1].mass,  0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5][2:].mass,  [0.2,0.3,0.4] | nbody_system.mass)
        

    def test12(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        self.assertAlmostRelativeEquals(converted[5:][:2].mass,  [0.5,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3][:1].mass,  0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5][2:].mass,  [0.2,0.3,0.4] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,0.2] | nbody_system.mass)
    

class TestParticlesIndexingWithBinding(amusetest.TestCase):
    class TestLegacyCode(object):
            
        def __init__(self, offset=2):
            self.masses = {}
            self.offset = offset
            
        def get_mass(self, id):
            masses = []
            errors = []
            
            for x in id:
                masses.append(self.masses[x-self.offset])
                errors.append(0)
            return ( masses, errors, )
        
        def set_mass(self, id, mass):
            for i,m in zip(id,mass):
                self.masses[i-self.offset] = m
                
            return ( [0] * len(id),)
            
        def new_particle(self, mass):
            ids = []
            errors = []
            
            for x in mass:
                id = len(self.masses)
                self.masses[len(self.masses)]  = x
                ids.append(id + self.offset)
                errors.append(0)
            return (ids, errors)
        
        def delete_particle(self, ids):
            errors = []
            for x in ids:
                del self.masses[x]
                errors.append(0)
            return errors
            
        def get_number_of_particles(self):
            return (len(self.masses), 0)
            
        set_state = set_mass
        get_state = get_mass
        
    class TestInterface(InCodeComponentImplementation):
        
        def __init__(self, offset = 2):
            InCodeComponentImplementation.__init__(self, TestParticlesIndexingWithBinding.TestLegacyCode(offset=offset))
        
        def define_methods(self, handler):
            handler.add_method('get_mass',(handler.NO_UNIT,), (units.g, handler.ERROR_CODE))
            handler.add_method('set_mass',(handler.NO_UNIT, units.g,), (handler.ERROR_CODE,))
            handler.add_method('new_particle',(units.g,), (handler.INDEX, handler.ERROR_CODE))
            handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
            handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        def define_particle_sets(self, handler):
            handler.define_set('particles', 'id')
            handler.set_new('particles', 'new_particle')
            handler.set_delete('particles', 'delete_particle')
            handler.add_setter('particles', 'set_mass')
            handler.add_getter('particles', 'get_mass', names = ('mass',))
        
    def test1(self):
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        self.assertAlmostRelativeEquals(particles[5].mass, 5 | units.kg)
        self.assertAlmostRelativeEquals(particles[[1,3,2,6]].mass, [1,3,2,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,2,4,6,8] | units.kg)

    def test2(self):
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        self.assertAlmostRelativeEquals(particles[5:].mass,  [5,6,7,8,9] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3].mass,  [1,2] | units.kg)
        self.assertAlmostRelativeEquals(particles[:5].mass,  [0,1,2,3,4] | units.kg)

    def test3(self):
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        self.assertAlmostRelativeEquals(particles[5:][:2].mass,  [5,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3][:1].mass,  1 | units.kg)
        self.assertAlmostRelativeEquals(particles[:5][2:].mass,  [2,3,4] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,2] | units.kg)
        
    def test4(self):
        
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5].mass, 0.5 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[[1,3,2,6]].mass, [0.1,0.3,0.2,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,.2,.4,.6,.8] | nbody_system.mass)
   
    def test5(self):
        
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1|  nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5:].mass,  [0.5,0.6,0.7,0.8,0.9] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3].mass,  [0.1,0.2] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5].mass,  [0,0.1,0.2,0.3,0.4] | nbody_system.mass)
        
    def test6(self):
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        self.assertAlmostRelativeEquals(converted[5:][:2].mass, [0.5,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3][:1].mass,  0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5][2:].mass,  [0.2,0.3,0.4] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,0.2] | nbody_system.mass)
    
    def test7(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        self.assertAlmostRelativeEquals(particles[5].mass, 5 | units.kg)
        self.assertAlmostRelativeEquals(particles[[1,3,2,6]].mass, [1,3,2,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,2,4,6,8] | units.kg)


    def test8(self):
        
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])

        self.assertAlmostRelativeEquals(particles[5:].mass,  [5,6,7,8,9] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3].mass,  [1,2] | units.kg)
        self.assertAlmostRelativeEquals(particles[:5].mass,  [0,1,2,3,4] | units.kg)

    def test9(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])

        
        self.assertAlmostRelativeEquals(particles[5:][0].mass,  5 | units.kg)
        self.assertAlmostRelativeEquals(particles[5:][:2].mass,  [5,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3][:1].mass,  1 | units.kg)
        self.assertAlmostRelativeEquals(particles[:5][2:].mass,  [2,3,4] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,2] | units.kg)
           
    def test10(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5].mass, 0.5 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[[1,2,3]].mass, [0.1,0.2,0.3] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,.2,.4,.6,.8] | nbody_system.mass)
 
    def test11(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        self.assertAlmostRelativeEquals(converted[5:][:2].mass, [0.5,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3][:1].mass,  0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5][2:].mass,  [0.2,0.3,0.4] | nbody_system.mass)
        

    def test12(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        self.assertAlmostRelativeEquals(converted[5:][:2].mass,  [0.5,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3][:1].mass,  0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5][2:].mass,  [0.2,0.3,0.4] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,0.2] | nbody_system.mass)
    
    def test13(self):
        set1 = self.TestInterface().particles
        set2 = self.TestInterface().particles
        self.assertEquals(set1.mass,  [])
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        self.assertEquals(particles.mass,  [])
        
class TestParticlesIndexingWithSet(amusetest.TestCase):

    #
    # normal set
    #
    def test1(self):
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        particles.multidimensional = numpy.arange(60).reshape(10,2,3)
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        self.assertAlmostRelativeEquals(particles.multidimensional, numpy.arange(60).reshape(10,2,3))
        particles[5].mass = 15 | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,6,7,8,9) | units.kg )
        particles[5].multidimensional = 15 
        self.assertAlmostRelativeEquals(particles.multidimensional[5], [[15,15,15],[15,15,15]] )
        particles.mass = numpy.arange(10) | units.kg
        particles[[1,3,2,6]].mass = (11, 13, 12, 16) | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,12,13,4,5,16,7,8,9) | units.kg )
        particles.mass = numpy.arange(10) | units.kg
        particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass = [10,12,14,16,18] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,1,12,3,14,5,16,7,18,9) | units.kg )

    def test2(self):
        particles = datamodel.Particles(10)
        
        particles.mass = numpy.arange(10) | units.kg
        particles[5:].mass = [15,16,17,18,19] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,16,17,18,19) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[1:3].mass = [11,12] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,12,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[:5].mass = [10,11,12,13,14] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,11,12,13,14,5,6,7,8,9) | units.kg )

    def test3(self):
        particles = datamodel.Particles(10)
        
        particles.mass = numpy.arange(10) | units.kg
        particles[5:][:2].mass = [15,16] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,16,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[1:3][:1].mass =  11 | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,2,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[:5][2:].mass = [12,13,14] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,12,13,14,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass = [10,12] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,1,12,3,4,5,6,7,8,9) | units.kg )
    
    
    #
    # Converted unit sets
    #
    def test4(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        converted[5].mass = 15 | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[[1,3,2,6]].mass = (11, 13, 12, 16) | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,120,130,4,5,160,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass = [10,12,14,16,18] | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,1,120,3,140,5,160,7,180,9) | units.kg )


    def test5(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        
        particles.mass = numpy.arange(10) | units.kg
        converted[5:].mass = [15,16,17,18,19] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,160,170,180,190) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[1:3].mass = [11,12] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,120,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[:5].mass = [10,11,12,13,14] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,110,120,130,140,5,6,7,8,9) | units.kg )
        
    def test6(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        
        particles.mass = numpy.arange(10) | units.kg
        converted[5:][:2].mass = [15,16] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,160,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[1:3][:1].mass =  11 |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,2,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[:5][2:].mass = [12,13,14] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,120,130,140,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass = [10,12] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,1,120,3,4,5,6,7,8,9) | units.kg )
        
        
    
    #
    # Particles super sets
    #
    def test7(self):
        
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        particles[5].mass = 15 | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[[1,3,2,6]].mass = (11, 13, 12, 16) | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,12,13,4,5,16,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass = [10,12,14,16,18] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,1,12,3,14,5,16,7,18,9) | units.kg )

    def test8(self):
        
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])

        particles.mass = numpy.arange(10) | units.kg
        particles[5:].mass = [15,16,17,18,19] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,16,17,18,19) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[1:3].mass = [11,12] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,12,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[:5].mass = [10,11,12,13,14] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,11,12,13,14,5,6,7,8,9) | units.kg )


    def test9(self):
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])

        particles.mass = numpy.arange(10) | units.kg
        particles[5:][0].mass = [15] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[5:][:2].mass = [15,16] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,16,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[1:3][:1].mass =  11 | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,2,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[:5][2:].mass = [12,13,14] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,12,13,14,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass = [10,12] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,1,12,3,4,5,6,7,8,9) | units.kg )
          
    def test10(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5].mass, 0.5 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[[1,2,3]].mass, [0.1,0.2,0.3] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,.2,.4,.6,.8] | nbody_system.mass)
 
    #
    # superset with conversion
    #
    def test11(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[5:].mass = [15,16,17,18,19] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,160,170,180,190) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[1:3].mass = [11,12] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,120,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[:5].mass = [10,11,12,13,14] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,110,120,130,140,5,6,7,8,9) | units.kg )

    def test12(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[5:][0].mass = [15] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[5:][:2].mass = [15,16] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,160,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[1:3][:1].mass =  11 |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,2,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[:5][2:].mass = [12,13,14] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,120,130,140,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass = [10,12] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,1,120,3,4,5,6,7,8,9) | units.kg )
        
class TestParticlesIndexingWithBindingAndSet(amusetest.TestCase):
    class TestLegacyCode(object):
            
        def __init__(self, offset=2):
            self.masses = {}
            self.offset = offset
            
        def get_mass(self, id):
            masses = []
            errors = []
            
            for x in id:
                masses.append(self.masses[x-self.offset])
                errors.append(0)
            return ( masses, errors, )
        
        def set_mass(self, id, mass):
            try:
                
                for i,m in zip(id,mass):
                    self.masses[i-self.offset] = m
            except:
                if len(id) == 1:
                    self.masses[id[0]-self.offset] = mass
                    
            return ( [0] * len(id),)
            
        def new_particle(self, mass):
            ids = []
            errors = []
            
            for x in mass:
                id = len(self.masses)
                self.masses[len(self.masses)]  = x
                ids.append(id + self.offset)
                errors.append(0)
            return (ids, errors)
        
        def delete_particle(self, ids):
            errors = []
            for x in ids:
                del self.masses[x]
                errors.append(0)
            return errors
            
        def get_number_of_particles(self):
            return (len(self.masses), 0)
            
        set_state = set_mass
        get_state = get_mass
        
    class TestInterface(InCodeComponentImplementation):
        
        def __init__(self, offset = 2):
            InCodeComponentImplementation.__init__(self, TestParticlesIndexingWithBindingAndSet.TestLegacyCode(offset=offset))
        
        def define_methods(self, handler):
            handler.add_method('get_mass',(handler.NO_UNIT,), (units.g, handler.ERROR_CODE))
            handler.add_method('set_mass',(handler.NO_UNIT, units.g,), (handler.ERROR_CODE,))
            handler.add_method('new_particle',(units.g,), (handler.INDEX, handler.ERROR_CODE))
            handler.add_method('delete_particle',(handler.NO_UNIT,), (handler.ERROR_CODE,))
            handler.add_method('get_number_of_particles',(), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        def define_particle_sets(self, handler):
            handler.define_set('particles', 'id')
            handler.set_new('particles', 'new_particle')
            handler.set_delete('particles', 'delete_particle')
            handler.add_setter('particles', 'set_mass')
            handler.add_getter('particles', 'get_mass', names = ('mass',))
        
    #
    # normal set
    #
    def test1(self):
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        particles[5].mass = 15 | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,6,7,8,9) | units.kg )
        particles.mass = numpy.arange(10) | units.kg
        particles[[1,3,2,6]].mass = (11, 13, 12, 16) | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,12,13,4,5,16,7,8,9) | units.kg )
        particles.mass = numpy.arange(10) | units.kg
        particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass = [10,12,14,16,18] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,1,12,3,14,5,16,7,18,9) | units.kg )

    def test2(self):
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[5:].mass = [15,16,17,18,19] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,16,17,18,19) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[1:3].mass = [11,12] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,12,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[:5].mass = [10,11,12,13,14] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,11,12,13,14,5,6,7,8,9) | units.kg )

    def test3(self):
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[5:][:2].mass = [15,16] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,16,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[1:3][:1].mass =  11 | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,2,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[:5][2:].mass = [12,13,14] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,12,13,14,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass = [10,12] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,1,12,3,4,5,6,7,8,9) | units.kg )
    
    #
    # units converted
    #
    def test4(self):
        
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        converted[5].mass = 15 | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[[1,3,2,6]].mass = (11, 13, 12, 16) | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,120,130,4,5,160,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass = [10,12,14,16,18] | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,1,120,3,140,5,160,7,180,9) | units.kg )

    def test5(self):
        
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[5:].mass = [15,16,17,18,19] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,160,170,180,190) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[1:3].mass = [11,12] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,120,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[:5].mass = [10,11,12,13,14] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,110,120,130,140,5,6,7,8,9) | units.kg )
        
    def test6(self):
        particles = self.TestInterface().particles
        particles.add_particles_to_store(
             numpy.arange(10),
            ["mass"],
            [ numpy.arange(10) | units.kg]
        )
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[5:][0].mass = [15] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[5:][:2].mass = [15,16] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,160,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[1:3][:1].mass =  11 |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,2,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[:5][2:].mass = [12,13,14] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,120,130,140,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass = [10,12] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,1,120,3,4,5,6,7,8,9) | units.kg )
    
    #
    # supersets
    #
    def test7(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        
        particles[5].mass = 15 | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[[1,3,2,6]].mass = (11, 13, 12, 16) | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,12,13,4,5,16,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass = [10,12,14,16,18] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,1,12,3,14,5,16,7,18,9) | units.kg )

    def test8(self):
        
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])

        particles.mass = numpy.arange(10) | units.kg
        particles[5:].mass = [15,16,17,18,19] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,16,17,18,19) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[1:3].mass = [11,12] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,12,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[:5].mass = [10,11,12,13,14] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,11,12,13,14,5,6,7,8,9) | units.kg )

    def test9(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])

        
        particles.mass = numpy.arange(10) | units.kg
        particles[5:][0].mass = [15] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[5:][:2].mass = [15,16] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,15,16,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[1:3][:1].mass =  11 | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,11,2,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[:5][2:].mass = [12,13,14] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (0,1,12,13,14,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass = [10,12] | units.kg
        self.assertAlmostRelativeEquals(particles.mass, (10,1,12,3,4,5,6,7,8,9) | units.kg )
        
    #
    # superset and unit conversion
    #
    def test10(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        converted[5].mass = 15 | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[[1,3,2,6]].mass = (11, 13, 12, 16) | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,120,130,4,5,160,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass = [10,12,14,16,18] | nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,1,120,3,140,5,160,7,180,9) | units.kg )
 
    def test11(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        particles.mass = numpy.arange(10) | units.kg
        converted[5:].mass = [15,16,17,18,19] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,160,170,180,190) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[1:3].mass = [11,12] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,120,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[:5].mass = [10,11,12,13,14] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,110,120,130,140,5,6,7,8,9) | units.kg )

    def test12(self):
        set1 = self.TestInterface().particles
        set1.add_particles_to_store(
             numpy.arange(5),
            ["mass"],
            [ numpy.arange(5) | units.kg]
        )
        set2 = self.TestInterface().particles
        set2.add_particles_to_store(
             numpy.arange(5)+5,
            ["mass"],
            [ numpy.arange(5) +5 | units.kg]
        )
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        particles.mass = numpy.arange(10) | units.kg
        converted[5:][0].mass = [15] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[5:][:2].mass = [15,16] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,2,3,4,150,160,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[1:3][:1].mass =  11 |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,110,2,3,4,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[:5][2:].mass = [12,13,14] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (0,1,120,130,140,5,6,7,8,9) | units.kg )
        
        particles.mass = numpy.arange(10) | units.kg
        converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass = [10,12] |  nbody_system.mass
        self.assertAlmostRelativeEquals(particles.mass, (100,1,120,3,4,5,6,7,8,9) | units.kg )




class TestParticlesIndexingAfterPickle(amusetest.TestCase):


    def copy_with_pickle(self, set):
        return pickle.loads(pickle.dumps(set))
        
    def test1(self):
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        
        particles = self.copy_with_pickle(particles)
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        self.assertAlmostRelativeEquals(particles[5].mass, 5 | units.kg)
        self.assertAlmostRelativeEquals(particles[[1,3,2,6]].mass, [1,3,2,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,2,4,6,8] | units.kg)

    def test2(self):
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        
        particles = self.copy_with_pickle(particles)
        
        self.assertAlmostRelativeEquals(particles[5:].mass,  [5,6,7,8,9] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3].mass,  [1,2] | units.kg)
        self.assertAlmostRelativeEquals(particles[:5].mass,  [0,1,2,3,4] | units.kg)

    def test3(self):
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        
        particles = self.copy_with_pickle(particles)
        
        self.assertAlmostRelativeEquals(particles[5:][:2].mass,  [5,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3][:1].mass,  1 | units.kg)
        self.assertAlmostRelativeEquals(particles[:5][2:].mass,  [2,3,4] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,2] | units.kg)
        
    def test4(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        converted = self.copy_with_pickle(converted)
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5].mass, 0.5 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[[1,3,2,6]].mass, [0.1,0.3,0.2,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,.2,.4,.6,.8] | nbody_system.mass)
   
    def test5(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        converted = self.copy_with_pickle(converted)
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1|  nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5:].mass,  [0.5,0.6,0.7,0.8,0.9] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3].mass,  [0.1,0.2] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5].mass,  [0,0.1,0.2,0.3,0.4] | nbody_system.mass)
        
    def test6(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        converted = self.copy_with_pickle(converted)
        self.assertAlmostRelativeEquals(converted[5:][:2].mass, [0.5,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3][:1].mass,  0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5][2:].mass,  [0.2,0.3,0.4] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,0.2] | nbody_system.mass)
    
    def test7(self):
        
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        particles = self.copy_with_pickle(particles)
        
        self.assertAlmostRelativeEquals(particles.mass, numpy.arange(10) | units.kg)
        self.assertAlmostRelativeEquals(particles[5].mass, 5 | units.kg)
        self.assertAlmostRelativeEquals(particles[[1,3,2,6]].mass, [1,3,2,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,2,4,6,8] | units.kg)


    def test8(self):
        
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])

        particles = self.copy_with_pickle(particles)
        
        self.assertAlmostRelativeEquals(particles[5:].mass,  [5,6,7,8,9] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3].mass,  [1,2] | units.kg)
        self.assertAlmostRelativeEquals(particles[:5].mass,  [0,1,2,3,4] | units.kg)

    def test9(self):
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])

        particles = self.copy_with_pickle(particles)
        

        
        self.assertAlmostRelativeEquals(particles[5:][0].mass,  5 | units.kg)
        self.assertAlmostRelativeEquals(particles[5:][:2].mass,  [5,6] | units.kg)
        self.assertAlmostRelativeEquals(particles[1:3][:1].mass,  1 | units.kg)
        self.assertAlmostRelativeEquals(particles[:5][2:].mass,  [2,3,4] | units.kg)
        self.assertAlmostRelativeEquals(particles[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,2] | units.kg)
           
    def test10(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
        
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])
        
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )
        
        self.assertAlmostRelativeEquals(converted.mass, numpy.arange(10) * 0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[5].mass, 0.5 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[[1,2,3]].mass, [0.1,0.2,0.3] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])].mass, [0,.2,.4,.6,.8] | nbody_system.mass)
 
    def test11(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        set1 = datamodel.Particles(5)
        set1.mass = numpy.arange(5) | units.kg
        set2 = datamodel.Particles(5)
        set2.mass = numpy.arange(5) + 5 | units.kg
        particles = datamodel.ParticlesSuperset([set1, set2])
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )

        converted = self.copy_with_pickle(converted)
        
        self.assertAlmostRelativeEquals(converted[5:][:2].mass, [0.5,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3][:1].mass,  0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5][2:].mass,  [0.2,0.3,0.4] | nbody_system.mass)
        

    def test12(self):
        converter = nbody_system.nbody_to_si(10 | units.kg, 5 | units.m )
         
        particles = datamodel.Particles(10)
        particles.mass = numpy.arange(10) | units.kg
        converted = datamodel.ParticlesWithUnitsConverted(
                particles,
                converter.as_converter_from_generic_to_si()
        )

        converted = self.copy_with_pickle(converted)
        
        self.assertAlmostRelativeEquals(converted[5:][:2].mass,  [0.5,0.6] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[1:3][:1].mass,  0.1 | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[:5][2:].mass,  [0.2,0.3,0.4] | nbody_system.mass)
        self.assertAlmostRelativeEquals(converted[numpy.array([True,False,True,False,True,False,True,False,True,False])][:2].mass,  [0,0.2] | nbody_system.mass)
    
class TestParticlesOverlay(amusetest.TestCase):
    
    def test1(self):
        set1 = datamodel.Particles(2)
        set1.x = [1.0, 2.0] | units.kg
        set2 = datamodel.ParticlesOverlay(set1)
        set2.y = [4.0, 5.0] | units.m
        self.assertAlmostRelativeEquals(set2.x, [1.0, 2.0] | units.kg)
        self.assertAlmostRelativeEquals(set2.y, [4.0, 5.0] | units.m)
        self.assertTrue('x' in dir(set1))
        self.assertFalse('y' in dir(set1))
        
    def test2(self):
        set1 = datamodel.Particles(2)
        set1.x = [1.0, 2.0] | units.kg
        set2 = datamodel.ParticlesOverlay(set1)
        set2.y = [4.0, 5.0] | units.m
        set2.add_particle(datamodel.Particle(x=3.0 | units.kg, y = 6.0 | units.m))
        self.assertAlmostRelativeEquals(set2.x, [1.0, 2.0, 3.0] | units.kg)
        self.assertAlmostRelativeEquals(set2.y, [4.0, 5.0, 6.0] | units.m)
        self.assertAlmostRelativeEquals(set1.x, [1.0, 2.0, 3.0] | units.kg)
    
    def test3(self):
        set1 = datamodel.Particles(2)
        set1.x = [1.0, 2.0] | units.kg
        set2 = datamodel.ParticlesOverlay(set1)
        set2.y = [4.0, 5.0] | units.m
        particle = set2[0]
        self.assertAlmostRelativeEquals(particle.x, 1.0 | units.kg)
        self.assertAlmostRelativeEquals(particle.y, 4.0 | units.m)
        
        
    def test4(self):
        set1 = datamodel.Particles(2)
        set1.x = [1.0, 2.0] | units.kg
        set2 = datamodel.ParticlesOverlay(set1)
        set2.y = [4.0, 5.0] | units.m
        set2.remove_particle(set2[1])
        self.assertEquals(len(set2), 1)
        self.assertAlmostRelativeEquals(set2.x, [1.0] | units.kg)
        self.assertAlmostRelativeEquals(set2.y, [4.0] | units.m)
        self.assertAlmostRelativeEquals(set1.x, [1.0] | units.kg)

    def test5(self):
        set1 = datamodel.Particles(3)
        set1.x = [1.0, 2.0, 3.0] | units.kg
        set2 = datamodel.ParticlesOverlay(set1)
        set2.y = [4.0, 5.0, 6.0] | units.m
        subset = set2[1:]
        self.assertEquals(len(subset), 2)
        self.assertAlmostRelativeEquals(subset.x, [2.0, 3.0] | units.kg)
        self.assertAlmostRelativeEquals(subset.y, [5.0, 6.0] | units.m)
        xy = subset.get_values_in_store(subset.get_all_keys_in_store(), ['x','y'])
        self.assertAlmostRelativeEquals(xy[0], [2.0, 3.0] | units.kg)
        self.assertAlmostRelativeEquals(xy[1], [5.0, 6.0] | units.m)

    def test6(self):
        set1 = datamodel.Particles(3)
        set1.x = [1.0, 2.0, 3.0] | units.kg
        set2 = datamodel.ParticlesOverlay(set1)
        set2.y = [4.0, 5.0, 6.0] | units.m
        set2.x = [7.0,8.0,9.0] | units.kg
    
        self.assertAlmostRelativeEquals(set2.x, [7.0,8.0,9.0] | units.kg)
        self.assertAlmostRelativeEquals(set1.x, [7.0,8.0,9.0] | units.kg)

