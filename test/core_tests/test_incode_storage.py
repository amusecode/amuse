from amuse.test import amusetest

from amuse.support.units import units
from amuse.support.units import constants
from amuse.support.units import nbody_system

from amuse.support.data.incode_storage import *


import numpy

class TestParticles(amusetest.TestCase):
    
    def test1(self):
        class Code(object):
            def __init__(self):
                # x,y,z,mass
                self.data = []
                self.get_position_called = False
                self.set_position_called = False
                
            def get_number_of_particles(self):
                return  0 if not self.data else len(self.data[0])
                
            def get_position(self,index):
                self.get_position_called = True
                data_to_return = [(self.data[0][i], self.data[1][i], self.data[2][i]) for i in index]
                data_to_return = numpy.asarray(data_to_return).reshape(3,-1)
                return [units.m(x) for x in data_to_return]
                
            def set_position(self,index,x,y,z):
                self.set_position_called = True
                pass
                
            def new_particle(self, x, y, z):
                x = x.value_in(units.m)
                y = y.value_in(units.m)
                z = z.value_in(units.m)
                self.data = [x,y,z]
                return [i for i in range(len(x))]
                
        code = Code()
        storage = InCodeAttributeStorage(
            code,
            NewParticleMethod(code.new_particle,("x","y","z")),
            None,
            code.get_number_of_particles,
            [],
            [ParticleGetAttributesMethod(code.get_position,("x","y","z")),],
            name_of_the_index = "index"
        )
        
        self.assertEquals(len(storage), 0)
        self.assertEquals(storage.get_defined_attribute_names(), ["x","y","z"])
        
        self.assertFalse(code.get_position_called)
        storage.get_values_in_store([],["x","y","z"])
        self.assertFalse(code.get_position_called)
        
        storage.add_particles_to_store(
            [1,2,3,4],
            ["x","y","z"],
            [
                units.m([1,2,3,4]),
                units.m([2,3,4,5]),
                units.m([3,4,5,6])
            ]
        )
        
        self.assertEquals(len(storage), 4)
    
    def test2(self):
        class Code(object):
            def __init__(self):
                # x,y,z,mass
                self.data = []
                self.get_position_called = False
                self.set_position_called = False
                self.get_mass_called = False
                self.set_mass_called = False
                
            def get_number_of_particles(self):
                return  0 if not self.data else len(self.data[0])
                
            def get_position(self,index):
                self.get_position_called = True
                data_to_return = [(self.data[0][i], self.data[1][i], self.data[2][i]) for i in index]
                data_to_return = numpy.asarray(data_to_return).reshape(3,-1)
                return [units.m(x) for x in data_to_return]
            
            def get_mass(self,index):
                self.get_mass_called = True
                data_to_return = [self.data[3][i] for i in index]
                return units.kg(data_to_return)
                
            def set_position(self,index,x,y,z):
                self.set_position_called = True
                pass
                
            def set_mass(self,index,mass):
                self.set_mass_called = True
                pass
                
            def new_particle(self, x, y, z, mass):
                x = x.value_in(units.m)
                y = y.value_in(units.m)
                z = z.value_in(units.m)
                mass = mass.value_in(units.kg)
                self.data = [x,y,z, mass]
                return [i for i in range(len(x))]
                
        code = Code()
        storage = InCodeAttributeStorage(
            code,
            NewParticleMethod(code.new_particle,("x","y","z","mass")),
            None,
            code.get_number_of_particles,
            [],
            [
                ParticleGetAttributesMethod(code.get_position,("x","y","z")),
                ParticleGetAttributesMethod(code.get_mass,("mass",)),
            ],
            name_of_the_index = "index"
        )
        
        storage.add_particles_to_store(
            [1,2,3,4],
            ["x","y","z", "mass"],
            [
                units.m([1,2,3,4]),
                units.m([2,3,4,5]),
                units.m([3,4,5,6]),
                units.kg([13,14,15,16]),
            ]
        )
        
        self.assertEquals(len(storage), 4)
        
        self.assertEquals(storage.get_defined_attribute_names(), [ "mass", "x","y","z"])
        
        self.assertFalse(code.get_position_called)
        self.assertFalse(code.get_mass_called)
        x,y,mass = storage.get_values_in_store([2,3],["x","y","mass"])
        self.assertTrue(code.get_position_called)
        self.assertTrue(code.get_mass_called)
        self.assertEquals(x[1], 3 | units.m)
        self.assertEquals(mass[1], 15 | units.kg)
        

    def test3(self):
        class Code(object):
            def __init__(self):
                # mass
                self.data = []
                self.get_mass_called = False
                self.set_mass_called = False
                
            def get_number_of_particles(self):
                return  0 if not self.data else len(self.data[0])
                
            def get_mass(self,index):
                self.get_mass_called = True
                data_to_return = [self.data[0][i] for i in index]
                return units.kg(data_to_return)
                
            def set_mass(self,index,mass):
                self.set_mass_called = True
                pass
                
            def new_particle(self, mass):
                mass = mass.value_in(units.kg)
                self.data = [mass]
                return [i for i in range(len(mass))]
                
        code = Code()
        storage = InCodeAttributeStorage(
            code,
            NewParticleMethod(code.new_particle,("mass",)),
            None,
            code.get_number_of_particles,
            [],
            [
                ParticleGetAttributesMethod(code.get_mass,("mass",)),
            ],
            name_of_the_index = "index"
        )
        
        storage.add_particles_to_store(
            [1,2,3,4],
            ["mass"],
            [
                units.kg([1,2,3,4]),
            ]
        )
        
        self.assertEquals(len(storage), 4)
        
        self.assertEquals(storage.get_defined_attribute_names(), ["mass",])
        
        index,mass = storage.get_values_in_store([2,3],["index_in_code","mass"])
        self.assertTrue(code.get_mass_called)
        print index, mass
        self.assertEquals(index[0], 1)
        self.assertEquals(mass[0],  2 | units.kg)
        self.assertEquals(index[1], 2)
        self.assertEquals(mass[1],  3 | units.kg)
    
    def test4(self):
        class Code(object):
            def __init__(self):
                # mass
                self.data = []
                self.get_mass_called = False
                self.set_mass_called = False
                self.number_of_particles = 0
                
            def get_number_of_particles(self):
                return  self.number_of_particles
                
            def get_mass(self,index):
                self.get_mass_called = True
                data_to_return = [self.data[i] for i in index]
                return units.kg(data_to_return)
                
            def set_mass(self,index,mass):
                self.set_mass_called = True
                pass
                
            def new_particle(self, mass):
                mass = mass.value_in(units.kg)
                self.data = mass
                self.number_of_particles = len(self.data)
                return [i for i in range(len(mass))]
                
        code = Code()
        storage = InCodeAttributeStorage(
            code,
            NewParticleMethod(code.new_particle,("mass",)),
            None,
            code.get_number_of_particles,
            [],
            [
                ParticleGetAttributesMethod(code.get_mass,("mass",)),
            ],
            name_of_the_index = "index"
        )
        
        storage.add_particles_to_store(
            numpy.asarray([1,2,3,4], dtype='uint64'),
            ["mass"],
            [
                units.kg([1,2,3,4]),
            ]
        )
        
        self.assertEquals(len(storage), 4)
        
        storage._remove_indices([1,2,])
        code.number_of_particles = 2
        index,mass = storage.get_values_in_store([1,4],["index_in_code","mass"])
        
        self.assertEquals(index[0], 0)
        self.assertEquals(index[1], 3)
        self.assertEquals(mass[0],  1 | units.kg)
        self.assertEquals(mass[1],  4 | units.kg)
        
        
        self.assertEquals(len(storage), 2)
        
        storage._add_indices([4,5])
        code.data = numpy.concatenate((code.data, [5, 6]))
        
        code.number_of_particles = 4
        self.assertEquals(len(storage), 4)
        
        mass, = storage.get_values_in_store(storage.particle_keys,["mass"])
        
        self.assertEquals(mass[0],  1 | units.kg)
        self.assertEquals(mass[1],  4 | units.kg)
        self.assertEquals(mass[2],  5 | units.kg)
        self.assertEquals(mass[3],  6 | units.kg)
        
        storage._remove_indices([4,])
        code.number_of_particles = 3
        self.assertEquals(len(storage), 3)
        
        mass, = storage.get_values_in_store(storage.particle_keys,["mass"])
        
        self.assertEquals(mass[0],  1 | units.kg)
        self.assertEquals(mass[1],  4 | units.kg)
        self.assertEquals(mass[2],  6 | units.kg)
    
    def test5(self):
        class Code(object):
            def __init__(self):
                self.data = []
                self.number_of_particles = 0
                
            def get_number_of_particles(self):
                return  self.number_of_particles
                
            def get_mass(self,index):
                data_to_return = [self.data[i][0] for i in index]
                return units.kg(data_to_return)
            
            def get_children(self,index):
                return [(self.data[i][1]) for i in index], [(self.data[i][2]) for i in index]
            
            def new_particle(self, mass):
                mass = mass.value_in(units.kg)
                self.data = [[x,-1,-1] for x in mass]
                self.number_of_particles = len(self.data)
                return [i for i in range(len(mass))]
            
        code = Code()
    
        children_getter = ParticleGetAttributesMethod(
                    code.get_children,
                    ('child1', 'child2',)
        )
        children_getter.index_output_attributes = set(['child1','child2'])
     
        storage = InCodeAttributeStorage(
            code,
            NewParticleMethod(code.new_particle,("mass",)),
            None,
            code.get_number_of_particles,
            [],
            [
                ParticleGetAttributesMethod(code.get_mass,("mass",)),
                children_getter
            ],
            name_of_the_index = "index"
        )
        
        storage.add_particles_to_store(
            numpy.asarray([100,200,300,400], dtype='uint64'),
            ["mass"],
            [
                units.kg([1,2,3,4]),
            ]
        )
        
        self.assertEquals(len(storage), 4)
    
        mass = storage.get_values_in_store([100,400],["mass",])[0]
        print mass
        self.assertEquals(mass[0], 1.0 | units.kg)
        self.assertEquals(mass[1], 4.0 | units.kg)
    
        code.data[0][1] = 1
        code.data[0][2] = 2
    
        child1,child2 = storage.get_values_in_store([100],['child1', 'child2'])
        print child1
    
        self.assertEquals(child1[0].number, 200)
        self.assertEquals(child2[0].number, 300)
    
        x = Particles(storage  = storage)
    
        self.assertEquals(x[0].mass, 1.0 | units.kg)
        self.assertEquals(x[0].child1.mass, 2.0 | units.kg)
        self.assertEquals(x[0].child2.mass, 3.0 | units.kg)
        self.assertEquals(x[1].child1, None)
        self.assertEquals(x[1].child2, None)
    
    
        code.data[1][1] = 3
        code.data[1][2] = 2
    
        self.assertEquals(x[0].child1, x[1])
        self.assertEquals(x[0].child1.child1.mass, 4.0 | units.kg)
        self.assertEquals(x[0].child1.child2.mass, 3.0 | units.kg)

    def test6(self):
        class Code(object):
            def __init__(self):
                self.data = []
                self.number_of_particles = 0
                
            def get_number_of_particles(self):
                return  self.number_of_particles
                
            def get_mass(self,index):
                data_to_return = [self.data[i][0] for i in index]
                return units.kg(data_to_return)
            
            def get_children(self,index):
                return [(self.data[i][1]) for i in index], [(self.data[i][2]) for i in index]
            
            def new_particle(self, mass):
                mass = mass.value_in(units.kg)
                self.data = [[x,-1,-1] for x in mass]
                self.number_of_particles = len(self.data)
                return [i for i in range(len(mass))]
            
        code = Code()
    
        children_getter = ParticleGetAttributesMethod(
                    code.get_children,
                    ('child1', 'child2',)
        )
        children_getter.index_output_attributes = set(['child1','child2'])
     
        storage = InCodeAttributeStorage(
            code,
            NewParticleMethod(code.new_particle,("mass",)),
            None,
            code.get_number_of_particles,
            [],
            [
                ParticleGetAttributesMethod(code.get_mass,("mass",)),
                children_getter
            ],
            name_of_the_index = "index"
        )
        
        
        code_particles = Particles(storage  = storage)
    
        memory_particles = Particles(keys = 100 * (1 + numpy.arange(10)) )
        memory_particles.mass = range(10) | units.kg
    
        code_particles.add_particles(memory_particles)
    
        self.assertEquals(len(code_particles), 10)
    
        code.data[0][1] = 1
        code.data[0][2] = 2
        code.data[1][1] = 3
        code.data[1][2] = 4
    
        self.assertEquals(code_particles[0].child1, code_particles[1])
        self.assertEquals(code_particles[0].child1.mass, 1.0 | units.kg)
        self.assertEquals(code_particles[0].child2.mass, 2.0 | units.kg)
        self.assertEquals(code_particles[0].child1.key, 200)
        self.assertEquals(code_particles[0].child2.key, 300)
        self.assertEquals(code_particles[0].child1.child1.mass, 3.0 | units.kg)
        self.assertEquals(code_particles[0].child1.child2.mass, 4.0 | units.kg)
    
        channel = code_particles.new_channel_to(memory_particles)
        channel.copy()
    
        self.assertEquals(memory_particles[0].child1, memory_particles[1])
        self.assertEquals(memory_particles[0].child1.mass, 1.0 | units.kg)
        self.assertEquals(memory_particles[0].child2.mass, 2.0 | units.kg)
        self.assertEquals(memory_particles[0].child1.child1.mass, 3.0 | units.kg)
        self.assertEquals(memory_particles[0].child1.child2.mass, 4.0 | units.kg)

class TestGrids(amusetest.TestCase):
    
    def test1(self):
        class Code(object):
            def get_range(self):
                return (1,10,2,5,3,6)
                
            def get_ijk(self,i,j,k):
                return units.m(i), units.m(j), units.m(k)
        
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [],
            [ParticleGetAttributesMethod(code.get_ijk,("i","j","k")),],
        )
        
        self.assertEquals(storage.storage_shape(), (10, 4, 4))
        self.assertEquals(storage.get_defined_attribute_names(), ["i","j","k"])
        
        values = storage.get_values_in_store((0,1,1), ("i",))
        self.assertEquals(len(values), 1)
        print values
        self.assertEquals(values[0], 1 | units.m)
        
        values = storage.get_values_in_store((0,1,1), ("k","j","i",))
        self.assertEquals(values[0], 4 | units.m)
        self.assertEquals(values[1], 3 | units.m)
        self.assertEquals(values[2], 1 | units.m)
        
    
    def test2(self):
        class Code(object):
            def get_range(self):
                return (1,10,2,5,3,6)
                
            def get_ijk(self,i,j,k):
                return units.m(i), units.m(j), units.m(k)
        
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [],
            [ParticleGetAttributesMethod(code.get_ijk,("i","j","k")),],
        )
        values = storage.get_values_in_store(numpy.s_[0:2], ("i",))
        print values
        self.assertEquals(len(values), 1)
        self.assertEquals(len(values[0]), 2)
        self.assertEquals(values[0].number.shape, (2,4,4))
        self.assertEquals(values[0][0][0][0], 1 | units.m)
        self.assertEquals(values[0][1][0][0], 2 | units.m)
        
    
    def test3(self):
        
        shape = (11,5,5)
            
        class Code(object):
            
            def __init__(self):
                self.storage = numpy.arange(shape[0]*shape[1]*shape[2]).reshape(shape)
                
            def get_range(self):
                return (0,shape[0]-1,0,shape[1]-1,0,shape[2]-1)
                
            def get_a(self,i_s,j_s,k_s):
                return units.m.new_quantity(numpy.asarray([(self.storage[i][j][k]) for i,j,k in zip(i_s, j_s, k_s)]))
                
            def set_a(self, i_s, j_s, k_s, values):
                print i_s, j_s, k_s
                print "VALUES:", values
                index = 0
                for i,j,k in zip(i_s, j_s, k_s):
                    self.storage[i][j][k] = values[index].value_in(units.m)
                    index += 1
                    print index
                    
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [ParticleSetAttributesMethod(code.set_a,("a",)),],
            [ParticleGetAttributesMethod(code.get_a,("a",)),],
        )
        
        values = storage.get_values_in_store(None, ("a",))
        print values[0].value_in(units.m)
        self.assertTrue(numpy.all(values[0].value_in(units.m) == code.storage))
        #self.assertTrue(False)
        values = storage.get_values_in_store((0,0,0), ("a",))
        self.assertEquals(values[0], 0 | units.m)
        storage.set_values_in_store((0,0,0), ("a",), [11.0 | units.m,])
        values = storage.get_values_in_store((0,0,0), ("a",))
        self.assertEquals(values[0], 11.0 | units.m)
        values = storage.get_values_in_store((0,0), ("a",))
        storage.set_values_in_store((0,0), ("a",), [[11.0, 12.0, 13.0, 14.0, 15.0]| units.m,])
        print code.storage[0][0]
        self.assertTrue(numpy.all(code.storage[0][0] == [11.0, 12.0, 13.0, 14.0, 15.0]))
        

    def test4(self):
        class Code(object):
            def get_range(self, d, l):
                return (1,10,2,5,3,6)
                
            def get_ijk(self,i,j,k, d, l):
                return units.m(d), units.m(l), units.m(k)
        
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [],
            [ParticleGetAttributesMethod(code.get_ijk,("i","j","k")),],
            extra_keyword_arguments_for_getters_and_setters = {'d':1, 'l':2},
        )
        
        self.assertEquals(storage.storage_shape(), (10, 4, 4))
        self.assertEquals(storage.get_defined_attribute_names(), ["i","j","k"])
        
        values = storage.get_values_in_store((0,1,1), ("i",))
        self.assertEquals(len(values), 1)
        self.assertEquals(values[0], 1 | units.m)
        
        values = storage.get_values_in_store((0,1,1), ("k","j","i",))
        self.assertEquals(values[0], 4 | units.m)
        self.assertEquals(values[1], 2 | units.m)
        self.assertEquals(values[2], 1 | units.m)
    
    

    def test5(self):
        class Code(object):
            def get_range(self):
                return (1,10,2,5,3,6)
                
            def get_ijk(self,i,j,k):
                return units.m(i), units.m(j), units.m(k)
        
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [],
            [ParticleGetAttributesMethod(code.get_ijk,("i","j","k")),],
        )
        
        self.assertEquals(storage.storage_shape(), (10, 4, 4))
        self.assertEquals(storage.get_defined_attribute_names(), ["i","j","k"])
        
        values = storage.get_values_in_store(None, ("i",))
        self.assertEquals(len(values), 1)
        print values
        self.assertEquals(values[0].number.ndim, 3)
    
    
