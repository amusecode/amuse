from amuse.test import amusetest
from amuse.datamodel.incode_storage import *

import numpy
import time

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system

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
        
        self.assertEqual(len(storage), 0)
        self.assertEqual(storage.get_defined_attribute_names(), ["x","y","z"])
        
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
        
        self.assertEqual(len(storage), 4)
    
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
        
        self.assertEqual(len(storage), 4)
        
        self.assertEqual(storage.get_defined_attribute_names(), [ "mass", "x","y","z"])
        
        self.assertFalse(code.get_position_called)
        self.assertFalse(code.get_mass_called)
        indices = storage.get_indices_of([2,3])
        x,y,mass = storage.get_values_in_store(indices,["x","y","mass"])
        self.assertTrue(code.get_position_called)
        self.assertTrue(code.get_mass_called)
        self.assertEqual(x[1], 3 | units.m)
        self.assertEqual(mass[1], 15 | units.kg)
        

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
        
        self.assertEqual(len(storage), 4)
        
        self.assertEqual(storage.get_defined_attribute_names(), ["mass",])
        
        indices = storage.get_indices_of([2,3])
        index,mass = storage.get_values_in_store(indices,["index_in_code","mass"])
        self.assertTrue(code.get_mass_called)
        self.assertEqual(index[0], 1)
        self.assertEqual(mass[0],  2 | units.kg)
        self.assertEqual(index[1], 2)
        self.assertEqual(mass[1],  3 | units.kg)
    
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
        
        self.assertEqual(len(storage), 4)
        
        storage._remove_indices([1,2,])
        code.number_of_particles = 2
        indices = storage.get_indices_of([1,4])
        index,mass = storage.get_values_in_store(indices,["index_in_code","mass"])
        
        self.assertEqual(index[0], 0)
        self.assertEqual(index[1], 3)
        self.assertEqual(mass[0],  1 | units.kg)
        self.assertEqual(mass[1],  4 | units.kg)
        
        
        self.assertEqual(len(storage), 2)
        
        storage._add_indices([4,5])
        code.data = numpy.concatenate((code.data, [5, 6]))
        
        code.number_of_particles = 4
        self.assertEqual(len(storage), 4)
        
        indices = storage.get_indices_of(storage.particle_keys)
        mass, = storage.get_values_in_store(indices,["mass"])
        
        self.assertEqual(mass[0],  1 | units.kg)
        self.assertEqual(mass[1],  4 | units.kg)
        self.assertEqual(mass[2],  5 | units.kg)
        self.assertEqual(mass[3],  6 | units.kg)
        
        storage._remove_indices([4,])
        code.number_of_particles = 3
        self.assertEqual(len(storage), 3)
        
        indices = storage.get_indices_of(storage.particle_keys)
        mass, = storage.get_values_in_store(indices,["mass"])
        
        self.assertEqual(mass[0],  1 | units.kg)
        self.assertEqual(mass[1],  4 | units.kg)
        self.assertEqual(mass[2],  6 | units.kg)
    
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
        
        self.assertEqual(len(storage), 4)
    
        indices = storage.get_indices_of([100,400])
        mass = storage.get_values_in_store(indices,["mass",])[0]
        self.assertEqual(mass[0], 1.0 | units.kg)
        self.assertEqual(mass[1], 4.0 | units.kg)
    
        code.data[0][1] = 1
        code.data[0][2] = 2
    
        indices = storage.get_indices_of([100])
        child1,child2 = storage.get_values_in_store(indices,['child1', 'child2'])
    
        self.assertEqual(child1[0].number, 200)
        self.assertEqual(child2[0].number, 300)
      
    def test7(self):
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
                return data_to_return
                
            def set_position(self,index,x,y,z):
                self.set_position_called = True
                pass
                
            def set_mass(self,index,mass):
                self.set_mass_called = True
                for i,j  in enumerate(index):
                    self.data[3][j] = mass[i]
                return [0 for i in range(len(index))]
                
            def new_particle(self, x, y, z, mass):
                x = x.value_in(units.m)
                y = y.value_in(units.m)
                z = z.value_in(units.m)
                mass = mass 
                self.data = [x,y,z,mass]
                return [i for i in range(len(x))]
                
        code = Code()
        storage = InCodeAttributeStorage(
            code,
            NewParticleMethod(code.new_particle,("x","y","z","mass")),
            None,
            code.get_number_of_particles,
            [
                ParticleSetAttributesMethod(code.set_position,("x","y","z")),
                ParticleSetAttributesMethod(code.set_mass,("mass",)),
            ],
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
                numpy.asarray([13.0,14.0,15,16]),
            ]
        )
        
        self.assertEqual(len(storage), 4)
        
        self.assertEqual(storage.get_defined_attribute_names(), [ "mass", "x","y","z"])
        
        self.assertFalse(code.get_position_called)
        self.assertFalse(code.get_mass_called)
        indices = storage.get_indices_of([2,3])
        x,y,mass = storage.get_values_in_store(indices,["x","y","mass"])
        self.assertTrue(code.get_position_called)
        self.assertTrue(code.get_mass_called)
        self.assertEqual(x[1], 3 | units.m)
        self.assertEqual(mass[1], 15 )
        self.assertEqual(mass[0], 14 )
        storage.set_values_in_store(indices,["x","y", "z", "mass"], [[10,11] | units.m , [12,14] | units.m, [12,14] | units.m, [40.0, 50.0]])
        x,y,mass = storage.get_values_in_store(indices,["x","y","mass"])
        self.assertEqual(mass[1], 50 )
        self.assertEqual(mass[0], 40 )
        
    
        
        

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
        
        self.assertEqual(storage.storage_shape(), (10, 4, 4))
        self.assertEqual(storage.get_defined_attribute_names(), ["i","j","k"])
        
        values = storage.get_values_in_store((0,1,1), ("i",))
        self.assertEqual(len(values), 1)
        self.assertEqual(values[0], 1 | units.m)
        
        values = storage.get_values_in_store((0,1,1), ("k","j","i",))
        self.assertEqual(values[0], 4 | units.m)
        self.assertEqual(values[1], 3 | units.m)
        self.assertEqual(values[2], 1 | units.m)
    
    
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
        self.assertEqual(len(values), 1)
        self.assertEqual(len(values[0]), 2)
        self.assertEqual(values[0].number.shape, (2,4,4))
        self.assertEqual(values[0][0][0][0], 1 | units.m)
        self.assertEqual(values[0][1][0][0], 2 | units.m)
    
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
                #~ print i_s, j_s, k_s
                #~ print "VALUES:", values
                index = 0
                for i,j,k in zip(i_s, j_s, k_s):
                    self.storage[i][j][k] = values[index].value_in(units.m)
                    index += 1
                    #~ print index
                    
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [ParticleSetAttributesMethod(code.set_a,("a",)),],
            [ParticleGetAttributesMethod(code.get_a,("a",)),],
        )
        
        values = storage.get_values_in_store(None, ("a",))
        self.assertTrue(numpy.all(values[0].value_in(units.m) == code.storage))
        #self.assertTrue(False)
        values = storage.get_values_in_store((0,0,0), ("a",))
        self.assertEqual(values[0], 0 | units.m)
        storage.set_values_in_store((0,0,0), ("a",), [11.0 | units.m,])
        values = storage.get_values_in_store((0,0,0), ("a",))
        self.assertEqual(values[0], 11.0 | units.m)
        values = storage.get_values_in_store((0,0), ("a",))
        storage.set_values_in_store((0,0), ("a",), [[11.0, 12.0, 13.0, 14.0, 15.0]| units.m,])
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
        
        self.assertEqual(storage.storage_shape(), (10, 4, 4))
        self.assertEqual(storage.get_defined_attribute_names(), ["i","j","k"])
        
        values = storage.get_values_in_store((0,1,1), ("i",))
        self.assertEqual(len(values), 1)
        self.assertEqual(values[0], 1 | units.m)
        
        values = storage.get_values_in_store((0,1,1), ("k","j","i",))
        self.assertEqual(values[0], 4 | units.m)
        self.assertEqual(values[1], 2 | units.m)
        self.assertEqual(values[2], 1 | units.m)
    
    
        

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
        
        self.assertEqual(storage.storage_shape(), (10, 4, 4))
        self.assertEqual(storage.get_defined_attribute_names(), ["i","j","k"])
        
        values = storage.get_values_in_store(None, ("i",))
        self.assertEqual(len(values), 1)
        self.assertEqual(values[0].number.ndim, 3)
        
    
    def test6(self):
        
        shape = (11,5,5)
            
        class Code(object):
            
            def __init__(self):
                self.storage = numpy.arange(shape[0]*shape[1]*shape[2]).reshape(shape)
                
            def get_range(self):
                return (0,shape[0]-1,0,shape[1]-1,0,shape[2]-1)
                
            def get_a(self,i_s,j_s,k_s):
                return numpy.asarray([(self.storage[i][j][k]) for i,j,k in zip(i_s, j_s, k_s)])
                
            def set_a(self, i_s, j_s, k_s, values):
                #~ print i_s, j_s, k_s
                #~ print "VALUES:", values
                index = 0
                for i,j,k in zip(i_s, j_s, k_s):
                    self.storage[i][j][k] = values[index]
                    index += 1
                    #~ print index
                    
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [ParticleSetAttributesMethod(code.set_a,("a",)),],
            [ParticleGetAttributesMethod(code.get_a,("a",)),],
        )
        
        values = storage.get_values_in_store(None, ("a",))
        self.assertTrue(numpy.all(values[0] == code.storage))
        values = storage.get_values_in_store((0,0,0), ("a",))
        self.assertEqual(values[0], 0)
        storage.set_values_in_store((0,0,0), ("a",), [11.0,])
        values = storage.get_values_in_store((0,0,0), ("a",))
        self.assertEqual(values[0], 11.0)
        values = storage.get_values_in_store((0,0), ("a",))[0]
        self.assertTrue(numpy.all(values == [11.0, 1.0, 2.0, 3.0, 4.0]))
        storage.set_values_in_store((0,0), ("a",), [[11.0, 12.0, 13.0, 14.0, 15.0],])
        self.assertTrue(numpy.all(code.storage[0][0] == [11.0, 12.0, 13.0, 14.0, 15.0]))
    
    def test7(self):
        
        shape = (11,5,5)
            
        class Code(object):
            
            def __init__(self):
                self.storage = numpy.arange(shape[0]*shape[1]*shape[2]).reshape(shape)
                
            def get_range(self):
                return (0,shape[0]-1,0,shape[1]-1,0,shape[2]-1)
                
            def get_a(self,i_s,j_s,k_s):
                return numpy.asarray([(self.storage[i][j][k]) for i,j,k in zip(i_s, j_s, k_s)])
                
            def set_a(self, i_s, j_s, k_s, values):
                index = 0
                for i,j,k in zip(i_s, j_s, k_s):
                    self.storage[i][j][k] = values[index]
                    index += 1
                    
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [ParticleSetAttributesMethod(code.set_a,("a",)),],
            [ParticleGetAttributesMethod(code.get_a,("a",)),],
        )
        
        values = storage.get_values_in_store((), ())
        self.assertTrue(values==[])
        values = storage.get_values_in_store((0,0,1,), ("a",))
        self.assertTrue(values[0]==1)

    def test8(self):
        class Code(object):
            def __init__(self):
                self.storage = 1. | units.m

            def get_range(self):
                return ()
                
            def get_a(self):
                return self.storage
                
            def set_a(self, value):
                self.storage=value
                

        
        code = Code()
        
        storage = InCodeGridAttributeStorage(
            code,
            code.get_range,
            [ParticleSetAttributesMethod(code.set_a,("a",)),],
            [ParticleGetAttributesMethod(code.get_a,("a",)),],
        )
        
        self.assertEqual(storage.storage_shape(), ())
        self.assertEqual(storage.get_defined_attribute_names(), ['a'])
        
        values = storage.get_values_in_store((), ("a",))
        self.assertEqual(len(values), 1)
        print(values,"<")
        self.assertEqual(values[0], 1 | units.m)
        
    
