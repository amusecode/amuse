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
        self.assertEquals(storage._get_attribute_names(), set(["x","y","z"]))
        
        self.assertFalse(code.get_position_called)
        storage._get_values([],["x","y","z"])
        self.assertFalse(code.get_position_called)
        
        storage._add_particles(
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
        
        storage._add_particles(
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
        
        self.assertEquals(storage._get_attribute_names(), set(["x","y","z", "mass"]))
        
        self.assertFalse(code.get_position_called)
        self.assertFalse(code.get_mass_called)
        x,y,mass = storage._get_values([2,3],["x","y","mass"])
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
        
        storage._add_particles(
            [1,2,3,4],
            ["mass"],
            [
                units.kg([1,2,3,4]),
            ]
        )
        
        self.assertEquals(len(storage), 4)
        
        self.assertEquals(storage._get_attribute_names(), set(["mass",]))
        
        index,mass = storage._get_values([2,3],["index_in_code","mass"])
        self.assertTrue(code.get_mass_called)
        print index, mass
        self.assertEquals(index[0], 1)
        self.assertEquals(mass[0],  2 | units.kg)
        self.assertEquals(index[1], 2)
        self.assertEquals(mass[1],  3 | units.kg)
    
    
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
        self.assertEquals(storage._get_attribute_names(), set(["i","j","k"]))
        
        values = storage._get_values((0,1,1), ("i",))
        self.assertEquals(len(values), 1)
        print values
        self.assertEquals(values[0], 1 | units.m)
        
        values = storage._get_values((0,1,1), ("k","j","i",))
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
        values = storage._get_values(numpy.s_[0:2], ("i",))
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
        
        values = storage._get_values(None, ("a",))
        print values[0].value_in(units.m)
        self.assertTrue(numpy.all(values[0].value_in(units.m) == code.storage))
        #self.assertTrue(False)
        values = storage._get_values((0,0,0), ("a",))
        self.assertEquals(values[0], 0 | units.m)
        storage._set_values((0,0,0), ("a",), [11.0 | units.m,])
        values = storage._get_values((0,0,0), ("a",))
        self.assertEquals(values[0], 11.0 | units.m)
        values = storage._get_values((0,0), ("a",))
        storage._set_values((0,0), ("a",), [[11.0, 12.0, 13.0, 14.0, 15.0]| units.m,])
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
        self.assertEquals(storage._get_attribute_names(), set(["i","j","k"]))
        
        values = storage._get_values((0,1,1), ("i",))
        self.assertEquals(len(values), 1)
        self.assertEquals(values[0], 1 | units.m)
        
        values = storage._get_values((0,1,1), ("k","j","i",))
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
        self.assertEquals(storage._get_attribute_names(), set(["i","j","k"]))
        
        values = storage._get_values(None, ("i",))
        self.assertEquals(len(values), 1)
        print values
        self.assertEquals(values[0].number.ndim, 3)
    
    
