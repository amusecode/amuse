import numpy
import time

from amuse.test import amusetest
from amuse.datamodel.memory_storage import InMemoryAttributeStorageUseDictionaryForKeySet
from amuse.datamodel.memory_storage import InMemoryAttributeStorageUseSortedKeys
from amuse.datamodel.memory_storage import get_in_memory_attribute_storage_factory
from amuse.datamodel.memory_storage import InMemoryVectorQuantityAttribute
from amuse.datamodel.incode_storage import *


class _AbstractTestInMemoryAttributeStorage:

    def new_inmemory_storage(self, is_with_units = True):
        raise NotImplementedError()

    def test1(self):
        keys = 4, 5, 6
        attributes = "a", "b"    
        values = [
            units.m.new_quantity([1.0, 2.0, 3.0]), 
            units.g.new_quantity([4.0, 5.0, 6.0])
        ]

        instance = self.new_inmemory_storage()

        indices = instance.add_particles_to_store(keys, attributes, values)

        self.assertEqual(len(indices), len(keys))
        self.assertEqual(2.0 | units.m, instance.get_value_in_store(indices[1], "a"))
        self.assertEqual(2.0 | units.m, instance.get_value_of(indices[1], "a"))


    def test2(self):
        keys =  4, 5, 6
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0, 3.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0, 6.0]))
        ]

        instance = self.new_inmemory_storage()
        indices_stored = instance.add_particles_to_store(keys, attributes, values)

        indices = instance.get_indices_of([
            keys[2], 
            keys[1], 
            keys[0]
        ])

        for index, wanted in zip(indices, [indices_stored[2], indices_stored[1], indices_stored[0]]):
            self.assertEqual(index, wanted)

    def ntest2b(self):
        keys =  4, 5, 6
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0, 3.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0, 6.0]))
        ]

        instance = self.new_inmemory_storage()
        print(instance.add_particles_to_store(keys, attributes, values))
        indices = instance.get_indices_of([
            keys[2],
            keys[0]
        ])
        all_values = instance.get_values_in_store(
            indices,
            ["b", "a"]
        )

        self.assertEqual(all_values[0][0], 6.0 | units.g)
        self.assertEqual(all_values[0][1], 4.0 | units.g)
        self.assertEqual(all_values[1][0], 0.003 | units.km)
        self.assertEqual(all_values[1][1], 0.001 | units.km)


    def test3(self):
        keys =  4, 5, 6
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0, 3.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0, 6.0]))
        ]

        instance = self.new_inmemory_storage()
        instance.add_particles_to_store(keys, attributes, values)
        self.assertEqual(values[0][0], 1.0 | units.m)

        indices = instance.get_indices_of([
            keys[0],
            keys[2]
        ])
        instance.set_values_in_store(
            indices,
            ["b", "a"],
            [
                units.kg.new_quantity(numpy.array([9.0, 11.0])),  
                units.km.new_quantity(numpy.array([1.0, 2.0]))
            ]
        )
        values = instance.get_values_in_store(None, ["a", "b"])
        self.assertEqual(values[0][0], 1000.0 | units.m)
        self.assertEqual(values[0][2], 2000.0 | units.m)
        self.assertEqual(values[1][0], 9.0 | units.kg)
        self.assertEqual(values[1][2], 11.0 | units.kg)


    def test4(self):
        keys = 10, 5, 6, 7
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0, 3.0, 4.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0, 6.0, 7.0]))
        ]

        instance = self.new_inmemory_storage()
        instance.add_particles_to_store(keys, attributes, values)

        indices = instance.get_indices_of([
            keys[2],
            keys[0]
        ])
        print(indices)
        instance.remove_particles_from_store(indices)

        self.assertEqual(len(instance), 2)

        indices = instance.get_indices_of([
            keys[1],
            keys[3]
        ])
        all_values = instance.get_values_in_store(
            indices,
            ["a", "b"]
        )
        self.assertEqual(all_values[0][0], 0.002 | units.km)
        self.assertEqual(all_values[0][1], 0.004 | units.km)
        self.assertEqual(all_values[1][0], 5.0 | units.g)
        self.assertEqual(all_values[1][1], 7.0 | units.g)


    def test5(self):
        keys = 10, 5, 6, 7
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0, 3.0, 4.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0, 6.0, 7.0]))
        ]

        instance = self.new_inmemory_storage()
        indices1 = instance.add_particles_to_store(keys, attributes, values)

        self.assertEqual(len(instance), 4)

        keys = 20, 21, 22, 23
        indices2 = instance.add_particles_to_store(keys, attributes, values)

        self.assertEqual(len(instance), 8)

        indices = instance.get_indices_of([5, 21])
        all_values = instance.get_values_in_store(indices, ["a"])

        self.assertEqual(all_values[0][0], 2.0 | units.m)
        self.assertEqual(all_values[0][1], 2.0 | units.m)

        indices = instance.get_indices_of([10, 20])
        instance.remove_particles_from_store(indices)

        self.assertEqual(len(instance), 6)

        indices = instance.get_indices_of([5, 21])
        all_values = instance.get_values_in_store(indices, ["a"])

        self.assertEqual(all_values[0][0], 2.0 | units.m)
        self.assertEqual(all_values[0][1], 2.0 | units.m)

    def test6(self):
        keys = 10, 5, 6, 7
        attributes = "a", "b"
        values = [
            numpy.array([1.0, 2.0, 3.0, 4.0]), 
            numpy.array([4.0, 5.0, 6.0, 7.0])
        ]

        instance = self.new_inmemory_storage(False)
        instance.add_particles_to_store(keys, attributes, values)

        self.assertEqual(len(instance), 4)

        particles = [4, 5, 6, 7]
        keys = 20, 21, 22, 23
        instance.add_particles_to_store(keys, attributes, values)

        self.assertEqual(len(instance), 8)

        indices = instance.get_indices_of([5, 21])
        all_values = instance.get_values_in_store(indices, ["a"])

        self.assertEqual(all_values[0][0], 2.0)
        self.assertEqual(all_values[0][1], 2.0)

        indices = instance.get_indices_of([10, 20])
        instance.remove_particles_from_store(indices)

        self.assertEqual(len(instance), 6)

        indices = instance.get_indices_of([5, 21])
        all_values = instance.get_values_in_store(indices, ["a"])

        self.assertEqual(all_values[0][0], 2.0)
        self.assertEqual(all_values[0][1], 2.0)
        instance.set_values_in_store(indices, ["a", "b"], [[4.0, 5.0], [77, 88]])
        all_values = instance.get_values_in_store(indices, ["a"])[0]
        self.assertEqual(all_values[0], 4.0)
        self.assertEqual(all_values[1], 5.0)

        if hasattr(instance, "copy"):
            instance2 = instance.copy()
            indices = instance2.get_indices_of([5, 21])
            instance2.set_values_in_store(indices, ["a"], [[3.0, 1.0]])
            all_values = instance2.get_values_in_store(indices, ["a"])[0]
            self.assertEqual(all_values[0], 3.0)
            self.assertEqual(all_values[1], 1.0)
            all_values = instance.get_values_in_store(indices, ["a"])[0]
            self.assertEqual(all_values[0], 4.0)
            self.assertEqual(all_values[1], 5.0)

    def test7(self):

        keys = 10, 5, 6, 7
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0, 3.0, 4.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0, 6.0, 7.0]))
        ]
        instance = self.new_inmemory_storage()
        indices = instance.add_particles_to_store(keys, attributes, values)
        self.assertTrue(instance.has_key_in_store(5))
        self.assertFalse(instance.has_key_in_store(1))
        self.assertFalse(instance.has_key_in_store(8))
        self.assertFalse(instance.has_key_in_store(11))

    def test8(self):

        keys = 10, 5, 6, 7
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0, 3.0, 4.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0, 6.0, 7.0]))
        ]

        instance = self.new_inmemory_storage()
        indices = instance.add_particles_to_store(keys, attributes, values)

        self.assertEqual(instance.get_all_indices_in_store(), indices)
        self.assertEqual(instance.get_indices_of(keys), indices) 

    def test9(self):

        keys = 10, 5, 6, 7
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0, 3.0, 4.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0, 6.0, 7.0]))
        ]     
        instance = self.new_inmemory_storage()
        indices = instance.add_particles_to_store(keys, attributes, values)

        self.assertRaises(Exception, instance.get_indices_of, [1], expected_message = "Key not found in storage: 1")
        self.assertRaises(Exception, instance.get_indices_of, [8], expected_message = "Key not found in storage: 8")
        self.assertRaises(Exception, instance.get_indices_of, [11], expected_message = "Key not found in storage: 11")
        self.assertRaises(Exception, instance.get_indices_of, [30], expected_message = "Key not found in storage: 30")
        self.assertRaises(Exception, instance.get_indices_of, [0], expected_message = "Key not found in storage: 0")

    def test10(self):

        keys = 10, 5, 6, 7
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0, 3.0, 4.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0, 6.0, 7.0]))
        ]     
        instance = self.new_inmemory_storage()
        indices = instance.add_particles_to_store(keys, attributes, values)

        self.assertRaises(Exception, instance.get_indices_of, [1, 8, 11, 0], expected_message = "Keys not found in storage: [ 1  8 11  0]")
        self.assertRaises(Exception, instance.get_indices_of, [5, 1, 6, 8, 7, 11, 10, 0], expected_message = "Keys not found in storage: [ 1  8 11  0]")
        self.assertRaises(Exception, instance.get_indices_of, [5, 1, 6], expected_message = "Key not found in storage: 1")
        self.assertRaises(Exception, instance.get_indices_of, [1, 5, 6], expected_message = "Key not found in storage: 1")
        self.assertRaises(Exception, instance.get_indices_of, [5, 6, 1], expected_message = "Key not found in storage: 1")


    def test11(self):
        keys = 10, 5
        attributes = "a", "b"
        values = [
            units.m.new_quantity(numpy.array([1.0, 2.0])), 
            units.g.new_quantity(numpy.array([4.0, 5.0]))
        ]

        instance = self.new_inmemory_storage()
        instance.add_particles_to_store(keys, attributes, values)

        self.assertEqual(len(instance), 2)

        keys = 7, 6
        instance.add_particles_to_store(keys, attributes, values)

        self.assertEqual(len(instance), 4)

        indices = instance.get_all_indices_in_store()[0]
        instance.remove_particles_from_store([indices])

        self.assertEqual(len(instance), 3)

        self.assertEqual(instance.get_all_keys_in_store(), [5, 7, 6])

        indices = instance.get_all_indices_in_store()[1]
        instance.remove_particles_from_store([indices])
        self.assertEqual(instance.get_all_keys_in_store(), [5, 6])

class TestSortedKeysInMemoryAttributeStorage(amusetest.TestCase, _AbstractTestInMemoryAttributeStorage):

    def new_inmemory_storage(self, is_with_units = True):
        return InMemoryAttributeStorageUseSortedKeys()

class TestDictionaryKeysInMemoryAttributeStorage(amusetest.TestCase, _AbstractTestInMemoryAttributeStorage):

    def new_inmemory_storage(self, is_with_units = True):
        return InMemoryAttributeStorageUseDictionaryForKeySet()

class _Code(object):
    def __init__(self, is_with_units):
        self.data = {}
        self.index = 100
        self.is_with_units = is_with_units

    def get_number_of_particles(self):
        return len(self.data)

    def get_ab(self, index):
        a = [(self.data[i]['a']) for i in index]
        b = [(self.data[i]['b']) for i in index]
        if self.is_with_units:
            return [units.m(a),  units.g(b)]
        else:
            return [numpy.asarray(a), numpy.asarray(b)]

    def set_ab(self, index, a, b):
        if self.is_with_units:
            a = a.value_in(units.m)
            b = b.value_in(units.g)
        for j, i in enumerate(index):
            self.data[i]['a'] = a[j]
            self.data[i]['b'] = b[j]

    def new_particle(self, a, b):
        if self.is_with_units:
            a = a.value_in(units.m)
            b = b.value_in(units.g)
        result = []
        for i in range(len(a)):
            self.data[self.index] = {}
            self.data[self.index]['a'] = a[i]
            self.data[self.index]['b'] = b[i]
            result.append(self.index)
            self.index += 1
        return result

    def delete_particle(self, index):
        for i in index:
            del self.data[i]

class TestInCodeAttributeStorage(amusetest.TestCase, _AbstractTestInMemoryAttributeStorage):

    def new_inmemory_storage(self, is_with_units = True):
        code = _Code(is_with_units)
        return InCodeAttributeStorage(
            code,
            NewParticleMethod(code.new_particle, ("a", "b")),
            code.delete_particle,
            code.get_number_of_particles,
            [ParticleSetAttributesMethod(code.set_ab, ("a", "b")), ],
            [ParticleGetAttributesMethod(code.get_ab, ("a", "b")), ],
            name_of_the_index = "index"
        )




