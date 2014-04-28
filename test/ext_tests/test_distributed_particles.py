from amuse.datamodel import AbstractSet
from amuse.test.amusetest import TestWithMPI
from amuse.ext import concurrent
from amuse import datamodel

from amuse.rfi.core import *
from amuse.ic.plummer import new_plummer_sphere
from amuse.units import nbody_system
from amuse.units import units

from amuse.community import *
from amuse.support.interface import InCodeComponentImplementation

from collections import namedtuple
import numpy
import pickle
try:
    from mpi4py import MPI
except ImportError:
    MPI = None
    
def dump_and_encode(x):
  return pickle.dumps(x,0) # -1 does not work with sockets channel
def decode_and_load(x):
  return pickle.loads(x.encode("latin-1"))

class DistributedParticlesInterface(PythonCodeInterface):
    
    def __init__(self, implementation_factory, **options):
        PythonCodeInterface.__init__(self, implementation_factory = implementation_factory, **options)
        
    @legacy_function
    def get_length():
        function = LegacyFunctionSpecification()  
        function.addParameter('reference_in', dtype='int64', direction=function.IN)
        function.addParameter('output', dtype='int64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def do_init():
        function = LegacyFunctionSpecification()  
        function.addParameter('input', dtype='int64', direction=function.IN)
        function.addParameter('reference_out', dtype='int64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def do_set_attribute():
        function = LegacyFunctionSpecification()  
        function.addParameter('reference_in', dtype='int64', direction=function.IN)
        function.addParameter('name_of_the_attribute', dtype='string', direction=function.IN)
        function.addParameter('pickled_value', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def do_get_attribute():
        function = LegacyFunctionSpecification()  
        function.addParameter('reference_in', dtype='int64', direction=function.IN)
        function.addParameter('name_of_the_attribute', dtype='string', direction=function.IN)
        function.addParameter('pickled_value', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def do_getitem():
        function = LegacyFunctionSpecification()  
        function.addParameter('reference_in', dtype='int64', direction=function.IN)
        function.addParameter('pickled_index', dtype='string', direction=function.IN)
        function.addParameter('is_particle', dtype='bool', direction=function.OUT)
        function.addParameter('reference', dtype='int64', direction=function.OUT)
        function.result_type = 'int32'
        return function

    @legacy_function
    def set_from_generator():
        function = LegacyFunctionSpecification()  
        function.addParameter('reference_in', dtype='int64', direction=function.IN)
        function.addParameter('pickled_generator', dtype='string', direction=function.IN)
        function.addParameter('pickled_args', dtype='string', direction=function.IN)
        function.addParameter('pickled_kwargs', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function        

    @legacy_function
    def select_array():
        function = LegacyFunctionSpecification()  
        function.addParameter('reference_in', dtype='int64', direction=function.IN)
        function.addParameter('pickled_function', dtype='string', direction=function.IN)
        function.addParameter('pickled_attributes', dtype='string', direction=function.IN)
        function.addParameter('reference_out', dtype='int64', direction=function.OUT)
        function.result_type = 'int32'
        return function
  
class DistributedParticlesCode(InCodeComponentImplementation):
    
    def __init__(self, **options):
        interface = DistributedParticlesInterface(
            implementation_factory = DistributedParticlesImplementation,
            **options
        )
        
        InCodeComponentImplementation.__init__(
            self,
            interface,
            **options
        )
        
class DistributedParticles(object):
    
    def __init__(self, size = 0, reference = 0, code = None, **options):
        if code is None:
            code = DistributedParticlesCode(**options)
            code.do_init(size)
            reference = 0
        
        object.__setattr__(self, "code", code)
        object.__setattr__(self, "reference", reference)
        
    def __len__(self):
        return self.code.get_length(self.reference)
        
    def __setattr__(self, name_of_the_attribute, value):
        self.code.do_set_attribute(
            self.reference,
            name_of_the_attribute, 
            dump_and_encode(value)
        )
        
    def __getattr__(self, name_of_the_attribute):
        return decode_and_load(self.code.do_get_attribute(self.reference,name_of_the_attribute)[0])
        
    def __getitem__(self, index):
        is_particle, reference = self.code.do_getitem(
            self.reference,
            dump_and_encode(index)
        )
        if is_particle:
            return DistributedParticle(
                code = self.code,
                reference = reference
            )
        else:
            return DistributedParticles(
                code = self.code,
                reference = reference
            )

    def set_from_generator(self,func, args=(),kwargs={}):        
        self.code.set_from_generator(self.reference, dump_and_encode(func),
                dump_and_encode(args),dump_and_encode(kwargs))

    def select_array(self, function, attributes):
        
        reference=self.code.select_array(self.reference,
                                         dump_and_encode(function),
                                         dump_and_encode(attributes)
        )

        return DistributedParticles(
            code = self.code,
            reference = reference
        )
            
class DistributedParticle(object):
    
    def __init__(self, reference, code, **options):
        object.__setattr__(self, "code", code)
        object.__setattr__(self, "reference", reference)
        
    def __setattr__(self, name_of_the_attribute, value):
        self.code.do_set_attribute(
            self.reference,
            name_of_the_attribute, 
            dump_and_encode(value)
        )
        
    def __getattr__(self, name_of_the_attribute):
        return decode_and_load(self.code.do_get_attribute(self.reference,name_of_the_attribute)[0])
        
ReferencedParticles = namedtuple('ReferencedParticles', ['particles', 'local_offset', 'local_size', 'global_size'])

class DistributedParticlesImplementation(object):
    EmptyReference = ReferencedParticles(
        None,
        0,
        0,
        0
    )
    
    def do_init(self, size, reference_out):
        
        self.mpi_comm = MPI.COMM_WORLD
        self.rank =  MPI.COMM_WORLD.Get_rank()
        self.number_of_processes = MPI.COMM_WORLD.Get_size()
        local_size  = size / self.number_of_processes 
        left_over = size - (local_size * self.number_of_processes)
        if self.rank == 0:
            local_size += left_over
            local_offset = 0
        else:
            local_offset = (local_size * self.rank) + left_over
        self.reference_counter = 0
        self.references_to_particles = {}
        real_particles = datamodel.Particles(local_size)
        self.references_to_particles[self.reference_counter] = ReferencedParticles(
            real_particles,
            local_offset,
            local_size,
            size
        )
        reference_out.value = self.reference_counter
        self.reference_counter += 1
        return 0
    
    def get_length(self, reference, len_out):
        particles_len = self.references_to_particles[reference].local_size
        
        input = numpy.zeros(1,  dtype='int64')
        output = numpy.zeros(1,  dtype='int64')
        
        input[0] = particles_len
            
        self.mpi_comm.Reduce(
            [input, MPI.INTEGER], 
            [output, MPI.INTEGER],
            op=MPI.SUM, 
            root=0
        )
        len_out.value = output[0]
        return 0
        
    def do_set_attribute(self, reference, name_of_the_attribute, pickled_value):
        real_particles,local_offset, local_size, global_size = self.references_to_particles[reference]
        value = decode_and_load(pickled_value)
        if not real_particles is None:
            if (is_quantity(value) and not value.is_scalar()) or hasattr(value, '__iter__'):
                setattr(real_particles, name_of_the_attribute, value[local_offset:local_offset + local_size])
            else:
                setattr(real_particles, name_of_the_attribute, value)
        return 0
        
    def _merge_results(self, quantities): 
        values = []
        value = None
        quantities = [x for x in quantities if x[-1] is not None]
        
        if len(quantities) == 0:
            return None
            
        if len(quantities) == 1:
            local_size, local_offset, quantity = quantities[0]
            if local_size == 0:
                return quantity
        
        total_size = 0
        for local_size, local_offset, quantity in quantities:
            total_size += local_size
            
        first_quantity = quantities[0][-1]
        result = numpy.zeros(total_size, dtype = first_quantity.dtype)
        if is_quantity(first_quantity):
            result = result | first_quantity.unit
            
        for local_size, local_offset, quantity in quantities:
            result[local_offset:local_offset+local_size] = quantity
        return result
        
    def do_get_attribute(self, reference, name_of_the_attribute, output):
        real_particles,local_offset, local_size, global_size = self.references_to_particles[reference]
        
        if not real_particles is None:
            quantity = [local_size, local_offset, getattr(real_particles, name_of_the_attribute)]
        else:
            quantity = [0, 0, None]
        
        quantities = self.mpi_comm.gather(quantity, root = 0)
        
        if self.rank  == 0:
           value = self._merge_results(quantities)
        else:
            value = None
        output.value = dump_and_encode(value)
        return 0
        
    def do_getitem(self, reference_in, pickled_index, is_particle_out, reference_out):
        real_particles, local_offset, local_size, global_size  = self.references_to_particles[reference_in]
        index = decode_and_load(pickled_index)
        self.references_to_particles[self.reference_counter] = self.EmptyReference
        if isinstance(index, int):
            if index >= local_offset and index < local_offset + local_size:
                output_particles = real_particles[index - local_offset]
                self.references_to_particles[self.reference_counter] = ReferencedParticles(
                    output_particles,
                    0,
                    0,
                    0
                )
            is_particle_out.value = True
        elif isinstance(index, slice):
            start, stop, step = index.indices(global_size)
            start -= local_offset
            stop -= local_offset
            start = min(max(0, start), local_size)
            stop = min(max(0, stop), local_size)
            output_particles = real_particles[start:stop]
            input = numpy.zeros(1,  dtype='int64')
            output = numpy.zeros(self.number_of_processes,  dtype='int64')
            
            input[0] = len(output_particles)
                
            self.mpi_comm.Allgather(
                [input, MPI.INTEGER], 
                [output, MPI.INTEGER]
            )
            total_size = 0
            local_offset = 0
            for i,current_size in enumerate(output):
                if i < self.rank:
                    local_offset += current_size
                total_size += current_size
            referenced_particles = ReferencedParticles(
                output_particles,
                local_offset,
                len(output_particles),
                total_size
            )
            self.references_to_particles[self.reference_counter] = referenced_particles
            
            is_particle_out.value = False # check, a slice could still result in 1 particle!!!
        else:
            raise Exception("need to parse index and do smart things here!!")
        
        reference_out.value = self.reference_counter
        self.reference_counter += 1
        return 0
        
    def set_from_generator(self,reference, pickled_generator, pickled_args=None,pickled_kwargs=None):
        generator=decode_and_load(pickled_generator)
        args=()
        if pickled_args is not None:
          args=decode_and_load(pickled_args)
        kwargs={}
        if pickled_kwargs is not None:
          kwargs=decode_and_load(pickled_kwargs)

        real_particles,local_offset, local_size, global_size = self.references_to_particles[reference]
        particles=generator(local_offset,local_offset+local_size,*args,**kwargs)
        keys = real_particles.get_all_keys_in_store()
        real_particles.remove_particles(real_particles[:])
        attributes = particles.get_attribute_names_defined_in_store()
        indices = particles.get_all_indices_in_store()
        values =  particles.get_values_in_store(indices, attributes)
        result = datamodel.Particles()
        result.add_particles_to_store(keys, attributes, values)
        real_particles.add_particles(result)
        return 0

    def select_array(self, reference_in, pickled_function, pickled_attributes, reference_out):
        real_particles, local_offset, local_size, global_size  = self.references_to_particles[reference_in]
        function = decode_and_load(pickled_function)
        attributes = decode_and_load(pickled_attributes)
        
        output_particles = real_particles.select_array(function, attributes)
        
        input = numpy.zeros(1,  dtype='int64')
        output = numpy.zeros(self.number_of_processes,  dtype='int64')
            
        input[0] = len(output_particles)
                
        self.mpi_comm.Allgather(
                [input, MPI.INTEGER], 
                [output, MPI.INTEGER]
        )
        total_size = 0
        local_offset = 0
        for i,current_size in enumerate(output):
            if i < self.rank:
                local_offset += current_size
            total_size += current_size
        referenced_particles = ReferencedParticles(
            output_particles,
            local_offset,
            len(output_particles),
            total_size
        )
        self.references_to_particles[self.reference_counter] = referenced_particles
        reference_out.value = self.reference_counter
        self.reference_counter += 1
        return 0


def generate_set_example_function(start,end,*args,**kwargs):
    from amuse.datamodel import Particles
    p=Particles(end-start)
    p.index=range(start,end)
    return p

def select_example_function(x):
    return x > 5
    
class TestDistributedParticles(TestWithMPI):

    def setUp(self):
        if MPI is None or CodeInterface(must_start_worker=False).channel_type != 'mpi':
            self.skip("test needs mpi")
        
    def test1(self):
        x = DistributedParticles(
            size = 3,
            number_of_workers = 1
        )
        self.assertEquals(len(x) , 3)
        x.mass = 10 | units.MSun
        self.assertEquals(x.mass, [10, 10, 10]| units.MSun)
        

    def test2(self):
        x = DistributedParticles(
            size = 3,
            number_of_workers = 1
        )
        self.assertEquals(len(x) , 3)
        x.mass = 10 | units.MSun
        y = x[0:2]
        self.assertEquals(len(y) , 2)
        self.assertEquals(y.mass, [10, 10]| units.MSun)
        

    def test3(self):
        x = DistributedParticles(
            size = 3,
            number_of_workers = 1
        )
        self.assertEquals(len(x) , 3)
        x.mass = [1,2,3]| units.MSun
        y = x[0:2]
        z = x[1:]
        self.assertEquals(len(y) , 2)
        self.assertEquals(y.mass, [1,2]| units.MSun)
        self.assertEquals(len(z) , 2)
        self.assertEquals(z.mass, [2,3]| units.MSun)
        z.mass = [4,5] | units.MSun
        self.assertEquals(y.mass, [1,4]| units.MSun)
        self.assertEquals(z.mass, [4,5]| units.MSun)
        self.assertEquals(x.mass, [1,4,5]| units.MSun)
        
        
    def test4(self):
        x = DistributedParticles(
            size = 3,
            number_of_workers = 1
        )
        x.mass = [1,2,3]| units.MSun
        y = x[1]
        self.assertEquals(y.mass, 2 | units.MSun)
        y.mass = 10 | units.MSun
        self.assertEquals(y.mass, 10 | units.MSun)
        self.assertEquals(x.mass, [1,10,3]| units.MSun)
        
    def test5(self):
        x = DistributedParticles(
            size = 8,
            number_of_workers = 2
        )
        self.assertEquals(len(x) , 8)
        x.mass = [1,2,3,4,5,6,7,8] | units.MSun
        for index in range(len(x)):
            self.assertEquals(x[index].mass, (index+1)| units.MSun)
            
    
    def test6(self):
        x = DistributedParticles(
            size = 9,
            number_of_workers = 2
        )
        self.assertEquals(len(x) , 9)
        x.mass = [1,2,3,4,5,6,7,8,9] | units.MSun
        for index in range(len(x)):
            self.assertEquals(x[index].mass, (index+1)| units.MSun)
        
    def test7(self):
        x = DistributedParticles(
            size = 8,
            number_of_workers = 2
        )
        self.assertEquals(len(x) , 8)
        x.mass = [1,2,3,4,5,6,7,8] | units.MSun
        self.assertEquals(len(x[3:7]), 4)
        x[3:7].mass = [10,11,12,13] | units.MSun
        expected = [1,2,3,10,11,12,13,8]| units.MSun
        for index in range(len(x)):
            self.assertEquals(x[index].mass, expected[index] )
            
    def test8(self):
        x = DistributedParticles(
            size = 8,
            number_of_workers = 2
        )
        self.assertEquals(len(x) , 8)
        x.mass = [1,2,3,4,5,6,7,8] | units.MSun
        self.assertEquals(len(x[3:7]), 4)
        x[3:7].mass = [10,11,12,13] | units.MSun
        
        self.assertEquals(x[2:6].mass, [3,10,11,12]| units.MSun)
            
    def test9(self):
        x = DistributedParticles(
            size = 8,
            number_of_workers = 2
        )
        self.assertEquals(len(x) , 8)
        x.mass = [1,2,3,4,5,6,7,8] 
        self.assertEquals(len(x[3:7]), 4)
        x[3:7].mass = [10,11,12,13] 
        
        self.assertEquals(x[2:6].mass, [3,10,11,12])
        
    def test10(self):
        x = DistributedParticles(
            size = 40,
            number_of_workers = 4
        )
        self.assertEquals(len(x) , 40)
        x.mass = range(40) 
        self.assertEquals(len(x[15:25]), 10)
        self.assertEquals(x[15:25].mass, range(15,25))
        x[15:25].mass = range(10)
        self.assertEquals(x[15:25].mass, range(10))
            
    def test11(self):
        from test_distributed_particles import generate_set_example_function
        y=generate_set_example_function(0,10)
        
        x = DistributedParticles(
            size = 10,
            number_of_workers = 2
        )
        x.set_from_generator(generate_set_example_function)
        self.assertEqual(y.index,x.index)

    def test12(self):
        from test_distributed_particles import generate_set_example_function
        from test_distributed_particles import select_example_function
        y=generate_set_example_function(0,10)
        
        x = DistributedParticles(
            size = 10,
            number_of_workers = 2
        )
        x.set_from_generator(generate_set_example_function)
        self.assertEqual(y.index,x.index)        

        highy=y.select_array( select_example_function, ("index",))
        highx=x.select_array( select_example_function, ("index",))
        self.assertEqual(highy.index,highx.index)        
  
    def test13(self):
        x = DistributedParticles(
            size = 8,
            number_of_workers = 2
        )
        x.index= [0,10,10,10,0,0,10,0]
        x.mass = [1, 2, 3, 4,5,6, 7,8] | units.MSun
        highx=x.select_array( select_example_function, ("index",))
        self.assertEqual(highx.mass, [2,3,4,7] | units.MSun)
