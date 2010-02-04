from amuse.legacy.support.core import legacy_function, legacy_global, LegacyFunctionSpecification
from amuse.support.core import late, OrderedDictionary

from mpi4py import MPI
import numpy

class Message(object):
    
    def __init__(self, tag = -1, length= 1):
        self.tag = tag
        self.length = length
        self.ints = []
        self.doubles = []
        self.floats = []
        self.strings = []
        
    def recieve(self, comm):
        header = numpy.zeros(6,  dtype='i')
        
        comm.Bcast([header, MPI.INT], root = 0)
        
        self.tag = header[0]
        self.length = header[1]
        
        number_of_doubles = header[2]
        number_of_ints = header[3]
        number_of_floats = header[4]
        number_of_strings = header[5]
        
        self.doubles = self.recieve_doubles(comm, self.length, number_of_doubles)
        self.ints = self.recieve_ints(comm, self.length, number_of_ints)
        self.floats = self.recieve_floats(comm, self.length, number_of_floats)
        self.strings = self.recieve_strings(comm, self.length, number_of_strings)
    
    def recieve_doubles(self, comm, length, total):
        if total > 0:
            result = numpy.empty(total * length,  dtype='d')
            comm.Bcast([result, MPI.DOUBLE], root = 0)
            
            return result
        else:
            return []
            
    def recieve_ints(self, comm, length, total):
        if total > 0:
            result = numpy.empty(total * length,  dtype='i')
            comm.Bcast([result, MPI.INT], root = 0)
            
            return result
        else:
            return []
            
    def recieve_floats(self, comm, length, total):
        if total > 0:
            result = numpy.empty(total * length,  dtype='f')
            comm.Bcast([result, MPI.FLOAT], root = 0)
            
            return result
        else:
            return []
            
    def recieve_strings(self, comm, length, total):
        if total > 0:
            raise Exception("not implemented strings yet!")
        else:
            return []
        
    
    def send(self, comm):
        header = numpy.array([
            self.tag, 
            self.length, 
            len(self.doubles) / self.length, 
            len(self.ints) / self.length, 
            len(self.floats) / self.length, 
            len(self.strings) / self.length
        ], dtype='i')
        
        comm.Send([header, MPI.INT], dest=0, tag=999)
        
        
        self.send_doubles(comm, self.doubles)
        self.send_ints(comm, self.ints)
        self.send_floats(comm, self.floats)
        self.send_strings(comm, self.strings)
        
    
    def send_doubles(self, comm, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='d')
            comm.Send([buffer, MPI.DOUBLE], dest=0, tag = 999)
            
    def send_ints(self, comm, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='i')
            comm.Send([buffer, MPI.INT], dest=0, tag = 999)
            
    def send_floats(self, comm, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='f')
            comm.Send([buffer, MPI.FLOAT], dest=0, tag = 999)
            
    def send_strings(self, comm, array):
        if len(array) > 0:
            raise Exception("Sending strings not supported!")
            

def pack_array(array, length,  dtype):
    result = numpy.empty(length * len(array), dtype = dtype)
    for i in range(len(array)):
        offset = i * length
        result[offset:offset+length] = array[i]
    return result
    

def unpack_array(array, length, dtype):
    result = []
    total = len(array) / length
    for i in range(total):
        offset = i * length
        result.append(array[offset:offset+length])
    return result

class ValueHolder(object):
    
    def __init__(self, value = None):
        self.value = None


class PythonImplementation(object):
    dtype_to_message_attribute = { 
        'int32' : 'ints',
        'float64' : 'doubles',
        'float32' : 'floats',
        'string' : 'strings', 
    }
    
    def __init__(self, implementation, interface):
        self.implementation = implementation
        self.interface = interface
        self.must_run = False
        
    def start(self):
        parent = MPI.Comm.Get_parent()
        parent = MPI.Comm.Get_parent()
        
        rank = MPI.COMM_WORLD.rank
        
        self.must_run = True
        while self.must_run:
            
            message = Message()
            message.recieve(parent)
                
            result_message = Message(message.tag, message.length)
            
            if message.tag == 0:
                self.must_run = False
            else:
                self.handle_message(message, result_message)
            
            if rank == 0:
                result_message.send(parent)
        
        parent.Disconnect()
        
    def handle_message(self, input_message, output_message):
        legacy_function = self.mapping_from_tag_to_legacy_function[input_message.tag]
        specification = legacy_function.specification
        
        
        dtype_to_count = self.get_dtype_to_count(specification)
        
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            count = dtype_to_count.get(type,0)
            for x in range(count):
                getattr(output_message, attribute).append(numpy.empty(output_message.length, dtype=type))
        
        
        method = getattr(self.implementation, specification.name)
        
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            array = getattr(input_message, attribute)
            unpacked = unpack_array(array, input_message.length, type)
            setattr(input_message,attribute, unpacked)
            
        for index in range(input_message.length):
            keyword_arguments = self.new_keyword_arguments_from_message(input_message, index,  specification)
            
            result = method(**keyword_arguments)
            
            self.fill_output_message(output_message, index , result, keyword_arguments, specification)
        
        for type, attribute in self.dtype_to_message_attribute.iteritems():
            array = getattr(output_message, attribute)
            packed = pack_array(array, input_message.length, type)
            setattr(output_message, attribute, packed)
    
    def new_keyword_arguments_from_message(self, input_message, index, specification):
        keyword_arguments = OrderedDictionary()
        for parameter in specification.parameters:
            attribute = self.dtype_to_message_attribute[parameter.datatype]
            argument_value = None
            if parameter.direction == LegacyFunctionSpecification.IN:
                argument_value = getattr(input_message, attribute)[parameter.input_index][index]
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                argument_value = ValueHolder(getattr(input_message, attribute)[parameter.input_index][index])
            if parameter.direction == LegacyFunctionSpecification.OUT:
                argument_value = ValueHolder(None)
        
            keyword_arguments[parameter.name] = argument_value
        return keyword_arguments
        
    def fill_output_message(self, output_message, index, result, keyword_arguments, specification):
        
            
        if not specification.result_type is None:
            attribute = self.dtype_to_message_attribute[specification.result_type]
            getattr(output_message, attribute)[0][index] = result
                
        for parameter in specification.parameters:
            attribute = self.dtype_to_message_attribute[parameter.datatype]
            argument_value = None
            if parameter.direction == LegacyFunctionSpecification.INOUT:
                argument_value = keyword_arguments[parameter.name] 
                getattr(output_message, attribute)[parameter.output_index][index] = argument_value.value
            if parameter.direction == LegacyFunctionSpecification.OUT:
                argument_value = keyword_arguments[parameter.name]
                getattr(output_message, attribute)[parameter.output_index][index] = argument_value.value
    
    def get_dtype_to_count(self, specification):
        dtype_to_count = {}
        
        for parameter in specification.output_parameters:
            count = dtype_to_count.get(parameter.datatype, 0)
            dtype_to_count[parameter.datatype] = count + 1
                
        if not specification.result_type is None:
            count = dtype_to_count.get(specification.result_type, 0)
            dtype_to_count[specification.result_type] = count + 1
        
        return dtype_to_count
        
    @late
    def mapping_from_tag_to_legacy_function(self):
        result = {}
        for x in self.legacy_functions:
            result[x.specification.id] = x
        return result
        
    @late
    def legacy_functions(self):
        attribute_names = dir(self.interface)
        legacy_functions = []
        for x in attribute_names:
            if x.startswith('__'):
                continue
            value = getattr(self.interface, x)
            if isinstance(value, legacy_function):
                legacy_functions.append(value)
        
        legacy_functions.sort(key= lambda x: x.specification.id)
        return legacy_functions
        
