"""
"""
from amuse.support.core import late
from amuse.support.core import print_out
from amuse.support.core import OrderedDictionary

from zlib import crc32

import numpy

class LegacyCall(object):
    """A legacy_call implements the runtime call to the remote process.
    """
    def __init__(self, interface, owner, specification):
        self.interface = interface
        self.owner = owner
        self.specification = specification
    
    def __call__(self, *arguments_list, **keyword_arguments):
        dtype_to_values_and_keyword = {
            'd' : [[],'doubles_in',0],
            'i' : [[],'ints_in',0]
        }
        
        for name, dtype, direction in self.specification.parameters:
            if direction == RemoteFunction.IN or direction == RemoteFunction.INOUT:
                values = dtype_to_values_and_keyword[dtype]
                values[0].append(None)
            
        names_in_argument_list = set([])
        for index, argument in enumerate(arguments_list):
            name, dtype, direction = self.specification.parameters[index]
            names_in_argument_list.add(name)
            values = dtype_to_values_and_keyword[dtype]
            values[0][values[2]] = argument
            values[2] += 1
        
        for parameters in self.specification.dtype_to_input_parameters.values():
            for index, parameter in enumerate(parameters):
                name, dtype, direction = parameter
                if name in keyword_arguments:
                    values = dtype_to_values_and_keyword[dtype]
                    values[0][index] = keyword_arguments[name]
               
        call_keyword_arguments = {}
        for values, keyword, count in dtype_to_values_and_keyword.values():
            call_keyword_arguments[keyword] = values
            
        return self.do_call(self.specification.id, **call_keyword_arguments)
        
    def do_call(self, id, **keyword_arguments):
        self.interface.channel.send_message(id , **keyword_arguments)
        (doubles, ints) = self.interface.channel.recv_message(id)
        
        
        
        number_of_outputs = 0
        for name, dtype, direction in self.specification.parameters:
            if direction == RemoteFunction.OUT or direction == RemoteFunction.INOUT:
                number_of_outputs += 1
       
        
        if number_of_outputs == 0:
            if self.specification.result_type == 'i':
                return ints[0]       
            if self.specification.result_type == 'd':
                return doubles[0] 
        
        if number_of_outputs == 1:
            if len(ints) == 1:
                return ints[0]
            if len(doubles) == 1:
                return doubles[0]
        
        result = OrderedDictionary()
        dtype_to_array = {
            'd' : list(reversed(doubles)),
            'i' : list(reversed(ints))
        }
        for name, dtype, direction in self.specification.parameters:
            if direction == RemoteFunction.OUT or direction == RemoteFunction.INOUT:
                result[name] = dtype_to_array[dtype].pop()
                
        return result
       

class legacy_function(object):
    """The meta information for a function call to a code
    """
    def __init__(self, specification_function):
        self.specification_function = specification_function
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return LegacyCall(instance, owner, self.specification)
        
    def __set__(self, instance, value):
        return
        
    @late
    def specification(self):
        result = self.specification_function()
        if result.name is None:
            result.name = self.specification_function.__name__
        if result.id is None:
            result.id = crc32(result.name)
        return result
        
        
class legacy_global(object):
    """The meta information for globals of the code
    """
    def __init__(self, name , id = None, dtype = 'i'):
        self.name = name
        self.id = id
        self.dtype = dtype
        
        if self.id is None:
            self.id = crc32(self.name)
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return LegacyCall(instance, owner, self.get_specification)()
        
    def __set__(self, instance, value):
        return LegacyCall(instance, None, self.set_specification)(value)
        
        
    
    def to_c_string(self):
        uc = MakeACStringOfALegacyGlobalSpecification()
        uc.specification = self
        return uc.result 
        
    @late
    def set_specification(self):
        result = RemoteFunction()
        result.id = self.id
        result.name = self.name
        result.addParameter('value', dtype=self.dtype, direction=RemoteFunction.IN)
        return result
    @late
    def get_specification(self):
        result = RemoteFunction()
        result.id = self.id
        result.name = self.name
        result.result_type = self.dtype
        return result
        
class RemoteFunction(object):
    IN = object()
    OUT = object()
    INOUT = object()
    
    def __init__(self):
        self.parameters = []
        self.name = None
        self.id = None
        self.result_type = None
        
    def addParameter(self, name, dtype = 'i', direction = IN):
        self.parameters.append((name, dtype, direction))
    
    @property
    def dtype_to_input_parameters(self):
        result = {}
        for name, dtype, direction in self.parameters:
            parameters = result.get(dtype, [])
            parameters.append((name, dtype, direction))
            result[dtype] = parameters
        return result
                
    
    
class MpiChannel(object):
    from mpi4py import MPI
    
    def __init__(self, intercomm):
        self.intercomm = intercomm
        
    def send_message(self, tag, id=0, int_arg1=0, int_arg2=0, doubles_in=[], ints_in=[]):
        header = numpy.array([tag,len(doubles_in), len(ints_in)], dtype='i')
        self.intercomm.Send([header, self.MPI.INT], dest=0, tag=0)
        if doubles_in:
            doubles = numpy.array(doubles_in, dtype='d')
            self.intercomm.Send([doubles, self.MPI.DOUBLE], dest=0, tag=0)
        if ints_in:
            ints = numpy.array(ints_in, dtype='i')
            self.intercomm.Send([ints, self.MPI.INT], dest=0, tag=0)
            
    def recv_message(self, tag):
        header = numpy.empty(3,  dtype='i')
        self.intercomm.Recv([header, self.MPI.INT], source=0, tag=999)
        n_doubles = header[1]
        n_ints = header[2]
        if n_doubles > 0:
            doubles_result = numpy.empty(n_doubles,  dtype='d')
            self.intercomm.Recv([doubles_result, self.MPI.DOUBLE], source=0, tag=999)
        else:
            doubles_result = []
        if n_ints > 0:
            ints_result = numpy.empty(n_ints,  dtype='i')
            self.intercomm.Recv([ints_result, self.MPI.INT], source=0, tag=999)
        else:
            ints_result = []
        if header[0] < 0:
            raise Exception("Not a valid message!")
        return (doubles_result, ints_result)
    


        
                
        
            
