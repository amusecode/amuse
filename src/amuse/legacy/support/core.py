"""
"""
from amuse.support.core import late
from amuse.support.core import print_out
from amuse.support.core import OrderedDictionary

from zlib import crc32
from mpi4py import MPI
import os.path

import inspect
import numpy

class DataType(object):
    pass
    
class IntDataType(DataType):
    short_key = 'i'
    mpi_type = MPI.INT
    pass
    
class DoubleDataType(DataType):
    short_key = 'd'
    mpi_type = MPI.DOUBLE
    pass
    
class FloatDataType(DataType):
    short_key = 'f'
    mpi_type = MPI.FLOAT
    pass

    
class LegacyCall(object):
    """A legacy_call implements the runtime call to the remote process.
    """
    def __init__(self, interface, owner, specification):
        self.interface = interface
        self.owner = owner
        self.specification = specification
    
    def __call__(self, *arguments_list, **keyword_arguments):
        
        
        dtype_to_values = self.specification.new_dtype_to_values()
        
        names_in_argument_list = set([])
        for index, argument in enumerate(arguments_list):
            parameter = self.specification.input_parameters[index]
            names_in_argument_list.add(parameter.name)
            
            values = dtype_to_values[parameter.dtype]
            values[parameter.input_index] = argument
        
        for index, parameter in enumerate(self.specification.input_parameters):
                if parameter.name in keyword_arguments:
                    values = dtype_to_values[parameter.dtype]
                    values[parameter.input_index] = keyword_arguments[parameter.name]
        
        dtype_to_keyword = {
            'd' : 'doubles_in',
            'i' : 'ints_in'
        }       
        call_keyword_arguments = {}
        for dtype, values in dtype_to_values.iteritems():
            keyword = dtype_to_keyword[dtype]
            call_keyword_arguments[keyword] = values
            
        return self.do_call(self.specification.id, **call_keyword_arguments)
        
    def do_call(self, id, **keyword_arguments):
        self.interface.channel.send_message(id , **keyword_arguments)
        (doubles, ints) = self.interface.channel.recv_message(id)
        floats = []
        
        
        number_of_outputs = len(self.specification.output_parameters)
        
        if number_of_outputs == 0:
            if self.specification.result_type == 'i':
                return ints[0]       
            if self.specification.result_type == 'd':
                return doubles[0] 
            if self.specification.result_type == 'f':
                return floats[0] 
        
        if number_of_outputs == 1:
            if len(ints) == 1:
                return ints[0]
            if len(doubles) == 1:
                return doubles[0]
            if len(floats) == 1:
                return floats[0]
        
        result = OrderedDictionary()
        dtype_to_array = {
            'd' : list(reversed(doubles)),
            'i' : list(reversed(ints))
        }
        for parameter in self.specification.output_parameters:
            result[parameter.name] = dtype_to_array[parameter.dtype].pop()
                
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
            result.id = abs(crc32(result.name))
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
     
class Parameter(object):
    def __init__(self, name, dtype, direction):
        self.name = name
        self.dtype = dtype
        self.direction = direction
        self.input_index = -1
        self.output_index = -1
        
    def is_input(self):
        return ( self.direction == RemoteFunction.IN 
            or self.direction == RemoteFunction.INOUT)
            
    
    def is_output(self):
        return ( self.direction == RemoteFunction.OUT 
            or self.direction == RemoteFunction.INOUT)
                    
class RemoteFunction(object):
    IN = object()
    OUT = object()
    INOUT = object()
    
    def __init__(self):
        self.parameters = []
        self.name = None
        self.id = None
        self.result_type = None
        self.input_parameters = []
        self.output_parameters = []
        self.dtype_to_input_parameters = {}
        self.dtype_to_output_parameters = {}
        self.can_handle_array = False
        
    def addParameter(self, name, dtype = 'i', direction = IN):
        parameter = Parameter(name, dtype, direction)
        self.parameters.append(parameter)
        
        if parameter.is_input():
            self.add_input_parameter(parameter)
        if parameter.is_output():
            self.add_output_parameter(parameter)
            
    def add_input_parameter(self, parameter):
        self.input_parameters.append(parameter)
        
        parameters = self.dtype_to_input_parameters.get(parameter.dtype, [])
        parameters.append(parameter)
        parameter.input_index = len(parameters) - 1
        self.dtype_to_input_parameters[parameter.dtype] = parameters
   
    def add_output_parameter(self, parameter):
        self.output_parameters.append(parameter)
        
        parameters = self.dtype_to_output_parameters.get(parameter.dtype, [])
        parameters.append(parameter)
        parameter.output_index = len(parameters) - 1
        self.dtype_to_output_parameters[parameter.dtype] = parameters
   
    def new_dtype_to_values(self):
        result = {}
        for dtype, parameters in self.dtype_to_input_parameters.iteritems():
            result[dtype] =  [None] * len(parameters)   
        return result
    
    def prepare_output_parameters(self):
        for dtype, parameters in self.dtype_to_output_parameters.iteritems():
            if dtype == self.result_type:
                offset = 1
            else:
                offset = 0
            for index, parameter in enumerate(parameters):
                parameter.output_index = offset + index
    
class MpiChannel(object):
    
    def __init__(self, intercomm):
        self.intercomm = intercomm
        
    def send_message(self, tag, id=0, int_arg1=0, int_arg2=0, doubles_in=[], ints_in=[], floats_in=[], length = 1):
        if doubles_in:
            try:
                length = len(doubles_in[0])
            except:
                pass
        if ints_in:
            try:
                length = len(ints_in[0])
            except:
                pass
        header = numpy.array([tag, length, len(doubles_in), len(ints_in), len(floats_in)], dtype='i')
        
        self.intercomm.Bcast([header, MPI.INT], root=MPI.ROOT)
        if doubles_in:
            doubles = numpy.zeros(length * len(doubles_in), dtype='d')
            for i in range(len(doubles_in)):
                offset = i * length
                doubles[offset:offset+length] = doubles_in[i]
            
            self.intercomm.Bcast([doubles, MPI.DOUBLE], root=MPI.ROOT)
        if ints_in:
            ints = numpy.zeros(length * len(ints_in), dtype='i')
            for i in range(len(ints_in)):
                offset = i * length
                ints[offset:offset+length] = ints_in[i]
            self.intercomm.Bcast([ints, MPI.INT], root=MPI.ROOT)
        if floats_in:
            floats = numpy.array(floats_in, dtype='f')
            self.intercomm.Bcast([floats, MPI.FLOAT], root=MPI.ROOT)
            
    def recv_message(self, tag):
        header = numpy.empty(5,  dtype='i')
        self.intercomm.Recv([header, MPI.INT], source=0, tag=999)
        length = header[1]
        n_doubles = header[2]
        n_ints = header[3]
        n_floats = header[4]
        if length > 1:
            print "l:", length
        if n_doubles > 0:
            doubles_mpi = numpy.empty(n_doubles * length,  dtype='d')
            self.intercomm.Recv([doubles_mpi, MPI.DOUBLE], source=0, tag=999)
            doubles_result = []
            if length > 1:
                for i in range(n_doubles):
                    offset = i * length
                    doubles_result.append(doubles_mpi[offset:offset+length])
            else:
                doubles_result = doubles_mpi
        else:
            doubles_result = []
        if n_ints > 0:
            ints_mpi = numpy.empty(n_ints * length,  dtype='i')
            self.intercomm.Recv([ints_mpi, MPI.INT], source=0, tag=999)
            ints_result = []
            if length > 1:
                for i in range(n_ints):
                    offset = i * length
                    ints_result.append(ints_mpi[offset:offset+length])
            else:
                ints_result = ints_mpi
                
        else:
            ints_result = []
        if n_floats > 0:
            floats_result = numpy.empty(n_floats * length,  dtype='f')
            self.intercomm.Recv([floats_result, MPI.FLOAT], source=0, tag=999)
        else:
            floats_result = []
        if header[0] < 0:
            raise Exception("Not a valid message!")
        return (doubles_result, ints_result)
        

class LegacyInterface(object):
    
    def __init__(self, name_of_the_worker = 'muse_worker', number_of_workers = 1):
        self._start_worker(name_of_the_worker, number_of_workers)
        
    def _start_worker(self,  name_of_the_worker, number_of_workers):
        directory_of_this_module = os.path.dirname(inspect.getfile(type(self)))
        full_name_of_the_worker = os.path.join(directory_of_this_module , name_of_the_worker)
        if not os.path.exists(full_name_of_the_worker):
            raise Exception("The worker application does not exists, it should be at: %s".format(full_name_of_the_worker))
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None, number_of_workers)
        self.channel = MpiChannel(self.intercomm)

    def __del__(self):
        self._stop_worker()
        
    @legacy_function
    def _stop_worker():
        function = RemoteFunction()  
        function.id = 0
        return function   



        
                
        
            
