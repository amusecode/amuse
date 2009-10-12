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
import weakref
import atexit

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



def _typecode_to_datatype(typecode):
    if typecode is None:
        return None
    
    mapping = {
        'd':'float64',
        'i':'int32',
        'f':'float32',
        's':'string',
    }
    if typecode in mapping:
        return mapping[typecode]
    
    values = mapping.values()
    if typecode in values:
        return typecode
    
    print typecode
    raise Exception("{0} is not a valid typecode".format(typecode))
    
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
            
            values = dtype_to_values[parameter.datatype]
            values[parameter.input_index] = argument
        
        for index, parameter in enumerate(self.specification.input_parameters):
                if parameter.name in keyword_arguments:
                    values = dtype_to_values[parameter.datatype]
                    values[parameter.input_index] = keyword_arguments[parameter.name]
        
        dtype_to_keyword = {
            'float64' : 'doubles_in',
            'int32'  : 'ints_in',
            'string'  : 'chars_in',
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
            if self.specification.result_type is None:
                return None
                
            if self.specification.result_type == 'int32':
                return ints[0]       
            if self.specification.result_type == 'float64':
                return doubles[0] 
            if self.specification.result_type == 'float32':
                return floats[0] 
        
        if number_of_outputs == 1 \
            and self.specification.result_type is None:
            if len(ints) == 1:
                return ints[0]
            if len(doubles) == 1:
                return doubles[0]
            if len(floats) == 1:
                return floats[0]
        
        result = OrderedDictionary()
        dtype_to_array = {
            'float64' : list(reversed(doubles)),
            'int32' : list(reversed(ints)),
            'float32' : list(reversed(floats))
        }
        for parameter in self.specification.output_parameters:
            result[parameter.name] = dtype_to_array[parameter.datatype].pop()
        
        if not self.specification.result_type is None:
            result["__result"] =  dtype_to_array[self.specification.result_type].pop()
        
        return result
       
    def __str__(self):
        return str(self.specification)
        
        
        
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
        self.datatype = _typecode_to_datatype(dtype)
        
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
        result.addParameter('value', dtype=self.datatype, direction=RemoteFunction.IN)
        return result
        
    @late
    def get_specification(self):
        result = RemoteFunction()
        result.id = self.id
        result.name = self.name
        result.result_type = self.datatype
        return result
     
class Parameter(object):
    
    def __init__(self, name, dtype, direction):
        self.name = name
        self.direction = direction
        self.input_index = -1
        self.output_index = -1
        self.datatype = _typecode_to_datatype(dtype)
        
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
        
        parameters = self.dtype_to_input_parameters.get(parameter.datatype, [])
        parameters.append(parameter)
        parameter.input_index = len(parameters) - 1
        self.dtype_to_input_parameters[parameter.datatype] = parameters
   
    def add_output_parameter(self, parameter):
        self.output_parameters.append(parameter)
        
        parameters = self.dtype_to_output_parameters.get(parameter.datatype, [])
        parameters.append(parameter)
        parameter.output_index = len(parameters) - 1
        self.dtype_to_output_parameters[parameter.datatype] = parameters
   
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
                
    
    def __str__(self):
        typecode_to_name = {'i':'int', 'd':'double', 'f':'float'}
        p = print_out()
        p + 'function: '
        if self.result_type is None:
            p + 'void'
        else:
            p + typecode_to_name[self.result_type]
        p + ' '
        p + self.name
        p + '('
        first = True
        for x in self.input_parameters:
            if first:
                first = False
            else:
                p + ', '
            p + typecode_to_name[x.dtype]
            p + ' '
            p + x.name
        p + ')'
        if self.output_parameters:
            p + '\n'
            p + 'output: '
            first = True
            for x in self.output_parameters:
                if first:
                    first = False
                else:
                    p + ', '
                p + typecode_to_name[x.dtype]
                p + ' '
                p + x.name
            if not self.result_type is None:
                p + ', '
                p + typecode_to_name[self.result_type]
                p + ' '
                p + '__result'
        return p.string
        
    
    def get_result_type(self):
        return self._result_type
    
    def set_result_type(self, value):
        self._result_type = _typecode_to_datatype(value)
        
    result_type = property(get_result_type, set_result_type);
    
class MpiChannel(object):
    
    def __init__(self, name_of_the_worker, number_of_workers, legacy_interface_type):
        self.name_of_the_worker = name_of_the_worker
        self.number_of_workers = number_of_workers
        self.legacy_interface_type = legacy_interface_type
        
    def start(self):
        tried_workers = []
        found = False
        current_type=self.legacy_interface_type
        while not found:
            directory_of_this_module = os.path.dirname(inspect.getfile(current_type))
            full_name_of_the_worker = os.path.join(directory_of_this_module , self.name_of_the_worker)
            found = os.path.exists(full_name_of_the_worker)
            if not found:
            
                tried_workers.append(full_name_of_the_worker)
                current_type = current_type.__bases__[0]
                if current_type is LegacyInterface:
                    raise Exception("The worker application does not exists, it should be at: {0}".format(tried_workers))
            else:
                found = True
        self.intercomm = MPI.COMM_SELF.Spawn(full_name_of_the_worker, None, self.number_of_workers)

    def stop(self):
        self.intercomm.Free()
        
    def send_message(self, tag, id=0, int_arg1=0, int_arg2=0, doubles_in=[], ints_in=[], floats_in=[], chars_in=[], length = 1):
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
        
        if len(chars_in) == 1:
            number_of_characters = len(chars_in[0])+1
        else:
            number_of_characters = 0
        header = numpy.array([tag, length, len(doubles_in), len(ints_in), len(floats_in), number_of_characters], dtype='i')
        
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

        if chars_in:
            as_array_of_bytes = [ord(ch) for ch in chars_in[0]]
            as_array_of_bytes.append(0)
            chars = numpy.array(as_array_of_bytes, dtype=numpy.uint8)
            self.intercomm.Bcast([chars, MPI.CHARACTER], root=MPI.ROOT)
            
    def recv_message(self, tag):
        header = numpy.empty(6,  dtype='i')
        self.intercomm.Recv([header, MPI.INT], source=0, tag=999)
        length = header[1]
        n_doubles = header[2]
        n_ints = header[3]
        n_floats = header[4]
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
        

def stop_interfaces():
    for reference in LegacyInterface.instances:
        x = reference()
        if not x is None:
            x._stop()
        

class LegacyInterface(object):
    instances = []
    atexit.register(stop_interfaces)
    
    def __init__(self, name_of_the_worker = 'muse_worker', number_of_workers = 1):
        self.channel = MpiChannel(name_of_the_worker, number_of_workers, type(self))
        self.channel.start()
        self.instances.append(weakref.ref(self))
        
    def __del__(self):
        self._stop()
    
    def _stop(self):
        if hasattr(self, 'channel'):
            self._stop_worker()
            self.channel.stop()
            del self.channel
        
    @legacy_function
    def _stop_worker():
        function = RemoteFunction()  
        function.id = 0
        return function   



        
                
        
            
