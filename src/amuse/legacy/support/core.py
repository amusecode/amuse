"""
This module implements the code to the define interfaces between python
code and C++ or Fortran legacy codes. It provides the abstract base
class for all legacy codes.
"""

import inspect
import numpy
import weakref
import atexit
import os
import os.path
import cPickle as pickle


import sys
import struct

from zlib import crc32
from mpi4py import MPI
from subprocess import Popen, PIPE

from amuse.support.core import late
from amuse.support.core import print_out
from amuse.support.core import OrderedDictionary
from amuse.legacy.interface import LegacyDocStringProperty


def is_mpd_running():
    """
    Determine if the MPD daemon process is running.
    
    
    Needed for installations of AMUSE in a MPICH2 environment using
    the default MPD daemon. The MPD deamon must be
    running before the first MPI_COMN_SPAWN call is made.
    Returns True for other MPI vendors (OpenMPI)
    
    :returns: Boolean result of check whether MPD daemon is running.
    :rtype: bool
    
    
    >>> is_mpd_running()
    True
    
        
    """
    name_of_the_vendor, version = MPI.get_vendor()
    if name_of_the_vendor == 'MPICH2':
        process = Popen(['mpdtrace'], stdout = PIPE, stderr = PIPE)
        (output_string, error_string) = process.communicate()
        return not (process.returncode == 255)
    else:
        return True


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
    
    __doc__ = LegacyDocStringProperty()
   
    def __init__(self, interface, owner, specification):
        """
        Implementation of the runtime call to the remote process.

        Performs the encoding of python arguments into lists
        of values, sends a message over an MPI channel and
        waits for a result message, decodes this message and
        returns.
        """
        self.interface = interface
        self.owner = owner
        self.specification = specification
    
    def __call__(self, *arguments_list, **keyword_arguments):
        keyword_arguments_for_the_mpi_channel = self.converted_keyword_and_list_arguments( arguments_list, keyword_arguments)
        self.interface.channel.send_message(self.specification.id , **keyword_arguments_for_the_mpi_channel)
        (doubles, ints) = self.interface.channel.recv_message(self.specification.id)
        return self.converted_results(doubles, ints)
        
    """
    Convert results from an MPI message to a return value.
    """
    def converted_results(self, doubles, ints):
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
       
    """
    Convert keyword arguments and list arguments to an MPI message
    """
    def converted_keyword_and_list_arguments(self, arguments_list, keyword_arguments):
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

        return call_keyword_arguments
        
    def __str__(self):
        return str(self.specification)
        
        
class legacy_function(object):
    
    __doc__ = LegacyDocStringProperty()
    
    def __init__(self, specification_function):
        """Decorator for legacy functions.
    
        The decorated function cannot have any arguments. This
        means the decorated function must not have a ``self``
        argument.
        
        The decorated function must return 
        a LegacyFunctionSpecification.
        
            
        >>> class LegacyExample(object):
        ...     @legacy_function
        ...     def evolve():
        ...          specification = LegacyFunctionSpecification()
        ...          return specification
        ...
                    
        :argument specification_function: The function to be decorated
                    
        """
        self.specification_function = specification_function
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return LegacyCall(instance, owner, self.specification)
        
    def __set__(self, instance, value):
        return
        
    @late
    def specification(self):
        """
        The legacy function specification of the call.
        """
        result = self.specification_function()
        if result.name is None:
            result.name = self.specification_function.__name__
        if result.id is None:
            result.id = abs(crc32(result.name))
        if result.description is None:
            import pydoc
            result.description = pydoc.getdoc(self.specification_function)
        return result
    
        
class legacy_global(object):
  
    def __init__(self, name , id = None, dtype = 'i'):
        """
        Decorator for legacy globals.
        
        *to be removed*
        """
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
        result = LegacyFunctionSpecification()
        result.id = self.id
        result.name = self.name
        result.addParameter('value', dtype=self.datatype, direction=LegacyFunctionSpecification.IN)
        return result
        
    @late
    def get_specification(self):
        result = LegacyFunctionSpecification()
        result.id = self.id
        result.name = self.name
        result.result_type = self.datatype
        return result
     
class ParameterSpecification(object):
    def __init__(self, name, dtype, direction, description):
        """Specification of a parameter of a legacy function 
        """
        self.name = name
        self.direction = direction
        self.input_index = -1
        self.output_index = -1
        self.description = description
        self.datatype = _typecode_to_datatype(dtype)
        
    def is_input(self):
        return ( self.direction == LegacyFunctionSpecification.IN 
            or self.direction == LegacyFunctionSpecification.INOUT)
            
    
    def is_output(self):
        return ( self.direction == LegacyFunctionSpecification.OUT 
            or self.direction == LegacyFunctionSpecification.INOUT)
                    
class LegacyFunctionSpecification(object):
    """
    Specification of a legacy function.
    Describes the name, result type and parameters of a
    legacy function. 
    
    The legacy functions are implemented by legacy codes. 
    The implementation of legacy functions is in C/C++ or Fortran.
    To interact with these functions a specification of the
    legacy function is needed.
    This specification is used to determine how to encode
    and decode the parameters and results of the function.
    Objects of this class describe the specification of one
    function.
    
    >>> specification = LegacyFunctionSpecification()
    >>> specification.name = "test"
    >>> specification.addParameter("one", dtype="int32", direction = specification.IN)
    >>> specification.result_type = "int32"
    >>> print specification
    function: int test(int one)
    
    """
    
    IN = object()
    """Used to specify that a parameter is used as an input parameter, passed by value"""
    
    OUT = object()
    """Used to specify that a parameter is used as an output parameter, passed by reference"""
    
    INOUT = object()
    """Used to specify that a parameter is used as an input and an outpur parameter, passed by reference"""
    
    def __init__(self):
        self.parameters = []
        self.name = None
        self.id = None
        self.result_type = None
        self.description = None
        self.input_parameters = []
        self.output_parameters = []
        self.dtype_to_input_parameters = {}
        self.dtype_to_output_parameters = {}
        self.can_handle_array = False
        self.result_doc = ''
        
    def addParameter(self, name, dtype = 'i', direction = IN, description = ""):
        """
        Extend the specification with a new parameter.
        
        The sequence of calls to addParameter is important. The first
        call will be interpreted as the first argument, the second
        call as the second argument etc.
        
        :argument name: Name of the parameter, used in documentation and function generation
        :argument dtype: Datatype specification string
        :argument direction: Direction of the argument, can be IN, OUT or INOUT
        :argument description: Description of the argument, for documenting purposes
        """
        parameter = ParameterSpecification(name, dtype, direction, description)
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
        typecode_to_name = {'int32':'int', 'float64':'double', 'float32':'float'}
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
            p + typecode_to_name[x.datatype]
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
        
    
    def _get_result_type(self):
        return self._result_type
    
    def _set_result_type(self, value):
        self._result_type = _typecode_to_datatype(value)
        
    result_type = property(_get_result_type, _set_result_type);
    
class MessageChannel(object):
    
    def get_full_name_of_the_worker(self, type):
        tried_workers = []
        found = False
        current_type=type
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
        return full_name_of_the_worker


import time
import ctypes

clib_library = ctypes.CDLL("libc.so.6")

memcpy = clib_library.memcpy
memcpy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t]

class MpiChannel(MessageChannel):
    
    def __init__(self, name_of_the_worker, number_of_workers, legacy_interface_type = None, debug_with_gdb = False):
        self.name_of_the_worker = name_of_the_worker
        self.number_of_workers = number_of_workers
        if not legacy_interface_type is None:
            self.full_name_of_the_worker = self.get_full_name_of_the_worker( legacy_interface_type)
        else:
            self.full_name_of_the_worker = self.name_of_the_worker
        self.debug_with_gdb = debug_with_gdb
        
        self.cached = None
        
    def start(self):
        if self.debug_with_gdb:
            if not 'DISPLAY' in os.environ:
                arguments = None
                command = self.full_name_of_the_worker
            else:
                arguments = ['-hold', '-display', os.environ['DISPLAY'], '-e', 'gdb',  self.full_name_of_the_worker]
                command = 'xterm'
        else:
            arguments = None
            command = self.full_name_of_the_worker
        self.intercomm = MPI.COMM_SELF.Spawn(command, arguments, self.number_of_workers)

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
        if chars_in:
            try:
                if not isinstance(chars_in[0], str):
                    length = len(chars_in[0])
                    print "length:", length
            except:
                pass
                    
        
        #if len(chars_in) == 1:
        #    number_of_characters = len(chars_in[0])+1
        #else:
        #    number_of_characters = 0
       
        header = numpy.array([tag, length, len(doubles_in), len(ints_in), len(floats_in), len(chars_in)], dtype='i')
        self.intercomm.Bcast([header, MPI.INT], root=MPI.ROOT)
        
        if doubles_in:
            #t0 = time.time()
            if self.cached  is None:
                doubles = numpy.empty(length * len(doubles_in), dtype='d')
                self.cached = doubles
            else:
                if self.cached.size == length * len(doubles_in):
                    doubles = self.cached
                else:
                    doubles = numpy.empty(length * len(doubles_in), dtype='d')
                    self.cached = doubles
            
            dst_pointer = doubles.__array_interface__['data'][0]
            for i in range(len(doubles_in)):
                offset = i * length
                doubles[offset:offset+length] = doubles_in[i]
                #print  (8.0 * length) / 1024.0 / 1024.0
                #memcpy(ctypes.c_void_p(dst_pointer + (offset * 8)), ctypes.c_void_p( doubles_in[i].__array_interface__['data'][0]), 8 * length)
            
            self.intercomm.Bcast([doubles, MPI.DOUBLE], root=MPI.ROOT)
            #t1 = time.time()
            #print "D", t1 - t0
        if ints_in:
            ints = numpy.empty(length * len(ints_in), dtype='i')
            for i in range(len(ints_in)):
                offset = i * length
                ints[offset:offset+length] = ints_in[i]
            
            self.intercomm.Bcast([ints, MPI.INT], root=MPI.ROOT)
            
        if floats_in:
            floats = numpy.array(floats_in, dtype='f')
            self.intercomm.Bcast([floats, MPI.FLOAT], root=MPI.ROOT)

        if chars_in:
            print "C"
            offsets = numpy.zeros(length * len(chars_in), dtype='i')
            offset = 0
            index = 0
            for strings in chars_in:
                if length == 1:
                    offset += len(strings)
                    offsets[index] = offset
                    offset += 1
                    index += 1
                else:
                    for string in strings:
                        offset += len(string)
                        offsets[index] = offset
                        offset += 1
                        index += 1
                
            print offsets
            self.intercomm.Bcast([offsets, MPI.INT], root=MPI.ROOT)    
            bytes = []
            for strings in chars_in:
                if length == 1:
                    bytes.extend([ord(ch) for ch in strings])
                    bytes.append(0)
                    print "str:", strings
                else:
                    for string in strings:
                        bytes.extend([ord(ch) for ch in string])
                        bytes.append(0)
              
            print bytes
            chars = numpy.array(bytes, dtype=numpy.uint8)
            self.intercomm.Bcast([chars, MPI.CHARACTER], root=MPI.ROOT)
       
        
        
    def recv_message(self, tag):
        header = numpy.zeros(6,  dtype='i')
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
        


class MultiprocessingMPIChannel(MessageChannel):
    
    def __init__(self, name_of_the_worker, number_of_workers, legacy_interface_type, debug_with_gdb = False):
        self.name_of_the_worker = name_of_the_worker
        self.number_of_workers = number_of_workers
        if not legacy_interface_type is None:
            self.full_name_of_the_worker = self.get_full_name_of_the_worker( legacy_interface_type)
        else:
            self.full_name_of_the_worker = self.name_of_the_worker
        self.debug_with_gdb = debug_with_gdb
    
    def start(self):
        name_of_dir = "/tmp/amuse_"+os.getenv('USER')
        self.name_of_the_socket, self.server_socket = self._createAServerUNIXSocket(name_of_dir)
        environment = os.environ.copy()
        if 'PYTHONPATH' in environment:
            environment['PYTHONPATH'] = environment['PYTHONPATH'] + ':' +  self._extra_path_item(__file__)
        else:
            environment['PYTHONPATH'] =  self._extra_path_item(__file__)
        template = "from amuse.legacy.support import core\nm = core.MultiprocessingMPIChannel('{0}',{1},None,{2})\nm.run_mpi_channel('{3}')"
        code_string = template.format(self.full_name_of_the_worker, self.number_of_workers, self.debug_with_gdb, self.name_of_the_socket)
        self.process =  Popen([sys.executable, "-c", code_string], env = environment)
        self.client_socket, undef = self.server_socket.accept()
        
    def stop(self):
        self._send(self.client_socket,('stop',(),))
        result = self._recv(self.client_socket)    
        self.process.wait()
        self.client_socket.close()
        self.server_socket.close()
        self._remove_socket(self.name_of_the_socket)
        
    def run_mpi_channel(self, name_of_the_socket):
        channel = MpiChannel(self.full_name_of_the_worker, self.number_of_workers, None, self.debug_with_gdb)
        channel.start()
        socket = self._createAClientUNIXSocket(name_of_the_socket)
        try:
            is_running = True
            while is_running:
                message, args = self._recv(socket)
                result = None
                if message == 'stop':
                    channel.stop()
                    is_running = False
                if message == 'send_message':
                    result = channel.send_message(*args)
                if message == 'recv_message':
                    result = channel.recv_message(*args)
                self._send(socket, result)
        finally:
            socket.close()
         
    def send_message(self, tag, id=0, int_arg1=0, int_arg2=0, doubles_in=[], ints_in=[], floats_in=[], chars_in=[], length = 1):
        self._send(self.client_socket, ('send_message',(tag, id, int_arg1, int_arg2, doubles_in, ints_in, floats_in, chars_in, length),))
        result = self._recv(self.client_socket)    
        return result
            
    def recv_message(self, tag):
        self._send(self.client_socket, ('recv_message',(tag,),))
        result = self._recv(self.client_socket)    
        return result
    
    def _send(self, client_socket, message):
        message_string = pickle.dumps(message)
        header = struct.pack("i", len(message_string))
        client_socket.sendall(header)
        client_socket.sendall(message_string)
        
    def _recv(self, client_socket):
        header = self._recvall(client_socket, 4)
        length = struct.unpack("i", header)
        message_string = self._recvall(client_socket, length[0])
        return pickle.loads(message_string)
        
    def _recvall(self, client_socket, number_of_bytes):
        block_size = 4096
        bytes_left = number_of_bytes
        blocks = []
        while bytes_left > 0:
            if bytes_left < block_size:
                block_size = bytes_left
            block = client_socket.recv(block_size)
            blocks.append(block)
            bytes_left -= len(block)
        return ''.join(blocks)
        
            
    def _createAServerUNIXSocket(self, name_of_the_directory, name_of_the_socket = None):
        import uuid
        import socket
        
        if name_of_the_socket == None:
            name_of_the_socket = os.path.join(name_of_the_directory,str(uuid.uuid1()))
            
        if not os.path.exists(name_of_the_directory):
            os.makedirs(name_of_the_directory)
            
        server_socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        self._remove_socket(name_of_the_socket)
        server_socket.bind(name_of_the_socket)
        server_socket.listen(5)
        return (name_of_the_socket, server_socket)

    def _createAClientUNIXSocket(self, name_of_the_socket):
        import socket
        client_socket = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        #client_socket.settimeout(0)
        client_socket.connect(name_of_the_socket)
        return client_socket
                
    def _remove_socket(self, name_of_the_socket):
        try:
            os.remove(name_of_the_socket)
        except OSError:
            pass
            
    def _extra_path_item(self, path_of_the_module):
        result = ''
        for x in sys.path:
            if path_of_the_module.startswith(x):
                if len(x) > len(result):
                    result = x
        return result
        
            

def stop_interfaces():
    for reference in LegacyInterface.instances:
        x = reference()
        if not x is None:
            x._stop()
        

class LegacyInterface(object):
    instances = []
    atexit.register(stop_interfaces)
    channel_factory = MpiChannel
    
    def __init__(self, name_of_the_worker = 'muse_worker', number_of_workers = 1, debug_with_gdb = False):
        self.channel = self.channel_factory(name_of_the_worker, number_of_workers, type(self), debug_with_gdb)
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
        function = LegacyFunctionSpecification()  
        function.id = 0
        return function   

