import weakref
import atexit
import os.path
import sys
import logging
import pydoc
import traceback
import random
import sys


from subprocess import Popen, PIPE

from amuse.support import exceptions
from amuse.support.core import late
from amuse.support.core import print_out
from amuse.support.core import OrderedDictionary
from amuse.support.options import OptionalAttributes, option
from amuse.rfi.tools.create_definition import CodeDocStringProperty
from amuse.rfi.channel import MpiChannel
from amuse.rfi.channel import MultiprocessingMPIChannel
from amuse.rfi.channel import DistributedChannel
from amuse.rfi.channel import SocketChannel
from amuse.rfi.channel import is_mpd_running

CODE_LOG = logging.getLogger("code")
if CODE_LOG.level == logging.NOTSET:
    CODE_LOG.setLevel(logging.WARN)

"""
This module implements the code to the define interfaces between python
code and C++ or Fortran codes. It provides the abstract base
class for all community codes.
"""

def ensure_mpd_is_running():
    from mpi4py import MPI
    if not is_mpd_running():
        name_of_the_vendor, version = MPI.get_vendor()
        if name_of_the_vendor == 'MPICH2':
            process = Popen(['nohup','mpd'])

def _typecode_to_datatype(typecode):
    if typecode is None:
        return None
    
    mapping = {
        'd':'float64',
        'i':'int32',
        'f':'float32',
        's':'string',
        'b':'bool',
        'l':'int64',
    }
    if typecode in mapping:
        return mapping[typecode]
    
    values = mapping.values()
    if typecode in values:
        return typecode
    
    raise exceptions.AmuseException("{0} is not a valid typecode".format(typecode))
    


    
class CodeFunction(object):
    
    __doc__ = CodeDocStringProperty()
   
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
        dtype_to_values = self.converted_keyword_and_list_arguments( arguments_list, keyword_arguments)
        
        handle_as_array = self.must_handle_as_array(dtype_to_values)
        
        if not self.owner is None:
            logging.getLogger("code").info("start call '%s.%s'",self.owner.__name__, self.specification.name)
        
        call_id = random.randint(0, 1000)
        
        try:
            self.interface.channel.send_message(call_id, self.specification.id, dtype_to_arguments = dtype_to_values)
            
            dtype_to_result = self.interface.channel.recv_message(call_id, self.specification.id, handle_as_array)
        except Exception, ex:
            CODE_LOG.info("Exception when calling function '{0}', of code '{1}', exception was '{2}'".format(self.specification.name, type(self.interface).__name__, ex))
            raise exceptions.CodeException("Exception when calling function '{0}', of code '{1}', exception was '{2}'".format(self.specification.name, type(self.interface).__name__, ex))
        
        result = self.converted_results(dtype_to_result, handle_as_array)
        
        if not self.owner is None:
            logging.getLogger("code").info("end call '%s.%s'",self.owner.__name__, self.specification.name)
        
        return result
    
    def async(self, *arguments_list, **keyword_arguments):
        dtype_to_values = self.converted_keyword_and_list_arguments( arguments_list, keyword_arguments)
        
        handle_as_array = self.must_handle_as_array(dtype_to_values)
        
        call_id = random.randint(0, 1000)
              
        self.interface.channel.send_message(call_id, self.specification.id, dtype_to_arguments = dtype_to_values)
        
        request = self.interface.channel.nonblocking_recv_message(call_id, self.specification.id, handle_as_array)
        
        def handle_result(function):
            try:
                dtype_to_result = function()
            except Exception, ex:
                raise exceptions.CodeException("Exception when calling legacy code '{0}', exception was '{1}'".format(self.specification.name, ex))
            return self.converted_results(dtype_to_result, handle_as_array)
            
        request.add_result_handler(handle_result)
        return request
        
        
    
    def must_handle_as_array(self, keyword_arguments):
        result = False
        for x in keyword_arguments.values():
            if x and hasattr(x[0],"__len__"):
                result = len(x[0]) > 0
                break
        return result
        
    """
    Convert results from an MPI message to a return value.
    """
    def converted_results(self, dtype_to_result, must_handle_as_array):
        
        number_of_outputs = len(self.specification.output_parameters)
        
        if number_of_outputs == 0:
            if self.specification.result_type is None:
                return None
            return dtype_to_result[self.specification.result_type][0]
            
        if number_of_outputs == 1 \
            and self.specification.result_type is None:
            
            for value in dtype_to_result.values():
                if len(value) == 1:
                    if must_handle_as_array:
                        return value
                    else:
                        return value[0]
            
        result = OrderedDictionary()
        dtype_to_array = {}
        
        for key, value in dtype_to_result.iteritems():
            dtype_to_array[key] = list(reversed(value))
        
        if not self.specification.result_type is None:
            return_value =  dtype_to_array[self.specification.result_type].pop()
        
        for parameter in self.specification.output_parameters:
            result[parameter.name] = dtype_to_array[parameter.datatype].pop()
        
        if not self.specification.result_type is None:
            result["__result"] =  return_value
        
        return result
       
    """
    Convert keyword arguments and list arguments to an MPI message
    """
    def converted_keyword_and_list_arguments(self, arguments_list, keyword_arguments):
        dtype_to_values = self.specification.new_dtype_to_values()
        
        input_parameters_seen = set(map(lambda x : x.name, self.specification.input_parameters))
        names_in_argument_list = set([])
        for index, argument in enumerate(arguments_list):
            parameter = self.specification.input_parameters[index]
            names_in_argument_list.add(parameter.name)
            
            values = dtype_to_values[parameter.datatype]
            values[parameter.input_index] = argument
            input_parameters_seen.remove(parameter.name)
        
        for index, parameter in enumerate(self.specification.input_parameters):
            if parameter.name in keyword_arguments:
                values = dtype_to_values[parameter.datatype]
                values[parameter.input_index] = keyword_arguments[parameter.name]
                input_parameters_seen.remove(parameter.name)
        
        for parameter in self.specification.input_parameters:
            if (parameter.name in input_parameters_seen) and parameter.has_default_value():
                values = dtype_to_values[parameter.datatype]
                values[parameter.input_index] = parameter.default
                input_parameters_seen.remove(parameter.name)
                
        if input_parameters_seen:
            raise exceptions.CodeException("Not enough parameters in call, missing " + str(sorted(input_parameters_seen)))
         
        return dtype_to_values
        
    def __str__(self):
        return str(self.specification)
        
        
class legacy_function(object):
    
    __doc__ = CodeDocStringProperty()
    
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
        >>> x = LegacyExample()
        >>> x.evolve.specification #doctest: +ELLIPSIS
        <amuse.rfi.core.LegacyFunctionSpecification object at 0x...>
        >>> LegacyExample.evolve #doctest: +ELLIPSIS
        <amuse.rfi.core.legacy_function object at 0x...>
        >>> x.evolve #doctest: +ELLIPSIS
        <amuse.rfi.core.CodeFunction object at 0x...>
        
                    
        :argument specification_function: The function to be decorated
                    
        """
        self.specification_function = specification_function
        
    def __get__(self, instance, owner):
        if instance is None:
            return self
        return CodeFunction(instance, owner, self.specification)
        
    def __set__(self, instance, value):
        return
        
    @late
    def specification(self):
        """
        Returns the specification for the call.
        """
        result = self.specification_function()
        if result.name is None:
            result.name = self.specification_function.__name__
        if result.id is None:
            result.id = abs(self.crc32(result.name))
        if result.description is None:
            result.description = pydoc.getdoc(self.specification_function)
        return result
    
    def is_compiled_file_up_to_date(self, time_of_the_compiled_file):
        name_of_defining_file = self.specification_function.func_code.co_filename
        if os.path.exists(name_of_defining_file):
            time_of_defining_file = os.stat(name_of_defining_file).st_mtime
            return time_of_defining_file <= time_of_the_compiled_file
        return True
        
    @late
    def crc32(self):
        try:
        
            from zlib import crc32
            try:
                if crc32('amuse')&0xffffffff == 0xc0cc9367:
                    return crc32
            except Exception:
                #python 3, crc32 needs bytes...
                def python3_crc32(x):
                    x = crc32(bytes(x, 'ascii'))
                    return x - ((x & 0x80000000) <<1)
                if python3_crc32('amuse')&0xffffffff == 0xc0cc9367:
                    return python3_crc32
        except Exception:
            pass
        try:
            from binascii import crc32
            try:
                if crc32('amuse')&0xffffffff == 0xc0cc9367:
                    return crc32
            except Exception:
                #python 3, crc32 needs bytes...
                def python3_crc32(x):
                    x = crc32(bytes(x, 'ascii'))
                    return x - ((x & 0x80000000) <<1)
                if python3_crc32('amuse')&0xffffffff == 0xc0cc9367:
                    return python3_crc32
        except Exception:
            pass
        
        raise Exception("No working crc32 implementation found!")
     
class ParameterSpecification(object):
    def __init__(self, name, dtype, direction, description, default = None, unit = None):
        """Specification of a parameter of a legacy function 
        """
        self.name = name
        self.direction = direction
        self.input_index = -1
        self.output_index = -1
        self.description = description
        self.datatype = _typecode_to_datatype(dtype)
        self.default = default
        self.unit = unit
        
    def is_input(self):
        return ( self.direction == LegacyFunctionSpecification.IN 
            or self.direction == LegacyFunctionSpecification.INOUT)
            
    
    def is_output(self):
        return ( self.direction == LegacyFunctionSpecification.OUT 
            or self.direction == LegacyFunctionSpecification.INOUT)
                    

    def has_default_value(self):
        return not self.default is None
    
    
    
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
    >>> specification.addParameter("two", dtype="float64", direction = specification.OUT)
    >>> specification.result_type = "int32"
    >>> print specification
    function: int test(int one)
    output: double two, int __result
    
    """
    
    IN = object()
    """Used to specify that a parameter is used as an input parameter, passed by value"""
    
    OUT = object()
    """Used to specify that a parameter is used as an output parameter, passed by reference"""
    
    INOUT = object()
    """Used to specify that a parameter is used as an input and an outpur parameter, passed by reference"""
    
    LENGTH = object()
    """Used to specify that a parameter is used as the length parameter for the other parameters"""
    
    def __init__(self):
        self.parameters = []
        self.name = None
        self.id = None
        self.result_type = None
        self.result_unit = None
        self.description = None
        self.input_parameters = []
        self.output_parameters = []
        self.dtype_to_input_parameters = {}
        self.dtype_to_output_parameters = {}
        self.can_handle_array = False
        self.must_handle_array = False
        self.result_doc = ''
        
    def addParameter(self, name, dtype = 'i', direction = IN, description = "", default = None, unit = None):
        """
        Extend the specification with a new parameter.
        
        The sequence of calls to addParameter is important. The first
        call will be interpreted as the first argument, the second
        call as the second argument etc.
        
        :argument name: Name of the parameter, used in documentation and function generation
        :argument dtype: Datatype specification string
        :argument direction: Direction of the argument, can be IN, OUT or INOUT
        :argument description: Description of the argument, for documenting purposes
        :argument default: An optional default value for the parameter
        """
        parameter = ParameterSpecification(name, dtype, direction, description, default, unit)
        self.parameters.append(parameter)
        
        if parameter.is_input():
            self.add_input_parameter(parameter)
        if parameter.is_output():
            self.add_output_parameter(parameter)
            
    def add_input_parameter(self, parameter):
        has_default_parameters = any(map(lambda x : x.has_default_value(), self.input_parameters))
        if has_default_parameters and not parameter.has_default_value():
            raise exceptions.AmuseException("non default argument '{0}' follows default argument".format(parameter.name))
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
        typecode_to_name = {'int32':'int', 'float64':'double', 'float32':'float', 'string':'string', 'int64':'long', 'bool':'bool' }
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
                p + typecode_to_name[x.datatype]
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
    
    def iter_optional_input_parameters(self):
        for x in self.input_parameters:
            if x.has_default_value():
                yield x
                
    result_type = property(_get_result_type, _set_result_type);



def stop_interfaces(exceptions = []):
    """
    Stop the workers of all instantiated interfaces.
    
    All instantiated interfaces will become unstable
    after this call!
    """
    for reference in CodeInterface.instances:
        x = reference()
        if not x is None and x.__class__.__name__ not in exceptions:
            try:
                x._stop()
            except:
                pass

class CodeInterface(OptionalAttributes):
    """
    Abstract base class for all interfaces to legacy codes.
    
    When a subclass is instantiated, a number of subprocesses
    will be started. These subprocesses are called workers
    as they implement the interface and do the actual work
    of the instantiated object.
    """
    instances = []
    is_stop_interfaces_registered = False
    
    def __init__(self, name_of_the_worker = 'worker_code', **options):
        """
        Instantiates an object, starting the worker.
        
        Deleting the instance, with ``del``, will stop
        the worker.
        
        The worker can be started with a gdb session. The code
        will start gdb in a new xterm window. To enable this
        debugging support, the ``DISPLAY`` environment variable must be
        set and the X display must be accessable, ``xhost +``.
        
        :argument name_of_the_worker: The filename of the application to start
        :argument number_of_workers: Number of applications to start. The application must have parallel MPI support if this is more than 1.
        :argument debug_with_gdb: Start the worker(s) in a gdb session in a separate xterm
        :argument hostname: Start the worker on the node with this name
        """
        OptionalAttributes.__init__(self, **options)
           
       
        self.instances.append(weakref.ref(self))
        #
        #ave: no more redirection in the code
        #1) does not seem to work in fortran correctly
        #2) seems to break on hydra
        #
        #if not 'debugger' in options:
        #    self._redirect_outputs(*self.redirection_filenames)
        #
        
        if self.must_start_worker:
            self._start(name_of_the_worker = name_of_the_worker, **options)
        
    def __del__(self):
        self._stop()
    
    def _check_if_worker_is_up_to_date(self):
        self.channel.check_if_worker_is_up_to_date(self)
        
    
    def _start(self, name_of_the_worker = 'worker_code', interpreter_executable = None, **options):
        self.channel = self.channel_factory(name_of_the_worker, type(self), interpreter_executable = interpreter_executable, **options)
        
        self._check_if_worker_is_up_to_date()
        
        self.channel.redirect_stdout_file = self.redirection_filenames[0]
        self.channel.redirect_stderr_file = self.redirection_filenames[1]
        self.channel.polling_interval_in_milliseconds = self.polling_interval_in_milliseconds
        
        self.channel.start()
        
        # must register stop interfaces after channel start
        # if done before, the mpi atexit will be registered 
        # incorrectly
        self.ensure_stop_interface_at_exit()
        
        if self.channel.is_polling_supported():
            if self.polling_interval_in_milliseconds > 0:
                self.internal__set_message_polling_interval(int(self.polling_interval_in_milliseconds * 1000))
        
    @option(type="int", sections=("channel",))
    def polling_interval_in_milliseconds(self):
        return 0
        
    @classmethod
    def ensure_stop_interface_at_exit(cls):
        if not cls.is_stop_interfaces_registered:
            atexit.register(stop_interfaces)
            cls.is_stop_interfaces_registered = True
        
    def _stop(self):
        if hasattr(self, 'channel'):
            if not self.channel is None and self.channel.is_active():
                self._stop_worker()
                self.channel.stop()
                self.channel = None
            del self.channel
        
    
        
    @legacy_function
    def _stop_worker():
        function = LegacyFunctionSpecification()  
        function.id = 0
        return function
        
    
        
    @legacy_function
    def internal__get_message_polling_interval():
        """Gets the message polling interval for MPI 
        header messages, in microseconds"""
        
        function = LegacyFunctionSpecification() 
        function.addParameter('polling_interval', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def internal__set_message_polling_interval():
        function = LegacyFunctionSpecification()  
        function.addParameter('polling_interval', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
        
        
    def stop(self):
        self._stop()
    
    def get_data_directory(self):
        """
        Returns the root name of the directory for the 
        application data files
        """
        pass
    
    def get_output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        pass
        


    @option(choices=['mpi','remote','distributed', 'sockets'], sections=("channel",))
    def channel_type(self):
        return 'mpi'
        
    @option(choices=("none","null","file"), sections=("channel",))
    def redirection(self):
        """Redirect the output of the code to null, standard streams or file"""
        return "null"
        
    @late
    def redirection_filenames(self):
        return {
            "none":("none", "none"),
            "null":("/dev/null", "/dev/null"), 
            "file":(self.redirect_stdout_file, self.redirect_stderr_file),
        }[self.redirection]
        
    @option(sections=("channel",))
    def redirect_stdout_file(self):
        return self.redirect_file
        
    @option(sections=("channel",))
    def redirect_stderr_file(self):
        return self.redirect_file

    @option(sections=("channel",))
    def redirect_file(self):
        return "code.out"
        
    
    @option(type='boolean', sections=("channel",))
    def must_start_worker(self):
        return True
    
    @late
    def channel_factory(self):
        if self.channel_type == 'mpi':
            if  MpiChannel.is_supported():
                return MpiChannel
            else:   
                return SocketChannel
                
        elif self.channel_type == 'remote':
            return MultiprocessingMPIChannel
        elif self.channel_type == 'distributed':
            return DistributedChannel
        elif self.channel_type == 'sockets':
            return SocketChannel
        else:
            raise exceptions.AmuseException("Cannot create a channel with type {0!r}, type is not supported".format(self.channel_type))
    
    def before_get_parameter(self):
        """
        Called everytime just before a parameter is retrieved in using::
            instance.parameter.name
        """
        pass
        
    def before_set_parameter(self):
        """
        Called everytime just before a parameter is updated in using::
            instance.parameter.name = newvalue
        """
        pass
        
class PythonCodeInterface(CodeInterface):
    """
    Base class for codes having a python implementation
    
    :argument implementation_factory: Class of the python implementation
    """
    
    def __init__(self, implementation_factory = None, name_of_the_worker = None, **options):
        self.implementation_factory = implementation_factory
        
        CodeInterface.__init__(self, name_of_the_worker, **options)
    
    def _start(self, name_of_the_worker = 'worker_code', **options):
        if name_of_the_worker is None:
            if self.implementation_factory is None:
                raise exceptions.CodeException("Must provide the name of a worker script or the implementation_factory class")
            name_of_the_worker = self.make_executable_script_for(self.implementation_factory)
        
        if self.use_python_interpreter:
            CodeInterface._start(self, name_of_the_worker = name_of_the_worker, interpreter_executable = self.python_interpreter, **options)
        else:
            CodeInterface._start(self, name_of_the_worker = name_of_the_worker, **options)
        
    def _check_if_worker_is_up_to_date(self):
        pass
        
    def make_executable_script_for(self, implementation_factory):
        from amuse.rfi.tools.create_python_worker import CreateAPythonWorker
        
        x = CreateAPythonWorker()
        x.channel_type = self.channel_type
        x.interface_class = type(self)
        x.implementation_factory = implementation_factory
        if self.channel_factory.is_root():
            x.start()
        return x.worker_name
        
        
    @classmethod
    def new_executable_script_string_for(cls, implementation_factory, channel_type = 'mpi'):
        raise Exception("tracing use")
            
    
    @option(type='boolean', sections=("channel",))
    def use_python_interpreter(self):
        return False

        
    @option(type='string', sections=("channel",))
    def python_interpreter(self):
        return sys.executable

        


