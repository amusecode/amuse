"""
This module implements the code to the define interfaces between python
code and C++ or Fortran legacy codes. It provides the abstract base
class for all legacy codes.
"""

import weakref
import atexit
import os.path


from zlib import crc32
from mpi4py import MPI
from subprocess import Popen, PIPE

from amuse.support.core import late
from amuse.support.core import print_out
from amuse.support.core import OrderedDictionary
from amuse.legacy.support.create_definition import LegacyDocStringProperty
from amuse.legacy.support.channel import MpiChannel, MultiprocessingMPIChannel


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
        
def ensure_mpd_is_running():
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
    }
    if typecode in mapping:
        return mapping[typecode]
    
    values = mapping.values()
    if typecode in values:
        return typecode
    
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
        
        handle_as_array = False
        for x in keyword_arguments_for_the_mpi_channel.values():
           if x and hasattr(x[0],"__len__"):
              handle_as_array = len(x[0]) > 0
              break
                
        self.interface.channel.send_message(self.specification.id , **keyword_arguments_for_the_mpi_channel)
        (doubles, ints, floats, strings) = self.interface.channel.recv_message(self.specification.id, handle_as_array)
       
            
        return self.converted_results(doubles, ints, floats, strings)
        
    """
    Convert results from an MPI message to a return value.
    """
    def converted_results(self, doubles, ints, floats, strings):
        
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
            if self.specification.result_type == 'string':
                return strings[0] 
        
        if number_of_outputs == 1 \
            and self.specification.result_type is None:
            if len(ints) == 1:
                return ints[0]
            if len(doubles) == 1:
                return doubles[0]
            if len(floats) == 1:
                return floats[0]
            if len(strings) == 1:
                return strings[0]
        
        result = OrderedDictionary()
        dtype_to_array = {
            'float64' : list(reversed(doubles)),
            'int32' : list(reversed(ints)),
            'float32' : list(reversed(floats)),
            'string' : list(reversed(strings)),
        }
        
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
        
        if input_parameters_seen:
            raise Exception("Not enough parameters in call, missing " + str(sorted(input_parameters_seen)))
            
        dtype_to_keyword = {
            'float64' : 'doubles_in',
            'float32' : 'floats_in',
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
        >>> x = LegacyExample()
        >>> x.evolve.specification #doctest: +ELLIPSIS
        <amuse.legacy.support.core.LegacyFunctionSpecification object at 0x...>
        >>> LegacyExample.evolve #doctest: +ELLIPSIS
        <amuse.legacy.support.core.legacy_function object at 0x...>
        >>> x.evolve #doctest: +ELLIPSIS
        <amuse.legacy.support.core.LegacyCall object at 0x...>
        
                    
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
    
    def is_compiled_file_up_to_date(self, time_of_the_compiled_file):
        name_of_defining_file = self.specification_function.func_code.co_filename
        if os.path.exists(name_of_defining_file):
            time_of_defining_file = os.stat(name_of_defining_file).st_mtime
            return time_of_defining_file <= time_of_the_compiled_file
        return True
        

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


def stop_interfaces():
    """
    Stop the workers of all instantiated interfaces.
    
    All instantiated interfaces will become unstable
    after this call!
    """
    for reference in LegacyInterface.instances:
        x = reference()
        if not x is None:
            x._stop()

atexit.register(stop_interfaces)

class LegacyInterface(object):
    """
    Abstract base class for all interfaces to legacy codes.
    
    When a subclass is instantiated, a number of subprocesses
    will be started. These subprocesses are called workers
    as they implement the interface and do the actual work
    of the instantiated object.
    """
    instances = []
    channel_factory = MpiChannel
    
    def __init__(self, name_of_the_worker = 'worker_code', number_of_workers = 1, debug_with_gdb = False, hostname = None):
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
        self.channel = self.channel_factory(name_of_the_worker, number_of_workers, type(self), debug_with_gdb, hostname)
        self._check_if_worker_is_up_to_date()
        self.channel.start()
        self.instances.append(weakref.ref(self))
        
   
        
    def __del__(self):
        self._stop()
    
    def _stop(self):
        if hasattr(self, 'channel'):
            if self.channel.is_active():
                self._stop_worker()
                self.channel.stop()
            del self.channel
        
    def _check_if_worker_is_up_to_date(self):
        name_of_the_compiled_file = self.channel.full_name_of_the_worker
        modificationtime_of_worker = os.stat(name_of_the_compiled_file).st_mtime
        my_class = type(self)
        for x in dir(my_class):
            if x.startswith('__'):
                continue
            value = getattr(my_class, x)
            if isinstance(value, legacy_function):
                is_up_to_date = value.is_compiled_file_up_to_date(modificationtime_of_worker)
                if not is_up_to_date:
                    raise Exception("""The worker code of the '{0}' interface class is not up to date.
Please do a 'make clean; make' in the root directory.
""".format(type(self).__name__))
        
    @legacy_function
    def _stop_worker():
        function = LegacyFunctionSpecification()  
        function.id = 0
        return function   
        

class LegacyPythonInterface(LegacyInterface):
    """
    Base class for codes having a python implementation
    
    :argument implementation_factory: Class of the python implementation
    """
    
    def __init__(self, implementation_factory = None, name_of_the_worker = None , number_of_workers = 1, debug_with_gdb = False, hostname = None):
        if name_of_the_worker is None:
            if implementation_factory is None:
                raise Exception("Must provide the name of a worker script or the implementation_factory class")
            name_of_the_worker = self.make_executable_script_for(implementation_factory)
        
        LegacyInterface.__init__(self, name_of_the_worker, number_of_workers, False, hostname)
        
    def _check_if_worker_is_up_to_date(self):
        pass
        
    def make_executable_script_for(self, implementation_factory):
        import inspect
        import stat
        
        string = self.new_executable_script_string_for(implementation_factory)
        
        path = inspect.getfile(implementation_factory)
        path = os.path.dirname(path)
        path = os.path.join(path, 'worker-script')
        
        with open(path, 'w') as f:
            f.write(string)
        
        os.chmod(path, 0777)
        
        return path
        
    @classmethod
    def new_executable_script_string_for(cls, implementation_factory):
        import inspect
        import sys
        
        path = os.path.dirname(__file__)
        path = os.path.join(path, 'legacy_script.template')
        with open(path, "r") as f:
            template_string = f.read()
        

        return template_string.format(
            executable = sys.executable,
            syspath = ','.join(map(repr, sys.path)),
            factory_module = inspect.getmodule(implementation_factory).__name__,
            factory = implementation_factory.__name__,
            interface_module = inspect.getmodule(cls).__name__,
            interface = cls.__name__,
        )
            
            
            

        
        


