import weakref
import atexit
import errno
import os
import sys
import logging
import pydoc
import traceback
import random
import sys
import warnings

import inspect
import functools

# from collections import OrderedDict

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
from amuse.rfi.async_request import DependentASyncRequest

try:
    from amuse import config
except ImportError as ex:
    class config(object):
        is_mpi_enabled = False

CODE_LOG = logging.getLogger("code")
if CODE_LOG.level == logging.NOTSET:
    CODE_LOG.setLevel(logging.WARN)

"""
This module implements the code to the define interfaces between python
code and C++ or Fortran codes. It provides the abstract base
class for all community codes.
"""

import numpy

from amuse.rfi.channel import LocalChannel

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
        if self.interface.async_request:
            try:
                self.interface.async_request.wait()
            except Exception,ex:
                warnings.warn("Ignored exception in async call: " + str(ex))
            
        dtype_to_values = self.converted_keyword_and_list_arguments( arguments_list, keyword_arguments)
        
        handle_as_array = self.must_handle_as_array(dtype_to_values)
        
        if not self.owner is None:
            CODE_LOG.info("start call '%s.%s'",self.owner.__name__, self.specification.name)
        
        call_id = random.randint(0, 1000)
        
        try:
            self.interface.channel.send_message(call_id, self.specification.id, dtype_to_arguments = dtype_to_values)
            
            dtype_to_result = self.interface.channel.recv_message(call_id, self.specification.id, handle_as_array)
        except Exception, ex:
            CODE_LOG.info("Exception when calling function '{0}', of code '{1}', exception was '{2}'".format(self.specification.name, type(self.interface).__name__, ex))
            raise exceptions.CodeException("Exception when calling function '{0}', of code '{1}', exception was '{2}'".format(self.specification.name, type(self.interface).__name__, ex))
        
        result = self.converted_results(dtype_to_result, handle_as_array)
        
        if not self.owner is None:
            CODE_LOG.info("end call '%s.%s'",self.owner.__name__, self.specification.name)
        
        return result
    
    def _async_request(self, *arguments_list, **keyword_arguments):
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
            result=self.converted_results(dtype_to_result, handle_as_array)
            return result
            
        request.add_result_handler(handle_result)
        
        return request

    def asynchronous(self, *arguments_list, **keyword_arguments):
        if self.interface.async_request is not None:
            def factory():
              return self._async_request(*arguments_list, **keyword_arguments)
            request=DependentASyncRequest( self.interface.async_request, factory) 
        else:
            request=self._async_request(*arguments_list, **keyword_arguments)

        request._result_index=self.result_index()

        def handle_result(function):

            result=function()
            if  self.interface.async_request==request:
                self.interface.async_request=None
            return result

        request.add_result_handler(handle_result)

        self.interface.async_request=request

        return request
    
    def must_handle_as_array(self, keyword_arguments):
        for argument_type, argument_values in keyword_arguments.items():
            if argument_values:
                count = 0
                for argument_value in argument_values:
                    try:
                        if not isinstance(argument_value, basestring):
                            count = max(count, len(argument_value))
                    except:
                        count = max(count, 0)
                if count > 0:
                    return True
        return False
    
    """
    Get list of result keys
    """
    def result_index(self):
        index=[]
        for parameter in self.specification.output_parameters:
            index.append(parameter.name)
        if not self.specification.result_type is None:
            index.append("__result")
        return index
        
    """
    Convert results from an MPI message to a return value.
    """
    def converted_results(self, dtype_to_result, must_handle_as_array):
        
        number_of_outputs = len(self.specification.output_parameters)
        result_type = self.specification.result_type
        if number_of_outputs == 0:
            if result_type is None:
                return None
            return dtype_to_result[result_type][0]
            
        if number_of_outputs == 1 \
            and result_type is None:
            
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
        
        if not result_type is None:
            return_value =  dtype_to_array[result_type].pop()
        
        for parameter in self.specification.output_parameters:
            result[parameter.name] = dtype_to_array[parameter.datatype].pop()
        
        if not result_type is None:
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
        if self.specification.has_units:
            return CodeFunctionWithUnits(instance, owner, self.specification)
        else:
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

def derive_dtype_unit_and_default(value):
    if value is None:
        return None,None,None
    try: 
        unit=value.unit
        number=value.number
    except:
        unit=None
        number=value
    try:
        dtype=number.dtype.__str__()
        default=number
    except:
        if number in [ 'd','float64','i','int32','f','float32',
            's','string','b','bool','l','int64']:
            dtype=number
            default=None
        else:
          if isinstance(number,type):
              number=number()
              default=None
          else:
              default=number
          if isinstance(number, bool):
              dtype="b"
          elif isinstance(number,(int,long)):
              dtype="i"
          elif isinstance(number,(float,)):
              dtype="d"
          elif isinstance(number,(str,unicode)):
              dtype="s"
          else:
              raise Exception("undetectable type")
    return dtype,unit,default


def get_function_specification(name,in_arg,out_arg,must_handle_array=False,
                                   can_handle_array=False,length_arguments=None):
    function=LegacyFunctionSpecification()
    function.name=name
    function.must_handle_array=must_handle_array
    function.can_handle_array=can_handle_array
    if "__result" in out_arg:
        result=out_arg.pop("__result")
        dtype,unit,dummy=derive_dtype_unit_and_default(result)
        function.result_type=dtype
        function.result_unit=unit
    else:
        function.result_type='i'
        function.result_unit=None            
    inout_arg=dict()
    for arg in in_arg.keys():
        if arg in out_arg:
            inout_arg[arg]=in_arg.pop(arg)
            out_arg.pop(arg)     
    for arg,value in in_arg.items():
        dtype,unit,default=derive_dtype_unit_and_default(value)
        function.addParameter(arg, dtype=dtype, direction=function.IN ,unit=unit,default=default)
    for arg,value in inout_arg.items():
        dtype,unit,default=derive_dtype_unit_and_default(value)
        function.addParameter(arg, dtype=dtype, direction=function.INOUT ,unit=unit,default=default)
    for arg,value in out_arg.items():
        dtype,unit,default=derive_dtype_unit_and_default(value)
        function.addParameter(arg, dtype=dtype, direction=function.OUT ,unit=unit,default=default)
    if function.must_handle_array:
        if length_arguments:
            name=length_arguments[0]
        else:
            name="N"
        function.addParameter(name, dtype='i', direction=function.LENGTH)
    return function
  
def simplified_function_specification(must_handle_array=False,can_handle_array=False):
    def wrapper(f):
        argspec=inspect.getargspec(f)
        nkw=len(argspec.defaults) if argspec.defaults else 0
        defaults=argspec.defaults if argspec.defaults else []
        length_arguments=argspec.args[0:-nkw]
        kwargs=argspec.args[-nkw:]
        in_arg=OrderedDictionary()
        for x,y in zip(kwargs,defaults):
          in_arg[x]=y
        
        out_arg=[]
        flatsrc=inspect.getsource(f).replace("\n","").replace(" ","")
        def returns(**kwargs):
            start=flatsrc.find("returns(")
            order=lambda k: flatsrc.find(k[0]+"=",start)
            out_arg.extend(sorted(kwargs.items(),key=order))
        f.func_globals['returns']=returns
        f(*argspec.args)
        out_arg_mapping=OrderedDictionary()
        for x in out_arg:
            out_arg_mapping[x[0]] = x[1]
            
        function=get_function_specification(f.func_name,in_arg,out_arg_mapping,
            must_handle_array,can_handle_array,length_arguments)
        def g():
            return function
        return g
    return wrapper

def remote_function(f=None,must_handle_array=False,can_handle_array=False):
    # If called without method, we've been called with optional arguments.
    # We return a decorator with the optional arguments filled in.
    # Next time round we'll be decorating method.
    if f is None:
        return functools.partial(remote_function,must_handle_array=must_handle_array,can_handle_array=can_handle_array)
    return legacy_function(simplified_function_specification(must_handle_array=must_handle_array,can_handle_array=can_handle_array)(f))

     
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
        self.has_units = False
        self.result_doc = ''

    def set_name(self, name):
        self.name = name
        
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
        parameter.index_in_input = len(self.input_parameters) - 1
        parameters = self.dtype_to_input_parameters.get(parameter.datatype, [])
        parameters.append(parameter)
        parameter.input_index = len(parameters) - 1
        self.dtype_to_input_parameters[parameter.datatype] = parameters
   
    def add_output_parameter(self, parameter):
        self.output_parameters.append(parameter)
        parameter.index_in_output = len(self.output_parameters) - 1
        
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
                
    result_type = property(_get_result_type, _set_result_type)



def stop_interfaces(exceptions = []):
    """
    Stop the workers of all instantiated interfaces.
    
    All instantiated interfaces will become unstable
    after this call!
    """
    for reference in reversed(CodeInterface.instances):
        x = reference()
        if not x is None and x.__class__.__name__ not in exceptions:
            try:
                x._stop()
            except:
                pass
    for x in CodeInterface.classes:
        x.stop_reusable_channels()

class CodeInterface(OptionalAttributes):
    """
    Abstract base class for all interfaces to legacy codes.
    
    When a subclass is instantiated, a number of subprocesses
    will be started. These subprocesses are called workers
    as they implement the interface and do the actual work
    of the instantiated object.
    """
    instances = []
    classes = set([])
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
        
        self.async_request=None

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
        if self.reuse_worker:
            channel = self.retrieve_reusable_channel()
            if channel is not None:
                self.channel = channel
                return

        if interpreter_executable is None and self.use_interpreter:
            interpreter_executable = self.interpreter

        self.channel = self.channel_factory(name_of_the_worker, type(self), interpreter_executable = interpreter_executable, **options)
        
        self._check_if_worker_is_up_to_date()

        self.channel.redirect_stdout_file = self.redirection_filenames[0]
        self.channel.redirect_stderr_file = self.redirection_filenames[1]
        self.channel.polling_interval_in_milliseconds = self.polling_interval_in_milliseconds
        self.channel.initialize_mpi = self.initialize_mpi
        
        self.channel.start()
        
        # must register stop interfaces after channel start
        # if done before, the mpi atexit will be registered 
        # incorrectly
        self.ensure_stop_interface_at_exit()
        
        if self.channel.is_polling_supported():
            if self.polling_interval_in_milliseconds > 0:
                self.internal__set_message_polling_interval(int(self.polling_interval_in_milliseconds * 1000))
        
    def wait(self):
        if self.async_request is not None:
            self.async_request.wait()

    @option(type="int", sections=("channel",))
    def polling_interval_in_milliseconds(self):
        return 0
        
    @classmethod
    def ensure_stop_interface_at_exit(cls):
        if not cls.is_stop_interfaces_registered:
            atexit.register(stop_interfaces)
            cls.is_stop_interfaces_registered = True

    @classmethod
    def retrieve_reusable_channel(cls):
        if not 'REUSE_INSTANCE' in cls.__dict__:
            cls.REUSE_INSTANCE = set([])
        s = cls.REUSE_INSTANCE
        if len(s) > 0:
            return s.pop()
        else:
            return None

    @classmethod
    def store_reusable_channel(cls, instance):
        if not 'REUSE_INSTANCE' in cls.__dict__:
            cls.REUSE_INSTANCE = set([])
        s = cls.REUSE_INSTANCE
        s.add(instance)
        cls.classes.add(cls)
       
        
    @classmethod
    def stop_reusable_channels(cls):
        if not 'REUSE_INSTANCE' in cls.__dict__:
            cls.REUSE_INSTANCE = set([])
        s = cls.REUSE_INSTANCE
        while len(s) > 0:
            x = s.pop()
            call_id = random.randint(0, 1000)
            # do the _stop_worker call with low level send
            # (id == 0, no arguments)
            x.send_message(call_id, 0 , dtype_to_arguments = {})
            dtype_to_result = x.recv_message(call_id, 0, False)
            x.stop()

    def _stop(self):
        if hasattr(self, 'channel'):
            if not self.channel is None and self.channel.is_active():
                if self.reuse_worker:
                    self.store_reusable_channel(self.channel)
                    self.channel = None
                else:
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
        
    @legacy_function
    def internal__open_port():
        function = LegacyFunctionSpecification()  
        function.addParameter('port_identifier', dtype='string', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def internal__accept_on_port():
        function = LegacyFunctionSpecification()  
        function.addParameter('port_identifier', dtype='string', direction=function.IN)
        function.addParameter('comm_identifier', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def internal__connect_to_port():
        function = LegacyFunctionSpecification()  
        function.addParameter('port_identifier', dtype='string', direction=function.IN)
        function.addParameter('comm_identifier', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def internal__activate_communicator():
        function = LegacyFunctionSpecification()  
        function.addParameter('comm_identifier', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
        
        
    def stop(self):
        self._stop()
    
    @option(choices=['mpi','remote','distributed', 'sockets', 'local'], sections=("channel",))
    def channel_type(self):
        return 'mpi'
    


    @option(type="boolean", sections=("channel",))
    def initialize_mpi(self):
        """Is MPI initialized in the code or not. Defaults to True if MPI is available"""
        return config.mpi.is_enabled
        
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
        elif self.channel_type == 'local':
            return LocalChannel
        else:
            raise exceptions.AmuseException("Cannot create a channel with type {0!r}, type is not supported".format(self.channel_type))

    @option(type="boolean", sections=("channel",))
    def reuse_worker(self):
        """Do not stop a worker, re-use an existing one"""
        return False
    

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
        
    def before_set_interface_parameter(self):
        """
        Called everytime just before a interface parameter is updated in using::
            instance.parameter.name = newvalue
        """
        pass

    def before_new_set_instance(self):
        """
        (Can be) called everytime just before a new set is created
        """
        pass    

    def before_get_data_store_names(self):
        """
        called before getting data store names (for state model) - should eventually 
        not be necessary
        """
        pass    

    @option(type='string', sections=("channel",))
    def interpreter(self):
        return sys.executable




    @option(type='boolean', sections=("channel",))
    def use_interpreter(self):
        return False

        

    @legacy_function
    def internal__become_code():
        function = LegacyFunctionSpecification()                      
        function.addParameter('number_of_workers', dtype='int32', direction=function.IN)
        function.addParameter('modulename', dtype='string', direction=function.IN)
        function.addParameter('classname', dtype='string', direction=function.IN)
        function.result_type = 'int32'
        return function
        


class CodeWithDataDirectories(object):
    
    
    def __init__(self):
        if not self.channel_type == 'distributed':
            self.ensure_data_directory_exists(self.get_output_directory())
    
    def ensure_data_directory_exists(self, directory):
        directory = os.path.expanduser(directory)
        directory = os.path.expandvars(directory)
        
        try:
            os.makedirs(directory)
        except OSError as ex:
            if ex.errno == errno.EEXIST and os.path.isdir(directory):
                pass
            else:
                raise ex

    @property
    def module_name(self):
        return self.__module__.split('.')[-2]
    
    @property
    def data_directory(self):
        return self.get_data_directory()
    
    @property
    def output_directory(self):
        return self.get_output_directory()
    
    def get_data_directory(self):
        """
        Returns the root name of the directory for the 
        application data files.
        """
        return os.path.join(self.input_data_root_directory, self.module_name, 'input')
    
    def get_output_directory(self):
        """
        Returns the root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(self.output_data_root_directory, self.module_name, 'output')
    
    @option(type="string", sections=('data',))
    def amuse_root_directory(self):
        """
        The root directory of AMUSE, used as default root for all data directories
        """
        return self.channel.get_amuse_root_directory()
        
    @option(type="string", sections=('data',))
    def input_data_root_directory(self):
        """
        The root directory of the input data, read only directories
        """
        return os.path.join(self.amuse_root_directory, 'data')
        
    @option(type="string", sections=('data',))
    def output_data_root_directory(self):
        """
        The root directory of the output data,
        read - write directory
        """
        return os.path.join(self.amuse_root_directory, 'data')
    
    def get_code_src_directory(self):
        """
        Returns the root name of the application's source code directory.
        """
        return os.path.join(self.amuse_root_directory, 'src', 'amuse', 'community', self.module_name, 'src')
        
class PythonCodeInterface(CodeInterface):
    """
    Base class for codes having a python implementation
    
    :argument implementation_factory: Class of the python implementation
    """
    
    def __init__(self, implementation_factory = None, name_of_the_worker = None, **options):
        if self.channel_type == 'distributed':
            print "Warning! Distributed channel not fully supported by PythonCodeInterface yet"
        self.implementation_factory = implementation_factory
        self.worker_dir=options.get("worker_dir",None)
        
        CodeInterface.__init__(self, name_of_the_worker, **options)        
    
    def _start(self, name_of_the_worker = 'worker_code', **options):

        if name_of_the_worker is None:
            if self.implementation_factory is None:
                raise exceptions.CodeException("Must provide the name of a worker script or the implementation_factory class")
            name_of_the_worker = self.make_executable_script_for(self.implementation_factory)
            if not options.setdefault("dynamic_python_code",True):
                raise exceptions.CodeException("dynamic code set to false, but python code generated")
        
        if self.use_python_interpreter:
            CodeInterface._start(self, name_of_the_worker = name_of_the_worker, interpreter_executable = self.python_interpreter, **options)
        else:
            CodeInterface._start(self, name_of_the_worker = name_of_the_worker, **options)

    def _check_if_worker_is_up_to_date(self):
        pass
        
    def make_executable_script_for(self, implementation_factory):
        from amuse.rfi.tools.create_python_worker import CreateAPythonWorker
        
        x = CreateAPythonWorker()
        if self.worker_dir:
            x.worker_dir=self.worker_dir
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



class CodeFunctionWithUnits(CodeFunction):
   
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
        self.number_of_output_parameters = len(self.specification.output_parameters)
    
    def __call__(self, *arguments_list, **keyword_arguments):
        dtype_to_values, units = self.converted_keyword_and_list_arguments( arguments_list, keyword_arguments)
        encoded_units = self.convert_input_units_to_floats(units)
        
        handle_as_array = self.must_handle_as_array(dtype_to_values)
        
        if not self.owner is None:
            CODE_LOG.info("start call '%s.%s'",self.owner.__name__, self.specification.name)
        
        call_id = random.randint(0, 1000)
        try:
            self.interface.channel.send_message(call_id, self.specification.id, dtype_to_arguments = dtype_to_values, encoded_units = encoded_units)
            
            dtype_to_result , output_encoded_units = self.interface.channel.recv_message(call_id, self.specification.id, handle_as_array, has_units = True)
        except Exception, ex:
            CODE_LOG.info("Exception when calling function '{0}', of code '{1}', exception was '{2}'".format(self.specification.name, type(self.interface).__name__, ex))
            raise exceptions.CodeException("Exception when calling function '{0}', of code '{1}', exception was '{2}'".format(self.specification.name, type(self.interface).__name__, ex))
        
        output_units = self.convert_floats_to_units(output_encoded_units)
        result = self.converted_results(dtype_to_result, handle_as_array, output_units)
        
        if not self.owner is None:
            CODE_LOG.info("end call '%s.%s'",self.owner.__name__, self.specification.name)
        
        return result
    
    def _async_request(self, *arguments_list, **keyword_arguments):
        dtype_to_values, units = self.converted_keyword_and_list_arguments( arguments_list, keyword_arguments)
        encoded_units = self.convert_input_units_to_floats(units)
        
        handle_as_array = self.must_handle_as_array(dtype_to_values)
        
        call_id = random.randint(0, 1000)
              
        self.interface.channel.send_message(call_id, self.specification.id, dtype_to_arguments = dtype_to_values, encoded_units = encoded_units)
        
        request = self.interface.channel.nonblocking_recv_message(call_id, self.specification.id, handle_as_array, has_units = True)
        
        def handle_result(function):
            try:
                dtype_to_result, output_encoded_units = function()
            except Exception, ex:
                raise exceptions.CodeException("Exception when calling legacy code '{0}', exception was '{1}'".format(self.specification.name, ex))
            output_units = self.convert_floats_to_units(output_encoded_units)
            return self.converted_results(dtype_to_result, handle_as_array, output_units)
            
        request.add_result_handler(handle_result)
        return request
        
        
    

    def must_handle_as_array(self, keyword_arguments):
        for argument_type, argument_values in keyword_arguments.items():
            if argument_values:
                count = 0
                for argument_value in argument_values:
                    try:
                        if not isinstance(argument_value, basestring):
                            count = max(count, len(argument_value))
                    except:
                        count = max(count, 0)
                if count > 0:
                    return True
        return False
        
    """
    Convert results from an MPI message to a return value.
    """
    def converted_results(self, dtype_to_result, must_handle_as_array, units):
        
        number_of_outputs = self.number_of_output_parameters
        result_type = self.specification.result_type         
        
        if number_of_outputs == 0:
            if result_type is None:
                return None
            return dtype_to_result[result_type][0]
            
        if number_of_outputs == 1 \
            and result_type is None:
            
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
        
        if not result_type is None:
            return_value =  dtype_to_array[result_type].pop()
        
        for parameter in self.specification.output_parameters:
            result[parameter.name] = dtype_to_array[parameter.datatype].pop()
            if self.specification.has_units and not units[parameter.index_in_output] is None:
                result[parameter.name] = result[parameter.name] | units[parameter.index_in_output]
                
        
        if not result_type is None:
            result["__result"] =  return_value
        
        return result
       
    """
    Convert keyword arguments and list arguments to an MPI message
    """
    def converted_keyword_and_list_arguments(self, arguments_list, keyword_arguments):
        from  amuse.units import quantities
        dtype_to_values = self.specification.new_dtype_to_values()
        units = [None] * len(self.specification.input_parameters)
        
        input_parameters_seen = set(map(lambda x : x.name, self.specification.input_parameters))
        names_in_argument_list = set([])
        for index, argument in enumerate(arguments_list):
            parameter = self.specification.input_parameters[index]
            names_in_argument_list.add(parameter.name)
            
            if quantities.is_quantity(argument):
                units[parameter.index_in_input] = argument.unit
                argument = argument.number
                
            values = dtype_to_values[parameter.datatype]
            values[parameter.input_index] = argument
            input_parameters_seen.remove(parameter.name)
        
        for index, parameter in enumerate(self.specification.input_parameters):
            if parameter.name in keyword_arguments:
                argument = keyword_arguments[parameter.name]
                if quantities.is_quantity(argument):
                    units[parameter.index_in_input] = argument.unit
                    argument = argument.number
                
                values = dtype_to_values[parameter.datatype]
                values[parameter.input_index] = argument
                input_parameters_seen.remove(parameter.name)
        
        for parameter in self.specification.input_parameters:
            if (parameter.name in input_parameters_seen) and parameter.has_default_value():
                argument = parameter.default
                if quantities.is_quantity(argument):
                    units[parameter.index_in_input] = argument.unit
                    argument = argument.number
                    
                values = dtype_to_values[parameter.datatype]
                values[parameter.input_index] = argument
                input_parameters_seen.remove(parameter.name)
                
        if input_parameters_seen:
            raise exceptions.CodeException("Not enough parameters in call, missing " + str(sorted(input_parameters_seen)))
         
        return dtype_to_values, units
        
    def __str__(self):
        return str(self.specification)
        
        
    def convert_unit_to_floats(self, unit):
        if unit is None:
            return numpy.zeros(9, dtype=numpy.float64)
        else:
            return unit.to_array_of_floats()
        

    def convert_input_units_to_floats(self, units):
        result = numpy.zeros(len(units) * 9, dtype = numpy.float64)
        for index, unit in enumerate(units):
            offset = index*9
            result[offset:offset+9] = self.convert_unit_to_floats(unit)
        return result
        
    def convert_floats_to_units(self, floats):
        result = []
        for index in range(len(floats) // 9):
            offset = index*9
            unit_floats = floats[offset:offset+9]
            unit = self.convert_float_to_unit(unit_floats)
            result.append(unit)
        return result
        

    def convert_float_to_unit(self, floats):
        from amuse.units import core
        from amuse.units import units
        if numpy.all(floats == 0):
            return None
        factor = floats[0]
        result = factor
        system_index = floats[1]
        unit_system = None
        for x in core.system.ALL.values():
            if x.index == system_index:
                unit_system = x
                break
        for x in unit_system.bases:
            power = floats[x.index + 2]
            if not power == 0.0:
                result = result * (x ** power)
        return result
        






