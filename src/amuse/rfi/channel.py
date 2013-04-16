import inspect
import numpy
import os.path
import cPickle as pickle


import sys
import struct
import threading
import select
import atexit
import time

import socket
import array

import logging

try:
    from mpi4py import rc
    rc.initialize = False
    from mpi4py import MPI
except ImportError:
    MPI = None
    
from subprocess import Popen, PIPE

from amuse.support.options import OptionalAttributes, option, GlobalOptions
from amuse.support.core import late
from amuse.support import exceptions
from amuse.rfi import run_command_redirected

class ASyncRequest(object):
        
    def __init__(self, request, message,  comm, header):
        self.request = request
        self.message = message
        self.header = header
        self.comm = comm
        self.is_finished = False
        self.is_set = False
        self._result = None
        self.result_handlers = []

    def wait(self):
        if self.is_finished:
            return
    
        self.request.Wait()
        self.is_finished = True
    
    def is_result_available(self):
        if self.is_finished:
            return True
            
        self.is_finished = self.request.Test()
        return self.is_finished
        
    def add_result_handler(self, function):
        self.result_handlers.append(function)
    
    def get_message(self):
        return self.message
        
    def _set_result(self):
        class CallingChain(object):
            def __init__(self, outer, inner):
                self.outer = outer
                self.inner = inner
                
            def __call__(self):
                return self.outer(self.inner)
                
        self.message.receive_content(self.comm, self.header)
        
        current = self.get_message
        for x in self.result_handlers:
            current = CallingChain(x, current)
        
        self._result = current()
        
        self.is_set = True
        
    def result(self):
        self.wait()
        
        if not self.is_set:
            self._set_result()
        
        return self._result
        
    def _new_handler(self, result_handler, args = (), kwargs = {}):
        return
        
    def is_mpi_request(self):
        return True


class ASyncSocketRequest(object):
        
    def __init__(self, message,  socket):
        self.message = message
        self.socket = socket
        
        self.is_finished = False
        self.is_set = False
        self._result = None
        self.result_handlers = []

    def wait(self):
        if self.is_finished:
            return
    
        while True:
            readables, _r, _x = select.select([self.socket],[],[])
            if len(readables) == 1:
                break
        
        self.is_finished = True
    
    def is_result_available(self):
        if self.is_finished:
            return True
            
        readables, _r, _x = select.select([self.socket],[],[], 0.001)
        
        self.is_finished = len(readables) == 1
        return self.is_finished
        
    def add_result_handler(self, function):
        self.result_handlers.append(function)
    
    def get_message(self):
        return self.message
        
    def _set_result(self):
        class CallingChain(object):
            def __init__(self, outer, inner):
                self.outer = outer
                self.inner = inner
                
            def __call__(self):
                return self.outer(self.inner)
                
        self.message.receive(self.socket)
        
        current = self.get_message
        for x in self.result_handlers:
            current = CallingChain(x, current)
        
        self._result = current()
        
        self.is_set = True
        
    def result(self):
        self.wait()
        
        if not self.is_set:
            self._set_result()
        
        return self._result
    
    def is_mpi_request(self):
        return False

        
        
class AsyncRequestWithHandler(object):
    
    def __init__(self, async_request,  result_handler, args = (), kwargs = {}):
        self.async_request = async_request
        self.result_handler = result_handler
        self.args = args
        self.kwargs = kwargs
    
    def run(self):
        self.result_handler(self.async_request, *self.args, **self.kwargs)
        
class AsyncRequestsPool(object):
    
    def __init__(self):
        self.requests_and_handlers = []
        self.registered_requests = set([])
        
    def add_request(self, async_request, result_handler, args = (), kwargs = {}):
        if async_request in self.registered_requests:
            raise Exception("Request is already registered, cannot register a request more than once")
            
        self.registered_requests.add(async_request)
        
        self.requests_and_handlers.append( 
            AsyncRequestWithHandler(
                async_request,
                result_handler,
                args,
                kwargs
            )
        )
    
        
    def wait(self):
        # TODO need to cleanup this code
        #
        requests = [x.async_request.request for x in self.requests_and_handlers if x.async_request.is_mpi_request()]
        indices = [i for i, x in enumerate(self.requests_and_handlers) if x.async_request.is_mpi_request()]
        
        if len(requests) > 0:
            index = MPI.Request.Waitany(requests)
              
            index = indices[index]
            
            request_and_handler = self.requests_and_handlers[index]
            
            self.registered_requests.remove(request_and_handler.async_request)
            
            del self.requests_and_handlers[index]
            
            request_and_handler.async_request.wait() #will set the finished flag
            
            request_and_handler.run()
            
        
        sockets = [x.async_request.socket for x in self.requests_and_handlers if not x.async_request.is_mpi_request()]
        indices = [i for i, x in enumerate(self.requests_and_handlers) if not x.async_request.is_mpi_request()]
        if len(sockets) > 0:
            readable, _, _ = select.select(sockets, [], [])
            indices_to_delete =[]
            for read_socket in readable:
                
                index = sockets.index(read_socket)
                
                index = indices[index]
                
                request_and_handler = self.requests_and_handlers[index]
                
                self.registered_requests.remove(request_and_handler.async_request)
                
                indices_to_delete.append(index)
                
                request_and_handler.async_request.wait() #will set the finished flag
                
                request_and_handler.run()
            
            for x in reversed(list(sorted(indices_to_delete))):
                
                del self.requests_and_handlers[x]
            
    def __len__(self):
        return len(self.requests_and_handlers)
        
        
class AbstractMessage(object):
    
    def __init__(self,
        call_id = 0, function_id = -1, call_count = 1, 
        dtype_to_arguments = {}, 
        error = False, 
        big_endian = (sys.byteorder.lower() == 'big'),
        polling_interval = 0
    ):
        self.polling_interval = polling_interval
        
        #flags
        self.big_endian = big_endian
        self.error = error
        
        #header
        self.call_id = call_id
        self.function_id = function_id
        self.call_count = call_count
        
        #data (numpy arrays)
        self.ints = []
        self.longs = []
        self.floats = []
        self.doubles = []
        self.strings = []
        self.booleans = []
        
        self.pack_data(dtype_to_arguments)
        
    def pack_data(self,dtype_to_arguments):
        for dtype, attrname in self.dtype_to_message_attribute():
            if dtype in dtype_to_arguments:
                array = pack_array( dtype_to_arguments[dtype], self.call_count, dtype)
                setattr(self, attrname, array)
    
    def to_result(self, handle_as_array = False):
        dtype_to_result = {}
        for dtype, attrname in self.dtype_to_message_attribute():
            result = getattr(self, attrname)
            if self.call_count > 1 or handle_as_array:
                dtype_to_result[dtype] = unpack_array(result , self.call_count, dtype)
            else:
                dtype_to_result[dtype] = result
                    
        return dtype_to_result
    
    def dtype_to_message_attribute(self):
        return (
            ('int32', 'ints'),
            ('int64', 'longs'),
            ('float32', 'floats'),
            ('float64', 'doubles'),
            ('bool', 'booleans'),
            ('string', 'strings'),
        )
    
    def receive(self, comm):
        raise NotImplementedError
        
    def send(self, comm):
        raise NotImplementedError
        
    
class MPIMessage(AbstractMessage):
        
    def receive(self, comm):
        header = self.receive_header(comm)
        self.receive_content(comm, header)
        
    def receive_header(self, comm):
        header = numpy.zeros(10,  dtype='i')
        self.mpi_receive(comm, [header, MPI.INT])
        return header
        
    def receive_content(self, comm, header):
        #4 flags as 8bit booleans in 1st 4 bytes of header
        #endiannes(not supported by MPI channel), error, unused, unused 
        flags = header.view(dtype='bool8')
        self.big_endian = flags[0]
        self.error = flags[1]

        self.call_id = header[1]
        self.function_id = header[2]
        self.call_count = header[3]

        number_of_ints = header[4]
        number_of_longs = header[5]
        number_of_floats = header[6]
        number_of_doubles = header[7]
        number_of_booleans = header[8]
        number_of_strings = header[9]
        
        self.ints = self.receive_ints(comm, number_of_ints)
        self.longs = self.receive_longs(comm, number_of_longs)
        self.floats = self.receive_floats(comm, number_of_floats)
        self.doubles = self.receive_doubles(comm, number_of_doubles)
        self.booleans = self.receive_booleans(comm, number_of_booleans)
        self.strings = self.receive_strings(comm, number_of_strings)
        
    def nonblocking_receive(self, comm):
        header = numpy.zeros(10,  dtype='i')
        request = self.mpi_nonblocking_receive(comm, [header, MPI.INT])
        return ASyncRequest(request, self, comm,  header)
    
    def receive_doubles(self, comm, total):
        if total > 0:
            result = numpy.empty(total, dtype='d')
            self.mpi_receive(comm,[result, MPI.DOUBLE])
            return result
        else:
            return []
            
    def receive_ints(self, comm, total):
        if total > 0:
            result = numpy.empty(total, dtype='i')
            self.mpi_receive(comm,[result, MPI.INT])
            return result
        else:
            return []
            
    def receive_longs(self, comm, total):
        if total > 0:
            result = numpy.empty(total, dtype='int64')
            self.mpi_receive(comm,[result, MPI.INTEGER8])
            return result
        else:
            return []
            
    def receive_floats(self, comm, total):
        if total > 0:
            result = numpy.empty(total, dtype='f')
            self.mpi_receive(comm,[result, MPI.FLOAT])
            return result
        else:
            return []
            
    
    def receive_booleans(self, comm, total):
        if total > 0:
            result = numpy.empty(total,  dtype='int32')
            self.mpi_receive(comm,[result, MPI.LOGICAL])
            return numpy.logical_not(result == 0)
        else:
            return []
    
            
    def receive_strings(self, comm, total):
        if total > 0:
            sizes = numpy.empty(total, dtype='i')
            
            self.mpi_receive(comm,[sizes, MPI.INT])
            
            logging.getLogger("channel").debug("got %d strings of size %s", total, sizes)
            
            byte_size = 0
            for size in sizes:
                byte_size = byte_size + size + 1
                
            data_bytes = numpy.empty(byte_size, dtype=numpy.uint8)
            self.mpi_receive(comm,[data_bytes,  MPI.CHARACTER])
            
            strings = []
            begin = 0
            for size in sizes:
                strings.append(data_bytes[begin:begin+size].tostring())
                begin = begin + size + 1
                
            logging.getLogger("channel").debug("got %d strings of size %s, data = %s", total, sizes, strings)
            return strings
        else:
            return []
        
    
    def send(self, comm):
        header = numpy.array([
            0,
            self.call_id,
            self.function_id, 
            self.call_count,
            len(self.ints) , 
            len(self.longs) ,
            len(self.floats) , 
            len(self.doubles) ,
            len(self.booleans) ,
            len(self.strings) ,
        ], dtype='i')
        
        
        flags = header.view(dtype='bool8')
        flags[0] = self.big_endian
        flags[1] = self.error
        self.send_header(comm, header)
        self.send_content(comm)
    
    def send_header(self, comm, header):
        self.mpi_send(comm,[header, MPI.INT])
        
    def send_content(self, comm):
        self.send_ints(comm, self.ints)
        self.send_longs(comm, self.longs)
        self.send_floats(comm, self.floats)
        self.send_doubles(comm, self.doubles)
        self.send_booleans(comm, self.booleans)
        self.send_strings(comm, self.strings)
        
    def send_ints(self, comm, array):
        if len(array) > 0:
            sendbuffer = numpy.array(array,  dtype='int32')
            self.mpi_send(comm,[sendbuffer, MPI.INT])
            
    def send_longs(self, comm, array):
        if len(array) > 0:
            sendbuffer = numpy.array(array,  dtype='int64')
            self.mpi_send(comm,[sendbuffer, MPI.INTEGER8])    
        
    def send_doubles(self, comm, array):
        if len(array) > 0:
            sendbuffer = numpy.array(array,  dtype='d')
            self.mpi_send(comm,[sendbuffer, MPI.DOUBLE])
            
    def send_floats(self, comm, array):
        if len(array) > 0:
            sendbuffer = numpy.array(array,  dtype='f')
            self.mpi_send(comm, [sendbuffer, MPI.FLOAT])
            
    def send_strings(self, comm, array):
        if len(array) == 0:
            return
            
        lengths = self.string_lengths(array)
        self.mpi_send(comm, [lengths, MPI.INT])
        chars=""
        for string in array:
            chars=string.join((chars,chr(0)))
            
        chars=numpy.fromstring(chars,dtype='uint8')
        self.mpi_send(comm, [chars, MPI.CHARACTER])
        
    def send_booleans(self, comm, array):
        if len(array) > 0:
            sendbuffer = numpy.array(array,  dtype='int32')
            self.mpi_send(comm, [sendbuffer, MPI.LOGICAL])
    
    def string_lengths(self, array):
        lengths = numpy.zeros(len(array), dtype='i')
        index = 0
        
        for string in array:
            lengths[index] = len(string)
            index += 1
        
        return lengths
    
    def set_error(self, message):
        self.strings = [message]
        self.error = True
        
    def mpi_nonblocking_receive(self, comm, array):
        raise NotImplementedError()
    
    def mpi_receive(self, comm, array):
        raise NotImplementedError()
        
    def mpi_send(self, comm, array):
        raise NotImplementedError()
    
    
class ServerSideMPIMessage(MPIMessage):
    
    def mpi_receive(self, comm, array):
        request = comm.Irecv(array,  source=0, tag=999)
        request.Wait()
        
    def mpi_send(self, comm, array):
        comm.Bcast(array, root=MPI.ROOT)
        
    def send_header(self, comm, array):
        requests = []
        for rank in range(comm.Get_remote_size()):
            request = comm.Isend(array, dest=rank, tag=989)
            requests.append(request)
        MPI.Request.Waitall(requests)
    
        
    def mpi_nonblocking_receive(self, comm, array):
        return comm.Irecv(array,  source=0, tag=999)

    def receive_header(self, comm):
        header = numpy.zeros(10,  dtype='i')
        request = comm.Irecv([header, MPI.INT],  source=0, tag=999)
        if self.polling_interval  > 0:
            is_finished = request.Test()
            while not is_finished:
                time.sleep(self.polling_interval / 1000000.)
                is_finished = request.Test()
            request.Wait()
        else:
            request.Wait()
            
        return header
        
    
class ClientSideMPIMessage(MPIMessage):
    
    def mpi_receive(self, comm, array):
        comm.Bcast(array, root = 0)
        
    def mpi_send(self, comm, array):
        comm.Send(array, dest=0, tag = 999)

    def mpi_nonblocking_receive(self, comm, array):
        return comm.Irecv(array,  source=0, tag=999)

    def receive_header(self, comm):
        header = numpy.zeros(10,  dtype='i')
        request = comm.Irecv([header, MPI.INT],  source=0, tag=989)
        if self.polling_interval  > 0:
            is_finished = request.Test()
            while not is_finished:
                time.sleep(self.polling_interval / 1000000.)
                is_finished = request.Test()
            request.Wait()
        else:
            request.Wait()
        return header

MAPPING = {}

def pack_array(array, length,  dtype):
    if dtype == 'string':
        if length == 1 and len(array) > 0 and isinstance(array[0], str):
            return array
        else:
            result = []
            for x in array:
                result.extend(x)
            return result
    else:
        total_length = length * len(array)
        if dtype in MAPPING:
            result = MAPPING.dtype
            if len(result) != total_length:
                result = numpy.empty(length * len(array), dtype = dtype)
        else:        
            result = numpy.empty(length * len(array), dtype = dtype)
        
        for i in range(len(array)):
            offset = i * length
            result[offset:offset+length] = array[i]
        return result
    

def unpack_array(array, length, dtype = None):
    result = []
    total = len(array) / length
    for i in range(total):
        offset = i * length
        result.append(array[offset:offset+length])
    return result

class AbstractMessageChannel(OptionalAttributes):
    """
    Abstract base class of all message channel.
    
    A message channel is used to send and retrieve messages from
    a remote party. A message channel can also setup the remote
    party. For example starting an instance of an application
    using MPI calls.
    
    The messages are encoded as arguments to the send and retrieve
    methods. Each message has an id and and optional list of doubles,
    integers, floats and/or strings.
    
    """
    
    def __init__(self,   **options):
        OptionalAttributes.__init__(self, **options)
    
    @classmethod
    def GDB(cls, full_name_of_the_worker, channel, interpreter_executable = None):
        arguments = ['-hold', '-display', os.environ['DISPLAY'], '-e', 'gdb', '--args']
        
        if not interpreter_executable is None:
            arguments.append(interpreter_executable)
        
        arguments.append(full_name_of_the_worker)
        
        command = 'xterm'
        return command, arguments
        
    @classmethod
    def DDD(cls, full_name_of_the_worker, channel, interpreter_executable = None):
        arguments = ['-display', os.environ['DISPLAY'], '-e', 'ddd', '--args']
        
        if not interpreter_executable is None:
            arguments.append(interpreter_executable)
        
        arguments.append(full_name_of_the_worker)
        
        command = 'xterm'
        return command, arguments
        
    @classmethod
    def VALGRIND(cls, full_name_of_the_worker, channel, interpreter_executable = None):
        #arguments = ['-hold', '-display', os.environ['DISPLAY'], '-e', 'valgrind',  full_name_of_the_worker]
        arguments = []
        
        if not interpreter_executable is None:
            arguments.append(interpreter_executable)
        
        arguments.append(full_name_of_the_worker)
        command = 'valgrind'
        return command, arguments
        
        
    @classmethod
    def XTERM(cls, full_name_of_the_worker, channel, interpreter_executable = None):
        arguments = ['-hold', '-display', os.environ['DISPLAY'], '-e']
        
        if not interpreter_executable is None:
            arguments.append(interpreter_executable)
        
        arguments.append(full_name_of_the_worker)
        
        command = 'xterm'
        return command, arguments
        

    @classmethod
    def REDIRECT(cls, full_name_of_the_worker, stdoutname, stderrname, command = None, interpreter_executable = None, **options):
        
        fname = run_command_redirected.__file__                
        arguments = [fname , stdoutname, stderrname]
        
        if not interpreter_executable is None:
            arguments.append(interpreter_executable)
            
        arguments.append(full_name_of_the_worker)
        
        if command is None :
            command = sys.executable
        
            
        return command, arguments
    
    @classmethod
    def GDBR(cls, full_name_of_the_worker, channel, interpreter_executable = None):
        "remote gdb, can run without xterm"
        
        arguments = ['localhost:{0}'.format(channel.debugger_port)]
        
        if not interpreter_executable is None:
            arguments.append(interpreter_executable)
            
        arguments.append(full_name_of_the_worker)
        
        command = 'gdbserver'
        return command, arguments
        
    @classmethod
    def NODEBUGGER(cls, full_name_of_the_worker, channel, interpreter_executable = None):
        if not interpreter_executable is None:
            return interpreter_executable, [full_name_of_the_worker]
        else:
            return full_name_of_the_worker, []
    
    @option(type='string', sections=("channel",))
    def worker_code_suffix(self):
        return ''
        
    @option(type='string', sections=("channel",))
    def worker_code_prefix(self):
        return ''
        
    @option(type='string', sections=("channel",))
    def worker_code_directory(self):
        return ''
        
    @option(type='boolean', sections=("channel",))
    def must_check_if_worker_is_up_to_date(self):
        return True
    
    @option(type="int", sections=("channel",))
    def number_of_workers(self):
        return 1
        
    def check_if_worker_is_up_to_date(self, object):
        if not self.must_check_if_worker_is_up_to_date:
            return
            
        name_of_the_compiled_file = self.full_name_of_the_worker
        modificationtime_of_worker = os.stat(name_of_the_compiled_file).st_mtime
        my_class = type(object)
        for x in dir(my_class):
            if x.startswith('__'):
                continue
            value = getattr(my_class, x)
            if hasattr(value, 'crc32'):
                is_up_to_date = value.is_compiled_file_up_to_date(modificationtime_of_worker)
                if not is_up_to_date:
                    raise exceptions.CodeException("""The worker code of the '{0}' interface class is not up to date.
Please do a 'make clean; make' in the root directory.
""".format(type(object).__name__))

    def get_full_name_of_the_worker(self, type):
        
        if os.path.isabs(self.name_of_the_worker):
            if os.path.exists(self.name_of_the_worker):
                if not os.access(self.name_of_the_worker, os.X_OK):
                    raise exceptions.CodeException("The worker application exists, but it is not executable.\n{0}".format(self.name_of_the_worker))
       
                return self.name_of_the_worker
        
        exe_name = self.worker_code_prefix + self.name_of_the_worker + self.worker_code_suffix
        
        tried_workers = []
        found = False
        
        if len(self.worker_code_directory) > 0 and os.path.exists(self.worker_code_directory):
            full_name_of_the_worker = os.path.join(self.worker_code_directory, exe_name)
            full_name_of_the_worker = os.path.normpath(os.path.abspath(full_name_of_the_worker))
            found = os.path.exists(full_name_of_the_worker)
            if not found:
                tried_workers.append(full_name_of_the_worker)
                
        current_type=type
        while not found:
            directory_of_this_module = os.path.dirname(inspect.getfile(current_type))
            full_name_of_the_worker = os.path.join(directory_of_this_module, exe_name)
            full_name_of_the_worker = os.path.normpath(os.path.abspath(full_name_of_the_worker))
            found = os.path.exists(full_name_of_the_worker)
            if not found:
                tried_workers.append(full_name_of_the_worker)
                current_type = current_type.__bases__[0]
                if current_type.__bases__[0] is object:
                    break
            else:
                found = True
        
        if not found:
            directory_of_this_module = os.path.dirname(os.path.dirname(__file__))
            full_name_of_the_worker = os.path.join(directory_of_this_module, '_workers', exe_name)
            full_name_of_the_worker = os.path.normpath(os.path.abspath(full_name_of_the_worker))
            
            found = os.path.exists(full_name_of_the_worker)
            if not found:
                raise exceptions.CodeException("The worker application does not exists, it should be at: \n{0}".format('\n'.join(tried_workers)))
            else:
                found = True
            
        return full_name_of_the_worker
    
    def send_message(self, call_id = 0, function_id = -1, dtype_to_arguments = {}):
        pass
        
    def recv_message(self, call_id = 0, function_id = -1, handle_as_array = False):
        pass
    
    def nonblocking_recv_message(self, call_id = 0, function_id = -1, handle_as_array = False):
        pass
        
    def start(self):
        pass
        
    def stop(self):
        pass

    def is_active(self):
        return True
        
    @classmethod
    def is_root(self):
        return True
    
    def is_polling_supported(self):
        return False
 
AbstractMessageChannel.DEBUGGERS = {
    "none":None,
    "gdb":AbstractMessageChannel.GDB, 
    "ddd":AbstractMessageChannel.DDD, 
    "xterm":AbstractMessageChannel.XTERM,
    "gdb-remote":AbstractMessageChannel.GDBR, 
    "valgrind":AbstractMessageChannel.VALGRIND,
}

#import time
#import ctypes
#clib_library = ctypes.CDLL("libc.so.6")
#memcpy = clib_library.memcpy
#memcpy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t]


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
    if not MpiChannel.is_supported():
        return True
        
    name_of_the_vendor, version = MPI.get_vendor()
    if name_of_the_vendor == 'MPICH2':
        must_check_mpd = True
        if 'AMUSE_MPD_CHECK' in os.environ:
            must_check_mpd = os.environ['AMUSE_MPD_CHECK'] == '1'
        if 'PMI_PORT' in os.environ:
            must_check_mpd = False
        if 'PMI_RANK' in os.environ:
            must_check_mpd = False
        if 'HYDRA_CONTROL_FD' in os.environ:
            must_check_mpd = False
        
        if not must_check_mpd:
            return True
        try:
            process = Popen(['mpdtrace'], stdout = PIPE, stderr = PIPE)
            (output_string, error_string) = process.communicate()
            return not (process.returncode == 255)
        except OSError as ex:
            return True
    else:
        return True


class MpiChannel(AbstractMessageChannel):
    """
    Message channel based on MPI calls to send and recv the messages
    
    :argument name_of_the_worker: Name of the application to start
    :argument number_of_workers: Number of parallel processes
    :argument legacy_interface_type: Type of the legacy interface
    :argument debug_with_gdb: If True opens an xterm with a gdb to debug the remote process
    :argument hostname: Name of the node to run the application on
    """
    _mpi_is_broken_after_possible_code_crash = False
    _intercomms_to_disconnect = []
    _is_registered = False

    
    def __init__(self, name_of_the_worker, legacy_interface_type = None, interpreter_executable = None,  **options):
        AbstractMessageChannel.__init__(self, **options)
        
        logging.basicConfig(level=logging.WARN)
        #logging.getLogger("channel").setLevel(logging.DEBUG)
        #logging.getLogger("code").setLevel(logging.DEBUG)
        
        self.ensure_mpi_initialized()
        
        self.name_of_the_worker = name_of_the_worker
        self.interpreter_executable = interpreter_executable
                
        if not legacy_interface_type is None:
            self.full_name_of_the_worker = self.get_full_name_of_the_worker( legacy_interface_type)
        else:
            self.full_name_of_the_worker = self.name_of_the_worker
        
        if self.check_mpi:
            if not is_mpd_running():
                raise exceptions.CodeException("The mpd daemon is not running, please make sure it is started before starting this code")
        
        if self._mpi_is_broken_after_possible_code_crash:
            raise exceptions.CodeException("Another code has crashed, cannot spawn a new code, please stop the script and retry")
        
        if not self.hostname is None:
            self.info = MPI.Info.Create()
            self.info['host'] = self.hostname
        else:
            self.info = MPI.INFO_NULL
            
        self.cached = None
        self.intercomm = None
        self._is_inuse = False
        self._communicated_splitted_message = False
    
    
    @classmethod
    def ensure_mpi_initialized(cls):
        if not MPI.Is_initialized():
            if rc.threaded:
                MPI.Init_thread(MPI.THREAD_MULTIPLE)
            else:
                MPI.Init()
        cls.register_finalize_code()

    @classmethod
    def register_finalize_code(cls):
        if not cls._is_registered:
            atexit.register(cls.finialize_mpi_atexit)
            cls._is_registered = True
    
    @classmethod
    def finialize_mpi_atexit(cls):
        if not MPI.Is_initialized():
            return
        if MPI.Is_finalized():
            return
        try:
            for x in cls._intercomms_to_disconnect:
                x.Disconnect()
                
            MPI.Finalize()
        except MPI.Exception as ex:
            return
        
    @option(type="boolean", sections=("channel",))
    def check_mpi(self):
        return True
        
    @option(type="boolean", sections=("channel",))
    def debug_with_gdb(self):
        return False
        
    @option(sections=("channel",))
    def hostname(self):
        return None
    
        
    @option(type="int", sections=("channel",))
    def debugger_port(self):
        return 4343
        
    @option(sections=("channel",))
    def python_exe_for_redirection(self):
        return None
        
    @option(choices=AbstractMessageChannel.DEBUGGERS.keys(), sections=("channel",))
    def debugger(self):
        """Name of the debugger to use when starting the code"""
        return "none"
        
    
    @option(type="int", sections=("channel",))
    def max_message_length(self):
        """
        For calls to functions that can handle arrays, MPI messages may get too long for large N.
        The MPI channel will split long messages into blocks of size max_message_length.
        """ 
        return 1000000
    

    @late
    def redirect_stdout_file(self):
        return "/dev/null"
        
    @late
    def redirect_stderr_file(self):
        return "/dev/null"
        
    @late
    def debugger_method(self):
        return self.DEBUGGERS[self.debugger]
    
    @classmethod
    def is_supported(cls):
        return not MPI is None
    
    @option(type="boolean", sections=("channel",))
    def can_redirect_output(self):
        name_of_the_vendor, version = MPI.get_vendor()
        if name_of_the_vendor == 'MPICH2':
            if 'MPISPAWN_ARGV_0' in os.environ:
                return False
        return True
        
    
    @option(type="boolean", sections=("channel",))
    def must_disconnect_on_stop(self):
        name_of_the_vendor, version = MPI.get_vendor()
        if name_of_the_vendor == 'MPICH2':
            if 'MPISPAWN_ARGV_0' in os.environ:
                return False
        return True
    
    @option(type="int", sections=("channel",))
    def polling_interval_in_milliseconds(self):
        return 0
        
    @classmethod
    def is_root(self):
        return MPI.COMM_WORLD.rank == 0
        
    def start(self):
        if not self.debugger_method is None:
            command, arguments = self.debugger_method(self.full_name_of_the_worker, self, interpreter_executable = self.interpreter_executable)
        else:
            if not self.can_redirect_output or (self.redirect_stdout_file == 'none' and self.redirect_stderr_file == 'none'):
                
                if self.interpreter_executable is None:
                    command = self.full_name_of_the_worker
                    arguments = None
                else:
                    command = self.interpreter_executable
                    arguments = [self.full_name_of_the_worker]
            else:
                command, arguments = self.REDIRECT(self.full_name_of_the_worker, self.redirect_stdout_file, self.redirect_stderr_file, command = self.python_exe_for_redirection, interpreter_executable = self.interpreter_executable)
                
        #print arguments
        #print command
        self.intercomm = MPI.COMM_SELF.Spawn(command, arguments, self.number_of_workers, info = self.info)
            
        
        
    def stop(self):
        if not self.intercomm is None:
            try:
                if self.must_disconnect_on_stop:
                    self.intercomm.Disconnect()
                else:
                    self._intercomms_to_disconnect.append(self.intercomm)
            except MPI.Exception as ex:
                if ex.error_class == MPI.ERR_OTHER:
                    type(self)._mpi_is_broken_after_possible_code_crash = True
                
            self.intercomm = None
    
    def determine_length_from_data(self, dtype_to_arguments):
        def get_length(x):
            if x:
                try:
                    if not isinstance(x[0], str):
                        return len(x[0])
                except:
                    return 1
               
               
        
        lengths = map(get_length, dtype_to_arguments.values())
        if len(lengths) == 0:
            return 1
            
        return max(1, max(lengths))
        
        
    def send_message(self, call_id, function_id, dtype_to_arguments = {}):

        if self.is_inuse():
            raise exceptions.CodeException("You've tried to send a message to a code that is already handling a message, this is not correct")
        if self.intercomm is None:
            raise exceptions.CodeException("You've tried to send a message to a code that is not running")
        
        call_count = self.determine_length_from_data(dtype_to_arguments)
        
        if call_count > self.max_message_length:
            self.split_message(call_id, function_id, call_count, dtype_to_arguments)
            return
        
        message = ServerSideMPIMessage(
            call_id, function_id, 
            call_count, dtype_to_arguments
        )
        message.send(self.intercomm)
        
        self._is_inuse = True
        
    def split_message(self, call_id, function_id, call_count, dtype_to_arguments):
        
        def split_input_array(i, arr_in):
            if call_count == 1:
                return [tmp[i*self.max_message_length] for tmp in arr_in]
            else:
                result = []
                for x in arr_in:
                    if hasattr(x, '__iter__'):
                        result.append(x[i*self.max_message_length:(i+1)*self.max_message_length])
                    else:
                        result.append(x)
                return result
        
        dtype_to_result = {}
        
        for i in range(1+(call_count-1)/self.max_message_length):
            split_dtype_to_argument = {}
            for key, value in dtype_to_arguments.iteritems():
                split_dtype_to_argument[key] = split_input_array(i, value)
                
            self.send_message(
                call_id,
                function_id,
                #min(call_count,self.max_message_length),
                split_dtype_to_argument
            )
            
            partial_dtype_to_result = self.recv_message(call_id, function_id, True)
            for datatype, value in partial_dtype_to_result.iteritems():
                if not datatype in dtype_to_result:
                    dtype_to_result[datatype] = [] 
                    for j, element in enumerate(value):
                        if datatype == 'string':
                            dtype_to_result[datatype].append([])
                        else:
                            dtype_to_result[datatype].append(numpy.zeros((call_count,), dtype=datatype))
                            
                for j, element in enumerate(value):
                    if datatype == 'string':
                        dtype_to_result[datatype][j].extend(element)
                    else:
                        dtype_to_result[datatype][j][i*self.max_message_length:(i+1)*self.max_message_length] = element
                
            #print partial_dtype_to_result
            call_count -= self.max_message_length
        
        self._communicated_splitted_message = True
        self._merged_results_splitted_message = dtype_to_result
    
    def recv_message(self, call_id, function_id, handle_as_array):
        
        self._is_inuse = False
        
        if self._communicated_splitted_message:
            x = self._merged_results_splitted_message
            self._communicated_splitted_message = False
            del self._merged_results_splitted_message
            return x
        
        message = ServerSideMPIMessage(
            polling_interval = self.polling_interval_in_milliseconds * 1000
        )
        try:
            message.receive(self.intercomm)
        except MPI.Exception as ex:
            self.stop()
            raise ex
        
        if message.call_id != call_id:
            self.stop()
            raise exceptions.CodeException('Received reply for call id {0} but expected {1}'.format(message.call_id, call_id))
        if message.function_id != function_id:
            self.stop()
            raise exceptions.CodeException('Received reply for function id {0} but expected {1}'.format(message.function_id, function_id))
        
        if message.error:
                logging.getLogger("channel").info("error message!")
                raise exceptions.CodeException("Error in code: " + message.strings[0])
#        if message.tag == -1:
#            raise exceptions.CodeException("Not a valid message, message is not understood by legacy code")
#        elif message.tag == -2:
#            self.stop()
#            raise exceptions.CodeException("Fatal error in code, code has exited")
        
        return message.to_result(handle_as_array)
        
    def nonblocking_recv_message(self, call_id, function_id, handle_as_array):
        request = ServerSideMPIMessage().nonblocking_receive(self.intercomm)
        def handle_result(function):
            self._is_inuse = False
        
            message = function()
            
            if message.call_id != call_id:
                self.stop()
                raise exceptions.CodeException('Received reply for call id {0} but expected {1}'.format(message.call_id, call_id))
        
            if message.function_id != function_id:
                self.stop()
                raise exceptions.CodeException('Received reply for function id {0} but expected {1}'.format(message.function_id, function_id))
            
            if message.error:
                raise exceptions.CodeException("Error in (asynchronous) communication with worker: " + message.strings[0])
                
            return message.to_result(handle_as_array)
    
        request.add_result_handler(handle_result)
        
        return request
        
    def is_active(self):
        return self.intercomm is not None
        
    def is_inuse(self):
        return self._is_inuse
    
    def is_polling_supported(self):
        return True
        


class MultiprocessingMPIChannel(AbstractMessageChannel):
    """
    Message channel based on JSON messages. 
    
    The remote party functions as a message forwarder.
    Each message is forwarded to a real application using MPI.
    This is message channel is a lot slower than the MPI message
    channel. But, it is useful during testing with
    the MPICH2 nemesis channel. As the tests will run as one
    application on one node they will cause oversaturation 
    of the processor(s) on the node. Each legacy code
    will call the MPI_FINALIZE call and this call will wait
    for the MPI_FINALIZE call of the main test process. During
    this wait it will consume about 10% of the processor power.
    To mitigate this problem, we can use objects of this class
    instead of the normal MPIChannel. Then, part of the
    test is performed in a separate application (at least
    as MPI sees it) and this part can be stopped after each
    sub-test, thus removing unneeded applications. 
    """
    def __init__(self, name_of_the_worker, legacy_interface_type = None,  interpreter_executable = None, **options):
        AbstractMessageChannel.__init__(self, **options)
        
        self.name_of_the_worker = name_of_the_worker
        self.interpreter_executable = interpreter_executable
        
        if not legacy_interface_type is None:
            self.full_name_of_the_worker = self.get_full_name_of_the_worker(legacy_interface_type)
        else:
            self.full_name_of_the_worker = self.name_of_the_worker
            
        self.process = None
    
    @option(type="boolean")
    def debug_with_gdb(self):
        return False
        
    @option
    def hostname(self):
        return None
    
    def start(self):
        name_of_dir = "/tmp/amuse_"+os.getenv('USER')
        self.name_of_the_socket, self.server_socket = self._createAServerUNIXSocket(name_of_dir)
        environment = os.environ.copy()
        
        if 'PYTHONPATH' in environment:
            environment['PYTHONPATH'] = environment['PYTHONPATH'] + ':' +  self._extra_path_item(__file__)
        else:
            environment['PYTHONPATH'] =  self._extra_path_item(__file__)
         
         
        all_options = {}
        for x in self.iter_options():
            all_options[x.name] = getattr(self, x.name)
        
          
        template = """from {3} import {4}
o = {1!r}
m = channel.MultiprocessingMPIChannel('{0}',**o)
m.run_mpi_channel('{2}')"""
        modulename = type(self).__module__
        packagagename, thismodulename = modulename.rsplit('.', 1)
        
        code_string = template.format(
            self.full_name_of_the_worker, 
            all_options, 
            self.name_of_the_socket,
            packagagename,
            thismodulename,
        )
        self.process =  Popen([sys.executable, "-c", code_string], env = environment)
        self.client_socket, undef = self.server_socket.accept()
    
    def is_active(self):
        return self.process is not None
        
    def stop(self):
        self._send(self.client_socket,('stop',(),))
        result = self._recv(self.client_socket)    
        self.process.wait()
        self.client_socket.close()
        self.server_socket.close()
        self._remove_socket(self.name_of_the_socket)
        self.process = None
        
    def run_mpi_channel(self, name_of_the_socket):
        channel = MpiChannel(self.full_name_of_the_worker, **self._local_options)
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
    
    def send_message(self, call_id = 0, function_id = -1, dtype_to_arguments = {}):
        self._send(self.client_socket, ('send_message',(call_id, function_id, dtype_to_arguments),))
        result = self._recv(self.client_socket)
        return result

    def recv_message(self, call_id = 0, function_id = -1, handle_as_array = False):
        self._send(self.client_socket, ('recv_message',(call_id, function_id, handle_as_array),))
        result = self._recv(self.client_socket)    
        return result
    
    def _send(self, client_socket, message):
        message_string = pickle.dumps(message)
        header = struct.pack("i", len(message_string))
        client_socket.sendall(header)
        client_socket.sendall(message_string)
        
    def _recv(self, client_socket):
        header = self._receive_all(client_socket, 4)
        length = struct.unpack("i", header)
        message_string = self._receive_all(client_socket, length[0])
        return pickle.loads(message_string)
        
    def _receive_all(self, client_socket, number_of_bytes):
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
        #client_socket.settimeout(0)header 
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
        
            


    @option(choices=AbstractMessageChannel.DEBUGGERS.keys(), sections=("channel",))
    def debugger(self):
        """Name of the debugger to use when starting the code"""
        return "none"
    
    

    @option(type="boolean")
    def check_mpi(self):
        return True

class SocketMessage(AbstractMessage):
      
    def _receive_all(self, nbytes, thesocket):

        #logging.getLogger("channel").debug("receiving %d bytes", nbytes)
        
        result = []
        
        while nbytes > 0:
            chunk = min(nbytes, 10240)
            data_bytes = thesocket.recv(chunk, socket.MSG_WAITALL)
            
            if len(data_bytes) == 0:
                raise exceptions.CodeException("lost connection to code")
            
            result.append(data_bytes)
            nbytes -= len(data_bytes)
            #logging.getLogger("channel").debug("got %d bytes, result length = %d", len(data_bytes), len(result))
            
        return "".join(result)
     
    def receive(self, socket):
        
        #logging.getLogger("channel").debug("receiving message")
        
        header_bytes = self._receive_all(40, socket)

        flags = numpy.frombuffer(header_bytes, dtype="b", count=4, offset=0)
        
        if flags[0] != self.big_endian:
            raise exceptions.CodeException("endianness in message does not match native endianness")
        
        if flags[1]:
            self.error = True
        else:
            self.error = False
        
        header = numpy.copy(numpy.frombuffer(header_bytes, dtype="i", offset=0))
        
        #logging.getLogger("channel").debug("receiving message with flags %s and header %s", flags, header)

        #id of this call
        self.call_id = header[1]
        
        #function ID
        self.function_id = header[2]
        
        #number of calls in this message
        self.call_count = header[3]
        
        #number of X's in TOTAL
        number_of_ints = header[4]
        number_of_longs = header[5]
        number_of_floats = header[6]
        number_of_doubles = header[7]
        number_of_booleans = header[8]
        number_of_strings = header[9]

        self.ints = self.receive_ints(socket, number_of_ints)
        self.longs = self.receive_longs(socket, number_of_longs)
        self.floats = self.receive_floats(socket, number_of_floats)
        self.doubles = self.receive_doubles(socket, number_of_doubles)
        self.booleans = self.receive_booleans(socket, number_of_booleans)
        self.strings = self.receive_strings(socket, number_of_strings)
        
        #logging.getLogger("channel").debug("message received")
        
    def receive_ints(self, socket, count):
        if count > 0:
            nbytes = count * 4 # size of int
            
            data_bytes = self._receive_all(nbytes, socket)
            
            result = numpy.copy(numpy.frombuffer(data_bytes, dtype='int32'))
            
            return result
        else:
            return []        
            
    def receive_longs(self, socket, count):
        if count > 0:
            nbytes = count * 8 # size of long
            
            data_bytes = self._receive_all(nbytes, socket)
            
            result = numpy.copy(numpy.frombuffer(data_bytes, dtype='int64'))
            
            return result
        else:
            return []
 
        
    def receive_floats(self, socket, count):
        if count > 0:
            nbytes = count * 4 # size of float
            
            data_bytes = self._receive_all(nbytes, socket)
            
            result = numpy.copy(numpy.frombuffer(data_bytes, dtype='f4'))
            
            return result
        else:
            return []
    
          
    def receive_doubles(self, socket, count):
        if count > 0:
            nbytes = count * 8 # size of double
            
            data_bytes = self._receive_all(nbytes, socket)
            
            result = numpy.copy(numpy.frombuffer(data_bytes, dtype='f8'))
            
            return result
        else:
            return []
        

    def receive_booleans(self, socket, count):
        if count > 0:
            nbytes = count * 1 # size of boolean/byte
            
            data_bytes = self._receive_all(nbytes, socket)
            
            result = numpy.copy(numpy.frombuffer(data_bytes, dtype='b'))
            
            return result
        else:
            return []
    
            
    def receive_strings(self, socket, count):
        if count > 0:
            lengths = self.receive_ints(socket, count)
            
            strings = []
            
            for i in range(count):
                data_bytes = self._receive_all(lengths[i], socket)
                strings.append(str(data_bytes.decode('utf-8')))
            
            return strings
        else:
            return []
            
    def nonblocking_receive(self, socket):
        return ASyncSocketRequest(self, socket)
    
    
    def send(self, socket):
        
        flags = numpy.array([self.big_endian, False, False, False], dtype="b")
        
        header = numpy.array([
            self.call_id,
            self.function_id,
            self.call_count,
            len(self.ints),
            len(self.longs),
            len(self.floats),
            len(self.doubles),
            len(self.booleans),
            len(self.strings),
        ], dtype='i')
        
        #logging.getLogger("channel").debug("sending message with flags %s and header %s", flags, header)
        
        socket.sendall(flags.tostring())
        
        socket.sendall(header.tostring())

        self.send_ints(socket, self.ints)
        self.send_longs(socket, self.longs)
        self.send_floats(socket, self.floats)
        self.send_doubles(socket, self.doubles)
        self.send_booleans(socket, self.booleans)
        self.send_strings(socket, self.strings)
        
        #logging.getLogger("channel").debug("message send")
    
    def send_doubles(self, socket, array):
        if len(array) > 0:
            data_buffer = numpy.array(array, dtype='f8')
            socket.sendall(data_buffer.tostring())
            
    def send_ints(self, socket, array):
        if len(array) > 0:
            data_buffer = numpy.array(array, dtype='int32')
            socket.sendall(data_buffer.tostring())
            
    def send_floats(self, socket, array):
        if len(array) > 0:
            data_buffer = numpy.array(array, dtype='f4')
            socket.sendall(data_buffer.tostring())
            
    def send_strings(self, socket, array):
        header = []
        data_bytes = []
        
        for i in range(len(array)):
            
            logging.getLogger("channel").info("sending string %s", array[i])
            
            utf8_string = array[i].encode('utf-8')
            header.append(len(utf8_string))
            data_bytes.append(utf8_string)
  
        self.send_ints(socket, header);
        
        for i in range(len(data_bytes)):
            socket.sendall(data_bytes[i])
        
    def send_booleans(self, socket, array):
        if len(array) > 0:
            data_buffer = numpy.array(array, dtype='b')
            socket.sendall(data_buffer.tostring())

    def send_longs(self, socket, array):
        if len(array) > 0:
            data_buffer = numpy.array(array, dtype='int64')
            socket.sendall(data_buffer.tostring())
        
class SocketChannel(AbstractMessageChannel):
    
    def __init__(self, name_of_the_worker, legacy_interface_type=None, interpreter_executable = None,**options):
        AbstractMessageChannel.__init__(self, **options)
        
        #logging.basicConfig(level=logging.DEBUG)
        
        logging.getLogger("channel").debug("initializing SocketChannel with options %s", options)
       
        #self.name_of_the_worker = name_of_the_worker + "_sockets"
        self.name_of_the_worker = name_of_the_worker

        self.interpreter_executable = interpreter_executable
        
        if self.hostname != None and self.hostname != 'localhost':
            raise exceptions.CodeException("can only run codes on local machine using SocketChannel, not on %s", self.hostname)
            
#        if self.number_of_workers != 0 and self.number_of_workers != 1:
#            raise exceptions.CodeException("can only a single worker for each code using Socket Channel, not " + str(self.number_of_workers))
            
        self.id = 0
        
        if not legacy_interface_type is None:
            self.full_name_of_the_worker = self.get_full_name_of_the_worker(legacy_interface_type)
        else:
            self.full_name_of_the_worker = self.name_of_the_worker
            
        #logging.getLogger("channel").debug("full name of worker is %s", self.full_name_of_the_worker)
        
        self._is_inuse = False
        self.socket = None
    

    @late
    def debugger_method(self):
        return self.DEBUGGERS[self.debugger]
        
    def start(self):
        server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        server_socket.bind(('', 0))
        server_socket.listen(1)
        
        #logging.getLogger("channel").debug("starting socket worker process")
    
        
        #set arguments to name of the worker, and port number we listen on 
    
    
        self.stdout = None
        self.stderr = None
        
        arguments = []
        
        if not self.debugger_method is None:
            command, arguments = self.debugger_method(self.full_name_of_the_worker, self, interpreter_executable = self.interpreter_executable)
        else:
            if self.redirect_stdout_file is None  or self.redirect_stdout_file == "none":
                self.stdout = None
            else:
                self.stdout = open(self.redirect_stdout_file, "w")
    
            if self.redirect_stderr_file is None or self.redirect_stderr_file == "none":
                self.stderr = None
            elif self.redirect_stderr_file == self.redirect_stdout_file:
                #stderr same file as stdout, do not open file twice
                self.stderr = self.stdout
            else:
                self.stderr = open(self.redirect_stderr_file, "w")
            
            if not self.interpreter_executable is None:
                command = self.interpreter_executable
                arguments = [self.full_name_of_the_worker]
            else:
                command = self.full_name_of_the_worker
            
        arguments.insert(0, command)        
        arguments.append(str(server_socket.getsockname()[1]))
        
        if self.number_of_workers > 1:
            logging.getLogger("channel").warn("multiple workers instances for socket worker not properly tested yet")
            command = "mpiexec"
            #prepend with mpiexec and arguments back to front
            arguments.insert(0, str(self.number_of_workers))
            arguments.insert(0, "-np")
            arguments.insert(0, "/usr/bin/mpiexec")

        logging.getLogger("channel").warn("starting process with arguments %s", arguments)

        
        self.process = Popen(arguments, executable = command, stdout = self.stdout, stderr = self.stderr)
        logging.getLogger("channel").warn("waiting for connection from worker")
     
        self.socket, address = server_socket.accept()
        
        self.socket.setblocking(1)
        
        self.socket.setsockopt(socket.SOL_TCP, socket.TCP_NODELAY, 1)
        
        server_socket.close()
        
        #logging.getLogger("channel").debug("got connection from %s", address)
        
        #logging.getLogger("channel").info("worker %s initialized", self.name_of_the_worker)
        
    @option(choices=AbstractMessageChannel.DEBUGGERS.keys(), sections=("channel",))
    def debugger(self):
        """Name of the debugger to use when starting the code"""
        return "none"
        
    @option(sections=("channel",))
    def hostname(self):
        return None
       
    def stop(self):
        if (self.socket == None):
            return
        
        logging.getLogger("channel").info("stopping socket worker %s", self.name_of_the_worker)
        self.socket.close()
        
        self.socket = None
        
        # should lookinto using poll with a timeout or some other mechanism
        # when debugger method is on, no killing
        count = 0
        while(count < 5):
            returncode = self.process.poll()
            if not returncode is None:
                break
            time.sleep(0.2)
            count += 1  
                 
        if not self.stdout is None:
            self.stdout.close()
            
        if not self.stderr is None:
            self.stderr.close()

    def is_active(self):
        return self.socket is not None
    
    def is_inuse(self):
        return self._is_inuse
    
    def determine_length_from_data(self, dtype_to_arguments):
        def get_length(x):
            if x:
                try:
                    if not isinstance(x[0], str):
                        return len(x[0])
                except:
                    return 1
               
               
        
        lengths = map(get_length, dtype_to_arguments.values())
        if len(lengths) == 0:
            return 1
            
        return max(1, max(lengths))
    
    def send_message(self, call_id, function_id, dtype_to_arguments = {}):
        
        call_count = self.determine_length_from_data(dtype_to_arguments)
        
        #logging.getLogger("channel").info("sending message for call id %d, function %d, length %d", id, tag, length)
        
        if self.is_inuse():
            raise exceptions.CodeException("You've tried to send a message to a code that is already handling a message, this is not correct")
        if self.socket is None:
            raise exceptions.CodeException("You've tried to send a message to a code that is not running")
        
        message = SocketMessage(call_id, function_id, call_count, dtype_to_arguments)
        message.send(self.socket)
        
        self._is_inuse = True
        
    def recv_message(self, call_id, function_id, handle_as_array):
           
        self._is_inuse = False
        
        message = SocketMessage()
        
        message.receive(self.socket)

        if message.call_id != call_id:
            self.stop()
            raise exceptions.CodeException('Received reply for call id {0} but expected {1}'.format(message.call_id, call_id))
        if message.function_id != function_id:
            self.stop()
            raise exceptions.CodeException('Received reply for function id {0} but expected {1}'.format(message.function_id, function_id))
        
        if message.error:
            logging.getLogger("channel").info("error message!")
            raise exceptions.CodeException("Error in code: " + message.strings[0])

        return message.to_result(handle_as_array)
        
    def nonblocking_recv_message(self, call_id, function_id, handle_as_array):
        request = SocketMessage().nonblocking_receive(self.socket)
    
        def handle_result(function):
            self._is_inuse = False
    
            message = function()
        
            if message.call_id != call_id:
                self.stop()
                raise exceptions.CodeException('Received reply for call id {0} but expected {1}'.format(message.call_id, call_id))
    
            if message.function_id != function_id:
                self.stop()
                raise exceptions.CodeException('Received reply for function id {0} but expected {1}'.format(message.function_id, function_id))
        
            if message.error:
                raise exceptions.CodeException("Error in (asynchronous) communication with worker: " + message.strings[0])
        
            return message.to_result(handle_as_array)

        request.add_result_handler(handle_result)
    
        return request

class OutputHandler(threading.Thread):
    
    def __init__(self, stream, port):
        threading.Thread.__init__(self)
        self.stream = stream

        logging.getLogger("channel").debug("output handler connecting to daemon")
        
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        
        address = ('localhost', 61575)
        
        try:
            self.socket.connect(address)
        except:
            raise exceptions.CodeException("Could not connect to Ibis Daemon at " + str(address))
        
        self.socket.setsockopt(socket.SOL_TCP, socket.TCP_NODELAY, 1)
        
        self.socket.sendall('TYPE_OUTPUT'.encode('utf-8'))

        #fetch ID of this connection
        
        result = SocketMessage()
        result.receive(self.socket)
        
        self.id = result.strings[0]
        
        self.daemon = True
        self.start()
        
    def run(self):
        
        while True:
            logging.getLogger("channel").debug("receiving data for output")
            data = self.socket.recv(1024)
            
            if len(data) == 0:
                logging.getLogger("channel").debug("end of output", len(data))
                return
            
            logging.getLogger("channel").debug("got %d bytes", len(data))
            
            self.stream.write(data)

class IbisChannel(AbstractMessageChannel):
    
    stdoutHandler = None
    
    stderrHandler = None
    
    socket = None
    
    @classmethod
    def getStdoutID(cls, port):
        if IbisChannel.stdoutHandler is None:
            IbisChannel.stdoutHandler = OutputHandler(sys.stdout, port)
            
        return IbisChannel.stdoutHandler.id
    
    @classmethod
    def getStderrID(cls, port):
        if IbisChannel.stderrHandler is None:
            IbisChannel.stderrHandler = OutputHandler(sys.stderr, port)
            
        return IbisChannel.stderrHandler.id
    
    def __init__(self, name_of_the_worker, legacy_interface_type=None, interpreter_executable = None, **options):
        AbstractMessageChannel.__init__(self, **options)
        
        logging.basicConfig(level=logging.WARN)
        #logging.getLogger("channel").setLevel(logging.DEBUG)
        
        #logging.getLogger("channel").debug("initializing IbisChannel with options %s", options)
       
        self.name_of_the_worker = name_of_the_worker
        self.interpreter_executable = interpreter_executable
        
        if self.hostname == None or self.hostname == 'local':
            self.hostname = 'localhost'
            
        if self.number_of_workers == 0:
            self.number_of_workers = 1
            
        if self.node_label == None:
            self.node_label = "default"
            
        logging.getLogger("channel").debug("number of workers is %d, number of nodes is %s", self.number_of_workers, self.number_of_nodes)
        
        self.daemon_host = 'localhost'    # Ibis deamon always running on the local machine
        self.daemon_port = self.port          # A random-but-fixed port number for the Ibis daemon

        logging.getLogger("channel").debug("port is %d", self.daemon_port)
        
        self.id = 0
        
        if not legacy_interface_type is None:
            #worker specified by type. Figure out where this file is
            #mostly (only?) used by dynamic python codes
            directory_of_this_module = os.path.dirname(inspect.getfile(legacy_interface_type))
            worker_path = os.path.join(directory_of_this_module, self.name_of_the_worker)
            self.full_name_of_the_worker = os.path.normpath(os.path.abspath(worker_path))
           
            self.name_of_the_worker = os.path.basename(self.full_name_of_the_worker)
            
        else:
            #worker specified by executable (usually already absolute)
            self.full_name_of_the_worker = os.path.normpath(os.path.abspath(self.name_of_the_worker))
        
        global_options = GlobalOptions()
        
        relative_worker_path = os.path.relpath(self.full_name_of_the_worker, global_options.amuse_rootdirectory)
            
        self.worker_dir = os.path.dirname(relative_worker_path)
            
        logging.getLogger("channel").debug("name of the worker is %s", self.name_of_the_worker)
        logging.getLogger("channel").debug("worker dir is %s", self.worker_dir)
            
        self._is_inuse = False
      

    def check_if_worker_is_up_to_date(self, object):
        if self.hostname != 'localhost':
            return
        
        logging.getLogger("channel").debug("hostname = %s, checking for worker", self.hostname)
        
        AbstractMessageChannel.check_if_worker_is_up_to_date(self, object)
   
    def start(self):
        logging.getLogger("channel").debug("connecting to daemon")
        
        #if redirect = none, set output file to console stdout stream ID, otherwise make absolute
        if (self.redirect_stdout_file == 'none'):
            self.redirect_stdout_file = IbisChannel.getStdoutID(self.port)
        else:
            self.redirect_stdout_file = os.path.abspath(self.redirect_stdout_file)

        #if redirect = none, set error file to console stderr stream ID, otherwise make absolute
        if (self.redirect_stderr_file == 'none'):
            self.redirect_stderr_file = IbisChannel.getStderrID(self.port)
        else:
            self.redirect_stderr_file = os.path.abspath(self.redirect_stderr_file)
        
        logging.getLogger("channel").debug("output send to = " + self.redirect_stdout_file)
        
        logging.getLogger("channel").debug("error send to = " + self.redirect_stderr_file)
        
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            self.socket.connect((self.daemon_host, self.daemon_port))
        except:
            self.socket = None
            raise exceptions.CodeException("Could not connect to Ibis Daemon at " + str(self.daemon_port))
        
        self.socket.setblocking(1)
        
        self.socket.setsockopt(socket.SOL_TCP, socket.TCP_NODELAY, 1)
        
        self.socket.sendall('TYPE_WORKER'.encode('utf-8'))
        
        arguments = {'string': [self.name_of_the_worker, self.worker_dir, self.hostname, self.redirect_stdout_file, self.redirect_stderr_file, self.node_label], 'int32': [self.number_of_workers, self.number_of_nodes, self.number_of_threads], 'bool': [self.copy_worker_code]}
        
        message = SocketMessage(call_id=1, function_id=10101010, call_count=1, dtype_to_arguments=arguments);

        message.send(self.socket)
        
        logging.getLogger("channel").info("waiting for worker %s to be initialized", self.name_of_the_worker)

        result = SocketMessage()
        result.receive(self.socket)
        
        if result.error:
            logging.getLogger("channel").error("Could not start worker: %s", result.strings[0])
            self.stop()
            raise exceptions.CodeException("Could not start worker for " + self.name_of_the_worker + ": " + result.strings[0])
        
        logging.getLogger("channel").info("worker %s initialized", self.name_of_the_worker)
        
    @option(choices=AbstractMessageChannel.DEBUGGERS.keys(), sections=("channel",))
    def debugger(self):
        """Name of the debugger to use when starting the code"""
        return "none"
        
    @option(sections=("channel",))
    def hostname(self):
        return None
    
    @option(type="int", sections=("channel",))
    def number_of_nodes(self):
        return 1
    
    @option(type="int", sections=("channel",))
    def port(self):
        return 61575
    
    @option(type="int", sections=("channel",))
    def number_of_threads(self):
        return 0
    
    @option(type="string", sections=("channel",))
    def node_label(self):
        return None
    
    @option(type="boolean", sections=("channel",))
    def copy_worker_code(self):
        return False
       
    def stop(self):
        if self.socket is not None:
            logging.getLogger("channel").info("stopping worker %s", self.name_of_the_worker)
            self.socket.close()
            self.socket = None

    def is_active(self):
        return self.socket is not None 
    
    def is_inuse(self):
        return self._is_inuse
    
    def determine_length_from_data(self, dtype_to_arguments):
        def get_length(x):
            if x:
                try:
                    if not isinstance(x[0], str):
                        return len(x[0])
                except:
                    return 1
               
               
        
        lengths = map(get_length, dtype_to_arguments.values())
        if len(lengths) == 0:
            return 1
            
        return max(1, max(lengths))
    
    def send_message(self, call_id, function_id, dtype_to_arguments = {}):
        
        call_count = self.determine_length_from_data(dtype_to_arguments)
        
        logging.getLogger("channel").info("sending message for call id %d, function %d, length %d", id, function_id, call_count)
        
        if self.is_inuse():
            raise exceptions.CodeException("You've tried to send a message to a code that is already handling a message, this is not correct")
        if self.socket is None:
            raise exceptions.CodeException("You've tried to send a message to a code that is not running")
        
        message = SocketMessage(call_id, function_id, call_count, dtype_to_arguments, False, False)
        message.send(self.socket)
        
        self._is_inuse = True
        
    def recv_message(self, call_id, function_id, handle_as_array):
           
        self._is_inuse = False
        
        message = SocketMessage()
        
        message.receive(self.socket)
        
        if message.error:
            raise exceptions.CodeException("Error in worker: " + message.strings[0])
        
        return message.to_result(handle_as_array)
    
    def nonblocking_recv_message(self, tag, handle_as_array):
        #       raise exceptions.CodeException("Nonblocking receive not supported by IbisChannel")
        request = SocketMessage().nonblocking_receive(self.socket)
        
        def handle_result(function):
            self._is_inuse = False
        
            message = function()
            
            if message.error:
                raise exceptions.CodeException("Error in worker: " + message.strings[0])
                
            return message.to_result(handle_as_array)
    
        request.add_result_handler(handle_result)
        
        return request
    
