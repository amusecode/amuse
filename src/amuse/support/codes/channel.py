import inspect
import numpy
import os
import os.path
import cPickle as pickle


import sys
import struct
import threading
import select
import tempfile

import socket
import array

import logging

from mpi4py import MPI
from subprocess import Popen, PIPE

from amuse.support.options import OptionalAttributes, option
from amuse.support.core import late
from amuse.support import exceptions

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
                
        self.message.recieve_content(self.comm, self.header)
        
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
        
class AbstractMessage(object):
    
    def __init__(self, tag = -1, length= 1, dtype_to_arguments = {}):
        self.tag = tag
        self.length = length
        self.ints = []
        self.doubles = []
        self.floats = []
        self.strings = []
        self.booleans = []
        self.longs = []
        
        self.pack_data(dtype_to_arguments)
        
    def pack_data(self,dtype_to_arguments):
        for dtype, attrname in self.dtype_to_message_attribute():
            if dtype in dtype_to_arguments:
                array = pack_array( dtype_to_arguments[dtype], self.length, dtype)
                setattr(self, attrname, array)
    
    def to_result(self, handle_as_array = False):
        dtype_to_result = {}
        for dtype, attrname in self.dtype_to_message_attribute():
            result = getattr(self, attrname)
            if self.length > 1 or handle_as_array:
                dtype_to_result[dtype] = unpack_array(result , self.length, dtype)
            else:
                dtype_to_result[dtype] = result
                    
        return dtype_to_result
    
    def dtype_to_message_attribute(self):
        return (
            ('float64', 'doubles'),
            ('float32', 'floats'),
            ('int32', 'ints'),
            ('bool', 'booleans'),
            ('string', 'strings'),
            ('int64', 'longs'),
        )
    
    def recieve(self, comm):
        raise NotImplementedError
        
    def send(self, comm):
        raise NotImplementedError
        
    
class MPIMessage(AbstractMessage):
    
        
    def recieve(self, comm):
        header = numpy.zeros(8,  dtype='i')
        
        self.mpi_recieve(comm, [header, MPI.INT])
    
        self.recieve_content(comm, header)
        
    def recieve_content(self, comm, header):
        self.tag = header[0]
        self.length = header[1]
        
        number_of_doubles = header[2]
        number_of_ints = header[3]
        number_of_floats = header[4]
        number_of_strings = header[5]
        number_of_booleans = header[6]
        number_of_longs = header[7]
        
        self.doubles = self.recieve_doubles(comm, self.length, number_of_doubles)
        self.ints = self.recieve_ints(comm, self.length, number_of_ints)
        self.floats = self.recieve_floats(comm, self.length, number_of_floats)
        self.strings = self.recieve_strings(comm, self.length, number_of_strings)
        self.booleans = self.recieve_booleans(comm, self.length, number_of_booleans)
        self.longs = self.recieve_longs(comm, self.length, number_of_longs)
        
    def nonblocking_recieve(self, comm):
        header = numpy.zeros(8,  dtype='i')
        request = self.mpi_nonblocking_recieve(comm, [header, MPI.INT])
        return ASyncRequest(request, self, comm,  header)
    
    def recieve_doubles(self, comm, length, total):
        if total > 0:
            result = numpy.empty(total * length,  dtype='d')
            self.mpi_recieve(comm,[result, MPI.DOUBLE])
            return result
        else:
            return []
            
    def recieve_ints(self, comm, length, total):
        if total > 0:
            result = numpy.empty(total * length,  dtype='i')
            self.mpi_recieve(comm,[result, MPI.INT])
            return result
        else:
            return []
            
    def recieve_longs(self, comm, length, total):
        if total > 0:
            result = numpy.empty(total * length,  dtype='int64')
            self.mpi_recieve(comm,[result, MPI.INTEGER8])
            return result
        else:
            return []
            
    def recieve_floats(self, comm, length, total):
        if total > 0:
            result = numpy.empty(total * length,  dtype='f')
            self.mpi_recieve(comm,[result, MPI.FLOAT])
            return result
        else:
            return []
            
    
    def recieve_booleans(self, comm, length, total):
        if total > 0:
            result = numpy.empty(total * length,  dtype='int32')
            self.mpi_recieve(comm,[result, MPI.LOGICAL])
            return result == 1
        else:
            return []
    
            
    def recieve_strings(self, comm, length, total):
        if total > 0:
            offsets = numpy.empty(length * total, dtype='i')
            self.mpi_recieve(comm,[offsets, MPI.INT])
            
            bytes = numpy.empty((offsets[-1] + 1), dtype=numpy.uint8)
            self.mpi_recieve(comm,[bytes,  MPI.CHARACTER])
            
            strings = []
            begin = 0
            
            for end in offsets:
                bytes_of_string = bytes[begin:end]
                string = ''.join(map(chr, bytes_of_string))
                begin = end + 1
                strings.append(string)
            return strings
        else:
            return []
        
    
    def send(self, comm):
        header = numpy.array([
            self.tag, 
            self.length, 
            len(self.doubles) / self.length, 
            len(self.ints) / self.length, 
            len(self.floats) / self.length, 
            len(self.strings) / self.length,
            len(self.booleans) / self.length,
            len(self.longs) / self.length,
        ], dtype='i')
        
        self.mpi_send(comm, [header, MPI.INT])
        
        
        self.send_doubles(comm, self.doubles)
        self.send_ints(comm, self.ints)
        self.send_floats(comm, self.floats)
        self.send_strings(comm, self.strings)
        self.send_booleans(comm, self.booleans)
        self.send_longs(comm, self.longs)
        
    
    def send_doubles(self, comm, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='d')
            self.mpi_send(comm,[buffer, MPI.DOUBLE])
            
    def send_ints(self, comm, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='int32')
            self.mpi_send(comm,[buffer, MPI.INT])
            
    def send_floats(self, comm, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='f')
            self.mpi_send(comm, [buffer, MPI.FLOAT])
            
    def send_strings(self, comm, array):
        if len(array) == 0:
            return
            
        offsets = self.string_offsets(array)
        self.mpi_send(comm, [offsets, MPI.INT])
        
        bytes = []
        for string in array:
            bytes.extend([ord(ch) for ch in string])
            bytes.append(0)
          
        chars = numpy.array(bytes, dtype=numpy.uint8)
        self.mpi_send(comm, [chars, MPI.CHARACTER])
        
    def send_booleans(self, comm, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='int32')
            self.mpi_send(comm, [buffer, MPI.LOGICAL])
            
        
    def send_longs(self, comm, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='int64')
            print buffer
            self.mpi_send(comm,[buffer, MPI.INTEGER8])    
    
    def string_offsets(self, array):
        offsets = numpy.zeros(len(array), dtype='i')
        offset = 0
        index = 0
        
        for string in array:
            offset += len(string)
            offsets[index] = offset
            offset += 1
            index += 1
        
        return offsets
        
    def mpi_nonblocking_recieve(self, comm, array):
        raise NotImplementedError()
    
    def mpi_recieve(self, comm, array):
        raise NotImplementedError()
        
    def mpi_send(self, comm, array):
        raise NotImplementedError()
    
    
class ServerSideMPIMessage(MPIMessage):
    
    def mpi_recieve(self, comm, array):
        comm.Recv(array,  source=0, tag=999)
        
    def mpi_send(self, comm, array):
        comm.Bcast(array, root=MPI.ROOT)
    
    def mpi_nonblocking_recieve(self, comm, array):
        return comm.Irecv(array,  source=0, tag=999)

    
class ClientSideMPIMessage(MPIMessage):
    
    def mpi_recieve(self, comm, array):
        comm.Bcast(array, root = 0)
        
    def mpi_send(self, comm, array):
        comm.Send(array, dest=0, tag = 999)

    def mpi_nonblocking_recieve(self, comm, array):
        return comm.Irecv(array,  source=0, tag=999)


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

class MessageChannel(OptionalAttributes):
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
    def GDB(cls, full_name_of_the_worker):
        arguments = ['-hold', '-display', os.environ['DISPLAY'], '-e', 'gdb',  full_name_of_the_worker]
        command = 'xterm'
        return command, arguments
        
    @classmethod
    def DDD(cls, full_name_of_the_worker):
        arguments = ['-display', os.environ['DISPLAY'], '-e', 'ddd',  full_name_of_the_worker]
        command = 'xterm'
        return command, arguments
        
    @classmethod
    def VALGRIND(cls, full_name_of_the_worker):
        #arguments = ['-hold', '-display', os.environ['DISPLAY'], '-e', 'valgrind',  full_name_of_the_worker]
        arguments = [full_name_of_the_worker]
        command = 'valgrind'
        return command, arguments
        
        
    @classmethod
    def XTERM(cls, full_name_of_the_worker):
        arguments = ['-hold', '-display', os.environ['DISPLAY'], '-e', full_name_of_the_worker]
        command = 'xterm'
        return command, arguments
        

    @classmethod
    def SAVE(cls, full_name_of_the_worker):
        code = """from mpi4py import rc
rc.initialize = False
from mpi4py import MPI
from subprocess import call
import sys
import time
import signal

def sighandler(s, frame):
    print "SIGNAL",s

signal.signal(signal.SIGINT, sighandler)
returncode = call(sys.argv[1:])
p = MPI.Comm.Get_parent()
p.Disconnect()"""
        fd, name = tempfile.mkstemp()
        with os.fdopen(fd, 'w') as f:
            f.write(code)
        
        arguments = [name , full_name_of_the_worker]
        #arguments = ['-hold', '-display', os.environ['DISPLAY'], '-e', sys.executable, '-c', code, full_name_of_the_worker]
        command = sys.executable
        #command = 'xterm'

        print arguments
        
        return command, arguments
        
    
        
    def get_full_name_of_the_worker(self, type):
        if os.path.isabs(self.name_of_the_worker):
            if os.path.exists(self.name_of_the_worker):
                return self.name_of_the_worker
            
        tried_workers = []
        found = False
        current_type=type
        while not found:
            directory_of_this_module = os.path.dirname(inspect.getfile(current_type))
            full_name_of_the_worker = os.path.join(directory_of_this_module, self.name_of_the_worker)
            full_name_of_the_worker = os.path.normpath(os.path.abspath(full_name_of_the_worker))
            found = os.path.exists(full_name_of_the_worker)
            if not found:
                tried_workers.append(full_name_of_the_worker)
                current_type = current_type.__bases__[0]
                if current_type.__bases__[0] is object:
                    raise exceptions.CodeException("The worker application does not exists, it should be at: {0}".format(tried_workers))
            else:
                found = True
        return full_name_of_the_worker

    def send_message(self, tag, id=0, dtype_to_arguments = {}, length = 1):
        pass
        
    def recv_message(self, tag, handle_as_array = False):
        pass
        
    def start(self):
        pass
        
    def stop(self):
        pass

    def is_active(self):
        return True
 
MessageChannel.DEBUGGERS = {
    "none":None,
    "gdb":MessageChannel.GDB, 
    "ddd":MessageChannel.DDD, 
    "xterm":MessageChannel.XTERM,
    "valgrind":MessageChannel.VALGRIND,
    "save": MessageChannel.SAVE,
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
    name_of_the_vendor, version = MPI.get_vendor()
    if name_of_the_vendor == 'MPICH2':
        process = Popen(['mpdtrace'], stdout = PIPE, stderr = PIPE)
        (output_string, error_string) = process.communicate()
        return not (process.returncode == 255)
    else:
        return True


class MpiChannel(MessageChannel):
    """
    Message channel based on MPI calls to send and recv the messages
    
    :argument name_of_the_worker: Name of the application to start
    :argument number_of_workers: Number of parallel processes
    :argument legacy_interface_type: Type of the legacy interface
    :argument debug_with_gdb: If True opens an xterm with a gdb to debug the remote process
    :argument hostname: Name of the node to run the application on
    """
    _mpi_is_broken_after_possible_code_crash = False
    
    
    def __init__(self, name_of_the_worker, legacy_interface_type = None,  **options):
        MessageChannel.__init__(self, **options)
               
        self.name_of_the_worker = name_of_the_worker
                
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
    
    @option(type="boolean")
    def check_mpi(self):
        return True
        
    @option(type="boolean")
    def debug_with_gdb(self):
        return False
        
    @option(sections=("channel",))
    def hostname(self):
        return None
    
    @option(type="int", sections=("channel",))
    def number_of_workers(self):
        return 1
        
    @option(choices=MessageChannel.DEBUGGERS.keys(), sections=("channel",))
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
    
    def start(self):
        
        must_close_std_streams = True
        if self.debug_with_gdb or (not self.debugger_method is None):
            if not 'DISPLAY' in os.environ:
                arguments = None
                command = self.full_name_of_the_worker
            else:
                must_close_std_streams = False
                if self.debugger is None:
                    command, arguments = self.GDB(self.full_name_of_the_worker)
                else:
                    command, arguments = self.debugger_method(self.full_name_of_the_worker)
        else:
            arguments = None
            command = self.full_name_of_the_worker
        

        fd_stdin = None
        fd_stdout = None
        fd_stderr = None
        if must_close_std_streams:
            fd_stdin = os.dup(0)
            zero = open('/dev/null','r')
            os.dup2(zero.fileno(), 0)
            if not self.redirect_stdout_file == "none":
                fd_stdout = os.dup(1)
                zero = open(self.redirect_stdout_file,'a')
                os.dup2(zero.fileno(), 1)
            if not self.redirect_stderr_file == "none":
                fd_stderr = os.dup(2)
                zero = open(self.redirect_stderr_file,'a')
                os.dup2(zero.fileno(), 2)
        try:
            self.intercomm = MPI.COMM_SELF.Spawn(command, arguments, self.number_of_workers, info = self.info)
        finally: 
            if must_close_std_streams:
                os.dup2(fd_stdin, 0)
                os.close(fd_stdin)
                if not fd_stdout == None:
                    os.dup2(fd_stdout, 1)
                    os.close(fd_stdout)
                if not fd_stderr == None:
                    os.dup2(fd_stderr, 2)
                    os.close(fd_stderr)
            
        
        
    def stop(self):
        if not self.intercomm is None:
            try:
                self.intercomm.Disconnect()
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
        
        
    def send_message(self, tag, id=0, dtype_to_arguments = {}, length = 1):

        if self.is_inuse():
            raise exceptions.CodeException("You've tried to send a message to a code that is already handling a message, this is not correct")
        if self.intercomm is None:
            raise exceptions.CodeException("You've tried to send a message to a code that is not running")
        
        length = self.determine_length_from_data(dtype_to_arguments)
        
        if length > self.max_message_length:
            self.split_message(tag, id, dtype_to_arguments, length)
            return
            
        message = ServerSideMPIMessage(tag, length, dtype_to_arguments)
        message.send(self.intercomm)
        
        self._is_inuse = True
        
    def split_message(self, tag, id, dtype_to_arguments, length):
        
        def split_input_array(i, arr_in):
            if length == 1:
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
        
        for i in range(1+(length-1)/self.max_message_length):
            split_dtype_to_argument = {}
            for key, value in dtype_to_arguments.iteritems():
                split_dtype_to_argument[key] = split_input_array(i, value)
                
            self.send_message(
                tag, 
                id, 
                split_dtype_to_argument, 
                min(length,self.max_message_length)
            )
            
            partial_dtype_to_result = self.recv_message(id, True)
            for datatype, value in partial_dtype_to_result.iteritems():
                if not datatype in dtype_to_result:
                    dtype_to_result[datatype] = [] 
                    for j, element in enumerate(value):
                        if datatype == 'string':
                            dtype_to_result[datatype].append([])
                        else:
                            dtype_to_result[datatype].append(numpy.zeros((length,), dtype=datatype))
                            
                for j, element in enumerate(value):
                    if datatype == 'string':
                        dtype_to_result[datatype][j].extend(element)
                    else:
                        dtype_to_result[datatype][j][i*self.max_message_length:(i+1)*self.max_message_length] = element
                
            #print partial_dtype_to_result
            length -= self.max_message_length
        
        self._communicated_splitted_message = True
        self._merged_results_splitted_message = dtype_to_result
    
    def recv_message(self, tag, handle_as_array):
        
        self._is_inuse = False
        
        if self._communicated_splitted_message:
            x = self._merged_results_splitted_message
            self._communicated_splitted_message = False
            del self._merged_results_splitted_message
            return x
        
        message = ServerSideMPIMessage()
        try:
            message.recieve(self.intercomm)
        except MPI.Exception as ex:
            self.stop()
            raise
        
        if message.tag == -1:
            raise exceptions.CodeException("Not a valid message, message is not understood by legacy code")
        elif message.tag == -2:
            self.stop()
            raise exceptions.CodeException("Fatal error in code, code has exited")
        
        return message.to_result(handle_as_array)
        
    def nonblocking_recv_message(self, tag, handle_as_array):
        request = ServerSideMPIMessage().nonblocking_recieve(self.intercomm)
        
        def handle_result(function):
            self._is_inuse = False
        
            message = function()
            
            if message.tag < 0:
                raise exceptions.CodeException("Not a valid message, message is not understood by legacy code")
                
            return message.to_result(handle_as_array)
    
        request.add_result_handler(handle_result)
        
        return request
        
    def is_active(self):
        return self.intercomm is not None
        
    def is_inuse(self):
        return self._is_inuse
        


class MultiprocessingMPIChannel(MessageChannel):
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
    def __init__(self, name_of_the_worker, legacy_interface_type = None,  **options):
        MessageChannel.__init__(self, **options)
        
        self.name_of_the_worker = name_of_the_worker
        
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
    
    @option(type="int")
    def number_of_workers(self):
        return 1
        
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
        print repr(all_options)
        
          
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
         
    def send_message(self, tag, id=0, dtype_to_arguments = {}, length = 1):
        self._send(self.client_socket, ('send_message',(tag, id, dtype_to_arguments, length),))
        result = self._recv(self.client_socket)
        return result
            
    def recv_message(self, tag, handle_as_array):
        self._send(self.client_socket, ('recv_message',(tag,handle_as_array),))
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
        
            


    @option(choices=MessageChannel.DEBUGGERS.keys(), sections=("channel",))
    def debugger(self):
        """Name of the debugger to use when starting the code"""
        return "none"
    
    

    @option(type="boolean")
    def check_mpi(self):
        return True
   
   
class IbisMessage(AbstractMessage):
  
    def _send_int(self, some_integer, socket):
        'Sends an integer'
        narray = numpy.empty(1, dtype=">i4");
        narray[0] = some_integer
        
        socket.sendall(narray.tostring())
        
    def _receive_int(self, socket):
        bytes = self.socket.recv(4);
        
        while len(bytes) < 4:
            bytes = bytes + socket.recv(4 - len(bytes))
            
        return numpy.frombuffer(bytes, dtype=">i4")[0]
    
    def _send_string(self, some_string, socket):
        'Sends a string through a socket by first sending a length header, then the string'
        
        logging.getLogger("ibis").debug("sending string '%s'",some_string)

        utf8_string = some_string.encode('utf-8')

        self._send_int(len(utf8_string), socket)

        socket.sendall(utf8_string)
#        
    def _receive_string(self, socket):
        'Receives a string from a socket by first receiving a length header, then the string'

        logging.getLogger("ibis").debug("receiving string");
        
        length = self._receive_int()

        bytes = socket.recv(length);
        
        while len(bytes) < length:
            bytes = bytes + socket.recv(length - len(bytes))

        result = bytes.decode('utf-8')
    
        logging.getLogger("ibis").debug("received string '%s'", result)
            
        return result
    

    def _receive_all(self, nbytes, socket):
        bytes = socket.recv(nbytes)

        while len(bytes) < nbytes:
            bytes = bytes + socket.recv(nbytes - len(bytes))

        return bytes
     
    def recieve(self, socket):
        
        header_bytes = self._receive_all(32, socket)
        
        header = numpy.frombuffer(header_bytes, dtype=">i4")
        
        #function ID
        self.tag = header[0]
        
        #number of calls in this message
        self.length = header[1]
        
        #number of X's per call
        number_of_doubles = header[2]
        number_of_ints = header[3]
        number_of_floats = header[4]
        number_of_strings = header[5]
        number_of_booleans = header[6]
        number_of_longs = header[7]
        
        self.doubles = self.recieve_doubles(socket, self.length, number_of_doubles)
        self.ints = self.recieve_ints(socket, self.length, number_of_ints)
        self.floats = self.recieve_floats(socket, self.length, number_of_floats)
        self.strings = self.recieve_strings(socket, self.length, number_of_strings)
        self.booleans = self.recieve_booleans(socket, self.length, number_of_booleans)
        self.longs = self.recieve_longs(socket, self.length, number_of_longs)
        
    def recieve_doubles(self, socket, length, count):
        if count > 0:
            nbytes = count * length * 8 # size of double
            
            bytes = self._receive_all(nbytes, socket)
            
            result = numpy.frombuffer(bytes, dtype='d')
            
            return result
        else:
            return []
        
    def recieve_ints(self, socket, length, count):
        if count > 0:
            nbytes = count * length * 4 # size of int
            
            bytes = self._receive_all(nbytes, socket)
            
            result = numpy.frombuffer(bytes, dtype='i')
            
            return result
        else:
            return []
            
    def recieve_longs(self, socket, length, count):
        if count > 0:
            nbytes = count * length * 8 # size of long
            
            bytes = self._receive_all(nbytes, socket)
            
            result = numpy.frombuffer(bytes, dtype='i')
            
            return result
        else:
            return []
            
    def recieve_floats(self, socket, length, count):
        if count > 0:
            nbytes = count * length * 8 # size of float
            
            bytes = self._receive_all(nbytes, socket)
            
            result = numpy.frombuffer(bytes, dtype='f')
            
            return result
        else:
            return []
    
    def recieve_booleans(self, socket, length, count):
        if count > 0:
            nbytes = count * length * 4 # size of boolean/int
            
            bytes = self._receive_all(nbytes, socket)
            
            result = numpy.frombuffer(bytes, dtype='int32')
            
            return result
        else:
            return []
    
            
    def recieve_strings(self, socket, length, count):
        if count > 0:
            strings = []
            
            for i in range(count * length):
                strings.append(self._receive_string(socket))
            
            return strings
        else:
            return []
    
    def send(self, socket):
        header = numpy.array([
            self.tag, 
            self.length, 
            len(self.doubles) / self.length, 
            len(self.ints) / self.length, 
            len(self.floats) / self.length, 
            len(self.strings) / self.length,
            len(self.booleans) / self.length,
            len(self.longs) / self.length,
        ], dtype='i')
        
        socket.sendall(header.tostring())
        
        self.send_doubles(socket, self.doubles)
        self.send_ints(socket, self.ints)
        self.send_floats(socket, self.floats)
        self.send_strings(socket, self.strings)
        self.send_booleans(socket, self.booleans)
        self.send_longs(socket, self.longs)
        
    
    def send_doubles(self, socket, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='d')
            socket.sendall(buffer.tostring())
            
    def send_ints(self, socket, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='int32')
            socket.sendall(buffer.tostring())
            
    def send_floats(self, socket, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='f')
            socket.sendall(buffer.tostring())
            
    def send_strings(self, socket, array):
        for i in range(len(array)):
            self._send_string(array[i], socket)
        
    def send_booleans(self, socket, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='int32')
            socket.sendall(buffer.tostring())
        
    def send_longs(self, socket, array):
        if len(array) > 0:
            buffer = numpy.array(array,  dtype='int64')
            socket.sendall(buffer.tostring())
        
class IbisChannel(MessageChannel):
    
    def __init__(self, name_of_the_worker, legacy_interface_type = None, hostname="localhost", number_of_workers=1, **options):
        MessageChannel.__init__(self, **options)
       
        self.name_of_the_worker = name_of_the_worker

        if hostname == None:
            hostname = 'localhost'
        
        self.hostname = hostname
        self.number_of_workers = number_of_workers
        
        self.daemon_host = 'localhost'    # Ibis deamon always running on the local machine
        self.daemon_port = 61575          # A random-but-fixed port number for the Ibis daemon
        
        logging.basicConfig(level=logging.DEBUG)
        
        if not legacy_interface_type is None:
            self.full_name_of_the_worker = self.get_full_name_of_the_worker( legacy_interface_type)
        else:
            self.full_name_of_the_worker = self.name_of_the_worker
            
        self._is_inuse = False
        self.socket = None
    

    def start(self):
        logging.getLogger("ibis").debug("connecting to daemon")
        
        self.socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.socket.connect((self.daemon_host, self.daemon_port))

        #send options, start  with magic string (daemon will expect this to be
        #send first
#        self._send_string('magic string')
#        self._send_string(self.name_of_the_worker)
#        self._send_string(self.hostname)
#        self._send_int(self.number_of_workers)
        
        arguments = {'strings': ['magic_string', self.name_of_the_worker, self.hostname], 'ints': [self.number_of_workers]}
        
        message = IbisMessage(10101010, 1, arguments)

        message.send(self.socket)
        
        logging.getLogger("ibis").info("waiting for worker %s to be initialized", self.name_of_the_worker)

        result = IbisMessage()
        result.recieve(self.socket)
        
        print result.to_result()

        if result != "OK":
            raise Exception('cannot start Ibis worker')
        
        logging.getLogger("ibis").info("worker %s initialized", self.name_of_the_worker)
        
    @option(choices=MessageChannel.DEBUGGERS.keys(), sections=("channel",))
    def debugger(self):
        """Name of the debugger to use when starting the code"""
        return "none"

        
    def stop(self):
        logging.getLogger("ibis").info("stopping worker %s", self.name_of_the_worker)
        self.socket.close()

    def is_active(self):
        return True    
    
    def is_inuse(self):
        return self._is_inuse
    
    def send_message(self, tag, id=0, dtype_to_arguments = {}, length = 1):
        
        print "tag" 
        print tag
        print "id" 
        print id
        print "type" 
        print dtype_to_arguments
        print "length"
        print length
        
        if self.is_inuse():
            raise exceptions.CodeException("You've tried to send a message to a code that is already handling a message, this is not correct")
        if self.socket is None:
            raise exceptions.CodeException("You've tried to send a message to a code that is not running")
        
        message = IbisMessage(tag, length, dtype_to_arguments)
        message.send(self.socket)
        
        self._is_inuse = True
        
        raise NotImplementedError
        
    def recv_message(self, tag, handle_as_array = False):
           
        self._is_inuse = False
        
        message = IbisMessage()
        
        message.recieve(self.socket)
        
        if message.tag == -1:
            raise exceptions.CodeException("Not a valid message, message is not understood by legacy code")
        elif message.tag == -2:
            self.stop()
            raise exceptions.CodeException("Fatal error in code, code has exited")
        
        return message.to_result(handle_as_array)
