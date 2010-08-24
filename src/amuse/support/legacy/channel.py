import inspect
import numpy
import os
import os.path
import cPickle as pickle


import sys
import struct
import threading
import select

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
        
        self.mpi_recieve(comm, [header, MPI.INT])
    
        self.recieve_content(comm, header)
        
    def recieve_content(self, comm, header):
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
        
    def nonblocking_recieve(self, comm):
        header = numpy.zeros(6,  dtype='i')
        
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
            
    def recieve_floats(self, comm, length, total):
        if total > 0:
            result = numpy.empty(total * length,  dtype='f')
            self.mpi_recieve(comm,[result, MPI.FLOAT])
            return result
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
            len(self.strings) / self.length
        ], dtype='i')
        
        self.mpi_send(comm, [header, MPI.INT])
        
        
        self.send_doubles(comm, self.doubles)
        self.send_ints(comm, self.ints)
        self.send_floats(comm, self.floats)
        self.send_strings(comm, self.strings)
        
    
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
            
            
    
    def mpi_recieve(self, comm, array):
        comm.Bcast(array, root = 0)
        
        
    def mpi_send(self, comm, array):
        comm.Send(array, dest=0, tag = 999)


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
    
    
class ServerSideMessage(Message):
    
    
    def mpi_recieve(self, comm, array):
        comm.Recv(array,  source=0, tag=999)
        
    def mpi_send(self, comm, array):
        comm.Bcast(array, root=MPI.ROOT)
    
    def mpi_nonblocking_recieve(self, comm, array):
        return comm.Irecv(array,  source=0, tag=999)


MAPPING = {}

def pack_array(array, length,  dtype):
    if dtype == 'string':
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
        
    
        
    def get_full_name_of_the_worker(self, type):
        if os.path.isabs(self.name_of_the_worker):
            if os.path.exists(self.name_of_the_worker):
                return self.name_of_the_worker
            
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
                if current_type.__bases__[0] is object:
                    raise exceptions.LegacyException("The worker application does not exists, it should be at: {0}".format(tried_workers))
            else:
                found = True
        return full_name_of_the_worker

    def send_message(self, tag, id=0, doubles_in=[], ints_in=[], floats_in=[], chars_in=[], length = 1):
        pass
        
    def recv_message(self, tag, handle_as_array = False):
        pass

 
MessageChannel.DEBUGGERS = {
    "none":None,
    "gdb":MessageChannel.GDB, 
    "ddd":MessageChannel.DDD, 
    "xterm":MessageChannel.XTERM,
    "valgrind":MessageChannel.VALGRIND,
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

class RunStandardStreamsRedirection(object):
    _INSTANCE = None
    
    @classmethod
    def instance(cls):
        if cls._INSTANCE is None:
            cls._INSTANCE = cls()
            cls._INSTANCE.start()
        return cls._INSTANCE
    
    @classmethod
    def clear(cls):
        if not cls._INSTANCE is None:
            cls._INSTANCE.stop()
            cls._INSTANCE = None
    
    def __init__(self):
        self.copied_fds = False
        self.restored_fds = False
        self.must_run = False
        
    def setup(self):
        self.fd_stdin_r, self.fd_stdin_w = os.pipe()
        self.fd_stdout_r , self.fd_stdout_w= os.pipe()
        self.fd_stderr_r, self.fd_stderr_w = os.pipe()
    
    def before_spawn(self):
        if not self.copied_fds:
            self.fd_stdin = os.dup(0)
            self.fd_stdout = os.dup(1)
            self.fd_stderr = os.dup(2)
            
            os.dup2(self.fd_stdin_r, 0)
            os.dup2(self.fd_stdout_w, 1)
            os.dup2(self.fd_stderr_w, 2)
            
            self.copied_fds = True
    
    def after_spawn(self):
        if not self.copied_fds:
            return
        if not self.restored_fds:
            os.dup2(self.fd_stdin, 0)
            os.dup2(self.fd_stdout, 1)
            os.dup2(self.fd_stderr, 2)
            
            os.close(self.fd_stdin)
            os.close(self.fd_stdout)
            os.close(self.fd_stderr)
        
            self.restored_fds = True
    
    def start(self):
        self.must_run = True
        self.setup()
        self.process = threading.Thread(target = self.runloop)
        self.process.daemon = True
        self.process.start()
    
    def stop(self):
        self.must_run = False
        
    def runloop(self):
        readers = (self.fd_stdout_r, self.fd_stderr_r)
        while self.must_run:
            ready_to_read, ignore1, ignore2 = select.select(readers, [], [])
            for x in ready_to_read:
                if x == self.fd_stdout_r:
                    bytes = os.read(self.fd_stdout_r, 1024)
                    sys.stdout.write(bytes)
                elif x == self.fd_stderr_r:
                    bytes = os.read(self.fd_stderr_r, 1024)
                    sys.stderr.write(bytes)
                else:
                    pass
        
        
        os.close(self.fd_stdout_r)
        os.close(self.fd_stderr_r)
        os.close(self.fd_stdin_r)
        os.close(self.fd_stdout_w)
        os.close(self.fd_stderr_w)
        os.close(self.fd_stdin_w)
            
                

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
                raise exceptions.LegacyException("The mpd daemon is not running, please make sure it is started before starting this code")
        
        if self._mpi_is_broken_after_possible_code_crash:
            raise exceptions.LegacyException("Another code has crashed, cannot spawn a new code, please stop the script and retry")
        
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
        
    @option(choices=("none","null","file"), sections=("channel",))
    def redirection(self):
        """Redirect the output of the code to null, standard streams or file"""
        return "null"
        
    @option(type="int", sections=("channel",))
    def max_message_length(self):
        """
        For calls to functions that can handle arrays, MPI messages may get too long for large N.
        The MPI channel will split long messages into blocks of size max_message_length.
        """
        return 1000000
    
    @option(sections=("channel",))
    def redirect_file(self):
        return "code.out"
        
    @late
    def redirection_filenames(self):
        return {
            "none":None,
            "null":("/dev/null", "/dev/null", "/dev/null"), 
            "file":("/dev/null", self.redirect_file, self.redirect_file),
        }[self.redirection]
        
    @late
    def debugger_method(self):
        return self.DEBUGGERS[self.debugger]
    
    def start(self):
        fd_stdin = None
        
        if self.debug_with_gdb or (not self.debugger_method is None):
            if not 'DISPLAY' in os.environ:
                arguments = None
                command = self.full_name_of_the_worker
            else:
                if self.debugger is None:
                    command, arguments = self.GDB(self.full_name_of_the_worker)
                else:
                    command, arguments = self.debugger_method(self.full_name_of_the_worker)
        else:
            arguments = None
            command = self.full_name_of_the_worker
            
            RunStandardStreamsRedirection.instance().before_spawn() 
            if False:
                if not self.redirection_filenames is None:
                    input_filename, output_filename, error_filename = self.redirection_filenames
                    print "redirecting...", self.redirection_filenames
                    fd_stdin = os.dup(0)
                    print fd_stdin
                    zero1 = os.pipe()
                    os.dup2(zero1.fileno(), 0)
                    
                    fd_stdout = os.dup(1)
                    zero2 = open(output_filename,'a')
                    os.dup2(zero2.fileno(), 1)
                    
                    fd_stderr = os.dup(2)
                    zero3 = open(error_filename,'a')
                    os.dup2(zero3.fileno(), 2)
            
        try:
            self.intercomm = MPI.COMM_SELF.Spawn(command, arguments, self.number_of_workers, info = self.info)
        finally:
            RunStandardStreamsRedirection.instance().after_spawn() 
            if False:
                if not self.redirection_filenames is None and not fd_stdin is None:
                    os.dup2(fd_stdin, 0)
                    os.dup2(fd_stdout, 1)
                    os.dup2(fd_stderr, 2)
                    os.close(fd_stdin)
                    os.close(fd_stdout)
                    os.close(fd_stderr)
                    zero1.close()
                    zero2.close()
                    zero3.close()
            
    def stop(self):
        if not self.intercomm is None:
            try:
                self.intercomm.Disconnect()
            except MPI.Exception as ex:
                if ex.error_class == MPI.ERR_OTHER:
                    type(self)._mpi_is_broken_after_possible_code_crash = True
                
            self.intercomm = None
    
    def determine_length_from_data(self, doubles_in, ints_in, floats_in, chars_in):
        def get_length(x):
            if x:
                try:
                    if not isinstance(x[0], str):
                        return len(x[0])
                except:
                    return -1
               
        lengths = map(get_length, (doubles_in, ints_in, floats_in, chars_in))
        
        return max(1, max(lengths))
        
        
    def send_message(self, tag, id=0, doubles_in=[], ints_in=[], floats_in=[], chars_in=[], length = 1):
    
        
        if self.is_inuse():
            raise exceptions.LegacyException("You've tried to send a message to a code that is already handling a message, this is not correct")
        if self.intercomm is None:
            raise exceptions.LegacyException("You've tried to send a message to a code that is not running")
        
        length = self.determine_length_from_data(doubles_in, ints_in, floats_in, chars_in)
        
        if length > self.max_message_length:
            self.split_message(tag, id, doubles_in, ints_in, floats_in, chars_in, length)
            return
        message = ServerSideMessage(tag,length)
        message.doubles = pack_array( doubles_in, message.length, 'float64')
        message.floats = pack_array(floats_in, message.length, 'float32')
        message.ints = pack_array(ints_in, message.length, 'int32')
        
        if message.length > 1:
            message.strings = pack_array(chars_in, message.length, 'string')
        else:
            message.strings = chars_in
    
        message.send(self.intercomm)
        self._is_inuse = True
        
    def split_message(self, tag, id, doubles_in, ints_in, floats_in, chars_in, length):
        doubles_result=[]
        ints_result=[]
        floats_result=[]
        strings_result=[]
        
        def append_results(arr, values, datatype):
            for j, subarr in enumerate(values):
                if datatype == 'string':
                    if not i: arr.append(subarr)
                    else: arr[j].extend(subarr)
                else:
                    if not i: arr.append(numpy.empty([length], dtype=datatype))
                    arr[j][i*self.max_message_length:(i+1)*self.max_message_length]=subarr
        
        def split_input_array(arr_in):
            if length == 1:
                return [tmp[i*self.max_message_length] for tmp in arr_in]
            else:
                return [tmp[i*self.max_message_length:(i+1)*self.max_message_length] for tmp in arr_in]
        
        for i in range(1+(length-1)/self.max_message_length):
            (doubles, ints, floats, chars) = map(split_input_array,[doubles_in, ints_in, floats_in, chars_in])
            self.send_message(tag, id, doubles, ints, floats, chars, min(length,self.max_message_length))
            next_part = self.recv_message(id, True)
            map(append_results,(doubles_result,ints_result,floats_result,strings_result),next_part,
                ('float64','int32','float32','string'))
            length -= self.max_message_length
        
        self._communicated_splitted_message = True
        self._merged_results_splitted_message = (doubles_result,ints_result,floats_result,strings_result)
    
    def recv_message(self, tag, handle_as_array):
        
        self._is_inuse = False
        
        if self._communicated_splitted_message:
            x = self._merged_results_splitted_message
            self._communicated_splitted_message = False
            del self._merged_results_splitted_message
            return x
        
        message = ServerSideMessage()
        try:
            message.recieve(self.intercomm)
        except MPI.Exception as ex:
            self.stop()
            raise
            
        if message.tag == -1:
            raise exceptions.LegacyException("Not a valid message, message is not understood by legacy code")
        elif message.tag == -2:
            self.stop()
            raise exceptions.LegacyException("Fatal error in code, code has exited")
                
            
            
        if message.length > 1 or handle_as_array:
            doubles_result = unpack_array(message.doubles, message.length, 'float64')
            floats_result = unpack_array(message.floats, message.length, 'float32')
            ints_result = unpack_array(message.ints, message.length, 'int32')
            strings_result = unpack_array(message.strings, message.length, 'string')
        else:
            doubles_result = message.doubles
            floats_result = message.floats
            ints_result = message.ints
            strings_result = message.strings
        
        return (doubles_result, ints_result, floats_result, strings_result)
        
    def nonblocking_recv_message(self, tag, handle_as_array):
        request = ServerSideMessage().nonblocking_recieve(self.intercomm)
        
        def handle_result(function):
            self._is_inuse = False
        
            message = function()
            
            if message.tag < 0:
                raise exceptions.LegacyException("Not a valid message, message is not understood by legacy code")
                
            if message.length > 1 or handle_as_array:
                doubles_result = unpack_array(message.doubles, message.length, 'float64')
                floats_result = unpack_array(message.floats, message.length, 'float32')
                ints_result = unpack_array(message.ints, message.length, 'int32')
                strings_result = unpack_array(message.strings, message.length, 'string')
            else:
                doubles_result = message.doubles
                floats_result = message.floats
                ints_result = message.ints
                strings_result = message.strings
            
            return (doubles_result, ints_result, floats_result, strings_result)
    
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
            
        template = """from amuse.support.legacy import channel
m = channel.MultiprocessingMPIChannel('{0}',number_of_workers = {1})
m.run_mpi_channel('{2}')"""
        code_string = template.format(
            self.full_name_of_the_worker, 
            self.number_of_workers, 
            self.name_of_the_socket
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
        channel = MpiChannel(self.full_name_of_the_worker, number_of_workers = self.number_of_workers)
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
         
    def send_message(self, tag, id=0, doubles_in=[], ints_in=[], floats_in=[], chars_in=[], length = 1):
        self._send(self.client_socket, ('send_message',(tag, id, doubles_in, ints_in, floats_in, chars_in, length),))
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
        
            

