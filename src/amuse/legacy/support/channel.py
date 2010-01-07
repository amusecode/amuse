import inspect
import numpy
import os
import os.path
import cPickle as pickle


import sys
import struct

from mpi4py import MPI
from subprocess import Popen


class MessageChannel(object):
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
    DEBUGGER = None
    
    
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
    def XTERM(cls, full_name_of_the_worker):
        arguments = ['-hold', '-display', os.environ['DISPLAY'], '-e', full_name_of_the_worker]
        command = 'xterm'
        return command, arguments
        
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
                if current_type.__bases__[0] is object:
                    raise Exception("The worker application does not exists, it should be at: {0}".format(tried_workers))
            else:
                found = True
        return full_name_of_the_worker

    def send_message(self, tag, id=0, int_arg1=0, int_arg2=0, doubles_in=[], ints_in=[], floats_in=[], chars_in=[], length = 1):
        pass

    
    def recv_message(self, tag, handle_as_array = False):
        pass
        
#import time
#import ctypes
#clib_library = ctypes.CDLL("libc.so.6")
#memcpy = clib_library.memcpy
#memcpy.argtypes = [ctypes.c_void_p, ctypes.c_void_p, ctypes.c_size_t]

class MpiChannel(MessageChannel):
    """
    Message channel based on MPI calls to send and recv the messages
    
    :argument name_of_the_worker: Name of the application to start
    :argument number_of_workers: Number of parallel processes
    :argument legacy_interface_type: Type of the legacy interface
    :argument debug_with_gdb: If True opens an xterm with a gdb to debug the remote process
    :argument hostname: Name of the node to run the application on
    """
    def __init__(self, name_of_the_worker, number_of_workers, legacy_interface_type = None, debug_with_gdb = False, hostname = None):
        self.name_of_the_worker = name_of_the_worker
        self.number_of_workers = number_of_workers
        
        if not legacy_interface_type is None:
            self.full_name_of_the_worker = self.get_full_name_of_the_worker( legacy_interface_type)
        else:
            self.full_name_of_the_worker = self.name_of_the_worker
            
        self.debug_with_gdb = debug_with_gdb
            
        if not hostname is None:
            self.info = MPI.Info.Create()
            self.info['host'] = hostname
        else:
            self.info = MPI.INFO_NULL
            
        self.cached = None
        self.intercomm = None
        
    def start(self):
        if self.debug_with_gdb or (not self.DEBUGGER is None):
            if not 'DISPLAY' in os.environ:
                arguments = None
                command = self.full_name_of_the_worker
            else:
                if self.DEBUGGER is None:
                    command, arguments = self.GDB(self.full_name_of_the_worker)
                else:
                    command, arguments = self.DEBUGGER(self.full_name_of_the_worker)
        else:
            arguments = None
            command = self.full_name_of_the_worker
        self.intercomm = MPI.COMM_SELF.Spawn(command, arguments, self.number_of_workers, info = self.info)

    def stop(self):
        if not self.intercomm is None:
            self.intercomm.Free()
            self.intercomm = None
        
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
                
            self.intercomm.Bcast([offsets, MPI.INT], root=MPI.ROOT)    
            bytes = []
            for strings in chars_in:
                if length == 1:
                    bytes.extend([ord(ch) for ch in strings])
                    bytes.append(0)
                else:
                    for string in strings:
                        bytes.extend([ord(ch) for ch in string])
                        bytes.append(0)
              
            chars = numpy.array(bytes, dtype=numpy.uint8)
            self.intercomm.Bcast([chars, MPI.CHARACTER], root=MPI.ROOT)
       
        
        
    def recv_message(self, tag, handle_as_array):
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
            if length > 1 or handle_as_array:
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
            if length > 1 or handle_as_array:
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
        
    def is_active(self):
        return self.intercomm is not None
        


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
    def __init__(self, name_of_the_worker, number_of_workers, legacy_interface_type, debug_with_gdb = False, hostname = None):
        self.name_of_the_worker = name_of_the_worker
        self.number_of_workers = number_of_workers
        if not legacy_interface_type is None:
            self.full_name_of_the_worker = self.get_full_name_of_the_worker( legacy_interface_type)
        else:
            self.full_name_of_the_worker = self.name_of_the_worker
        self.debug_with_gdb = debug_with_gdb
        self.process = None
    
    def start(self):
        name_of_dir = "/tmp/amuse_"+os.getenv('USER')
        self.name_of_the_socket, self.server_socket = self._createAServerUNIXSocket(name_of_dir)
        environment = os.environ.copy()
        
        if 'PYTHONPATH' in environment:
            environment['PYTHONPATH'] = environment['PYTHONPATH'] + ':' +  self._extra_path_item(__file__)
        else:
            environment['PYTHONPATH'] =  self._extra_path_item(__file__)
            
        template = """from amuse.legacy.support import core
m = core.MultiprocessingMPIChannel('{0}',{1},None,{2})
m.run_mpi_channel('{3}')"""
        code_string = template.format(
            self.full_name_of_the_worker, 
            self.number_of_workers, 
            self.debug_with_gdb, 
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
        
            

