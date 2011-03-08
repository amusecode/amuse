# Echo client program
import socket
import struct
import array
import numpy

HOST = 'localhost'    # The remote host
PORT = 61575              # The same port as used by the server
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.connect((HOST, PORT))
#s.send('Hello, world')
list = [1,2,3,4]
array1 = array.array('I',list)
numpyarray = numpy.array(list, dtype='int32')
s.send(numpyarray)
data = s.recv(1024)
s.close()

result = array.array("I")
result.fromstring(data)

print 'Received', result
