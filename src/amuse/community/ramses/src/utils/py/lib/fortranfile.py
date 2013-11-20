import numpy

class FortranFile(file):
	"""
	Class for reading Fortran unformatted files.
	"""
	def __init__(self, fname, mode='rb', buf=0):
		file.__init__(self, fname, mode, buf)

	def read_fortran_record(self, dtype, endian="="):
		""" Read a FORTRAN record of given numpy dtype """

		mytype = numpy.dtype(dtype).newbyteorder(endian)
		mint32 = numpy.dtype('i4').newbyteorder(endian)

		nbytes = numpy.array(0,dtype=mytype).itemsize

		n1 = numpy.fromfile(self, dtype=mint32, count=1)[0]/nbytes
		data = numpy.fromfile(self, dtype=mytype, count=n1)
		n2 = numpy.fromfile(self, dtype=mint32, count=1)[0]/nbytes


		if(n1 != n2):
			raise IOError("Error reading FORTRAN binary record")

		return data

