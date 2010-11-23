from amuse.support.core import late
from amuse.support import exceptions
import os.path

import textwrap
import struct
import numpy



registered_fileformat_processors = {}

class IoException(exceptions.CoreException):
    formatstring = "IO exception: {0}"
    
class UnsupportedFormatException(IoException):
    """Raised when the given format is not supported by AMUSE.
    """
    formatstring = "You tried to load or save a file with fileformat '{0}', but this format is not in the supported formats list"

class CannotSaveException(IoException):
    """Raised when the given format cannot save data (only reading of data is supported for the format)
    """
    formatstring = "You tried to save a file with fileformat '{0}', but this format is not supported for writing files"

class CannotLoadException(IoException):
    """Raised when the given format cannot read data (only saving of data is supported for the format)
    """
    formatstring = "You tried to load a file with fileformat '{0}', but this format is not supported for reading files"

    
class format_option(late):
    
    def __init__(self, initializer):
        late.__init__(self, initializer)
        self.__doc__ = initializer.__doc__
        
    def get_name(self):
        return self.initializer.__name__
       

def _get_processor_factor(format):
    if isinstance(format, basestring):
        if not format in registered_fileformat_processors:
            raise UnsupportedFormatException(format)
        processor_factory = registered_fileformat_processors[format]
    else:
        processor_factory = format
    
    return processor_factory
    
def write_set_to_file(
    set,
    filename,
    format = 'csv', 
    **format_specific_keyword_arguments):
    """
    Write a set to the given file in the given format.
    
    :argument filename: name of the file to write the data to
    :argument format: name of a registered format or
        a :class:`FileFormatProcessor` subclass (must be a 
        class and not an instance)
    
    All other keywords are set as attributes on the fileformat processor. To
    determine the supported options for a processor call 
    :func:`get_options_for_format`
    """
        
    processor_factory = _get_processor_factor(format)
    
    processor = processor_factory(filename, set=set, format = format)
    processor.set_options(format_specific_keyword_arguments)
    processor.store()
    

def read_set_from_file(
    filename,
    format = 'csv', 
    **format_specific_keyword_arguments):
    """
    Read a set from the given file in the given format.
    
    :argument filename: name of the file to read the data from
    :argument format: name of a registered format or
        a :class:`FileFormatProcessor` subclass (must be a 
        class and not an instance)
    
    All other keywords are set as attributes on the fileformat processor. To
    determine the supported options for a processor call 
    :func:`get_options_for_format`
    """
    if not os.path.exists(filename):
        raise IoException("Error: file '{0}' does not exist.".format(filename))
        
    
    processor_factory = _get_processor_factor(format)
    
    processor = processor_factory(filename, format = format)
    processor.set_options(format_specific_keyword_arguments)
    return processor.load()

def get_options_for_format(
        format = 'csv', 
    ):
    """Retuns a list of tuples, each tuple contains the
    name of the option, a description of the option and
    the default values.
    
    :argument format: name of a registered format or
        a :class:`FileFormatProcessor` subclass (must be a 
        class and not an instance)
    """
    
    processor_factory = _get_processor_factor(format)
    
    processor = processor_factory(format = format)
    
    return list(processor.get_description_of_options())



def add_fileformat_processor(class_of_the_format):
    """
    Register the specified class, so that it can be used
    by the :func:`write_set_to_file` and :func:`read_set_from_file`
    functions.
    
    Do not call this method directly, instead use :func:`FileFormatProcessor.register`
    """
    for x in class_of_the_format.provided_formats:
        registered_fileformat_processors[x] = class_of_the_format
    _update_documentation_strings()
    

    
def _update_documentation_strings():
    for methodname in ['write_set_to_file', 'read_set_from_file']:
        method = globals()[methodname]
        if not hasattr(method, '_original_doc'):
            method._original_doc = method.__doc__
        
        new_doc = method._original_doc
        
        new_doc += "\n    Registered file formats:\n\n"
        sorted_formatnames = sorted(registered_fileformat_processors.keys())
        for x in sorted_formatnames:
            processor = registered_fileformat_processors[x]
            processor_doc = processor.__doc__
            if processor_doc is None or len(processor_doc) == 0:
                continue
            processor_doc = processor_doc.strip()
            line = processor_doc.splitlines()[0]
            line = '    **' + x + '**,\n      ' + line + '\n'
            new_doc += line
        method.__doc__ = new_doc
    
class FileFormatProcessor(object):
    """
    Abstract base class of all fileformat processors
    
    All classes providing loading or storing of files should be
    subclasses of this base class.
    
    Every subclass must support the *filename*, *set* and
    *format* arguments. The arguments must all be optional.
    
    :argument filename: name of the file the read the data from
    :argument set: set (of particles or entities) to store in the file
    :argument format: format of the file, will be a string or class
    
    :attribute provided_formats: list of strings of the formats provided
        by the processor
    """
    
    provided_formats = []
    
    
    def __init__(self, filename = None, set = None, format = None):
        self.filename = filename
        self.set = set
        self.format = format
        
    @classmethod
    def get_options(cls):
        attribute_names = dir(cls)
        result = {}
        for x in attribute_names:
            if x.startswith('_'):
                continue
                
            attribute_value = getattr(cls, x)
            if isinstance(attribute_value, format_option):
                result[x] = attribute_value
        return result
        
    @classmethod
    def register(cls):
        """
        Register this class, so that it can be found by name
        int the :func:`write_set_to_file` and :func:`read_set_from_file`
        functions.
        """
        add_fileformat_processor(cls)
        
    def set_options(self, dictionary):
        supported_options = self.get_options()
        for key, value in dictionary.iteritems():
            if key in supported_options:
                setattr(self, key, value)
                
    def store(self):
        """
        Stores the set in the file.
        The set and the file are both properties
        of the processor.
        """
        raise CannotSaveException(self.format)
                
    def load(self):
        """
        Loads the set from the file and returns
        the set.
        """
        raise CannotLoadException(self.format)
        
    def store_string(self):
        raise CannotSaveException(self.format)
                
    def load_string(self, string):
        raise CannotLoadException(self.format)
    
    def get_description_of_options(self):
        """Yields tuples, each tuple contains the
        name of the option, a description of the option and
        the default values
        """
        for option, method in self.get_options().iteritems():
            default_value = getattr(self, option)
            description = textwrap.dedent(method.__doc__)
            yield (option, description, default_value)
            
class FullTextFileFormatProcessor(FileFormatProcessor):
    """
    Abstract base class of all fileformat processors that process
    their data by first reading the complete text string
    
    Subclasses need to implement the
    :func:`store_string` and :func:`load_string` methods.
    
    """

    def store(self):
        with open(self.filename, 'w') as f:
            f.write(self.store_string())
                
    def load(self):
        with open(self.filename, 'r') as f:
            return self.load_string(f.read())
        
    def store_string(self):
        """Return a string representation of the particle set"""
        raise CannotSaveException(self.format)
                
    def load_string(self, string):
        """Return a particle set, read from the string"""
        raise CannotLoadException(self.format)
        

class BinaryFileFormatProcessor(FileFormatProcessor):
    """
    Abstract base class of all fileformat processors that process
    their data by first reading the complete text string
    
    Subclasses need to implement the
    :func:`store_file` and / or :func:`load_file` methods.
    
    """

    def store(self):
        with open(self.filename, 'wb') as f:
            self.store_file(f)
                
    def load(self):
        with open(self.filename, 'rb') as f:
            return self.load_file(f)
        
    def store_file(self, file):
        """Store the data on the opened file"""
        raise CannotSaveException(self.format)
                
    def load_file(self, string):
        """Return a particle set, read from the binary file"""
        raise CannotLoadException(self.format)
    

class FortranFileFormatProcessor(BinaryFileFormatProcessor):
    """
    Abstract base class of all fileformat processors that process
    their data by first reading fortran blocks
    
    Subclasses need to implement the
    :func:`store_file` and / or :func:`load_file` methods.
    
    """
    
    @format_option
    def endianness(self):
        return '@' #native
        
    @late
    def float_type(self):
        result = numpy.dtype(numpy.float32)
        if self.endianness == '@':
            return result
        else:
            return result.newbyteorder(self.endianness)
    
    @late
    def uint_type(self):
        result = numpy.dtype(numpy.uint32)
        if self.endianness == '@':
            return result
        else:
            return result.newbyteorder(self.endianness)
    
    @late
    def int_type(self):
        result = numpy.dtype(numpy.int32)
        if self.endianness == '@':
            return result
        else:
            return result.newbyteorder(self.endianness)
            
    def read_fortran_block(self, file):
        """Returns one block read from file. Checks if the 
        block is consistant. Result is an array of bytes
        """
        format = self.endianness+'I'
        bytes = file.read(4)
        if not bytes:
            return None
        length_of_block = struct.unpack(format, bytes)[0]
        result = file.read(length_of_block)
        bytes = file.read(4)
        length_of_block_after = struct.unpack(format, bytes)[0]
        if(length_of_block_after != length_of_block):
            raise base.IoException("Block is mangled sizes don't match before: {0}, after: {1}".format(length_of_block, length_of_block_after))
        return result
        
    def read_fortran_block_floats(self, file):
        bytes = self.read_fortran_block(file)
        return numpy.frombuffer(bytes, dtype=self.float_type)
        
    def read_fortran_block_uints(self, file):
        bytes = self.read_fortran_block(file)
        return numpy.frombuffer(bytes, dtype=self.uint_type)
        
    def read_fortran_block_ints(self, file):
        bytes = self.read_fortran_block(file)
        return numpy.frombuffer(bytes, dtype=self.int_type)
        
    def read_fortran_block_float_vectors(self, file, size = 3):
        result = self.read_fortran_block_floats(file)
        return result.reshape(len(result)/size,size)
        
    def write_fortran_block(self, file, bytes):
        format = self.endianness+'I'
        length_of_block = len(bytes)
        file.write(struct.pack(format, length_of_block))
        file.write(bytes)
        file.write(struct.pack(format, length_of_block))
    
    
