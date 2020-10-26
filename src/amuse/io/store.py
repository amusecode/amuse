
import numpy
import pickle
import os.path

from amuse.io import base
from amuse.units import si
from amuse.units import units
from amuse.units import core
from amuse.units.quantities import is_quantity
from amuse.support import exceptions

from amuse.datamodel import Particles
from amuse.datamodel import AttributeStorage

from amuse.io import store_v1
from amuse.io import store_v2

StoreHDF = store_v1.StoreHDF

class _FileContext(object):
    
    def __init__(self, processor, result):
        self.processor = processor
        self.result = result
    
    def __enter__(self):
        return self.result
    
    def __exit__(self, type, value, traceback):
        self.processor.close()
        
class HDF5FileFormatProcessor(base.FileFormatProcessor):
    """
    Process an HDF5 file
    """
    
    provided_formats = ['hdf5', 'amuse']
    
    def __init__(self, filename = None, set = None, format = None):
        base.FileFormatProcessor.__init__(self, filename, set, format)
    
    def load(self):
        processor, result = self.load_base()
        if self.return_context:
            return _FileContext(processor, result)
        else:
            return result
            
    def load_base(self):
        processor = store_v2.StoreHDF(
                self.filename, 
                open_for_writing = False, 
                append_to_file = self.append_to_file or self.allow_writing, 
                copy_history = self.copy_history,
                overwrite_file = self.overwrite_file )
        if not processor.is_correct_version():                
            processor.close()
            processor = store_v1.StoreHDF(
                    self.filename, 
                    open_for_writing = False, 
                    append_to_file = self.append_to_file or self.allow_writing, 
                    copy_history = self.copy_history,
                    overwrite_file = self.overwrite_file )

        if len(self.names) > 0:
            result = processor.load_sets(self.names)
            if self.close_file:
                if not self.copy_history:
                    for part in result:
                        part._private.previous = None
                processor.close()
        else:
            result = processor.load()
            if self.close_file:
                if not self.copy_history:
                    result = result.copy()
                    result._private.previous = None
                processor.close()
        return processor, result
        

    def store(self):
        
        if self.version == '1.0':
            processor = store_v1.StoreHDF(
                self.filename, 
                append_to_file = self.append_to_file, 
                open_for_writing = True,
                overwrite_file = self.overwrite_file
            )
            
            if not processor.is_correct_version():
                raise Exception("You are trying to append to a file that was not written in version 1.0 format")
        else:
            processor = store_v2.StoreHDF(
                self.filename, 
                append_to_file = self.append_to_file, 
                open_for_writing = True,
                overwrite_file = self.overwrite_file            
            )
        
            if not processor.is_correct_version():
                raise Exception("You are trying to append to a file that was written in version 1.0 format")
        try:
            if len(self.names) > 0:
                return processor.store_sets(self.set, self.names, self.extra_attributes)
            else:
                return processor.store(self.set, self.extra_attributes)
        finally:
            processor.close()
    
    @base.format_option
    def append_to_file(self):
        """If set to True, new data is appended to HDF5 files. 
        If set to False, the existing file is removed and overwritten.
        Only relevant for write set to file. (default: False)"""
        return False
    
    @base.format_option
    def close_file(self):
        """If set to True, the file is closed after reading, unless you
        set copy_history to True no previous savepoints will be returned"""
        return True
        
    @base.format_option
    def copy_history(self):
        """If set to True, the savepoint history is read from file into memory."""
        return True
        
    @base.format_option
    def names(self):
        """A list of names to save the data under. If filled this 
        will load the sets or grids with the given names and
        return this as a list. When saving the names will be used to
        save each set, a list of sets is needed in the write_set_to_file 
        function (default: [])"""
        return []
    

    @base.format_option
    def version(self):
        """AMUSE storage version to use, needs to be >= '2.0' if you want
        to store links between particles and grids (default: '2.0')"""
        return '2.0'

    @base.format_option
    def return_context(self):
        """If set to True, will return a context manager instead of
        the loaded set(s). This context manager will take care of properly 
        closing any connected resources. To access the set, use the with 
        statement.
        
        .. code-block:: python
            
            with load_set_from_file("example.h5", "amuse") as particles:
                print particles
        
        Usefull for cases where close_file == False. (default: False)"""
        return False
        
    @base.format_option
    def allow_writing(self):
        """If set to True, data can be written to the file, equivalent of append_to_file 
        even if read_set_from_file is used"""
        return False

    @base.format_option
    def overwrite_file(self):
        """If set to True, overwrite file if it exists, otherwise writing an existing
        file generates an exception"""
        return False
