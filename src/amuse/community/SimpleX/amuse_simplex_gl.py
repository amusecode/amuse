import numpy

from amuse.community import *
from amuse_simplex import SimpleX as SimpleX_basic


class SimpleX(SimpleX_basic):            
    include_headers=['worker_gl.h']
    def __init__(self, **kwargs):
        LegacyInterface.__init__(self,debug_with_gdb=False,
           name_of_the_worker = 'simplex_worker_gl', **kwargs)
        
    @legacy_function
    def start_viewer():
        function = LegacyFunctionSpecification()  
        function.result_type = 'i'
        return function
    @legacy_function
    def stop_viewer():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'          
        return function
