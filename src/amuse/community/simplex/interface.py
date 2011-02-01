import os.path
from amuse.community.interface.common import CommonCodeInterface, CommonCode
from amuse.community import *
from amuse.support.options import option


class SimpleXInterface(LegacyInterface, CommonCodeInterface, LiteratureRefs):
    """
    SimpleX(2) is a method for radiative transfer on an unstructured Delaunay 
    grid. The grid samples the medium through which photons are transported in 
    an optimal way for fast radiative transfer calculations.
    
    The relevant references are:
        .. [#] Paardekooper, J.-P., Kruip, C.J.H., Icke, V. 2010, A&A, 515, A79 (SimpleX2)
        .. [#] Ritzerveld, J., & Icke, V. 2006, Phys. Rev. E, 74, 26704 (SimpleX)
    """
    include_headers=['worker.h']
    
    def __init__(self, **kwargs):
        LegacyInterface.__init__(self, name_of_the_worker = "simplex_worker", **kwargs)
        LiteratureRefs.__init__(self)
    
    @option(type="string")
    def data_directory(self):
        """
        The root name of the directory for the SimpleX
        application data files.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'simplex', 'input')
    
    @option(type="string")
    def output_directory(self):
        """
        The root name of the directory to use by the 
        application to store it's output / temporary files in.
        """
        return os.path.join(get_amuse_root_dir(), 'data', 'simplex', 'output')
    
    @legacy_function
    def set_output_directory():
        function = LegacyFunctionSpecification()
        function.addParameter('output_path', dtype='string', direction=function.IN,
            description = "Name of the output directory")
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def commit_particles():
        """
        Let the code perform initialization actions
        after all particles have been loaded in the model.
        Should be called before the first evolve call and
        after the last new_particle call.
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention
        """
        return function

    @legacy_function
    def recommit_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def new_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT)
        for x in ['x','y','z','rho','flux','xion']:
            function.addParameter(x, dtype='float64', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def delete_particle():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN)
        for x in ['x','y','z','rho','flux','xion']:
            function.addParameter(x, dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        return function
        
    @legacy_function
    def evolve():
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN)
        function.addParameter('synchronize', dtype='int32', direction=function.IN)
        function.result_type = 'int32'
        return function
    


class GlSimpleXInterface(SimpleXInterface):
    include_headers=['worker_gl.h']
    
    def __init__(self, **kwargs):
        LegacyInterface.__init__(self,debug_with_gdb=False,
           name_of_the_worker = 'simplex_worker_gl', **kwargs)
    
    @legacy_function
    def start_viewer():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    
    @legacy_function
    def stop_viewer():
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        return function
    

class SimpleXDoc(object):

    def __get__(self, instance, owner):
        return instance.legacy_interface.__doc__+"\n\n"+instance.parameters.__doc__


class SimpleX(CommonCode):
    
    __doc__ = SimpleXDoc()
    
    def __init__(self, **options):
        CodeInterface.__init__(self, SimpleXInterface(), **options)
        self.set_output_directory(self.output_directory)
