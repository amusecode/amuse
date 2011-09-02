from amuse.community import *
from amuse.units import units

class AarsethZareInterface(CodeInterface, LiteratureReferencesMixIn):
    """
    Interface to the regularized Burlish-Stoer integrator of
    Aarseth & Zare

    The relevant references are:
        .. [#] Aarseth, S. & Zare, K., 1974, Celestial Mechanics 10, 185.
        .. [#] Aarseth, S. & Zare, K., 1974, Celestial Mechanics 10, 516.
    """    
    include_headers = ['worker_code.h']

    def __init__(self, **keyword_arguments):
        CodeInterface.__init__(self, name_of_the_worker = 'aarsethzare_worker', **keyword_arguments)
        LiteratureReferencesMixIn.__init__(self)

    @legacy_function
    def evolve_triple():
        function = LegacyFunctionSpecification()  
        function.addParameter('time', dtype='float64', direction=function.INOUT)
        function.addParameter('masses', dtype='float64', direction=function.IN)
        function.addParameter('x', dtype='float64', direction=function.INOUT)
        function.addParameter('y', dtype='float64', direction=function.INOUT)
        function.addParameter('z', dtype='float64', direction=function.INOUT)
        function.addParameter('vx', dtype='float64', direction=function.INOUT)
        function.addParameter('vy', dtype='float64', direction=function.INOUT)
        function.addParameter('vz', dtype='float64', direction=function.INOUT)
        function.addParameter('tend', dtype='float64', direction=function.IN)
        function.addParameter('nl', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

    @legacy_function
    def construct_orbital_elements() :
        function = LegacyFunctionSpecification()  
        function.addParameter('m', dtype='float64', direction=function.IN)
        function.addParameter('r', dtype='float64', direction=function.IN)
        function.addParameter('v', dtype='float64', direction=function.IN)
        function.addParameter('e1', dtype='float64', direction=function.OUT)
        function.addParameter('e2', dtype='float64', direction=function.OUT)
        function.addParameter('nl', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function



    @legacy_function
    def construct_orbital_coordinates() :
        function = LegacyFunctionSpecification()  
        function.addParameter('m', dtype='float64', direction=function.IN)
        function.addParameter('e1', dtype='float64', direction=function.IN)
        function.addParameter('e2', dtype='float64', direction=function.IN)
        function.addParameter('r', dtype='float64', direction=function.OUT)
        function.addParameter('v', dtype='float64', direction=function.OUT)
        function.addParameter('nl', dtype='int32', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function

class AarsethZare(InCodeComponentImplementation):

    def __init__(self, unit_converter=None,**keyword_arguments):
        self.unit_converter = unit_converter
        InCodeComponentImplementation.__init__(self,  AarsethZareInterface(**keyword_arguments))
        
        
    def define_converter(self, object):
        if not self.unit_converter is None:
            object.set_converter(self.unit_converter.as_converter_from_si_to_generic())
            
    def define_particle_sets(self, object):
        object.define_inmemory_set('particles')
        
    def define_methods(self, handler):
        handler.add_method(
            "evolve_triple",
            (
                nbody_system.time,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.time
            ),
            (
                nbody_system.time,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "construct_orbital_elements",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.speed
            ),
            (
                units.none,
                units.none,
                handler.ERROR_CODE
            )
        )

        handler.add_method(
            "construct_orbital_coordinates",
            (
                nbody_system.mass,
                handler.NO_UNIT,
                handler.NO_UNIT
            ),
            (
                nbody_system.length,
                nbody_system.speed,
                handler.ERROR_CODE
            )
        )

    def evolve_model(self, tend):
        mass = self.particles.mass
        x = self.particles.x
        y = self.particles.y
        z = self.particles.z
        vx = self.particles.vx
        vy = self.particles.vy
        vz = self.particles.vz
        
        if hasattr(self.particles,"time"):
            time = self.particles.time
        else:
            time = tend.as_vector_with_length(len(self.particles)).aszeros()
        
        print time,mass,x,y,z,vx,vy,vz,tend.as_vector_with_length(len(self.particles))
        time, x, y, z, vx, vy, vz = self.evolve_triple(time,mass,x,y,z,vx,vy,vz,tend.as_vector_with_length(len(self.particles)))
        
        self.particles.time = time
        self.particles.x= x
        self.particles.y= y
        self.particles.z=z
        self.particles.vx=vx
        self.particles.vy=vy
        self.particles.vz=vz
        
