import numpy
import os.path
from amuse.community import *
from amuse.datamodel import Particles
from amuse import datamodel
from amuse.units import nbody_system
from amuse.community.interface.common import CommonCodeInterface, CommonCode

class FractalClusterInterface(CodeInterface,  LiteratureReferencesMixIn):
    """
    makes fractal of nstar particles of dimension fdim, using ndiv 
    subunits forcing the number of cells if force=.true.
    
    reference:
        .. [#] Simon Goodwin & Ant Whitworth (2004, A&A, 413, 929) [2004A&A...413..929G]

    """

    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker = self.name_of_the_worker(), **options)
        LiteratureReferencesMixIn.__init__(self)

    def name_of_the_worker(self):
        return 'fractal_worker'

    @legacy_function
    def get_state():
        function = LegacyFunctionSpecification()
        function.must_handle_array = True
        function.addParameter('id', dtype='i', direction=function.IN)
        for x in ['x','y','z','vx','vy','vz']:
            function.addParameter(x, dtype='d', direction=function.OUT)
        function.addParameter('nstar', dtype='i', direction=function.LENGTH)
        function.result_type = 'i'
        return function

    @legacy_function
    def generate_particles():
        function = LegacyFunctionSpecification()
        function.result_type = 'i'
        return function

    @legacy_function
    def get_fractal_dimension():
        function = LegacyFunctionSpecification()
        function.addParameter('fdim', dtype='d', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function
    def set_fractal_dimension():
        function = LegacyFunctionSpecification()
        function.addParameter('fdim', dtype='d', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_random_seed():
        function = LegacyFunctionSpecification()
        function.addParameter('seed', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function
    def set_random_seed():
        function = LegacyFunctionSpecification()
        function.addParameter('seed', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function

    @legacy_function
    def get_nstar():
        function = LegacyFunctionSpecification()
        function.addParameter('nstar', dtype='i', direction=function.OUT)
        function.result_type = 'i'
        return function
    @legacy_function
    def set_nstar():
        function = LegacyFunctionSpecification()
        function.addParameter('nstar', dtype='i', direction=function.IN)
        function.result_type = 'i'
        return function
    
    new_particle = None
    
    def delete_particle(self, index_of_the_particle):
        return 0
    
    @legacy_function
    def get_number_of_particles_updated():
        """
        Return the number of particles added during the last generate_particles.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        return function
    


class FractalCluster(CommonCode):
    
    def __init__(self, unit_converter=None, **options):
        self.unit_converter = unit_converter
        InCodeComponentImplementation.__init__(self, FractalClusterInterface(**options), **options)
    
    def initialize_code(self):
        pass
    
    def cleanup_code(self):
        pass
    
    def commit_parameters(self):
        pass
    
    def recommit_parameters(self):
        pass
    
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_nstar",
            "set_nstar",
            "number_of_particles",
            "the number of particles to be generated in the model",
            default_value = 0
        )
        
        handler.add_method_parameter(
            "get_fractal_dimension",
            "set_fractal_dimension",
            "fractal_dimension",
            "the fractal dimension of the spatial particle distribution",
            default_value = 1.6
        )
        
        handler.add_method_parameter(
            "get_random_seed",
            "set_random_seed",
            "random_seed",
            "the initial seed to be used by the random number generator",
            default_value = 1234321
        )
    
    def define_methods(self, handler):
        CommonCode.define_methods(self, handler)
        handler.add_method("generate_particles", (), (handler.ERROR_CODE,))
        handler.add_method("get_number_of_particles_updated", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        
        handler.add_method("get_state", (handler.INDEX,), 
            [nbody_system.length]*3 + [nbody_system.speed]*3 + [handler.ERROR_CODE]
        )
        
        handler.add_method("get_target_number_of_particles", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_target_number_of_particles", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_fractal_dimension", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_fractal_dimension", (handler.NO_UNIT, ), (handler.ERROR_CODE,))
        
        handler.add_method("get_random_seed", (), (handler.NO_UNIT, handler.ERROR_CODE,))
        handler.add_method("set_random_seed", (handler.NO_UNIT, ), (handler.ERROR_CODE,))

    def define_converter(self, handler):
        if not self.unit_converter is None:
            handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())
    
    def define_particle_sets(self, handler):
        handler.define_set('particles', 'index_of_the_particle')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_getter('particles', 'get_state')
    
    def define_state(self, handler):
        CommonCode.define_state(self, handler)
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('EDIT', 'before_get_parameter')
        
        handler.add_transition('EDIT', 'RUN', 'generate_particles', False)
        handler.add_transition('RUN', 'EDIT', 'clear_particle_set')
        handler.add_method('EDIT', 'get_number_of_particles_updated')
        handler.add_method('RUN', 'get_number_of_particles_updated')
        handler.add_method('RUN', 'get_state')
    
    def generate_particles(self):
        result = self.overridden().generate_particles()
        self.update_particle_set()
    
    def update_particle_set(self):
        """
        update the particle set after changes in the code
        
        this implementation needs to move to the
        amuse.datamodel.incode_storage module, as
        it uses a lot of internal methods and info!
        """
        number_of_updated_particles = self.get_number_of_particles_updated()
        if number_of_updated_particles:
            self.particles._private.attribute_storage._add_indices(
                list(range(1, number_of_updated_particles+1))
            )
    
    def clear_particle_set(self):
        if len(self.particles):
            self.particles.remove_particles(self.particles)
    


class MakeFractalCluster(object):
    
    def __init__(self, N=None, convert_nbody=None, masses=None, do_scale=True,
            random_seed=None, fractal_dimension=1.6, virial_ratio=0.5, verbose=False, match_N=True):
        if masses is None:
            if N is None:
                raise exceptions.AmuseException("Either keyword argument 'N' (number of particles) or "
                    "'masses' (vector quantity with mass of each particle) is required.")
            self.masses = numpy.ones(N) / N | nbody_system.mass
            self.N = N
        else:
            if not N is None and len(masses) != N:
                print("warning: provided mass array not equal to masses")
            self.masses = masses / masses.sum() | nbody_system.mass
            self.N = len(masses)
        
        self.convert_nbody=convert_nbody
        self.do_scale=do_scale
        self.random_seed=random_seed
        self.fractal_dimension=fractal_dimension
        self.virial_ratio=virial_ratio
        self.verbose=verbose
        self.match_N=match_N

    def new_model(self):
        generator = FractalCluster(redirection=("none" if self.verbose else "null"))
        generator.parameters.number_of_particles = self.N
        generator.parameters.fractal_dimension = self.fractal_dimension
        if self.random_seed is not None:
            generator.parameters.random_seed = self.random_seed
        generator.generate_particles()
        if self.match_N:
            while len(generator.particles)<self.N:
                generator.generate_particles()
        result = generator.particles.copy()
        generator.stop()
        result.mass = self.masses
        result.radius = 0 | nbody_system.length
        return result

    @property
    def result(self):
        particles = self.new_model()
        particles.move_to_center()
        if self.do_scale:
            particles.scale_to_standard(virial_ratio=self.virial_ratio)
        
        if not self.convert_nbody is None:
            return datamodel.ParticlesWithUnitsConverted(particles, self.convert_nbody.as_converter_from_si_to_generic()).copy()
        return particles

def new_fractal_cluster_model(*list_arguments, **keyword_arguments):
    """
    create a fractal stellar distribution

    :argument N: Number of particles, if None then masses must be provided [None]
    :argument masses: optional masses, if not provided masses will be 1/N [None]
    :argument fractal_dimension: fractal dimension of distribution, between 1.6 and 3 [1.6]
    :argument random_seed: random seed [None]
    :argument convert_nbody:  When given will convert the resulting set to SI units
    :argument do_scale: scale the result to exact nbody units (G=1, M=1, U=-0.5)
    :argument virial_ratio: ratio of kinetic to potential energy [0.5]
    """
    uc = MakeFractalCluster(*list_arguments, **keyword_arguments)
    return uc.result


Fractalcluster = FractalCluster
