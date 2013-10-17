from amuse.community import *
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSubset
from amuse.community.interface import se

class SeBaInterface(CodeInterface, se.StellarEvolutionInterface, LiteratureReferencesMixIn):
    
    """
    Stellar evolution is performed by the rapid single-star evolution
    and binary evolution using SeBa.This is a package of
    semi-analytical formulae which covers all phases of evolution from
    the zero-age main-sequence up to and including remnant phases. It
    is valid for masses in the range 0.01-1000 Msun with variable
    metallicity.  SeBa includes prescriptions for mass loss by stellar
    winds, supernova and supports binary evolution.
    
        .. [#] ** Portegies Zwart S.F. & Verbunt F., 1996, A&A, 309, 179:
        .. [#] ... "Population synthesis of high-mass binaries"
        .. [#] Toonen, S., Nelemans, G., Portegies Zwart S.F., 2012 submitted to A&A (arXiv 1101.2787)
        .. [#] ... "New population synthesis model: Preliminary results for close double white dwarf populations"
    """

    include_headers = ['worker_code.h']
    
    def __init__(self, **options):
        CodeInterface.__init__(self, name_of_the_worker="seba_worker", **options)
        LiteratureReferencesMixIn.__init__(self)

    
    @legacy_function
    def evolve_star():
        function = LegacyFunctionSpecification()  
        function.addParameter('mass', dtype='float64', direction=function.IN)
        function.addParameter('endtime', dtype='float64', direction=function.IN)
        function.addParameter('metal', dtype='float64', direction=function.IN)
        function.addParameter('resulttime', dtype='float64', direction=function.OUT)
        function.addParameter('end_mass', dtype='float64', direction=function.OUT)
        function.addParameter('end_radius', dtype='float64', direction=function.OUT)
        function.addParameter('end_luminosity', dtype='float64', direction=function.OUT)
        function.addParameter('end_temperature', dtype='float64', direction=function.OUT)
        function.addParameter('time_step', dtype='float64', direction=function.OUT)
        function.addParameter('stellar_type', dtype='int32', direction=function.OUT)
        function.result_type = 'int32'
        function.can_handle_array = True
        return function

    @legacy_function   
    def new_binary():
        """
        Define a new star in the code. The star will start with the given mass.
        """
        function = LegacyFunctionSpecification()  
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.OUT
            , description="The new index for the star. This index can be used to refer to this star in other functions")
        function.addParameter('semi_major_axis', dtype='float64', direction=function.IN
            , description="The eccentricity of the orbit")
        function.addParameter('eccentricity', dtype='float64', direction=function.IN
            , description="The eccentricity of the orbit")
        function.addParameter('child1', dtype='int32', direction=function.IN
            , description="The index of the first child, as returned by new_particle")
        function.addParameter('child2', dtype='int32', direction=function.IN
            , description="The index of the second child, as returned by new_particle")
            
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            New star was loaded and the index_of_the_star parameter set.
        -1 - ERROR
            New star could not be created.
        """
        return function
        
    @legacy_function
    def delete_binary():
        """
        Remove the definition of binary from the code. After calling this function the particle is
        no longer part of the model evolution. It's children are still a part of particles model.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the binary to be removed. This index must have been returned by an earlier call to :meth:`new_binary`")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            binary was removed from the model
        -1 - ERROR
            binary could not be found
        -2 - ERROR
            not yet implemented
        """
        return function
        
    @legacy_function
    def evolve_system():
        """
        Evolve the model until the given time, or until a stopping condition is set.
        Need to call this evolve_system as evolve_model is overriden in se.StellarEvolution
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN,
            description = "Model time to evolve the code to. The model will be "
                "evolved until this time is reached exactly or just after.")
        function.result_type = 'int32'
        return function  
    
    @legacy_function   
    def get_eccentricity():
        """
        Retrieve the current eccentricity of the binary star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The current eccentricity.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function
    
    @legacy_function   
    def get_semi_major_axis():
        """
        Retrieve the current semi major axis of the elliptical orbit of the parts in the binary star.
        """
        function = LegacyFunctionSpecification() 
        function.can_handle_array = True 
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The current semi major axis.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function
        
    @legacy_function
    def get_children_of_binary():
        """
        Return the indices of both children
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_star', dtype='int32',
                              direction=function.IN, 
                 description = 'index of the parent particle',
                 unit = INDEX)
        function.addParameter('child1', dtype='int32', direction=function.OUT,
                description = 'index of the first child particle, -1 if none',
                unit = LINK('particles') )
        function.addParameter('child2', dtype='int32', direction=function.OUT,
                unit = LINK('particles'))
        function.can_handle_array = True
        function.result_type = 'int32'
        return function
        
    
    @legacy_function
    def get_is_logging_of_evolve_enabled():
        """
        If True log the star state before and after evolve
        in starev.data
        """
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was retrieved
        -1 - ERROR
            could not retrieve parameter
        """
        return function
        
    @legacy_function
    def set_is_logging_of_evolve_enabled():
        """
        If True log the star state before and after evolve
        in starev.data
        """
        function = LegacyFunctionSpecification()
        function.addParameter('value', dtype='bool', direction=function.IN)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            the parameter was set
        -1 - ERROR
            could not set parameter
        """
        return function

    @legacy_function
    def get_supernova_kick_velocity():
        """
        Retrieve the current value of the supernova kick velocity (in kms). 
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('v_disp', dtype='float64', direction=function.OUT,
            description = "The current value of the kick velocity dispersion")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the metallicity was retrieved
        -1 - ERROR
            The code does not have support for retrieving the metallicity
        """
        return function
        
    
    @legacy_function
    def set_supernova_kick_velocity():
        """
        Update the value of the kick velocity dispersion.
        """
        function = LegacyFunctionSpecification()  
        function.addParameter('v_disp', dtype='float64', direction=function.IN,
            description = "The new value of the supernova kick velocity dispersion.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the metallicity was set
        -1 - ERROR
            The code does not have support for updating the metallicity
        """
        return function

    def evolve_model(self, time):
        return self.evolve_system(time)

class SeBa(se.StellarEvolution):

    def __init__(self, **options):
        se.StellarEvolution.__init__(self,  SeBaInterface(**options), **options)
    
        self.model_time = 0.0 | units.yr

    
    def evolve_model(self, end_time):
        return self.evolve_system(end_time);

    def define_methods(self, object):
        se.StellarEvolution.define_methods(self, object)

        object.add_method(
            "evolve_for",
            (object.INDEX, units.Myr),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "evolve_star",
            (units.MSun, units.Myr, units.none),
            (units.Myr, units.MSun, units.RSun, units.LSun, units.K, units.Myr,units.stellar_type, object.ERROR_CODE)
        )
        object.add_method(
            "evolve_system",
            (units.Myr,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "new_binary",
            (units.RSun, object.NO_UNIT, object.LINK('particles'), object.LINK('particles')),
            (object.INDEX, object.ERROR_CODE,)
        )
        object.add_method(
            "delete_binary",
            (object.INDEX,),
            (object.ERROR_CODE,)
        )
        object.add_method(
            "get_eccentricity",
            (object.INDEX,),
            (object.NO_UNIT, object.ERROR_CODE,)
        )
        object.add_method(
            "get_semi_major_axis",
            (object.INDEX,),
            (units.RSun, object.ERROR_CODE,)
        )
        object.add_method(
            "get_age", 
            (object.INDEX,), 
            (units.Myr, object.ERROR_CODE,)
        )
        object.add_method(
            "get_time_step", 
            (object.INDEX,), 
            (units.Myr, object.ERROR_CODE,)
        )
        object.add_method(
            "get_supernova_kick_velocity", 
            (), 
            (units.kms, object.ERROR_CODE,)
        )
        object.add_method(
            "set_supernova_kick_velocity", 
            (units.kms,), 
            (object.ERROR_CODE,)
        )
    
    def update_time_steps(self):
        pass
    
    def define_parameters(self, object):
        object.add_method_parameter(
            "get_metallicity",
            "set_metallicity",
            "metallicity", 
            "Metallicity of all stats", 
            default_value = 0.02
        )
        object.add_method_parameter(
            "get_supernova_kick_velocity",
            "set_supernova_kick_velocity",
            "supernova_kick_velocity", 
            "Kick velocity to compact object formed in supernova", 
            default_value = 600 | units.kms
        )
        
        object.add_method_parameter(
            "get_is_logging_of_evolve_enabled",
            "set_is_logging_of_evolve_enabled",
            "is_logging_of_evolve_enabled", 
            "if True will log star state before and after evolve in starev.data", 
            default_value = False
        )
        
        
    def define_particle_sets(self, object):
       
        object.define_set('particles', 'index_of_the_star')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_star')
        
        object.add_getter('particles', 'get_radius', names = ('radius',))
        object.add_getter('particles', 'get_stellar_type', names = ('stellar_type',))
        object.add_getter('particles', 'get_mass', names = ('mass',))
        object.add_getter('particles', 'get_age', names = ('age',))
        object.add_getter('particles', 'get_time_step', names = ('time_step',))
        #object.add_getter('particles', 'get_spin', names = ('spin',))
        object.add_getter('particles', 'get_luminosity', names = ('luminosity',))
        object.add_getter('particles', 'get_temperature', names = ('temperature',))
        object.add_method('particles', 'evolve_one_step')
        object.add_method('particles', 'evolve_for')
    

        object.define_set('binaries', 'index_of_the_star')
        object.set_new('binaries', 'new_binary')
        object.set_delete('binaries', 'delete_binary')
        
        object.add_getter('binaries', 'get_semi_major_axis', names = ('semi_major_axis',))
        object.add_getter('binaries', 'get_eccentricity', names = ('eccentricity',))
        object.add_getter('binaries', 'get_mass', names = ('mass',))
        object.add_getter('binaries', 'get_time_step', names = ('time_step',))
        object.add_getter('binaries', 'get_age', names = ('age',))
        object.add_getter("binaries", 'get_children_of_binary')
