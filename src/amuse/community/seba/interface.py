from amuse.community import *
from amuse.datamodel import Particles
from amuse.datamodel import ParticlesSubset
from amuse.community.interface import se
from amuse.support import code

class SeBaInterface(CodeInterface, se.StellarEvolutionInterface, LiteratureReferencesMixIn, StoppingConditionInterface):

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
        .. [#] Toonen, S., Nelemans, G., Portegies Zwart S.F., 2012, A&A, 546A, 70T
        .. [#] ... "Supernova Type Ia progenitors from merging double white dwarfs. Using a new population synthesis model"
    """

    include_headers = ['worker_code.h', 'stopcond.h']

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
    def get_core_mass():
        """
        Retrieve the current core mass of a star
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The current core_mass.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function

    @legacy_function
    def get_COcore_mass():
        """
        Retrieve the current CO core mass of a star
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The current core_mass.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function

    @legacy_function
    def change_mass():
        """
        Set the add_mass increase (positive) or decrease (negative) of a star
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('value', dtype='float64', direction=function.IN
            , description="The current change_mass.")
        function.addParameter('dt', dtype='float64', direction=function.IN
            , description="The time in which the mass was changed.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function


    @legacy_function
    def merge_the_binary():
        """
        Merges the stars in the binary, creates a 'merged' binary. 
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_binary', dtype='int32', direction=function.IN
            , description="The index of the binary to set the value of")
        function.addParameter('child1', dtype='int32', direction=function.IN
            , description="The index of the consumer")
        function.addParameter('child2', dtype='int32', direction=function.IN
            , description="The index of the dinner")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function

    @legacy_function
    def merge_with_other_star():
        """
        Merges the star with another star, companion star remains unaltered
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('child1', dtype='int32', direction=function.IN
            , description="The index of the consumer")
        function.addParameter('child2', dtype='int32', direction=function.IN
            , description="The index of the dinner")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function



    @legacy_function
    def refresh_memory():
        """
        Refresh the memory of SeBa. Update previous parameters in SeBa to current values. 
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function


    @legacy_function
    def recall_memory_one_step():
        """
        Recall the memory of SeBa one time_step ago. Update current parameters in SeBa to previous values. 
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function



    @legacy_function
    def get_envelope_mass():
        """
        Retrieve the current envelope mass of a star
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The current envelope_mass.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function

    @legacy_function
    def get_core_radius():
        """
        Retrieve the current radius of the core
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('value', dtype='float64', direction=function.OUT
            , description="The current core_radius.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function


    @legacy_function
    def set_semi_major_axis():
        """
        Update the current semi major axis of the elliptical orbit of the parts in the binary star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('value', dtype='float64', direction=function.IN
            , description="The new semi major axis.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            The value has been set.
        -1 - ERROR
            A binary with the given index was not found.
        """
        return function

    @legacy_function
    def set_eccentricity():
        """
        Update the current eccentricity of the elliptical orbit of the parts in the binary star.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to get the value of")
        function.addParameter('value', dtype='float64', direction=function.IN
            , description="The new eccentricity.")
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
            Current value of the supernova kick velocity was retrieved
        -1 - ERROR
            The code does not have support for retrieving the supernova kick velocity
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
            Current value of the supernova kick velocity was set
        -1 - ERROR
            The code does not have support for updating the supernova kick velocity
        """
        return function


    @legacy_function
    def get_gyration_radius_sq():
        """
        Retrieve the current value of the square of the gyration radius (no units).
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('gyration_radius_sq', dtype='float64', direction=function.OUT,
            description = "The current value of the square of the gyration radius")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the square of the gyration radius was retrieved
        -1 - ERROR
            The code does not have support for retrieving the square of the gyration radius
        """
        return function

    @legacy_function
    def get_relative_age():
        """
        Retrieve the current value of the square of the relative age (Myr).
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('relative_age', dtype='float64', direction=function.OUT,
            description = "The current value of the square of the relative age")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the square of the relative age was retrieved
        -1 - ERROR
            The code does not have support for retrieving the square of the relative age
        """
        return function

    @legacy_function
    def get_natal_kick_velocity():
        """
        Retrieve the current value of the square of the natal kick velocity (Myr).
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('natal_kick_velocity_x', dtype='float64', direction=function.OUT)
        function.addParameter('natal_kick_velocity_y', dtype='float64', direction=function.OUT)
        function.addParameter('natal_kick_velocity_z', dtype='float64', direction=function.OUT)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the natal_kick_velocity was retrieved
        -1 - ERROR
            The code does not have support for retrieving the square of the natal_kick_velocity
        """
        return function

    @legacy_function
    def get_relative_mass():
        """
        Retrieve the current value of the square of the relative mass (MSun).
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('relative_mass', dtype='float64', direction=function.OUT,
            description = "The current value of the square of the relative mass")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the square of the relative mass was retrieved
        -1 - ERROR
            The code does not have support for retrieving the square of the relative mass
        """
        return function

    @legacy_function
    def get_effective_radius():
        """
        Retrieve the current value of the effective radius (Rsun).
        This can be different from the (equilibrium) radius due to accretion or mass loss.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('effective_radius', dtype='float64', direction=function.OUT,
            description = "The current value of the effective radius")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the effective_radius was retrieved
        -1 - ERROR
            The code does not have support for retrieving the effective_radius
        """
        return function


    @legacy_function
    def get_convective_envelope_mass():
        """
        Retrieve the current value of the mass of the part of the envelope that is convective (MSun).
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('convective_envelope_mass', dtype='float64', direction=function.OUT,
            description = "The current value of the mass of the part of the envelope that is convective")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the convective envelope mass was retrieved
        -1 - ERROR
            The code does not have support for retrieving the square of the relative mass
        """
        return function




    @legacy_function
    def get_convective_envelope_radius():
        """
        Retrieve the current value of the radius of the part of the envelope that is convective (RSun).
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('convective_envelope_radius', dtype='float64', direction=function.OUT,
            description = "The current value of the radius of the part of the envelope that is convective")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the convective envelope radius was retrieved
        -1 - ERROR
            The code does not have support for retrieving the square of the relative mass
        """
        return function




    @legacy_function
    def get_time():
        """
        Retrieve the model time. This time should be close to the end time specified
        in the evolve code. Or, when a collision was detected, it will be the
        model time of the collision.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT,
            description = "The current model time")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the time was retrieved
        -1 - ERROR
            The code does not have support for querying the time
        """
        return function

    @legacy_function
    def get_wind_mass_loss_rate():
        """
        Retrieve the current value of the wind_mass_loss_rate (Msun/yr).
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_star', dtype='int32', direction=function.IN
            , description="The index of the star to set the value of")
        function.addParameter('wind_mass_loss_rate', dtype='float64', direction=function.OUT,
            description = "The current value of the wind_mass_loss_rate")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the wind_mass_loss_rate was retrieved
        -1 - ERROR
            The code does not have support for retrieving the wind_mass_loss_rate
        """
        return function


    # SPZ&SLWMcM
    # No stopping conditions in this version for now.
    def evolve_model(self, end_time=None, keep_synchronous=True):
        if not keep_synchronous:
            raise Exception("non_synchronous evolution not implemented")
        if end_time is None:
            end_time = self.model_time + min(self.particles.time_step)
        return self.evolve_system(end_time)

class SeBa(se.StellarEvolution):
    
    __interface__ = SeBaInterface

    def __init__(self, **options):
        self.stopping_conditions = StoppingConditions(self)
        self.stopping_conditions.supernova_detection = code.StoppingCondition('supernova_detection')

        se.StellarEvolution.__init__(self,  SeBaInterface(**options), **options)

    def evolve_model(self, end_time=None, keep_synchronous=True):
        self.stopping_conditions.supernova_detection.unset()
        if not keep_synchronous:
            raise Exception("non_synchronous evolution not implemented")
        evolve_a_success = 0
        if end_time is None:
            end_time = self.model_time + min(self.particles.time_step)
#        print "t=", self.model_time, end_time
        if not self.stopping_conditions.supernova_detection.is_enabled():
            evolve_a_success =  self.evolve_system(end_time)        
        else:
            new_end_time = min(end_time, self.model_time + min(self.particles.time_step))
            while self.model_time<new_end_time:
                old_particles = self.particles.copy()
                evolve_a_success =  self.evolve_system(new_end_time)
                psn = self.particles[numpy.logical_and(self.particles.stellar_type >= 13|units.stellar_type, self.particles.stellar_type <= 15|units.stellar_type)]
                psn -= self.particles[numpy.logical_and(old_particles.stellar_type >= 13|units.stellar_type, old_particles.stellar_type <= 15|units.stellar_type)]
    
                if len(psn)>0:
#                    print "Supernova at time:", psn
                    self.stopping_conditions.supernova_detection.set(psn)
                    break
                new_end_time = min(end_time, self.model_time + min(self.particles.time_step))

            """
            to be used if stopping condition is implemented in c
        old_particles = self.particles.copy()
        evolve_a_success =  self.evolve_system(end_time)
        if self.stopping_conditions.supernova_detection.is_set():
            print "XXStellar supernova stopping condition is set"
            print self.particles
            print "age=", self.particles.age, self.particles.relative_age

            self.reset_all_stars(old_particles)
            evolve_a_success =  self.evolve_system(new_end_time)
            """

        return evolve_a_success

    def reset_all_stars(self, old_particles):
        psn = self.particles[self.particles.stellar_type==14|units.stellar_type]
        psn -= self.particles[old_particles.stellar_type==14|units.stellar_type]

        # print "t=", psn.relative_age[0]
        tsn = self.model_time - psn.get_relative_age()[0]
#        channel_from_old_to_new_star = old_particles.new_channel_to(self.particles)
#        channel_from_old_to_new_star.copy_attributes(["relative_age", "relative_mass"])
#        self.evolve_model(tsn)

        self.particles.remove_particles(self.particles)
        self.particles.add_particles(old_particles)
        self.model_time = tsn
        self.evolve_model(tsn)

    def define_properties(self, handler):
        se.StellarEvolution.define_properties(self, handler)
        handler.add_property('get_time', public_name = "model_time")


    def define_methods(self, handler):
        se.StellarEvolution.define_methods(self, handler)

        handler.add_method(
            "evolve_for",
            (handler.INDEX, units.Myr),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "evolve_star",
            (units.MSun, units.Myr, units.none),
            (units.Myr, units.MSun, units.RSun, units.LSun, units.K, units.Myr,units.stellar_type, handler.ERROR_CODE)
        )
        handler.add_method(
            "evolve_system",
            (units.Myr,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "new_binary",
            (units.RSun, handler.NO_UNIT, handler.LINK('particles'), handler.LINK('particles')),
            (handler.INDEX, handler.ERROR_CODE,)
        )
        handler.add_method(
            "delete_binary",
            (handler.INDEX,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_eccentricity",
            (handler.INDEX,),
            (handler.NO_UNIT, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_semi_major_axis",
            (handler.INDEX,),
            (units.RSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_core_mass",
            (handler.INDEX,),
            (units.MSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_COcore_mass",
            (handler.INDEX,),
            (units.MSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "change_mass",
            (handler.INDEX,units.MSun,units.Myr),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "merge_the_binary",
            (handler.INDEX,handler.LINK('particles'),handler.LINK('particles')),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "merge_with_other_star",
            (handler.INDEX,handler.LINK('particles')),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "refresh_memory",
            (handler.INDEX),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "recall_memory_one_step",
            (handler.INDEX),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_envelope_mass",
            (handler.INDEX,),
            (units.MSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_core_radius",
            (handler.INDEX,),
            (units.RSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_semi_major_axis",
            (handler.INDEX, units.RSun,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_age",
            (handler.INDEX,),
            (units.Myr, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_time_step",
            (handler.INDEX,),
            (units.Myr, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_supernova_kick_velocity",
            (),
            (units.kms, handler.ERROR_CODE,)
        )
        handler.add_method(
            "set_supernova_kick_velocity",
            (units.kms,),
            (handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_gyration_radius_sq",
            (handler.INDEX,),
            (units.none, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_relative_age",
            (handler.INDEX,),
            (units.Myr, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_natal_kick_velocity",
            (handler.INDEX,),
            (units.kms, units.kms, units.kms, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_relative_mass",
            (handler.INDEX,),
            (units.MSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_effective_radius",
            (handler.INDEX,),
            (units.RSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_convective_envelope_mass",
            (handler.INDEX,),
            (units.MSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_convective_envelope_radius",
            (handler.INDEX,),
            (units.RSun, handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_time",
            (),
            (units.Myr,handler.ERROR_CODE,)
        )
        handler.add_method(
            "get_wind_mass_loss_rate",
            (handler.INDEX,),
            (units.MSun/units.yr, handler.ERROR_CODE,)
        )
        self.stopping_conditions.define_methods(handler)

    def update_time_steps(self):
        pass

    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_metallicity",
            "set_metallicity",
            "metallicity",
            "Metallicity of all stats",
            default_value = 0.02
        )
        handler.add_method_parameter(
            "get_supernova_kick_velocity",
            "set_supernova_kick_velocity",
            "supernova_kick_velocity",
            "Kick velocity to compact handler formed in supernova",
            default_value = 600 | units.kms
        )

        handler.add_method_parameter(
            "get_is_logging_of_evolve_enabled",
            "set_is_logging_of_evolve_enabled",
            "is_logging_of_evolve_enabled",
            "if True will log star state before and after evolve in starev.data",
            default_value = False
        )

        self.stopping_conditions.define_parameters(handler)


    def define_particle_sets(self, handler):

        handler.define_set('particles', 'index_of_the_star')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_star')

        handler.add_getter('particles', 'get_radius', names = ('radius',))
        handler.add_getter('particles', 'get_stellar_type', names = ('stellar_type',))
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        handler.add_getter('particles', 'get_core_mass', names = ('core_mass',))
        handler.add_getter('particles', 'get_COcore_mass', names = ('CO_core_mass',))
        handler.add_getter('particles', 'get_envelope_mass', names = ('envelope_mass',))
        handler.add_getter('particles', 'get_core_radius', names = ('core_radius',))
        handler.add_getter('particles', 'get_age', names = ('age',))
        handler.add_getter('particles', 'get_time_step', names = ('time_step',))
        #handler.add_getter('particles', 'get_spin', names = ('spin',))
        handler.add_getter('particles', 'get_luminosity', names = ('luminosity',))
        handler.add_getter('particles', 'get_temperature', names = ('temperature',))
        handler.add_getter('particles', 'get_natal_kick_velocity', names = ('natal_kick_x','natal_kick_y','natal_kick_z'))
        handler.add_getter('particles', 'get_convective_envelope_mass', names = ('convective_envelope_mass',))
        handler.add_getter('particles', 'get_convective_envelope_radius', names = ('convective_envelope_radius',))
        handler.add_getter('particles', 'get_gyration_radius_sq', names = ('gyration_radius_sq',))
        handler.add_getter('particles', 'get_relative_age', names = ('relative_age',))
        handler.add_getter('particles', 'get_relative_mass', names = ('relative_mass',))
        handler.add_getter('particles', 'get_wind_mass_loss_rate', names = ('wind_mass_loss_rate',))
        handler.add_getter('particles', 'get_effective_radius', names = ('effective_radius',))

        handler.add_method('particles', 'evolve_one_step')
        handler.add_method('particles', 'evolve_for')
        handler.add_method('particles', 'change_mass')
        handler.add_method('particles', 'refresh_memory')
        handler.add_method('particles', 'recall_memory_one_step')
        handler.add_method('particles', 'merge_with_other_star')

        handler.define_set('binaries', 'index_of_the_star')
        handler.set_new('binaries', 'new_binary')
        handler.set_delete('binaries', 'delete_binary')

        handler.add_getter('binaries', 'get_semi_major_axis', names = ('semi_major_axis',))
        handler.add_getter('binaries', 'get_eccentricity', names = ('eccentricity',))
        handler.add_getter('binaries', 'get_mass', names = ('mass',))
        handler.add_getter('binaries', 'get_time_step', names = ('time_step',))
        handler.add_getter('binaries', 'get_age', names = ('age',))
        handler.add_getter("binaries", 'get_children_of_binary')
        handler.add_setter('binaries', 'set_semi_major_axis', names = ('semi_major_axis',))
        handler.add_setter('binaries', 'set_eccentricity', names = ('eccentricity',))
        handler.add_method('binaries', 'merge_the_binary')
    def define_state(self, handler):
        se.StellarEvolution.define_state(self, handler)
        
        self.stopping_conditions.define_state(handler)



Seba = SeBa
