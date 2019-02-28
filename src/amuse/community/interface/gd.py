"""
Stellar Dynamics Interface Defintion
"""

from amuse.support.interface import InCodeComponentImplementation
from amuse.units import nbody_system
from amuse.units import generic_unit_converter
from amuse.community.interface import common

from amuse.rfi.core import legacy_function
from amuse.rfi.core import LegacyFunctionSpecification


class GravitationalDynamicsInterface(common.CommonCodeInterface):

    @legacy_function
    def new_particle():
        """
        Define a new particle in the stellar dynamics code. The particle is initialized with the provided
        mass, radius, position and velocity. This function returns an index that can be used to refer
        to this particle.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT, description =
            """
            An index assigned to the newly created particle.
            This index is supposed to be a local index for the code
            (and not valid in other instances of the code or in other codes)
            """
            )

        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of the particle", default = 0)
        function.result_type = 'int32'
        function.result_doc = """ 0 - OK
            particle was created and added to the model
        -1 - ERROR
            particle could not be created"""
        return function

    @legacy_function
    def delete_particle():
        """
        Remove the definition of particle from the code. After calling this function the particle is
        no longer part of the model evolution. It is up to the code if the index will be reused.
        This function is optional.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to be removed. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be removed
        -2 - ERROR
            not yet implemented
        """
        return function
        
    @legacy_function
    def get_state():
        """
        Retrieve the current state of a particle. The *minimal* information of a stellar
        dynamics particle (mass, radius, position and velocity) is returned.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The current radius of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_state():
        """
        Update the current state of a particle. The *minimal* information of a stellar
        dynamics particle (mass, radius, position and velocity) is updated.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The new radius of the particle", default = 0)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def get_mass():
        """
        Retrieve the mass of a particle. Mass is a scalar property of a particle,
        this function has one OUT argument.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_mass():
        """
        Update the mass of a particle. Mass is a scalar property of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        """
        return function

    @legacy_function
    def get_radius():
        """
        Retrieve the radius of a particle. Radius is a scalar property of a particle,
        this function has one OUT argument.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the radius of. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The current radius of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was retreived
        -1 - ERROR
            particle could not be found
        """
        return function


    @legacy_function
    def set_radius():
        """
        Set the radius of a particle. Radius is a scalar property of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the radius of. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The new radius of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was retreived
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def get_position():
        """
        Retrieve the position vector of a particle. Position is a vector property,
        this function has 3 OUT arguments.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def set_position():
        """
        Update the position of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        """
        return function

    @legacy_function
    def get_velocity():
        """
        Retrieve the velocity vector of a particle. Position is a vector property,
        this function has 3 OUT arguments.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the velocity from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current x component of the position vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current y component of the position vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current z component of the position vector of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function


    @legacy_function
    def set_velocity():
        """
        Set the velocity vector of a particle.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The current x component of the velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The current y component of the velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The current z component of the velocity vector of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function


    @legacy_function
    def get_acceleration():
        """
        Retrieve the acceleration vector of a particle. Second time derivative of the position.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('ax', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('ay', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('az', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def set_acceleration():
        """
        Update the acceleration of a particle.
        *Defined for symetry with the get_acceleration function.*
        *Should be removed if physaccily unsound*
        *Maybe moved to snapshot support functionality*
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('ax', dtype='float64', direction=function.IN, description = "The new acceleration vector of the particle")
        function.addParameter('ay', dtype='float64', direction=function.IN, description = "The new acceleration vector of the particle")
        function.addParameter('az', dtype='float64', direction=function.IN, description = "The new acceleration vector of the particle")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle
        -3 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def get_potential():
        """
        Retrieve the potential at a particle position, for retrieving the potential anywhere in
        the field use get_potential_at_point.
        """
        function = LegacyFunctionSpecification()
        
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('potential', dtype='float64', direction=function.OUT, description = "The current scalar potential...")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def evolve_model():
        """
        Evolve the model until the given time, or until a stopping condition is set.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN,
            description = "Model time to evolve the code to. The model will be "
                "evolved until this time is reached exactly or just after.")
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
    def synchronize_model():
        """
        After an evolve the particles may be at different simulation
        times. Synchronize the particles to a consistent stat
        at the current simulation time
        """
        function = LegacyFunctionSpecification()
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK

        """
        return function

    @legacy_function
    def recommit_particles():
        """
        Let the code perform initialization actions
        after the number of particles have been updated
        or particle attributes have been updated from
        the script.
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
    def get_eps2():
        """
        Retrieve the current value of the squared smoothing parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared', dtype='float64', direction=function.OUT,
            description = "The current value of the smooting parameter, squared.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the smoothing parameter was set
        -1 - ERROR
            The code does not have support for a smoothing parameter
        """
        return function


    @legacy_function
    def set_eps2():
        """
        Update the value of the squared smoothing parameter.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('epsilon_squared', dtype='float64', direction=function.IN,
            description = "The new value of the smooting parameter, squared.")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the smoothing parameter was set
        -1 - ERROR
            The code does not have support for a smoothing parameter
        """
        return function



    @legacy_function
    def get_kinetic_energy():
        """
        Retrieve the current kinetic energy of the model
        """
        function = LegacyFunctionSpecification()
        function.addParameter('kinetic_energy', dtype='float64', direction=function.OUT,
            description = "The kinetic energy of the model")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the kinetic energy was set
        -1 - ERROR
            Kinetic energy could not be provided
        """
        return function

    @legacy_function
    def get_potential_energy():
        """
        Retrieve the current potential energy of the model
        """
        function = LegacyFunctionSpecification()
        function.addParameter('potential_energy', dtype='float64', direction=function.OUT,
            description = "The potential energy of the model")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the potential energy was set
        -1 - ERROR
            Kinetic potential could not be provided
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
    def get_begin_time():
        """
        Retrieve the model time to start the evolution at.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.OUT,
            description = "The begin time", unit = nbody_system.time)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the time was retrieved
        -2 - ERROR
            The code does not have support for querying the begin time
        """
        return function
    
    @legacy_function
    def set_begin_time():
        """
        Set the model time to start the evolution at. This is an offset for
        all further calculations in the code.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time', dtype='float64', direction=function.IN,
            description = "The model time to start at", unit = nbody_system.time)
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Time value was changed
        -2 - ERROR
            The code does not support setting the begin time
        """
        return function
        
    @legacy_function
    def get_time_step():
        """
        Retrieve the model timestep.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time_step', dtype='float64', direction=function.OUT,
            description = "The current model timestep")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the time step was retrieved
        -1 - ERROR
            The code does not have support for querying the time
        """
        return function

    @legacy_function
    def get_total_mass():
        """
        Retrieve the sum of the masses of all particles.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('mass', dtype='float64', direction=function.OUT,
            description = "The total mass of the model")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the kinetic mass was retrieved
        -1 - ERROR
            Total mass could not be provided
        -2 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def get_center_of_mass_position():
        """
        Retrieve the center of mass (a point in space) of all particles.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('x', dtype='float64', direction=function.OUT,
            description = "The center of mass of the model")
        function.addParameter('y', dtype='float64', direction=function.OUT,
            description = "The center of mass of the model")
        function.addParameter('z', dtype='float64', direction=function.OUT,
            description = "The center of mass of the model")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the center was retrieved
        -1 - ERROR
            The mass could not be provided
        -2 - ERROR
            not yet implemented
        """
        return function

    @legacy_function
    def get_center_of_mass_velocity():
        """
        Retrieve the velocity of the center of mass of all particles. This
        velocity is mass weighted mean of the velocity of all particles.
        """
        function = LegacyFunctionSpecification()
        function.can_handle_array = True
        function.addParameter('vx', dtype='float64', direction=function.OUT,
            description = "The mean velocity of the model")
        function.addParameter('vy', dtype='float64', direction=function.OUT,
            description = "The mean velocity of the model")
        function.addParameter('vz', dtype='float64', direction=function.OUT,
            description = "The mean velocity  of the model")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the center of mass velocity was retrieved
        -1 - ERROR
            The value could not be provided
        """
        return function


    @legacy_function
    def get_total_radius():
        """
        Return the radius of the sphere, centered on the center of mass that
        contains all the particles. *get_size?*
        """
        function = LegacyFunctionSpecification()
        function.addParameter('radius', dtype='float64', direction=function.OUT,
            description = "The maximum distance from a star to the center of mass of the model")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the radius was retrieved
        -1 - ERROR
            The value could not be provided
        """
        return function

   
    @legacy_function
    def get_number_of_particles():
        """
        Retrieve the total number of particles define  d in the code
        """
        function = LegacyFunctionSpecification()
        function.addParameter('number_of_particles', dtype='int32', direction=function.OUT,
            description = "Count of the particles in the code")
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            Count could be determined
         -1 - ERROR
            Unable to determine the count
        """
        return function

    @legacy_function
    def get_index_of_first_particle():
        """
        Retrieve the index of first particle. When this index is used as the
        starting index for the ``get_index_of_next_particle`` method, all
        particles can be iterated over::

            error, first_index = instance.get_index_of_first_particle()
            current_index = first_index
            while error == 0:
                status, mass = instance.get_mass(current_index)
                print mass
                error, current_index = instance.get_index_of_next_particle(current_index)
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT,
            description = "Index of the first particle")
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            Index was set
         1 - ERROR
            Code has no particles, or cannot set the index
        """
        return function

    @legacy_function
    def get_index_of_next_particle():
        """
        Retrieve the index of the particle following the provided index. The
        iteration direction is determined by the code.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle")
        function.addParameter('index_of_the_next_particle', dtype='int32', direction=function.OUT,
            description = "Index of the particle following the particle with the provided index")
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            Index was set
         1 - STATE
            Last index was reached
         -1 - ERROR
            Particle could not be found
        """
        return function
    
        
class GravityFieldInterface(object):
    """
    Codes implementing the gravity field interface provide functions to
    calculate the force and potential energy fields at any point.
    """
    
    @legacy_function    
    def get_gravity_at_point():
        """
        Get the gravitational acceleration at the given points. To calculate the force on
        bodies at those points, multiply with the mass of the bodies
        """
        function = LegacyFunctionSpecification()  
        for x in ['eps','x','y','z']:
            function.addParameter(
              x, 
              dtype='float64', 
              direction=function.IN,
              unit=nbody_system.length
            )
        for x in ['ax','ay','az']:
            function.addParameter(
                x, 
                dtype='float64', 
                direction=function.OUT,     
                unit=nbody_system.acceleration
            )
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32' 
        function.must_handle_array = True
        return function
        
    @legacy_function    
    def get_potential_at_point():
        """
        Determine the gravitational potential on any given point
        """
        function = LegacyFunctionSpecification()  
        for x in ['eps','x','y','z']:
            function.addParameter(
                x, 
                dtype='float64', 
                direction=function.IN,
                unit=nbody_system.length
            )
        for x in ['phi']:
            function.addParameter(
                x, 
                dtype='float64', 
                direction=function.OUT,
                unit=nbody_system.potential
            )
        function.addParameter('npoints', dtype='i', direction=function.LENGTH)
        function.result_type = 'int32'
        function.must_handle_array = True
        return function


class SinglePointGravityFieldInterface(object):
    """
    Codes implementing the gravity field interface provide functions to
    calculate the force and potential energy fields at any point.
    """
    
    @legacy_function    
    def get_gravity_at_point():
        """
        Get the gravitational acceleration at the given points. To calculate the force on
        bodies at those points, multiply with the mass of the bodies
        """
        function = LegacyFunctionSpecification()  
        for x in ['eps','x','y','z']:
            function.addParameter(
              x, 
              dtype='float64', 
              direction=function.IN,
              unit=nbody_system.length
            )
        for x in ['ax','ay','az']:
            function.addParameter(
                x, 
                dtype='float64', 
                direction=function.OUT,     
                unit=nbody_system.acceleration
            )
        function.result_type = 'int32' 
        function.can_handle_array = True
        return function
        
    @legacy_function    
    def get_potential_at_point():
        """
        Determine the gravitational potential on any given point
        """
        function = LegacyFunctionSpecification()  
        for x in ['eps','x','y','z']:
            function.addParameter(
                x, 
                dtype='float64', 
                direction=function.IN,
                unit=nbody_system.length
            )
        for x in ['phi']:
            function.addParameter(
                x, 
                dtype='float64', 
                direction=function.OUT,
                unit=nbody_system.potential
            )
        function.result_type = 'int32'
        function.can_handle_array = True
        return function


class GravitationalDynamicsDocumentation(object):

    def __get__(self, instance, owner):

        string = ""

        string = instance.parameters.__doc__
        return string

class GravitationalDynamics(common.CommonCode):
    NBODY = object()

    __doc__ = GravitationalDynamicsDocumentation()

    def __init__(self, legacy_interface, unit_converter = None,  **options):
        self.unit_converter = unit_converter
        
        common.CommonCode.__init__(self, legacy_interface, **options)

    def define_properties(self, handler):
        handler.add_property("get_kinetic_energy")
        handler.add_property("get_potential_energy")
        handler.add_property("get_total_radius")
        handler.add_property("get_center_of_mass_position")
        handler.add_property("get_center_of_mass_velocity")
        handler.add_property("get_total_mass")
        handler.add_property('get_time', public_name = "model_time")

    def define_state(self, handler): 
        common.CommonCode.define_state(self, handler)   
        handler.add_transition('END', 'INITIALIZED', 'initialize_code', False)    
        
        handler.add_transition('INITIALIZED','EDIT','commit_parameters')
        handler.add_transition('RUN','CHANGE_PARAMETERS_RUN','before_set_parameter', False)
        handler.add_transition('EDIT','CHANGE_PARAMETERS_EDIT','before_set_parameter', False)
        handler.add_transition('UPDATE','CHANGE_PARAMETERS_UPDATE','before_set_parameter', False)
        handler.add_transition('CHANGE_PARAMETERS_RUN','RUN','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_EDIT','EDIT','recommit_parameters')
        handler.add_transition('CHANGE_PARAMETERS_UPDATE','UPDATE','recommit_parameters')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_set_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_set_parameter')
        
        handler.add_method('CHANGE_PARAMETERS_RUN', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_EDIT', 'before_get_parameter')
        handler.add_method('CHANGE_PARAMETERS_UPDATE','before_get_parameter')
        handler.add_method('RUN', 'before_get_parameter')
        handler.add_method('EDIT', 'before_get_parameter')
        handler.add_method('UPDATE','before_get_parameter')
        handler.add_method('EVOLVED','before_get_parameter')
        
        
        handler.add_method('EDIT', 'new_particle')
        handler.add_method('EDIT', 'delete_particle')
        handler.add_method('UPDATE', 'new_particle')
        handler.add_method('UPDATE', 'delete_particle')
        handler.add_transition('EDIT', 'RUN', 'commit_particles')
        handler.add_transition('RUN', 'UPDATE', 'new_particle', False)
        handler.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        handler.add_transition('UPDATE', 'RUN', 'recommit_particles')
        handler.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        handler.add_method('EVOLVED', 'evolve_model')
        handler.add_transition('EVOLVED','RUN', 'synchronize_model')
        handler.add_method('RUN', 'synchronize_model')
        handler.add_method('RUN', 'get_state')
        handler.add_method('RUN', 'get_mass')
        handler.add_method('RUN', 'get_position')
        handler.add_method('RUN', 'get_velocity')
        handler.add_method('RUN', 'get_potential')
        handler.add_method('RUN', 'get_potential_energy')
        handler.add_method('RUN', 'get_kinetic_energy')
        
        
    def define_parameters(self, handler):
        handler.add_method_parameter(
            "get_time_step",
            None,
            "timestep",
            "constant timestep for iteration",
            default_value = 0.7 | nbody_system.time
        )
        
        handler.add_method_parameter(
            "get_begin_time",
            "set_begin_time",
            "begin_time",
            "model time to start the simulation at",
            default_value = 0.0 | nbody_system.time
        )

    def define_methods(self, handler):
        common.CommonCode.define_methods(self, handler)
        handler.add_method(
            'evolve_model',
            (
                nbody_system.time,
            ),
            (
                handler.ERROR_CODE,
            )
        )

        handler.add_method(
            "new_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
            ),
            (
                handler.INDEX,
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "delete_particle",
            (
                handler.NO_UNIT,
            ),
            (
                handler.ERROR_CODE,
            )
        )
        handler.add_method(
            "get_state",
            (
                handler.NO_UNIT,
            ),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_state",
            (
                handler.NO_UNIT,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_mass",
            (
                handler.NO_UNIT,
                nbody_system.mass,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_mass",
            (
                handler.NO_UNIT,
            ),
            (
                nbody_system.mass,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_radius",
            (
                handler.NO_UNIT,
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_radius",
            (
                handler.NO_UNIT,
            ),
            (
                nbody_system.length,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_position",
            (
                handler.NO_UNIT,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_position",
            (
                handler.NO_UNIT,
            ),
            (
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "set_velocity",
            (
                handler.NO_UNIT,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ),
            (
                handler.ERROR_CODE
            )
        )
        handler.add_method(
            "get_velocity",
            (
                handler.NO_UNIT,
            ),
            (
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                handler.ERROR_CODE
            )
        )
        
        handler.add_method(
            "get_potential",
            (
                handler.NO_UNIT,
            ),
            (nbody_system.length ** 2  * nbody_system.time ** -2, handler.ERROR_CODE,)
        )


        handler.add_method(
            'get_indices_of_colliding_particles',
            (),
            (
                handler.NO_UNIT,
                handler.NO_UNIT,
                handler.ERROR_CODE,
            )
        )

       
        handler.add_method(
            'commit_particles',
            (),
            (handler.ERROR_CODE)
        )
        
        handler.add_method(
            'recommit_particles',
            (),
            (handler.ERROR_CODE)
        )
        
        handler.add_method(
            'synchronize_model',
            (),
            (handler.ERROR_CODE)
        )


        handler.add_method(
            "get_time_step",
            (),
            (nbody_system.time, handler.ERROR_CODE,)
        )


        handler.add_method(
            "get_kinetic_energy",
            (),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, handler.ERROR_CODE,)
        )


        handler.add_method(
            "get_potential_energy",
            (),
            (nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2, handler.ERROR_CODE,)
        )


        handler.add_method(
            "get_total_radius",
            (),
            (nbody_system.length, handler.ERROR_CODE,)
        )


        handler.add_method(
            "get_center_of_mass_position",
            (),
            (nbody_system.length,nbody_system.length,nbody_system.length, handler.ERROR_CODE,)
        )


        handler.add_method(
            "get_center_of_mass_velocity",
            (),
            (nbody_system.length / nbody_system.time,nbody_system.length / nbody_system.time,nbody_system.length / nbody_system.time, handler.ERROR_CODE,)
        )


        handler.add_method(
            "get_total_mass",
            (),
            (nbody_system.mass, handler.ERROR_CODE,)
        )


        handler.add_method(
            'get_time',
            (),
            (nbody_system.time, handler.ERROR_CODE,)
        )


    def define_particle_sets(self, handler):
        handler.define_set('particles', 'index_of_the_particle')
        handler.set_new('particles', 'new_particle')
        handler.set_delete('particles', 'delete_particle')
        handler.add_setter('particles', 'set_state')
        handler.add_getter('particles', 'get_state')
        handler.add_setter('particles', 'set_mass')
        handler.add_getter('particles', 'get_mass', names = ('mass',))
        handler.add_setter('particles', 'set_position')
        handler.add_getter('particles', 'get_position')
        handler.add_setter('particles', 'set_velocity')
        handler.add_getter('particles', 'get_velocity')
        handler.add_setter('particles', 'set_radius')
        handler.add_getter('particles', 'get_radius')
        handler.add_query('particles', 'get_indices_of_colliding_particles', public_name = 'select_colliding_particles')

    def get_colliding_particles(self):
        subset = self.colliding_particles_method._run(self, self.particles)
        return subset

    def define_converter(self, handler):
        if not self.unit_converter is None:
            handler.set_converter(self.unit_converter.as_converter_from_si_to_generic())
            
    def commit_parameters(self):
        self.parameters.send_not_set_parameters_to_code()
        self.parameters.send_cached_parameters_to_code()
        self.overridden().commit_parameters()
        
    def cleanup_code(self):
        self.overridden().cleanup_code()
    
        handler = self.get_handler('PARTICLES')
        handler._cleanup_instances()

    def reset(self):
        parameters = self.parameters.copy()
        self.cleanup_code()
        self.initialize_code()
        self.parameters.reset_from_memento(parameters)
        
    def get_total_energy(self):
        return self.get_potential_energy() + self.get_kinetic_energy()
        

class GravityFieldCode(object):
        
    def define_state(self, handler): 
        handler.add_method('RUN', 'get_gravity_at_point')
        handler.add_method('RUN', 'get_potential_at_point')
        
        
