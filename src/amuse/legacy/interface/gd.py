"""
Stellar Dynamics Interface Defintion
"""

from amuse.legacy.support.core import legacy_function, LegacyFunctionSpecification
from amuse.support.interface import CodeInterface
from amuse.support.units import nbody_system
from amuse.legacy.interface import common

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
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The radius of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The initial position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The initial velocity vector of the particle")
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
        function.addParameter('radius', dtype='float64', direction=function.OUT, description = "The current radius of the particle")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.OUT, description = "The current velocity vector of the particle")
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
        function.addParameter('radius', dtype='float64', direction=function.IN, description = "The new radius of the particle")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('vx', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vy', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
        function.addParameter('vz', dtype='float64', direction=function.IN, description = "The new velocity vector of the particle")
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
        Retrieve the potential at a (particle) position (vector).
        
        *Need better description of use and relation to get_acceleration and get_gravity*
        """
        function = LegacyFunctionSpecification()
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The current position vector of the particle")
        function.addParameter('V', dtype='float64', direction=function.OUT, description = "The current scalar potential...")
        function.result_type = 'int32'
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
    def evolve():
        """
        Evolve the model until the given time or until a collision happens.
        """
        function = LegacyFunctionSpecification()
        function.addParameter('time',
            dtype='float64',
            direction=function.IN,
            description =
                "Model time to evolve the code to. "
                "The model will be evolved until "
                "this time is reached exactly or just before")
        """
        cello, probably gonna kill this one...........
        function.addParameter('collision_flag',
            dtype='float64',
            direction=function.IN,
            description =
                "(1) Stop evolving the model when a collision is detected by the code "
                "(0) Continue evolving the code, ignore collisions")
        """
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Evolved until time, no collision happened
        1 - COLLISION DETECTED
            Stopped after a collision
        -1 - ERROR
            Model did not converge
        -2 - ERROR
            Procedure get_kinetic_energy FAILED
        """
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
    def get_indices_of_colliding_particles():
        """
        Retrieve the two indices of the colliding particles
        """
        function = LegacyFunctionSpecification()
        function.addParameter('index_of_particle1', dtype='int32', direction=function.OUT,
            description = "Index of the first colliding partner")
        function.addParameter('index_of_particle2', dtype='int32', direction=function.OUT,
            description = "Index of the second colliding partner")
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            The indices of the collision partners were set
         -1 - ERROR
            No collision detected during evolve
         -2 - ERROR
            not yet implemented
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
    def get_gravity_at_point():
        """
        Determine the gravitational force on a given point
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eps', dtype='float64', direction=function.IN,
            description = "The smoothing parameter")
        function.addParameter('x', dtype='float64', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('y', dtype='float64', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('z', dtype='float64', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('forcex', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.addParameter('forcey', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.addParameter('forcez', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.result_type = 'int32'
        function.can_handle_array = True
        function.result_doc = """
         0 - OK
            Force could be calculated
        -1 - ERROR
            No force calculation supported
        """
        return function


    @legacy_function
    def get_potential_at_point():
        """
        Determine the potential on a given point
        """
        function = LegacyFunctionSpecification()
        function.addParameter('eps', dtype='float64', direction=function.IN,
         description = "The smoothing factor, may be ignored by the code")
        function.addParameter('x', dtype='float64', direction=function.IN)
        function.addParameter('y', dtype='float64', direction=function.IN)
        function.addParameter('z', dtype='float64', direction=function.IN)
        function.addParameter('phi', dtype='float64', direction=function.OUT)
        function.can_handle_array = True
        function.result_type = 'int32'
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

class GdAutoDoc(object):

    def __get__(self, instance, owner):

        string = ""

        string = instance.parameters.__doc__
        return string

class GravitationalDynamics(common.CommonCode):
    NBODY = object()

    __doc__ = GdAutoDoc()

    def __init__(self, legacy_interface, convert_nbody = None,  **options):
        self.convert_nbody = convert_nbody

        CodeInterface.__init__(self, legacy_interface, **options)

    def define_properties(self, object):
        object.add_property("get_kinetic_energy", nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2)
        object.add_property("get_potential_energy", nbody_system.mass * nbody_system.length ** 2  * nbody_system.time ** -2)
        object.add_property("get_total_radius", nbody_system.length)
        object.add_property("get_center_of_mass_position", nbody_system.length)
        object.add_property("get_center_of_mass_velocity", nbody_system.length / nbody_system.time)
        object.add_property("get_total_mass", nbody_system.mass)
        object.add_property('get_time', nbody_system.time, "model_time")

    def define_state(self, object):
        common.CommonCode.define_state(self, object)

        object.add_transition('INITIALIZED','EDIT','commit_parameters')
        object.add_method('EDIT', 'new_particle')
        object.add_method('EDIT', 'delete_particle')
        object.add_transition('EDIT', 'RUN', 'commit_particles')
        object.add_transition('RUN', 'UPDATE', 'new_particle', False)
        object.add_transition('RUN', 'UPDATE', 'delete_particle', False)
        object.add_transition('UPDATE', 'RUN', 'recommit_particles')
        object.add_transition('RUN', 'EVOLVED', 'evolve_model', False)
        object.add_method('EVOLVED', 'evolve_model')
        object.add_transition('EVOLVED','RUN', 'synchronize_model')
        object.add_method('RUN', 'synchronize_model')
        object.add_method('RUN', 'get_state')
        object.add_method('RUN', 'get_mass')
        object.add_method('RUN', 'get_position')
        object.add_method('RUN', 'get_gravity_at_point')
        object.add_method('RUN', 'get_potential_at_point')

    def define_parameters(self, object):
        object.add_method_parameter(
            "get_time_step",
            None,
            "timestep",
            "constant timestep for iteration",
            nbody_system.time,
            0.7 | nbody_system.time
            )

    def define_methods(self, object):

        object.add_method(
            'evolve',
            (nbody_system.time,),
            public_name = 'evolve_model'
        )

        object.add_method(
            "new_particle",
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ),
            (
                object.INDEX,
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "delete_particle",
            (
                object.NO_UNIT,
            ),
            (
                object.ERROR_CODE,
            )
        )
        object.add_method(
            "get_state",
            (
                object.NO_UNIT,
            ),
            (
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_state",
            (
                object.NO_UNIT,
                nbody_system.mass,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                nbody_system.speed,
                nbody_system.speed,
                nbody_system.speed,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_mass",
            (
                object.NO_UNIT,
                nbody_system.mass,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_mass",
            (
                object.NO_UNIT,
            ),
            (
                nbody_system.mass,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_radius",
            (
                object.NO_UNIT,
                nbody_system.length,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_radius",
            (
                object.NO_UNIT,
            ),
            (
                nbody_system.length,
                object.ERROR_CODE
            )
        )
        object.add_method(
            "set_position",
            (
                object.NO_UNIT,
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
            ),
            (
                object.ERROR_CODE
            )
        )
        object.add_method(
            "get_position",
            (
                object.NO_UNIT,
            ),
            (
                nbody_system.length,
                nbody_system.length,
                nbody_system.length,
                object.ERROR_CODE
            )
        )

        #object.add_method(
        #    "get_time_step",
        #    (
        #        object.NO_UNIT,
        #    ),
        #    (
        #        nbody_system.time,
        #        object.ERROR_CODE
        #    )
        #)
        #
        object.add_method(
            'get_indices_of_colliding_particles',
            (),
            (
                object.NO_UNIT,
                object.NO_UNIT,
                object.ERROR_CODE,
            )
        )

        object.add_method(
            'get_center_of_mass_position',
            (),
            (nbody_system.length, nbody_system.length, nbody_system.length, object.ERROR_CODE)
        )

        object.add_method(
            'get_center_of_mass_velocity',
            (),
            (nbody_system.speed, nbody_system.speed, nbody_system.speed, object.ERROR_CODE)
        )

        object.add_method(
            'get_gravity_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.acceleration, nbody_system.acceleration, nbody_system.acceleration, object.ERROR_CODE)
        )

        object.add_method(
            'get_potential_at_point',
            (nbody_system.length, nbody_system.length, nbody_system.length, nbody_system.length),
            (nbody_system.potential, object.ERROR_CODE)
        )


    def define_particle_sets(self, object):
        object.define_set('particles', 'index_of_the_particle')
        object.set_new('particles', 'new_particle')
        object.set_delete('particles', 'delete_particle')
        object.add_setter('particles', 'set_state')
        object.add_getter('particles', 'get_state')
        object.add_setter('particles', 'set_mass')
        object.add_getter('particles', 'get_mass', names = ('mass',))
        object.add_setter('particles', 'set_position')
        object.add_getter('particles', 'get_position')
        object.add_query('particles', 'get_indices_of_colliding_particles', public_name = 'select_colliding_particles')

    def get_colliding_particles(self):
        subset = self.colliding_particles_method._run(self, self.particles)
        return subset

    def define_converter(self, object):
        if not self.convert_nbody is self.NBODY:
            object.set_nbody_converter(self.convert_nbody)
