"""
Stellar Dynamics Interface Defintion
"""

from amuse.legacy.support.core import legacy_function, RemoteFunction

class GravitationalDynamics(object):

    @legacy_function
    def new_particle():
        """
        Define a new particle in the stellar dynamics code. The particle is initialized with the provided
        mass, radius, position and velocity. This function returns an index that can be used to refer
        to this particle.
        """
        function = RemoteFunction()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.OUT
            , description = 
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
        function = RemoteFunction()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to be removed. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was removed from the model
        -1 - ERROR
            particle could not be removed
        """
        return function

    @legacy_function
    def get_state():
        """
        Retrieve the current state of a particle. The *minimal* information of a stellar 
        dynamics particle (mass, radius, position and velocity) is returned.
        """
        function = RemoteFunction()  
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
        function = RemoteFunction()  
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
        """
        return function    
        

    @legacy_function
    def get_mass():
        """
        Retrieve the mass of a particle. Mass is a scalar property of a particle,
        this function has one OUT argument.
        """
        function = RemoteFunction()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.OUT, description = "The current mass of the particle")
        function.result_type = 'int32'
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
        Update the mass of a particle. Mass is a scalar property of a particle
        """
        function = RemoteFunction()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('mass', dtype='float64', direction=function.IN, description = "The new mass of the particle")
        function.result_type = 'int32'
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
    def get_position():
        """
        Retrieve the position vector of a particle. Position is a vector property,
        this function has 3 OUT arguments.
        """
        function = RemoteFunction()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function

    @legacy_function
    def set_position():
        """
        Update the position of a particle.
        """
        function = RemoteFunction()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new position vector of the particle")
        function.result_type = 'int32'
        function.result_type = 'int32'
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
    def get_acceleration():
        """
        Retrieve the acceleration vector of a particle. Second time derivative of the position.
        """
        function = RemoteFunction()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current position vector of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
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
        function = RemoteFunction()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle for which the state is to be updated. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.IN, description = "The new acceleration vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.IN, description = "The new acceleration vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.IN, description = "The new acceleration vector of the particle")
        function.result_type = 'int32'
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            particle was found in the model and the information was set
        -1 - ERROR
            particle could not be found
        -2 - ERROR
            code does not support updating of a particle 
        """
    
    @legacy_function
    def get_potential():
        """
        Retrieve the potential vector of a particle. 
        *Need better description of use and relation to get_acceleration and get_gravity*
        """
        function = RemoteFunction()  
        function.addParameter('index_of_the_particle', dtype='int32', direction=function.IN,
            description = "Index of the particle to get the state from. This index must have been returned by an earlier call to :meth:`new_particle`")
        function.addParameter('x', dtype='float64', direction=function.OUT, description = "The current potential vector of the particle")
        function.addParameter('y', dtype='float64', direction=function.OUT, description = "The current potential vector of the particle")
        function.addParameter('z', dtype='float64', direction=function.OUT, description = "The current potential vector of the particle")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            current value was retrieved
        -1 - ERROR
            particle could not be found
        """
        return function
    
        
    @legacy_function
    def evolve():
        """
        Evolve the model until the given time or until a collision happens.
        """
        function = RemoteFunction()  
        function.addParameter('time', dtype='float64', direction=function.IN,
            description = "Model time to evolve the code to. The model will be evolved until this time is reached exactly or just before")
        function.addParameter('collision_flag', dtype='float64', direction=function.IN, 
            description = "(1) Stop evolving the model when a collision is detected by the code (0) Continue evolving the code, ignore collisions")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Evolved until time, no collision happened
        1 - COLLISION DETECTED
            Stopped after a collision
        -1 - ERROR
            Model did not converge 
        """
        return function  
    
    @legacy_function
    def initialize_particles():
        """
        Let the code perform initialization actions after all particles have loaded
        in the model. Called before the first evolve call and after the last new_particle call.
        """
        function = RemoteFunction()  
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Model is initialized and evolution can start
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention 
        """
        return function  
    
    @legacy_function
    def initialize_code():
        """
        Let the code perform initialization actions after all parameters have been set.
        """
        function = RemoteFunction()  
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Code is initialized
         -1 - ERROR
            Error happened during initialization, this error needs to be further specified by every code implemention 
        """
        return function  
        
    
    @legacy_function
    def get_eps2():
        """
        Retrieve the current value of the squared smoothing parameter.
        """
        function = RemoteFunction()  
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
        function = RemoteFunction()  
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
        function = RemoteFunction()  
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
        function = RemoteFunction()  
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
        function = RemoteFunction()  
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
        """
        return function  
    
    @legacy_function
    def get_time():
        """
        Retrieve the model time. This time should be close to the end time specified
        in the evolve code. Or, when a collision was detected, it will be the
        model time of the collision.
        """
        function = RemoteFunction()  
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
    def get_total_mass():
        """
        Retrieve the sum of the masses of all particles.
        """
        function = RemoteFunction()  
        function.addParameter('mass', dtype='float64', direction=function.OUT,
            description = "The total mass of the model")
        function.result_type = 'int32'
        function.result_doc = """
        0 - OK
            Current value of the kinetic mass was retrieved 
        -1 - ERROR
            Total mass could not be provided
        """
        return function
    
    @legacy_function
    def get_center_of_mass_position():
        """
        Retrieve the center of mass (a point in space) of all particles.
        """
        function = RemoteFunction()  
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
        """
        return function
    
    @legacy_function
    def get_center_of_mass_velocity():
        """
        Retrieve the velocity of the center of mass of all particles. This 
        velocity is mass weighted mean of the velocity of all particles.
        """
        function = RemoteFunction()  
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
        Return the radius of the sphere centered on the center of mass that
        contains all the particles. *get_size?*
        """
        function = RemoteFunction()  
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
        function = RemoteFunction()  
        function.addParameter('x', dtype='int32', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('y', dtype='int32', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('z', dtype='int32', direction=function.IN,
            description = "The position vector of the point")
        function.addParameter('forcex', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.addParameter('forcey', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.addParameter('forcez', dtype='float64', direction=function.OUT,
            description = "Force created by the particles in the code at the given position")
        function.result_type = 'int32'
        function.result_doc = """
         0 - OK
            Force could be calculated
        -1 - ERROR
            No force calculation supported
        """
        return function  
    


    @legacy_function
    def get_number_of_particles():
        """
        Retrieve the total number of particles defined in the code
        """
        function = RemoteFunction()  
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
        function = RemoteFunction()  
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
        function = RemoteFunction()  
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

    