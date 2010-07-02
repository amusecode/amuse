from amuse.support.units import nbody_system
from amuse.support.data import base
from amuse.support.data.particles import ParticlesWithUnitsConverted, AbstractParticleSet
from amuse.support.data.values import zero
from amuse.support.data import values
from amuse.support.data.values import Quantity, new_quantity, zero
from amuse.support.units import constants
from amuse.support.units import units

def move_to_center(particles):
    """
    Move the particle positions to the center of mass, also
    moves the particle velocity to the center of mass velocity.
    """
    particles.position -= particles.center_of_mass()
    particles.velocity -= particles.center_of_mass_velocity()
    
def scale_to_standard(particles, convert_nbody = None, smoothing_length_squared = zero):
    """
    Scale the particles to a standard NBODY model with
    **total mass=1**, **kinetic energy=0.25** and 
    **potential_energy=0.5** (or **viridial_radius=1.0**)
    
    :argument convert_nbody: the scaling is in nbody units and
        when the particles are in si units a convert_nbody is needed
    :argument smoothing_length_squared: needed for calculating
        the potential energy
    """
    if not convert_nbody is None:
        particles = ParticlesWithUnitsConverted(particles, convert_nbody.as_converter_from_nbody_to_si())
        if not smoothing_length_squared is zero:
            smoothing_length_squared = convert_nbody.to_nbody(smoothing_length_squared)
    
    total_mass = particles.mass.sum()
    particles.mass *= ((1 | total_mass.unit) / total_mass)
    
    kinetic_energy = particles.kinetic_energy()
    target_energy =  0.25 | nbody_system.energy
    scale_factor = (target_energy / kinetic_energy).sqrt()
    
    particles.velocity *= scale_factor
    
    potential_energy = particles.potential_energy(G=nbody_system.G)
    target_energy = -0.5 | nbody_system.energy
    scale_factor = (potential_energy / target_energy)
    
    particles.position *= scale_factor 
    


def center_of_mass(particles):
    """
    Returns the center of mass of the particles set.
    The center of mass is defined as the average 
    of the positions of the particles, weighted by their masses.
    
    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.x = [-1.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.center_of_mass()
    quantity<[0.0, 0.0, 0.0] m>
    """

    masses = particles.mass
    x_values = particles.x
    y_values = particles.y
    z_values = particles.z
    
    total_mass = masses.sum()
    massx = (masses * x_values).sum()
    massy = (masses * y_values).sum()
    massz = (masses * z_values).sum()

    return values.VectorQuantity.new_from_scalar_quantities(
        massx/total_mass,
        massy/total_mass,
        massz/total_mass
    )

def center_of_mass_velocity(particles):
    """
    Returns the center of mass velocity of the particles set.
    The center of mass velocity is defined as the average 
    of the velocities of the particles, weighted by their masses.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.vx = [-1.0, 1.0] | units.ms
    >>> particles.vy = [0.0, 0.0] | units.ms
    >>> particles.vz = [0.0, 0.0] | units.ms
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.center_of_mass_velocity()
    quantity<[0.0, 0.0, 0.0] m * s**-1>
    """


    masses = particles.mass
    x_values = particles.vx
    y_values = particles.vy
    z_values = particles.vz
    
    total_mass = masses.sum()
    massx = (masses * x_values).sum()
    massy = (masses * y_values).sum()
    massz = (masses * z_values).sum()

    return values.VectorQuantity.new_from_scalar_quantities(
        massx/total_mass,
        massy/total_mass,
        massz/total_mass
    )
    
def kinetic_energy(particles):
    """
    Returns the total kinetic energy of the
    particles in the particles set.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.vx = [-1.0, 1.0] | units.ms
    >>> particles.vy = [0.0, 0.0] | units.ms
    >>> particles.vz = [0.0, 0.0] | units.ms
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.kinetic_energy()
    quantity<1.0 m**2 * kg * s**-2>
    """

    mass = particles.mass
    vx = particles.vx
    vy = particles.vy
    vz = particles.vz
    v_squared = (vx * vx) + (vy * vy) + (vz * vz)
    m_v_squared = mass * v_squared
    return 0.5 * m_v_squared.sum()
    

def potential_energy(particles, smoothing_length_squared = zero, G = constants.G):
    """
    Returns the total potential energy of the particles in the particles set.
    
    :argument smooting_length_squared: the smoothing length is added to every distance.
    :argument G: gravitational constant, need to be changed for particles in different units systems
    
    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.x = [0.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.potential_energy()
    quantity<-6.67428e-11 m**2 * kg * s**-2>
    """

    mass = particles.mass
    x_vector = particles.x
    y_vector = particles.y
    z_vector = particles.z
           
    sum_of_energies = zero
    
    for i in range(len(particles)):
       x = x_vector[i]
       y = y_vector[i]
       z = z_vector[i]
       dx = x - x_vector[i+1:]
       dy = y - y_vector[i+1:]
       dz = z - z_vector[i+1:]
       dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
       dr = (dr_squared+smoothing_length_squared).sqrt()
       m_m = mass[i] * mass[i+1:]
       
       energy_of_this_particle = (m_m / dr).sum()
       sum_of_energies -= energy_of_this_particle
        
    return G * sum_of_energies 



def particle_specific_kinetic_energy(set, particle):
    """
    Returns the kinetic energy of a particle.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.vx = [0.0, 1.0] | units.m
    >>> particles.vy = [0.0, 0.0] | units.m
    >>> particles.vz = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles[1].specific_kinetic_energy()
    quantity<0.5 m**2>
    """

    return 0.5*(particle.velocity**2).sum()

def particle_potential(set, particle, smoothing_length_squared = zero, gravitationalConstant = constants.G):
    """
    Returns the potential energy of a particle.

    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.x = [0.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles[1].potential()
    quantity<-6.67428e-11 m**2 * s**-2>
    """

    particles = set - particle
    dx = particle.x - particles.x
    dy = particle.y - particles.y
    dz = particle.z - particles.z 
    dr_squared = (dx * dx) + (dy * dy) + (dz * dz)
    dr = (dr_squared+smoothing_length_squared).sqrt()
    return - gravitationalConstant * (particles.mass / dr).sum()

def virial_radius(particles):
    """
    Returns the virial radius of the particles set.
    The virial radius is the inverse of the average inverse 
    distance between particles, weighted by their masses.
    
    >>> from amuse.support.data.core import Particles
    >>> particles = Particles(2)
    >>> particles.x = [-1.0, 1.0] | units.m
    >>> particles.y = [0.0, 0.0] | units.m
    >>> particles.z = [0.0, 0.0] | units.m
    >>> particles.mass = [1.0, 1.0] | units.kg
    >>> particles.virial_radius()
    quantity<4.0 m>
    """
    if len(particles) < 2:
        raise Exception("Cannot calculate virial radius for a particles set with fewer than 2 particles.")
    partial_sum = zero
    length_unit = particles.position.unit
    for i, particle_1 in enumerate(particles-particles[-1]):
        distance_vecs = [particle_1.position.number - particle_2.position.number 
            for particle_2 in particles[i+1:]] | length_unit
        partial_sum += ((particle_1.mass * particles[i+1:].mass) / distance_vecs.lengths()).sum()        
    return (particles.mass.sum()**2) / (2*partial_sum)


    
AbstractParticleSet.add_global_function_attribute("center_of_mass", center_of_mass)
AbstractParticleSet.add_global_function_attribute("center_of_mass_velocity", center_of_mass_velocity)
AbstractParticleSet.add_global_function_attribute("kinetic_energy", kinetic_energy)
AbstractParticleSet.add_global_function_attribute("potential_energy", potential_energy)
AbstractParticleSet.add_global_function_attribute("virial_radius", virial_radius)

AbstractParticleSet.add_global_vector_attribute("position", ["x","y","z"])
AbstractParticleSet.add_global_vector_attribute("velocity", ["vx","vy","vz"])

AbstractParticleSet.add_global_function_attribute("specific_kinetic_energy", None, particle_specific_kinetic_energy)
AbstractParticleSet.add_global_function_attribute("potential", None, particle_potential)


AbstractParticleSet.add_global_function_attribute("move_to_center", move_to_center)
AbstractParticleSet.add_global_function_attribute("scale_to_standard", scale_to_standard)
    
    
    
    
