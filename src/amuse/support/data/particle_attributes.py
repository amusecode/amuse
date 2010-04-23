from amuse.support.units import nbody_system
from amuse.support.data import core
from amuse.support.data.values import zero

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
        particles = core.ParticlesWithUnitsConverted(particles, convert_nbody.as_converter_from_nbody_to_si())
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
    

core.AbstractParticleSet.add_global_function_attribute("move_to_center", move_to_center)
core.AbstractParticleSet.add_global_function_attribute("scale_to_standard", scale_to_standard)
    
    
    
    
