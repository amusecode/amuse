from amuse.support.data.core import ParticlesWithUnitsConverted
from amuse.community.halogen.interface import Halogen

def new_halogen_model(number_of_particles, unit_converter = None, do_scale = False, **keyword_arguments):
    """
    Create an alpha-beta-gamma-model using Halogen with the given number of 
    particles. Returns a set of equal-mass particles self-consistently sampled 
    from the spherically symmetric density distribution defined by the alpha, 
    beta, and gamma parameters. The model is centered around the origin. 
    Positions and velocities are optionally scaled such that the kinetic and 
    potential energies are 0.25 and -0.5 in nbody-units, respectively.
    
    The alpha, beta, and gamma parameters are (of course) required, but all 
    other Halogen parameters can be used too, e.g. 
    new_halogen_model(..., black_hole_mass = 1.0e6 | units.MSun)
    will set halogen.parameters.black_hole_mass to this value. See 
    help(Halogen().parameters) for an overview of the Halogen parameters.

    :argument number_of_particles: Number of particles to generate in the model
    :argument unit_converter:  When given will convert the resulting set to SI units
    :argument do_scale: scale the result to exact nbody units (M=1, K=0.25, U=-0.5)
    :argument alpha: alpha parameter in density profile (see amuse/community/halogen/src/doc for details)
    :argument beta:  beta parameter in density profile (see amuse/community/halogen/src/doc for details)
    :argument gamma: gamma parameter in density profile (see amuse/community/halogen/src/doc for details)
    """
    instance = Halogen(unit_converter = unit_converter)
    instance.parameters.number_of_particles = number_of_particles
    for (key, value) in keyword_arguments.iteritems():
        setattr(instance.parameters, key, value)
    
    instance.generate_particles()
    result = instance.particles.copy()
    instance.stop()
    
    result.move_to_center()
    if do_scale:
        result.scale_to_standard()
    
    if not unit_converter is None:
        result = ParticlesWithUnitsConverted(result, unit_converter.as_converter_from_si_to_nbody())
        result = result.copy_to_memory()
    return result

