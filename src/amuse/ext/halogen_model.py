from amuse.community.halogen.interface import Halogen
from amuse.datamodel import ParticlesWithUnitsConverted

def new_halogen_model(number_of_particles, convert_nbody = None, do_scale = False, 
        redirection = 'null', **keyword_arguments):
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
    :argument convert_nbody:  When given will convert the resulting set to SI units
    :argument do_scale: scale the result to exact nbody units (M=1, K=0.25, U=-0.5)
    :argument alpha: alpha parameter in density profile (see amuse/community/halogen/src/doc for details)
    :argument beta:  beta parameter in density profile (see amuse/community/halogen/src/doc for details)
    :argument gamma: gamma parameter in density profile (see amuse/community/halogen/src/doc for details)
    """
    instance = Halogen(unit_converter=convert_nbody, redirection=redirection)
    instance.parameters.number_of_particles = number_of_particles
    for (key, value) in keyword_arguments.items():
        setattr(instance.parameters, key, value)
    
    instance.generate_particles()
    result = instance.particles.copy()
    instance.stop()
    
    result.move_to_center()
    if do_scale:
        result.scale_to_standard()
    
    if not convert_nbody is None:
        result = ParticlesWithUnitsConverted(result, convert_nbody.as_converter_from_si_to_generic())
        result = result.copy()
    return result

