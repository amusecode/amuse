from amuse.support.units import si

class Domain(object):
    pass
    
class Gravity(Domain):
    mass = 0.0 | si.g , "the mass of a star"
    position = [0.0, 0.0, 0.0] | si.m , "the position vector of a star"
    velocity = [0.0, 0.0, 0.0] | si.m / si.s , "the velocity vector of a star"
    radius = 0.0 | si.m , "the radius of a star"
    acceleration = [0.0, 0.0, 0.0] | si.m / (si.s ** 2), "the acceleraration vector of a star"
    
class StellarEvolution(Domain):
    mass = 0.0 | si.g , "the mass of a star"
    radius = 0.0 | si.m , "the radius of a star"
    age = 0.0 | si.s , "the age of a star, time evolved since star formation"


class SseCode(StellarEvolution):
    zams_mass = 0.0 | si.g , "the mass of a star after formation"
    type = 0 | si.no_unit, "stars evolve through typical stages, during each stage one can classify a star as belonging to a specific type"
    luminosity = 0.0 | si.cd / (si.m ** 2), "brightness of a star"
    radius = 0.0 | si.m, "total radius of a star"
    core_mass = 0.0 | si.g, "mass of the innermost layer of a star"
    core_radius = 0.0 | si.m, "radius of the innermost layer of a star"
    envelope_mass = 0.0 | si.g, "mass of the radiative / convective envelope around the core of a star"
    envelope_radius = 0.0 | si.m, "radius of the radiative / convective envelope around the core of a star"
    spin = 0.0 | si.m / si.s, "speed of rotation around the central axis of a star"
    epoch = 0.0 | si.s, "set when a star changes type"
    physical_time = 0.0 | si.s, "age of a star relative to last change of type"
