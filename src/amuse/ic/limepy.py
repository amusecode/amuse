from amuse.support.exceptions import AmuseWarning
from amuse.datamodel import Particles, ParticlesWithUnitsConverted
from amuse.units import nbody_system
from amuse.community import LiteratureReferencesMixIn

try:
    from ._limepy import limepy, sample
    scipy_imported = True
except ImportError:
    raise
    scipy_imported = False

    AmuseWarning("import limepy failed, maybe scipy is not installed")


class Limepy(LiteratureReferencesMixIn):
    """
    LIMEPY : Lowered Isothermal Model Explorer in PYthon
    for help:
    print help(limepy.limepy)
    print help(limepy.sample)

    Relevant references:
        .. [#] Gieles & Zocchi 2015, MNRAS, 454,576
    """

    def __init__(self, *args, **kwargs):
        LiteratureReferencesMixIn.__init__(self)
        kwargs["M"] = 1
        kwargs["G"] = 1
        kwargs["rv"] = 1

        self.model = limepy(*args, **kwargs)
        self.kwargs = kwargs



    @property
    def result(self):
        stars = sample(self.model, **self.kwargs)
        self.sample = stars
        p = Particles(stars.N)
        p.mass = stars.m | nbody_system.mass
        p.x = stars.x | nbody_system.length
        p.y = stars.y | nbody_system.length
        p.z = stars.z | nbody_system.length
        p.vx = stars.vx | nbody_system.length / nbody_system.time
        p.vy = stars.vy | nbody_system.length / nbody_system.time
        p.vz = stars.vz | nbody_system.length / nbody_system.time
        return p


def new_limepy_model(*args, **kwargs):
    conv = kwargs.pop("converter", None)

    l = Limepy(*args, **kwargs)
    p = l.result

    if conv is not None:
        p = ParticlesWithUnitsConverted(
            p, conv.as_converter_from_si_to_generic())

    return p
