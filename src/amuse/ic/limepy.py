from amuse.support.exceptions import AmuseWarning
from amuse.datamodel import Particles, ParticlesWithUnitsConverted
from amuse.units import nbody_system

try:
    from _limepy import limepy, sample
    scipy_imported = True
except ImportError:
    scipy_imported = False

    raise AmuseWarning("import limepy failed, maybe scipy is not installed")


class Limepy(object):
    __doc__ = limepy.__init__.__doc__ + sample.__init__.__doc__

    def __init__(self, *args, **kwargs):
        kwargs["scale"] = True
        kwargs["MS"] = 1
        kwargs["GS"] = 1
        kwargs["RS"] = 1
        kwargs["scale_radius"] = "rv"

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
