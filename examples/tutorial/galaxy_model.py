import numpy
from amuse.lab import (
        # AdaptingVectorQuantity,
        constants, units,
        set_printing_strategy,
        )


class MilkyWay_galaxy(object):
    def get_gravity_at_point(self, eps, x, y, z):
        phi_0 = self.get_potential_at_point(eps, x, y, z)
        # grav = AdaptingVectorQuantity()
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0, x+dpos, y, z) - phi_0
        phi_dy = self.get_potential_at_point(0, x, y+dpos, z) - phi_0
        phi_dz = self.get_potential_at_point(0, x, y, z+dpos) - phi_0
        return phi_dx/dpos, phi_dy/dpos, phi_dz/dpos

    def disk_and_bulge_potentials(self, x, y, z, a, b, mass):
        r = (x**2+y**2).sqrt()
        return constants.G * mass /\
            (r**2 + (a + (z**2 + b**2).sqrt())**2).sqrt()

    def halo_potential(
            self, x, y, z, Mc=5.0E+10 | units.MSun, Rc=1.0 | units.kpc**2):
        r = (x**2+y**2+z**2).sqrt()
        rr = (r/Rc)
        return (
                -constants.G * (Mc/Rc)*(
                    0.5*numpy.log(1 + rr**2) + numpy.arctan(rr)/rr
                    )
                )

    def get_potential_at_point(self, eps, x, y, z):
        pot_disk = self.disk_and_bulge_potentials(
                x, y, z,
                0.0 | units.kpc, 0.277 | units.kpc, 1.12E+10 | units.MSun)
        pot_bulge = self.disk_and_bulge_potentials(
                x, y, z,
                3.7 | units.kpc, 0.20 | units.kpc, 8.07E+10 | units.MSun)
        pot_halo = self.halo_potential(
                x, y, z,
                Mc=5.0E+10 | units.MSun, Rc=6.0 | units.kpc)
        return pot_disk + pot_bulge + pot_halo


if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom",  # nbody_converter = converter,
                          preferred_units=[units.MSun, units.RSun, units.yr],
                          precision=4, prefix="",
                          separator=" [", suffix="]")
    mwg = MilkyWay_galaxy()
    sun_pos = [8.5, 0, 0] | units.kpc
    eps = 1 | units.AU
    print "Milky Way Galaxy:"
    print "gravity at solar location:", \
        mwg.get_gravity_at_point(eps, sun_pos[0], sun_pos[1], sun_pos[2])
    print "potential at solar location:", \
        mwg.get_potential_at_point(eps, sun_pos[0], sun_pos[1], sun_pos[2])
