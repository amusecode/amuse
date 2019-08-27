# -*- coding: ascii -*-
"""
Integrates a stellar orbit in the galactic potential

This example illustrates the use of a simple external potential and simple
integrator, no
amuse community code is used.
"""
from __future__ import print_function
# import numpy
from amuse.units.optparse import OptionParser
from math import atan
from math import log
from amuse.units import units, constants
from amuse.units.quantities import AdaptingVectorQuantity
from amuse.datamodel import Particle
from matplotlib import pyplot
from amuse.plot import plot


class MilkyWay_galaxy(object):

    def __init__(self, potential="point_particle", M=1.6e10 | units.MSun):
        print("pot=", potential)
        if potential.find("Milky") >= 0:
            print("Milky Way Potential")
            self.potential = self.Milky_Way_potential
        else:
            self.potential = self.point_particle_potential
        self.M = M

    # Generic function to get gravity at a point given the potential
    def get_gravity_at_point(self, pos):
        phi_0 = self.get_potential_at_point(pos)
        grav = AdaptingVectorQuantity()
        dpos = 0.001 * pos.length()
        for ii in range(len(pos)):
            ppos = 1.0 * pos
            ppos[ii] += dpos
            phi_1 = self.get_potential_at_point(ppos)
            grav.append((phi_1 - phi_0) / dpos)
        return grav

    def get_potential_at_point(self, pos):
        phi = self.potential(pos)
        return phi

    def point_particle_potential(self, pos):
        return self.Kepler_potential(self.M, pos)

    def Kepler_potential(self, mass, pos, eps2=0.0 | units.kpc**2):
        return self.Power_law_potential(mass, pos, eps2=eps2, eta=1)

    def Power_law_potential(self, mass, pos, eps2=0.0 | units.kpc**2, eta=1):
        return constants.G * mass / (pos.length()**2 + eps2)**(eta / 2.)

    def disk_and_bulge_potentials(self, pos, a, b, mass):
        r = (pos.x**2 + pos.y**2).sqrt()
        return constants.G * mass /\
            (r**2 + (a + (pos.z**2 + b**2).sqrt())**2).sqrt()

    def halo_potential(
            self, pos, Mc=5.0E+10 | units.MSun, Rc=1.0 | units.kpc**2):
        r = pos.length()
        rr = (r / Rc)
        return -constants.G * (Mc / Rc) * \
            (0.5 * log(1 + rr**2) + atan(rr) / rr)

    # 1990ApJ...348..485P
    def Milky_Way_potential(self, pos):
        pot_disk = self.disk_and_bulge_potentials(pos,
                                                  0.0 | units.kpc,
                                                  0.277 | units.kpc,
                                                  1.12E+10 | units.MSun)
        pot_bulge = self.disk_and_bulge_potentials(pos, 3.7 | units.kpc,
                                                   0.20 | units.kpc,
                                                   8.07E+10 | units.MSun)
        pot_halo = self.halo_potential(pos, Mc=5.0E+10 | units.MSun,
                                       Rc=6.0 | units.kpc)
        return pot_disk + pot_bulge + pot_halo


def new_single_star(mass, pos, vel):
    single_star = Particle()
    single_star.mass = mass
    single_star.position = pos
    single_star.velocity = vel
    single_star.radius = 1.0 | units.RSun
    return single_star


def evolve_particle_in_potential(single_star, potential, t_end):
    time = 0 | units.Myr
    dt_min = 0.1 * t_end
    while time < t_end:
        acc = potential.get_gravity_at_point(single_star.position)
        dt = min(dt_min, 0.1 * single_star.velocity.length() / acc.length())
        single_star.velocity += acc * dt
        single_star.position += single_star.velocity * dt
        time += dt


def evolve_particle_trajectory_in_potential(
        single_star, potential, dt_diag, t_end):
    time = 0 | units.Myr
    x = [] | size_unit
    y = [] | size_unit
    z = [] | size_unit
    while time < t_end:
        evolve_particle_in_potential(single_star, potential, dt_diag)
        time += dt_diag
        print(
            "time=", time,
            single_star.position.length().as_quantity_in(units.AU))
        x.append(single_star.x)
        y.append(single_star.y)
        z.append(single_star.y)
    return x, y, z


def plot_orbit(x, y):
    pyplot.figure(figsize=(10, 10))
    plot(x, y)
    pyplot.show()


def new_option_parser():
    result = OptionParser()
    result.add_option("-t", dest="t_end", type="float", default=250,
                      unit=units.Myr, help="end time [%unit]")
    result.add_option("-d", dest="dt_diag", type="float", default=10,
                      unit=units.Myr, help="diagnostic timestep [%unit]")
    result.add_option("-P", dest="potential", default="MilkyWay",
                      help="name of potential")
    result.add_option("-M", dest="mass", type="float", default=1.e+11,
                      unit=units.MSun, help="mass of the galaxy [%unit]")
    result.add_option("-e", dest="eps", type="float", default=0.0,
                      unit=units.parsec,
                      help="softening of the potential [%unit]")
    result.add_option("-x", dest="x", type="float", default=8500,
                      unit=units.parsec, help="x-position [%unit]")
    result.add_option("-y", dest="y", type="float", default=0,
                      unit=units.parsec, help="y-position [%unit]")
    result.add_option("-z", dest="z", type="float", default=0,
                      unit=units.parsec, help="z-position [%unit]")
    result.add_option("--vx", dest="vx", type="float", default=0,
                      unit=units.km / units.s, help="x-velocity [%unit]")
    result.add_option("--vy", dest="vy", type="float", default=220,
                      unit=units.km / units.s, help="y-velocity [%unit]")
    result.add_option("--vz", dest="vz", type="float", default=0,
                      unit=units.km / units.s, help="z-velocity [%unit]")
    return result


if __name__ == "__main__":
    o, arguments = new_option_parser().parse_args()

    t_end = o.t_end
    dt_diag = min(o.dt_diag, 0.1 * o.t_end)
    size_unit = units.parsec
    mass = 1.0
    pos = [o.x, o.y, o.z]
    vel = [o.vx, o.vy, o.vz]
    single_star = new_single_star(mass, pos, vel)
    galaxy = MilkyWay_galaxy(o.potential, M=o.mass)
    x, y, z = evolve_particle_trajectory_in_potential(
        single_star, galaxy, dt_diag, t_end)
    plot_orbit(x, y)
