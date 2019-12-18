from __future__ import print_function
import numpy
from amuse.lab import *
from amuse.couple import bridge
from matplotlib import pyplot

class MilkyWay_galaxy(object):
    def get_gravity_at_point(self, eps, x,y,z):
        phi_0 = self.get_potential_at_point(eps, x,y,z)
        grav = AdaptingVectorQuantity()
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0,x+dpos,y,z) - phi_0
        phi_dy = self.get_potential_at_point(0,x,y+dpos,z) - phi_0
        phi_dz = self.get_potential_at_point(0,x,y, z+dpos) - phi_0
        return phi_dx/dpos, phi_dy/dpos, phi_dz/dpos

    def disk_and_bulge_potentials(self, x,y,z, a, b, mass):
        r = (x**2+y**2).sqrt()
        return constants.G * mass /\
            (r**2 + (a + (z**2 + b**2).sqrt())**2).sqrt()

    def halo_potential(self, x,y,z, Mc=5.0E+10|units.MSun, Rc=1.0|units.kpc**2):
        r=(x**2+y**2+z**2).sqrt()
        rr = (r/Rc)
        return -constants.G * (Mc/Rc)*(0.5*numpy.log(1 +rr**2) + numpy.arctan(rr)/rr)

    def get_potential_at_point(self, eps, x, y, z):
        pot_disk = self.disk_and_bulge_potentials(x,y,z, 
            0.0|units.kpc, 0.277|units.kpc, 1.12E+10|units.MSun) 
        pot_bulge = self.disk_and_bulge_potentials(x,y,z, 
            3.7|units.kpc, 0.20|units.kpc, 8.07E+10|units.MSun) 
        pot_halo = self.halo_potential(x,y,z, 
            Mc=5.0E+10|units.MSun, Rc=6.0|units.kpc)
        return pot_disk + pot_bulge + pot_halo

def movie(time, sun_and_planets):

    R = [] | units.kpc
    for sp in sun_and_planets:
        R.append(sp.position.length())
    # - sun_and_planets.z
    print(R)
    pyplot.subplot(2,2,1)
    pyplot.scatter(sun_and_planets.x.value_in(units.kpc), 
                   sun_and_planets.y.value_in(units.kpc), 
                   c=['k', 'r'], s=10, lw=0)
    pyplot.subplot(2,2,2)
    pyplot.scatter(R.value_in(units.kpc), 
                   sun_and_planets.z.value_in(units.kpc), 
                   c=['k', 'r'], s=10, lw=0)
    pyplot.xlabel("R [kpc]")
    pyplot.ylabel("Z [kpc]")

    R = [] | units.kpc
    R.append((sun_and_planets[1].position-sun_and_planets[0].position).length())
    pyplot.subplot(2,2,3)
    pyplot.scatter(-time.value_in(units.Gyr), 
                   R.value_in(units.kpc), 
                   c=['k', 'r'], s=10, lw=0)
    pyplot.xlabel("t [Myr]")
    pyplot.ylabel("r [kpc]")

    pyplot.draw()

def main(t_end, filename):

    bodies = Particles(2)
    Sun = bodies[0]
    Sun.mass = 1|units.MSun
    Sun.position = (8.4, 0.0, 0.0) | units.kpc
    Sun.velocity = (-10.1, 235.5, 7.5) | units.kms # SPZ2009

    M67 = bodies[1]
    M67.mass = 50000 | units.MSun
    M67.position = Sun.position + ((0.766, 0.0, 0.49) |units.kpc) 
    M67.velocity = Sun.velocity + ((31.92, -21.66, -8.71) |units.kms)

    converter = nbody_system.nbody_to_si(bodies.mass.sum(), Sun.x)
    sunandm67 = Huayno(converter)
    sunandm67.particles.add_particle(bodies)
    channel_from_sunandm67 = sunandm67.particles.new_channel_to(bodies)

    gravity = bridge.Bridge()
    gravity.add_system(sunandm67, (MilkyWay_galaxy(),) )
    dt = 1|units.Myr
    gravity.timestep = dt 

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    t_end = 4.5|units.Gyr
    time = 0 * t_end

    if filename:
        write_set_to_file(bodies.savepoint(0.0 | t_end.unit), 
                          filename, "hdf5", 
                          append_to_file=False)
        pyplot.draw()
    else:
        R = [] | units.kpc
        for bi in bodies:
            R.append(bi.position.length())
        pyplot.ion()
        pyplot.scatter(R.value_in(units.kpc), 
                       bodies.z.value_in(units.kpc), 
                       c=['k', 'r'], s=10, lw=0)
        pyplot.xlabel("R [kpc]")
        pyplot.ylabel("Z [kpc]")
    while time < t_end:
        time += dt
        gravity.evolve_model(time)
        channel_from_sunandm67.copy()

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print("T=", time, "M=", bodies.mass.sum(), end=' ') 
        print("E= ", Etot, "Q= ", Ekin/Epot, end=' ')
        print("dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot) 
        Etot_prev = Etot
        if filename:
            write_set_to_file(bodies.savepoint(time), filename, "hdf5")
        else:
            R = [] | units.kpc
            for bi in bodies:
                R.append(bi.position.length())
            pyplot.scatter(R.value_in(units.kpc), 
                           bodies.z.value_in(units.kpc), 
                           c=['k', 'r'], s=10, lw=0)
            pyplot.draw()
    gravity.stop()
    
def new_option_parser():
   
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-t", unit=units.Gyr,
                      dest="t_end", type="float", default = 4.5|units.Gyr,
                      help="end time of the simulation [%default]")
    result.add_option("-f", dest="filename", default = "",
                      help="output filename")

    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

