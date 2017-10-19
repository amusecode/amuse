import numpy
from amuse.lab import *
from amuse.couple import bridge

class MilkyWay_galaxy(object):
    def get_gravity_at_point(self, eps, x,y,z):
        phi_0 = self.get_potential_at_point(eps, x,y,z)
        grav = AdaptingVectorQuantity()
        dpos = 0.001*(x**2+y**2+z**2).sqrt()
        phi_dx = self.get_potential_at_point(0,x+dpos,y,z) - phi_0
        phi_dy = self.get_potential_at_point(0,x,y+dpos,z) - phi_0
        phi_dz = self.get_potential_at_point(0,x,y, z+dpos) - phi_0
        return phi_dx/dpos, phi_dy/dpos, phi_dz/dpos

    def disk_or_bulge_potentials(self, x, y, z, a, b, mass):
        r2 = x**2 + y**2
        b2 = (a + (z**2 + b**2).sqrt())**2
        return constants.G * mass / (r2 + b2).sqrt()

    def halo_potential(self, x,y,z, Mc=5.0E+10|units.MSun, Rc=1.0|units.kpc**2):
        r=(x**2 + y**2 + z**2).sqrt()
        rr = (r/Rc)
        return -constants.G * (Mc/Rc)*(0.5*numpy.log(1 +rr**2) + numpy.arctan(rr)/rr)

    def get_potential_at_point(self, eps, x, y, z):
        pot_disk = self.disk_or_bulge_potentials(x,y,z, 
            0.0|units.kpc, 0.277|units.kpc, 1.12E+10|units.MSun) 
        pot_bulge = self.disk_or_bulge_potentials(x,y,z, 
            3.7|units.kpc, 0.20|units.kpc, 8.07E+10|units.MSun) 
        pot_halo = self.halo_potential(x,y,z, 
            Mc=5.0E+10|units.MSun, Rc=6.0|units.kpc)
        return pot_disk + pot_bulge + pot_halo

def main(Ncl, mcl, rcl, W0, Rgal, vgal, t_end, n_steps):

    bodies = Particles(2)
    Sun = bodies[0]
    v_LSR = (-10, 5.2, 7.2) | units.kms
    Sun.mass = 1|units.MSun
    Sun.position = (8.4, 0.0, 0.0) | units.kpc
    Sun.velocity = (-10.1, 235.5, 7.5) | units.kms # SPZ2009

    M67 = bodies[1]
    M67.mass = 50000 | units.MSun
    M67.position = Sun.position + ((0.766, 0.0, 0.49) |units.kpc) 
    M67.velocity = Sun.velocity + ((31.92, -21.66, -8.71) |units.kms)

    Sun.velocity *= -1
    M67.velocity *= -1

    converter = nbody_system.nbody_to_si(bodies.mass.sum(), Sun.x)

    sunandm67 = ph4(converter)
    sunandm67.particles.add_particle(bodies)
    channel_from_sunandm67 = sunandm67.particles.new_channel_to(bodies)

    gravity = bridge.Bridge()
    gravity.add_system(sunandm67, (MilkyWay_galaxy(),) )
    dt = 0.1|units.Myr
    gravity.timestep = dt 
    
    filename="sunandM67.hdf5"
    write_set_to_file(bodies.savepoint(0.0 | t_end.unit), filename, "hdf5", append_to_file=False)

    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    t_end = 4.56|units.Gyr
    time = 0 | units.yr
    while time < t_end:
        time += 10*dt

        gravity.evolve_model(time)
        channel_from_sunandm67.copy()
        write_set_to_file(bodies.savepoint(time), filename, "hdf5")

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        print "T=", time, "M=", bodies.mass.sum(), 
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot 
        Etot_prev = Etot
    gravity.stop()
    
    
def new_option_parser():
   
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="Ncl", type="int",default = 100,
                      help="number of stars [%default]")
    result.add_option("-t", unit=units.Gyr,
                      dest="t_end", type="float", default = 4.5|units.Gyr,
                      help="end time of the simulation [%default]")
    result.add_option("-n", dest="n_steps", type="float", default = 300,
                      help="number of diagnostics output steps [%default]")
    result.add_option("-m", unit=units.parsec,
                      dest="mcl", type="float", default = 10**7|units.MSun,
                      help="cluster mass [%default]")
    result.add_option("-r", unit=units.parsec,
                      dest="rcl", type="float", default = 100|units.parsec,
                      help="cluster half-mass radius [%default]")
    result.add_option("-R", unit=units.kpc,
                      dest="Rgal", type="float", default = 8.5|units.kpc,
                      help="distance to the GC [%default]")
    result.add_option("-v", unit=units.parsec,
                      dest="vgal", type="float", default = 100|units.kms,
                      help="orbital velocity of the CDG [%default]")
    result.add_option("-W", dest="W0", type="float", default = 7.0,
                      help="Dimension-less depth of the King potential (W0) [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.RSun, units.yr], 
                          precision = 4, prefix = "", 
                          separator = " [", suffix = "]")
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

