"""
   Simulate the radiative and hydrodynamial evolution of a disk with 
   a bump around a single star
"""
from time import time, localtime
from amuse.lab import *
from amuse.ext.molecular_cloud import ism_cube
from amuse.community.simplex.interface import SimpleXInterface, SimpleX, \
     SimpleXSplitSet
from amuse.ext.protodisk import ProtoPlanetaryDisk

set_printing_strategy("custom", 
                      preferred_units = [units.MSun, units.AU, units.Myr], 
                      precision = 12, prefix = "",
                      separator = " [", suffix = "]")

def mu(X = None, Y = 0.25, Z = 0.02, x_ion = 0.1):
    """
    Compute the mean molecular weight in kg (the average weight of
    particles in a gas) X, Y, and Z are the mass fractions of
    Hydrogen, of Helium, and of metals, respectively.  x_ion is the
    ionisation fraction (0 < x_ion < 1), 1 means fully ionised.
    """
    if X is None:
        X = 1.0 - Y - Z
    elif abs(X + Y + Z - 1.0) > 1e-6:
        raise Exception("Error in calculating mu: mass "
                         + "fractions do not sum to 1.0")
    return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0
                                     + Z*x_ion/2.0)

class RadHydro:
    def __init__(self, rad, hydro, star, disk):

        self.time = 0|units.day
        self.star = star
        self.disk = disk

        disk.r2 = disk.x**2 + disk.y**2
        disk = disk.sorted_by_attributes("r2")
        Rmax = disk.r2.max().sqrt()
        print "MaxR=", Rmax.in_(units.AU)

        self.hydro = hydro(nbody_system.nbody_to_si(self.disk.mass.sum(), Rmax),
                           number_of_workers=4)
        self.hydro.parameters.epsilon_squared = (10|units.AU)**2
        
        self.hydro.gas_particles.add_particles(self.disk)
        self.hydro.gas_particles.new_channel_to(self.disk)
        self.hydro.dm_particles.add_particles(self.star)
        self.hydro.dm_particles.new_channel_to(self.star)

        self.hydro_to_star = self.hydro.dm_particles.new_channel_to(self.star)
        self.hydro_to_disk = self.hydro.gas_particles.new_channel_to(self.disk)
        self.star_to_hydro = self.star.new_channel_to(self.hydro.dm_particles)
        self.disk_to_hydro = self.disk.new_channel_to(self.hydro.gas_particles)
            
        self.hydro.evolve_model(1|units.s)
        self.hydro_to_star.copy()
        self.hydro_to_disk.copy()

        self.rad = rad()
        for si in self.star:
            if si.mass>=5|units.MSun:
                self.rad.src_particles.add_particle(si)
        self.rad.gas_particles.add_particles(self.disk)
        self.rad.parameters.box_size=2.01*Rmax
        self.rad.parameters.timestep=1|units.day
        self.rad.set_source_Teff(star.temperature)

        self.rad_to_disk = self.rad.gas_particles.new_channel_to(self.disk,
                                                    attributes=["xion", "u"])
        self.star_to_rad = self.star.new_channel_to(self.rad.src_particles,
                                                    attributes=["x", "y", "z"])
        self.disk_to_rad = self.disk.new_channel_to(self.rad.gas_particles,
                                                    attributes=["x", "y", "z"])

        self.rad.stop()
        self.index = 0

    def write_file(self):
        self.index += 1
        filename = "hydro_disk_with_bump_i{0:04}.amuse".format(self.index)
        write_set_to_file(self.star, filename, "amuse", append_to_file=False)
        write_set_to_file(self.disk, filename, "amuse")
        
    def evolve_model(self, model_time):
        dt = model_time - self.time
        self.old_time = self.time
        self.time += dt/2.

        #self.disk_to_rad.copy()
        #self.star_to_rad.copy()
        #self.rad.evolve_model(self.time)
        #self.rad_to_disk.copy()
        
        self.time += dt/2.

        self.disk_to_hydro.copy()
        self.star_to_hydro.copy()
        self.hydro.evolve_model(self.time)
        self.hydro_to_disk.copy()
        self.hydro_to_star.copy()

        print "RT done at time:", self.time.in_(units.day)

    def print_diagnostics(self):
        umin = self.disk.u.min()
        umean = self.disk.u.mean()
        umax = self.disk.u.max()
        Tmin =  mu() / constants.kB * umax
        Tmean =  mu() / constants.kB * umean
        Tmax =  mu() / constants.kB * umin

        print "Time=", self.time.in_(units.day)
        print "Ionization:", self.disk.xion.min(), self.disk.xion.mean(), \
              self.disk.xion.max()
        print "Intenal energy:", umin, umean, umax
        print "Temperature:", Tmin, Tmean, Tmax
        print "Density:", self.disk.density.min().in_(units.amu/units.cm**3), \
              self.disk.density.mean().in_(units.amu/units.cm**3), \
              self.disk.density.max().in_(units.amu/units.cm**3)
        print "scaleheight:", abs(self.disk.z.value_in(units.AU)).mean()
        
    def stop(self):
        self.hydro.stop()

###BOOKLISTSTART1###
def new_disk_with_bump(Mstar = 10|units.MSun,
                       Ndisk=100, Mdisk=1.0|units.MSun, 
                       Rmin=1.0|units.AU, Rmax=100.0|units.AU, 
                       Mbump=0.1|units.MSun,Rbump=5.0|units.AU, 
                       abump=10|units.AU):

    converter=nbody_system.nbody_to_si(Mstar, Rmin)
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter,
                              densitypower=1.5, Rmin=1, Rmax=Rmax/Rmin,
                              q_out=1.0, discfraction=Mdisk/Mstar).result
    com = disk.center_of_mass()

    # determine bump's local velocity
    
    inner_particles = disk.select(lambda r: (com-r).length()
                                                   < abump,["position"])
    M_inner = Mstar + inner_particles.mass.sum() 
    v_circ = (constants.G*M_inner*(2./abump - 1./abump)) \
                 .sqrt().value_in(units.kms)

    # initialize bump
    
    Nbump = int(Ndisk*Mbump/Mdisk)
    bump = new_plummer_gas_model(Nbump,
                                 convert_nbody=nbody_system.nbody_to_si(Mbump,
                                                                        Rbump))
    bump.x += abump
    bump.velocity += [0, v_circ, 0] | units.kms

    disk.add_particles(bump)
    disk.move_to_center()
    return disk
###BOOKLISTSTOP1###

def evolve_star(Mstar, tstar):
    stars = Particles(1)
    stars = Particles(2)
    stars[0].mass = Mstar
    stars[1].mass = 0.1*Mstar
    stellar = SeBa()
    stellar.particles.add_particle(stars)
    stellar.evolve_model(tstar)
    stars.mass = stellar.particles.mass
    stars.position = (0, 0, 0) |units.AU
    stars.velocity = (0, 0, 0) |units.kms
    stars.luminosity = stellar.particles.luminosity/(20. | units.eV)
    stars.temperature = stellar.particles.temperature
    stars.flux = stars.luminosity
    stars.rho = 1.0|(units.g/units.cm**3)
    stars.xion = 0.0 #ionization_fraction
    stars.u = (9. |units.kms)**2 #internal_energy
    print stars

    if len(stars)>1:
        stars[1].x = 50|units.AU
        vc = 1.0*(constants.G*stars.mass.sum()/(100.|units.AU)).sqrt()
        stars[1].vy += vc
    
    stellar.stop()
    return stars

def hydro_disk_with_bump(Mstar = 10|units.MSun,
                         Ndisk=100,
                         Mdisk=1.0|units.MSun,
                         Rmin=1.0|units.AU,
                         Rmax=100.0|units.AU,
                         Mbump=0.1|units.MSun,
                         Rbump=5.0|units.AU,
                         abump=10|units.AU,
                         t_end=10|units.yr,
                         n_steps=10):

    star = evolve_star(Mstar, t_end)

    disk = new_disk_with_bump(star[0].mass, Ndisk, Mdisk, Rmin, Rmax,
                              Mbump, Rbump, abump)

    radhydro = RadHydro(SimpleXSplitSet, Gadget2, star, disk)
    radhydro.write_file()

    dt = t_end/float(n_steps)
    print "dt=", dt.in_(units.day)
    time = 0 | units.day
    while time<t_end:
        time += dt
        radhydro.evolve_model(time)
        radhydro.print_diagnostics()
        radhydro.write_file()
    radhydro.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="Ndisk", type="int",default = 10000,
                      help="number of stars [%default]")
    result.add_option("-t", unit=units.yr,
                      dest="t_end", type="float", default = 2000.|units.yr,
                      help="radiation time [%default]")
    result.add_option("-n", dest="n_steps", type="int",default = 100,
                      help="number of steps [%default]")
    result.add_option("--Mstar", unit=units.MSun,
                      dest="Mstar", type="float", default = 10|units.MSun,
                      help="Mass of the central star [%default]")
    result.add_option("--Mdisk", unit=units.MSun, 
                      dest="Mdisk", type="float", default = 1|units.MSun,
                      help="Mass of the disk [%default]")
    result.add_option("-r", unit=units.AU,
                      dest="Rmin", type="float", default = 10 |units.AU,
                      help="inner disk radius [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Rmax", type="float", default = 100 | units.AU,
                      help="outer disk radius [%default]")
    result.add_option("--Mbump", unit=units.MSun,
                      dest="Mbump", type="float", default = 0.5 | units.MSun,
                      help="bump mass [%default]")
    result.add_option("--Rbump", unit=units.AU,
                      dest="Rbump", type="float", default = 5 | units.AU,
                      help="bump radius [%default]")
    result.add_option("-a", unit=units.AU,
                      dest="abump", type="float", default = 50 | units.AU,
                      help="distance of bump from star [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    hydro_disk_with_bump(**o.__dict__)

