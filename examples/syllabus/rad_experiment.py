"""
   Simulate the radiative and hydrodynamial evolution of a disk with 
   a single bump around a single star
"""
from time import time, localtime
from amuse.lab import *
from amuse.ext.molecular_cloud import ism_cube
from amuse.community.simplex.interface import SimpleXInterface, SimpleX,SimpleXSplitSet

set_printing_strategy("custom", 
                      preferred_units = [units.MSun, units.AU, units.Myr], 
                      precision = 5, prefix = "", separator = " [", suffix = "]"
)

def new_disk_with_bump(Mstar = 1|units.MSun,
                       Ndisk=100, Mdisk=0.9|units.MSun, 
                       Rmin=1.0|units.AU, Rmax=100.0|units.AU, 
                       Mbump=0.1|units.MSun,Rbump=10.0|units.AU, 
                       abump=10|units.AU):

    converter=nbody_system.nbody_to_si(Mdisk, Rmin)
    from amuse.ext.protodisk import ProtoPlanetaryDisk
    bodies = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter, 
                                densitypower=1.5, 
                                Rmin=1,
                                Rmax=Rmax/Rmin,
                                q_out=1.0,
                                discfraction=1.0).result
    Mdisk = bodies.mass.sum()
    bodies.move_to_center()
    com = bodies.center_of_mass()

    mm = Mdisk/float(Ndisk)
    Nbump = Mbump/mm
    bump = new_plummer_gas_model(Nbump, convert_nbody=nbody_system.nbody_to_si(Mbump, Rbump))

    bump.x += abump
    r_bump = abump 
    inner_particles = bodies.select(lambda r: (com-r).length()<abump,["position"])
    M_inner = inner_particles.mass.sum() + Mstar

    v_circ = (constants.G*M_inner*(2./r_bump - 1./abump)).sqrt().value_in(units.kms)
    bump.velocity += [0, v_circ, 0] | units.kms
    bodies.add_particles(bump)

    bodies = bodies.select(lambda r: (com-r).length()>Rmin,["position"])
    bodies = bodies.select(lambda r: (com-r).length()<Rmax,["position"])

    print "Nbump=", Ndisk, Mbump, Rbump, Nbump
    print "Mass =", Mstar, Mdisk, bodies.mass.sum().in_(units.MSun), bodies.mass.sum()/Mstar

    print "X-min/max:", min(bodies.x), max(bodies.x)
    bodies = bodies.select(lambda r: r.length()<1.0*Rmax,["position"])
    print "X-min/max:", min(bodies.x), max(bodies.x)

    return bodies

def read_disk_with_bump(image_id=1, filename="hydro.hdf5"):

    star = read_set_from_file(filename, "hdf5")
    ism = read_set_from_file(filename, "hdf5")
    snapshot_id = 0
    for si in ism.history:
        snapshot_id += 1
        print "reading snapshot_id=", snapshot_id
        if image_id<0 or image_id == snapshot_id:
            com = si.center_of_mass()
            return si

def irradiate_disk_with_bump(Mstar = 10|units.MSun, tstar = 20|units.Myr, Ndisk=100, Mdisk=0.9|units.MSun,  Rmin=1.0|units.AU, Rmax=100.0|units.AU,  Mbump=0.1|units.MSun,Rbump=10.0|units.AU, abump=10|units.AU, t_end=10|units.yr, n_steps=10):
    dt = t_end/float(n_steps)

    stellar = SeBa()
    stellar.particles.add_particle(Particle(mass=Mstar))
    stellar.evolve_model(tstar)

    source=Particles(1)
    source.mass = stellar.particles[0].mass
    source.position = (0, 0, 0) |units.AU
    source.velocity = (0, 0, 0) |units.kms
    source.luminosity = stellar.particles[0].luminosity/(20. | units.eV)
    source.temperature = stellar.particles[0].temperature
    Teff = stellar.particles[0].temperature
    source.flux = source.luminosity
    source.rho = 1.0|(units.g/units.cm**3)
    source.xion = 0.0 #ionization_fraction
    source.u = (9. |units.kms)**2 #internal_energy
    stellar.stop()

    ism = new_disk_with_bump(source[0].mass, Ndisk, Mdisk, Rmin, Rmax, Mbump, Rbump, abump)
    ism.flux = 0 | units.s**-1
    ism.xion = 0.0 #ionization_fraction

    hydro = Gadget2(nbody_system.nbody_to_si(Mdisk, Rmax))
    hydro.gas_particles.add_particles(ism)
    hydro.dm_particles.add_particles(source)
    hydro.evolve_model(1|units.day)
    hydro.gas_particles.new_channel_to(ism).copy()
    hydro.stop()
    
    rad = SimpleXSplitSet(redirect="none",numer_of_workers=4)
    rad.parameters.box_size=2.01*Rmax
    rad.parameters.timestep=0.1*dt
    rad.set_source_Teff(Teff)
    rad.src_particles.add_particle(source)
    rad.gas_particles.add_particles(ism)

    rad_to_framework = rad.gas_particles.new_channel_to(ism)
    particles = ParticlesSuperset([source, ism])
    write_set_to_file(particles, "rad.hdf5", 'hdf5')

    while rad.model_time<t_end:
        rad.evolve_model(rad.model_time + dt)
        rad_to_framework.copy_attributes(["x","y", "z", "xion"])
        write_set_to_file(particles, "rad.hdf5", 'hdf5')
        print "Time=", rad.model_time, "Ionization (min, mean, max):", ism.xion.min(), ism.xion.mean(), ism.xion.max()
    rad.stop()

def _irradiate_disk_with_pump(Mstar = 10|units.MSun,
         tstar = 20|units.Myr,
         Ndisk=100, Mdisk=0.9|units.MSun, 
         Rmin=1.0|units.AU, Rmax=100.0|units.AU, 
         Mbump=0.1|units.MSun,Rbump=10.0|units.AU, abump=10|units.AU,
         t_end=10|units.yr, n_steps=10, filename = None, image_id=1):
    model_time = 0.0 | t_end.unit
    dt = t_end/float(n_steps)

    ionization_fraction = 0.0
    internal_energy = (9. |units.kms)**2

    stellar = SeBa()
    stellar.particles.add_particle(Particle(mass=Mstar))
    stellar.evolve_model(tstar)

    print "L=", stellar.particles[0].luminosity.in_(units.LSun), stellar.particles[0].temperature

    source=Particles(1)
    source.mass = stellar.particles[0].mass
    source.position = (0, 0, 0) |units.AU
    source.velocity = (0, 0, 0) |units.kms
    source.luminosity = stellar.particles[0].luminosity/(20. | units.eV)
    source.temperature = stellar.particles[0].temperature
    Teff = stellar.particles[0].temperature
    source.flux = source.luminosity
    source.rho = 1.0|(units.g/units.cm**3)
    source.xion = ionization_fraction
    source.u = internal_energy
    stellar.stop()

    Mstar = source[0].mass
    if filename ==None:
        ism = new_disk_with_bump(Mstar = Mstar,
                                 Ndisk=Ndisk, Mdisk=Mdisk,
                                 Rmin=Rmin, Rmax=Rmax,
                                 Mbump=Mbump, Rbump=Rbump,
                                 abump=abump)
    else:
        ism = read_disk_with_bump(image_id=image_id, filename=filename)#, Rmax=Rmax)
    ism = ism.select(lambda r: r.length()<1.0*Rmax,["position"])
    ism.flux = 0 | units.s**-1
    ism.xion = ionization_fraction 

    if filename ==None:
        converter=nbody_system.nbody_to_si(Mdisk, Rmax)
        hydro = Gadget2(converter)
        hydro.gas_particles.add_particles(ism)
        hydro.dm_particles.add_particles(source)
        hydro.evolve_model(1|units.day)
        hydro.gas_particles.new_channel_to(ism).copy()
        hydro.stop()
    
    rad = SimpleXSplitSet(redirect="none",numer_of_workers=4)
    rad.parameters.number_of_freq_bins=5
    rad.parameters.thermal_evolution_flag=1
    rad.parameters.blackbody_spectrum_flag=1
    rad.parameters.metal_cooling_flag=0
    rad.parameters.box_size=2.01*Rmax
    if isinstance(rad, SPHRay):
        ism.h_smooth = 0.1 | (units.RSun)

    else:
        rad.parameters.timestep=0.1*dt
        rad.set_source_Teff(Teff)

    rad.src_particles.add_particle(source)
    rad.gas_particles.add_particles(ism)

    channel_from_rad_to_framework = rad.gas_particles.new_channel_to(ism)

    particles = ParticlesSuperset([source, ism])
    write_set_to_file(particles, "rad.hdf5", 'hdf5')

    tCPU = time()
    while model_time<t_end:
        model_time += dt

        rad.evolve_model(model_time)
        channel_from_rad_to_framework.copy_attributes(["x","y", "z", "xion"])
        write_set_to_file(particles, "rad.hdf5", 'hdf5')

        print "Date:"+str(localtime()[2])+"."+str(localtime()[1])+"."+str(localtime()[1]), "at", str(localtime()[3])+"h", str(localtime()[4])+"m"
        print "Time=", model_time, "dt_CPU=", time()-tCPU
        print "min ionization:", ism.xion.min()
        print "average Xion:", ism.xion.mean()
        print "max ionization:", ism.xion.max()
        tCPU = time()

    rad.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="Ndisk", type="int",default = 1000,
                      help="number of stars [%default]")
    result.add_option("-t", unit=units.yr,
                      dest="t_end", type="float", default = 1|units.yr,
                      help="radiation time [%default]")
    result.add_option("-f",
                      dest="filename", default = None,
                      help="input filename [%default]")
    result.add_option("--tstar", unit=units.Myr,
                      dest="tstar", type="float", default = 20|units.Myr,
                      help="age of the star [%default]")
    result.add_option("-n", dest="n_steps", type="int",default = 36500,
                      help="number of steps [%default]")
    result.add_option("-i", dest="image_id", type="int",default = 1,
                      help="id of the input (filename) snapshot [%default]")
    result.add_option("--Mstar", unit=units.MSun,
                      dest="Mstar", type="float", default = 10|units.MSun,
                      help="Mass of the central star [%default]")
    result.add_option("--Mdisk", unit=units.MSun, 
                      dest="Mdisk", type="float", default = 0.0002|units.MSun,
                      help="Mass of the disk [%default]")
    result.add_option("-r", unit=units.AU,
                      dest="Rmin", type="float", default = 10 |units.AU,
                      help="inner disk radius [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Rmax", type="float", default = 100 | units.AU,
                      help="outer disk radius [%default]")
    result.add_option("--Mbump", unit=units.MSun,
                      dest="Mbump", type="float", default = 0.0001 | units.MSun,
                      help="bump mass [%default]")
    result.add_option("--Rbump", unit=units.AU,
                      dest="Rbump", type="float", default = 5 | units.AU,
                      help="bump radius [%default]")
    result.add_option("-a", unit=units.AU,
                      dest="abump", type="float", default = 10 | units.AU,
                      help="distance of bump from star [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    irradiate_disk_with_bump(**o.__dict__)

