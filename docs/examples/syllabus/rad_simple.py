from amuse.lab import *
from amuse.ext.molecular_cloud import ism_cube

from amuse.ext.protodisk import ProtoPlanetaryDisk
from amuse.community.simplex.interface import SimpleXInterface, SimpleX,SimpleXSplitSet
def new_disk_with_bump(Mstar = 1|units.MSun,
                       Ndisk=100, Mdisk=0.9|units.MSun, 
                       Rmin=1.0|units.AU, Rmax=100.0|units.AU, 
                       Mbump=0.1|units.MSun,Rbump=10.0|units.AU, 
                       abump=10|units.AU):
    converter=nbody_system.nbody_to_si(Mdisk, Rmin)
    bodies = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter, 
                                densitypower=1.5, 
                                Rmin=1.0, 
                                Rmax=Rmax/Rmin,
                                q_out=1.0,
                                discfraction=1.0).result
    Mdisk = bodies.mass.sum()
    bodies.move_to_center()
    com = bodies.center_of_mass()

    mm = Mdisk/float(Ndisk)
    Nbump = Mbump/mm
    print "Nbump=", Mbump, Rbump, Nbump
    print "Mass =", Mstar, Mdisk, bodies.mass.sum().in_(units.MSun), bodies.mass.sum()/Mstar
    bump = new_plummer_gas_model(Nbump, convert_nbody=nbody_system.nbody_to_si(Mbump, Rbump))

    bump.x += abump
    r_bump = abump 
    inner_particles = bodies.select(lambda r: (com-r).length()<abump,["position"])
    M_inner = inner_particles.mass.sum() + Mstar

    v_circ = (constants.G*M_inner*(2./r_bump - 1./abump)).sqrt().value_in(units.kms)
    bump.velocity += [0, v_circ, 0] | units.kms
    bodies.add_particles(bump)
    return bodies

def main(N=1000, Lstar=100| units.LSun, boxsize=10| units.parsec, 
         rho=1.0| (units.amu/units.cm**3), t_end=0.1|units.Myr, n_steps=10):

    ionization_fraction = 0.0
    internal_energy = (9. |units.kms)**2

    source=Particles(1)
    source.position = (0, 0, 0) |units.parsec
    source.flux = Lstar/(20. | units.eV)
    source.luminosity = Lstar/(20. | units.eV)
    source.rho = rho
    source.xion = ionization_fraction
    source.u = internal_energy

    converter=nbody_system.nbody_to_si(1|units.MSun, boxsize)
    ism = ProtoPlanetaryDisk(N, convert_nbody=converter, 
                                densitypower=1.5, 
                                Rmin=0.1, 
                                Rmax=1,
                                q_out=1.0,
                                discfraction=1.0).result
    ism = ism.select(lambda r: r.length()<0.5*boxsize,["position"])
    gamma=5./3.
    mu=1.| units.amu
    Tinit = 10000|units.K
    ism.u = 1/(gamma-1)*constants.kB * Tinit/mu
    ism.rho = rho
    ism.flux = 0. | units.s**-1
    ism.xion = ionization_fraction 
    ism.h_smooth = 0 | units.AU

    rad = SimpleXSplitSet(redirect="none")
#    rad = SimpleX()
#    rad = SPHRay()
    rad.parameters.box_size=1.001*boxsize    
    rad.parameters.timestep=0.001 | units.Myr

    rad.src_particles.add_particle(source)
    rad.gas_particles.add_particles(ism)

    channel_to_local_gas = rad.gas_particles.new_channel_to(ism)
    write_set_to_file(ism, "rad.hdf5", 'hdf5')

    time = 0.0 | t_end.unit
    dt = t_end/float(n_steps)
    while time<t_end:
        time += dt
        rad.evolve_model(time)
        channel_to_local_gas.copy_attributes(["xion",])
        write_set_to_file(ism, "rad.hdf5", 'hdf5')

        print "Time=", time
        print "min ionization:", rad.gas_particles.xion.min()
        print "average Xion:", rad.gas_particles.xion.mean()
        print "max ionization:", rad.gas_particles.xion.max()
    rad.stop()
    
def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-N", dest="N", type="int",default = 1000,
                      help="number of stars [10]")
    result.add_option("-t", unit=units.Myr,
                      dest="t_end", default = 0.1|units.Myr,
                      help="radiation time [%default]")
    result.add_option("-L", unit=units.LSun,
                      dest="Lstar", default = 100|units.LSun,
                      help="luminosity of ionizing source [%default]")
    result.add_option("-p", unit=units.amu/units.cm**3, 
                      dest="rho", default = 1|units.amu/units.cm**3,
                      help="interstellar density [%default] amu/cm^3")
    result.add_option("-d", unit=units.parsec,
                      dest="boxsize", default = 100|units.parsec,
                      help="size of the density box [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)

