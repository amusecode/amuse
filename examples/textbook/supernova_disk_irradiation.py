import numpy
from amuse.lab import *
from amuse.ext.protodisk import ProtoPlanetaryDisk
from matplotlib import pyplot
from supernova_IIp_Lightcurve import Supernova_IIp
    
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
        raise Exception(
            "Error in calculating mu: mass fractions do not sum to 1.0")
    return constants.proton_mass / (X*(1.0+x_ion) + Y*(1.0+2.0*x_ion)/4.0
                                     + Z*x_ion/2.0)

def update_source_particle(rad, time, old_source, efficiency_factor,
                           supernova_IIp):
    source=Particle()
    source.position = old_source.position
    source.luminosity = efficiency_factor**2 \
                         * supernova_IIp.luminosity_at_time(time)/(20.|units.eV)
    source.flux = source.luminosity
    source.xion = old_source.xion
    source.u = old_source.u

    rad.src_particles.remove_particle(old_source)
    rad.src_particles.add_particle(source)
    return source

def print_diagnostics(time, supernova, disk):
    umin = disk.u.min()
    umean = disk.u.mean()
    umax = disk.u.max()
    Tmin =  mu() / constants.kB * umax
    Tmean =  mu() / constants.kB * umean
    Tmax =  mu() / constants.kB * umin

    print("Time=", time.in_(units.day))
    print("Supernova luminosity:", \
          (supernova.luminosity*(20.|units.eV)).in_(units.LSun))
    print("Ionization:", disk.xion.min(), disk.xion.mean(), disk.xion.max())
    print("Intenal energy:", umin, umean, umax)
    print("Temperature:", Tmin, Tmean, Tmax)
    print("Density:", disk.density.min().in_(units.amu/units.cm**3), \
          disk.density.mean().in_(units.amu/units.cm**3), \
          disk.density.max().in_(units.amu/units.cm**3))
    print("scaleheight:", abs(disk.z.value_in(units.AU)).mean())

def main(Ndisk, Mstar, Mdisk, Rin, Rout, t_end, Nray, x, y, z):

    time = 0 | units.Myr
    supernova_IIp = Supernova_IIp(10|units.day)
    
    efficiency_factor = 0.1
    Rsn = efficiency_factor * (x**2 + y**2 + z**2)**0.5

    supernova = Particle()
    supernova.position = (x.value_in(units.parsec),
                          y.value_in(units.parsec),
                          z.value_in(units.parsec)) |units.parsec
    supernova.position *= efficiency_factor
    supernova_IIp.particles.add_particle(supernova)
    supernova_IIp.evolve_model(time)
    supernova.luminosity = efficiency_factor**2 * supernova.luminosity
    supernova.xion = 0.0
    supernova.u = (10**51 | units.erg)/(10|units.MSun)
    
    stellar = SeBa()
    star = Particle()
    star.mass = Mstar
    star.position = (0,0,0) | units.AU
    star.velocity = (0,0,0) | units.kms
    stellar.particles.add_particle(star)
    stellar.evolve_model(1|units.Myr)
    star.luminosity = stellar.particles[0].luminosity/(20. | units.eV)
    star.temperature = stellar.particles[0].temperature
    stellar.stop()
    star.u = (9. |units.kms)**2
    star.xion = 0.0
    print(star)

    print("M=", Mdisk/Mstar)
    converter=nbody_system.nbody_to_si(Mstar, 1 | units.AU)
    disk = ProtoPlanetaryDisk(Ndisk, convert_nbody=converter,
                              Rmin=Rin.value_in(units.AU), 
                              Rmax=Rout.value_in(units.AU),
                              q_out=25.0, discfraction=Mdisk/Mstar).result
####                              q_out=2.0, discfraction=Mdisk/Mstar).result
    print(disk.x.max().in_(units.AU))
    print(disk.mass.sum().in_(units.MSun))
    print(disk.u.max().in_(units.kms**2))
    print(disk.mass.min().in_(units.MSun))
    disk.flux = 0. | units.s**-1
    disk.xion = 0.0

    dt = t_end/1000.
    
    hydro = Fi(converter)
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.radiation_flag=False
    hydro.parameters.self_gravity_flag=True
    
    hydro.parameters.gamma=1.
    hydro.parameters.isothermal_flag=True

# Non-isothermal: lowest temperature remains constat at 20K
#    hydro.parameters.gamma=1.66667
#    hydro.parameters.isothermal_flag=False
    
    hydro.parameters.integrate_entropy_flag=False
    hydro.parameters.timestep=0.5 | units.hour #0.125 | units.day
#    hydro.parameters.verbosity=99
    hydro.parameters.epsilon_squared=0.1 | units.AU**2
    hydro.parameters.courant=0.2    
    hydro.parameters.artificial_viscosity_alpha = 0.1 
    
    hydro.gas_particles.add_particles(disk)
    hydro.dm_particles.add_particle(star)
    hydro_to_disk = hydro.gas_particles.new_channel_to(disk)
    hydro_to_star = hydro.dm_particles.new_channel_to(star.as_set())
    disk_to_hydro = disk.new_channel_to(hydro.gas_particles)
    star_to_hydro = star.as_set().new_channel_to(hydro.dm_particles)
    hydro.evolve_model(1|units.hour)
    hydro_to_disk.copy()
    hydro_to_star.copy()
    
    radiative = SPHRay(redirection="file",
                       number_of_workers=4)#, debugger="gdb")
    radiative.parameters.number_of_rays=Nray/dt
    print(dt.in_(units.yr))
    radiative.parameters.default_spectral_type=-3.
    radiative.parameters.box_size=10000. | units.AU
    radiative.parameters.ionization_temperature_solver=2
    print(radiative.parameters)

#    radiative.src_particles.add_particle(star)
    radiative.src_particles.add_particle(supernova)
    radiative.gas_particles.add_particles(disk)

    gas_to_rad = disk.new_channel_to(radiative.gas_particles)
    rad_to_gas = radiative.gas_particles.new_channel_to(disk)
    
    print("Before")
    print("Luminosity:", radiative.src_particles.luminosity)
    print("min ionization:", radiative.gas_particles.xion.min())
    print("average Xion:", radiative.gas_particles.xion.mean())
    print("max ionization:", radiative.gas_particles.xion.max())
    print("min u:", radiative.gas_particles.u.min())
    print("average u:", radiative.gas_particles.u.mean())
    print("max u:", radiative.gas_particles.u.max())

    Tmean = [] | units.K
    Tmin = [] | units.K
    Tmax = [] | units.K
    t = [] | units.day
    while radiative.model_time<t_end:

        supernova = update_source_particle(radiative, time+0.5*dt, supernova,
                                           efficiency_factor, supernova_IIp)

        radiative.evolve_model(time+0.5*dt)
        print("RT done at time:", time.in_(units.day))
        rad_to_gas.copy()

        disk_to_hydro.copy()
        star_to_hydro.copy()
        hydro.evolve_model(time + dt)
        hydro_to_disk.copy()
        hydro_to_star.copy()

        supernova = update_source_particle(radiative, time+dt, supernova,
                                           efficiency_factor, supernova_IIp)

        radiative.evolve_model(time+dt)
        print("RT done at time:", time.in_(units.day))
        rad_to_gas.copy()

        time += dt
        print_diagnostics(time, supernova, disk)
        
        Temperature =  mu() / constants.kB * disk.u
        t.append(time)
        Tmean.append(Temperature.mean())
        Tmin.append(Temperature.min())
        Tmax.append(Temperature.max())

        #write_set_to_file(disk, "disk_irradiation.amuse", "amuse")
        
        print("timescale:", (disk.mass.sum().value_in(units.amu) \
                              / ((Rout/Rsn)**2*supernova.luminosity)).in_(units.yr))
        print("scaleheight:", abs(disk.z.value_in(units.AU)).mean())
        
        #pyplot.hist2d(abs(disk.x.value_in(units.AU)), abs(numpy.log10(Temperature.value_in(units.K))), bins=200)
       # pyplot.hist2d(abs(disk.x.value_in(units.AU)), abs(disk.z.value_in(units.AU)), bins=200)
        #,norm=LogNorm())
       #, Temperature.in_(units.K))
#       pyplot.tripcolor(abs(disk.x.value_in(units.AU)), abs(disk.y.value_in(units.AU)), Temperature.in_(units.K))
        #pyplot.hist(abs(disk.x.value_in(units.AU)), 100)
       # pyplot.show()

    print_diagnostics(time, supernova, disk)

    #plot_ionization_fraction(disk.position, disk.xion)
#    plot_ionization_fraction(disk.z, disk.u.value_in(units.kms**2))
#    plot_ionization_fraction(disk.x, disk.u.value_in(units.kms**2))
    radiative.stop()
    plot_temperature(t, Tmin, Tmean, Tmax)

#~ from prepare_figure import single_frame, figure_frame
#~ from distinct_colours import get_distinct

def plot_temperature(t, tmin, tmean, tmax):

    x_label = "t [day]"
    y_label = 'T [K]'
    figure = single_frame(x_label, y_label, logx=False, logy=False,
                          xsize=14, ysize=8)
    pyplot.plot(t.value_in(units.day), tmean.value_in(units.K), c='k')
    pyplot.plot(t.value_in(units.day), tmin.value_in(units.K), c='r')
    pyplot.plot(t.value_in(units.day), tmax.value_in(units.K), c='b')
    pyplot.show()

def plot_ionization_fraction(pos, xion):
    r = []
    x = []
    for pi, xi in zip(pos, xion):
        #r.append(pi.length())
        r.append(numpy.log10(pi.value_in(units.AU)+0.000001))
        r.append(pi.value_in(units.AU))
        x.append(numpy.log10(xi+0.000001))
    r, x = list(zip(*sorted(zip(r, x))))
    
    from matplotlib import pyplot
    x_label = "r [pc]"
    y_label = r'$\xi_{\rm ion}$'
    figure = single_frame(x_label, y_label, logx=False, logy=False,
                          xsize=14, ysize=8)
    pyplot.scatter(r, x, c=get_distinct(1), lw=0, s=100)
#    pyplot.xlim(-1, 1)
#    pyplot.ylim(-0.04, 1.19)
    #pyplot.savefig("fig_ionization_of_GMC")
    pyplot.show()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("--Ndisk", dest="Ndisk", type="int", default = 10000,
                      help="number of disk particles [%default]")
    result.add_option("-t", unit=units.yr,
                      dest="t_end", type="float", default = 2|units.yr,
                      help="radiation time [%default]")
    result.add_option("-r", unit=units.AU, type="float",
                      dest="Rin", default = 1.|units.AU,
                      help="inner disk radius [%default]")
    result.add_option("-R", unit=units.AU, type="float",
                      dest="Rout", default = 100.|units.AU,
                      help="outer disk radius [%default]")
    result.add_option("--Mstar", unit=units.MSun, type="float",
                      dest="Mstar", default = 1|units.MSun,
                      help="stellar mass")
    result.add_option("--Mdisk", unit=units.MSun, type="float",
                      dest="Mdisk", default = 0.01|units.MSun,
                      help="disk mass")
    result.add_option("-x", unit=units.parsec, type="float",
                      dest="x", default = 0.1|units.parsec,
                      help="supnova x-position")
    result.add_option("-y", unit=units.parsec, type="float",
                      dest="y", default = 0|units.parsec,
                      help="supnova y-position")
    result.add_option("-z", unit=units.parsec, type="float",
                      dest="z", default = 0.01|units.parsec,
                      help="supnova z-position")
    result.add_option("--Nray", type="int",
                      dest="Nray", default = 10**7,
                      help="number of rays [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    main(**o.__dict__)


