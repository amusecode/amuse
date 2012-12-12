"""
   example code for bridging a gravity solver with a hydrodynamics solver
"""
import numpy
from amuse.lab import *
from amuse.couple import bridge
from amuse import datamodel
from amuse.ext.evrard_test import uniform_unit_sphere

def new_sph_particles_from_stellar_wind(stars, mgas):

    new_sph=datamodel.Particles(0)
    for si in stars:
        p = si.position
        v = si.velocity
        Ngas = int(-si.Mwind/mgas)
        print "new Ngas=", si.mass, Ngas,
        if Ngas==0:
          continue 
        Mgas = mgas*Ngas
        si.Mwind += Mgas
#        Ngas = 10
#        mgas = Mgas/10. 
        print "new Ngas=", Ngas, mgas
        add=datamodel.Particles(Ngas)
        add.mass = mgas
        add.h_smooth=0. | units.parsec
        dx,dy,dz=uniform_unit_sphere(Ngas).make_xyz()
        
        add.x=si.x+(dx * si.radius)
        add.y=si.y+(dy * si.radius)
        add.z=si.z+(dz * si.radius)
        for ri in range(len(add)):
          r = add[ri].position-p
          r = r/r.length()
          v_wind = (constants.G*si.mass/(add[ri].position-p).length()).sqrt()
          add.u= 0.5 * (v_wind)**2
          v_wind = si.terminal_wind_velocity
          add.vx=v.x + r[0]*v_wind
          add.vy=v.y + r[1]*v_wind
          add.vz=v.z + r[2]*v_wind
        new_sph.add_particles(add)  
    return new_sph

def v_terminal_teff(star):
  t4=numpy.log10(star.temperature.value_in(units.K))-4.
  t4=t4.clip(0.,1.)
  return (30 | units.km/units.s) + ((4000 | units.km/units.s)*t4)

def get_kepler_elements(model_time, bh, star, converter):
    kep = Kepler(converter)
    kep.initialize_code()
    pos = bh.position - star.position
    vel = bh.velocity - star.velocity
    print "Kep:", bh.mass + star.mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]
    kep.initialize_from_dyn(bh.mass + star.mass, pos[0], pos[1], pos[2], vel[0], vel[1], vel[2])
    a,e = kep.get_elements()
    kep.stop()
    return a, e

def gravity_hydro_bridge(a, ecc, t_end, n_steps, Rgas, Mgas, Ngas):

    stars = Particles(3)
    stars.mass = [5.0, 9.9, 10.0] | units.MSun
    stellar = SeBa()
    stellar.particles.add_particles(stars)
    stellar_to_framework = stellar.particles.new_channel_to(stars)
    stellar.evolve_model(26|units.Myr)
    stellar_to_framework.copy_attributes(["mass","radius","temperature"])
    print "stars=", stars
    stellar.evolve_model(26.1|units.Myr)
    stars.dmdt = (stellar.particles.mass-stars.mass)/(0.1|units.Myr)
    stars.Mwind = 0 | units.MSun
    stars.terminal_wind_velocity=v_terminal_teff(stars)
    stellar.stop()
    print "dmdt=", stars.dmdt
    dt = 0.1|units.day
    mgas =  0.1*abs(stars.dmdt.sum()*dt)
    print "mgas=", mgas, stars.dmdt/mgas

    vc = constants.G*stars.mass.sum()/a
    Porb = 2*numpy.pi*(a**3/(constants.G*stars.mass.sum())).sqrt()
    stars[0].position = (0,0,0) | units.AU
    stars[0].velocity = (0,0,0) | units.kms
    vc = (constants.G*stars[:2].mass.sum()/(a*(1+ecc))).sqrt()
    vc *= numpy.sqrt((1-ecc)/(1+ecc)) 
    stars[1].position = (a.value_in(units.AU),0,0) | units.AU
    stars[1].velocity = (0,vc.value_in(units.kms),0) | units.kms
    stars[:2].move_to_center()
    ecc = 0.2
    vc = (constants.G*stars.mass.sum()/(10*a*(1+ecc))).sqrt()
    vc *= numpy.sqrt((1-ecc)/(1+ecc)) 
    stars[2].position = (10*a.value_in(units.AU),0,0) | units.AU
    stars[2].velocity = (0,vc.value_in(units.kms),0) | units.kms
    stars.move_to_center()
    stars.radius = 0.2*a

    #define for printing
#    stars.h_smooth= 0.0*a
#    stars.u = 0 | units.kms**2

    converter=nbody_system.nbody_to_si(stars.mass.sum(), a)
    gravity = ph4(converter, redirection="none")
    gravity.particles.add_particles(stars)
    gravity.parameters.epsilon_squared = (10|units.RSun)**2
    Ed0_tot = gravity.kinetic_energy + gravity.potential_energy

    channel_from_gravity = gravity.particles.new_channel_to(stars)
    channel_from_to_gravity = stars.new_channel_to(gravity.particles)

    dt = t_end/float(n_steps)
    converter=nbody_system.nbody_to_si(1.0|units.MSun, a)
    ism = Particles(0)
    ism.mass = mgas
    ism.position = (0,0,0)|units.AU
    ism.velocity = (0,0,0)|units.kms
    ism.u = 0 | units.m**2 * units.s**-2 
    ism.h_smooth= 0.01*a

    hydro = Fi(converter, redirection="none")
    hydro.parameters.timestep = dt/8.
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.radiation_flag=False
    hydro.parameters.self_gravity_flag=True
    hydro.parameters.integrate_entropy_flag=False
    hydro.parameters.gamma=1.
    hydro.parameters.isothermal_flag=True
    hydro.parameters.epsilon_squared = (10|units.RSun)**2
    if len(ism)>0:
        hydro.gas_particles.add_particles(ism)
    Eh0_tot = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
    hydro.parameters.periodic_box_size = 10000*a


    channel_from_hydro = hydro.gas_particles.new_channel_to(ism)
    channel_from_to_hydro = ism.new_channel_to(hydro.gas_particles)

    moving_bodies = ParticlesSuperset([stars, ism])

    model_time = 0 | units.Myr
    filename = "stellargravhydro.hdf5"
    if len(ism)>0:
        write_set_to_file(moving_bodies, filename, 'hdf5')

    gravhydro = bridge.Bridge(use_threading=False)
    gravhydro.add_system(gravity, (hydro,) )
    gravhydro.add_system(hydro, (gravity,) )
    gravhydro.timestep = min(dt, 2*hydro.parameters.timestep)

    istep = 0
    while model_time < t_end:
        model_time += dt
        a, e = get_kepler_elements(gravity.model_time, stars[0], stars[1], converter) 
        print "AB: time=", model_time, a, e
        com_star = Particles(1)
        com_star.mass = stars[:2].mass.sum()
        com_star.position = stars[:2].center_of_mass()
        com_star.velocity = stars[:2].center_of_mass_velocity()
        a, e = get_kepler_elements(gravity.model_time, com_star[0], stars[2], converter) 
        print "(AB)C: time=", model_time, a, e

        stars.Mwind += stars.dmdt*dt
        print "Mw=", stars.Mwind, stars.Mwind/mgas
        new_sph = new_sph_particles_from_stellar_wind(stars, mgas)
        print "Ngas=", len(new_sph), len(ism), len(hydro.gas_particles)

        if len(new_sph)>0: # and len(bodies)<4000:
            ism.add_particles(new_sph)
            ism.synchronize_to(hydro.gas_particles)

        if len(ism)>100:
            print "t=", hydro.model_time, dt
            gravhydro.evolve_model(model_time)
            channel_from_gravity.copy()
            channel_from_hydro.copy()
            channel_from_hydro.copy_attributes(["u"])
            print "N=", len(hydro.particles)
            Ed_tot = gravity.kinetic_energy + gravity.potential_energy
            Eh_tot = hydro.kinetic_energy + hydro.potential_energy + hydro.thermal_energy
            print "Energies:", Ed_tot/Ed0_tot, Eh_tot/Eh0_tot

            if istep%10==0:
                write_set_to_file(moving_bodies, filename, 'hdf5')
            istep+=1

    gravity.stop()
    hydro.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-n", dest="n_steps", type="int", default = 1000,
                      help="number of diagnostics time steps [%default]")
    result.add_option("-N", dest="Ngas", type="int", default = 1024,
                      help="number of gas particles [%default]")
    result.add_option("-M", unit=units.MSun,
                      dest="Mgas", type="float", default = 1|units.MSun,
                      help="Mass of the gas [%default]")
    result.add_option("-R", unit=units.AU,
                      dest="Rgas", type="float", default = 1|units.AU,
                      help="Size of the gas distribution [%default]")
    result.add_option("-a", unit=units.AU,
                      dest="a", type="float", default = 0.2|units.AU,
                      help="initial orbital separation [%default]")
    result.add_option("-e", dest="ecc", type="float", default = 0.0,
                      help="initial orbital eccentricity [%default]")
    result.add_option("-t", unit=units.yr, 
                      dest="t_end", type="float", default = 10|units.yr,
                      help="end time of the simulation [%default]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    gravity_hydro_bridge(**o.__dict__)
