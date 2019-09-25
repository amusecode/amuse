import numpy
from amuse.lab import *
from amuse import datamodel
from amuse.ext.evrard_test import uniform_unit_sphere

set_printing_strategy("custom", #nbody_converter = converter, 
                      preferred_units = [units.MSun, units.AU, units.Myr], 
                      precision = 5, prefix = "", separator = " [", suffix = "]"
)

def new_sph_particles_from_stellar_wind(stars, mgas):
    new_sph=datamodel.Particles(0)
    for si in stars:
        Ngas = int(-si.Mwind/mgas)
        if Ngas==0:
            continue 
        Mgas = mgas*Ngas
        si.Mwind += Mgas
        add=datamodel.Particles(Ngas)
        add.mass = mgas
        add.h_smooth=0. | units.parsec

        dx,dy,dz=uniform_unit_sphere(Ngas).make_xyz()
        add.x=si.x+(dx * si.radius)
        add.y=si.y+(dy * si.radius)
        add.z=si.z+(dz * si.radius)
        for ri in range(len(add)):
            r = add[ri].position-si.position
            r = r/r.length()
            v_wind = (constants.G*si.mass/(add[ri].position-si.position).length()).sqrt()
            add.u= 0.5 * (v_wind)**2
            add.vx=si.vx + r[0]*si.terminal_wind_velocity
            add.vy=si.vy + r[1]*si.terminal_wind_velocity
            add.vz=si.vz + r[2]*si.terminal_wind_velocity
        new_sph.add_particles(add)  
    return new_sph

def v_terminal_teff(star):
  t4=(numpy.log10(star.temperature.value_in(units.K))-4.).clip(0.,1.)
  return (30 | units.km/units.s) + ((4000 | units.km/units.s)*t4)

def main():
    stars = Particles(2)
    stars.mass = (9.5, 10) | units.MSun
    stars[0].position = (1, 0, 0) | units.AU
    stars[0].velocity = (0, 0, 0) | units.kms
    stars[1].position = (0, 0, 0) | units.AU
    stars[1].velocity = (0, 0, 0) | units.kms
    stars.move_to_center()

    a = stars.position.length().amax()
    vc = constants.G*stars.mass.sum()/a
    stellar = SeBa()
    stellar.particles.add_particles(stars)
    stellar_to_framework = stellar.particles.new_channel_to(stars)
    stellar.evolve_model(26|units.Myr)
    stellar_to_framework.copy_attributes(["mass","radius","temperature"])
    dt = 0.1|units.Myr
    stellar.evolve_model((26|units.Myr)+dt)
    stars.dmdt = (stellar.particles.mass-stars.mass)/dt
    stars.Mwind = 0 | units.MSun
    stars.terminal_wind_velocity=v_terminal_teff(stars)
    stellar.stop()
    dt = 0.1|units.day
    mgas =  0.1*abs(stars.dmdt.sum()*dt)

    converter=nbody_system.nbody_to_si(1|units.MSun, a)
    bodies = Particles(0)
    bodies.mass = mgas
    bodies.position = (0,0,0)|units.AU
    bodies.velocity = (0,0,0)|units.kms
    bodies.u = 0 | units.m**2 * units.s**-2 
    bodies.h_smooth= 0.01*a

    hydro = Fi(converter, redirection="none")
    if len(bodies)>0:
        hydro.gas_particles.add_particles(bodies)
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.timestep=dt
    hydro.parameters.periodic_box_size = 1000*a
    hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)

    moving_bodies = ParticlesSuperset([stars, bodies])
    filename = "hydro_outflow.hdf5"
    istep = 0
    while hydro.model_time < 10|units.yr:
        stars.Mwind += stars.dmdt*dt
        new_sph = new_sph_particles_from_stellar_wind(stars, mgas)

        if len(new_sph)>0: 
            bodies.add_particles(new_sph)
            bodies.synchronize_to(hydro.gas_particles)
        print "time=", hydro.model_time, "Ngas=", len(bodies), mgas*len(bodies)
        if len(bodies)>100:
            hydro.evolve_model(hydro.model_time+dt)
            hydro_to_framework.copy()
            if istep%10 == 0:
                write_set_to_file(moving_bodies, filename, 'hdf5')
            istep += 1
    hydro.stop()

if __name__ in ('__main__', '__plot__'):
    main()

