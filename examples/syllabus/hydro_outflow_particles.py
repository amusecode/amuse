import os.path
import math
import numpy
from amuse.lab import *
from optparse import OptionParser
import time
from amuse import datamodel

from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.units.optparse import OptionParser

set_printing_strategy("custom", #nbody_converter = converter, 
                      preferred_units = [units.MSun, units.AU, units.Myr], 
                      precision = 5, prefix = "", separator = " [", suffix = "]"
)

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
        Ngas = 10
        mgas = Mgas/10. 
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
          """
          add.u= 0.5 * (si.terminal_wind_velocity)**2
          v_wind = (constants.G*si.mass/(add[ri].position-p).length()).sqrt()
          add.vx=v.x + r[0]*v_wind
          add.vy=v.y + r[1]*v_wind
          add.vz=v.z + r[2]*v_wind
          """
#          add.u= 0.5 * (si.terminal_wind_velocity)**2
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

def main():

    stars = Particles(2)
    stars.mass = (9.5, 10) | units.MSun
#    stars.radius = (100,100)|units.RSun
    stars[0].position = (1, 0, 0) | units.AU
    stars[0].velocity = (0, 0, 0) | units.kms
    stars[1].position = (0, 0, 0) | units.AU
    stars[1].velocity = (0, 0, 0) | units.kms
    stars.move_to_center()

    a = stars.position.length().amax()
    print "a=", a
    vc = constants.G*stars.mass.sum()/a
    stellar = SeBa()
    stellar.particles.add_particles(stars)
    stellar_to_framework = stellar.particles.new_channel_to(stars)
    stellar.evolve_model(26|units.Myr)
    stellar_to_framework.copy_attributes(["mass","radius","temperature"])
    print "stars=", stars
    dt = 0.1|units.Myr
    stellar.evolve_model((26|units.Myr)+dt)
    stars.dmdt = (stellar.particles.mass-stars.mass)/dt
    stars.Mwind = 0 | units.MSun
    stars.terminal_wind_velocity=v_terminal_teff(stars)
    stellar.stop()
    print "dmdt=", stars.dmdt
    dt = 0.1|units.day
    mgas =  0.01*abs(stars.dmdt.sum()*dt)
    print "mgas=", mgas, stars.dmdt/mgas

#    converter=nbody_system.nbody_to_si(0.001 | units.MSun, 0.01*a)
    converter=nbody_system.nbody_to_si(1|units.MSun, a)
    a2 = (10*a)**2

    """
    bodies = new_plummer_gas_model(200, convert_nbody=converter)
    print "N=", len(bodies),
    bodies=bodies.select_array(lambda x,y,z: (x)**2+(y)**2+(z)**2 < a2,['x','y','z'])  
    print len(bodies)
    """

    bodies = Particles(0)
    bodies.mass = mgas
    bodies.position = (0,0,0)|units.AU
    bodies.velocity = (0,0,0)|units.kms
    bodies.u = 0 | units.m**2 * units.s**-2 
    bodies.h_smooth= 0.01*a

#    hydro = Gadget2(converter)
    hydro = Fi(converter, redirection="none")
    if len(bodies)>0:
        hydro.gas_particles.add_particles(bodies)
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.radiation_flag=False
    hydro.parameters.self_gravity_flag=True
    hydro.parameters.integrate_entropy_flag=False
#    hydro.parameters.time_max = 100|units.day
    print "parameters=", hydro.parameters
    hydro.parameters.gamma=1.
    hydro.parameters.isothermal_flag=True
    hydro.parameters.timestep=dt
    hydro.parameters.periodic_box_size = 1000*a
    hydro_to_framework = hydro.gas_particles.new_channel_to(bodies)

    moving_bodies = ParticlesSuperset([stars, bodies])

    filename = "hydro_outflow.hdf5"
    model_time = 0 |units.yr
    istep = 0
    while model_time < 10|units.yr:
        stars.Mwind += stars.dmdt*dt
        print "Mw=", stars.Mwind, stars.Mwind/mgas
        new_sph = new_sph_particles_from_stellar_wind(stars, mgas)
        print "Ngas=", len(new_sph), len(bodies), len(hydro.gas_particles)
        if len(new_sph)>0: # and len(bodies)<4000:
#            new_sph=new_sph.select_array(lambda x,y,z: (x)**2+(y)**2+(z)**2 < a2,['x','y','z'])  
            bodies.add_particles(new_sph)
            bodies.synchronize_to(hydro.gas_particles)

        if len(bodies)>100:
            print "t=", hydro.model_time, dt
            hydro.evolve_model(hydro.model_time+dt)
            hydro_to_framework.copy()
            print "N=", len(hydro.particles)
            if istep%10 == 0:
                write_set_to_file(moving_bodies, filename, 'hdf5')
            istep += 1
#            bodies=bodies.select_array(lambda x,y,z: (x)**2+(y)**2+(z)**2 < a2,['x','y','z'])  
#            bodies.synchronize_to(hydro.gas_particles)
       
    hydro.stop()

if __name__ in ('__main__', '__plot__'):
    main()

