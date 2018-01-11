import numpy
from amuse.lab import *
from amuse import datamodel
from amuse.ext.evrard_test import uniform_unit_sphere

def hydro_sink_particles(sinks, gas):
    removed_particles = Particles()
    for s in sinks:
        xs,ys,zs=s.x,s.y,s.z
        radius_squared = s.sink_radius**2
        #print "R=", s.key, numpy.sqrt(radius_squared).in_(units.AU)
        insink=gas.select_array(lambda x,y,z: (x-xs)**2+(y-ys)**2+(z-zs)**2 < radius_squared,['x','y','z'])  
        if len(insink)>0:
            cm=s.position*s.mass
            p=s.velocity*s.mass
            s.mass+=insink.total_mass()
            s.position=(cm+insink.center_of_mass()*insink.total_mass())/s.mass
            s.velocity=(p+insink.total_momentum())/s.mass
            removed_particles.add_particles(insink)
    return removed_particles

def new_sph_particles_from_stellar_wind(stars, mgas):
    new_sph=datamodel.Particles(0)
    for si in stars:
        Ngas = int(-si.Mwind/mgas)
        if Ngas==0:
            continue 
        Mgas = mgas*Ngas
        si.Mwind += Mgas
        si.mass -= Mgas
        add=datamodel.Particles(Ngas)
        add.mass = mgas
        add.h_smooth=0. | units.parsec

        dx,dy,dz=uniform_unit_sphere(Ngas).make_xyz()
        add.x=si.x+(dx * si.wind_radius)
        add.y=si.y+(dy * si.wind_radius)
        add.z=si.z+(dz * si.wind_radius)
        for ri in range(len(add)):
            r = add[ri].position-si.position
            r = r/r.length()
            v_wind = (constants.G*si.mass/(add[ri].position-si.position).length()).sqrt()
            add[ri].u= 0.5 * (v_wind)**2
            add[ri].vx=si.vx + r[0]*si.terminal_wind_velocity
            add[ri].vy=si.vy + r[1]*si.terminal_wind_velocity
            add[ri].vz=si.vz + r[2]*si.terminal_wind_velocity
        new_sph.add_particles(add)  
    return new_sph

def v_terminal_teff(star):
  t4=(numpy.log10(star.temperature.value_in(units.K))-4.).clip(0.,1.)
  return (30 | units.km/units.s) + ((4000 | units.km/units.s)*t4)

def new_hydro_code(hydro_code, dt, converter):
    hydro = hydro_code(converter, redirection="none")
    hydro.parameters.use_hydro_flag=True
    hydro.parameters.timestep=dt
    return hydro

def remove_gas(gas, rremove):
    escaped_gas = Particles()
    if len(gas)>=1:
#        escaped_gas = gas.select_array(lambda r: r.length()>rremove,["position"])
#        escaped_gas = gas[gas.position.length()>rremove].as_set()
        radius_squared = rremove**2
        escaped_gas=gas.select_array(lambda x,y,z: (x**2+y**2+z**2) >= radius_squared,['x','y','z'])  
        
#        print len(escaped_gas), escaped_gas.position.length()
    return escaped_gas

def main(filename):
    stars = Particles(2)
    stars.mass = (1.924785833858, 1.0) | units.MSun  # age=1494.4Myr
    separation = 10|units.AU
    vc = numpy.sqrt(constants.G*stars.mass.sum()/separation)
    stars[0].position = (1, 0, 0) * separation
    stars[0].velocity = (0, 1, 0) * vc
    stars[1].position = (0, 0, 0) | units.AU
    stars[1].velocity = (0, 0, 0) | units.kms
    stars.move_to_center()

    dt = 1|units.day
    stars.dmdt = [-0.11, 0.0] |units.MSun/units.Myr
    stars.temperature = [3309.6, 5672.25] | units.K
    stars.terminal_wind_velocity=v_terminal_teff(stars)
    stars.Mwind = 0.0|units.MSun
    stars.radius = [0.992, 0.00434665] |units.AU
    stars.wind_radius = [1, 0] * stars.radius
    R_BH = 2*constants.G*stars[1].mass/stars[0].terminal_wind_velocity**2
    print "Bondi Hoyle accretion radius:", R_BH
    stars.sink_radius = [0, 1] * R_BH
    #print "stars R=", stars.key, stars.sink_radius
    mgas =  0.1*abs(stars.dmdt.sum()*dt)

    converter=nbody_system.nbody_to_si(1|units.MSun, 1|units.AU)
    gas = Particles(0)
    gas.mass = mgas
    gas.position = (0,0,0)|units.AU
    gas.velocity = (0,0,0)|units.kms
    gas.u = 0 | units.m**2 * units.s**-2 
    gas.h_smooth= 1*units.RSun

    hydro = new_hydro_code(Fi, dt, converter)
    hydro.gas_particles.add_particles(gas)
    hydro.dm_particles.add_particles(stars)
    hydro_to_framework = hydro.gas_particles.new_channel_to(gas, 
        attributes=["x", "y", "z", "vx", "vy", "vz",
                    "mass", "u", "rho", "h_smooth"]) 

    generated_gas = Particles()
    accreted_gas = Particles()
    escaped_gas = Particles()
    dt_snap = 1|units.yr
    t_snap = dt_snap
    while hydro.model_time < 10000|units.yr:
        time = hydro.model_time
        stars.Mwind += stars.dmdt*dt

        accreted = hydro_sink_particles(stars, gas)
        if len(accreted)>0:
            print "N accreted:", time.in_(units.yr), len(accreted), \
                  "m=", accreted.mass.sum().in_(units.MSun)
            accreted_gas.add_particles(accreted.copy())
            gas.remove_particles(accreted)
            hydro.gas_particles.remove_particles(accreted)

        escaped = remove_gas(gas, 5*separation)
        if len(escaped):
            print "N escaped:", len(escaped)
            escaped_gas.add_particles(escaped.copy())
            gas.remove_particles(escaped)
            hydro.gas_particles.remove_particles(escaped)

        new_sph = new_sph_particles_from_stellar_wind(stars, mgas)
        if len(new_sph)>0: 
            #print "N wind:", len(new_sph)
            generated_gas.add_particles(new_sph.copy())
            gas.add_particles(new_sph)
            hydro.gas_particles.add_particles(new_sph)
            
        if len(gas)>100:
            #gas.synchronize_to(hydro.gas_particles)
            hydro.evolve_model(hydro.model_time+dt)
            hydro_to_framework.copy()
            if time>t_snap:
                t_snap += dt_snap
                write_set_to_file(gas, filename, 'hdf5')
                print "time=", hydro.model_time, "Ngas=", len(gas), \
                      mgas*len(gas)
                print "T=", time, "M=", stars[0].mass, stars[1].mass
                print "Gas = ", len(generated_gas), len(accreted_gas), \
                      len(escaped_gas)
        
    hydro.stop()

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser()
    result.add_option("-f", dest="filename", default = "hydro_outflow.h5",
                      help="output filename [hydro.hdf5]")
    return result

if __name__ in ('__main__', '__plot__'):
    o, arguments  = new_option_parser().parse_args()
    set_printing_strategy("custom", #nbody_converter = converter, 
                          preferred_units = [units.MSun, units.AU, units.Myr], 
                          precision = 14, prefix = "", separator = " [", suffix = "]")
    main(o.filename)

