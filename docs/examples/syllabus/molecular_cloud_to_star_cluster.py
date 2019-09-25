"""
  example of molecular cloud evolution with explictly 
  split SPH and grav evolution

  Initial condition is a smooth spherical cloud with random velocities
  as in Bonnell et al. (2003)  
  
"""  

import numpy
  
from matplotlib import pyplot 

from amuse.lab import *
from amuse.units import nbody_system
from amuse.units import units
from amuse.units import constants

from amuse.community.fi.interface import Fi

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube

from amuse.ext.derived_grav_systems import copycat
from amuse.ext.bridge import bridge

from amuse.io import write_set_to_file
from amuse.support.data import ParticlesWithUnitsConverted
from amuse.datamodel import Particles

def make_map(sph,N=100,L=1):

    x,y=numpy.indices( ( N+1,N+1 ))

    x=L*(x.flatten()-N/2.)/N
    y=L*(y.flatten()-N/2.)/N
    z=x*0.
    vx=0.*x
    vy=0.*x
    vz=0.*x

    x=units.parsec(x)
    y=units.parsec(y)
    z=units.parsec(z)
    vx=units.kms(vx)
    vy=units.kms(vy)
    vz=units.kms(vz)

    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rho.reshape((N+1,N+1))

    return rho
    
'''
def write_output(filename, parts, conv):
    
    output= file(filename, 'w')
    for i in range (0,len(parts)):
        #print i
        rho = conv.to_nbody(parts.rho[i])
        mass= conv.to_nbody(parts.mass[i])
        x= conv.to_nbody(parts.x[i])
        y= conv.to_nbody(parts.y[i])
        z= conv.to_nbody(parts.z[i])
        vx= conv.to_nbody(parts.vx[i])
        vy= conv.to_nbody(parts.vy[i])
        vz= conv.to_nbody(parts.vz[i])
        print>> output, rho.value_in(nbody_system.mass/nbody_system.length**3), mass.value_in(nbody_system.mass), x.value_in(nbody_system.length), y.value_in(nbody_system.length),  z.value_in(nbody_system.length), vx.value_in(nbody_system.length/nbody_system.time), vy.value_in(nbody_system.length/nbody_system.time),vz.value_in(nbody_system.length/nbody_system.time)
    output.close()
    return 0
'''

def write_output(filename, parts, conv):

    particles_nbody = ParticlesWithUnitsConverted(parts, conv.as_converter_from_nbody_to_si())
    #print  particles_nbody

    write_set_to_file(particles_nbody, filename, "txt", attribute_names= ('rho', 'mass', 'x', 'y', 'z','vx', 'vy', 'vz'))

    return 0

def plot_stars(time, stars, i, L=6.):
    fig=pyplot.figure(figsize=(12,12))
    m =  100.0*stars.mass/max(stars.mass)
    x = -stars.x.value_in(units.parsec)
    y = stars.y.value_in(units.parsec)
    pyplot.scatter(x, y, s=m)
    pyplot.title("Star cluster at"+time.as_string_in(units.Myr))
    pyplot.xlim(-L/2., L/2.)
    pyplot.ylim(-L/2., L/2.)
    pyplot.xlabel("x [pc]")
    pyplot.ylabel("y [pc]")
    pyplot.savefig("SC_"+str(i)+".png")

def plot_hydro(time, sph, i, L=10):
    fig=pyplot.figure(figsize=(12,12))
    rho=make_map(sph,N=200,L=L)
    pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)), extent=[-L/2,L/2,-L/2,L/2],vmin=1,vmax=5)
#    subplot.set_title("GMC at zero age")
    pyplot.title("Molecular cloud at time="+time.as_string_in(units.Myr))
    pyplot.xlabel("x [pc]")
    pyplot.ylabel("x [pc]")
    pyplot.title("GMC at time="+time.as_string_in(units.Myr))
    pyplot.savefig("GMC_"+str(i)+".png")

def plot_hydro_and_stars(time, sph, L=10):
    fig=pyplot.figure(figsize=(12,12))
    rho=make_map(sph,N=200,L=L)
    pyplot.imshow(numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)), extent=[-L/2,L/2,-L/2,L/2],vmin=1,vmax=5)
#    subplot.set_title("GMC at zero age")

    stars = get_stars_from_molecular_clous(sph.gas_particles)
    m =  100.0*stars.mass/max(stars.mass)
    x = -stars.x.value_in(units.parsec)
    y = stars.y.value_in(units.parsec)
    pyplot.scatter(x, y, s=m)
    pyplot.xlim(-L/2., L/2.)
    pyplot.ylim(-L/2., L/2.)
    pyplot.title("Molecular cloud at time="+time.as_string_in(units.Myr))
    pyplot.xlabel("x [pc]")
    pyplot.ylabel("x [pc]")
    pyplot.title("GMC at time="+time.as_string_in(units.Myr))
    pyplot.savefig("GMC_SC.png")

def run_molecular_cloud(N=100, Mcloud=100. | units.MSun, Rcloud=1. | units.parsec):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)

    rho_cloud = 3.*Mcloud/(4.*numpy.pi*Rcloud**3)
    print rho_cloud
    tff = 0.5427/numpy.sqrt(constants.G*rho_cloud)
    print "t_ff=", tff.value_in(units.Myr), 'Myr'

    dt = 5.e-2 | units.Myr
    tend=1.0 | units.Myr
#    tend=2.0 | units.Myr

    parts=molecular_cloud(targetN=N,convert_nbody=conv,
            base_grid=body_centered_grid_unit_cube, seed=100).result

    sph=Fi(conv, number_of_workers=3)
    sph.parameters.use_hydro_flag=True
    sph.parameters.radiation_flag=False
    sph.parameters.gamma=1
    sph.parameters.isothermal_flag=True
    sph.parameters.integrate_entropy_flag=False
    sph.parameters.timestep=dt  
    sph.parameters.verbosity = 0
    sph.parameters.eps_is_h_flag = False# h_smooth is constant
    eps = 0.1 | units.parsec
    sph.parameters.gas_epsilon = eps
    sph.parameters.sph_h_const = eps
    parts.h_smooth= eps

    print 'eps-h flag', sph.get_eps_is_h(), sph.get_consthsm()

    expected_dt = 0.2*numpy.pi*numpy.power(eps, 1.5)/numpy.sqrt(constants.G*Mcloud/N)

    print "dt_exp=", expected_dt.value_in(units.Myr)
    print "dt=", dt
    print "eps=", sph.parameters.gas_epsilon.in_(units.parsec)

    sph.gas_particles.add_particles(parts)

    #grav=copycat(Fi, sph, conv)
    #sys=bridge(verbose=False)
    #sys.add_system(sph,(grav,),False)
    channel_from_sph_to_parts= sph.gas_particles.new_channel_to(parts)
    channel_from_parts_to_sph= parts.new_channel_to(sph.gas_particles)

    i=0
    L=6
    E0 = 0.0
    ttarget = 0.0 | units.Myr

    plot_hydro(ttarget, sph, i, L)

    while ttarget < tend:
        ttarget=float(i)*dt
        print ttarget
        sph.evolve_model(ttarget, timestep=dt)
        E = sph.gas_particles.kinetic_energy()+sph.gas_particles.potential_energy() + sph.gas_particles.thermal_energy()
        E_th = sph.gas_particles.thermal_energy()
        if i==0:
            E0 = E
        Eerr = (E-E0)/E0
        print 'energy=', E, 'energy_error=', Eerr, 'e_th=', E_th
        channel_from_sph_to_parts.copy()
        """
        filename = 'm400k_r10pc_e01_'+ str(i).zfill(2) + '.dat'
        print filename
        parts_sorted = parts.sorted_by_attribute('rho')
        write_output(filename, parts_sorted, conv)        
        """
        plot_hydro(ttarget, sph, i, L)
        i=i+1

    plot_hydro_and_stars(ttarget, sph, L)
    
    sph.stop()
    return parts

def make_stars(cluster_particle):
    sfe = 0.3
    mmean = 1.0|units.MSun
    N = int(sfe*cluster_particle.mass/mmean)
    stars = Particles(0)
    print "N_cluster=", N
    if N>0:
        masses = new_salpeter_mass_distribution(N, 0.3|units.MSun, min(100|units.MSun, cluster_particle.mass))

        r = cluster_particle.h_smooth  
        converter=nbody_system.nbody_to_si(masses.sum(),r)
        stars = new_plummer_model(N, convert_nbody=converter)
        stars.mass = masses
        stars.position += cluster_particle.position
        stars.velocity += cluster_particle.velocity
    return stars

def get_stars_from_molecular_clous(parts):
    cutoff_density = 10000 | units.amu/units.cm**3
    stars = Particles(0)
    for ip in parts:
        if ip.rho>cutoff_density:
            local_stars = make_stars(ip)
            stars.add_particles(local_stars)
    return stars

def run_dynamics(bodies, t_end, nsteps):

    Mtot_init = bodies.mass.sum()
    stellar = SeBa()
    stellar.parameters.metallicity = 0.02
    stellar.particles.add_particle(bodies)
    Rvir = 1|units.parsec
    converter=nbody_system.nbody_to_si(stars.mass.sum(),Rvir)

    gravity = ph4(converter)
#    gravity = bhtree(converter)
    gravity.parameters.timestep_parameter = 0.01

    gravity.particles.add_particles(bodies)

    channel_from_se_to_framework = stellar.particles.new_channel_to(bodies)
    channel_from_gd_to_framework = gravity.particles.new_channel_to(bodies)
    channel_from_framework_to_gd = bodies.new_channel_to(gravity.particles)
    channel_from_se_to_framework.copy_attributes(["mass","radius","luminosity"])

    bodies.scale_to_standard(convert_nbody=converter)
    
    filename = "GMC_stars.hdf5"
    write_set_to_file(bodies.savepoint(0|units.Myr), filename, 'hdf5')
    
    Etot_init = gravity.kinetic_energy + gravity.potential_energy
    Etot_prev = Etot_init

    time = 0.0 | t_end.unit
    dt = t_end/nsteps
    i = 0
    plot_stars(time, bodies, i)
    while time < t_end:
        time += dt
        gravity.evolve_model(time)
        Etot_prev_se = gravity.kinetic_energy + gravity.potential_energy

        stellar.evolve_model(time)

        channel_from_gd_to_framework.copy()
        channel_from_se_to_framework.copy_attributes(["mass","radius","luminosity", "temperature"])
        channel_from_framework_to_gd.copy_attributes(["mass"])

        write_set_to_file(bodies.savepoint(time), filename, 'hdf5')

        Ekin = gravity.kinetic_energy 
        Epot = gravity.potential_energy
        Etot = Ekin + Epot
        dE = Etot_prev-Etot
        dE_se = Etot_prev_se-Etot
        Mtot = bodies.mass.sum()
        print "T=", time, 
        print "M=", Mtot, "(dM[SE]=", Mtot/Mtot_init, ")",
        print "E= ", Etot, "Q= ", Ekin/Epot,
        print "dE=", (Etot_init-Etot)/Etot, "ddE=", (Etot_prev-Etot)/Etot, 
        print "(dE[SE]=", dE_se/Etot, ")"
        Etot_init -= dE
        Etot_prev = Etot

        plot_stars(time, bodies, i)
        i+=1

    gravity.stop()
    stellar.stop()

    plot_stars(time, bodies, i, L=20)
  
if __name__ in ("__main__","__plot__"):
#    parts = run_molecular_cloud(4000, Mcloud=400000. | units.MSun, Rcloud=10. | units.parsec)
    parts = run_molecular_cloud(1000, Mcloud=10000. | units.MSun, Rcloud=3. | units.parsec)
    stars = get_stars_from_molecular_clous(parts)
#    write_set_to_file(stars, "stars.hdf5", 'hdf5')
    run_dynamics(stars, 10.0|units.Myr, 10)

    
