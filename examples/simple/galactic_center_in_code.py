"""
Evolve a cluster orbiting a massive central particle.
"""

import numpy

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.datamodel import particles

from amuse.community.phiGRAPE.interface import PhiGRAPE

from matplotlib import pyplot

from amuse.ic.kingmodel import new_king_model
        
        
def new_galactic_center(mass):
    result = particles.Particle()
    result.mass = mass
    result.position = [0,0,0] | units.parsec
    result.velocity = [0,0,0] | units.kms
    result.radius = 0 | units.parsec
    return result
    
def circular_velocity_for_stable_orbit(radius, central_mass ):  
    return (constants.G*central_mass/radius)**0.5
    
def new_cluster(
        number_of_particles,
        W0,
        converter
    ):
            
    particles=new_king_model(N,W0,convert_nbody=converter)
    particles.radius=0.0| units.parsec
    return particles
    
def new_code(converter, timestep = 0.0025 | units.Myr):
    result=PhiGRAPE(converter)
    result.parameters.epsilon_squared = (0.01 | units.parsec)**2
    return result
    
def shift_particles(particles,dx,dy,dz,dvx,dvy,dvz):
    particles.x += dx
    particles.y += dy
    particles.z += dz
    particles.vx += dvx
    particles.vy += dvy
    particles.vz += dvz
    
if __name__ in ('__main__', '__plot__'):
    # parameter setup:
    N=256
    W0=3
    Rinit=50. | units.parsec
    timestep=0.01 | units.Myr
    Mcluster=4.e4 | units.MSun
    Rcluster=0.7 | units.parsec
    central_mass = 1.6e7 | units.MSun
    
    converter=nbody_system.nbody_to_si(Mcluster,Rcluster)
    cluster = new_cluster(number_of_particles = N, W0 = W0, converter = converter)
    code = new_code(converter, timestep = timestep)
    galactic_center = new_galactic_center(central_mass)
    
    # shift the cluster to an orbit around GC
    # (note the systems share the same coordinate frame, although units may differ)    
    vcirc=circular_velocity_for_stable_orbit(Rinit, central_mass)
    shift_particles(
        cluster,
        Rinit,0| units.parsec,0.|units.parsec,
        0| units.kms,0.8*vcirc,0| units.kms
    )

    code.particles.add_particles(cluster)
    
    # add galactic center to the code
    galactic_center = code.particles.add_particle(galactic_center)


    # evolve and make plots
    times=units.Myr([0.,0.2,0.4,0.6])
    f=pyplot.figure(figsize=(8,8))

    for i,t in enumerate(times):
        code.evolve_model(t,timestep=timestep)

        x=code.particles.x.value_in(units.parsec)
        y=code.particles.y.value_in(units.parsec)

        subplot=f.add_subplot(2,2,i+1)
        subplot.plot(x,y,'r .')
        subplot.plot([0.],[0.],'b +')
        subplot.set_xlim(-60,60)
        subplot.set_ylim(-60,60)
        subplot.set_title(t)
        if i==7:
            subplot.set_xlabel('parsec')
        
    code.stop()
    pyplot.show()
    
