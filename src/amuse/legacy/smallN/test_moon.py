import numpy
from amuse.legacy.smallN.muse_dynamics_mpi import SmallN
from amuse.support.units import nbody_system
from amuse.support.units.units import *
from amuse.support.data.core import Particle

if __name__=='__main__':
    myunits = nbody_system.nbody_to_si(1 | MSun, 1 | AU)

    # Orbital data
    pos_earth=numpy.array([ 8.418982185410142E-01,\
                            5.355823303978186E-01,\
                            2.327960005926782E-05]) | AU 
    vel_earth=numpy.array([-9.488931818313919E-03,\
                            1.447515189957170E-02,\
                            3.617712172296458E-07]) | AUd
    pos_moon=numpy.array([8.426656721530955E-01,\
                          5.331110650437484E-01,\
                         -6.837900390288286E-05]) | AU
    vel_moon=numpy.array([-8.943352544499154E-03,\
                           1.467795416516487E-02,\
                           4.840393580601162E-05]) | AUd

    earth = Particle()
    earth.mass = 5.9742e24 | kg
    earth.radius = 6371 | km 
    (earth.x, earth.y, earth.z) = pos_earth 
    (earth.vx, earth.vy, earth.vz) = vel_earth 

    moon = Particle()
    moon.mass = 7.3477e22 | kg 
    moon.radius = 1737.1 | km
    (moon.x, moon.y, moon.z) = pos_moon 
    (moon.vx, moon.vy, moon.vz) = vel_moon 


    cmpos=(pos_earth*earth.mass+pos_moon*moon.mass)/(earth.mass+moon.mass)
    cmvel=(vel_earth*earth.mass+vel_moon*moon.mass)/(earth.mass+moon.mass)
    
    pos_earth=pos_earth-cmpos
    vel_earth=vel_earth-cmvel
    pos_moon=pos_moon-cmpos
    vel_moon=vel_moon-cmvel

    nb=SmallN(myunits)
    
    earth = nb.particles.add_particle(earth)
    moon = nb.particles.add_particle(moon)

    nb.evolve()  
    
    print "moon:"
    print moon.mass
    print moon.position
    print moon.velocity
    
    print "earth:"
    print earth.mass
    print earth.position
    print earth.velocity
