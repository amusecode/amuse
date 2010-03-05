#1BPY
"""
Name:  'sun-earth'
Blender: 249
Group:'Add'
Tooltip: 'Amuse example'
"""
import Blender
import bpy
import pylab as pl
from Blender import Mesh

from amuse.legacy.hermite0.interface import HermiteInterface, Hermite
from amuse.support.data import core
from amuse.support.units import nbody_system
from amuse.support.units import units
from amuse.ext.blender import blender
import numpy as np

class SunEarth(object):
    
    def new_system_of_sun_and_earth(self):
        stars = core.Stars(2)
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(np.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(np.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = units.km(np.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(np.array((0.0,29800,0.0)))
        
        return stars
    
    def goaround(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        hermite = Hermite(convert_nbody)
        hermite.parameters.epsilon_squared = 0.0 | units.AU**2
        hermite.setup_module()
        hermite.dt_dia = 5000
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        sun = stars[0]
        Earth = blender.Primitives.sphere(10,10,0.1)
        Sun = blender.Primitives.sphere(32,32,1)

        hermite.setup_particles(stars)

        for i in range(1*365):
            hermite.evolve_model(i | units.day)
            hermite.update_particles(stars)
            Earth.loc = (1*earth.position.value_in(units.AU)[0],1*earth.position.value_in(units.AU)[1],earth.position.value_in(units.AU)[2])
            Sun.loc = (1*sun.position.value_in(units.AU)[0],1*sun.position.value_in(units.AU)[1],sun.position.value_in(units.AU)[2])
            #Earth.RotZ+=2*3.141592/24 	            
            Blender.Redraw() 

        del hermite
        
I = SunEarth()
I.goaround()

