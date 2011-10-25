from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.bonsai.interface import BonsaiInterface

import os
import sys 
import numpy

from amuse.community.octgrav.interface import OctgravInterface, Octgrav

from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
from amuse.rfi import channel
from amuse.ic.plummer import *

default_options = dict()
#default_options = dict(redirection="none")


class TestBonsaiInterface(TestWithMPI):
    
    def test0(self):
        print "Instantie aanmaken"
        instance = self.new_instance_of_an_optional_code(BonsaiInterface, **default_options)
        print "aangemaakt"
        result,error = instance.echo_int(12)
        print "call instance done"
        self.assertEquals(error, 0)
        self.assertEquals(result, 12)
        instance.stop()
    
    def test1(self):
        plummer_size = 500
#        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
        plummer =  new_plummer_sphere(plummer_size)
#        stars.radius = range(1, plummer_size+1)|units.km
        mass=plummer.mass.number
        radius=plummer.radius.number
        x=plummer.x.number
        y=plummer.y.number
        z=plummer.z.number
        vx=plummer.vx.number
        vy=plummer.vy.number
        vz=plummer.vz.number


#        instance = self.new_instance_of_an_optional_code(Octgrav, convert_nbody)
        instance = self.new_instance_of_an_optional_code(BonsaiInterface, **default_options)
        instance.initialize_code()
        #ids,err = instance.new_particles(mass,radius,x,y,z,vx,vy,vz)
        for i in range(0,plummer_size):
            ids,err = instance.new_particle(mass[i],radius[i],x[i],y[i],z[i],vx[i],vy[i],vz[i])
            #print ids,err
        instance.commit_particles()
        instance.evolve_model(1)

        #energy_total_init = instance.potential_energy + instance.kinetic_energy
        #instance.evolve_model(100 | units.day)
        #energy_total_final = instance.potential_energy + instance.kinetic_energy

        #self.assertAlmostRelativeEqual(energy_total_init, energy_total_final, 2)
        mass,radius,x,y,z,vx,vy,vz,err=instance.get_state(ids)

        instance.stop()

