import os
import sys
import numpy

from numpy import random
from amuse.test import amusetest
from amuse.test.amusetest import get_path_to_results



from amuse.community.fi import interface as interface
from amuse.community.hop.interface import HopInterface, Hop
from amuse.units import units
from amuse.units import nbody_system
from amuse.datamodel.particles import Particles
from amuse.rfi.core import is_mpd_running
from amuse.ic.plummer import new_plummer_model

class TestHopInterface(amusetest.TestCase):
    def test1(self):
        print "First test: adding particles, setting and getting."
        hop = HopInterface()
        n, err = hop.get_number_of_particles()
        self.assertEquals(n, 0)
        self.assertEquals(err, 0)
        
        for i in range(6):
            id, err = hop.new_particle(0, i*i, 0, 0)
            n, err = hop.get_number_of_particles()
            self.assertEquals(n, i+1)
            self.assertEquals(err, 0)
        
        for i in range(6):
            x, y, z, err = hop.get_position(i)
            self.assertEquals(x, i*i)
            self.assertEquals(y, 0)
            self.assertEquals(z, 0)
            self.assertEquals(err, 0)
            
            hop.set_position(i, x, i*i, 0)
        
        for i in range(6):
            x, y, z, err  = hop.get_position(i)
            self.assertEquals(x, y)
            self.assertEquals(z, 0)
            self.assertEquals(err, 0)
        
        hop.stop()
            
    def test2(self):
        random.seed(1001)
        
        hop = HopInterface()
        
        particles = new_plummer_model(1000)
        ids, errors = hop.new_particle(
            particles.mass.value_in(nbody_system.mass),
            particles.x.value_in(nbody_system.length),
            particles.y.value_in(nbody_system.length),
            particles.z.value_in(nbody_system.length)
        )
        
        n, err = hop.get_number_of_particles()
        self.assertEquals(n, 1000)
        self.assertEquals(err, 0)
        
        #distance_to_center = (particles.position - particles.center_of_mass()).lengths()
        
        #print distance_to_center
        ds = {0: 0.482308834791, 1:0.4885137677192688, 2:0.27442726492881775}
        for method in [0,1,2]:
            hop.set_nDens(7)
            hop.set_density_method(method)
            hop.calculate_densities()
            
            d, err = hop.get_density(0)
            self.assertAlmostRelativeEquals(d,ds[method], 5)
    
    def test3(self): 
    
        random.seed(1001)
        
        print "Third test: densest neighbors and groups."
                
        hop = HopInterface()
        
        particles1 = new_plummer_model(10)
        particles2 = new_plummer_model(10)
        particles3 = new_plummer_model(10)
        
        particles2.position += (10,0,0) | nbody_system.length
        
        particles3.position += (0,20,0) | nbody_system.length
        
        ids1, errors = hop.new_particle(
            particles1.mass.value_in(nbody_system.mass),
            particles1.x.value_in(nbody_system.length),
            particles1.y.value_in(nbody_system.length),
            particles1.z.value_in(nbody_system.length)
        )
        
        ids2, errors = hop.new_particle(
            particles2.mass.value_in(nbody_system.mass),
            particles2.x.value_in(nbody_system.length),
            particles2.y.value_in(nbody_system.length),
            particles2.z.value_in(nbody_system.length)
        )
        
        ids3, errors = hop.new_particle(
            particles3.mass.value_in(nbody_system.mass),
            particles3.x.value_in(nbody_system.length),
            particles3.y.value_in(nbody_system.length),
            particles3.z.value_in(nbody_system.length)
        )
        
        
        hop.set_nDens(5)
        hop.calculate_densities()
        hop.do_hop()
        
        n, err = hop.get_group_id(ids1)
        self.assertEquals(n, 2)
        self.assertEquals(err, 0)
        n, err = hop.get_group_id(ids2)
        self.assertEquals(err, 0)
        n, err = hop.get_group_id(ids3)
        self.assertEquals(err, 0)
        
        n, err = hop.get_densest_particle_in_group(2)
        self.assertEquals(n, 7)
        for i in range(3):
            n, err = hop.get_number_of_particles_in_group(0)
            self.assertEquals(err, 0)
            self.assertEquals(n, 10)
            
        n, err = hop.get_number_of_groups()
        self.assertEquals(n, 3)
        
        n, err = hop.get_densest_neighbor(ids1)
        self.assertEquals(n, [7,7,12,0,7,7,7,7,12,7])
        

    def test4(self):
        hop = HopInterface()
        hop.initialize_code()
        value, error = hop.get_nDens()
        self.assertEquals(error,0)
        self.assertEquals(value,64)
        error = hop.set_nDens(7)
        self.assertEquals(error,0)
        value, error = hop.get_nDens()
        self.assertEquals(value,7)
        
        value, error = hop.get_nHop()
        self.assertEquals(error,0)
        self.assertEquals(value, -1)
        error = hop.set_nHop(7)
        self.assertEquals(error,0)
        value, error = hop.get_nHop()
        self.assertEquals(value,7)
        
        value, error = hop.get_nBucket()
        self.assertEquals(error,0)
        self.assertEquals(value, 16)
        error = hop.set_nHop(7)
        self.assertEquals(error,0)
        value, error = hop.get_nHop()
        self.assertEquals(value,7)
        
class TestHop(amusetest.TestCase):
    def test1(self):
        print "First test: adding particles, setting and getting."
        hop = Hop()
        particles = Particles(6)
        particles.mass = 1.0 | nbody_system.mass
        particles.x = [i*i for i in range(6)] | nbody_system.length
        particles.y = 0.0 | nbody_system.length
        particles.z = 0.0 | nbody_system.length
        
        hop.particles.add_particles(particles)
        
        positions = hop.particles.position
        for i in range(6):
            x, y, z = positions[i]
            self.assertEquals(x, i*i | nbody_system.length)
            self.assertEquals(y, 0 | nbody_system.length)
            self.assertEquals(z, 0 | nbody_system.length)
            
        hop.stop()
            
    def test2(self):
        random.seed(1001)
        
        hop = Hop()
        hop.initialize_code()
        hop.parameters.number_of_neighbors_for_local_density = 7
        hop.commit_parameters()
        
        particles = new_plummer_model(1000)
        hop.particles.add_particles(particles)
        
        
        #distance_to_center = (particles.position - particles.center_of_mass()).lengths()
        
        #print distance_to_center
        ds = {0: 0.482308834791, 1:0.4885137677192688, 2:0.27442726492881775}
        for method in [0,1,2]:
            hop.set_density_method(method)
            hop.calculate_densities()
            
            d = hop.particles[0].density
            
            self.assertAlmostRelativeEquals(d, ds[method] | nbody_system.density, 5)
    
    def test3(self): 
    
        random.seed(1001)
        
        print "Third test: densest neighbors and groups."
                
        hop = Hop()
        hop.parameters.number_of_neighbors_for_local_density = 5
        
        particles1 = new_plummer_model(10)
        particles2 = new_plummer_model(10)
        particles3 = new_plummer_model(10)
        
        particles2.position += (10,0,0) | nbody_system.length
        
        particles3.position += (0,20,0) | nbody_system.length
        
        hop.particles.add_particles(particles1)
        hop.particles.add_particles(particles2)
        hop.particles.add_particles(particles3)        
        
        hop.calculate_densities()
        hop.do_hop()
        
        print hop.particles.group_id
        
        groups = list(hop.groups())
        
        self.assertEquals(len(groups), 3)
        
        self.assertEquals(hop.get_number_of_particles_outside_groups(), 0)
        
        #densities = (0,0,0) | nbody_system.density
        for index, group in enumerate(groups):
            self.assertEquals(len(group), 10)
            self.assertEquals(group.id_of_group(), index)
            #self.assertEquals(group.get_density_of_group(), densities[index])

