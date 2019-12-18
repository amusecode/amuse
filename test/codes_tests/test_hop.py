import numpy
from numpy import random

from amuse.test.amusetest import get_path_to_results, TestCase
from amuse.units import units, nbody_system
from amuse.support.exceptions import AmuseException
from amuse.datamodel.particles import Particles
from amuse.ic.plummer import new_plummer_model
from amuse.community.hop.interface import HopInterface, Hop

class TestHopInterface(TestCase):
    def test1(self):
        print("First test: adding particles, setting and getting.")
        hop = HopInterface()
        hop.initialize_code()
        n, err = hop.get_number_of_particles()
        self.assertEqual(n, 0)
        self.assertEqual(err, 0)
        
        for i in range(6):
            id, err = hop.new_particle(0, i*i, 0, 0)
            n, err = hop.get_number_of_particles()
            self.assertEqual(n, i+1)
            self.assertEqual(err, 0)
        
        for i in range(6):
            x, y, z, err = hop.get_position(i)
            self.assertEqual(x, i*i)
            self.assertEqual(y, 0)
            self.assertEqual(z, 0)
            self.assertEqual(err, 0)
            
            hop.set_position(i, x, i*i, 0)
        
        for i in range(6):
            x, y, z, err  = hop.get_position(i)
            self.assertEqual(x, y)
            self.assertEqual(z, 0)
            self.assertEqual(err, 0)
        
        hop.stop()
    
    def test2(self):
        random.seed(1001)
        
        hop = HopInterface()
        hop.initialize_code()
        
        particles = new_plummer_model(1000)
        ids, errors = hop.new_particle(
            particles.mass.value_in(nbody_system.mass),
            particles.x.value_in(nbody_system.length),
            particles.y.value_in(nbody_system.length),
            particles.z.value_in(nbody_system.length)
        )
        
        n, err = hop.get_number_of_particles()
        self.assertEqual(n, 1000)
        self.assertEqual(err, 0)
        
        #distance_to_center = (particles.position - particles.center_of_mass()).lengths()
        
        #print distance_to_center
        ds = {0: 0.482308834791, 1:0.4885137677192688, 2:0.27442726492881775}
        for method in [0,1,2]:
            hop.set_nDens(7)
            hop.set_density_method(method)
            hop.calculate_densities()
            
            d, err = hop.get_density(0)
            self.assertAlmostRelativeEquals(d,ds[method], 5)
        hop.stop()
    
    def test3(self): 
    
        random.seed(1001)
        
        print("Third test: densest neighbors and groups.")
                
        hop = HopInterface()
        hop.initialize_code()
        
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
        self.assertEqual(n, 2)
        self.assertEqual(err, 0)
        n, err = hop.get_group_id(ids2)
        self.assertEqual(err, 0)
        n, err = hop.get_group_id(ids3)
        self.assertEqual(err, 0)
        
        n, err = hop.get_densest_particle_in_group(2)
        self.assertEqual(n, 7)
        for i in range(3):
            n, err = hop.get_number_of_particles_in_group(0)
            self.assertEqual(err, 0)
            self.assertEqual(n, 10)
            
        n, err = hop.get_number_of_groups()
        self.assertEqual(n, 3)
        
        n, err = hop.get_densest_neighbor(ids1)
        self.assertEqual(n, [7,7,12,0,7,7,7,7,12,7])
        hop.stop()
    
    def test4(self):
        hop = HopInterface()
        hop.initialize_code()
        value, error = hop.get_nDens()
        self.assertEqual(error,0)
        self.assertEqual(value,64)
        error = hop.set_nDens(7)
        self.assertEqual(error,0)
        value, error = hop.get_nDens()
        self.assertEqual(value,7)
        
        value, error = hop.get_nHop()
        self.assertEqual(error,0)
        self.assertEqual(value, -1)
        error = hop.set_nHop(7)
        self.assertEqual(error,0)
        value, error = hop.get_nHop()
        self.assertEqual(value,7)
        
        value, error = hop.get_nBucket()
        self.assertEqual(error,0)
        self.assertEqual(value, 16)
        error = hop.set_nHop(7)
        self.assertEqual(error,0)
        value, error = hop.get_nHop()
        self.assertEqual(value,7)
        hop.stop()
    
class TestHop(TestCase):
    def test1(self):
        print("First test: adding particles, setting and getting.")
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
            self.assertEqual(x, i*i | nbody_system.length)
            self.assertEqual(y, 0 | nbody_system.length)
            self.assertEqual(z, 0 | nbody_system.length)
            
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
        hop.stop()
    
    def test3(self): 
    
        random.seed(1001)
        
        print("Third test: densest neighbors and groups.")
                
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
        
        print(hop.particles.group_id)
        
        groups = list(hop.groups())
        
        self.assertEqual(len(groups), 3)
        
        self.assertEqual(hop.get_number_of_particles_outside_groups(), 0)
        
        #densities = (0,0,0) | nbody_system.density
        for index, group in enumerate(groups):
            self.assertEqual(len(group), 10)
            self.assertEqual(group.id_of_group(), index)
            #self.assertEquals(group.get_density_of_group(), densities[index])
        hop.stop()
    
    def test4(self):
        random.seed(1001)
        print("Test 4: complicated density field.")
        
        # A separate group below peak_density_threshold -> should be dropped
        particles0 = new_plummer_model(90, convert_nbody=nbody_system.nbody_to_si(0.9 | units.MSun, 1.0 | units.RSun))
        
        # A nearby group below peak_density_threshold -> should be attached to proper group
        particles1 = new_plummer_model(80, convert_nbody=nbody_system.nbody_to_si(0.8 | units.MSun, 1.0 | units.RSun))
        particles1.x += 10 | units.RSun
        
        # A proper group very nearby other proper group -> groups should merge
        particles2a = new_plummer_model(200, convert_nbody=nbody_system.nbody_to_si(2.0 | units.MSun, 1.0 | units.RSun))
        particles2b = new_plummer_model(300, convert_nbody=nbody_system.nbody_to_si(3.0 | units.MSun, 1.0 | units.RSun))
        particles2a.x += 11.0 | units.RSun
        particles2b.x += 11.2 | units.RSun
        
        # A separate proper group other proper group -> groups should be preserved
        particles3 = new_plummer_model(400, convert_nbody=nbody_system.nbody_to_si(4.0 | units.MSun, 1.0 | units.RSun))
        particles3.x += 20 | units.RSun
        
        hop = Hop(unit_converter=nbody_system.nbody_to_si(10.7 | units.MSun, 1.0 | units.RSun))
        hop.parameters.number_of_neighbors_for_local_density = 100
        hop.parameters.saddle_density_threshold_factor = 0.5
        hop.parameters.relative_saddle_density_threshold = True
        hop.commit_parameters()
        
        for set in [particles0, particles1, particles2a, particles2b, particles3]:
            hop.particles.add_particles(set)
        
        hop.calculate_densities()
        hop.parameters.outer_density_threshold = 0.1 * hop.particles.density.mean()
        hop.parameters.peak_density_threshold = hop.particles.density.amax() / 4.0
        hop.recommit_parameters()
        hop.do_hop()
        groups = list(hop.groups())
        
        self.assertEqual(len(hop.particles), 1070)
        self.assertEqual(len(groups), 2)
        self.assertEqual(hop.particles.select(lambda x: x < 5|units.RSun, "x").group_id, -1)
        self.assertEqual(hop.get_number_of_particles_outside_groups(), 299)
        self.assertEqual(1070 - len(groups[0]) - len(groups[1]), 299)
        
        expected_size = [477, 294] # Less than [580, 400], because particles below outer_density_threshold are excluded
        expected_average_x = [11.0, 20] | units.RSun
        for index, group in enumerate(groups):
            self.assertEqual(group.id_of_group(), index)
            self.assertAlmostEqual(group.center_of_mass()[0], expected_average_x[index], 1)
            self.assertEqual(len(group), expected_size[index])
        
        if False: # Make a plot
            original = hop.particles.copy()
            from amuse.plot import scatter, native_plot
            colors = ["r", "g", "b", "y", "k", "w"]*100
            for group, color in zip(hop.groups(), colors):
                scatter(group.x, group.y, c=color)
                original -= group
            scatter(original.x, original.y, c="m", marker="s")
            native_plot.show()
        
        hop.stop()
    
    def test5(self):
        print("Test error codes")
        unit_converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        particles = new_plummer_model(200, convert_nbody=unit_converter)
        hop = Hop(unit_converter=unit_converter)#, redirection="none")
        hop.parameters.number_of_neighbors_for_local_density = 100
        hop.particles.add_particles(particles[:99])
        self.assertRaises(AmuseException, hop.calculate_densities, expected_message=
            "Error when calling 'calculate_densities' of a 'Hop', errorcode is -5, error is 'Too few particles.'")
        hop.particles.add_particles(particles[99:101])
        hop.calculate_densities()
        hop.parameters.number_of_neighbors_for_hop = 200
        self.assertRaises(AmuseException, hop.calculate_densities, expected_message=
            "Error when calling 'calculate_densities' of a 'Hop', errorcode is -5, error is 'Too few particles.'")
        hop.particles.add_particles(particles[101:])
        hop.calculate_densities()
        
        self.assertRaises(AmuseException, hop.get_mass, 200, expected_message=
            "Error when calling 'get_mass' of a 'Hop', errorcode is -3, error is 'A particle with the given index was not found.'")
        hop.stop()
    
    def test6(self):
        print("Test with different masses")
        # Particles on a cubic grid with masses according to a gaussian density profile
        grid = numpy.mgrid[-1:1:21j, -1:1:21j, -1:1:21j] | units.m
        particles = Particles(9261, x=grid[0].flatten(), y=grid[1].flatten(), z=grid[2].flatten())
        peak_positions = [[0.2, -0.4, 0.3], [-0.6, 0.2, 0.7]] | units.m
        particles.mass = 2*numpy.exp(-(particles.position-peak_positions[0]).lengths_squared() / (0.1|units.m**2)) | units.kg
        particles.mass += numpy.exp(-(particles.position-peak_positions[1]).lengths_squared() / (0.1|units.m**2)) | units.kg
        self.assertAlmostEqual(particles.position[particles.mass.argmax()], peak_positions[0])
        self.assertAlmostEqual(particles[:4000].position[particles[:4000].mass.argmax()], peak_positions[1])
        
        hop = Hop(unit_converter=nbody_system.nbody_to_si(particles.mass.sum(), 1.0 | units.m))#, redirection="none")
        hop.parameters.density_method = 2
        hop.parameters.number_of_neighbors_for_local_density = 50
        hop.parameters.relative_saddle_density_threshold = True
        hop.commit_parameters()
        hop.particles.add_particles(particles)
        hop.calculate_densities()
        self.assertAlmostEqual(hop.particles.position[hop.particles.density.argmax()], peak_positions[0])
        self.assertAlmostEqual(hop.particles[:4000].position[hop.particles[:4000].density.argmax()], peak_positions[1])
        hop.do_hop()
        groups = list(hop.groups())
        self.assertEqual(len(groups), 2)
        for group, peak_position in zip(groups, peak_positions):
            self.assertAlmostEqual(group.center_of_mass(), peak_position, 1)
        hop.stop()
    
    def test7(self):
        print("Testing Hop states")
        unit_converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        particles = new_plummer_model(200, convert_nbody=unit_converter)
        
        print("First do everything manually:", end=' ')
        instance = Hop(unit_converter=unit_converter)
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.particles.add_particles(particles)
        instance.commit_particles()
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.cleanup_code()
        self.assertEqual(instance.get_name_of_current_state(), 'END')
        instance.stop()
        print("ok")

        print("initialize_code(), commit_parameters(), (re)commit_particles(), " \
            "and cleanup_code() should be called automatically:", end=' ')
        instance = Hop(unit_converter=unit_converter)
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.number_of_neighbors_for_local_density = 50
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particles(particles[:100])
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        mass = instance.particles[0].mass
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.particles.add_particles(particles[100:])
        self.assertEqual(instance.get_name_of_current_state(), 'UPDATE')
        mass = instance.particles[100].mass
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        instance.stop()
        self.assertEqual(instance.get_name_of_current_state(), 'STOPPED')
        print("ok")
    
    def test8(self):
        random.seed(1001)
        print("Test 8: SI vs nbody units.")
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.RSun)
        numpy.random.seed(1234)
        particles = new_plummer_model(1000, convert_nbody=converter)
        hop = Hop(unit_converter=converter)
        hop.particles.add_particles(particles)
        hop.calculate_densities()
        hop.parameters.outer_density_threshold = 0.1 | nbody_system.mass / nbody_system.length**3
        hop.do_hop()
        groups = list(hop.groups())
        
        self.assertEqual(len(hop.particles), 1000)
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(groups[0]), 511)
        self.assertEqual(len(hop.no_group()), 489)
        hop.stop()
        
        numpy.random.seed(1234)
        particles = new_plummer_model(1000)
        hop = Hop()
        hop.particles.add_particles(particles)
        hop.calculate_densities()
        hop.parameters.outer_density_threshold = 0.1 | nbody_system.mass / nbody_system.length**3
        hop.do_hop()
        groups = list(hop.groups())
        
        self.assertEqual(len(hop.particles), 1000)
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(groups[0]), 511)
        self.assertEqual(len(hop.no_group()), 489)
        hop.stop()
    


