import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.community.simplex.interface import SimpleXInterface, SimpleX
from amuse.units import units
from amuse.support.data import Particles
default_options = dict(number_of_workers=2)
default_options = dict(number_of_workers=2, redirection="none")

class TestSimpleXInterface(TestWithMPI):

    def test1(self):
        print "Test 1: initialization"
        instance = SimpleXInterface(**default_options)
        self.assertEqual(0, instance.set_output_directory(instance.output_directory))
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()
    
    def test2(self):
        print "Test 2: commit_particles, getters and setters"
        instance = SimpleXInterface(**default_options)
        self.assertEqual(0, instance.set_output_directory(instance.output_directory))
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        
        input_file = os.path.join(instance.data_directory, 'vertices_test2.txt')
        x, y, z, n_H, flux, X_ion = read_input_file(input_file)
        x=numpy.array(x)
        y=numpy.array(y)
        z=numpy.array(z)
        number_of_particles = len(x)
        indices, errors = instance.new_particle(x, y, z, n_H, flux, X_ion)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertEqual(indices, range(number_of_particles))
        self.assertEqual(0, instance.commit_particles())
        x_out, y_out, z_out, n_H_out, flux_out, X_ion_out, error = instance.get_state(indices)
        self.assertAlmostEqual((x-x_out)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((y-y_out)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((z-z_out)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual(flux,flux_out, 7)
        self.assertAlmostEqual(n_H,n_H_out, 7)
        self.assertAlmostEqual(X_ion,X_ion_out, 7)
        
        
#        for expected, received in zip([13200.*x, 13200.*y, 13200.*z, n_H, flux, X_ion, [0]*number_of_particles], 
#                [13200.*x_out, 13200.*y_out, 13200.*z_out, n_H_out, flux_out, X_ion_out, error]):
#            self.assertAlmostEqual(expected, received, 5)
        
        x, y, z, n_H, flux, X_ion, error = instance.get_state(0)
        for expected, received in zip([0, 0, 0, 0.001, 5.0, 0.0, 0], 
                [x, y, z, n_H, flux, X_ion, error]):
            self.assertAlmostRelativeEqual(expected, received,6)
        x, y, z, error1 = instance.get_position(0)
        n_H, error2     = instance.get_density(0)
        flux, error3    = instance.get_flux(0)
        X_ion, error4   = instance.get_ionisation(0)
        for expected, received in zip([0.,0.,0., 0.001, 5.0, 0.0, 0, 0, 0, 0], 
                [x, y, z, n_H, flux, X_ion, error1, error2, error3, error4]):
            self.assertAlmostRelativeEqual(expected, received, 5)
        
        self.assertEqual(0, instance.set_state(3, 1.0, 2.0, 3.0, 4.0, 5.0, 0.6))
        x, y, z, n_H, flux, X_ion, error = instance.get_state(3)
        for expected, received in zip([1.0, 2.0, 3.0, 4.0, 5.0, 0.6, 0], 
                [x, y, z, n_H, flux, X_ion, error]):
            self.assertAlmostRelativeEqual(expected, received, 5)
        self.assertEqual(0, instance.set_position(4, 3.0, 2.0, 1.0))
        self.assertEqual(0, instance.set_density(4, 0.6))
        self.assertEqual(0, instance.set_flux(4, 0.5))
        self.assertEqual(0, instance.set_ionisation(4, 0.4))
        x, y, z, n_H, flux, X_ion, error = instance.get_state(4)
        for expected, received in zip([3.0, 2.0, 1.0, 0.6, 0.5, 0.4, 0], 
                [x, y, z, n_H, flux, X_ion, error]):
            self.assertAlmostRelativeEqual(expected, received, 5)
        
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()
    
    def test3(self):
        print "Test 3: evolve"
        instance = SimpleXInterface(**default_options)
        self.assertEqual(0, instance.set_output_directory(instance.output_directory))
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        
        input_file = os.path.join(instance.data_directory, 'vertices_test2.txt')
        x, y, z, n_H, flux, X_ion = read_input_file(input_file)
        number_of_particles = len(x)
        indices, errors = instance.new_particle(x, y, z, n_H, flux, X_ion)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertEqual(indices, range(number_of_particles))
        self.assertEqual(0, instance.commit_particles())
        X_ion, errors = instance.get_ionisation(indices)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertAlmostEqual(X_ion.sum()/number_of_particles, 0.0)
        self.assertEqual(0, instance.evolve_model(0.5))
        X_ion, errors = instance.get_ionisation(indices)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertAlmostEqual(X_ion.sum()/number_of_particles, 0.000933205)
        
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test4(self):
        print "Test 4: set boxsize, hilbert_order, timestep"
        instance = SimpleXInterface(**default_options)
        self.assertEqual(0, instance.set_output_directory(instance.output_directory))
        self.assertEqual(0, instance.initialize_code())
        
        instance.set_box_size_parameter(16384.)
        instance.set_hilbert_order_parameter(1)
        instance.set_timestep_parameter(0.5)
                
        self.assertEqual(0, instance.commit_parameters())
        
        self.assertEqual(16384.,instance.get_box_size_parameter()['box_size'])
        self.assertEqual(1,instance.get_hilbert_order_parameter()['hilbert_order'])
        self.assertEqual(0.5,instance.get_timestep_parameter()['timestep'])

        input_file = os.path.join(instance.data_directory, 'vertices_test2.txt')
        x, y, z, n_H, flux, X_ion = read_input_file(input_file)
        number_of_particles = len(x)
        indices, errors = instance.new_particle(x, y, z, n_H, flux, X_ion)
        instance.commit_particles()

        self.assertEqual(16384.,instance.get_box_size_parameter()['box_size'])
        self.assertEqual(1,instance.get_hilbert_order_parameter()['hilbert_order'])
#        self.assertEqual(0.5,instance.get_timestep_parameter()['timestep'])

        self.assertEqual(0, instance.cleanup_code())
        instance.stop()
    

class TestSimpleX(TestWithMPI):

    def test1(self):
        print "Test 1: initialization"
        instance = SimpleX(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Test 2: commit_particles, getters and setters"
        instance = SimpleX(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        
        input_file = os.path.join(instance.data_directory, 'vertices_test2.txt')
        particles = particles_from_input_file(input_file)
        instance.particles.add_particles(particles)
        instance.commit_particles()
#        for attribute in ['position', 'rho', 'flux', 'xion']:
#            self.assertAlmostEqual(13200.*getattr(particles, attribute),
#                                   13200.*getattr(instance.particles, attribute), 5)
#            setattr(instance.particles, attribute, getattr(particles, attribute)/2.0)
#            self.assertAlmostEqual(13200.*getattr(particles, attribute)/2.0,
#                                   13200.*getattr(instance.particles, attribute), 5)
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):
        print "Test 3: evolve"
        instance = SimpleX(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        
        input_file = os.path.join(instance.data_directory, 'vertices_test2.txt')
        particles = particles_from_input_file(input_file)
        instance.particles.add_particles(particles)
        instance.commit_particles()
        self.assertAlmostEqual(instance.particles.xion.mean(), 0.0 | units.none)
        instance.evolve_model(0.5 | units.Myr)
        self.assertAlmostEqual(instance.particles.xion.mean(), 0.000933205 | units.none)
        instance.cleanup_code()
        instance.stop()
    
def read_input_file(input_file):
    file = open(input_file, 'r')
    lines = file.readlines()
    lines.pop(0)
    x, y, z, nh, flux, xion = [], [], [], [], [], []
    for line in lines:
        l = line.strip().split()
        if len(l) >= 7:
            x.append(float(l[1]))
            y.append(float(l[2]))
            z.append(float(l[3]))
            nh.append(float(l[4]))
            flux.append(float(l[5]))
            xion.append(float(l[6]))
    return x, y, z, nh, flux, xion

def particles_from_input_file(input_file):
    x, y, z, n_H, flux, X_ion = read_input_file(input_file)
    particles = Particles(len(x))
    particles.x = x | units.parsec
    particles.y = y | units.parsec
    particles.z = z | units.parsec
    particles.rho = n_H | units.amu / units.cm**3
    particles.flux = flux | 1.0e48 / units.s
    particles.xion = X_ion | units.none
    return particles

