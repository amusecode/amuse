import os.path
import numpy
from amuse.test.amusetest import TestWithMPI

from amuse.community.simplex.interface import SimpleXInterface, SimpleX,SimpleXSplitSet
from amuse.units import units
from amuse.datamodel import Particles
#default_options = dict(number_of_workers=2, redirection="none")
default_options = dict(number_of_workers=1)#,debugger='gdb')

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
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        x, y, z, n_H, flux, X_ion,u = read_input_file(input_file)
        x=numpy.array(x)
        y=numpy.array(y)
        z=numpy.array(z)
        number_of_particles = len(x)
        indices, errors = instance.new_particle(x, y, z, n_H, flux, X_ion,u)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertEqual(indices, range(number_of_particles))
        self.assertEqual(0, instance.commit_particles())
        x_out, y_out, z_out, n_H_out, flux_out, X_ion_out,u_out, metallicity_out, error = instance.get_state(indices)

        self.assertAlmostEqual((x-x_out)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((y-y_out)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual((z-z_out)/13200., numpy.zeros_like(x), 7)
        self.assertAlmostEqual(flux,flux_out, 7)
        self.assertAlmostEqual(n_H,n_H_out, 7)
        self.assertAlmostEqual(X_ion,X_ion_out, 7)
        self.assertAlmostRelativeEqual(u,u_out, 7)
        #self.assertAlmostRelativeEqual(metallicity,metallicity_out, 7)
        
        
        
        x, y, z, n_H, flux, X_ion,u,metallicity, error = instance.get_state(0)
        for expected, received in zip([0, 0, 0, 0.001, 5.0, 0.0, 831447247704,0], 
                [x, y, z, n_H, flux, X_ion, u,error]):
            self.assertAlmostRelativeEqual(expected, received,6)
        x, y, z, error1 = instance.get_position(0)
        n_H, error2     = instance.get_density(0)
        flux, error3    = instance.get_flux(0)
        X_ion, error4   = instance.get_ionisation(0)
        for expected, received in zip([0.,0.,0., 0.001, 5.0, 0.0, 0, 0, 0, 0], 
                [x, y, z, n_H, flux, X_ion, error1, error2, error3, error4]):
            self.assertAlmostRelativeEqual(expected, received, 5)
        
        self.assertEqual(0, instance.set_state(3, 1.0, 2.0, 3.0, 4.0, 5.0, 0.6,77.0))
        x, y, z, n_H, flux, X_ion, u,metallicity,error = instance.get_state(3)
        for expected, received in zip([1.0, 2.0, 3.0, 4.0, 5.0, 0.6, 77,0], 
                [x, y, z, n_H, flux, X_ion,u, error]):
            self.assertAlmostRelativeEqual(expected, received, 5)
        self.assertEqual(0, instance.set_position(4, 3.0, 2.0, 1.0))
        self.assertEqual(0, instance.set_density(4, 0.6))
        self.assertEqual(0, instance.set_flux(4, 0.5))
        self.assertEqual(0, instance.set_ionisation(4, 0.4))
        self.assertEqual(0, instance.set_internal_energy(4, 1234.))
        x, y, z, n_H, flux, X_ion,u,metallicity, error = instance.get_state(4)
        for expected, received in zip([3.0, 2.0, 1.0, 0.6, 0.5, 0.4,1234., 0], 
                [x, y, z, n_H, flux, X_ion,u, error]):
            self.assertAlmostRelativeEqual(expected, received, 5)


        self.assertEqual(0, instance.set_dinternal_energy_dt(4, 12345.))
        du_dt,err=instance.get_dinternal_energy_dt(4)
        self.assertEqual(12345, du_dt)
        
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()
    
    def test3(self):
        print "Test 3: evolve"
        instance = SimpleXInterface(**default_options)
        self.assertEqual(0, instance.set_output_directory(instance.output_directory))
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        x, y, z, n_H, flux, X_ion,u = read_input_file(input_file)

        number_of_particles = len(x)
        indices, errors = instance.new_particle(x, y, z, n_H, flux, X_ion,u)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertEqual(indices, range(number_of_particles))
        
        self.assertEqual(0, instance.commit_particles())
        
        X_ion, errors = instance.get_ionisation(indices)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertAlmostEqual(X_ion.sum()/number_of_particles, 0.0)
        
        density, errors = instance.get_density(indices)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertAlmostEqual(density.sum(), 1.0,6)
        
        self.assertEqual(0, instance.evolve_model(0.5))

        density, errors = instance.get_density(indices)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertAlmostEqual(density.sum(), 1.0,6)
        
        flux, errors = instance.get_flux(indices)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertEqual(flux.sum(), 5.0)
        
        X_ion, errors = instance.get_ionisation(indices)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertAlmostEqual(X_ion.sum()/number_of_particles, 0.000845247683257)
        
        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test4(self):
        print "Test 4: set boxsize, hilbert_order, timestep"
        instance = SimpleXInterface(**default_options)
        self.assertEqual(0, instance.set_output_directory(instance.output_directory))
        self.assertEqual(0, instance.initialize_code())
        
        instance.set_box_size(16384.)
        instance.set_hilbert_order(1)
        instance.set_timestep(0.5)
                
        self.assertEqual(0, instance.commit_parameters())
        
        self.assertEqual(16384.,instance.get_box_size()['box_size'])
        self.assertEqual(1,instance.get_hilbert_order()['hilbert_order'])
        self.assertEqual(0.5,instance.get_timestep()['timestep'])

        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        x, y, z, n_H, flux, X_ion,u = read_input_file(input_file)
        number_of_particles = len(x)
        indices, errors = instance.new_particle(x, y, z, n_H, flux, X_ion,u)
        instance.commit_particles()

        self.assertEqual(16384.,instance.get_box_size()['box_size'])
        self.assertEqual(1,instance.get_hilbert_order()['hilbert_order'])
        self.assertEqual(0.5,instance.get_timestep()['timestep'])

        self.assertEqual(0, instance.cleanup_code())
        instance.stop()

    def test5(self):
        print "Test 2: delete particles"
        instance = SimpleXInterface(**default_options)
        self.assertEqual(0, instance.set_output_directory(instance.output_directory))
        self.assertEqual(0, instance.initialize_code())
        self.assertEqual(0, instance.commit_parameters())
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        x, y, z, n_H, flux, X_ion,u = read_input_file(input_file)
        x=numpy.array(x)
        y=numpy.array(y)
        z=numpy.array(z)
        number_of_particles = len(x)
        indices, errors = instance.new_particle(x, y, z, n_H, flux, X_ion,u)
        self.assertEqual(errors, [0]*number_of_particles)
        self.assertEqual(indices, range(number_of_particles))
        error=instance.delete_particle(indices[0])
        self.assertEqual(error, -1)
        instance.commit_particles()
        error=instance.delete_particle(indices[0])
# this one I don't understand:
        self.assertEqual(error, -1)
        self.assertEqual(0, instance.evolve_model(0.125))
        error=instance.delete_particle(indices[0])
        self.assertEqual(error, 0)


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
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
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
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        particles = particles_from_input_file(input_file)
        particles.du_dt = particles.u/(10|units.Myr)
        instance.particles.add_particles(particles)
#        instance.particles.du_dt=particles.du_dt
#        instance.commit_particles()
        instance.particles.du_dt=particles.du_dt
        self.assertAlmostEqual(instance.particles.xion.mean(), 0.0)
        self.assertAlmostEqual(instance.particles.du_dt.mean().in_(units.cm**2/units.s**3),particles.du_dt.mean().in_(units.cm**2/units.s**3))
        instance.evolve_model(0.5 | units.Myr)
        self.assertAlmostEqual(instance.particles.du_dt.mean().in_(units.cm**2/units.s**3),particles.du_dt.mean().in_(units.cm**2/units.s**3))
        self.assertAlmostEqual(instance.particles.xion.mean(), 0.000845247683257)
        instance.cleanup_code()
        instance.stop()

    def test4(self):
        print "Test 4: default parameters"
        instance = SimpleX(**default_options)

        default=dict( timestep= 0.05| units.Myr, 
                  source_effective_T=  1.e5 | units.K,
                  hilbert_order= 1,
                  number_of_freq_bins= 1,
                  thermal_evolution_flag = 0,
                  blackbody_spectrum_flag = 0,
                  box_size=13200 | units.parsec,
                  metal_cooling_flag=0,
                  collisional_ionization_flag=0)
        for x in default:
            self.assertEqual(getattr(instance.parameters,x), default[x])
        instance.commit_parameters()
        for x in default:
            self.assertEqual(getattr(instance.parameters,x), default[x])

        tnow=instance.model_time
        self.assertEqual(tnow, 0. | units.Myr)    
        instance.model_time=321. | units.Myr
        tnow=instance.model_time
        self.assertEqual(tnow, 321. | units.Myr)    

    def test5(self):
        print "Test 4: default parameters"
        instance = SimpleX(**default_options)
        param=dict( timestep= 0.1| units.Myr, 
                  source_effective_T=  2.e5 | units.K,
                  hilbert_order= 3,
                  number_of_freq_bins= 4,
                  thermal_evolution_flag = 1,
                  blackbody_spectrum_flag = 1,
                  box_size=32100 | units.parsec,
                  metal_cooling_flag=1,
                  collisional_ionization_flag=1,
                  simplex_data_directory='.')
        for x in param:
            setattr(instance.parameters,x, param[x])
        for x in param:
            self.assertEqual(getattr(instance.parameters,x), param[x])

    def test6(self):
        print "Test 2: print parameters,data directory"
        instance = SimpleX(**default_options)
        print instance.parameters
        
        instance.parameters.simplex_data_directory="some/dir"
        self.assertEqual(instance.parameters.simplex_data_directory, "some/dir")

    def test7(self):
        print "Test 7: two step evolve"
        instance = SimpleX(**default_options)
        instance.parameters.recombination_radiation_flag=1
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        particles = particles_from_input_file(input_file)
        instance.particles.add_particles(particles)
        self.assertAlmostEqual(instance.particles.xion.mean(), 0.0)
        instance.evolve_model(0.25 | units.Myr)
        instance.evolve_model(0.5 | units.Myr)
        self.assertAlmostRelativeEqual(instance.particles.xion.mean(),0.00084660917243,3)
        self.assertEqual( instance.particles.flux.max().value_in(1.e48* units.s**-1), 5)
        instance.cleanup_code()
        instance.stop()

    def test8(self):
        print "Test 8: two step evolve"
        instance = SimpleX(**default_options)
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        particles = particles_from_input_file(input_file)
        instance.particles.add_particles(particles)
        self.assertAlmostEqual(instance.particles.xion.mean(), 0.0)
        instance.evolve_model(0.25 | units.Myr)
        instance.evolve_model(0.5 | units.Myr)
        self.assertEqual( instance.particles.flux.max().value_in(1.e48* units.s**-1), 5)
        self.assertAlmostRelativeEqual(instance.particles.xion.mean(),0.00084660917243,3)
        instance.cleanup_code()
        instance.stop()

    def test9(self):
        print "Test 9: add test"
        instance = SimpleX(number_of_workers=1)
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        particles = particles_from_input_file(input_file)
        N=len(particles)
        toadd=particles[0:10].copy()
        particles=particles[10:].copy()
        instance.particles.add_particles(particles)
        instance.evolve_model(0.25 | units.Myr)
#        instance.commit_particles()
        instance.particles.add_particles(toadd)
        self.assertEqual( len(instance.particles), N)

    def test10(self):
        print "Test 10: add test"
        instance = SimpleX(number_of_workers=2)
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        particles = particles_from_input_file(input_file)
        N=len(particles)
        toadd=particles[0:10].copy()
        particles=particles[10:].copy()
        instance.particles.add_particles(particles)
        instance.evolve_model(0.25 | units.Myr)
#        instance.commit_particles()
        instance.particles.add_particles(toadd)
        self.assertEqual( len(instance.particles), N)

    def test11(self):
        print "Test 11: add test"
        instance = SimpleX(number_of_workers=1)
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        particles = particles_from_input_file(input_file)
        N=len(particles)
        toadd=particles[0:10].copy()
        particles=particles[10:].copy()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        instance.particles.add_particles(toadd)
        self.assertEqual( len(instance.particles), N)


class TestSimpleXSplitSet(TestWithMPI):

    def test1(self):
        print "Test 1: initialization"
        instance = SimpleXSplitSet(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Test 2: commit_particles, getters and setters"
        instance = SimpleXSplitSet(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        particles,src_particles = splitset_from_input_file(input_file)
        instance.gas_particles.add_particles(particles)
        instance.src_particles.add_particles(src_particles)
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
        instance = SimpleXSplitSet(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        
        input_file = os.path.join(instance.data_directory, 'vertices_test3.txt')
        particles,src_particles = splitset_from_input_file(input_file)
        instance.src_particles.add_particles(src_particles)
        particles.du_dt = particles.u/(10|units.Myr)
        instance.gas_particles.add_particles(particles)

        self.assertAlmostEqual(instance.gas_particles.xion.mean(), 0.0)
        self.assertAlmostEqual(instance.gas_particles.du_dt.mean().in_(units.cm**2/units.s**3),particles.du_dt.mean().in_(units.cm**2/units.s**3))
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.evolve_model(0.5 | units.Myr)
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        self.assertAlmostEqual(instance.gas_particles.du_dt.mean().in_(units.cm**2/units.s**3),particles.du_dt.mean().in_(units.cm**2/units.s**3))
        self.assertAlmostEqual(instance.gas_particles.xion.mean(), 0.000845247683257)
        instance.gas_particles.remove_particles(particles[0:4])
# this is what we would like....
#        self.assertEquals(instance.get_name_of_current_state(), 'UPDATE')
#        instance.recommit_particles()
        instance.evolve_model(0.75 | units.Myr)
        self.assertEqual(len(instance.particles), len(particles)-4)
        instance.cleanup_code()
        instance.stop()

    def test4(self):
        print "Test 4: default parameters"
        instance = SimpleXSplitSet(**default_options)

        default=dict( timestep= 0.05| units.Myr, 
                  source_effective_T=  1.e5 | units.K,
                  hilbert_order= 1,
                  number_of_freq_bins= 1,
                  thermal_evolution_flag = 0,
                  blackbody_spectrum_flag = 0,
                  box_size=13200 | units.parsec,
                  metal_cooling_flag=0,
                  collisional_ionization_flag=0)
        for x in default:
            self.assertEqual(getattr(instance.parameters,x), default[x])
        instance.commit_parameters()
        for x in default:
            self.assertEqual(getattr(instance.parameters,x), default[x])

        tnow=instance.model_time
        self.assertEqual(tnow, 0. | units.Myr)    
        instance.model_time=321. | units.Myr
        tnow=instance.model_time
        self.assertEqual(tnow, 321. | units.Myr)    

    def test5(self):
        print "Test 4: default parameters"
        instance = SimpleXSplitSet(**default_options)
        param=dict( timestep= 0.1| units.Myr, 
                  source_effective_T=  2.e5 | units.K,
                  hilbert_order= 3,
                  number_of_freq_bins= 4,
                  thermal_evolution_flag = 1,
                  blackbody_spectrum_flag = 1,
                  box_size=32100 | units.parsec,
                  metal_cooling_flag=1,
                  collisional_ionization_flag=1)
        for x in param:
            setattr(instance.parameters,x, param[x])
        for x in param:
            self.assertEqual(getattr(instance.parameters,x), param[x])


    
def read_input_file(input_file):
    file = open(input_file, 'r')
    lines = file.readlines()
    lines.pop(0)
    x, y, z, nh, flux, xion,u = [], [], [], [], [], [],[]
    for line in lines:
        l = line.strip().split()
        if len(l) >= 7:
            x.append(float(l[1]))
            y.append(float(l[2]))
            z.append(float(l[3]))
            nh.append(float(l[4]))
            flux.append(float(l[5]))
            xion.append(float(l[6]))
            u.append(float(l[7]))
            
    return x, y, z, nh, flux, xion,u

def particles_from_input_file(input_file):
    x, y, z, n_H, flux, X_ion,u = read_input_file(input_file)
    particles = Particles(len(x))
    particles.x = x | units.parsec
    particles.y = y | units.parsec
    particles.z = z | units.parsec
    particles.rho = n_H | units.amu / units.cm**3
    particles.flux = flux | 1.0e48 / units.s
    particles.xion = X_ion | units.none
    particles.u = u | (units.cm**2/units.s**2)
    return particles

def splitset_from_input_file(input_file):
    x, y, z, n_H, flux, X_ion,u = read_input_file(input_file)
    
    particles = Particles(len(x))
    particles.x = x | units.parsec
    particles.y = y | units.parsec
    particles.z = z | units.parsec
    particles.rho = n_H | units.amu / units.cm**3
    particles.xion = X_ion | units.none
    particles.u = u | (units.cm**2/units.s**2)
    
    a=numpy.where(numpy.array(flux) > 0.)[0]
    src_particles=Particles(len(a))
    src_particles.x = x[a] | units.parsec
    src_particles.y = y[a] | units.parsec
    src_particles.z = z[a] | units.parsec
    src_particles.luminosity = flux[a] | 1.e48*units.s**-1    
    
    return particles,src_particles

