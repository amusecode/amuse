import os
import os.path
import numpy
from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.galactics.interface import GalactICsInterface, GalactICs

# Change the default for some GalactICs(-Interface) keyword arguments:
default_options = dict()
#default_options = dict(redirection = "none")

class GalactICsInterfaceTests(TestWithMPI):
    
    def test1(self):
        print "Testing GalactICsInterface initialization"
        instance = GalactICsInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.set_output_path(instance.get_output_directory()), 0)
        ensure_data_directory_exists(instance.get_output_directory())
        self.assertEqual(instance.set_generate_bulge_flag(False), 0)
        self.assertEqual(instance.set_generate_disk_flag(False), 0)
        self.assertEqual(instance.set_order_of_multipole_expansion(0), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    
    def test2(self):
        print "Testing GalactICsInterface parameters"
        instance = GalactICsInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.set_output_path(os.path.join(instance.get_output_directory(), "test")), 0)
        ensure_data_directory_exists(instance.get_output_directory())
        
        self.assertEquals(instance.set_generate_halo_flag(False), 0)
        self.assertEquals(instance.set_disk_do_center_flag(False), 0)
        self.assertEquals(instance.set_number_of_grid_intervals(50000), 0)
        self.assertEquals(instance.set_disk_random_seed(-1234), 0)
        self.assertEquals(instance.set_halo_outer_radius(250.0), 0)
        self.assertEquals(instance.set_bulge_streaming_fraction(0.4), 0)
        
        self.assertEquals([False, 0], instance.get_generate_halo_flag().values())
        self.assertEquals([False, 0], instance.get_disk_do_center_flag().values())
        self.assertEquals([50000, 0], instance.get_number_of_grid_intervals().values())
        self.assertEquals([-1234, 0], instance.get_disk_random_seed().values())
        self.assertEquals([250.0, 0], instance.get_halo_outer_radius().values())
        self.assertEquals([0.4, 0], instance.get_bulge_streaming_fraction().values())
        
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    
    def slowtest3(self):
        print "Testing GalactICsInterface generate_particles"
        n_particles_halo = 100
        n_particles_bulge = 100
        n_particles_disk = 100
        number_of_particles = n_particles_disk + n_particles_bulge + n_particles_halo
        
        instance = GalactICsInterface(**default_options)
        self.assertEquals(instance.initialize_code(), 0)
        self.assertEquals(instance.set_output_path(instance.get_output_directory()), 0)
        ensure_data_directory_exists(instance.get_output_directory())
        self.assertEquals(instance.set_halo_number_of_particles(n_particles_halo), 0)
        self.assertEquals(instance.set_bulge_number_of_particles(n_particles_bulge), 0)
        self.assertEquals(instance.set_disk_number_of_particles(n_particles_disk), 0)
        self.assertEquals(instance.commit_parameters(), 0)
        
        self.assertEquals(instance.get_number_of_particles_updated().values(), [0, 0])
        self.assertEquals(instance.generate_particles(), 0)
        self.assertEquals(instance.get_number_of_particles_updated().values(), [number_of_particles, 0])
        
        mass_disk, mass_bulge, mass_halo = 26.571852, 14.6317065, 1186.23991
        masses, errors = instance.get_mass(range(number_of_particles))
        self.assertEquals(errors, numpy.zeros(number_of_particles))
        self.assertAlmostRelativeEquals(masses, numpy.concatenate((
            numpy.ones(n_particles_disk)*mass_disk/n_particles_disk, 
            numpy.ones(n_particles_bulge)*mass_bulge/n_particles_bulge, 
            numpy.ones(n_particles_halo)*mass_halo/n_particles_halo
        )), 5)
        
        x_positions, y_positions, z_positions, errors = instance.get_position(range(number_of_particles))
        self.assertEquals(errors, numpy.zeros(number_of_particles))
        self.assertAlmostEquals(numpy.array([numpy.mean(x_positions), numpy.mean(y_positions), 
            numpy.mean(z_positions)]), numpy.array([0.0]*3), 5)
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_positions[:n_particles_disk])), 
            numpy.mean(abs(y_positions[:n_particles_disk])), 
            numpy.mean(abs(z_positions[:n_particles_disk]))]), numpy.array([7.61088372, 7.67542887, 0.38014104]), 5)
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_positions[n_particles_disk:n_particles_disk+n_particles_bulge])), 
            numpy.mean(abs(y_positions[n_particles_disk:n_particles_disk+n_particles_bulge])), 
            numpy.mean(abs(z_positions[n_particles_disk:n_particles_disk+n_particles_bulge]))]), numpy.array([1.09453949, 1.08477648, 0.91855720]), 5)
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_positions[-n_particles_halo:])), 
            numpy.mean(abs(y_positions[-n_particles_halo:])), 
            numpy.mean(abs(z_positions[-n_particles_halo:]))]), numpy.array([63.08369136, 82.63743181, 71.36220391]), 5)

        
        x_velocities, y_velocities, z_velocities, errors = instance.get_velocity(range(number_of_particles))
        self.assertEquals(errors, numpy.zeros(number_of_particles))
        self.assertAlmostEquals(numpy.array([numpy.mean(x_velocities), numpy.mean(y_velocities), 
            numpy.mean(z_velocities)]), numpy.array([0.0]*3))
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_velocities[:n_particles_disk])), 
            numpy.mean(abs(y_velocities[:n_particles_disk])), 
            numpy.mean(abs(z_velocities[:n_particles_disk]))]), numpy.array([1.55797144, 1.55624395, 0.19477551]), 5)
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_velocities[n_particles_disk:])), 
            numpy.mean(abs(y_velocities[n_particles_disk:])), 
            numpy.mean(abs(z_velocities[n_particles_disk:]))]), numpy.array([1.03654369, 1.05683125, 0.97451790]), 5)
        
        self.assertEquals(instance.cleanup_code(), 0)
        instance.stop()
    
    def test4(self):
        print "Testing GalactICsInterface generate_particles"
        number_of_particles_halo = 1000
        instance = GalactICsInterface(**default_options)
        self.assertEquals(instance.initialize_code(), 0)
        self.assertEquals(instance.set_output_path(instance.get_output_directory()), 0)
        ensure_data_directory_exists(instance.get_output_directory())
        self.assertEquals(instance.set_halo_number_of_particles(number_of_particles_halo), 0)
        self.assertEquals(instance.set_generate_bulge_flag(False), 0)
        self.assertEquals(instance.set_generate_disk_flag(False), 0)
        self.assertEquals(instance.set_order_of_multipole_expansion(0), 0)
        self.assertEquals(instance.commit_parameters(), 0)
        
        self.assertEquals(instance.get_number_of_particles_updated().values(), [0, 0])
        self.assertEquals(instance.generate_particles(), 0)
        self.assertEquals(instance.get_number_of_particles_updated().values(), [number_of_particles_halo, 0])

        
        
        masses, errors = instance.get_mass(range(number_of_particles_halo))
        self.assertEquals(errors, numpy.zeros(number_of_particles_halo))
        self.assertAlmostRelativeEquals(masses, numpy.ones(number_of_particles_halo)*masses[0])
        total_mass = masses.sum() 
        
        if ((total_mass-1179.03507)/1179.03507) < 1e-5:
            expected_mean_pos = numpy.array([73.71981114, 79.16010066, 76.47781355])
            expected_mean_vel = numpy.array([0.95465595, 0.90994227, 0.92280032])
        else:
            expected_mean_pos = numpy.array([73.73911318, 79.18082641, 76.49783728])
            expected_mean_vel = numpy.array([0.95453909, 0.90982887, 0.92268653])
            
        x_positions, y_positions, z_positions, errors = instance.get_position(range(number_of_particles_halo))
        self.assertEquals(errors, numpy.zeros(number_of_particles_halo))
        self.assertAlmostEquals(numpy.array([numpy.mean(x_positions), numpy.mean(y_positions), 
            numpy.mean(z_positions)]), numpy.array([0.0]*3), 5)
        self.assertAlmostRelativeEquals(numpy.array([numpy.mean(abs(x_positions)), numpy.mean(abs(y_positions)), 
            numpy.mean(abs(z_positions))]), expected_mean_pos, 3)
        
        x_velocities, y_velocities, z_velocities, errors = instance.get_velocity(range(number_of_particles_halo))
        self.assertEquals(errors, numpy.zeros(number_of_particles_halo))
        self.assertAlmostEquals(numpy.array([numpy.mean(x_velocities), numpy.mean(y_velocities), 
            numpy.mean(z_velocities)]), numpy.array([0.0]*3))
        self.assertAlmostRelativeEquals(numpy.array([numpy.mean(abs(x_velocities)), numpy.mean(abs(y_velocities)), 
            numpy.mean(abs(z_velocities))]), expected_mean_vel, 4)
        
        self.assertEquals(instance.cleanup_code(), 0)
        instance.stop()
    


class GalactICsTests(TestWithMPI):
     
    default_unit_converter = nbody_system.nbody_to_si(1.0 | units.kpc, 1.0e6 | units.MSun)
   
    def test1(self):
        print "Testing GalactICs initialization"
        instance = GalactICs(**default_options)
        instance.initialize_code()
        instance.parameters.generate_bulge_flag = False
        instance.parameters.generate_disk_flag = False
        instance.parameters.order_of_multipole_expansion = 0
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Testing GalactICs parameters (with unit converter)"
        instance = GalactICs(self.default_unit_converter, **default_options)
        instance.initialize_code()
        
        for par, value in [('generate_halo_flag', True), ('generate_disk_flag', True), 
                ('generate_bulge_flag', True), ('halo_do_center_flag', True), 
                ('bulge_do_center_flag', True), ('disk_do_center_flag', True)]:
            self.assertTrue(value is getattr(instance.parameters, par))
            setattr(instance.parameters, par, not value)
            self.assertFalse(value is getattr(instance.parameters, par))
        
        for par, value in [('number_of_grid_intervals', 90000), ('order_of_multipole_expansion', 10), 
                ('number_of_radial_steps_correction_fns_disk_df', 10),
                ('number_of_iterations_disk_df', 50), ('halo_number_of_particles', 200000), 
                ('bulge_number_of_particles', 50000), ('disk_number_of_particles', 100000), 
                ('halo_random_seed', -1), ('bulge_random_seed', -1), ('disk_random_seed', -1)]:
            self.assertEquals(value | units.none, getattr(instance.parameters, par))
            setattr(instance.parameters, par, 1 | units.none)
            self.assertEquals(1 | units.none, getattr(instance.parameters, par))
        
        for par, value in [('halo_outer_radius', 300.0 | nbody_system.length), 
                ('halo_scale_velocity', 3.26331115 | nbody_system.speed), 
                ('halo_scale_radius', 6.06699419 | nbody_system.length), 
                ('halo_truncation_width', 100.0 | nbody_system.length)]:
            self.assertEquals(instance.unit_converter.to_si(value), 
                getattr(instance.parameters, par))
            setattr(instance.parameters, par, 3.0 | value.unit)
            self.assertEquals(instance.unit_converter.to_si(3.0 | value.unit),
                getattr(instance.parameters, par))
        
        self.assertEquals(os.path.join(instance.get_output_directory()) | units.string, instance.parameters.output_directory)
        instance.parameters.output_directory = 'test' | units.string
        self.assertEquals("test" | units.string, instance.parameters.output_directory)
        
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):
        print "Testing GalactICs parameters (nbody units, no converter)"
        instance = GalactICs(**default_options)
        instance.initialize_code()
        
        for par, value in [('halo_outer_radius', 300.0 | nbody_system.length), 
                ('halo_scale_velocity', 3.26331115 | nbody_system.speed), 
                ('halo_scale_radius', 6.06699419 | nbody_system.length), 
                ('halo_truncation_width', 100.0 | nbody_system.length)]:
            self.assertEquals(value, getattr(instance.parameters, par))
            setattr(instance.parameters, par, 3.0 | value.unit)
            self.assertEquals(3.0 | value.unit, getattr(instance.parameters, par))
        instance.cleanup_code()
        instance.stop()
    
    def slowtest4(self):
        print "Testing GalactICs generate_particles"
        n_particles_halo = 100
        n_particles_bulge = 100
        n_particles_disk = 100
        number_of_particles = n_particles_disk + n_particles_bulge + n_particles_halo
        
        instance = GalactICs(**default_options)
        instance.initialize_code()
        instance.parameters.disk_number_of_particles = n_particles_disk
        instance.parameters.bulge_number_of_particles = n_particles_bulge
        instance.parameters.halo_number_of_particles = n_particles_halo
        instance.commit_parameters()
        instance.generate_particles()
        self.assertEquals(len(instance.particles), number_of_particles)
        self.assertAlmostRelativeEquals(instance.particles.total_mass(), 1227.4434685 | nbody_system.mass, 5)
        self.assertAlmostRelativeEquals(instance.particles.kinetic_energy(), 2912.27638811 | nbody_system.energy, 5)
        self.assertAlmostRelativeEquals(instance.particles.potential_energy(G = nbody_system.G), -6318.78987337 | nbody_system.energy, 5)
        self.assertAlmostRelativeEquals(instance.particles.virial_radius(), 119.217247175 | nbody_system.length, 5)
        
        instance.cleanup_code()
        instance.stop()
     
    def test5(self):
        print "Testing GalactICs generate_particles"
        instance = GalactICs(**default_options)
        instance.initialize_code()
        instance.parameters.halo_number_of_particles = 1000
        instance.parameters.generate_bulge_flag = False
        instance.parameters.generate_disk_flag = False
        instance.parameters.order_of_multipole_expansion = 0
        instance.commit_parameters()
        instance.generate_particles()
        self.assertEquals(len(instance.particles), 1000)
        accuracy = 3
        mass_halo = 1179 | nbody_system.mass
        expected_kinetic_energy = 2506 | nbody_system.energy
        self.assertAlmostRelativeEquals(instance.particles.total_mass(), mass_halo, accuracy)
        self.assertAlmostRelativeEquals(instance.particles.kinetic_energy(), expected_kinetic_energy, accuracy)
        instance.cleanup_code()
        instance.stop()
    
    def test6(self):
        print "Testing GalactICs generate_particles: generate multiple sets"
        number_of_particles = 1000
        instance = GalactICs(**default_options)
        instance.initialize_code()
        instance.parameters.halo_number_of_particles = number_of_particles
        instance.parameters.generate_bulge_flag = False
        instance.parameters.generate_disk_flag = False
        instance.parameters.order_of_multipole_expansion = 0
        instance.parameters.halo_random_seed = -1.0
        instance.commit_parameters()
        
        instance.generate_particles()
        set1 = instance.particles.copy()
        self.assertEquals(len(set1), number_of_particles)
        
        instance.generate_particles()
        set2 = instance.particles.copy()
        self.assertEquals(len(set2), number_of_particles)
        # GalactICs' random-number generator is re-seeded with 'halo_random_seed'
        # each time, and the result should be the same:
        for attribute in ["mass", "x", "y", "z", "vx", "vy", "vz"]:
            self.assertEquals(getattr(set1, attribute), getattr(set2, attribute))
        
        instance.parameters.halo_random_seed = -42.0
        instance.generate_particles()
        # halo_random_seed changed: draw a different random set of particles
        set3 = instance.particles.copy()
        self.assertEquals(len(set3), number_of_particles)
        self.assertEquals(set1.mass, set3.mass)
        self.assertRaises(self.failureException, self.assertEquals, set1.x, set3.x)
        self.assertAlmostRelativeEquals(abs(set1.x).median(), abs(set3.x).median(), 1)
        self.assertAlmostRelativeEquals(abs(set1.vy).median(), abs(set3.vy).median(), 1)
        
        instance.cleanup_code()
        instance.stop()
    
    def test7(self):
        print "Testing GalactICs state"
        number_of_particles = 1000

        print "First do everything manually:"
        instance = GalactICs(**default_options)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.parameters.halo_number_of_particles = number_of_particles
        instance.parameters.generate_bulge_flag = False
        instance.parameters.generate_disk_flag = False
        instance.parameters.order_of_multipole_expansion = 0
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.overridden().generate_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'UPDATE')
        instance.invoke_state_change_updated()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        self.assertEquals(len(instance.particles), number_of_particles)
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()

        print "initialize_code(), (re)commit_parameters(), update_particle_set(), " \
            "and cleanup_code() should be called automatically:"
        instance = GalactICs(**default_options)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.halo_number_of_particles = number_of_particles
        instance.parameters.generate_bulge_flag = False
        instance.parameters.generate_disk_flag = False
        instance.parameters.order_of_multipole_expansion = 0
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        self.assertEquals(instance.get_number_of_particles_updated(), 0)
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.parameters.halo_random_seed = -42.0
        self.assertEquals(instance.get_name_of_current_state(), 'PARAMETER_CHANGE_B')
        self.assertEquals(instance.get_number_of_particles_updated(), 0)
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.generate_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        self.assertEquals(len(instance.particles), number_of_particles)
        self.assertEquals(instance.get_number_of_particles_updated(), 0)
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
    
