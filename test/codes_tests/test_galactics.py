import os
import os.path
import numpy
import platform
from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.galactics.interface import GalactICsInterface, GalactICs

# Change the default for some GalactICs(-Interface) keyword arguments:
default_options = dict()
#default_options = dict(redirection = "none")

class GalactICsInterfaceTests(TestWithMPI):
    
    def test1(self):
        print("Testing GalactICsInterface initialization")
        instance = GalactICsInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.set_output_path(instance.get_output_directory()), 0)
        self.assertEqual(instance.set_generate_bulge_flag(False), 0)
        self.assertEqual(instance.set_generate_disk_flag(False), 0)
        self.assertEqual(instance.set_order_of_multipole_expansion(0), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    
    def test2(self):
        print("Testing GalactICsInterface parameters")
        instance = GalactICsInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.set_output_path(os.path.join(instance.get_output_directory(), "test")), 0)
        
        self.assertEqual(instance.set_generate_halo_flag(False), 0)
        self.assertEqual(instance.set_disk_do_center_flag(False), 0)
        self.assertEqual(instance.set_number_of_grid_intervals(50000), 0)
        self.assertEqual(instance.set_disk_random_seed(-1234), 0)
        self.assertEqual(instance.set_halo_outer_radius(250.0), 0)
        self.assertEqual(instance.set_bulge_streaming_fraction(0.4), 0)
        
        self.assertEqual([False, 0], list(instance.get_generate_halo_flag().values()))
        self.assertEqual([False, 0], list(instance.get_disk_do_center_flag().values()))
        self.assertEqual([50000, 0], list(instance.get_number_of_grid_intervals().values()))
        self.assertEqual([-1234, 0], list(instance.get_disk_random_seed().values()))
        self.assertEqual([250.0, 0], list(instance.get_halo_outer_radius().values()))
        self.assertEqual([0.4, 0], list(instance.get_bulge_streaming_fraction().values()))
        
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    
    def slowtest3(self):
        print("Testing GalactICsInterface generate_particles")
        n_particles_halo = 100
        n_particles_bulge = 100
        n_particles_disk = 100
        number_of_particles = n_particles_disk + n_particles_bulge + n_particles_halo
        
        instance = GalactICsInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.set_output_path(instance.get_output_directory()), 0)
        self.assertEqual(instance.set_halo_number_of_particles(n_particles_halo), 0)
        self.assertEqual(instance.set_bulge_number_of_particles(n_particles_bulge), 0)
        self.assertEqual(instance.set_disk_number_of_particles(n_particles_disk), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        
        self.assertEqual(list(instance.get_number_of_particles_updated().values()), [0, 0])
        self.assertEqual(instance.generate_particles(), 0)
        self.assertEqual(list(instance.get_number_of_particles_updated().values()), [number_of_particles, 0])
        
        mass_disk, mass_bulge, mass_halo = 26.578816771507263, 14.632800221443176, 1184.2350006103516
        masses, errors = instance.get_mass(list(range(number_of_particles)))
        self.assertEqual(errors, numpy.zeros(number_of_particles))
        self.assertAlmostRelativeEquals(masses, numpy.concatenate((
            numpy.ones(n_particles_disk)*mass_disk/n_particles_disk, 
            numpy.ones(n_particles_bulge)*mass_bulge/n_particles_bulge,
            numpy.ones(n_particles_halo)*mass_halo/n_particles_halo,
        )), 3)
        
        x_positions, y_positions, z_positions, errors = instance.get_position(list(range(number_of_particles)))
        self.assertEqual(errors, numpy.zeros(number_of_particles))
        self.assertAlmostEqual(numpy.array([numpy.mean(x_positions), numpy.mean(y_positions), 
            numpy.mean(z_positions)]), numpy.array([0.0]*3), 5)
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_positions[:n_particles_disk])), 
            numpy.mean(abs(y_positions[:n_particles_disk])), 
            numpy.mean(abs(z_positions[:n_particles_disk]))]), 
            numpy.array([7.3994484072923656, 7.1570298135280606, 0.33854196755215527]), 3)
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_positions[n_particles_disk:n_particles_disk+n_particles_bulge])), 
            numpy.mean(abs(y_positions[n_particles_disk:n_particles_disk+n_particles_bulge])), 
            numpy.mean(abs(z_positions[n_particles_disk:n_particles_disk+n_particles_bulge]))]), 
            numpy.array([1.244429082274437,1.1639373835548759 , 0.8550614269822836]), 3)
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_positions[-n_particles_halo:])), 
            numpy.mean(abs(y_positions[-n_particles_halo:])), 
            numpy.mean(abs(z_positions[-n_particles_halo:]))]), 
            numpy.array([94.242819476127622,88.41320479869843 , 85.234394512176507]), 3)

        
        x_velocities, y_velocities, z_velocities, errors = instance.get_velocity(list(range(number_of_particles)))
        self.assertEqual(errors, numpy.zeros(number_of_particles))
        self.assertAlmostEqual(numpy.array([numpy.mean(x_velocities), numpy.mean(y_velocities), 
            numpy.mean(z_velocities)]), numpy.array([0.0]*3))
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_velocities[:n_particles_disk])), 
            numpy.mean(abs(y_velocities[:n_particles_disk])), 
            numpy.mean(abs(z_velocities[:n_particles_disk]))]), 
            numpy.array([1.5026254250109197, 1.5649469271302223, 0.20230436498299242]), 5)
        self.assertAlmostRelativeEquals(numpy.array([
            numpy.mean(abs(x_velocities[n_particles_disk:])), 
            numpy.mean(abs(y_velocities[n_particles_disk:])), 
            numpy.mean(abs(z_velocities[n_particles_disk:]))]),
            numpy.array([0.99470628838986164,0.95913934175856408 , 0.9359876788407564]), 5)
        
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    
    def test4(self):
        print("Testing GalactICsInterface generate_particles")
        number_of_particles_halo = 1000
        instance = GalactICsInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.set_output_path(instance.get_output_directory()), 0)
        self.assertEqual(instance.set_halo_number_of_particles(number_of_particles_halo), 0)
        self.assertEqual(instance.set_generate_bulge_flag(False), 0)
        self.assertEqual(instance.set_generate_disk_flag(False), 0)
        self.assertEqual(instance.set_order_of_multipole_expansion(0), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        
        self.assertEqual(list(instance.get_number_of_particles_updated().values()), [0, 0])
        self.assertEqual(instance.generate_particles(), 0)
        self.assertEqual(list(instance.get_number_of_particles_updated().values()), [number_of_particles_halo, 0])

        
        
        masses, errors = instance.get_mass(list(range(number_of_particles_halo)))
        self.assertEqual(errors, numpy.zeros(number_of_particles_halo))
        self.assertAlmostRelativeEquals(masses, numpy.ones(number_of_particles_halo)*masses[0])
        total_mass = masses.sum() 
        
        if platform.processor() == 'ppc64le':
            # on ppc64le, the model generation has small differences from intel
            # change expected pos
            expected_mean_pos = numpy.array([73.5628, 76.251034, 75.53434])
        else:
            expected_mean_pos = numpy.array([73.768384103536604, 76.03533643054962, 75.176319462463255])

        expected_mean_vel = numpy.array([0.92904859858192501, 0.94953939936682585, 0.92897711758688095])
            
        x_positions, y_positions, z_positions, errors = instance.get_position(list(range(number_of_particles_halo)))
        self.assertEqual(errors, numpy.zeros(number_of_particles_halo))
        self.assertAlmostEqual(numpy.array([numpy.mean(x_positions), numpy.mean(y_positions), 
            numpy.mean(z_positions)]), numpy.array([0.0]*3), 5)
        self.assertAlmostRelativeEquals(numpy.array([numpy.mean(abs(x_positions)), numpy.mean(abs(y_positions)), 
            numpy.mean(abs(z_positions))]), expected_mean_pos, 3)
        
        x_velocities, y_velocities, z_velocities, errors = instance.get_velocity(list(range(number_of_particles_halo)))
        self.assertEqual(errors, numpy.zeros(number_of_particles_halo))
        self.assertAlmostEqual(numpy.array([numpy.mean(x_velocities), numpy.mean(y_velocities), 
            numpy.mean(z_velocities)]), numpy.array([0.0]*3))
        self.assertAlmostRelativeEquals(numpy.array([numpy.mean(abs(x_velocities)), numpy.mean(abs(y_velocities)), 
            numpy.mean(abs(z_velocities))]), expected_mean_vel, 2)
        
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    


class GalactICsTests(TestWithMPI):
     
    default_unit_converter = nbody_system.nbody_to_si(1.0 | units.kpc, 1.0e6 | units.MSun)
   
    def test1(self):
        print("Testing GalactICs initialization")
        instance = GalactICs(**default_options)
        instance.initialize_code()
        instance.parameters.generate_bulge_flag = False
        instance.parameters.generate_disk_flag = False
        instance.parameters.order_of_multipole_expansion = 0
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print("Testing GalactICs parameters (with unit converter)")
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
            self.assertEqual(value, getattr(instance.parameters, par))
            setattr(instance.parameters, par, 1)
            self.assertEqual(1, getattr(instance.parameters, par))
        
        for par, value in [('halo_outer_radius', 300.0 | nbody_system.length), 
                ('halo_scale_velocity', 3.26331115 | nbody_system.speed), 
                ('halo_scale_radius', 6.06699419 | nbody_system.length), 
                ('halo_truncation_width', 100.0 | nbody_system.length)]:
            self.assertEqual(instance.unit_converter.to_si(value), 
                getattr(instance.parameters, par))
            setattr(instance.parameters, par, 3.0 | value.unit)
            self.assertEqual(instance.unit_converter.to_si(3.0 | value.unit),
                getattr(instance.parameters, par))
        
        self.assertEqual(os.path.join(instance.get_output_directory()), instance.parameters.output_directory)
        instance.parameters.output_directory = 'test'
        self.assertEqual("test", instance.parameters.output_directory)
        
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):
        print("Testing GalactICs parameters (nbody units, no converter)")
        instance = GalactICs(**default_options)
        instance.initialize_code()
        
        for par, value in [('halo_outer_radius', 300.0 | nbody_system.length), 
                ('halo_scale_velocity', 3.26331115 | nbody_system.speed), 
                ('halo_scale_radius', 6.06699419 | nbody_system.length), 
                ('halo_truncation_width', 100.0 | nbody_system.length)]:
            self.assertEqual(value, getattr(instance.parameters, par))
            setattr(instance.parameters, par, 3.0 | value.unit)
            self.assertEqual(3.0 | value.unit, getattr(instance.parameters, par))
        instance.cleanup_code()
        instance.stop()
    
    def slowtest4(self):
        print("Testing GalactICs generate_particles")
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
        self.assertEqual(len(instance.particles), number_of_particles)
        self.assertAlmostRelativeEquals(instance.particles.total_mass(), 1225.4466176 | nbody_system.mass, 3)
        self.assertAlmostRelativeEquals(instance.particles.kinetic_energy(), 2564.69894361 | nbody_system.energy, 3)
        self.assertAlmostRelativeEquals(instance.particles.potential_energy(G = nbody_system.G), -4531.58416742 | nbody_system.energy, 3)
        self.assertAlmostRelativeEquals(instance.particles.virial_radius(), 165.694750127 | nbody_system.length, 3)
        
        instance.cleanup_code()
        instance.stop()
     
    def test5(self):
        print("Testing GalactICs generate_particles")
        instance = GalactICs(**default_options)
        instance.initialize_code()
        instance.parameters.halo_number_of_particles = 1000
        instance.parameters.generate_bulge_flag = False
        instance.parameters.generate_disk_flag = False
        instance.parameters.order_of_multipole_expansion = 0
        instance.commit_parameters()
        instance.generate_particles()
        self.assertEqual(len(instance.particles), 1000)
        accuracy = 3
        mass_halo = 1178.89297009 | nbody_system.mass
        expected_kinetic_energy = 2418.49730735 | nbody_system.energy
        self.assertAlmostRelativeEquals(instance.particles.total_mass(), mass_halo, accuracy)
        self.assertAlmostRelativeEquals(instance.particles.kinetic_energy(), expected_kinetic_energy, accuracy)
        
        self.assertEqual(len(instance.halo_particles),1000)
        self.assertEqual(len(instance.disk_particles),0)
        self.assertEqual(len(instance.bulge_particles),0)
        
        
        instance.cleanup_code()
        instance.stop()
    
    def test6(self):
        print("Testing GalactICs generate_particles: generate multiple sets")
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
        self.assertEqual(len(set1), number_of_particles)
        
        instance.generate_particles()
        set2 = instance.particles.copy()
        self.assertEqual(len(set2), number_of_particles)
        # GalactICs' random-number generator is re-seeded with 'halo_random_seed'
        # each time, and the result should be the same:
        for attribute in ["mass", "x", "y", "z", "vx", "vy", "vz"]:
            self.assertEqual(getattr(set1, attribute), getattr(set2, attribute))
        
        instance.parameters.halo_random_seed = -42.0
        instance.generate_particles()
        # halo_random_seed changed: draw a different random set of particles
        set3 = instance.particles.copy()
        self.assertEqual(len(set3), number_of_particles)
        self.assertEqual(set1.mass, set3.mass)
        self.assertRaises(self.failureException, self.assertEqual, set1.x, set3.x)
        self.assertAlmostRelativeEquals(abs(set1.x).median(), abs(set3.x).median(), 1)
        self.assertAlmostRelativeEquals(abs(set1.vy).median(), abs(set3.vy).median(), 1)
        
        instance.cleanup_code()
        instance.stop()
    
    def test7(self):
        print("Testing GalactICs state")
        number_of_particles = 1000

        print("First do everything manually:")
        instance = GalactICs(**default_options)
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.parameters.halo_number_of_particles = number_of_particles
        instance.parameters.generate_bulge_flag = False
        instance.parameters.generate_disk_flag = False
        instance.parameters.order_of_multipole_expansion = 0
        instance.commit_parameters()
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.overridden().generate_particles()
        self.assertEqual(instance.get_name_of_current_state(), 'UPDATE')
        instance.invoke_state_change_updated()
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        self.assertEqual(len(instance.particles), number_of_particles)
        instance.cleanup_code()
        self.assertEqual(instance.get_name_of_current_state(), 'END')
        instance.stop()

        print("initialize_code(), (re)commit_parameters(), update_particle_set(), " \
            "and cleanup_code() should be called automatically:")
        instance = GalactICs(**default_options)
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.halo_number_of_particles = number_of_particles
        instance.parameters.generate_bulge_flag = False
        instance.parameters.generate_disk_flag = False
        instance.parameters.order_of_multipole_expansion = 0
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        self.assertEqual(instance.get_number_of_particles_updated(), 0)
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.parameters.halo_random_seed = -42.0
        self.assertEqual(instance.get_name_of_current_state(), 'CHANGE_PARAMETERS_EDIT')
        self.assertEqual(instance.get_number_of_particles_updated(), 0)
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.generate_particles()
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        self.assertEqual(len(instance.particles), number_of_particles)
        self.assertEqual(instance.get_number_of_particles_updated(), 0)
        instance.stop()
        self.assertEqual(instance.get_name_of_current_state(), 'STOPPED')
    
