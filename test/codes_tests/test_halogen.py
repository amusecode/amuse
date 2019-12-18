import os.path
import numpy
import subprocess
from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from amuse.community.halogen.interface import HalogenInterface, Halogen
import amuse.community.halogen

# Change the default for some Halogen(-Interface) keyword arguments:
default_options = dict()
#default_options = dict(redirection = "none")

class HalogenInterfaceTests(TestWithMPI):
    
    def test1(self):
        print("Testing HalogenInterface initialization")
        instance = HalogenInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.set_model_alpha(2.0), 0)
        self.assertEqual(instance.set_model_beta(5.0), 0)
        self.assertEqual(instance.set_model_gamma(0.0), 0)
        self.assertEqual(instance.set_target_number_of_particles(1000), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    
    def test2(self):
        print("Testing HalogenInterface parameters")
        instance = HalogenInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        
        self.assertEqual(instance.set_model_alpha(2.0), 0)
        self.assertEqual(instance.set_model_beta(5.0),  0)
        self.assertEqual(instance.set_model_gamma(0.0), 0)
        self.assertEqual(instance.set_target_number_of_particles(100), 0)
        self.assertEqual(instance.set_random_seed(1), 0)
        self.assertEqual(instance.set_total_mass(9.0),  0)
        self.assertEqual(instance.set_scale_radius(9.0), 0)
        self.assertEqual(instance.set_cutoff_radius(9.0), 0)
        self.assertEqual(instance.set_black_hole_mass(9.0), 0)
        self.assertEqual(instance.set_do_exact_virial_radius_flag(1.0), 0)
        
        self.assertEqual([2.0, 0], list(instance.get_model_alpha().values()))
        self.assertEqual([5.0, 0], list(instance.get_model_beta().values()))
        self.assertEqual([0.0, 0], list(instance.get_model_gamma().values()))
        self.assertEqual([100, 0], list(instance.get_target_number_of_particles().values()))
        self.assertEqual([1, 0], list(instance.get_random_seed().values()))
        self.assertEqual([9.0, 0], list(instance.get_total_mass().values()))
        self.assertEqual([9.0, 0], list(instance.get_scale_radius().values()))
        self.assertEqual([9.0, 0], list(instance.get_cutoff_radius().values()))
        self.assertEqual([9.0, 0], list(instance.get_black_hole_mass().values()))
        self.assertEqual([1.0, 0], list(instance.get_do_exact_virial_radius_flag().values()))
        
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    
    def test3(self):
        print("Testing HalogenInterface generate_particles")
        number_of_particles = 100000
        instance = HalogenInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        self.assertEqual(instance.set_model_alpha(2.0), 0)
        self.assertEqual(instance.set_model_beta(5.0), 0)
        self.assertEqual(instance.set_model_gamma(0.0), 0)
        self.assertEqual(instance.set_target_number_of_particles(number_of_particles), 0)
        self.assertEqual(instance.set_random_seed(1), 0)
        self.assertEqual(instance.set_output_path(instance.get_output_directory()), 0)
        self.assertEqual(instance.commit_parameters(), 0)
        
        self.assertEqual(list(instance.get_number_of_particles_updated().values()), [0, 0])
        self.assertEqual(instance.generate_particles(), 0)
        self.assertEqual(list(instance.get_number_of_particles_updated().values()), [number_of_particles, 0])
        
        masses, errors = instance.get_mass(list(range(number_of_particles)))
        self.assertEqual(errors, numpy.zeros(number_of_particles))
        self.assertEqual(masses, numpy.ones(number_of_particles)/number_of_particles)
        
        x_positions, y_positions, z_positions, errors = instance.get_position(list(range(number_of_particles)))
        self.assertEqual(errors, numpy.zeros(number_of_particles))
        self.assertAlmostEqual(numpy.array([numpy.mean(x_positions), numpy.mean(y_positions), 
            numpy.mean(z_positions)]), numpy.array([0.0]*3))
        self.assertAlmostEqual(numpy.array([numpy.mean(abs(x_positions)), numpy.mean(abs(y_positions)), 
            numpy.mean(abs(z_positions))]), numpy.array([1.0]*3), 1)
        
        x_velocities, y_velocities, z_velocities, errors = instance.get_velocity(list(range(number_of_particles)))
        self.assertEqual(errors, numpy.zeros(number_of_particles))
        self.assertAlmostEqual(numpy.array([numpy.mean(x_velocities), numpy.mean(y_velocities), 
            numpy.mean(z_velocities)]), numpy.array([0.0]*3))
        self.assertAlmostEqual(numpy.array([numpy.mean(abs(x_velocities)), numpy.mean(abs(y_velocities)), 
            numpy.mean(abs(z_velocities))]), numpy.array([0.25]*3), 1)
        
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    
    def test4(self):
        print("Testing HalogenInterface output file name and directory")
        instance = HalogenInterface(**default_options)
        self.assertEqual(instance.initialize_code(), 0)
        
        self.assertEqual(instance.set_model_alpha(2.0), 0)
        self.assertEqual(instance.set_model_beta(5.0), 0)
        self.assertEqual(instance.set_model_gamma(0.0), 0)
        self.assertEqual(instance.set_target_number_of_particles(1000), 0)
        self.assertEqual(instance.set_random_seed(1), 0)
        self.assertEqual(instance.set_write_output_flag(1.0), 0)
        
        self.assertEqual(list(instance.get_output_basename().values()), ["halogen", 0])
        self.assertEqual(list(instance.get_output_path().values()), ["./", 0])
        
        self.assertEqual(instance.set_output_basename("oops_this_output_basename_is_way_too_long"*2), -1)
        self.assertEqual(list(instance.get_output_basename().values()), ["halogen", 0])
        self.assertEqual(instance.set_output_path("/oops/this/output/path/has/way/too/many/subdirs"*6), -1)
        self.assertEqual(list(instance.get_output_path().values()), ["./", 0])
        
        self.assertEqual(instance.set_output_basename("test"), 0)
        self.assertEqual(instance.set_output_path(instance.get_output_directory()), 0)
        
        self.assertEqual(list(instance.get_output_basename().values()), ["test", 0])
        self.assertEqual(list(instance.get_output_path().values()), 
            [os.path.join(instance.get_output_directory(), ""), 0])
        
        self.assertEqual(instance.commit_parameters(), 0)
        outputfile = os.path.join(instance.get_output_directory(), "test.IC.ascii")
        if os.path.exists(outputfile):
            os.remove(outputfile)
        self.assertEqual(instance.generate_particles(), 0)
        self.assertTrue(os.path.exists(outputfile))
        
        
        halogen4muse_path = os.path.join(os.path.dirname(amuse.community.halogen.__file__), 'src', 'halogen4muse')
        
        if not os.path.exists(halogen4muse_path) or not os.access(halogen4muse_path, os.X_OK):
            return
            
        process = subprocess.Popen([
            halogen4muse_path, 
            '-a', '2', 
            '-b', '5', 
            '-c', '0', 
            '-N', '1000', 
            '-name', 'test_stand_alone', 
            '-randomseed', '1'
            ]
            , cwd = instance.get_output_directory()
            , stdout = subprocess.PIPE
            , stderr = subprocess.PIPE
        )
        process.communicate()
        
        self.compare_files(
            os.path.join(instance.get_output_directory(), "test.IC.ascii"),
            os.path.join(instance.get_output_directory(), "test_stand_alone.IC.ascii")
        )
        stdoutput = subprocess.Popen(["diff", "test.out", "test_stand_alone.out"], 
            cwd = instance.get_output_directory(), stdout = subprocess.PIPE).communicate()[0]
        stdoutput=stdoutput.decode()
        self.assertTrue("< N/A (executed by AMUSE)" in stdoutput)
        self.assertTrue("halogen4muse -a 2 -b 5 -c 0 -N 1000 -name test_stand_alone -randomseed 1" in stdoutput)
        
        self.assertEqual(instance.cleanup_code(), 0)
        instance.stop()
    
    def compare_files(self, filename1, filename2):
        with open(filename1) as stream1:
            lines1 = stream1.readlines()
        with open(filename2) as stream2:
            lines2 = stream2.readlines()

        self.assertTrue(len(lines1), len(lines2))
        for line1, line2 in zip(lines1, lines2):
            columns1 = line1.split()
            columns2 = line2.split()
            self.assertTrue(len(columns1), len(columns2))
            if (len(columns1)) == 0:
                continue
            self.assertEqual(int(columns1[0]), int(columns2[0]))
            for column1, column2 in zip(columns1[1:], columns2[1:]):
                self.assertAlmostRelativeEquals(
                    float(column1),
                    float(column2),
                    8
                )
                
            
class HalogenTests(TestWithMPI):
    
    default_unit_converter = nbody_system.nbody_to_si(1.0 | units.kpc, 1.0e6 | units.MSun)
    
    def test1(self):
        print("Testing Halogen initialization")
        instance = Halogen(**default_options)
        instance.initialize_code()
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print("Testing Halogen parameters (with unit converter)")
        instance = Halogen(self.default_unit_converter, **default_options)
        instance.initialize_code()
        
        for par, value in [('do_exact_virial_radius_flag', False), 
                ('outputgridr_flag', False), ('outputgriddf_flag', False), 
                ('write_output_flag', False)]:
            self.assertTrue(value is getattr(instance.parameters, par))
            setattr(instance.parameters, par, not value)
            self.assertFalse(value is getattr(instance.parameters, par))
        
        for par, value in [('alpha', -1.0), ('beta', -1.0), ('gamma', -1.0),
                ('number_of_particles', -1), ('random_seed', 42)]:
            self.assertEqual(value, getattr(instance.parameters, par))
            setattr(instance.parameters, par, 1)
            self.assertEqual(1, getattr(instance.parameters, par))
        
        for par, value in [('total_mass', 1.0 | nbody_system.mass), 
                ('scale_radius',    1.0 | nbody_system.length), 
                ('cutoff_radius',  -1.0 | nbody_system.length), 
                ('black_hole_mass', 0.0 | nbody_system.mass)]:
            self.assertEqual(instance.unit_converter.to_si(value), 
                getattr(instance.parameters, par))
            setattr(instance.parameters, par, 3.0 | value.unit)
            self.assertEqual(instance.unit_converter.to_si(3.0 | value.unit),
                getattr(instance.parameters, par))
        
        for par, value in [('output_directory', os.path.join(instance.get_output_directory(), "")), 
                ('output_basename', "halogen")]:
            self.assertEqual(value, getattr(instance.parameters, par))
            setattr(instance.parameters, par, 'test/')
            self.assertEqual("test/", getattr(instance.parameters, par))
        
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):
        print("Testing Halogen parameters (nbody units, no converter)")
        instance = Halogen(**default_options)
        instance.initialize_code()
        
        for par, value in [('total_mass', 1.0 | nbody_system.mass), 
                ('scale_radius',    1.0 | nbody_system.length), 
                ('cutoff_radius',  -1.0 | nbody_system.length), 
                ('black_hole_mass', 0.0 | nbody_system.mass)]:
            self.assertEqual(value, getattr(instance.parameters, par))
            setattr(instance.parameters, par, 3.0 | value.unit)
            self.assertEqual(3.0 | value.unit, getattr(instance.parameters, par))
        instance.cleanup_code()
        instance.stop()
    
    def test4(self):
        print("Testing Halogen generate_particles")
        number_of_particles = 100
        instance = Halogen(**default_options)
        instance.initialize_code()
        instance.parameters.alpha = 2.0
        instance.parameters.beta  = 5.0
        instance.parameters.gamma = 0.0
        instance.parameters.number_of_particles = number_of_particles
        instance.parameters.random_seed = 1
        instance.commit_parameters()
        instance.generate_particles()
        self.assertEqual(len(instance.particles), number_of_particles)
        self.assertAlmostEqual(instance.particles.total_mass(), 1.0 | nbody_system.mass)
        self.assertAlmostEqual(instance.particles.kinetic_energy(), 
            0.17345836639 | nbody_system.energy)
        self.assertAlmostEqual(instance.particles.potential_energy(G = nbody_system.G), 
            -0.310395778644 | nbody_system.energy)
        self.assertAlmostEqual(instance.particles.virial_radius(), 
            1.61084664935 | nbody_system.length)
        
        instance.cleanup_code()
        instance.stop()
    
    def test5(self):
        print("Testing Halogen generate_particles: generate multiple sets")
        number_of_particles = 1000
        instance = Halogen(**default_options)
        instance.initialize_code()
        instance.parameters.alpha = 2.0
        instance.parameters.beta  = 5.0
        instance.parameters.gamma = 0.0
        instance.parameters.number_of_particles = number_of_particles
        instance.parameters.random_seed = 1
        instance.commit_parameters()
        
        instance.generate_particles()
        set1 = instance.particles.copy()
        self.assertEqual(len(set1), number_of_particles)
        
        instance.parameters.random_seed = 1
        instance.generate_particles()
        set2 = instance.particles.copy()
        self.assertEqual(len(set2), number_of_particles)
        # Since a (any would do!) parameter was changed, recommit_parameters was
        # called, re-seeding, and the result should be the same:
        for attribute in ["mass", "x", "y", "z", "vx", "vy", "vz"]:
            self.assertEqual(getattr(set1, attribute), getattr(set2, attribute))
        
        instance.generate_particles()
        # No parameter change: draw the next random set of particles
        set3 = instance.particles.copy()
        self.assertEqual(len(set3), number_of_particles)
        self.assertEqual(set1.mass, set3.mass)
        self.assertRaises(self.failureException, self.assertEqual, set1.x, set3.x)
        self.assertIsOfOrder(abs(set1.x).median(), abs(set3.x).median(), 1)
        self.assertAlmostEqual(abs(set1.vy).median(), abs(set3.vy).median(), 1)
        
        instance.cleanup_code()
        instance.stop()
    
    def test6(self):
        print("Testing Halogen state")
        number_of_particles = 1000

        print("First do everything manually:")
        instance = Halogen(**default_options)
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.parameters.alpha = 2.0
        instance.parameters.beta  = 5.0
        instance.parameters.gamma = 0.0
        instance.parameters.number_of_particles = number_of_particles
        instance.parameters.random_seed = 1
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
        instance = Halogen(**default_options)
        self.assertEqual(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.alpha = 2.0
        instance.parameters.beta  = 5.0
        instance.parameters.gamma = 0.0
        instance.parameters.number_of_particles = number_of_particles
        instance.parameters.random_seed = 1
        self.assertEqual(instance.get_name_of_current_state(), 'INITIALIZED')
        self.assertEqual(instance.get_number_of_particles_updated(), 0)
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.parameters.random_seed = 2
        self.assertEqual(instance.get_name_of_current_state(), 'CHANGE_PARAMETERS_EDIT')
        self.assertEqual(instance.get_number_of_particles_updated(), 0)
        self.assertEqual(instance.get_name_of_current_state(), 'EDIT')
        instance.generate_particles()
        self.assertEqual(instance.get_name_of_current_state(), 'RUN')
        self.assertEqual(len(instance.particles), number_of_particles)
        self.assertEqual(instance.get_number_of_particles_updated(), 0)
        instance.stop()
        self.assertEqual(instance.get_name_of_current_state(), 'STOPPED')
    
    def test7(self):
        print("Testing Halogen error handling")
        number_of_particles = 1000
        instance = Halogen(**default_options)
        instance.initialize_code()
        self.assertRaises(exceptions.AmuseException, instance.commit_parameters, expected_message = 
            "Error when calling 'commit_parameters' of a 'Halogen', errorcode is -2, error is "
            "'Missing or bad parameter for halo (see amuse/community/halogen/src/doc for details on required parameters).'")
        instance.parameters.alpha = 2.0
        instance.parameters.beta  = 5.0
        instance.parameters.gamma = 5.0
        instance.parameters.number_of_particles = number_of_particles
        instance.parameters.random_seed = 1
        self.assertRaises(exceptions.AmuseException, instance.commit_parameters, expected_message = 
            "Error when calling 'commit_parameters' of a 'Halogen', errorcode is -2, error is "
            "'Missing or bad parameter for halo (see amuse/community/halogen/src/doc for details on required parameters).'")
        instance.parameters.gamma = -0.5
        self.assertRaises(exceptions.AmuseException, instance.commit_parameters, expected_message = 
            "Error when calling 'commit_parameters' of a 'Halogen', errorcode is -2, error is "
            "'Missing or bad parameter for halo (see amuse/community/halogen/src/doc for details on required parameters).'")
        instance.parameters.gamma = 0.0
        instance.parameters.beta  = 2.0
        self.assertRaises(exceptions.AmuseException, instance.commit_parameters, expected_message = 
            "Error when calling 'commit_parameters' of a 'Halogen', errorcode is -2, error is "
            "'Missing or bad parameter for halo (see amuse/community/halogen/src/doc for details on required parameters).'")
        instance.parameters.beta  = 5.0
        instance.commit_parameters()
        
        instance.cleanup_code()
        instance.stop()
    

