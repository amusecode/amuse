import numpy
from amuse.test.amusetest import TestWithMPI
from amuse.community.gadget2.interface import Gadget2Interface, Gadget2
from amuse.ext.plummer import MakePlummerModel, new_plummer_sphere
from amuse.ext.evrard_test import MakeEvrardTest, new_evrard_gas_sphere
from amuse.ext.spherical_model import new_uniform_spherical_particle_distribution

from amuse.support.exceptions import AmuseException
from amuse.support.units import nbody_system
from amuse.support.units import generic_unit_converter
from amuse.support.units import generic_unit_system
from amuse.support.units import units
from amuse.support.data import core
from amuse.support.codes import channel

# Change the default for some Gadget2(-Interface) keyword arguments:
default_options = dict(number_of_workers=2)
#default_options = dict(number_of_workers=2, redirection="none")

# ... but never use (number_of_workers>1) for tests with only a few particles:
few_particles_default_options = dict()
#few_particles_default_options = dict(redirection="none")

class TestGadget2Interface(TestWithMPI):

    def test1(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test2(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals([1.989000e43, 0], instance.get_unit_mass().values())
        self.assertEquals([3.085678e21, 0], instance.get_unit_length().values())
        self.assertEquals([3.085678e16, 0], instance.get_unit_time().values())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test3(self):
        instance = Gadget2Interface(**few_particles_default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals([1, 0], instance.new_dm_particle(0.01,  1, 0, 0,  0, 1, 0).values())
        self.assertEquals([2, 0], instance.new_dm_particle(0.01, -1, 0, 0,  0,-1, 0).values())
        self.assertEquals(-1, instance.get_index_of_first_particle()['__result'])
        self.assertEquals(0, instance.commit_particles())
        self.assertEquals([1, 0], instance.get_index_of_first_particle().values())
        self.assertEquals([2, 1], instance.get_index_of_next_particle(1).values())
        self.assertEquals(-1, instance.get_index_of_next_particle(2)['__result'])
        self.assertEquals(0, instance.evolve(0.01))
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test4(self):
        instance = Gadget2Interface(**few_particles_default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals([1, 0], instance.new_dm_particle(0.01,  1, 0, 0,  0, 1, 0).values())
        self.assertEquals([2, 0], instance.new_dm_particle(0.02, -1, 0, 0,  0,-1, 0).values())
        self.assertEquals(-3, instance.get_mass(1)['__result'])
        self.assertEquals(0, instance.commit_particles())

        mass, result = instance.get_mass(1)
        self.assertAlmostEquals(0.01,mass)
        self.assertEquals(0,result)
        mass, result = instance.get_mass(2)
        self.assertAlmostEquals(0.02,mass)
        self.assertEquals(0,result)
        self.assertEquals(-3, instance.get_mass(3)['__result'])
        for result,expected in zip(instance.get_position(1),[1,0,0, 0]):
            self.assertAlmostEquals(result,expected)
        for result,expected in zip(instance.get_position(2),[-1,0,0, 0]):
            self.assertAlmostEquals(result,expected)
        for result,expected in zip(instance.get_velocity(1),[0,1,0, 0]):
            self.assertAlmostEquals(result,expected)
        for result,expected in zip(instance.get_velocity(2),[0,-1,0, 0]):
            self.assertAlmostEquals(result,expected)

        self.assertEquals(0, instance.evolve(0.01))
        mass, result = instance.get_mass(1)
        self.assertAlmostEquals(0.01,mass)
        self.assertEquals(0,result)
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test5(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        target_number_of_particles = 100
        evrard = MakeEvrardTest(target_number_of_particles)
        mass, x,y,z, vx,vy,vz, u = evrard.new_model()
        number_of_particles = len(mass)
        indices, results = instance.new_sph_particle(mass,x,y,z,vx,vy,vz,u)
        self.assertEquals([0 for i in range(number_of_particles)], list(results))
        self.assertEquals([i+1 for i in range(number_of_particles)], list(indices))
        self.assertEquals(0, instance.commit_particles())
        self.assertEquals(0, instance.evolve(0.00005))
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test6(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        target_number_of_particles = 100
        evrard = MakeEvrardTest(target_number_of_particles)
        mass, x,y,z, vx,vy,vz, u = evrard.new_model()
        number_of_particles = len(mass)
        first_half = number_of_particles/2
        indices, results = instance.new_sph_particle(mass[:first_half],x[:first_half],y[:first_half],z[:first_half],
            vx[:first_half],vy[:first_half],vz[:first_half],u[:first_half])
        self.assertEquals([0 for i in range(first_half)], list(results))
        self.assertEquals([i+1 for i in range(first_half)], list(indices))
        self.assertEquals([number_of_particles/2+1, 0], instance.new_dm_particle(0.01,  1, 0, 0,  0, 1, 0).values())
        self.assertEquals([number_of_particles/2+2, 0], instance.new_dm_particle(0.02, -1, 0, 0,  0,-1, 0).values())
        indices, results = instance.new_sph_particle(mass[first_half:],x[first_half:],y[first_half:],z[first_half:],
            vx[first_half:],vy[first_half:],vz[first_half:],u[first_half:])
        self.assertEquals([0 for i in range(number_of_particles-first_half)], list(results))
        self.assertEquals([first_half+i+1+2 for i in range(number_of_particles-first_half)], list(indices))
        self.assertEquals(0, instance.commit_particles())

        mass_list = [x for sublist in [mass[:first_half], [0.01, 0.02], mass[first_half:]] for x in sublist]
        first_index, result = instance.get_index_of_first_particle().values()
        self.assertEquals([1, 0], [first_index, result])
        for i in range(first_index, number_of_particles+2):
            index, result = instance.get_index_of_next_particle(i)
            self.assertEquals([i+1, 0 if (i < number_of_particles+1) else 1], [index, result])
            mass, result = instance.get_mass(index)
            self.assertEquals(0, result)
            self.assertAlmostEquals(mass_list[i], mass)

        self.assertEquals(0, instance.evolve(0.0001))
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test7(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        target_number_of_particles = 100
        evrard = MakeEvrardTest(target_number_of_particles)
        mass, x,y,z, vx,vy,vz, u = evrard.new_model()
        number_of_particles = len(mass)
        indices, results = instance.new_sph_particle(mass,x,y,z,vx,vy,vz,u)
        self.assertEquals([0 for i in range(number_of_particles)], list(results))
        self.assertEquals([i+1 for i in range(number_of_particles)], list(indices))
        self.assertEquals([number_of_particles+1, 0], instance.new_dm_particle(0.01,  1, 0, 0,  0, 1, 0).values())
        self.assertEquals([number_of_particles+2, 0], instance.new_dm_particle(0.02, -1, 0, 0,  0,-1, 0).values())
        self.assertEquals(0, instance.commit_particles())
        self.assertEquals(0, instance.evolve(0.00005))
        self.assertEquals(0, instance.delete_particle(number_of_particles-1))
        self.assertEquals(0, instance.delete_particle(number_of_particles+1))
        self.assertEquals(-3, instance.delete_particle(number_of_particles-1))
        indices, results = instance.new_sph_particle(mass,x,y,z,vx,vy,vz,u)
        self.assertEquals([2*number_of_particles+3, 0], instance.new_dm_particle(0.02, -1, 0, 0,  0,-1, 0).values())
        self.assertEquals(0, instance.recommit_particles())
        mass_list = [x for sublist in [mass[:number_of_particles-2], [mass[number_of_particles-1], 0.02],
            mass, [0.02]] for x in sublist]
        index_list = [x for sublist in [range(1,number_of_particles-1), [number_of_particles],
            range(number_of_particles+2,2*number_of_particles+4)] for x in sublist]
        index, result = instance.get_index_of_first_particle().values()
        self.assertEquals([index_list[0], 0], [index, result])
        for i in range(1, 2*number_of_particles+1):
            index, result = instance.get_index_of_next_particle(index)
            self.assertEquals([index_list[i], 0 if (i < 2*number_of_particles) else 1], [index, result])
            mass, result = instance.get_mass(index)
            self.assertEquals(0, result)
            self.assertAlmostEquals(mass_list[i], mass)

        self.assertEquals(0, instance.evolve(0.00008))
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test8(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        for i in range(1,10):
            self.assertEquals([i, 0], instance.new_dm_particle(i*0.01,  i, -i, 0,  -i*10, i*10, 0).values())
        self.assertEquals(0, instance.commit_particles())
        
        for i in range(1,10):
            mass, result = instance.get_mass(i)
            self.assertAlmostEquals(0.01*i,mass)
            self.assertEquals(0,result)
            
            x, y, z, result = instance.get_position(i)
            self.assertAlmostEquals(i, x)
            self.assertAlmostEquals(-i, y)
            self.assertAlmostEquals(0, z)
            self.assertEquals(0,result)
            
            vx, vy, vz, result = instance.get_velocity(i)
            self.assertAlmostEquals(-i*10, vx)
            self.assertAlmostEquals(i*10, vy)
            self.assertAlmostEquals(0, vz)
            self.assertEquals(0,result)
            
            self.assertEquals(0, instance.set_mass(i, i*0.1))
            mass, result = instance.get_mass(i)
            self.assertAlmostEquals(0.1*i,mass)
            self.assertEquals(0,result)
            
            self.assertEquals(0, instance.set_position(i, 2*i, -2*i, 0))
            x, y, z, result = instance.get_position(i)
            self.assertAlmostEquals(2*i, x)
            self.assertAlmostEquals(-2*i, y)
            self.assertAlmostEquals(0, z)
            self.assertEquals(0,result)
            
            self.assertEquals(0, instance.set_velocity(i, -i*20, i*20, 0))
            vx, vy, vz, result = instance.get_velocity(i)
            self.assertAlmostEquals(-i*20, vx)
            self.assertAlmostEquals(i*20, vy)
            self.assertAlmostEquals(0, vz)
            self.assertEquals(0,result)
        
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()


    def test9(self):
        instance = Gadget2Interface(mode = Gadget2Interface.MODE_PERIODIC_BOUNDARIES, 
            **few_particles_default_options)
        instance.initialize_code()
        instance.set_min_size_timestep(1.0)
        instance.set_gadget_output_directory(instance.get_output_directory())
        
        instance.set_box_size(2.)
        value, error = instance.get_box_size()
        self.assertEqual(error, 0)
        self.assertEqual(value, 2.)
        
        self.assertEqual(instance.get_periodic_boundaries_flag()['value'], False) # default, has to be changed for periodic runs
        self.assertEqual(instance.commit_parameters(), -4)
        instance.set_periodic_boundaries_flag(True)
        self.assertEqual(instance.get_periodic_boundaries_flag()['value'], True)
        self.assertEqual(instance.recommit_parameters(), 0)
        ids,err=instance.new_particle( 
           [1.0,1.0,1.0],
           [0.5,0.0,0.0],
           [0.0,-0.5,0.0],
           [0.0,0.0,0.5],
           [-1.0,0.0,0.0],
           [0.0,1.0,0.0],
           [0.0,0.0,-1.0])
        instance.commit_particles()
        m,x,y,z,vx,vy,vz,err=instance.get_state(ids)
        self.assertAlmostEqual(x, [0.5,0.,0.], places=6)
        self.assertAlmostEqual(y, [0.,1.5,0.], places=6)
        self.assertAlmostEqual(z, [0.,0.,0.5], places=6)
        instance.evolve(0.1)
        m,x,y,z,vx,vy,vz,err=instance.get_state(ids)
        self.assertAlmostEqual(x, [0.4,0.,0.], places=6)
        self.assertAlmostEqual(y, [0.,1.6,0.], places=6)
        self.assertAlmostEqual(z, [0.,0.,0.4], places=6)
        
        instance.evolve(1.0)
        m,x,y,z,vx,vy,vz,err=instance.get_state(ids)
        self.assertAlmostEqual(x, [1.5,0.,0.], places=6)
        self.assertAlmostEqual(y, [0.,0.5,0.], places=6)
        self.assertAlmostEqual(z, [0.,0.,1.5], places=6)
        instance.cleanup_code()
        instance.stop()
        
        instance = Gadget2Interface(mode = Gadget2Interface.MODE_NORMAL,
            **few_particles_default_options) # MODE_NORMAL is default: non-periodic
        instance.initialize_code()
        instance.set_gadget_output_directory(instance.get_output_directory())
        
        self.assertEqual(instance.get_periodic_boundaries_flag()['value'], False) # default, has to be changed for periodic runs
        self.assertEqual(instance.commit_parameters(), 0)
        instance.set_periodic_boundaries_flag(True)
        self.assertEqual(instance.get_periodic_boundaries_flag()['value'], True)
        self.assertEqual(instance.recommit_parameters(), -4)
        instance.stop()
        
    
    def test10(self):
        instance = Gadget2Interface(**few_particles_default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        
        indices_and_value = []
        for i in range(10):
            index, error = instance.new_sph_particle(0.1, i,0,0, 0,0,0, 0.0)
            self.assertEquals(0, error)
            indices_and_value.append((index, i))
            
        instance.commit_particles()
        for i, value in indices_and_value:
            x, y, z, error = instance.get_position(i)
            self.assertEquals(0, error)
            self.assertEquals(x, value)
            self.assertEquals(y, 0)
            self.assertEquals(z, 0)
        
        
        self.assertEquals(0, error)
        for i in range(10,20):
            index, error = instance.new_sph_particle(0.1, i,0,0, 0,0,0, 0.0)
            self.assertEquals(0, error)
            indices_and_value.append((index, i))
            
        error = instance.recommit_particles()
        self.assertEquals(0, error)
        
        for i, value in indices_and_value:
            x, y, z, error = instance.get_position(i)
            self.assertEquals(0, error)
            self.assertEquals(x, value)
            self.assertEquals(y, 0)
            self.assertEquals(z, 0)
            
        instance.stop()
    
    
class TestGadget2(TestWithMPI):

    UnitLength = 3.085678e21 | units.cm     # ~ 1.0 kpc
    UnitMass = 1.989e43 | units.g           # 1.0e10 solar masses
    UnitVelocity = 1e5 | units.cm / units.s # 1 km/sec
    default_converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(UnitLength, UnitMass, UnitVelocity)
    default_convert_nbody = nbody_system.nbody_to_si(UnitLength, UnitMass)
    
    three_particles_IC = core.Particles(3)
    three_particles_IC.position = [[0.5, 0.0, 0.0], [0.0,-0.5, 0.0], [0.0, 0.0, 0.5]] | units.kpc 
    three_particles_IC.velocity =[[-1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0,-1.0]] | units.km / units.s
    three_particles_IC.mass = 1.0e10 | units.MSun
    
    def test1(self):
        print "Testing Gadget initialization"
        instance = Gadget2(self.default_converter, **default_options)
        instance.initialize_code()
        self.assertTrue("/data/gadget2/output/" in str(instance.parameters.gadget_output_directory))
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()

    def test2(self):
        print "Testing Gadget parameters"
        instance = Gadget2(self.default_converter, **default_options)
        instance.initialize_code()
        self.assertAlmostEquals(instance.parameters.epsilon_squared, (0.01 | units.kpc)**2)
        self.assertAlmostRelativeEquals(instance.parameters.code_mass_unit, self.UnitMass, 7)
        self.assertAlmostRelativeEquals(instance.parameters.code_time_unit, 3.085678e16 | units.s, 7)
        self.assertAlmostRelativeEquals(instance.parameters.code_length_unit, self.UnitLength, 7)
        instance.parameters.epsilon_squared = 0.01 | units.kpc**2
        instance.parameters.opening_angle = 0.4 | units.none
        instance.commit_parameters()
        self.assertAlmostEquals(instance.parameters.epsilon_squared, 0.01 | units.kpc**2)
        self.assertAlmostEquals(instance.parameters.opening_angle, 0.4 | units.none)
        self.assertAlmostRelativeEquals(instance.parameters.code_mass_unit, 1.989e43 | units.g, 7)
        self.assertAlmostRelativeEquals(instance.parameters.code_time_unit, 3.085678e16 | units.s, 7)
        instance.stop()

    def test3(self):
        print "Testing Gadget, 2 nbody particles"
         # No default_options since this test fails for (number_of_workers > 1):
        instance = Gadget2(self.default_converter, **few_particles_default_options)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.01 | units.kpc**2

        dark = core.Particles(2)
        dark.mass = [0.4, 0.4] | generic_unit_system.mass
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.kpc
        dark.velocity = [[100.0,100.0,100.0], [1.0,1.0,1.0]] | units.km / units.s
        dark.radius = [0.0, 0.0] | units.RSun

        instance.dm_particles.add_particles(dark)
        instance.evolve_model(1.0 | units.Myr)
        self.assertAlmostEquals(instance.model_time, 1.0 | units.Myr, 3)
        instance.stop()

    def test4(self):
        print "Testing Gadget, 100 nbody particles"
        instance = Gadget2(self.default_converter, **default_options)
        instance.initialize_code()
        convert_nbody = nbody_system.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        dark = MakePlummerModel(100, convert_nbody).result
        instance.dm_particles.add_particles(dark)
        instance.evolve_model(1.0 | units.Myr)
        instance.stop()

    def test5(self):
        print "Test 5: testing SPH particles"
        target_number_sph_particles = 100
        gas = new_evrard_gas_sphere(target_number_sph_particles, self.default_convert_nbody, seed = 1234)
        instance = Gadget2(self.default_converter, **default_options)
        instance.initialize_code()
        instance.gas_particles.add_particles(gas)
        instance.evolve_model(0.0001 | generic_unit_system.time)
        instance.stop()

    def test6(self):
        print "Test 6: testing dark matter + SPH particles"
        target_number_sph_particles = 100
        gas = new_evrard_gas_sphere(target_number_sph_particles, self.default_convert_nbody, seed = 1234)

        dark = core.Particles(2)
        dark.mass = [0.4, 0.4] | generic_unit_system.mass
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.kpc
        dark.velocity = [[100.0,100.0,100.0], [1.0,1.0,1.0]] | units.km / units.s

        instance = Gadget2(self.default_converter, **default_options)
        instance.dm_particles.add_particles(dark)
        instance.gas_particles.add_particles(gas)
        instance.evolve_model(1.0 | units.Myr)
        instance.stop()

    def test7(self):
        print "Testing more Gadget parameters"
        instance = Gadget2(self.default_converter, **default_options)
        instance.initialize_code()
        
        for par, value in [('gadget_cell_opening_flag', True), 
                ('comoving_integration_flag', False), 
                ('periodic_boundaries_flag', False)]:
            self.assertTrue(value is eval("instance.parameters."+par))
            exec("instance.parameters."+par+" = not value")
            self.assertFalse(value is eval("instance.parameters."+par))
        
        for par, value in [('time_limit_cpu', 36000 | units.s), 
                ('hubble_param', 0.7 | 100 * units.km / units.s / units.Mpc),
                ('min_gas_temp', 0.0 | units.K)]:
            self.assertEquals(value, eval("instance.parameters."+par))
            exec("instance.parameters."+par+" = 2 * value")
            self.assertEquals(2 * value, eval("instance.parameters."+par))
        
        for par, value in [('n_smooth_tol',0.1), ('n_smooth',50), ('opening_angle',0.5),
                ('gadget_cell_opening_constant',0.005), ('artificial_viscosity_alpha',0.5),
                ('courant',0.3), ('type_of_timestep_criterion',0), ('omega_zero',0.0),
                ('omega_lambda',0.0), ('omega_baryon',0.0), ('min_gas_hsmooth_fractional',0.0),
                ('timestep_accuracy_parameter',0.025), ('tree_domain_update_frequency',0.05)]:
            self.assertEquals(value | units.none, eval("instance.parameters."+par))
            exec("instance.parameters."+par+" = 1 | units.none")
            self.assertEquals(1 | units.none, eval("instance.parameters."+par))
        
        for par, value in [('gas_epsilon', 0.01 | generic_unit_system.length), 
                ('time_begin', 0.0 | generic_unit_system.time), 
                ('time_max', 100.0 | generic_unit_system.time), 
                ('max_size_timestep', 0.01 | generic_unit_system.time), 
                ('min_size_timestep', 0.0 | generic_unit_system.time), 
                ('periodic_box_size', 1.0 | generic_unit_system.length),
                ('time_between_statistics', 0.1 | generic_unit_system.time),
                ('softening_gas_max_phys', 0.0 | generic_unit_system.length),
                ('softening_halo_max_phys', 0.0 | generic_unit_system.length)]:
            self.assertEquals(instance.unit_converter.to_si(value), 
                eval("instance.parameters."+par))
            exec("instance.parameters."+par+" = 3.0 | value.unit")
            self.assertEquals(instance.unit_converter.to_si(3.0 | value.unit),
                eval("instance.parameters."+par))
        
        for par, value in [('energy_file',"energy.txt"),('info_file',"info.txt"),
                ('timings_file',"timings.txt"),('cpu_file',"cpu.txt")]:
            self.assertEquals(value | units.string, eval("instance.parameters."+par))
            exec("instance.parameters."+par+" = 'test.txt' | units.string")
            self.assertEquals("test.txt" | units.string, eval("instance.parameters."+par))
        
        
        instance.stop()

    def test8(self):
        print "Testing read-only Gadget parameters"
        instance = Gadget2(self.default_converter, **default_options)
        instance.initialize_code()
        for par, value in [ ('no_gravity_flag',False),
                            ('isothermal_flag', False),
                            ('eps_is_h_flag',   False),
                            ('code_mass_unit',     self.default_converter.to_si(generic_unit_system.mass)),
                            ('code_length_unit',   self.default_converter.to_si(generic_unit_system.length)),
                            ('code_time_unit',     self.default_converter.to_si(generic_unit_system.time)),
                            ('code_velocity_unit', self.default_converter.to_si(generic_unit_system.speed)),
                            ('polytropic_index_gamma', (5.0/3) | units.none)]:
            self.assertEquals(value, eval("instance.parameters."+par))
            def try_set_parameter(par, value, instance):
                exec("instance.parameters."+par+" = value")
            self.assertRaises(AmuseException, try_set_parameter, par, value, instance,
                expected_message = "Could not set value for parameter '"+par+"' of a 'Gadget2' object, "
                    "parameter is read-only")
        instance.stop()

    def test9(self):
        print "Testing Gadget properties"
        target_number_of_particles = 100
        gas = new_evrard_gas_sphere(target_number_of_particles, self.default_convert_nbody, do_scale=True, seed = 1234)
        instance = Gadget2(self.default_converter, **default_options)
        instance.gas_particles.add_particles(gas)
        self.assertEquals(instance.model_time,                        0.0 | units.s)
        self.assertAlmostEquals(instance.potential_energy, -4.27843220393 | 1e+50*units.J)
        self.assertAlmostEquals(instance.kinetic_energy,              0.0 | 1e+49*units.J)
        self.assertAlmostEquals(instance.thermal_energy,    4.27851824913 | 1e+49*units.J)
        self.assertAlmostEquals(instance.total_radius,      3.96592921066 | 1e+19*units.m)
        self.assertAlmostEquals(instance.center_of_mass_position, [0,0,0] | 1e+19*units.m)
        self.assertAlmostEquals(instance.center_of_mass_velocity, [0,0,0] | units.m/units.s)
        self.assertAlmostEquals(instance.total_mass,                1.989 | 1e+40*units.kg)
        instance.stop()

    def test10(self):
        print "Testing Gadget states"
        target_number_sph_particles = 100
        gas = new_evrard_gas_sphere(target_number_sph_particles, self.default_convert_nbody, seed = 1234)
        dark = core.Particles(2)
        dark.mass = [0.4, 0.4] | generic_unit_system.mass
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.kpc
        dark.velocity = [[100.0,100.0,100.0], [1.0,1.0,1.0]] | units.km / units.s

        print "First do everything manually:"
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        instance.gas_particles.add_particles(gas)
        instance.commit_particles()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        mass = instance.gas_particles[0].mass
        instance.evolve_model(0.001 | generic_unit_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()

        print "commit_parameters(), (re)commit_particles(), and cleanup_code() should be called " \
            "automatically before new_xx_particle(), get_xx(), and stop():"
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.gas_particles.add_particles(gas)
        self.assertEquals(instance.get_name_of_current_state(), 'EDIT')
        mass = instance.gas_particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.dm_particles.add_particles(dark)
        self.assertEquals(instance.get_name_of_current_state(), 'UPDATE')
        mass = instance.gas_particles[0].mass
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.evolve_model(0.001 | generic_unit_system.time)
        self.assertEquals(instance.get_name_of_current_state(), 'EVOLVED')
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'END')

    def test11(self):
        particles = core.Particles(2)

        particles.x = [0.0,10.0] | generic_unit_system.length
        particles.y = 0 | generic_unit_system.length
        particles.z = 0 | generic_unit_system.length
        particles.radius = 0.005 | generic_unit_system.length
        particles.vx =  0 | generic_unit_system.speed
        particles.vy =  0 | generic_unit_system.speed
        particles.vz =  0 | generic_unit_system.speed
        particles.mass = 1.0 | generic_unit_system.mass

        instance = Gadget2(self.default_converter, **few_particles_default_options)
        instance.initialize_code()
        instance.parameters.stopping_conditions_number_of_steps = 2
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps, 2 | units.none)
        instance.parameters.epsilon_squared = (0.01 | generic_unit_system.length)**2
        instance.particles.add_particles(particles)
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(1.0e15 | units.s)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        self.assertTrue(instance.model_time < 1.0e15 | units.s)
        
        instance.stop()
    
    def test12(self):
        print "Testing Gadget get_hydro_state_at_point"
        number_sph_particles = 100
        gas = new_evrard_gas_sphere(number_sph_particles, self.default_convert_nbody, seed = 1234)
        instance = Gadget2(self.default_converter, **default_options)
        instance.parameters.n_smooth     =   64 | units.none
        instance.parameters.n_smooth_tol = 0.01 | units.none
        instance.gas_particles.add_particles(gas)
        
        coords = [0.0 | units.kpc]*3
        speeds = [0.0 | units.m / units.s]*3
        hydro_state = instance.get_hydro_state_at_point(*(coords+speeds))
        expected = [ 3.5469e-19 | units.kg * units.m**-3, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                     8.4252e-10 | units.kg * units.m**-1 * units.s**-2]
        for value, expect in zip(hydro_state, expected):
            self.assertAlmostRelativeEqual(value, expect, places=3)
        
        coords = [0.1 | units.kpc]*3
        speeds = [0.0 | units.m / units.s]*3
        hydro_state = instance.get_hydro_state_at_point(*(coords+speeds))
        expected = [ 4.1456e-19 | units.kg * units.m**-3, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                            0.0 | units.kg * units.m**-2 / units.s, 
                     1.0742e-09 | units.kg * units.m**-1 * units.s**-2]
        for value, expect in zip(hydro_state, expected):
            self.assertAlmostRelativeEqual(value, expect, places=3)
        
        instance.stop()
    
    def test13(self):
        print "Testing Gadget get_hydro_state_at_point II: uniform sphere"
        number_sph_particles = 1000
        gas = new_uniform_spherical_particle_distribution(number_sph_particles, self.UnitLength, self.UnitMass)
        gas.velocity = [10.0, 20.0, 30.0] | units.km / units.s
        gas.u = 0.05 | generic_unit_system.specific_energy
        density = (1.0e10 | units.MSun) / (4.0/3.0 * numpy.pi * (1.0 | units.kpc)**3)
        
        instance = Gadget2(self.default_converter, **default_options)
        instance.parameters.n_smooth     =   64 | units.none
        instance.parameters.n_smooth_tol = 0.01 | units.none
        instance.gas_particles.add_particles(gas)
        
        coords = [0.0 | units.kpc]*3
        speeds = [0.0 | units.m / units.s]*3
        rho, rhovx, rhovy, rhovz, rhoe = instance.get_hydro_state_at_point(*(coords + speeds))
        self.assertAlmostRelativeEqual(rho,   density,                              places=3)
        self.assertAlmostRelativeEqual(rhovx, density*instance.gas_particles[0].vx, places=3)
        self.assertAlmostRelativeEqual(rhovy, density*instance.gas_particles[0].vy, places=3)
        self.assertAlmostRelativeEqual(rhovz, density*instance.gas_particles[0].vz, places=3)
        self.assertAlmostRelativeEqual(rhoe,  density * (instance.gas_particles[0].u + 
            instance.gas_particles[0].velocity.length_squared()),  places=3)
        
        coords = [0.1 | units.kpc]*3
        rho, rhovx, rhovy, rhovz, rhoe = instance.get_hydro_state_at_point(*(coords + speeds))
        self.assertAlmostRelativeEqual(rho,   density,                              places=3)
        self.assertAlmostRelativeEqual(rhovx, density*instance.gas_particles[0].vx, places=3)
        self.assertAlmostRelativeEqual(rhovy, density*instance.gas_particles[0].vy, places=3)
        self.assertAlmostRelativeEqual(rhovz, density*instance.gas_particles[0].vz, places=3)
        self.assertAlmostRelativeEqual(rhoe,  density * (instance.gas_particles[0].u + 
            instance.gas_particles[0].velocity.length_squared()),  places=3)
        instance.stop()
    
    def test14(self):
        print "Testing Gadget SPH particle properties"
        number_sph_particles = 1000
        gas = new_evrard_gas_sphere(number_sph_particles, self.default_convert_nbody, seed = 1234)
        instance = Gadget2(self.default_converter, **default_options)
        instance.gas_particles.add_particles(gas)
        self.assertEqual(instance.gas_particles.num_neighbours >= 
            instance.parameters.n_smooth*(1-instance.parameters.n_smooth_tol), True)
        self.assertEqual(instance.gas_particles.num_neighbours <= 
            instance.parameters.n_smooth*(1+instance.parameters.n_smooth_tol), True)
        self.assertIsOfOrder(instance.gas_particles.h_smooth, 
            self.default_convert_nbody.to_si(1.0 | nbody_system.length) * 
            (instance.parameters.n_smooth*1.0/number_sph_particles)**(1.0/3))
        self.assertAlmostRelativeEqual(instance.gas_particles.u, 
            self.default_convert_nbody.to_si(0.05 | nbody_system.specific_energy))
        
        # the density of the cloud scales with 1/r:
        r_sort, rho_sort = instance.gas_particles.position.lengths().sorted_with(instance.gas_particles.rho)
        mean_density = self.default_convert_nbody.to_si(3.0/(4.0*numpy.pi) | nbody_system.density)
        select = slice(number_sph_particles/2) # select 50% particles closest to center to avoid boundaries
        self.assertIsOfOrder(rho_sort[select]/mean_density, r_sort.mean()/r_sort[select])
    
    def test15(self):
        instance = Gadget2(mode = Gadget2Interface.MODE_PERIODIC_BOUNDARIES, 
            **few_particles_default_options)
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        self.assertEqual(instance.parameters.periodic_boundaries_flag, False)
        # 'False' is default, will be changed for periodic runs in initialize_code():
        instance.parameters.periodic_box_size = 2.0 | generic_unit_system.length # implicitly calls initialize_code()...
        self.assertEqual(instance.parameters.periodic_boundaries_flag, True)
        self.assertAlmostEqual(instance.parameters.periodic_box_size, 2.0 | units.kpc, places=6)
        instance.parameters.min_size_timestep = 1.0 | generic_unit_system.time
        
        instance.dm_particles.add_particles(self.three_particles_IC)
        self.assertAlmostEqual(instance.dm_particles.x, [0.5,0.,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.y, [0.,1.5,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.z, [0.,0.,0.5] | units.kpc, places=6)
        
        instance.evolve_model(0.1 | generic_unit_system.time)
        self.assertAlmostEqual(instance.dm_particles.x, [0.4,0.,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.y, [0.,1.6,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.z, [0.,0.,0.4] | units.kpc, places=6)
        
        instance.evolve_model(1.0 | generic_unit_system.time)
        self.assertAlmostEqual(instance.dm_particles.x, [1.5,0.,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.y, [0.,0.5,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.z, [0.,0.,1.5] | units.kpc, places=6)
        instance.stop()
        
    def test16(self):
        instance = Gadget2(mode = Gadget2Interface.MODE_PERIODIC_BOUNDARIES,
            **few_particles_default_options)
        self.assertEqual(instance.parameters.periodic_boundaries_flag, False) # 'False' is default
        instance.initialize_code() # initialize_code() will set periodic_boundaries_flag correctly:
        self.assertEqual(instance.parameters.periodic_boundaries_flag, True)
        instance.parameters.periodic_boundaries_flag = False # but we are stubborn...
        self.assertEqual(instance.parameters.periodic_boundaries_flag, False)
        self.assertRaises(AmuseException, instance.commit_parameters, expected_message = "Error when "
            "calling 'commit_parameters' of a 'Gadget2', errorcode is -4, error is 'Parameter check failed.'")
        instance.parameters.periodic_boundaries_flag = True
        instance.commit_parameters() # Now it's ok again.
        instance.dm_particles.add_particles(self.three_particles_IC)
        instance.parameters.periodic_boundaries_flag = False # suppose we are REALLY stubborn...
        self.assertRaises(AmuseException, instance.evolve_model, 1.0 | generic_unit_system.time, expected_message = "Error when "
            "calling 'recommit_parameters' of a 'Gadget2', errorcode is -4, error is 'Parameter check failed.'")
        instance.stop()
    
    def test17(self):
        print "Testing Gadget parameters further"
        target_number_of_particles = 100
        gas = new_evrard_gas_sphere(target_number_of_particles, self.default_convert_nbody, seed = 1234)
        instance = Gadget2(self.default_converter, **default_options)
        instance.parameters.time_limit_cpu = 0.001 | units.s
        instance.gas_particles.add_particles(gas)
        self.assertRaises(AmuseException, instance.evolve_model, 1.0 | generic_unit_system.time, expected_message = 
            "Error when calling 'evolve_model' of a 'Gadget2', errorcode is -5, error is 'CPU-time limit reached.'")
        instance.parameters.time_limit_cpu = 10.0 | units.s
        instance.evolve_model(0.0001 | generic_unit_system.time)
        self.assertRaises(AmuseException, instance.evolve_model, 0.00005 | generic_unit_system.time, expected_message = 
            "Error when calling 'evolve_model' of a 'Gadget2', errorcode is -6, error is 'Can't evolve backwards in time.'")
        instance.stop()
    
    def test18(self):
        print "Testing Gadget parameters further"
        instance = Gadget2(mode = Gadget2Interface.MODE_PERIODIC_BOUNDARIES,
            **few_particles_default_options)
        instance.parameters.periodic_box_size = 2.0 | generic_unit_system.length # implicitly calls initialize_code()...
        instance.parameters.time_begin = 10.04 | generic_unit_system.time
        instance.parameters.time_max = 10.4 | generic_unit_system.time
        instance.parameters.min_size_timestep = 0.1 | generic_unit_system.time
        
        instance.dm_particles.add_particles(self.three_particles_IC)
        self.assertAlmostEqual(instance.dm_particles.x, [0.5,0.,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.y, [0.,1.5,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.z, [0.,0.,0.5] | units.kpc, places=6)
        
        instance.evolve_model(10.14 | generic_unit_system.time)
        self.assertAlmostEqual(instance.dm_particles.x, [0.4,0.,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.y, [0.,1.6,0.] | units.kpc, places=6)
        self.assertAlmostEqual(instance.dm_particles.z, [0.,0.,0.4] | units.kpc, places=6)
        instance.stop()
    
    def test19(self):
        particles = new_plummer_sphere(31)
       
        instance = Gadget2(self.default_converter, number_of_workers=1)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.01 | generic_unit_system.length ** 2
        instance.particles.add_particles(particles)
        instance.commit_particles()
        
        instance.evolve_model(0.001 | generic_unit_system.time)
        expected_positions = instance.particles.position
        instance.stop()
        positions_per_workers = []
        for n in [2,3,4,5]:
            instance = Gadget2(self.default_converter, number_of_workers=n)
            instance.initialize_code()
            instance.parameters.epsilon_squared = 0.01 | generic_unit_system.length ** 2
            instance.commit_parameters()
            instance.particles.add_particles(particles)
            instance.commit_particles()
            
            instance.evolve_model(0.001 | generic_unit_system.time)
            positions_per_workers.append(instance.particles.position)
            instance.stop()
         
         
        for index, n in enumerate([2,3,4,5]):
            self.assertAlmostRelativeEqual(expected_positions, positions_per_workers[index])
