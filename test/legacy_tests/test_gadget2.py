import numpy
from amuse.test.amusetest import TestWithMPI
from amuse.legacy.gadget2.interface import Gadget2Interface, Gadget2
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.evrard_test import MakeEvrardTest, new_evrard_gas_sphere
from amuse.ext.spherical_model import new_uniform_spherical_particle_distribution

from amuse.support.exceptions import AmuseException
from amuse.support.units import nbody_system
from amuse.support.units import generic_unit_converter
from amuse.support.units import generic_unit_system
from amuse.support.units import units
from amuse.support.data import core
from amuse.support.legacy import channel

default_options = dict()#redirection = "null")
#default_options = dict(debugger = "none")

class TestGadget2Interface(TestWithMPI):

    def test1(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.set_parameterfile_path(instance.default_path_to_parameterfile))
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test2(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.set_parameterfile_path(instance.default_path_to_parameterfile))
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals([1.989000e43, 0], instance.get_unit_mass().values())
        self.assertEquals([3.085678e21, 0], instance.get_unit_length().values())
        self.assertEquals([3.085678e16, 0], instance.get_unit_time().values())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test3(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.set_parameterfile_path(instance.default_path_to_parameterfile))
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
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.set_parameterfile_path(instance.default_path_to_parameterfile))
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_gadget_output_directory(instance.get_output_directory()))
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals([1, 0], instance.new_dm_particle(0.01,  1, 0, 0,  0, 1, 0).values())
        self.assertEquals([2, 0], instance.new_dm_particle(0.02, -1, 0, 0,  0,-1, 0).values())
        self.assertEquals(-1, instance.get_mass(1)['__result'])
        self.assertEquals(0, instance.commit_particles())

        mass, result = instance.get_mass(1)
        self.assertAlmostEquals(0.01,mass)
        self.assertEquals(0,result)
        mass, result = instance.get_mass(2)
        self.assertAlmostEquals(0.02,mass)
        self.assertEquals(0,result)
        self.assertEquals(-1, instance.get_mass(3)['__result'])
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
        self.assertEquals(0, instance.set_parameterfile_path(instance.default_path_to_parameterfile))
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
        self.assertEquals(0, instance.evolve(0.001))
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test6(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.set_parameterfile_path(instance.default_path_to_parameterfile))
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

        self.assertEquals(0, instance.evolve(0.001))
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()

    def test7(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.set_parameterfile_path(instance.default_path_to_parameterfile))
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
        self.assertEquals(0, instance.evolve(0.0001))
        self.assertEquals(0, instance.evolve(0.0002))
        self.assertEquals(0, instance.delete_particle(number_of_particles-1))
        self.assertEquals(0, instance.delete_particle(number_of_particles+1))
        self.assertEquals(-1, instance.delete_particle(number_of_particles-1))
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

        self.assertEquals(0, instance.evolve(0.0003))
        self.assertEquals(0, instance.evolve(0.0004))
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test8(self):
        instance = Gadget2Interface(**default_options)
        self.assertEquals(0, instance.set_parameterfile_path(instance.default_path_to_parameterfile))
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

class TestGadget2(TestWithMPI):

    UnitLength = 3.085678e21 | units.cm     # ~ 1.0 kpc
    UnitMass = 1.989e43 | units.g           # 1.0e10 solar masses
    UnitVelocity = 1e5 | units.cm / units.s # 1 km/sec
    default_converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(UnitLength, UnitMass, UnitVelocity)
    default_convert_nbody = nbody_system.nbody_to_si(UnitLength, UnitMass)

    def test1(self):
        print "Testing Gadget initialization"
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertTrue("/data/gadget2/output/" in str(instance.parameters.gadget_output_directory))
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()

    def test2(self):
        print "Testing Gadget parameters"
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(0, instance.initialize_code())
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
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(0, instance.initialize_code())
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
        self.assertEquals(0, instance.initialize_code())
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
        self.assertEquals(0, instance.initialize_code())
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
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(True, instance.parameters.gadget_cell_opening_flag)
        instance.parameters.gadget_cell_opening_flag = False
        self.assertEquals(False, instance.parameters.gadget_cell_opening_flag)

        for par, value in [('n_smooth_tol',0.1),('n_smooth',50),('opening_angle',0.5),('gadget_cell_opening_constant',0.005),
            ('artificial_viscosity_alpha',0.5),('courant',0.3)]:
            self.assertEquals(value | units.none, eval("instance.parameters."+par))
            exec("instance.parameters."+par+" = 1 | units.none")
            self.assertEquals(1 | units.none, eval("instance.parameters."+par))

        self.assertEquals(instance.unit_converter.to_si(0.01 | generic_unit_system.length),
            instance.parameters.gas_epsilon)
        instance.parameters.gas_epsilon = 0.1 | generic_unit_system.length
        self.assertEquals(instance.unit_converter.to_si(0.1 | generic_unit_system.length),
            instance.parameters.gas_epsilon)
        instance.stop()

    def test8(self):
        print "Testing read-only Gadget parameters"
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(0, instance.initialize_code())
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
        self.assertEquals(0, instance.initialize_code())
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
        self.assertEquals(0, instance.initialize_code())
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

        instance = Gadget2(self.default_converter, **default_options)
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
        gas = new_uniform_spherical_particle_distribution(number_sph_particles, self.UnitLength, self.UnitMass, seed = 1234)
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
    
