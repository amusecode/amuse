from amuse.test.amusetest import TestWithMPI
from amuse.legacy.gadget2.interface import Gadget2Interface, Gadget2
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.evrard_test import MakeEvrardTest, new_evrard_gas_sphere

from amuse.support.units import nbody_system
from amuse.support.units import generic_unit_converter
from amuse.support.units import generic_unit_system 
from amuse.support.units import units
from amuse.support.data import core
from amuse.legacy.support import channel

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
    

class TestGadget2(TestWithMPI):

    UnitLength = 3.085678e21 | units.cm     # 1.0 kpc
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
        
        for par, value in [('n_neighbour_tol',0.1),('nsmooth',50),('opening_angle',0.5),('gadget_cell_opening_constant',0.005),
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
            try:
                exec("instance.parameters."+par+" = value")
                self.fail("Should not be able to set read-only parameters.")
            except Exception as ex:
                self.assertEquals("Could not set value for parameter '"+par+"' of a 'Gadget2' object, "
                    "parameter is read-only", str(ex))
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
    
