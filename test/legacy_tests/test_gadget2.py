from amuse.test.amusetest import TestWithMPI
from amuse.legacy.gadget2.interface import Gadget2Interface, Gadget2
from amuse.ext.plummer import MakePlummerModel
from amuse.ext.evrard_test import MakeEvrardTest

from amuse.support.units import nbody_system as nbody
from amuse.support.units import generic_unit_converter as generic_system
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

    UnitLength_in_cm = 3.085678e21 | units.cm # 1.0 kpc
    UnitMass_in_g = 1.989e43 | units.g        # 1.0e10 solar masses
    UnitVelocity_in_cm_per_s = 1e5 | units.cm / units.s # 1 km/sec
    default_converter = generic_system.generic_to_si(UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s)
    
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
        self.assertAlmostRelativeEquals(instance.parameters.code_mass_unit, self.UnitMass_in_g, 7)
        self.assertAlmostRelativeEquals(instance.parameters.code_time_unit, 3.085678e16 | units.s, 7)
        self.assertAlmostRelativeEquals(instance.parameters.code_length_unit, self.UnitLength_in_cm, 7)
        instance.parameters.epsilon_squared = 0.01 | units.kpc**2
        instance.parameters.timestep = 0.1 | units.Myr
        instance.commit_parameters()
        self.assertAlmostEquals(instance.parameters.epsilon_squared, 0.01 | units.kpc**2)
        self.assertAlmostEquals(instance.parameters.timestep, 0.1 | units.Myr)
        self.assertAlmostRelativeEquals(instance.parameters.code_mass_unit, 1.989e43 | units.g, 7)
        self.assertAlmostRelativeEquals(instance.parameters.code_time_unit, 3.085678e16 | units.s, 7)
        instance.cleanup_code()
        instance.stop()
    
    def test3(self):
        print "Testing Gadget, 2 nbody particles"
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(0, instance.initialize_code())
        instance.parameters.epsilon_squared = 0.01 | units.kpc**2
        instance.commit_parameters()
        
        dark = core.Particles(2)
        dark.mass = [0.4, 0.4] | generic_system.mass
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.kpc
        dark.velocity = [[100.0,100.0,100.0], [1.0,1.0,1.0]] | units.km / units.s
        dark.radius = [0.0, 0.0] | units.RSun
        
        instance.dm_particles.add_particles(dark)
        
        instance.commit_particles()
        instance.evolve_model(1.0 | units.Myr)
        print instance.model_time
        self.assertAlmostEquals(instance.model_time, 1.0 | units.Myr, 3)
        instance.cleanup_code()
        instance.stop()
    
    def test4(self):
        print "Testing Gadget, 1000 nbody particles"
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(0, instance.initialize_code())
        instance.commit_parameters()
        convert_nbody = nbody.nbody_to_si(1.0e9 | units.MSun, 1.0 | units.kpc)
        dark = MakePlummerModel(1000, convert_nbody).result
        dark.radius = 0.0 | units.RSun
        instance.dm_particles.add_particles(dark)
        instance.commit_particles()
        instance.evolve_model(1.0 | units.Myr)
        instance.cleanup_code()
        instance.stop()
    
    def test5(self):
        print "Test 5: testing SPH particles"
        target_number_sph_particles = 1000
        mass,x,y,z,vx,vy,vz,u=MakeEvrardTest(target_number_sph_particles).new_model()
        number_sph_particles = len(mass)
        attributes = ['mass','x','y','z','vx','vy','vz','u']
        attr_units = [generic_system.mass, generic_system.length, generic_system.length, generic_system.length, 
            generic_system.speed, generic_system.speed, generic_system.speed, generic_system.energy]
        
        gas = core.Particles(number_sph_particles)
        for attribute, unit in zip(attributes,attr_units):
            setattr(gas, attribute, unit.new_quantity(eval(attribute)))
        
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(0, instance.initialize_code())
        instance.parameters.timestep = 0.005 | generic_system.time
        instance.commit_parameters()
        instance.gas_particles.add_particles(gas)
        
        instance.commit_particles()
        instance.evolve_model(0.0001 | generic_system.time)
        instance.cleanup_code()
        instance.stop()
    
    def test6(self):
        print "Test 6: testing dark matter + SPH particles"
        target_number_sph_particles = 100
        mass,x,y,z,vx,vy,vz,u=MakeEvrardTest(target_number_sph_particles).new_model()
        number_sph_particles = len(mass)
        attributes = ['mass','x','y','z','vx','vy','vz','u']
        attr_units = [generic_system.mass, generic_system.length, generic_system.length, generic_system.length, 
            generic_system.speed, generic_system.speed, generic_system.speed, generic_system.energy]
        
        gas = core.Particles(number_sph_particles)
        for attribute, unit in zip(attributes,attr_units):
            setattr(gas, attribute, unit.new_quantity(eval(attribute)))
        
        dark = core.Particles(2)
        dark.mass = [0.4, 0.4] | generic_system.mass
        dark.position = [[0.0,0.0,0.0], [1.0,0.0,0.0]] | units.kpc
        dark.velocity = [[100.0,100.0,100.0], [1.0,1.0,1.0]] | units.km / units.s
        
        instance = Gadget2(self.default_converter, **default_options)
        self.assertEquals(0, instance.initialize_code())
        instance.commit_parameters()
        instance.dm_particles.add_particles(dark)
        instance.gas_particles.add_particles(gas)
        
        instance.commit_particles()
        instance.evolve_model(1.0 | units.Myr)
        instance.cleanup_code()
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
        
        self.assertEquals(instance.unit_converter.to_si(0.01 | generic_system.length), 
            instance.parameters.gas_epsilon)
        instance.parameters.gas_epsilon = 0.1 | generic_system.length
        self.assertEquals(instance.unit_converter.to_si(0.1 | generic_system.length), 
            instance.parameters.gas_epsilon)
        instance.cleanup_code()
        instance.stop()
    
