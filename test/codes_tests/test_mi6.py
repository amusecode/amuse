import numpy
import math

from amuse.test.amusetest import TestWithMPI
from amuse.units import nbody_system, units, constants
from amuse.datamodel import Particles
from amuse.community.mi6.interface import MI6Interface, MI6

from amuse.ic.plummer import new_plummer_model


default_options = dict()

class TestMI6Interface(TestWithMPI):
    
    def test1(self):
        print "Test MI6Interface initialization"
        instance = MI6Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test2(self):
        print "Test MI6Interface new_particle / get_state"
        instance = MI6Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        
        id, error = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        self.assertEquals(0, error)
        self.assertEquals(1, id)
        id, error = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        self.assertEquals(0, error)
        self.assertEquals(2, id)
        self.assertEquals(0, instance.commit_particles())
        
        retrieved_state1 = instance.get_state(1)
        retrieved_state2 = instance.get_state(2)
        self.assertEquals(0,  retrieved_state1['__result'])
        self.assertEquals(0,  retrieved_state2['__result'])
        self.assertEquals(11.0,  retrieved_state1['mass'])
        self.assertEquals(21.0,  retrieved_state2['mass'])
        self.assertEquals( 0.0,  retrieved_state1['x'])
        self.assertEquals(10.0,  retrieved_state2['x'])
        
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def xtest3(self):
        print "Test MI6Interface get_index_of_first_particle, get_index_of_next_particle"
        instance = MI6Interface(**default_options)
        instance.initialize_code()
        instance.commit_parameters()
        for i in [1, 2, 3]:
            result = instance.new_particle(mass = i, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
            self.assertEquals(i, result['index_of_the_particle'])
        
        instance.commit_particles()
        self.assertEquals(1, instance.get_index_of_first_particle()['index_of_the_particle'])
        self.assertEquals(2, instance.get_index_of_next_particle(1)['index_of_the_next_particle']) 
        self.assertEquals(3, instance.get_index_of_next_particle(2)['index_of_the_next_particle'])
            
        instance.delete_particle(1)
      
        self.assertEquals(2, instance.get_number_of_particles()['number_of_particles'])
        
        #the deletion does a swap, so 3 is copied to 1, (overwriting old 1 and treesize -> treesize-1
        self.assertEquals(3, instance.get_index_of_first_particle()['index_of_the_particle'])
        
        self.assertEquals(1, instance.get_index_of_next_particle(2)['__result'])

        instance.cleanup_code()
        instance.stop()
    
    def test4(self):
        print "Test MI6Interface particle property getters/setters"
        instance = MI6Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals([1, 0], instance.new_particle(0.01,  1, 0, 0,  0, 1, 0, 0.1).values())
        self.assertEquals([2, 0], instance.new_particle(0.02, -1, 0, 0,  0,-1, 0, 0.1).values())
        self.assertEquals(-3, instance.get_mass(1)['__result']) # Have to commit first
        self.assertEquals(0, instance.commit_particles())
        
        # getters
        mass, result = instance.get_mass(1)
        self.assertAlmostEquals(0.01, mass)
        self.assertEquals(0,result)
        radius, result = instance.get_radius(2)
        self.assertAlmostEquals(0.1, radius)
        self.assertEquals(0,result)
        self.assertEquals(-3, instance.get_mass(3)['__result']) # Particle not found
        self.assertEquals([ 1, 0, 0,  0], instance.get_position(1).values())
        self.assertEquals([-1, 0, 0,  0], instance.get_position(2).values())
        self.assertEquals([ 0, 1, 0,  0], instance.get_velocity(1).values())
        self.assertEquals([ 0,-1, 0,  0], instance.get_velocity(2).values())
        
        # setters
        self.assertEquals(0, instance.set_state(1, 0.01, 1,2,3, 4,5,6, 0.1))
        self.assertEquals([0.01, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.1, 0], instance.get_state(1).values())
        self.assertEquals(0, instance.set_mass(1, 0.02))
        self.assertEquals([0.02, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.1, 0], instance.get_state(1).values())
        self.assertEquals(0, instance.set_radius(1, 0.2))
        self.assertEquals([0.02, 1.0,2.0,3.0, 4.0,5.0,6.0, 0.2, 0], instance.get_state(1).values())
        self.assertEquals(0, instance.set_position(1, 10,20,30))
        self.assertEquals([0.02, 10.0,20.0,30.0, 4.0,5.0,6.0, 0.2, 0], instance.get_state(1).values())
        self.assertEquals(0, instance.set_velocity(1, 40,50,60))
        self.assertEquals([0.02, 10.0,20.0,30.0, 40.0,50.0,60.0, 0.2, 0], instance.get_state(1).values())

        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test5(self):
        print "Test MI6Interface parameters"
        instance = MI6Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        
        # MI6 has separate epsilon_squared parameters for different interactions!
        self.assertEquals([0, 0], instance.get_eps2_fs_fs().values())
        self.assertEquals([0, 0], instance.get_eps2_fs_bh().values())
        self.assertEquals([0, 0], instance.get_eps2_bh_bh().values())
        self.assertEquals(-2, instance.get_eps2()['__result']) # Not implemented (would be ambiguous)
        
        self.assertEquals(0, instance.set_eps2_fs_fs(0.2))
        self.assertEquals([0.2, 0], instance.get_eps2_fs_fs().values())
        self.assertEquals(0, instance.set_eps2_fs_bh(0.3))
        self.assertEquals([0.3, 0], instance.get_eps2_fs_bh().values())
        self.assertEquals(0, instance.set_eps2_bh_bh(0.4))
        self.assertEquals([0.4, 0], instance.get_eps2_bh_bh().values())
        self.assertEquals(-2, instance.set_eps2(0.1)) # Not implemented (would be ambiguous)
        
        self.assertEquals([1.0e-4, 0], instance.get_eta_s().values())
        self.assertEquals([0.1, 0], instance.get_eta_fs().values())
        self.assertEquals([0.4, 0], instance.get_eta_smbh().values())
        self.assertEquals([0.4, 0], instance.get_eta_imbh().values())
        
        self.assertEquals(0, instance.set_eta_s(0.01))
        self.assertEquals([0.01, 0], instance.get_eta_s().values())
        self.assertEquals(0, instance.set_eta_fs(0.02))
        self.assertEquals([0.02, 0], instance.get_eta_fs().values())
        self.assertEquals(0, instance.set_eta_smbh(0.03))
        self.assertEquals([0.03, 0], instance.get_eta_smbh().values())
        self.assertEquals(0, instance.set_eta_imbh(0.04))
        self.assertEquals([0.04, 0], instance.get_eta_imbh().values())
        
        self.assertEquals(0, instance.set_max_relative_energy_error(1.0e-6))
        self.assertEquals([1.0e-6, 0], instance.get_max_relative_energy_error().values())
        self.assertEquals(0, instance.set_maximum_timestep(1.0e-6))
        self.assertEquals([1.0e-6, 0], instance.get_maximum_timestep().values())
        
        self.assertEquals(0, instance.set_smbh_mass(0.1))
        self.assertEquals([0.1, 0], instance.get_smbh_mass().values())
        
        self.assertEquals([0, 0], instance.get_include_smbh_flag().values())
        self.assertEquals(0, instance.set_include_smbh_flag(1))
        self.assertEquals([1, 0], instance.get_include_smbh_flag().values())
        self.assertEquals([1, 0], instance.get_calculate_postnewtonian().values())
        self.assertEquals(0, instance.set_calculate_postnewtonian(0))
        self.assertEquals([0, 0], instance.get_calculate_postnewtonian().values())
        
        self.assertEquals(0, instance.commit_parameters())
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test6(self):
        print "Test MI6Interface evolve_model, single particle (+SMBH)"
        instance = MI6Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_include_smbh_flag(1))
        self.assertEquals(0, instance.set_calculate_postnewtonian(0))
        self.assertEquals(0, instance.commit_parameters())
        
        # Set up a light particle on a circular orbit around the SMBH:
        self.assertEquals([1, 0], instance.new_particle(1.0e-4,  1.0, 0, 0,  0, 1.0, 0, 0.001).values())
        self.assertEquals(0, instance.commit_particles())
        
        self.assertEquals(0, instance.evolve_model(math.pi)) # half an orbit
        for result, expected in zip(instance.get_position(1).values(), [-1.0, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 4)
        
        self.assertEquals(0, instance.evolve_model(2 * math.pi)) # full orbit
        for result, expected in zip(instance.get_position(1).values(), [1.0, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 4)
        
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test7(self):
        print "Test MI6Interface evolve_model, 2 particles orbiting the SMBH"
        instance = MI6Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_include_smbh_flag(1))
        self.assertEquals(0, instance.set_calculate_postnewtonian(0))
        self.assertEquals(0, instance.commit_parameters())
        # Set up a light binary on circular orbits around the SMBH:
        mass = 1.0e-4
        dv = 0.5 * (math.sqrt(1.0 + 0.5*mass) - 1.0)
        self.assertEquals([1, 0], instance.new_particle(mass,  1.0, 0, 0,  0, 1.0+dv, 0, 0.001).values())
        self.assertEquals([2, 0], instance.new_particle(mass, -1.0, 0, 0,  0,-1.0-dv, 0, 0.001).values())
        self.assertEquals(0, instance.commit_particles())
        
        P = 2 * math.pi / (1 + dv)
        self.assertEquals(0, instance.evolve_model(P / 2)) # half an orbit
        for result, expected in zip(instance.get_position(1).values(), [-1.0, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 4)
        
        self.assertEquals(0, instance.evolve_model(P)) # full orbit
        for result, expected in zip(instance.get_position(1).values(), [1.0, 0.0, 0.0, 0]):
            self.assertAlmostEquals(result, expected, 4)
        
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    
    def test8(self):
        print "Test MI6Interface evolve_model, binary (+SMBH) --> accretion!?"
        instance = MI6Interface(**default_options)
        self.assertEquals(0, instance.initialize_code())
        self.assertEquals(0, instance.set_include_smbh_flag(1))
        self.assertEquals(0, instance.set_calculate_postnewtonian(0))
        self.assertEquals(0, instance.commit_parameters())
        # Set up a light binary on circular orbits around the SMBH:
        mass = 1.0e-6
        self.assertEquals([1, 0], instance.new_particle(mass,  1.0, 0, 0,  0, 1.0, 0, 0.001).values())
        self.assertEquals([2, 0], instance.new_particle(mass, -1.0, 0, 0,  0,-1.0, 0, 0.001).values())
        self.assertEquals(0, instance.commit_particles())
        
        self.assertEquals(0, instance.evolve_model(0.00001))
        result1 = instance.get_position(1).values()
        result2 = instance.get_position(2).values()
        print result1
        print result2
        
        self.assertEquals(0, instance.cleanup_code())
        instance.stop()
    


class TestMI6(TestWithMPI):
    
    default_converter = nbody_system.nbody_to_si(1.0e4 | units.MSun, 1.0 | units.AU)
    
    def new_sun_earth_system(self):
        particles = Particles(2)
        particles.mass = [1.0, 3.0037e-6] | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * particles.total_mass() / (1.0 | units.AU)).sqrt()
        return particles
    
    def test1(self):
        print "Testing MI6 initialization"
        instance = MI6(self.default_converter, **default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.cleanup_code()
        instance.stop()
    
    def test2(self):
        print "Testing MI6 parameters"
        instance = MI6(self.default_converter, **default_options)
        instance.initialize_code()
        
        self.assertEquals(instance.parameters.epsilon_squared, 
            instance.unit_converter.to_si(0.0 | nbody_system.length**2))
        self.assertEquals(instance.parameters.timestep_parameter, 0.1)
        
        for par, value in [('epsilon_squared_star_star', 0.0 | nbody_system.length**2), 
                ('epsilon_squared_star_blackhole', 0.0 | nbody_system.length**2), 
                ('epsilon_squared_blackhole_blackhole', 0.0 | nbody_system.length**2),
                ('initial_timestep_parameter', 1.0e-4),
                ('timestep_parameter_stars', 0.1),
                ('timestep_parameter_supermassive_black_holes', 0.4),
                ('timestep_parameter_intermediate_mass_black_holes', 0.4),
                ('max_relative_energy_error', 5.0e-5),
                ('maximum_timestep', 1.0/1024.0 | nbody_system.time),
                ('smbh_mass', 1.0 | nbody_system.mass),
                ('lightspeed', 1.0 | nbody_system.speed)]:
            self.assertEquals(instance.unit_converter.to_si(value), 
                getattr(instance.parameters, par))
                
            if hasattr(value, 'unit'):
                new_value = 3.0 | value.unit
            else:
                new_value = 3.0
                
            setattr(instance.parameters, par, new_value)
            self.assertEquals(instance.unit_converter.to_si(new_value),
                getattr(instance.parameters, par))
        
        # epsilon_squared is an alias for epsilon_squared_star_star, so epsilon_squared also has become 3:
        self.assertEquals(instance.parameters.epsilon_squared, 
            instance.unit_converter.to_si(3.0 | nbody_system.length**2))
        instance.parameters.epsilon_squared = 0.1 | nbody_system.length**2
        self.assertEquals(instance.parameters.epsilon_squared, 
            instance.unit_converter.to_si(0.1 | nbody_system.length**2))
        # timestep_parameter is an alias for timestep_parameter_stars, so timestep_parameter also has become 3:
        self.assertEquals(instance.parameters.timestep_parameter, 3.0)
        instance.parameters.timestep_parameter = 0.01
        self.assertEquals(instance.parameters.timestep_parameter, 0.01)
        
        self.assertEquals(instance.parameters.include_smbh, False)
        instance.parameters.include_smbh = True
        self.assertEquals(instance.parameters.include_smbh, True)
        self.assertEquals(instance.parameters.calculate_postnewtonian, True)
        instance.parameters.calculate_postnewtonian = False
        self.assertEquals(instance.parameters.calculate_postnewtonian, False)
        
        self.assertEquals(instance.parameters.drink, "Vodka martini. Shaken, not stirred.")
        
        instance.stop()
    
    def test3(self):
        print "Testing MI6 particles"
        instance = MI6(self.default_converter, **default_options)
        instance.initialize_code()
        instance.commit_parameters()
        instance.particles.add_particles(self.new_sun_earth_system())
        instance.commit_particles()
        
        self.assertAlmostEquals(instance.particles.mass, [1.0, 3.0037e-6] | units.MSun)
        self.assertAlmostEquals(instance.particles.radius, 1.0 | units.RSun)
        self.assertAlmostEquals(instance.particles.position, 
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU)
        self.assertAlmostEquals(instance.particles.velocity, 
            [[0.0, 0.0, 0.0], [0.0, 29.7885, 0.0]] | units.km / units.s, 3)
        
        instance.cleanup_code()
        instance.stop()
    
    def test4(self):
        print "Testing MI6 evolve_model, 2 particles orbiting the SMBH"
        particles = Particles(2)
        particles.mass = 1.0 | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = ((constants.G * (10001.0 | units.MSun) / (1.0 | units.AU)).sqrt() + 
                           (constants.G * (10000.0 | units.MSun) / (1.0 | units.AU)).sqrt())
        particles.move_to_center()
        print particles
        
        instance = MI6(self.default_converter, **default_options)
        instance.initialize_code()
        instance.parameters.include_smbh = True
        instance.parameters.lightspeed = constants.c
        instance.commit_parameters()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        primary = instance.particles[0]
        
        P = 2 * math.pi * primary.x / primary.vy
        P_corrected = (P) * (2.0 / (1.0 + math.sqrt(1.0001)))
        
        position_at_start = primary.position.x
        instance.evolve_model(P_corrected / 4.0)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.y, 3)
        
        instance.evolve_model(P_corrected / 2.0)
        self.assertAlmostRelativeEqual(position_at_start, -primary.position.x, 3)
        
        instance.evolve_model(P_corrected)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.x, 3)
        
        instance.cleanup_code()
        instance.stop()
    
    def test5(self):
        print "Testing MI6 evolve_model, 2 particles, no SMBH"
        particles = Particles(2)
        particles.mass = 1.0 | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * (2.0 | units.MSun) / (2.0 | units.AU)).sqrt()
        particles.move_to_center()
        print particles
        
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = MI6(converter, **default_options)
        instance.initialize_code()
        instance.parameters.smbh_mass = 0.0 | units.MSun
        instance.commit_parameters()
        instance.particles.add_particles(particles)
        instance.commit_particles()
        primary = instance.particles[0]
        
        P = 2 * math.pi * primary.x / primary.vy
        
        position_at_start = primary.position.x
        instance.evolve_model(P / 4.0)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.y, 3)
        
        instance.evolve_model(P / 2.0)
        self.assertAlmostRelativeEqual(position_at_start, -primary.position.x, 3)
        
        instance.evolve_model(P)
        self.assertAlmostRelativeEqual(position_at_start, primary.position.x, 3)
        
        instance.cleanup_code()
        instance.stop()
    
    def test6(self):
        print "Testing MI6 evolve_model, earth-sun system, no SMBH"
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        instance = MI6(converter, **default_options)
        instance.initialize_code()
        instance.parameters.smbh_mass = 0.0 | units.MSun
        instance.commit_parameters()
        instance.particles.add_particles(self.new_sun_earth_system())
        instance.commit_particles()
        earth = instance.particles[1]
        
        position_at_start = earth.position.x
        instance.evolve_model(0.25 | units.yr)
        self.assertAlmostRelativeEqual(position_at_start, earth.position.y, 3)
        
        instance.evolve_model(0.5 | units.yr)
        self.assertAlmostRelativeEqual(position_at_start, -earth.position.x, 3)
        
        instance.evolve_model(1.0 | units.yr)
        self.assertAlmostRelativeEqual(position_at_start, earth.position.x, 3)
        
        instance.cleanup_code()
        instance.stop()
    
    def test7(self):
        print "Testing effect of MI6 parameter epsilon_squared"
        converter = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        particles = self.new_sun_earth_system()
        particles.rotate(0.0, 0.0, -math.pi/4)
        particles.move_to_center()
        
        tan_initial_direction = particles[1].vy/particles[1].vx
        self.assertAlmostEquals(tan_initial_direction, math.tan(math.pi/4))
        tan_final_direction =  []
        for log_eps2 in range(-9,10,2):
            instance = MI6(converter, **default_options)
            instance.initialize_code()       
            instance.parameters.epsilon_squared = 10.0**log_eps2 | units.AU ** 2
            instance.parameters.smbh_mass = 0.0 | units.MSun
            instance.commit_parameters()
            instance.particles.add_particles(particles)
            instance.commit_particles()
            instance.evolve_model(0.25 | units.yr)
            tan_final_direction.append(instance.particles[1].velocity[1]/
                instance.particles[1].velocity[0])
            instance.cleanup_code()
            instance.stop()
        # Small values of epsilon_squared should result in normal earth-sun dynamics: rotation of 90 degrees
        self.assertAlmostEquals(tan_final_direction[0], math.tan(3 * math.pi / 4.0), 2)
        # Large values of epsilon_squared should result in ~ no interaction
        self.assertAlmostEquals(tan_final_direction[-1], tan_initial_direction, 2)
        # Outcome is most sensitive to epsilon_squared when epsilon_squared = d(earth, sun)^2
        delta = [abs(tan_final_direction[i+1]-tan_final_direction[i]) for i in range(len(tan_final_direction)-1)]
        self.assertEquals(delta[len(tan_final_direction)/2 -1], max(delta))
    
    def test8(self):
        print "Testing MI6 get_gravity_at_point and get_potential_at_point"
        instance = MI6(**default_options)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        instance.parameters.smbh_mass = 0.0 | nbody_system.mass
        
        particles = Particles(2)
        particles.mass = 1.0 | nbody_system.mass
        particles.radius =  0.0 | nbody_system.length
        particles.position = [[0.0,0.0,0.0], [2.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        
        zero = 0.0 | nbody_system.length
        fx, fy, fz = instance.get_gravity_at_point(zero, 1.0 | nbody_system.length, zero, zero)
        self.assertAlmostEqual(fx, 0.0 | nbody_system.acceleration)
        self.assertAlmostEqual(fy, 0.0 | nbody_system.acceleration)
        self.assertAlmostEqual(fz, 0.0 | nbody_system.acceleration)

        for x in (0.25, 0.5, 0.75):
            x0 = x | nbody_system.length
            x1 = (2.0 - x) | nbody_system.length
            potential0 = instance.get_potential_at_point(zero, x0, zero, zero)
            potential1 = instance.get_potential_at_point(zero, x1, zero, zero)
            fx0, fy0, fz0 = instance.get_gravity_at_point(zero, x0, zero, zero)
            fx1, fy1, fz1 = instance.get_gravity_at_point(zero, x1, zero, zero)
            
            self.assertAlmostEqual(fy0, 0.0 | nbody_system.acceleration)
            self.assertAlmostEqual(fz0, 0.0 | nbody_system.acceleration)
            self.assertAlmostEqual(fy1, 0.0 | nbody_system.acceleration)
            self.assertAlmostEqual(fz1, 0.0 | nbody_system.acceleration)
            
            self.assertAlmostEqual(fx0, -1.0 * fx1)
            fx = (-1.0 / (x0**2) + 1.0 / (x1**2)) * (1.0 | nbody_system.length ** 3 / nbody_system.time ** 2)
            self.assertAlmostEqual(fx, fx0)
            self.assertAlmostEqual(potential0, potential1)
        
        instance.stop()
    
    def test9(self):
        print "Testing MI6 evolve_model and getters energy, plummer sphere, no SMBH"
        converter = nbody_system.nbody_to_si(1.0e2 | units.MSun, 1.0 | units.parsec)
        instance = MI6(converter, **default_options)
        instance.initialize_code()
        instance.parameters.smbh_mass = 0.0 | units.MSun
        instance.commit_parameters()
        numpy.random.seed(987654321)
        instance.particles.add_particles(new_plummer_model(100, convert_nbody=converter))
        instance.commit_particles()
        
        kinetic_energy = instance.kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy, 2.12292810174e+37 | units.J, 10)
        self.assertAlmostRelativeEqual(potential_energy, -4.33511391248e+37 | units.J, 10)
        
        initial_total_energy = kinetic_energy + potential_energy
        instance.evolve_model(0.1 | nbody_system.time)
        kinetic_energy = instance.kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy, 2.13633848369e+37 | units.J, 5)
        self.assertAlmostRelativeEqual(potential_energy, -4.34851806763e+37 | units.J, 5)
        
        self.assertAlmostRelativeEqual(potential_energy + kinetic_energy, initial_total_energy, 5)
        
        instance.cleanup_code()
        instance.stop()
    
    def test10(self):
        print "Testing MI6 collision_detection"
        particles = Particles(7)
        particles.mass = 0.00000001 | nbody_system.mass
        particles.radius = 0.0001 | nbody_system.length
        particles.x = [-101.0, -100.0, -0.5, 0.5, 100.0, 101.0, 104.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [[2, 0, 0], [-2, 0, 0]]*3 + [[-4, 0, 0]] | nbody_system.speed
        
        instance = MI6(**default_options)
        instance.initialize_code()
        instance.parameters.set_defaults()
        instance.particles.add_particles(particles)
        collisions = instance.stopping_conditions.collision_detection
        collisions.enable()
        instance.evolve_model(1.0 | nbody_system.time)
        
        self.assertTrue(collisions.is_set())
        self.assertAlmostEquals(instance.model_time, -(particles[0].x-particles[1].x)/(particles[0].vx-particles[1].vx), 3)
        self.assertEquals(len(collisions.particles(0)), 3)
        self.assertEquals(len(collisions.particles(1)), 3)
        self.assertEquals(len(particles - collisions.particles(0) - collisions.particles(1)), 1)
        self.assertEquals(abs(collisions.particles(0).x - collisions.particles(1).x) < 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True, True, True])
        
        sticky_merged = Particles(len(collisions.particles(0)))
        sticky_merged.mass = collisions.particles(0).mass + collisions.particles(1).mass
        sticky_merged.radius = collisions.particles(0).radius
        for p1, p2, merged in zip(collisions.particles(0), collisions.particles(1), sticky_merged):
            merged.position = (p1 + p2).center_of_mass()
            merged.velocity = (p1 + p2).center_of_mass_velocity()
        
        print instance.model_time
        print instance.particles
        instance.particles.remove_particles(collisions.particles(0) + collisions.particles(1))
        instance.particles.add_particles(sticky_merged)
        
        instance.evolve_model(1.0 | nbody_system.time)
        print
        print instance.model_time
        print instance.particles
        self.assertTrue(collisions.is_set())
        self.assertAlmostEquals(instance.model_time, -(particles[6].x-0.5*(particles[4].x+particles[5].x))/particles[6].vx, 3)
        self.assertEquals(len(collisions.particles(0)), 1)
        self.assertEquals(len(collisions.particles(1)), 1)
        self.assertEquals(len(instance.particles - collisions.particles(0) - collisions.particles(1)), 2)
        self.assertEquals(abs(collisions.particles(0).x - collisions.particles(1).x) < 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True])
        instance.stop()
    
    def test11(self):
        print "Testing MI6 properties"
        numpy.random.seed(12345)
        particles = new_plummer_model(100, do_scale=True)
        particles.position += [1, 2, 3] | nbody_system.length
        cluster_velocity = [4, 5, 6] | nbody_system.speed
        particles.velocity += cluster_velocity
        external_kinetic_energy = (0.5 | nbody_system.mass) * cluster_velocity.length_squared()
        
        instance = MI6(**default_options)
        instance.particles.add_particles(particles)
        
        kinetic_energy = instance.kinetic_energy - external_kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy, 0.25 | nbody_system.energy, 10)
        self.assertAlmostRelativeEqual(potential_energy, -0.5 | nbody_system.energy, 10)
        self.assertAlmostRelativeEqual(instance.total_mass, 1.0 | nbody_system.mass, 10)
        self.assertAlmostRelativeEqual(instance.center_of_mass_position, 
            [1, 2, 3] | nbody_system.length, 10)
        self.assertAlmostRelativeEqual(instance.center_of_mass_velocity, 
            [4, 5, 6] | nbody_system.speed, 10)
        initial_total_energy = kinetic_energy + potential_energy
        
        instance.evolve_model(0.1 | nbody_system.time)
        self.assertAlmostRelativeEqual(instance.model_time, 0.1 | nbody_system.time)
        kinetic_energy = instance.kinetic_energy - external_kinetic_energy
        potential_energy = instance.potential_energy
        self.assertAlmostRelativeEqual(kinetic_energy+potential_energy, -0.25 | nbody_system.energy, 3)
        self.assertAlmostRelativeEqual(instance.total_mass, 1.0 | nbody_system.mass, 3)
        self.assertAlmostRelativeEqual(instance.center_of_mass_position, 
            [1.4, 2.5, 3.6] | nbody_system.length, 3)
        self.assertAlmostRelativeEqual(instance.center_of_mass_velocity, 
            [4, 5, 6] | nbody_system.speed, 3)
        
        instance.cleanup_code()
        instance.stop()
    
