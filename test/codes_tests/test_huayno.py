from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.community.huayno.interface import HuaynoInterface, Huayno

from amuse.units import nbody_system
from amuse.units import units, constants
from amuse import datamodel
from amuse.ic import plummer
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

def elements(starmass,x,y,z,vx,vy,vz,G=constants.G):
    mu=G*starmass
    r=(x**2+y**2+z**2)**0.5
    v2=(vx**2+vy**2+vz**2)

    e=v2/2-mu/r

    a=-mu/2/e

    hx=y*vz-z*vy
    hy=z*vx-x*vz
    hz=x*vy-y*vx

    rdotv=x*vx+y*vy+z*vz

    ex=v2*x/mu-rdotv*vx/mu-x/r
    ey=v2*y/mu-rdotv*vy/mu-y/r
    ez=v2*z/mu-rdotv*vz/mu-z/r

    h2=(hx**2+hy**2+hz**2)

    eps=(1-h2/mu/a)**0.5

    return a,eps

def energy(mass,eps,orbiters):
   return (0.5*orbiters.velocity.lengths()**2-constants.G*mass/(orbiters.position.lengths()**2+eps**2)**0.5)

def angular_momentum(orbiters):
   return orbiters.position.cross(orbiters.velocity).lengths()

class TestHuaynoInterface(TestWithMPI):
    
    def test1(self):
        instance = HuaynoInterface()
        instance.initialize_code()
    
        res1 = instance.new_particle(mass = 11.0, radius = 2.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        res2 = instance.new_particle(mass = 21.0, radius = 5.0, x = 10.0, y = 0.0, z = 0.0, vx = 10.0, vy = 0.0, vz = 0.0)
        
        self.assertEqual(0, res1['index_of_the_particle'])
        self.assertEqual(1, res2['index_of_the_particle'])
    
        retrieved_state1 = instance.get_state(0)
        retrieved_state2 = instance.get_state(1)
    
        self.assertEqual(11.0,  retrieved_state1['mass'])
        self.assertEqual(21.0,  retrieved_state2['mass'])
        self.assertEqual(0.0,  retrieved_state1['x'])
        self.assertEqual(10.0,  retrieved_state2['x'])
    
        instance.cleanup_code()
        instance.stop()

    def test2(self):
        instance = HuaynoInterface()
        instance.initialize_code()

        for i in [0, 1, 2]:
            temp_particle = instance.new_particle(mass = i, radius = 1.0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
            self.assertEqual(i, temp_particle['index_of_the_particle'])
            
        instance.delete_particle(1)
      
        self.assertEqual(2, instance.get_number_of_particles()['number_of_particles'])
        
        self.assertEqual(0, instance.get_index_of_first_particle()['index_of_the_particle'])
        
        self.assertEqual(2, instance.get_index_of_next_particle(0)['index_of_the_next_particle'])
        self.assertEqual(0, instance.get_index_of_next_particle(0)['__result'])
        self.assertEqual(-2, instance.get_index_of_next_particle(1)['__result'])
        self.assertEqual(1, instance.get_index_of_next_particle(2)['__result'])
        instance.cleanup_code()
        instance.stop()
        
    def test3(self):
        huayno = HuaynoInterface()
        huayno.set_eps2(0.101)
        self.assertEqual(0.101, huayno.get_eps2()['eps2'])
        huayno.set_eps2( 0.110)
        self.assertEqual(0.110, huayno.get_eps2()['eps2'])
        huayno.cleanup_code()
        huayno.stop()

    def test5(self):
        huayno = HuaynoInterface()
        huayno.initialize_code()
        
        huayno.new_particle([10,20],[0,0],[0,0], [0,0], [0,0], [0,0], [0,0],[1,1])
        retrieved_state = huayno.get_state(0)
        
        self.assertEqual(10.0,  retrieved_state['mass'])
        self.assertEqual(1, retrieved_state['radius'])
    
        retrieved_state = huayno.get_state([0,1])
        self.assertEqual(20.0,  retrieved_state['mass'][1])
        self.assertEqual(huayno.get_number_of_particles()['number_of_particles'], 2)
        huayno.cleanup_code() 
        huayno.stop()
       

class TestHuayno(TestWithMPI):
    def new_system_of_sun_and_earth(self):
        particles = datamodel.Particles(2)
        particles.mass = [1.0, 3.0037e-6] | units.MSun
        particles.radius = 1.0 | units.RSun
        particles.position = [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]] | units.AU
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.km / units.s
        particles[1].vy = (constants.G * particles.total_mass() / (1.0 | units.AU)).sqrt()
        particles.move_to_center()
        return particles
    
    def test1(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        huayno = Huayno(convert_nbody)
        huayno.initialize_code()
        huayno.parameters.epsilon_squared = 0.0 | units.AU**2
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        
        huayno.particles.add_particles(stars)
        
        huayno.evolve_model(365.0 | units.day)
        huayno.particles.copy_values_of_all_attributes_to(stars)
        
        position_at_start = earth.position.value_in(units.AU)[0]
        position_after_full_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(position_at_start, position_after_full_rotation, 6)
        
        huayno.evolve_model(365.0 + (365.0 / 2) | units.day)
        
        huayno.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 2)
                
        huayno.evolve_model(365.0 + (365.0 / 2) + (365.0 / 4)  | units.day)
        
        huayno.particles.copy_values_of_all_attributes_to(stars)
        position_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        self.assertAlmostEqual(-position_at_start, position_after_half_a_rotation, 3)
        
        huayno.cleanup_code()
        
        huayno.stop()
        

    def test2(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)
    
        instance = Huayno(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.0 | units.AU**2
        
        stars = self.new_system_of_sun_and_earth()
        earth = stars[1]
        instance.particles.add_particles(stars)
    
        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.particles.copy_values_of_all_attributes_to(stars)
            stars.savepoint()
        
        if HAS_MATPLOTLIB:
            figure = pyplot.figure()
            plot = figure.add_subplot(1,1,1)
            
            x_points = earth.get_timeline_of_attribute("x")
            y_points = earth.get_timeline_of_attribute("y")
            
            x_points_in_AU = [t_x[1].value_in(units.AU) for t_x in x_points]
            y_points_in_AU = [t_x1[1].value_in(units.AU) for t_x1 in y_points]
            
            plot.scatter(x_points_in_AU,y_points_in_AU, color = "b", marker = 'o')
            
            plot.set_xlim(-1.5, 1.5)
            plot.set_ylim(-1.5, 1.5)
               
            
            test_results_path = self.get_path_to_results()
            output_file = os.path.join(test_results_path, "huayno-earth-sun2.svg")
            figure.savefig(output_file)
              
        instance.cleanup_code()
        del instance

    def test3(self):
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 149.5e6 | units.km)

        instance = Huayno(convert_nbody)
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.00001 | units.AU**2
        
        stars = datamodel.Stars(2)
        star1 = stars[0]
        star2 = stars[1]

        star1.mass = units.MSun(1.0)
        star1.position = units.AU(numpy.array((-1.0,0.0,0.0)))
        star1.velocity = units.AUd(numpy.array((0.0,0.0,0.0)))
        star1.radius = units.RSun(1.0)

        star2.mass = units.MSun(1.0)
        star2.position = units.AU(numpy.array((1.0,0.0,0.0)))
        star2.velocity = units.AUd(numpy.array((0.0,0.0,0.0)))
        star2.radius = units.RSun(100.0)
        
        instance.particles.add_particles(stars)
    
        for x in range(1,2000,10):
            instance.evolve_model(x | units.day)
            instance.particles.copy_values_of_all_attributes_to(stars)
            stars.savepoint()
        
    
    def test4(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = Huayno(convert_nbody)
        instance.initialize_code()
        
        particles = datamodel.Particles(2)
        self.assertEqual(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 2)
        
        instance.particles.mass =  [17.0, 33.0] | units.kg
        
        
        self.assertEqual(instance.get_mass(0), 17.0| units.kg) 
        self.assertEqual(instance.get_mass(1), 33.0| units.kg)  

    def test5(self):
        convert_nbody = nbody_system.nbody_to_si(5.0 | units.kg, 10.0 | units.m)

        instance = Huayno(convert_nbody)
        instance.initialize_code()
        
        particles = datamodel.Particles(2)
        self.assertEqual(len(instance.particles), 0)
        
        particles.mass = [15.0, 30.0] | units.kg
        particles.radius =  [10.0, 20.0] | units.m
        particles.position = [[10.0, 20.0, 30.0], [20.0, 40.0, 60.0]] | units.m
        particles.velocity = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]] | units.m / units.s

        
        instance.particles.add_particles(particles)
        self.assertEqual(len(instance.particles), 2)
        
        instance.set_state(1, 16|units.kg, 20.0|units.m, 40.0|units.m, 60.0|units.m, 
                                 1.0|units.ms, 1.0|units.ms, 1.0|units.ms)
        
        curr_state =  instance.get_state(1)
        for expected, actural in zip((16|units.kg, 20.0|units.m, 40.0|units.m, 60.0|units.m, 
                                 1.0|units.ms, 1.0|units.ms, 1.0|units.ms, 0 | units.m), curr_state):
            self.assertAlmostRelativeEquals(actural,expected)
        
        instance.set_state(1, 16|units.kg, 20.0|units.m, 40.0|units.m, 60.0|units.m, 
                                 1.0|units.ms, 1.0|units.ms, 1.0|units.ms , 20.0|units.m)
        
        curr_state =  instance.get_state(1)
        for expected, actural in zip((16|units.kg, 20.0|units.m, 40.0|units.m, 60.0|units.m, 
                                 1.0|units.ms, 1.0|units.ms, 1.0|units.ms, 20 | units.m), curr_state):
            self.assertAlmostRelativeEquals(actural,expected)
        
        
    
    def test6(self):
        print("Test6: Testing Huayno parameters")
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.yr, 1.0 | units.AU)
        instance = Huayno(convert_nbody)
        
        (value, error) = instance.legacy_interface.get_eps2()
        self.assertEqual(0, error)
        self.assertEqual(0.0, value)
        self.assertAlmostEqual(0.0 | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
        for x in [0.01, 0.1, 0.2]:
            instance.parameters.epsilon_squared = x | units.AU**2
            self.assertAlmostEqual(x | units.AU**2, instance.parameters.epsilon_squared, in_units=units.AU**2)
                
        (value, error) = instance.legacy_interface.get_time()
        self.assertEqual(0, error)
        self.assertEqual(0.0, value)
        self.assertAlmostEqual(0.0 | units.yr, instance.parameters.begin_time, in_units=units.yr)
        for x in [1.0, 10.0, 100.0]:
            instance.parameters.begin_time = x | units.yr
            self.assertAlmostEqual(x | units.yr, instance.parameters.begin_time, in_units=units.yr)
        instance.stop()
    
    def test7(self):
        print("Test7: Testing effect of Huayno parameter epsilon_squared")
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        
        particles = datamodel.Particles(2)
        sun = particles[0]
        sun.mass = 1.0 | units.MSun
        sun.position = [0.0, 0.0, 0.0] | units.AU
        sun.velocity = [0.0, 0.0, 0.0] | units.AU / units.yr
        sun.radius = 1.0 | units.RSun

        earth = particles[1]
        earth.mass = 5.9736e24 | units.kg
        earth.radius = 6371.0 | units.km
        earth.position = [0.0, 1.0, 0.0] | units.AU
        earth.velocity = [2.0*numpy.pi, -0.0001, 0.0] | units.AU / units.yr
        
        initial_direction = math.atan((earth.velocity[0]/earth.velocity[1]))
        final_direction = []
        for log_eps2 in range(-9,10,2):
            instance = Huayno(convert_nbody)
            instance.initialize_code()
            instance.parameters.epsilon_squared = 10.0**log_eps2 | units.AU ** 2
            instance.particles.add_particles(particles)
            instance.commit_particles()
            instance.evolve_model(0.25 | units.yr)
            final_direction.append(math.atan((instance.particles[1].velocity[0]/
                instance.particles[1].velocity[1])))
            instance.stop()
        # Small values of epsilon_squared should result in normal earth-sun dynamics: rotation of 90 degrees
        self.assertAlmostEqual(abs(final_direction[0]), abs(initial_direction+math.pi/2.0), 2)
        # Large values of epsilon_squared should result in ~ no interaction
        self.assertAlmostEqual(final_direction[-1], initial_direction, 2)
        # Outcome is most sensitive to epsilon_squared when epsilon_squared = d(earth, sun)^2
        delta = [abs(final_direction[i+1]-final_direction[i]) for i in range(len(final_direction)-1)]
        self.assertEqual(delta[len(final_direction)//2 -1], max(delta))
        
        
    def test13(self):
        particles = plummer.new_plummer_model(31)
       
        instance = Huayno(number_of_workers=1)#, debugger="xterm")
        instance.initialize_code()
        instance.parameters.epsilon_squared = 0.01 | nbody_system.length ** 2
        instance.particles.add_particles(particles)
        
        instance.evolve_model(0.1 | nbody_system.time)
        instance.synchronize_model()
        expected_positions = instance.particles.position
        instance.stop()
    
    def test14(self):
        import hashlib

        numpy.random.seed(123456)
        particles = plummer.new_plummer_model(64)
        sha=hashlib.sha1()

        class inttypes(object):
            SHARED2=1
            EXTRAPOLATE=5
            PASS_KDK=2
            PASS_DKD=7
            HOLD_KDK=3
            HOLD_DKD=8
            PPASS_DKD=9
            BRIDGE_KDK=4
            BRIDGE_DKD=10
            CC=11
            CC_KEPLER=12
            OK=13
            KEPLER=14
            SHARED4=15
            SHARED6=18
            SHARED8=19
            SHARED10=20
            SHAREDBS=21
           
            @classmethod
            def _list(cls):
                  return set([x for x in cls.__dict__.keys() if not x.startswith('_')])

        for itype in sorted(inttypes._list()):
            if itype in ("KEPLER"): continue
            instance = Huayno()
            instance.parameters.inttype_parameter=getattr(Huayno.inttypes,itype)
            instance.particles.add_particles(particles)
            instance.evolve_model(0.125 | nbody_system.time)
            part_out= instance.particles.copy()
            position = part_out.position.number
            if hasattr(position,'tobytes'):
                as_bytes = position.tobytes()
            else:
                as_bytes = numpy.array(position.data, copy=True, order='C')
            sha.update(as_bytes)
            
            instance.stop()
        
        # this result is probably dependent on system architecture hence no good for assert
        print() 
        print(sha.hexdigest())
        print("8e71f9441578a43af4af927943577ad2c4130e4c")
        
    def test15(self):
        particles = plummer.new_plummer_model(512)
        expected_positions = None
        for mode in ["cpu", "openmp", "opencl"]:
            try:
                instance = Huayno(mode=mode, number_of_workers=1)#, debugger="xterm")
            except:
                print("Running huayno with mode=", mode, " was unsuccessful.")
                continue
            else:
                print("Running huayno with mode=", mode, "... ")
                
            instance.initialize_code()
            instance.parameters.epsilon_squared = 0.01 | nbody_system.length ** 2
            instance.particles.add_particles(particles)
            
            instance.evolve_model(0.2 | nbody_system.time)
            instance.synchronize_model()
            if expected_positions is None:
                expected_positions = instance.particles.position
            else:
                self.assertAlmostRelativeEquals(expected_positions, instance.particles.position, 8)
            instance.stop()

    def test16(self):
        instance = Huayno()
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        
        particles = datamodel.Particles(2)
        particles.mass = [1.] | nbody_system.mass
        particles.radius =  [0.0] | nbody_system.length
        particles.position = [[0.0,0.0,0.0],[1.0,0.0,0.0]] | nbody_system.length
        particles.velocity = [[0.0, 0.0, 0.0]] | nbody_system.speed
        instance.particles.add_particles(particles)
        
        zero = 0.0 | nbody_system.length
        
        for x in (0.25, 0.5, 0.75):
            x0 = x | nbody_system.length
            potential0 = instance.get_potential_at_point(zero, x0, zero, zero)
            fx0, fy0, fz0 = instance.get_gravity_at_point(zero, x0, zero, zero)
            
            self.assertAlmostEqual(fy0, 0.0 | nbody_system.acceleration,14)
            self.assertAlmostEqual(fz0, 0.0 | nbody_system.acceleration,14)
            
            fx = (-1.0 / (x0**2)+1.0 / (((1.0|nbody_system.length)-x0)**2)) * (1.0 | nbody_system.length ** 3 / nbody_system.time ** 2)
            self.assertAlmostEqual(fx, fx0,14)
            self.assertAlmostEqual(potential0, -nbody_system.G*(1.|nbody_system.mass)*(1./x0+1./((1.|nbody_system.length)-x0)),14)
        instance.stop()
    
    def _compare_integrator_with_collision_integrator(self, inttype_parameter1, inttype_parameter2):
        numpy.random.seed(12345)
        particles = plummer.new_plummer_model(101)
        instance = Huayno()
        instance.parameters.inttype_parameter = inttype_parameter1
        instance.particles.add_particles(particles)
        instance.evolve_model(0.2 | nbody_system.time)
        expected_position = instance.particles.position
        expected_velocity = instance.particles.velocity
        instance.reset()
        instance.parameters.inttype_parameter = inttype_parameter2
        instance.particles.add_particles(particles)
        instance.evolve_model(0.2 | nbody_system.time)
        self.assertAlmostRelativeEquals(expected_position, instance.particles.position, 8)
        self.assertAlmostRelativeEquals(expected_velocity, instance.particles.velocity, 8)
        instance.stop()
    
    def _run_collision_with_integrator(self, inttype_parameter):
        particles = datamodel.Particles(7)
        particles.mass = 0.001 | nbody_system.mass
        particles.radius = 0.01 | nbody_system.length
        particles.x = [-101.0, -100.0, -0.5, 0.5, 100.0, 101.0, 104.0] | nbody_system.length
        particles.y = 0 | nbody_system.length
        particles.z = 0 | nbody_system.length
        particles.velocity = [[2, 0, 0], [-2, 0, 0]]*3 + [[-4, 0, 0]] | nbody_system.speed
        
        instance = Huayno()
        instance.parameters.inttype_parameter = inttype_parameter
        instance.particles.add_particles(particles)
        collisions = instance.stopping_conditions.collision_detection
        collisions.enable()
        instance.evolve_model(1.0 | nbody_system.time)
        
        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 0.5 | nbody_system.time)
        self.assertEqual(len(collisions.particles(0)), 3)
        self.assertEqual(len(collisions.particles(1)), 3)
        self.assertEqual(len(particles - collisions.particles(0) - collisions.particles(1)), 1)
        self.assertEqual(abs(collisions.particles(0).x - collisions.particles(1).x) <= 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True, True, True])
        
        sticky_merged = datamodel.Particles(len(collisions.particles(0)))
        sticky_merged.mass = collisions.particles(0).mass + collisions.particles(1).mass
        sticky_merged.radius = collisions.particles(0).radius
        for p1, p2, merged in zip(collisions.particles(0), collisions.particles(1), sticky_merged):
            merged.position = (p1 + p2).center_of_mass()
            merged.velocity = (p1 + p2).center_of_mass_velocity()
        
        print(instance.model_time)
        print(instance.particles)
        instance.particles.remove_particles(collisions.particles(0) + collisions.particles(1))
        instance.particles.add_particles(sticky_merged)
        
        instance.evolve_model(1.0 | nbody_system.time)
        print()
        print(instance.model_time)
        print(instance.particles)
        self.assertTrue(collisions.is_set())
        self.assertTrue(instance.model_time < 1.0 | nbody_system.time)
        self.assertEqual(len(collisions.particles(0)), 1)
        self.assertEqual(len(collisions.particles(1)), 1)
        self.assertEqual(len(instance.particles - collisions.particles(0) - collisions.particles(1)), 2)
        self.assertEqual(abs(collisions.particles(0).x - collisions.particles(1).x) <= 
                (collisions.particles(0).radius + collisions.particles(1).radius),
                [True])
        instance.stop()
    
    def test17(self):
        print("Compare the SHARED2 integrator with the collision-detection enabled SHARED2_COLLISIONS integrator")
        self._compare_integrator_with_collision_integrator(Huayno.inttypes.SHARED2_COLLISIONS, Huayno.inttypes.SHARED2)
        print("Testing Huayno collision_detection with SHARED2_COLLISIONS")
        self._run_collision_with_integrator(Huayno.inttypes.SHARED2_COLLISIONS)
    
    def test18(self):
        print("Compare the SHARED4 integrator with the collision-detection enabled SHARED4_COLLISIONS integrator")
        self._compare_integrator_with_collision_integrator(Huayno.inttypes.SHARED4_COLLISIONS, Huayno.inttypes.SHARED4)
        print("Testing Huayno collision_detection with SHARED4_COLLISIONS")
        self._run_collision_with_integrator(Huayno.inttypes.SHARED4_COLLISIONS)
    
    def test19(self):
        print("Compare the SHARED6 integrator with the collision-detection enabled SHARED6_COLLISIONS integrator")
        self._compare_integrator_with_collision_integrator(Huayno.inttypes.SHARED6_COLLISIONS, Huayno.inttypes.SHARED6)
        print("Testing Huayno collision_detection with SHARED6_COLLISIONS")
        self._run_collision_with_integrator(Huayno.inttypes.SHARED6_COLLISIONS)
    
    def test20(self):
        print("Compare the SHARED8 integrator with the collision-detection enabled SHARED8_COLLISIONS integrator")
        self._compare_integrator_with_collision_integrator(Huayno.inttypes.SHARED8_COLLISIONS, Huayno.inttypes.SHARED8)
        print("Testing Huayno collision_detection with SHARED8_COLLISIONS")
        self._run_collision_with_integrator(Huayno.inttypes.SHARED8_COLLISIONS)
    
    def test21(self):
        print("Compare the SHARED10 integrator with the collision-detection enabled SHARED10_COLLISIONS integrator")
        self._compare_integrator_with_collision_integrator(Huayno.inttypes.SHARED10_COLLISIONS, Huayno.inttypes.SHARED10)
        print("Testing Huayno collision_detection with SHARED10_COLLISIONS")
        self._run_collision_with_integrator(Huayno.inttypes.SHARED10_COLLISIONS)
    
    def test22(self):
        print("Testing zero-mass test particles in Huayno, can be used for removing particles when inside recursive evolve loop")
        sun_and_earth = self.new_system_of_sun_and_earth()
        period = (4.0 * math.pi**2 * (1.0 | units.AU)**3 / (constants.G * sun_and_earth.total_mass())).sqrt()
        convert_nbody = nbody_system.nbody_to_si(1.0 | units.MSun, 1.0 | units.AU)
        huayno = Huayno(convert_nbody)
        huayno.parameters.epsilon_squared = 0.0 | units.AU**2
        huayno.parameters.inttype_parameter = huayno.inttypes.SHARED8
        
        test_particle = datamodel.Particle(mass=0|units.MSun, position=[4,0,0]|units.AU, velocity=[0,0,0]|units.kms)
        test_particle.vy = (constants.G * sun_and_earth.total_mass() / (4.0 | units.AU)).sqrt()
        sun_and_earth.add_particle(test_particle)
        huayno.particles.add_particles(sun_and_earth)
        huayno.evolve_model(period)
        self.assertAlmostRelativeEqual(huayno.particles[:2].x, sun_and_earth[:2].x, 13)
        huayno.evolve_model(1.25 * period)
        self.assertAlmostRelativeEqual(huayno.particles[1].y, sun_and_earth[1].x, 13)
        huayno.evolve_model(8.0 * period)
        self.assertAlmostRelativeEqual(huayno.particles.x, sun_and_earth.x, 8)
        huayno.stop()
    
    def test23(self):
        print("testing removing and adding particles repeatedly")
        N=1100
        p1=plummer.new_plummer_model(N)
        p2=plummer.new_plummer_model(N)
        h1=Huayno()
        h1.particles.add_particles(p1[:N//2])
        h1.particles.add_particles(p2[-N//2:])

        h2=Huayno()
        h2.particles.add_particles(p1)
        h2.particles.remove_particles(p1[N//2:])
        h2.particles.add_particles(p2)
        h2.particles.remove_particles(p2[:-N//2])
        
        self.assertEqual(len(h1.particles),len(h2.particles))
        self.assertAlmostEqual(h1.kinetic_energy,h2.kinetic_energy,15)
        self.assertAlmostEqual(h1.potential_energy,h2.potential_energy,15)
        
    def test24(self):
        print("test massless particles/ kepler integrator")
        N=20        
        tend=2.| units.yr
        numpy.random.seed(12345)
        conv=nbody_system.nbody_to_si(4.| units.MSun, 5.|units.AU)
        orbiters=plummer.new_plummer_model(N,conv)
        sun=datamodel.Particle(mass=1.|units.MSun)
        sun.position=[12.3,1.,-.2]|units.AU
        sun.velocity=[10,50.,-20.]|units.kms
        orbiters.mass*=0.
        
        a0,eps0=elements(sun.mass,orbiters.x,orbiters.y,orbiters.z,
                      orbiters.vx,orbiters.vy,orbiters.vz)

        orbiters.position+=sun.position
        orbiters.velocity+=sun.velocity

        pos=dict()
        for inttype in [20,14]:
          code=Huayno(conv)
          code.parameters.inttype_parameter=inttype
          code.parameters.timestep_parameter=0.1
          code.particles.add_particle(sun)
          orbiters2=code.particles.add_particles(orbiters).copy()
          orbiters2.position-=sun.position
          orbiters2.velocity-=sun.velocity
          code.evolve_model(tend)
          a,eps=elements(sun.mass,orbiters2.x,orbiters2.y,orbiters2.z,
                    orbiters2.vx,orbiters2.vy,orbiters2.vz)

          da=abs((a-a0)/a0)
          deps=abs(eps-eps0)/eps0

          dev=numpy.where(da > 1.e-12)[0]
          self.assertEqual( len(dev),0)
          dev=numpy.where(deps > 1.e-12)[0]
          self.assertEqual( len(dev),0)
          pos[inttype]=[orbiters2.x.value_in(units.AU),orbiters2.y.value_in(units.AU),orbiters2.z.value_in(units.AU)]
        self.assertAlmostEqual(pos[20][0],pos[14][0],12)
        self.assertAlmostEqual(pos[20][1],pos[14][1],12)
        self.assertAlmostEqual(pos[20][2],pos[14][2],12)

    def test25(self):
        print("test massless particles/ kepler integrator, smoothed")
        N=10        
        tend=20.| units.yr
        numpy.random.seed(12345)
        conv=nbody_system.nbody_to_si(4.| units.MSun, 5.|units.AU)
        orbiters=plummer.new_plummer_model(N,conv)
        sun=datamodel.Particle(mass=1.|units.MSun)
        sun.position=[0,0,0]|units.AU
        sun.velocity=[0,0,0]|units.kms
        orbiters.mass*=0.
        eps=(5. | units.AU)

        e0=energy(sun.mass,eps,orbiters)
        l0=angular_momentum(orbiters)
        
        pos=dict()
        for inttype in [20,14]:
          code=Huayno(conv)
          code.parameters.inttype_parameter=inttype
          code.parameters.timestep_parameter=0.1
          code.parameters.epsilon_squared=eps**2
          code.particles.add_particle(sun)
          orbiters2=code.particles.add_particles(orbiters)
          code.evolve_model(tend)

          e1=energy(sun.mass,eps,orbiters2)
          l1=angular_momentum(orbiters2)
          de,dl=abs((e1-e0)/e0).max(),abs((l1-l0)/l1).max()
          self.assertTrue( numpy.all(de< 1.e-8))
          self.assertTrue( numpy.all(dl< 1.e-8))

          pos[inttype]=[orbiters2.x.value_in(units.AU),orbiters2.y.value_in(units.AU),orbiters2.z.value_in(units.AU)]
        self.assertAlmostRelativeEqual(pos[20][0],pos[14][0],4) # still not clear why 4
        self.assertAlmostRelativeEqual(pos[20][1],pos[14][1],4)
        self.assertAlmostRelativeEqual(pos[20][2],pos[14][2],4)


    def test26(self):
        print("test massless particles (negative time)")
        N=10        
        tend=-5.| units.yr
        numpy.random.seed(12345)
        conv=nbody_system.nbody_to_si(4.| units.MSun, 5.|units.AU)
        orbiters=plummer.new_plummer_model(N,conv)
        sun=datamodel.Particle(mass=1.|units.MSun)
        sun.position=[0,0,0]|units.AU
        sun.velocity=[0,0,0]|units.kms
        orbiters.mass*=0.
        
        a0,eps0=elements(sun.mass,orbiters.x,orbiters.y,orbiters.z,
                      orbiters.vx,orbiters.vy,orbiters.vz)

        pos=dict()
        for inttype in [20,14]:
          code=Huayno(conv)
          code.parameters.inttype_parameter=inttype
          code.parameters.timestep_parameter=0.1
          code.particles.add_particle(sun)
          orbiters2=code.particles.add_particles(orbiters)
          code.evolve_model(tend)
          a,eps=elements(sun.mass,orbiters2.x,orbiters2.y,orbiters2.z,
                    orbiters2.vx,orbiters2.vy,orbiters2.vz)

          da=abs((a-a0)/a0)
          deps=abs(eps-eps0)/eps0

          dev=numpy.where(da > 1.e-12)[0]
          self.assertEqual( len(dev),0)
          dev=numpy.where(deps > 1.e-12)[0]
          self.assertEqual( len(dev),0)
          pos[inttype]=[orbiters.x.value_in(units.AU),orbiters.y.value_in(units.AU),orbiters.z.value_in(units.AU)]
        self.assertAlmostEqual(pos[20][0],pos[14][0],12)
        self.assertAlmostEqual(pos[20][1],pos[14][1],12)
        self.assertAlmostEqual(pos[20][2],pos[14][2],12)
        
    def test27(self):
        particles = plummer.new_plummer_model(31)

        tend=0.25| nbody_system.time

        instance = Huayno()
        instance.particles.add_particles(particles)        
        instance.evolve_model(tend)
        expected_positions = instance.particles.position
        self.assertEqual(instance.model_time, tend)
        instance.stop()
        
        particles2=particles.copy()
        particles2.velocity*=-1

        instance = Huayno()
        instance.particles.add_particles(particles2)        
        instance.evolve_model(-tend)
        positions = instance.particles.position
        self.assertEqual(instance.model_time, -tend)
        instance.stop()

        self.assertAlmostEqual(positions,expected_positions)

    def test28(self):
        particles = plummer.new_plummer_model(31)

        instance = Huayno()
        instance.particles.add_particles(particles)        
        self.assertAlmostEqual(particles.total_mass(),instance.total_mass)
        self.assertAlmostEqual(particles.center_of_mass(),instance.center_of_mass_position)
        self.assertAlmostEqual(particles.center_of_mass_velocity(),instance.center_of_mass_velocity)

    def test29(self):
        instance = Huayno()
        instance.parameters.epsilon_squared = 0.0 | nbody_system.length**2
        
        particles = datamodel.Particles(2)
        particles.position = ([0,0,0], [1,0,0] )| nbody_system.length
        particles.velocity = ([-1,0,0], [2,0,0] )| nbody_system.speed
        particles.radius = 0| nbody_system.length
        particles.mass = 0.1| nbody_system.mass
                
        instance.particles.add_particles(particles)
        instance.stopping_conditions.out_of_box_detection.enable()
        instance.parameters.stopping_conditions_out_of_box_size = 2 | nbody_system.length
        instance.parameters.stopping_conditions_out_of_box_use_center_of_mass = False
        instance.evolve_model(1 | nbody_system.time)
        self.assertTrue(instance.stopping_conditions.out_of_box_detection.is_set())
        self.assertEqual(len(instance.stopping_conditions.out_of_box_detection.particles(0)), 1)
        self.assertEqual(instance.stopping_conditions.out_of_box_detection.particles(0)[0].key, particles[1].key)
        instance.stop()
        

