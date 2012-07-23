import numpy
from numpy import pi
from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from amuse.community.mikkola.interface import MikkolaInterface
from amuse.community.mikkola.interface import Mikkola

class MikkolaInterfaceTests(TestWithMPI):
    
    def xtest0(self):
        instance = MikkolaInterface()
    
    def test1(self):
#        instance = MikkolaInterface(debugger="ddd")
        instance = MikkolaInterface()
        #instance = MikkolaInterface()
        instance.initialize_code()
        instance.set_lightspeed(1e4)
        instance.commit_parameters()

        res1 = instance.new_particle(mass = 1.0, radius = 0, x = 0.0, y = 0.0, z = 0.0, vx = 0.0, vy = 0.0, vz = 0.0)
        res2 = instance.new_particle(mass = 0.001, radius = 0, x = 1.0, y = 0.0, z = 0.0, vx = 0.0, vy = 1.0, vz = 0.0)
        
        self.assertEquals(1, res1['index_of_the_particle'])
        self.assertEquals(2, res2['index_of_the_particle'])
    
        retrieved_state1 = instance.get_state(1)
        retrieved_state2 = instance.get_state(2)
    
        end_time = 10.0 
        instance.evolve_model(end_time)
        instance.cleanup_code()
        instance.stop()
        
    def test2(self):
        instance = MikkolaInterface()
        instance.initialize_code()
        instance.set_lightspeed(1e4)
        instance.commit_parameters()

        index1, error = instance.new_particle(mass = 1.0, radius = 0, x = 1.0, y = 2.0, z = 3.0, vx = 1.1, vy = 2.2, vz = 3.3)
        self.assertEquals(error, 0)
        index2, error = instance.new_particle(mass = 0.001, radius = 0, x = 3.0, y = 4.0, z = 5.0, vx = 3.3, vy = 4.4, vz = 5.5)
        self.assertEquals(error, 0)
        
        self.assertEquals(1, index1)
        self.assertEquals(2, index2)
        
        mass,error = instance.get_mass(1)
        self.assertEquals(error, 0)
        self.assertEquals(mass, 1.0)
                
        mass,error = instance.get_mass(2)
        self.assertEquals(error, 0)
        self.assertEquals(mass, 0.001)
        
        
        x,y,z,error = instance.get_position(1)
        self.assertEquals(error, 0)
        self.assertEquals(x, 1.0)
        self.assertEquals(y, 2.0)
        self.assertEquals(z, 3.0)
                
        x,y,z,error = instance.get_position(2)
        self.assertEquals(error, 0)
        self.assertEquals(x, 3.0)
        self.assertEquals(y, 4.0)
        self.assertEquals(z, 5.0)
        
        vx,vy,vz,error = instance.get_velocity(1)
        self.assertEquals(error, 0)
        self.assertEquals(vx, 1.1)
        self.assertEquals(vy, 2.2)
        self.assertEquals(vz, 3.3)
    
        instance.cleanup_code()
        instance.stop()
        
    def test3(self):
        instance = MikkolaInterface()
        instance.initialize_code()
        instance.set_lightspeed(1e4)
        instance.commit_parameters()

        index1, error = instance.new_particle(mass = 1.0, radius = 0, x = 1.0, y = 2.0, z = 3.0, vx = 1.1, vy = 2.2, vz = 3.3)
        self.assertEquals(error, 0)
        index2, error = instance.new_particle(mass = 0.001, radius = 0, x = 3.0, y = 4.0, z = 5.0, vx = 3.3, vy = 4.4, vz = 5.5)
        self.assertEquals(error, 0)
        
        self.assertEquals(1, index1)
        self.assertEquals(2, index2)
        
        mass,error = instance.get_mass(1)
        self.assertEquals(error, 0)
        self.assertEquals(mass, 1.0)
                
        mass,error = instance.get_mass(2)
        self.assertEquals(error, 0)
        self.assertEquals(mass, 0.001)
        
        error  = instance.delete_particle(1)
        mass,error = instance.get_mass(1)
        self.assertEquals(error, -1)
        
        index1, error = instance.new_particle(mass = 2.0, radius = 0, x = 1.0, y = 2.0, z = 3.0, vx = 1.1, vy = 2.2, vz = 3.3)
        self.assertEquals(error, 0)
        self.assertEquals(1, index1)
        
        mass,error = instance.get_mass(1)
        self.assertEquals(error, 0)
        self.assertEquals(mass, 2.0)
                
        mass,error = instance.get_mass(2)
        self.assertEquals(error, 0)
        self.assertEquals(mass, 0.001)
        
    def test4(self):
        instance = MikkolaInterface()
        instance.initialize_code()
        instance.set_lightspeed(1e4)
        instance.commit_parameters()

        index1, error = instance.new_particle(mass = 1.0, radius = 0, x = 1.0, y = 2.0, z = 3.0, vx = 1.1, vy = 2.2, vz = 3.3)
        self.assertEquals(error, 0)
        self.assertEquals(index1, 1)
        
        instance.commit_particles()
        
        error = instance.cleanup_code()
        self.assertEquals(error, 0)
        
        mass,error = instance.get_mass(1)
        self.assertEquals(error, -1)
        instance.commit_parameters()
        index1, error = instance.new_particle(mass = 1.0, radius = 0, x = 1.0, y = 2.0, z = 3.0, vx = 1.1, vy = 2.2, vz = 3.3)
        self.assertEquals(error, 0)
        self.assertEquals(index1, 1)
        instance.commit_particles()
        error = instance.cleanup_code()
        self.assertEquals(error, 0)
        
        
        
# run with: 
# %>nosetests -v test_mikkola.py:TestMikkola.test1
class TestMikkola(TestWithMPI):
    def new_system_of_sun_and_earth(self):
        stars = datamodel.Stars(2)
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

        earth = stars[1]
        earth.mass = units.kg(5.9736e24)
        earth.radius = units.km(6371) 
        earth.position = units.km(numpy.array((149.5e6,0.0,0.0)))
        earth.velocity = units.ms(numpy.array((0.0,29800,0.0)))
        
        return stars

    def test1(self):
        convert_nbody=nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.yr/(2.0*pi))
        instance = Mikkola(convert_nbody)
        instance.parameters.timestep = 0.5 | units.day
        stars = self.new_system_of_sun_and_earth()
        instance.particles.add_particles(stars)
        Sun = stars[0]
        earth = stars[1]
        
        postion_at_start = earth.position.value_in(units.AU)[0]
        
#        instance.evolve_model(365.0 | units.day)

        instance.evolve_model(1.0 | units.yr)
        channel = instance.particles.new_channel_to(stars)
        channel.copy()
        self.assertAlmostRelativeEquals(instance.potential_energy * -0.5, instance.kinetic_energy, 3)
        self.assertAlmostRelativeEquals(instance.radiated_gravitational_energy, -6222456075.98| units.m**2 * units.kg * units.s**-2, 4)
        
        postion_after_full_rotation = earth.position.value_in(units.AU)[0]
       
        self.assertAlmostEqual(postion_at_start, instance.particles[1].position.value_in(units.AU)[0], 3)
        self.assertAlmostEqual(postion_at_start, postion_after_full_rotation, 3)
        
        instance.evolve_model(1.5 | units.yr)
        
        channel.copy()
        
        postion_after_half_a_rotation = earth.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        
        instance.evolve_model(1.75  | units.yr)
         
        channel.copy()
        
        postion_after_half_a_rotation = earth.position.value_in(units.AU)[1]
        
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        instance.cleanup_code()
        del instance

# run with: 
# %>nosetests -v test_mikkola.py:TestMikkola.test1
    def new_system_of_sun_and_mercury(self):
        stars = datamodel.Stars(2)
        sun = stars[0]
        sun.mass = units.MSun(1.0)
        sun.position = units.m(numpy.array((0.0,0.0,0.0)))
        sun.velocity = units.ms(numpy.array((0.0,0.0,0.0)))
        sun.radius = units.RSun(1.0)

        mercury = stars[1]
        mercury.mass = 3.3022e23 | units.kg
        mercury.radius = 0 | units.RSun
        mercury.position = [0.387098, 0, 0] | units.AU
        mercury.velocity = [0.0, 47.87, 0.0] | units.km/units.s
        return stars

    def test2(self):
        convert_nbody=nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.yr/(2.0*pi))
        instance = Mikkola(convert_nbody)
        instance.parameters.timestep = 1 | units.day
        stars = self.new_system_of_sun_and_mercury()
        instance.particles.add_particles(stars)
        Sun = stars[0]
        mercury = stars[1]
        
        postion_at_start = mercury.position.value_in(units.AU)[0]
        
        period_mercury = 87.9691 | units.day
        instance.evolve_model(period_mercury)
        channel = instance.particles.new_channel_to(stars)
        channel.copy()
        
        postion_after_full_rotation = mercury.position.value_in(units.AU)[0]
       
        self.assertAlmostEqual(postion_at_start, postion_after_full_rotation, 4)
        
        instance.evolve_model(1.5 * period_mercury)
        
        channel.copy()
        
        postion_after_half_a_rotation = mercury.position.value_in(units.AU)[0]
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        
        instance.evolve_model(1.75 *period_mercury)
         
        channel.copy()
        
        postion_after_half_a_rotation = mercury.position.value_in(units.AU)[1]
        
        self.assertAlmostEqual(-postion_at_start, postion_after_half_a_rotation, 3)
        instance.cleanup_code()
        del instance

# run with: 
# %>nosetests -v test_mikkola.py:TestMikkola.test3
    def new_system_of_Hulse_Taylor_pulsar(self):
        stars = datamodel.Stars(2)
        Hulse = stars[0]
        Hulse.mass = 1.441 | units.MSun
        Hulse.radius = 1.4e-5 | units.RSun
        Hulse.position = [-1576800.0, 0, 0] | units.km
#        Hulse.velocity = [0.0, -55.0, 0.0] | units.km/units.s
#        Hulse.position = [-1546691.3719943422, 0, 0] | units.km
        Hulse.velocity = [0.0, -110.0, 0.0] | units.km/units.s

        Taylor = stars[1]
        Taylor.mass = 1.387 | units.MSun
        Taylor.radius = 1.4e-5 | units.RSun
        Taylor.position = [1606908.6280056578, 0, 0] | units.km
        Taylor.velocity = [0.0, 114.28262436914201, 0.0] | units.km/units.s
#        Taylor.position = [1576800.0, 0, 0] | units.km
#        Taylor.velocity = [0.0, 55.0, 0.0] | units.km/units.s
        return stars

    def test3(self):
        convert_nbody=nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.yr/(2.0*pi))
        instance = Mikkola(convert_nbody)
        stars = self.new_system_of_Hulse_Taylor_pulsar()
        instance.particles.add_particles(stars)
        Hulse = stars[0]
        Taylor = stars[1]
        
        postion_at_start = Taylor.position.value_in(units.AU)[0]
        
        #orbital period
        #see http://www.johnstonsarchive.net/relativity/binpulsar.html
        period_HTpulsar = 7.75 | units.hour
        #period_HTpulsar = 77.51939106 | units.hour
        # period for abseidal motion
        #        period_HTpulsar = 85.0 | units.yr #4.2degrees/year
        #period_HTpulsar = 1.0 | units.yr 
        instance.evolve_model(period_HTpulsar)
        channel = instance.particles.new_channel_to(stars)
        channel.copy()

        postion_after_full_rotation = Taylor.position.value_in(units.AU)[0]
       
        print "Time=", instance.model_time, period_HTpulsar
        self.assertAlmostEqual(postion_at_start, postion_after_full_rotation, 3)
        
        instance.stop()
        del instance
        
    def test4(self):
        convert_nbody=nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.yr/(2.0*pi))
        instance = Mikkola(convert_nbody)

        stars = self.new_system_of_Hulse_Taylor_pulsar()
        instance.particles.add_particles(stars)
        
        instance.commit_particles()
        self.assertEquals(len(instance.particles), 2)
        instance.cleanup_code()
        self.assertEquals(len(instance.particles), 0)
        
        instance.initialize_code()
        instance.particles.add_particles(stars)
        self.assertEquals(len(instance.particles), 2)
        
        instance.commit_particles()
        self.assertEquals(len(instance.particles), 2)
        instance.cleanup_code()
        self.assertEquals(len(instance.particles), 0)

    def test5(self):
        convert_nbody=nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.yr/(2.0*pi))
        
        instance = Mikkola(convert_nbody)
        
        stars = datamodel.Stars(3)
        stars[0].mass = 1.0 | units.MSun
        stars[0].radius = 1.0| units.RSun
        stars[0].position = [0.0, 0, 0] | units.km
        stars[0].velocity = [0.0,0.0,0.0] | units.km/units.s

        stars[1].mass = 1.0 | units.MSun
        stars[1].radius = 1.0| units.RSun
        stars[1].position = [0.1, 0, 0] | units.RSun
        stars[1].velocity = [-0.4, 0.0, 0.0] | units.km/units.s
        
        
        stars[2].mass = 0.01 | units.MSun
        stars[2].radius = 1.0| units.RSun
        stars[2].position = [6000, 0, 0] | units.RSun
        stars[2].velocity = [0.0, -10, 0.0] | units.km/units.s
        
        
        instance.particles.add_particles(stars)
        
        instance.evolve_model(0.25 | units.yr)
        
        self.assertEquals(instance.get_number_of_particles_added(), 1)
        self.assertRaises(Exception, instance.get_id_of_added_particle, [2])
        self.assertEquals(instance.get_id_of_added_particle(0), 4)
        self.assertAlmostRelativeEquals(instance.get_mass(4), 2 | units.MSun)
        self.assertAlmostRelativeEquals(instance.get_mass(1), 1 | units.MSun)
        self.assertAlmostRelativeEquals(instance.get_mass(2), 1 | units.MSun)
        instance.update_particle_set()
        print instance.particles[0].mass.as_quantity_in(units.MSun)
        #self.assertEquals(len(instance.particles), 4)

    def test6(self):
        convert_nbody=nbody_system.nbody_to_si(1.0|units.MSun, 1.0|units.yr/(2.0*pi))
        
        instance = Mikkola(convert_nbody)
        
        stars = datamodel.Stars(4)
        stars[0].mass = 1.0 | units.MSun
        stars[0].radius = 1.0| units.RSun
        stars[0].position = [0.0, 0, 0] | units.km
        stars[0].velocity = [0.0,0.0,0.0] | units.km/units.s

        stars[1].mass = 1.0 | units.MSun
        stars[1].radius = 1.0| units.RSun
        stars[1].position = [0.1, 0, 0] | units.RSun
        stars[1].velocity = [-0.4, 0.0, 0.0] | units.km/units.s
        
        stars[2].mass = 1.0 | units.MSun
        stars[2].radius = 1.0| units.RSun
        stars[2].position = [0.3, 0, 0] | units.RSun
        stars[2].velocity = [-0.4, 0.0, 0.0] | units.km/units.s
        
        
        stars[3].mass = 0.01 | units.MSun
        stars[3].radius = 1.0| units.RSun
        stars[3].position = [6000, 0, 0] | units.RSun
        stars[3].velocity = [0.0, -10, 0.0] | units.km/units.s
        
        
        instance.particles.add_particles(stars)
        
        instance.evolve_model(0.25 | units.yr)
        
        self.assertEquals(instance.get_number_of_particles_added(), 2)
        self.assertRaises(Exception, instance.get_id_of_added_particle, [2])
        self.assertEquals(instance.get_id_of_added_particle(0), 5)
        self.assertEquals(instance.get_id_of_added_particle(1), 6)
        self.assertAlmostRelativeEquals(instance.get_mass(6), 3 | units.MSun)
        self.assertAlmostRelativeEquals(instance.get_mass(5), 2 | units.MSun)
        instance.update_particle_set()
        print instance.particles[0].mass.as_quantity_in(units.MSun)
        #self.assertEquals(len(instance.particles), 4)
