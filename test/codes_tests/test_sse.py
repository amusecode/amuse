import os

from amuse.community.sse.interface import SSEInterface, SSE

from amuse.support.data.core import Particles
from amuse.support.units import units

from amuse.test.amusetest import TestWithMPI

from amuse.support import io

class TestMPIInterface(TestWithMPI):
    
    class state(object):
        def __init__(self):
            self.stellar_type = 0.0
            self.zams_mass = 0.0
            self.mass = 0.0
            self.radius = 0.0
            self.luminosity  = 0.0
            self.core_mass = 0.0
            self.core_radius = 0.0
            self.envelope_mass = 0.0
            self.envelope_radius = 0.0
            self.spin = 0.0
            self.epoch = 0.0
            self.t_ms = 0.0
            self.sse_age = 0.0
        
    def initialize_module_with_default_parameters(self, sse):
        metallicity = 0.02
        
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        sigma =  190.0
        
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
    
    
        status = sse.initialize(metallicity,
            neta, bwind, hewind, sigma,
            ifflag, wdflag, bhflag, nsflag, mxns,
            pts1, pts2, pts3)
        
    def test1(self):
        sse = SSEInterface()
        
        metallicity = 0.02
        
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        sigma =  190.0
        
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
    
    
        status = sse.initialize(metallicity,
            neta, bwind, hewind, sigma,
            ifflag, wdflag, bhflag, nsflag, mxns,
            pts1, pts2, pts3)
            
        self.assertEqual(status,0)
        
        sse.stop()
        
    def test2(self):
        sse = SSEInterface()
        
        metallicity = 0.02
        
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        sigma =  190.0
        
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
    
    
        status = sse.initialize(metallicity,
            neta, bwind, hewind, sigma,
            ifflag, wdflag, bhflag, nsflag, mxns,
            pts1, pts2, pts3)
        self.assertEqual(status,0)
        new_state = self.state()
        new_state.mass = 1.0
        new_state.zams_mass = 1.0
        new_state.stellar_type = 1.0
        new_state.age = 1e-06
        result = sse.evolve(
            new_state.stellar_type, 
            new_state.zams_mass, new_state.mass, new_state.radius, 
            new_state.luminosity, new_state.core_mass, new_state.core_radius,
            new_state.envelope_mass, new_state.envelope_radius, new_state.spin,
            new_state.epoch, new_state.t_ms, new_state.sse_age, new_state.age
        )
        updated_state = self.state()
        (updated_state.stellar_type,updated_state.zams_mass, updated_state.mass, updated_state.radius, 
            updated_state.luminosity, updated_state.core_mass, updated_state.core_radius,
            updated_state.envelope_mass, updated_state.envelope_radius, updated_state.spin,
            updated_state.epoch, updated_state.t_ms, updated_state.sse_age, updated_state.age) = result
        attributes = ('stellar_type', 'zams_mass', 'mass', 'radius', 'luminosity', 'core_mass', 'core_radius',
            'envelope_mass', 'envelope_radius', 'spin', 'epoch', 't_ms', 'sse_age', 'age')
        
         
        expected = {
            'zams_mass': '0x1.0000000000000p+0',
            'mass': '0x1.0000000000000p+0',
            'radius': '0x1.c6c8a1c793bcep-1',
            'luminosity': '0x1.653b1b2d0333bp-1',
            'core_mass': '0x0.0p+0',
            'core_radius': '0x0.0p+0',
            'envelope_mass': '0x1.0d6fc100ab510p-5',
            'envelope_radius': '0x1.db27631ba0e5ap-3',
            'spin': '0x1.07413b0522d63p+10',
            'epoch': '0x0.0p+0',
            't_ms': '0x1.57d90abe54643p+13',
            'sse_age': '0x1.0c6f7a0b5ed8dp-20',
            'age': '0x1.0c6f7a0b5ed8dp-20',
        };    


        for x in expected:
            self.assertAlmostRelativeEqual(float.fromhex(expected[x]),getattr(updated_state, x))
            
        self.assertEquals(updated_state.age, 1e-06)
        dt = sse.get_time_step(updated_state.stellar_type,
            updated_state.zams_mass, 
            updated_state.age, 
            updated_state.mass, 
            updated_state.t_ms, 
            updated_state.epoch)
        self.assertAlmostEqual(dt, 550.1565, 2)
        sse.stop()
     
    def test3(self):
        sse = SSEInterface()
        self.initialize_module_with_default_parameters(sse)  
        types = [1,1,1]
        masses = [10,5,4]
        radii = [5.0, 2.0, 1.0]
        luminosity = core_mass = core_radius =  envelope_mass =\
        envelope_radius =  spin = epoch = t_ms = [0.0,0.0,0.0]
        sse_age = age = [1e-6, 1e-06, 1e-6]
        result = sse.evolve(
            types, 
            masses, 
            masses, 
            radii, 
            luminosity, 
            core_mass, 
            core_radius,
            envelope_mass,
            envelope_radius, 
            spin,
            epoch, 
            t_ms, 
            sse_age, 
            age
        )
        self.assertEquals(result['mass'][0], 10)
        self.assertEquals(result['mass'][1], 5)
        self.assertAlmostEqual(result['mass'][2], 4.0, 2)
        sse.stop()
        
    def test4(self):
        sse = SSEInterface()
        self.initialize_module_with_default_parameters(sse) 
        types = [1 for x in range(1,4000)]
        masses = [1.0 + ((x / 4000.0) * 10.0) for x in range(1,4000)]
        radii = [1.0 for x in range(1,4000)]
        luminosity = core_mass = core_radius =  envelope_mass =\
        envelope_radius =  spin = epoch =\
        t_ms = [0.0 for x in range(1,4000)]
        
        sse_age = age = [1e-06 for x in range(1,4000)]
        result = sse.evolve(
            types, 
            masses, 
            masses, 
            radii, 
            luminosity, 
            core_mass, 
            core_radius,
            envelope_mass,
            envelope_radius, 
            spin,
            epoch, 
            t_ms, 
            sse_age, 
            age
        )
        self.assertEquals(len(result['mass']), 3999)
        sse.stop()

        
class TestSSE(TestWithMPI):
    
    def test1(self):
        sse = SSE()
        sse.commit_parameters() 
        stars = Particles(1)
        star = stars[0]
        star.mass = 5 | units.MSun
        star.radius = 0.0 | units.RSun
        
        sse.particles.add_particles(stars)
        from_sse_to_model = sse.particles.new_channel_to(stars)
        from_sse_to_model.copy()
        
        previous_type = star.stellar_type
        results = []
        t0 = 0 | units.Myr
        current_time = t0
        
        while current_time < (125 | units.Myr):
            sse.update_time_steps()
            
            current_time = current_time + sse.particles[0].time_step
            
            sse.evolve_model(current_time)

            from_sse_to_model.copy()
            
            if not star.stellar_type == previous_type:
                results.append((star.age, star.mass, star.stellar_type))
                previous_type = star.stellar_type
                
        self.assertEqual(len(results), 6)
        
        times = ( 
            104.0 | units.Myr, 
            104.4 | units.Myr, 
            104.7 | units.Myr, 
            120.1 | units.Myr,
            120.9 | units.Myr,
            121.5 | units.Myr
        )
        for result, expected in zip(results, times):
            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 1)
            
        masses = ( 
            5.000 | units.MSun, 
            5.000 | units.MSun, 
            4.998 | units.MSun, 
            4.932 | units.MSun,
            4.895 | units.MSun,
            0.997 | units.MSun
        )
        for result, expected in zip(results, masses):
            self.assertAlmostEqual(result[1].value_in(units.MSun), expected.value_in(units.MSun), 3)
         
        types = (
            "Hertzsprung Gap",
            "First Giant Branch",
            "Core Helium Burning",
            "First Asymptotic Giant Branch",
            "Second Asymptotic Giant Branch",
            "Carbon/Oxygen White Dwarf",
        )
        
        for result, expected in zip(results, types):
            self.assertEquals(str(result[2]), expected)
        
        sse.stop()
            
    def test2(self):
        sse = SSE()
        sse.commit_parameters() 
        stars = Particles(1)
        
        star = stars[0]
        star.mass = 5 | units.MSun
        star.radius = 0.0 | units.RSun
        
        sse.particles.add_particles(stars)
        sse.evolve_model(120.1 | units.Myr)
                
        self.assertAlmostEqual(sse.particles[0].mass.value_in(units.MSun), 4.932, 3)
        self.assertAlmostEqual(sse.particles[0].temperature.value_in(units.K), 4221., 0)
         
        sse.stop()
        
    
    def test3(self):
        sse = SSE()
        sse.commit_parameters() 
        stars = Particles(1)
        
        star = stars[0]
        star.mass = 5 | units.MSun
        star.radius = 0.0 | units.RSun
        
        stars.synchronize_to(sse.particles)
        
        channel = sse.particles.new_channel_to(stars)
        channel.copy_attributes(sse.particles.get_attribute_names_defined_in_store())   
        
        previous_type = sse.particles.stellar_type
        results = []
        
        sse.evolve_model(121.5 | units.Myr)
        
        channel.copy_attributes(sse.particles.get_attribute_names_defined_in_store())   
        
        self.assertAlmostEqual(star.mass.value_in(units.MSun), 0.997, 3)
         
        sse.stop()
        
    
    def test5(self):
        sse = SSE()
        sse.commit_parameters() 
        stars = Particles(1)
        
        star = stars[0]
        star.mass = 35 | units.MSun
        star.radius = 0.0 | units.RSun
        
        stars.synchronize_to(sse.particles)
        
        channel = sse.particles.new_channel_to(stars)
        channel.copy_attributes(sse.particles.get_attribute_names_defined_in_store())   
        
        previous_type = star.stellar_type
        results = []
        
        dt = 1 | units.Myr
        t = 0 | units.Myr
        while t < 30 | units.Myr:
            t += dt
            sse.evolve_model(t)
                
        self.assertTrue(sse.particles[0].mass.value_in(units.MSun) < 10.6)
         
        sse.stop()


    def test6(self):
        print "Test whether a set of stars evolves synchronously..."
#       Create an array of stars with a range in stellar mass
        masses = [.5, 1., 2., 5., 10., 30.] | units.MSun
        number_of_stars = len(masses)
        stars = Particles(number_of_stars)
        stars.mass = masses

#       Initialize stellar evolution code
        instance = SSE()
        instance.commit_parameters() 
        instance.particles.add_particles(stars)
        instance.commit_particles()
        
        from_code_to_model = instance.particles.new_channel_to(stars)
        from_code_to_model.copy()
        
        instance.evolve_model(end_time = 125 | units.Myr)
        from_code_to_model.copy()
                
        end_types = (
            "deeply or fully convective low mass MS star",
            "Main Sequence star",
            "Main Sequence star",
            "Carbon/Oxygen White Dwarf",
            "Neutron Star",
            "Black Hole",
        )
        for i in range(number_of_stars):
            self.assertAlmostEquals(stars[i].age, 125.0 | units.Myr)
            self.assertTrue(stars[i].mass <= masses[i])
            self.assertEquals(str(stars[i].stellar_type), end_types[i])
        instance.stop()
    
    def test7(self):
        print "Test: evolve particles one at a time."
        print "Used to be problematic, since initial_mass of idle particle is set to zero."
        stars = Particles(2)
        stars.mass = 1.0 | units.MSun
        for star in stars:
            print star
            stellar_evolution = SSE()
            stellar_evolution.commit_parameters()
            stellar_evolution.particles.add_particles(star.as_set())
            stellar_evolution.commit_particles()
            from_stellar_evolution_to_model = stellar_evolution.particles.new_channel_to(star.as_set())
            stellar_evolution.evolve_model()
            from_stellar_evolution_to_model.copy()
            stellar_evolution.stop()
        self.assertEquals(stars[0].initial_mass, stars[1].initial_mass)
        self.assertEquals(stars[0].luminosity, stars[1].luminosity)
        self.assertEquals(stars[0].age, stars[1].age)
        print "Solved: SSE_muse_interface.f sets initial_mass to mass when necessary."
    
    def test8(self):
        instance = SSE()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5 | units.none)
        myvalue = 0.7 | units.none
        instance.parameters.reimers_mass_loss_coefficient = myvalue
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, myvalue)
        instance.commit_parameters()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, myvalue)
        instance.stop()
        
        instance = SSE()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5 | units.none)
        myvalue = 0.7 | units.none
        instance.parameters.reimers_mass_loss_coefficient = myvalue
        instance.parameters.set_defaults()
        instance.commit_parameters()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5 | units.none)
        instance.stop()
        
    def test9(self):
        print "Test: large number of particles"
        stellar_evolution = SSE(max_message_length=500)
        stellar_evolution.commit_parameters()
        number_of_particles = 10000
        print "Has been tested with up to a million particles!"
        print "Now using ", number_of_particles, "particles only, for speed."
        stars = Particles(number_of_particles)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.particles.add_particles(stars)
        self.assertEqual(len(stellar_evolution.particles), number_of_particles)
        stellar_evolution.stop()
    

    def test10(self):
        stellar_evolution = SSE()
        stellar_evolution.commit_parameters()
        stars = Particles(10)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.particles.add_particles(stars)
        self.assertEquals(stellar_evolution.particles._factory_for_new_collection(), Particles)
        if os.path.exists('test.h5'):
            os.remove('test.h5')
            
        io.write_set_to_file(stellar_evolution.particles, 'test.h5', 'hdf5')
        stored_stars = io.read_set_from_file('test.h5', 'hdf5')
        self.assertEquals(len(stars), len(stored_stars))
    
        self.assertAlmostRelativeEquals(stars.mass, stored_stars.mass)
    
    def test11(self):
        print "Test evolve_model optional arguments: end_time and keep_synchronous"
        stars = Particles(3)
        stars.mass = [1.0, 2.0, 3.0] | units.MSun
        instance = SSE()
        instance.commit_parameters()
        instance.particles.add_particles(stars)
        
        self.assertAlmostEqual(instance.particles.age, [0.0, 0.0, 0.0] | units.yr)
        self.assertAlmostEqual(instance.particles.time_step, [550.1565, 58.2081, 18.8768] | units.Myr, 3)
        
        print "evolve_model without arguments: use shared timestep = min(particles.time_step)"
        instance.evolve_model()
        self.assertAlmostEqual(instance.particles.age, [18.8768, 18.8768, 18.8768] | units.Myr, 3)
        self.assertAlmostEqual(instance.particles.time_step, [550.1565, 58.2081, 18.8768] | units.Myr, 3)
        self.assertAlmostEqual(instance.model_time, 18.8768 | units.Myr, 3)
        
        print "evolve_model with end_time: take timesteps, until end_time is reached exactly"
        instance.evolve_model(100 | units.Myr)
        self.assertAlmostEqual(instance.particles.age, [100.0, 100.0, 100.0] | units.Myr, 3)
        self.assertAlmostEqual(instance.particles.time_step, [550.1565, 58.2081, 18.8768] | units.Myr, 3)
        self.assertAlmostEqual(instance.model_time, 100.0 | units.Myr, 3)
        
        print "evolve_model with keep_synchronous: use non-shared timestep, particle ages will typically diverge"
        instance.evolve_model(keep_synchronous = False)
        self.assertAlmostEqual(instance.particles.age, (100 | units.Myr) + ([550.1565, 58.2081, 18.8768] | units.Myr), 3)
        self.assertAlmostEqual(instance.particles.time_step, [550.1565, 58.2081, 18.8768] | units.Myr, 3)
        self.assertAlmostEqual(instance.model_time, 100.0 | units.Myr, 3) # Unchanged!
        instance.stop()
        
    def test12(self):
        print "Testing adding and removing particles from stellar evolution code..."
        
        particles = Particles(3)
        particles.mass = 1.0 | units.MSun
        
        instance = SSE()
        instance.initialize_code()
        instance.commit_parameters()
        self.assertEquals(len(instance.particles), 0) # before creation
        instance.particles.add_particles(particles[:-1])
        instance.commit_particles()
        instance.evolve_model(1.0 | units.Myr)
        self.assertEquals(len(instance.particles), 2) # before remove
        self.assertAlmostEqual(instance.particles.age, 1.0 | units.Myr)
        
        instance.particles.remove_particle(particles[0])
        self.assertEquals(len(instance.particles), 1)
        instance.evolve_model(2.0 | units.Myr)
        self.assertAlmostEqual(instance.particles[0].age, 2.0 | units.Myr)
        
        instance.particles.add_particles(particles[::2])
        self.assertEquals(len(instance.particles), 3) # it's back...
        self.assertAlmostEqual(instance.particles[0].age, 2.0 | units.Myr)
        self.assertAlmostEqual(instance.particles[1].age, 0.0 | units.Myr)
        self.assertAlmostEqual(instance.particles[2].age, 0.0 | units.Myr) # ... and rejuvenated.
        
        instance.evolve_model(3.0 | units.Myr) # The young stars keep their age offset from the old star
        self.assertAlmostEqual(instance.particles.age, [3.0, 1.0, 1.0] | units.Myr)
        instance.evolve_model(4.0 | units.Myr)
        self.assertAlmostEqual(instance.particles.age, [4.0, 2.0, 2.0] | units.Myr)
        instance.stop()
    
    def test13(self):
        print "Testing SSE states"
        stars = Particles(1)
        stars.mass = 1.0 | units.MSun
        
        print "First do everything manually:",
        instance = SSE()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.initialize_code()
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.commit_parameters()
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.cleanup_code()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        instance.stop()
        print "ok"

        print "initialize_code(), commit_parameters(), " \
            "and cleanup_code() should be called automatically:",
        instance = SSE()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.reimers_mass_loss_coefficient = 0.5
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particles(stars)
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'END')
        print "ok"
    
