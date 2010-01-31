import platform

from amuse.legacy.bse.interface import BSE

from amuse.support.data import core
from amuse.support.units import units

from legacy_support import TestWithMPI

class TestMPIInterface(TestWithMPI):
    
    class state(object):
        def __init__(self):
            self.type1 = 0.0
            self.type2 = 0.0
            self.initial_mass1 = 0.0
            self.initial_mass2 = 0.0
            self.mass1 = 0.0
            self.mass2 = 0.0
            self.radius1 = 0.0
            self.radius2 = 0.0
            self.luminosity1  = 0.0
            self.luminosity2  = 0.0
            self.core_mass1 = 0.0
            self.core_mass2 = 0.0
            self.core_radius1 = 0.0
            self.core_radius2 = 0.0
            self.envelope_mass1 = 0.0
            self.envelope_mass2 = 0.0
            self.envelope_radius1 = 0.0
            self.envelope_radius2 = 0.0
            self.spin1 = 0.0
            self.spin2 = 0.0
            self.epoch1 = 0.0
            self.epoch2 = 0.0
            self.t_ms1 = 0.0
            self.t_ms2 = 0.0
            self.bse_age = 0.0
            self.orbital_period = 0.0
            self.eccentricity = 0.0
        
    def test1(self):
        instance = BSE()
        metallicity = 0.02
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        alpha1 = 1.0
        CElambda = 0.5
        ceflag = 0
        tflag = 1
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
        sigma =  190.0
        beta = 1.0/8.0
        xi = 1.0
        acc2 = 3.0/2.0
        epsnov = 0.001
        eddfac = 1.0
        gamma = -1.0

        status = instance.initialize(metallicity,
            neta, bwind, hewind, alpha1, CElambda,
            ceflag, tflag, ifflag, wdflag, bhflag,
            nsflag, mxns, pts1, pts2, pts3,
            sigma,beta,xi,acc2,epsnov,eddfac,gamma)
        self.assertEqual(status,0)
        del instance
        
    def test2(self):
        instance = BSE()
        metallicity = 0.02
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        alpha1 = 1.0
        CElambda = 0.5
        ceflag = 0
        tflag = 1
        ifflag = 0
        wdflag =  1
        bhflag =  0 
        nsflag =  1
        mxns =  3.0
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
        sigma =  190.0
        beta = 1.0/8.0
        xi = 1.0
        acc2 = 3.0/2.0
        epsnov = 0.001
        eddfac = 1.0
        gamma = -1.0

        status = instance.initialize(metallicity,
            neta, bwind, hewind, alpha1, CElambda,
            ceflag, tflag, ifflag, wdflag, bhflag,
            nsflag, mxns, pts1, pts2, pts3,
            sigma,beta,xi,acc2,epsnov,eddfac,gamma)
        self.assertEqual(status,0)
        new_state = self.state()
        new_state.mass1 = 3.0
        new_state.mass2 = 1.0
        new_state.initial_mass1 = 3.0
        new_state.initial_mass2 = 1.0
        new_state.type1 = 1.0
        new_state.type2 = 1.0
        new_state.end_time = 1e-06
        new_state.orbital_period = 200.0
        new_state.eccentricity = 0.5
        result = instance.evolve(
            new_state.type1,new_state.type2,new_state.initial_mass1,new_state.initial_mass2,
            new_state.mass1, new_state.mass2, new_state.radius1, new_state.radius2, 
            new_state.luminosity1, new_state.luminosity2, new_state.core_mass1, 
            new_state.core_mass2, new_state.core_radius1, new_state.core_radius2,
            new_state.envelope_mass1, new_state.envelope_mass2, new_state.envelope_radius1,
            new_state.envelope_radius2, new_state.spin1, new_state.spin2, new_state.epoch1,
            new_state.epoch2, new_state.t_ms1, new_state.t_ms2, new_state.bse_age,
            new_state.end_time, new_state.orbital_period, new_state.eccentricity
        )
        print result
        updated_state = self.state()
        (updated_state.type1,updated_state.type2,updated_state.initial_mass1,updated_state.initial_mass2,
            updated_state.mass1, updated_state.mass2, updated_state.radius1, updated_state.radius2, 
            updated_state.luminosity1, updated_state.luminosity2, updated_state.core_mass1,
            updated_state.core_mass2, updated_state.core_radius1, updated_state.core_radius2,
            updated_state.envelope_mass1,updated_state.envelope_mass2,updated_state.envelope_radius1,
            updated_state.envelope_radius2, updated_state.spin1, updated_state.spin2,
            updated_state.epoch1, updated_state.epoch2, updated_state.t_ms1, updated_state.t_ms2,
            updated_state.bse_age, updated_state.end_time,
            updated_state.orbital_period, updated_state.eccentricity) = result        
         
        expected = {
            'radius2' : '0x1.c6c8a1c793bcep-1',
            'luminosity2' : '0x1.653b1b2d0333bp-1',
            'core_mass2' : '0x0.0p+0',
            'bse_age' : '0x1.0c6f7a0b5ed8dp-20',
            'end_time' : '0x1.0c6f7a0b5ed8dp-20',
            'envelope_mass2' : '0x1.0d6fc100ab510p-5',
            'mass2' : '0x1.0000000000000p+0',
            'initial_mass2' : '0x1.0000000000000p+0',
            'envelope_radius2' : '0x1.db27631ba0e5ap-3',
            'core_radius2' : '0x0.0p+0',
            'epoch2' : '0x0.0p+0',
            't_ms2' : '0x1.57d90abe54643p+13',
            'spin2' : '0x1.07413b0522aebp+10',
        };    

        expected_32bit = {
            'core_radius2' : '0x0.0p+0',
            'initial_mass2' : '0x1.0000000000000p+0',
            'envelope_mass2' : '0x1.0d6fc100ab50fp-5',
            'radius2' : '0x1.c6c8a1c793bd0p-1',
            't_ms2' : '0x1.57d90abe54642p+13',
            'spin2' : '0x1.07413b0522d61p+10',
            'luminosity2' : '0x1.653b1b2d0333ap-1',
            'end_time' : '0x1.0c6f7a0b5ed8dp-20',
            'envelope_radius2' : '0x1.db27631ba0e5cp-3',
            'bse_age' : '0x1.0c6f7a0b5ed8dp-20',
            'epoch2' : '0x0.0p+0',
            'mass2' : '0x1.0000000000000p+0',
            'core_mass2' : '0x0.0p+0',
        };

        architecture, linkage_format = platform.architecture()
        if architecture == '32bit' and platform.system() != 'Darwin':
            for x in expected:
                #print "'%s' : '%s'," % (x, getattr(updated_state, x).hex())
                self.assertEqual(float.fromhex(expected_32bit[x]),getattr(updated_state, x))
        else:
            for x in expected:
                #print "'%s' : '%s'," % (x, getattr(updated_state, x).hex())
                self.assertEqual(float.fromhex(expected[x]),getattr(updated_state, x))
            
        self.assertEquals(updated_state.end_time, 1e-06)
        dt = instance.get_time_step(updated_state.type1, updated_state.type2,
            updated_state.initial_mass1, updated_state.initial_mass2, updated_state.mass1,
            updated_state.mass2, updated_state.t_ms1, updated_state.t_ms2,
            updated_state.epoch1, updated_state.epoch2, updated_state.bse_age)
        self.assertAlmostEqual(dt, 18.8768, 3)
        del instance
     
    def xtest3(self):
        instance = BSE()
        instance.initialize_module_with_default_parameters()  
        types = [1,1,1]
        masses = [10,5,4]
        radii = [5.0, 2.0, 1.0]
        luminosity = core_mass = core_radius =  envelope_mass =\
        envelope_radius =  spin = epoch = t_ms = [0.0,0.0,0.0]
        sse_age = age = [1e-6, 1e-06, 1e-6]
        result = instance.evolve(
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
        del instance
        
    def xtest4(self):
        instance = BSE()
        instance.initialize_module_with_default_parameters()  
        types = [1 for x in range(1,4000)]
        masses = [1.0 + ((x / 4000.0) * 10.0) for x in range(1,4000)]
        radii = [1.0 for x in range(1,4000)]
        luminosity = core_mass = core_radius =  envelope_mass =\
        envelope_radius =  spin = epoch =\
        t_ms = [0.0 for x in range(1,4000)]
        
        sse_age = age = [1e-06 for x in range(1,4000)]
        result = instance.evolve(
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
        del instance

    def xtest5(self):
        instance = BSE()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5 | units.none)
        myvalue = 0.7 | units.none
        instance.parameters.reimers_mass_loss_coefficient = myvalue
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, myvalue)
        instance.initialize_module_with_current_parameters()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, myvalue)
        del instance
        
        instance = BSE()
        myvalue = 0.7 | units.none
        instance.parameters.reimers_mass_loss_coefficient = myvalue
        instance.initialize_module_with_default_parameters()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5 | units.none)
        del instance
        
class TestSSE(TestWithMPI):
    
    def xtest1(self):
        instance = BSE()
        instance.initialize_module_with_default_parameters() 
        stars =  core.Stars(1)
        
        star = stars[0]
        star.mass = 5 | units.MSun
        star.radius = 0.0 | units.RSun
        
        instance.particles.add_particles(stars)
        from_sse_to_model = instance.particles.new_channel_to(stars)
        from_sse_to_model.copy()
        
        previous_type = star.type
        results = []
        t0 = 0 | units.Myr
        current_time = t0
        
        while current_time < (125 | units.Myr):
            instance.update_time_steps()
            
            current_time = current_time + instance.particles[0].time_step
            
            instance.evolve_model(current_time)

            from_sse_to_model.copy()
            
            if not star.type == previous_type:
                results.append((star.age, star.mass, star.type))
                previous_type = star.type
                
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
        
        del instance
            
    def xtest2(self):
        instance = BSE()
        instance.initialize_module_with_default_parameters() 
        stars =  core.Stars(1)
        
        star = stars[0]
        star.mass = 5 | units.MSun
        star.radius = 0.0 | units.RSun
        
        instance.particles.add_particles(stars)
        instance.evolve_model(120.1 | units.Myr)
                
        self.assertAlmostEqual(instance.particles[0].mass.value_in(units.MSun), 4.932, 3)
         
        del instance
        
    
    def xtest3(self):
        instance = BSE()
        instance.initialize_module_with_default_parameters() 
        stars =  core.Stars(1)
        
        star = stars[0]
        star.mass = 5 | units.MSun
        star.radius = 0.0 | units.RSun
        
        stars.synchronize_to(instance.particles)
        
        channel = instance.particles.new_channel_to(stars)
        channel.copy_attributes(instance.particles._get_attributes())   
        
        previous_type = instance.particles.type
        results = []
        
        instance.evolve_model(121.5 | units.Myr)
        
        channel.copy_attributes(instance.particles._get_attributes())   
        
        self.assertAlmostEqual(star.mass.value_in(units.MSun), 0.997, 3)
         
        del instance
        
    
    def xtest5(self):
        instance = BSE()
        instance.initialize_module_with_default_parameters() 
        stars =  core.Stars(1)
        
        star = stars[0]
        star.mass = 35 | units.MSun
        star.radius = 0.0 | units.RSun
        
        stars.synchronize_to(instance.particles)
        
        channel = instance.particles.new_channel_to(stars)
        channel.copy_attributes(instance.particles._get_attributes())   
        
        previous_type = star.type
        results = []
        
        dt = 1 | units.Myr
        t = 0 | units.Myr
        while t < 30 | units.Myr:
            t += dt
            instance.evolve_model(t)
                
        self.assertTrue(instance.particles[0].mass.value_in(units.MSun) < 10.6)
         
        del instance


    def xtest6(self):
#       Test whether a set of stars evolve synchronously
#       Create an array of stars with a range in stellar mass
        masses = [.5, 1., 2., 5., 10., 30.] | units.MSun
        number_of_stars = len(masses)
        stars = core.Stars(number_of_stars)
        for i, star in enumerate(stars):
            star.mass = masses[i]
            star.radius = 0.0 | units.RSun

#       Initialize stellar evolution code
        instance = BSE()
        instance.initialize_module_with_default_parameters() 
        instance.setup_particles(stars)
#       Let the code perform initialization actions after all particles have been created. 
#       Called before the first evolve call and after the last new_particle call.
        instance.initialize_stars()
        
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
            self.assertTrue(stars[i].age.value_in(units.Myr) > 125)
            self.assertTrue(stars[i].mass <= masses[i])
            self.assertEquals(str(stars[i].type), end_types[i])

        del instance

