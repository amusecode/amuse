from amuse.legacy.bse.interface import BSE

from amuse.support.data import core
from amuse.support.units import units
from amuse.test.amusetest import TestWithMPI

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
        print "Test initialization..."
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
        idum = 29769
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
            nsflag, mxns, idum, pts1, pts2, pts3,
            sigma,beta,xi,acc2,epsnov,eddfac,gamma)
        self.assertEqual(status,0)
        del instance
        
    def test2(self):
        print "Test basic operations (legacy functions evolve & get_time_step)..."
        instance = BSE()
        instance.initialize_module_with_default_parameters()
        
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
        
        result = instance.legacy_interface.evolve(
            new_state.type1,new_state.type2,new_state.initial_mass1,new_state.initial_mass2,
            new_state.mass1, new_state.mass2, new_state.radius1, new_state.radius2, 
            new_state.luminosity1, new_state.luminosity2, new_state.core_mass1, 
            new_state.core_mass2, new_state.core_radius1, new_state.core_radius2,
            new_state.envelope_mass1, new_state.envelope_mass2, new_state.envelope_radius1,
            new_state.envelope_radius2, new_state.spin1, new_state.spin2, new_state.epoch1,
            new_state.epoch2, new_state.t_ms1, new_state.t_ms2, new_state.bse_age,
            new_state.orbital_period, new_state.eccentricity, new_state.end_time
        )
        
        updated_state = self.state()
        (updated_state.type1,updated_state.type2,updated_state.initial_mass1,updated_state.initial_mass2,
            updated_state.mass1, updated_state.mass2, updated_state.radius1, updated_state.radius2, 
            updated_state.luminosity1, updated_state.luminosity2, updated_state.core_mass1,
            updated_state.core_mass2, updated_state.core_radius1, updated_state.core_radius2,
            updated_state.envelope_mass1,updated_state.envelope_mass2,updated_state.envelope_radius1,
            updated_state.envelope_radius2, updated_state.spin1, updated_state.spin2,
            updated_state.epoch1, updated_state.epoch2, updated_state.t_ms1, updated_state.t_ms2,
            updated_state.bse_age, updated_state.orbital_period,
            updated_state.eccentricity, updated_state.end_time) = result        
         
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

        for x in expected:
            #print "'%s' : '%s'," % (x, getattr(updated_state, x).hex())
            self.assertAlmostRelativeEqual(float.fromhex(expected[x]),getattr(updated_state, x))
            
        self.assertEquals(updated_state.end_time, 1e-06)
        dt = instance.legacy_interface.get_time_step(updated_state.type1, updated_state.type2,
            updated_state.initial_mass1, updated_state.initial_mass2, updated_state.mass1,
            updated_state.mass2, updated_state.t_ms1, updated_state.t_ms2,
            updated_state.epoch1, updated_state.epoch2, updated_state.bse_age)
        self.assertAlmostEqual(dt, 18.8768, 3)
        del instance
     
    def test3(self):
        print "Test whether the interface can handle arrays..."
        instance = BSE()
        instance.initialize_module_with_default_parameters()
        masses1 = [10.0,5.0,4.0]
        masses2 = [1.0,1.0,1.0]
        types1 = types2 = [1,1,1]
        orbital_periods = [100.0,200.0,300.0]
        eccentricities = [0.5,0.6,0.7]

        radii1 = luminosity1 = core_mass1 = core_radius1 =  envelope_mass1 =\
        envelope_radius1 = spin1 = epoch1 = t_ms1 = [0.0,0.0,0.0]
        radii2 = luminosity2 = core_mass2 = core_radius2 =  envelope_mass2 =\
        envelope_radius2 = spin2 = epoch2 = t_ms2 = [0.0,0.0,0.0]
        init_mass1 = masses1
        init_mass2 = masses2
        bse_age = [0.0,0.0,0.0]
        end_time = [10.0, 10.0, 10.0]
        result = instance.legacy_interface.evolve(
            types1, types2, init_mass1, init_mass2,
            masses1, masses2, radii1, radii2,
            luminosity1, luminosity2, core_mass1, core_mass2,
            core_radius1, core_radius2, envelope_mass1, envelope_mass2,
            envelope_radius1, envelope_radius2, spin1, spin2,
            epoch1, epoch2, t_ms1, t_ms2, 
            bse_age, orbital_periods, eccentricities, end_time
        )
        self.assertAlmostEqual(result['mass1'][0], 9.977, 2)
        self.assertAlmostEqual(result['mass1'][1], 5.0, 2)
        self.assertAlmostEqual(result['mass1'][2], 4.0, 2)
        del instance
        
    def test4(self):
        print "Test large number of particles..."
        number_of_particles = 2000
        instance = BSE()
        instance.initialize_module_with_default_parameters()  
        masses1 = [1.0 + ((x / 1.0*number_of_particles) * 10.0) for x in range(1,number_of_particles+1)]
        masses2 = [2.0 + ((x / 1.0*number_of_particles) * 5.0) for x in range(1,number_of_particles+1)]
        orbital_periods = [100.0 + ((x / 1.0*number_of_particles) * 900.0) for x in range(1,number_of_particles+1)]
        eccentricities = [0.5 + ((x / 1.0*number_of_particles) * 0.4) for x in range(1,number_of_particles+1)]
        
        types1 = types2 = [1 for x in range(1,number_of_particles+1)]
        radii1 = luminosity1 = core_mass1 = core_radius1 =  envelope_mass1 =\
        envelope_radius1 =  spin1 = epoch1 = t_ms1 =\
        radii2 = luminosity2 = core_mass2 = core_radius2 =  envelope_mass2 =\
        envelope_radius2 =  spin2 = epoch2 = t_ms2 =\
        bse_age = [0.0 for x in range(1,number_of_particles+1)]
        end_time = [1.0 for x in range(1,number_of_particles+1)]
        init_mass1 = masses1
        init_mass2 = masses2
        
        result = instance.legacy_interface.evolve(
            types1, types2, init_mass1, init_mass2,
            masses1, masses2, radii1, radii2,
            luminosity1, luminosity2, core_mass1, core_mass2,
            core_radius1, core_radius2, envelope_mass1, envelope_mass2,
            envelope_radius1, envelope_radius2, spin1, spin2,
            epoch1, epoch2, t_ms1, t_ms2, 
            bse_age, orbital_periods, eccentricities, end_time
        )
        self.assertEquals(len(result['mass1']), number_of_particles)
        del instance

        
class TestBSE(TestWithMPI):
    
    def test1(self):
        print "Testing evolution of a close binary system..."
        instance = BSE()
        instance.parameters.metallicity = 0.001 | units.none
        instance.parameters.common_envelope_efficiency = 3.0 | units.none
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0 | units.none
        instance.initialize_module_with_current_parameters()
        stars =  core.Stars(1)
        
        binary = stars[0]
        binary.mass1 = 3.0 | units.MSun
        binary.mass2 = 0.3 | units.MSun
        binary.orbital_period = 200.0 | units.day
        binary.eccentricity = 0.5 | units.none
        
        instance.particles.add_particles(stars)
        from_bse_to_model = instance.particles.new_channel_to(stars)
        from_bse_to_model.copy()
        
        previous_type1 = binary.type1
        results = []
        current_time = 0 | units.Myr
        
        while current_time < (480 | units.Myr):
            instance.update_time_steps()
            # The next line appears a bit weird, but saves time for this simple test.
            current_time = current_time + max(5.0*instance.particles[0].time_step, 0.3 | units.Myr)
            instance.evolve_model(current_time)
            from_bse_to_model.copy()
            if not binary.type1 == previous_type1:
                results.append((binary.age, binary.mass1, binary.type1))
                previous_type1 = binary.type1
            
        self.assertEqual(len(results), 6)
        
        types = (
            "Hertzsprung Gap",
            "First Giant Branch",
            "Core Helium Burning",
            "First Asymptotic Giant Branch",
            "Hertzsprung Gap Naked Helium star",
            "Carbon/Oxygen White Dwarf",
        )
        
        for result, expected in zip(results, types):
            self.assertEquals(str(result[2]), expected)
        
        times = ( 
            284.8516 | units.Myr, 
            287.0595 | units.Myr, 
            287.7848 | units.Myr, 
            331.1454 | units.Myr, 
            331.3983 | units.Myr, 
            332.2786 | units.Myr,
        )
        for result, expected in zip(results, times):
            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 0)
            
        masses = ( 
            3.000 | units.MSun, 
            3.000 | units.MSun, 
            2.999 | units.MSun, 
            2.956 | units.MSun,
            0.888 | units.MSun,
            0.707 | units.MSun,
        )
        for result, expected in zip(results, masses):
            self.assertAlmostEqual(result[1].value_in(units.MSun), expected.value_in(units.MSun), 2)
         
        del instance
            
    def test2(self):
        print "Testing evolution of a wide binary system."
        instance = BSE()
        instance.parameters.metallicity = 0.001 | units.none
        instance.parameters.common_envelope_efficiency = 3.0 | units.none
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0 | units.none
        instance.initialize_module_with_current_parameters()
        stars =  core.Stars(1)
        
        binary = stars[0]
        binary.mass1 = 3.0 | units.MSun
        binary.mass2 = 0.3 | units.MSun
        binary.orbital_period = 2.0e5 | units.day
        binary.eccentricity = 0.0 | units.none

        instance.particles.add_particles(stars)
        from_bse_to_model = instance.particles.new_channel_to(stars)
        from_bse_to_model.copy()
        
        previous_type1 = binary.type1
        results = []
        current_time = 0 | units.Myr
        while current_time < (335 | units.Myr):
            instance.update_time_steps()
            # The next line appears a bit weird, but saves time for this simple test.
            current_time = current_time + max(2.0*instance.particles[0].time_step, 0.04 | units.Myr)
            instance.evolve_model(current_time)
            from_bse_to_model.copy()
    
            if not binary.type1 == previous_type1:
                results.append((binary.age, binary.mass1, binary.type1))
                previous_type1 = binary.type1
            
        print results
        self.assertEqual(len(results), 6)
        
        times = ( 
            284.8516 | units.Myr, 
            287.0595 | units.Myr, 
            287.7848 | units.Myr, 
            331.1454 | units.Myr,
            332.7407 | units.Myr,
            333.4146 | units.Myr
        )
        for result, expected in zip(results, times):
            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 0)
            
        masses = ( 
            3.000 | units.MSun, 
            3.000 | units.MSun, 
            2.999 | units.MSun, 
            2.956 | units.MSun,
            2.919 | units.MSun,
            0.928 | units.MSun
        )
        for result, expected in zip(results, masses):
            self.assertAlmostEqual(result[1].value_in(units.MSun), expected.value_in(units.MSun), 2)
         
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
            
    def test3(self):
        print "Testing standard BSE example 2..."
        instance = BSE()
        instance.parameters.common_envelope_efficiency = 3.0 | units.none
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0 | units.none
        instance.initialize_module_with_current_parameters()
        stars =  core.Stars(1)
        
        binary = stars[0]
        binary.mass1 = 7.816 | units.MSun
        binary.mass2 = 4.387 | units.MSun
        binary.orbital_period = 1964.18453 | units.day
        binary.eccentricity = 0.0 | units.none

        instance.particles.add_particles(stars)
        from_bse_to_model = instance.particles.new_channel_to(stars)
        from_bse_to_model.copy()
        
        previous_type1 = binary.type1
        previous_type2 = binary.type2
        results = []
        current_time = 0 | units.Myr        
        while current_time < (170 | units.Myr):
            instance.update_time_steps()
            # The next line appears a bit weird, but saves time for this simple test.
            current_time = current_time + max(2.0*instance.particles[0].time_step, 0.04 | units.Myr)
            instance.evolve_model(current_time)
            from_bse_to_model.copy()
    
            if not (binary.type1 == previous_type1 and binary.type2 == previous_type2):
                results.append((binary.age, str(binary.type1)+" and "+str(binary.type2)))
                previous_type1 = binary.type1
                previous_type2 = binary.type2
            
        print '\n'.join(map(str, results))
        self.assertEqual(len(results), 13)
        times = ( 
            38.9708 | units.Myr, 
            39.0897 | units.Myr, 
            39.1213 | units.Myr, 
            43.8025 | units.Myr,
            43.9923 | units.Myr,
            44.0686 | units.Myr,
            141.7077 | units.Myr, 
            142.3448 | units.Myr, 
            142.7827 | units.Myr,
            166.1043 | units.Myr,
            166.5795 | units.Myr,
            166.9627 | units.Myr,
            166.9863 | units.Myr
        )
        for result, expected in zip(results, times):
            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 0)
            
        types = (
            "Hertzsprung Gap and Main Sequence star",
            "First Giant Branch and Main Sequence star",
            "Core Helium Burning and Main Sequence star",
            "First Asymptotic Giant Branch and Main Sequence star",
            "Second Asymptotic Giant Branch and Main Sequence star",
            "Oxygen/Neon White Dwarf and Main Sequence star",
            "Oxygen/Neon White Dwarf and Hertzsprung Gap",
            "Oxygen/Neon White Dwarf and First Giant Branch",
            "Oxygen/Neon White Dwarf and Core Helium Burning",
            "Oxygen/Neon White Dwarf and First Asymptotic Giant Branch",
            "Oxygen/Neon White Dwarf and Hertzsprung Gap Naked Helium star",
            "Neutron Star and Hertzsprung Gap Naked Helium star",
            "Neutron Star and Carbon/Oxygen White Dwarf",
        )
        
        for result, expected in zip(results, types):
            self.assertEquals(result[1], expected)
        
        self.assertAlmostEqual(binary.mass1.value_in(units.MSun), 1.304, 3)
        self.assertAlmostEqual(binary.mass2.value_in(units.MSun), 0.800, 3)
        
        del instance
        
    def test4(self):
        print "Quick testing standard BSE example 2..."
        instance = BSE()
        instance.parameters.common_envelope_efficiency = 3.0 | units.none
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0 | units.none
        instance.initialize_module_with_current_parameters()
        stars =  core.Stars(1)
        
        binary = stars[0]
        binary.mass1 = 7.816 | units.MSun
        binary.mass2 = 4.387 | units.MSun
        binary.orbital_period = 1964.18453 | units.day
        binary.eccentricity = 0.0 | units.none

        instance.particles.add_particles(stars)
        from_bse_to_model = instance.particles.new_channel_to(stars)
        from_bse_to_model.copy()
        
        instance.evolve_model(170 | units.Myr)
        from_bse_to_model.copy()

        self.assertAlmostEqual(binary.mass1.value_in(units.MSun), 1.304, 3)
        self.assertAlmostEqual(binary.mass2.value_in(units.MSun), 0.800, 3)
        self.assertEquals(str(binary.type1), "Neutron Star")
        self.assertEquals(str(binary.type2), "Carbon/Oxygen White Dwarf")

        del instance
    
    def test5(self):
        print "Testing stellar collision..."
        instance = BSE()
        instance.parameters.common_envelope_efficiency = 3.0 | units.none
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0 | units.none
        instance.initialize_module_with_current_parameters()
        stars =  core.Stars(1)
        
        binary = stars[0]
        binary.mass1 = 3.0 | units.MSun
        binary.mass2 = 0.3 | units.MSun
        binary.orbital_period = 200.0 | units.day
        binary.eccentricity = 0.99 | units.none

        instance.particles.add_particles(stars)
        from_bse_to_model = instance.particles.new_channel_to(stars)
        from_bse_to_model.copy()
        
        instance.evolve_model(170 | units.Myr)
        from_bse_to_model.copy()

        self.assertAlmostEqual(binary.mass1.value_in(units.MSun), 3.300, 3)
        self.assertAlmostEqual(binary.mass2.value_in(units.MSun), 0.000, 3)
        self.assertEquals(str(binary.type1), "Main Sequence star")
        self.assertEquals(str(binary.type2), "Massless Supernova")

        del instance
        
    
    def test6(self):
        print "Testing additional parameters for initialization..."
        instance = BSE()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5 | units.none)
        myvalue = 0.7 | units.none
        instance.parameters.reimers_mass_loss_coefficient = myvalue
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, myvalue)
        instance.initialize_module_with_current_parameters()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, myvalue)
        del instance
        
        instance = BSE()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5 | units.none)
        myvalue = 0.7 | units.none
        instance.parameters.reimers_mass_loss_coefficient = myvalue
        instance.initialize_module_with_default_parameters()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5 | units.none)
        del instance
