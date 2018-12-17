from amuse.community.mobse.interface import MOBSE, MOBSEInterface

from amuse.test.amusetest import TestWithMPI
from amuse.units import units
from amuse.units import constants
from amuse.datamodel import Particles

import numpy

class TestMOBSEInterface(TestWithMPI):
    
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
        instance = MOBSEInterface()
        metallicity = 0.02
        neta = 0.5
        bwind =  0.0
        hewind =  0.5
        alpha1 = 1.0
        CElambda = 0.1
        ceflag = 0
        tflag = 1
        ifflag = 0
        wdflag =  1
        bhflag =  1 
        nsflag =  3
        piflag =  1
        mxns =  3.0
        idum = 29769
        pts1 = 0.05
        pts2 = 0.01
        pts3 = 0.02
        sigma1 =  265.0
        sigma2 =  7.0
        beta = 1.0/8.0
        xi = 1.0
        acc2 = 3.0/2.0
        epsnov = 0.001
        eddfac = 1.0
        gamma = -1.0

        status = instance.initialize(metallicity,
            neta, bwind, hewind, alpha1, CElambda,
            ceflag, tflag, ifflag, wdflag, bhflag,
            nsflag, piflag, mxns, idum, pts1, pts2, pts3,
            sigma1,sigma2,beta,xi,acc2,epsnov,eddfac,gamma)
        self.assertEqual(status,0)
        instance.stop()
        
    def test2(self):
        print "Test basic operations (legacy functions evolve & get_time_step)..."
        instance = MOBSEInterface()
        status = instance.initialize(0.02, 0.5, 0.0, 0.5, 1.0, 0.1, 0, 1, 0, 1, 1, 3, 1, 3.0,
            29769, 0.05, 0.01, 0.02, 265.0, 7.0, 1.0/8.0, 1.0, 3.0/2.0, 0.001, 1.0, -1.0)
        
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
        
        result = instance.evolve_binary(
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
            'spin2' : '0x1.07413b0522dabp+10',
        };    

        for x in expected:
            print "'%s' : '%s'," % (x, getattr(updated_state, x).hex())
            self.assertAlmostRelativeEqual(float.fromhex(expected[x]),getattr(updated_state, x))
            
        self.assertEquals(updated_state.end_time, 1e-06)
        dt = instance.get_time_step(updated_state.type1, updated_state.type2,
            updated_state.initial_mass1, updated_state.initial_mass2, updated_state.mass1,
            updated_state.mass2, updated_state.t_ms1, updated_state.t_ms2,
            updated_state.epoch1, updated_state.epoch2, updated_state.bse_age)
        self.assertAlmostEqual(dt, 18.8768, 3)
        instance.stop()
     
    def test3(self):
        print "Test whether the interface can handle arrays..."
        instance = MOBSEInterface()
        status = instance.initialize(0.02, 0.5, 0.0, 0.5, 1.0, 0.1, 0, 1, 0, 1, 1, 3, 1, 3.0,
            29769, 0.05, 0.01, 0.02, 265.0, 7.0, 1.0/8.0, 1.0, 3.0/2.0, 0.001, 1.0, -1.0)
        masses1 = [10.0,5.0,4.0,40.0,130.0]
        masses2 = [1.0,1.0,1.0,10.0,60.0]
        types1 = types2 = [1,1,1,1,1]
        orbital_periods = [100.0,200.0,300.0,500.0,600.0]
        eccentricities = [0.5,0.6,0.7,0.3,0.2]

        radii1 = luminosity1 = core_mass1 = core_radius1 =  envelope_mass1 =\
        envelope_radius1 = spin1 = epoch1 = t_ms1 = [0.0,0.0,0.0,0.0,0.0]
        radii2 = luminosity2 = core_mass2 = core_radius2 =  envelope_mass2 =\
        envelope_radius2 = spin2 = epoch2 = t_ms2 = [0.0,0.0,0.0,0.0,0.0]
        init_mass1 = masses1
        init_mass2 = masses2
        bse_age = [0.0,0.0,0.0,0.0,0.0]
        end_time = [10.0, 10.0, 10.0, 10.0, 10.0]
        result = instance.evolve_binary(
            types1, types2, init_mass1, init_mass2,
            masses1, masses2, radii1, radii2,
            luminosity1, luminosity2, core_mass1, core_mass2,
            core_radius1, core_radius2, envelope_mass1, envelope_mass2,
            envelope_radius1, envelope_radius2, spin1, spin2,
            epoch1, epoch2, t_ms1, t_ms2, 
            bse_age, orbital_periods, eccentricities, end_time
        )
        self.assertAlmostEqual(result['mass1'][0], 9.99356, 2)
        self.assertAlmostEqual(result['mass1'][1], 4.99956, 2)
        self.assertAlmostEqual(result['mass1'][2], 3.99992, 2)
        self.assertAlmostEqual(result['mass1'][3], 3.07374, 2)
        self.assertAlmostEqual(result['mass1'][4], 13.61644, 2)
        instance.stop()
        
    def test4(self):
        print "Test large number of particles..."
        number_of_particles = 2000
        instance = MOBSEInterface()
        status = instance.initialize(0.02, 0.5, 0.0, 0.5, 1.0, 0.1, 0, 1, 0, 1, 1, 3, 1, 3.0,
            29769, 0.05, 0.01, 0.02, 265.0, 7.0, 1.0/8.0, 1.0, 3.0/2.0, 0.001, 1.0, -1.0)
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
        
        result = instance.evolve_binary(
            types1, types2, init_mass1, init_mass2,
            masses1, masses2, radii1, radii2,
            luminosity1, luminosity2, core_mass1, core_mass2,
            core_radius1, core_radius2, envelope_mass1, envelope_mass2,
            envelope_radius1, envelope_radius2, spin1, spin2,
            epoch1, epoch2, t_ms1, t_ms2, 
            bse_age, orbital_periods, eccentricities, end_time
        )
        self.assertEquals(len(result['mass1']), number_of_particles)
        instance.stop()

        
class TestMOBSE(TestWithMPI):
    
    def test1(self):
        print "Testing evolution of a close binary system..."
        instance = MOBSE()
        instance.initialize_code()
        instance.parameters.metallicity = 0.001
        instance.parameters.common_envelope_efficiency = 3.0
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0
        instance.commit_parameters()
        stars =  Particles(2)
        stars[0].mass = 3.0 | units.MSun
        stars[1].mass = 0.3 | units.MSun
        
        orbital_period = 200.0 | units.day
        semi_major_axis = instance.orbital_period_to_semi_major_axis(orbital_period,  stars[0].mass , stars[1].mass)
        
        instance.particles.add_particles(stars)
        
        binaries =  Particles(1)
        
        binary = binaries[0]
        binary.semi_major_axis = semi_major_axis
        binary.eccentricity = 0.5
        binary.child1 = stars[0]
        binary.child2 = stars[1]
        
        instance.binaries.add_particles(binaries)
        
        from_mobse_to_model = instance.particles.new_channel_to(stars)
        from_mobse_to_model.copy()

        from_mobse_to_model_binaries = instance.binaries.new_channel_to(binaries)
        from_mobse_to_model_binaries.copy()
        
        previous_type = binary.child1.stellar_type
        results = []
        current_time = 0 | units.Myr
        
        while current_time < (480 | units.Myr):
            instance.update_time_steps()
            # The next line appears a bit weird, but saves time for this simple test.
            current_time = current_time + max(5.0*instance.binaries[0].time_step, 0.3 | units.Myr)
            instance.evolve_model(current_time)
            from_mobse_to_model.copy()
            from_mobse_to_model_binaries.copy()
            if not binary.child1.stellar_type == previous_type:
                results.append((binary.age, binary.child1.mass, binary.child1.stellar_type))
                previous_type = binary.child1.stellar_type
            
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
            284.8632 | units.Myr, 
            287.0713 | units.Myr, 
            287.7967 | units.Myr, 
            331.1631 | units.Myr, 
            331.4164 | units.Myr, 
            332.2864 | units.Myr,
        )
        for result, expected in zip(results, times):
            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 0)
            
        masses = ( 
            3.000 | units.MSun, 
            3.000 | units.MSun, 
            2.999 | units.MSun, 
            2.956 | units.MSun,
            0.888 | units.MSun,
            0.701 | units.MSun,
        )
        for result, expected in zip(results, masses):
            self.assertAlmostEqual(result[1].value_in(units.MSun), expected.value_in(units.MSun), 2)
         
        instance.stop()
            
    def test2(self):
        print "Testing evolution of a wide binary system."
        instance = MOBSE()
        instance.parameters.metallicity = 0.001
        instance.parameters.common_envelope_efficiency = 3.0
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0
        instance.commit_parameters()
        
        stars =  Particles(2)
        stars[0].mass = 4.0 | units.MSun
        stars[1].mass = 0.3 | units.MSun
        orbital_period =  2.0e5 | units.day
        semi_major_axis = instance.orbital_period_to_semi_major_axis(orbital_period,  stars[0].mass , stars[1].mass)
        
        instance.particles.add_particles(stars)
        
        binaries =  Particles(1)
        
        binary = binaries[0]
        binary.semi_major_axis = semi_major_axis
        binary.eccentricity = 0.5
        binary.child1 = stars[0]
        binary.child2 = stars[1]
        
        instance.binaries.add_particles(binaries)
        
        from_mobse_to_model = instance.particles.new_channel_to(stars)
        from_mobse_to_model.copy()

        from_mobse_to_model_binaries = instance.binaries.new_channel_to(binaries)
        from_mobse_to_model_binaries.copy()
        
        previous_type = binary.child1.stellar_type
        results = []
        current_time = 0 | units.Myr
        
        while current_time < (170 | units.Myr):
            instance.update_time_steps()
            # The next line appears a bit weird, but saves time for this simple test.
            current_time = current_time + max(2.0*instance.binaries[0].time_step, 0.04 | units.Myr)
            instance.evolve_model(current_time)
            from_mobse_to_model.copy()
            from_mobse_to_model_binaries.copy()
            if not binary.child1.stellar_type == previous_type:
                results.append((binary.age, binary.child1.mass, binary.child1.stellar_type))
                previous_type = binary.child1.stellar_type
        print results
        self.assertEqual(len(results), 6)
        
        times = ( 
            147.1282 | units.Myr, 
            148.0345 | units.Myr, 
            148.2282 | units.Myr, 
            167.2811 | units.Myr,
            168.0344 | units.Myr,
            168.7475 | units.Myr
        )
        for result, expected in zip(results, times):
            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 0)
            
        masses = ( 
            4.000 | units.MSun, 
            3.999 | units.MSun, 
            3.999 | units.MSun, 
            3.942 | units.MSun,
            3.906 | units.MSun,
            1.016 | units.MSun
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
        
        instance.stop()
            
    def test3(self):
        print "Testing standard MOBSE example 2..."
        instance = MOBSE()
        instance.parameters.common_envelope_efficiency = 3.0
        instance.parameters.common_envelope_binding_energy_factor= 0.5
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0
        instance.commit_parameters()
        
        stars =  Particles(2)
        stars[0].mass = 7.816 | units.MSun
        stars[1].mass = 4.387 | units.MSun
        
        orbital_period =   1964.18453 | units.day
        semi_major_axis = instance.orbital_period_to_semi_major_axis(orbital_period,  stars[0].mass , stars[1].mass)
        instance.particles.add_particles(stars)
        
        binaries =  Particles(1)
        
        binary = binaries[0]
        binary.semi_major_axis = semi_major_axis
        binary.eccentricity = 0.0
        binary.child1 = stars[0]
        binary.child2 = stars[1]
        
        instance.binaries.add_particles(binaries)
        
        from_mobse_to_model = instance.particles.new_channel_to(stars)
        from_mobse_to_model.copy()

        from_mobse_to_model_binaries = instance.binaries.new_channel_to(binaries)
        from_mobse_to_model_binaries.copy()
        
        previous_type1 = binary.child1.stellar_type
        previous_type2 = binary.child2.stellar_type
        results = []
        current_time = 0 | units.Myr
        
        while current_time < (170 | units.Myr):
            instance.update_time_steps()
            # The next line appears a bit weird, but saves time for this simple test.
            current_time = current_time + max(2.0*instance.binaries[0].time_step, 0.04 | units.Myr)
            instance.evolve_model(current_time)
            from_mobse_to_model.copy()
            from_mobse_to_model_binaries.copy()        
            if not (binary.child1.stellar_type  == previous_type1 and binary.child2.stellar_type == previous_type2):
                results.append((binary.age, str(binary.child1.stellar_type)+" and "+str(binary.child2.stellar_type)))
                previous_type1 = binary.child1.stellar_type
                previous_type2 = binary.child2.stellar_type
        
            
        print '\n'.join(map(str, results))
        self.assertEqual(len(results), 12)
        times = ( 
            39.1037 | units.Myr, 
            39.2242 | units.Myr, 
            39.2565 | units.Myr, 
            43.9911 | units.Myr,
            44.1842 | units.Myr,
            44.2644 | units.Myr,
            141.8444 | units.Myr, 
            142.4835 | units.Myr, 
            142.9234 | units.Myr,
            166.3238 | units.Myr,
            166.8385 | units.Myr,
            167.1731 | units.Myr
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
            "Neutron Star and Carbon/Oxygen White Dwarf",
        )
        
        for result, expected in zip(results, types):
            self.assertEquals(result[1], expected)
        
        self.assertAlmostEqual(binary.child1.mass.value_in(units.MSun), 1.26079, 3)
        self.assertAlmostEqual(binary.child2.mass.value_in(units.MSun), 0.76080, 3)
        
        instance.stop()
        
    def test4(self):
        print "Quick testing standard MOBSE example 2..."
        instance = MOBSE()
        instance.parameters.common_envelope_efficiency = 3.0
        instance.parameters.common_envelope_binding_energy_factor= 0.5
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0
        instance.commit_parameters()
        
        stars =  Particles(2)
        stars[0].mass = 7.816 | units.MSun
        stars[1].mass = 4.387 | units.MSun
        
        instance.particles.add_particles(stars)
        
        binaries =  Particles(1)
        
        binary = binaries[0]
        orbital_period =   1964.18453 | units.day
        semi_major_axis = instance.orbital_period_to_semi_major_axis(orbital_period,  stars[0].mass , stars[1].mass)
        binary.semi_major_axis = semi_major_axis
        binary.eccentricity = 0.0
        binary.child1 = stars[0]
        binary.child2 = stars[1]
        
        instance.binaries.add_particles(binaries)
        
        from_mobse_to_model = instance.particles.new_channel_to(stars)
        from_mobse_to_model.copy()

        from_mobse_to_model_binaries = instance.binaries.new_channel_to(binaries)
        from_mobse_to_model_binaries.copy()
        
        instance.evolve_model(170 | units.Myr)
        from_mobse_to_model.copy()
        from_mobse_to_model_binaries.copy()

        self.assertAlmostEqual(binary.child1.mass.value_in(units.MSun), 1.26079, 3)
        self.assertAlmostEqual(binary.child2.mass.value_in(units.MSun), 0.76080, 3)
        self.assertEquals(str(binary.child1.stellar_type), "Neutron Star")
        self.assertEquals(str(binary.child2.stellar_type), "Carbon/Oxygen White Dwarf")

        instance.stop()
    
    def test5(self):
        print "Testing stellar collision..."
        instance = MOBSE()
        instance.parameters.common_envelope_efficiency = 3.0
        instance.parameters.Eddington_mass_transfer_limit_factor = 10.0
        instance.commit_parameters()

        stars =  Particles(2)
        stars[0].mass = 130.0  | units.MSun
        stars[1].mass = 50 | units.MSun
        
        instance.particles.add_particles(stars)
        
        binaries =  Particles(1)
        
        binary = binaries[0]
        orbital_period =   300.0 | units.day
        semi_major_axis = instance.orbital_period_to_semi_major_axis(orbital_period,  stars[0].mass , stars[1].mass)
        binary.semi_major_axis = semi_major_axis
        binary.eccentricity = 0.99
        binary.child1 = stars[0]
        binary.child2 = stars[1]
        
        instance.binaries.add_particles(binaries)
        
        from_mobse_to_model = instance.particles.new_channel_to(stars)
        from_mobse_to_model.copy()

        from_mobse_to_model_binaries = instance.binaries.new_channel_to(binaries)
        from_mobse_to_model_binaries.copy()
        
        instance.evolve_model(170 | units.Myr)
        
        from_mobse_to_model.copy()
        from_mobse_to_model_binaries.copy()
        print binaries
        self.assertAlmostEqual(binary.child1.mass.value_in(units.MSun), 180.00, 3)
        self.assertAlmostEqual(binary.child2.mass.value_in(units.MSun), 0.000, 3)
        self.assertEquals(str(binary.child1.stellar_type), "Main Sequence star")
        self.assertEquals(str(binary.child2.stellar_type), "Massless Supernova")

        instance.stop()

    def test6(self):
        print "Testing additional parameters for initialization..."
        instance = MOBSE()
        instance.initialize_code()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5)
        myvalue = 0.7
        instance.parameters.reimers_mass_loss_coefficient = myvalue
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, myvalue)
        instance.commit_parameters()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, myvalue)
        instance.stop()
        
        instance = MOBSE()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5)
        myvalue = 0.7
        instance.parameters.reimers_mass_loss_coefficient = myvalue
        instance.parameters.set_defaults()
        instance.commit_parameters()
        self.assertEqual(instance.parameters.reimers_mass_loss_coefficient, 0.5)
        instance.stop()
    
    def test7(self):
        print "Test evolve_model optional arguments: end_time and keep_synchronous"

        instance = MOBSE()
        instance.commit_parameters()
        
        stars =  Particles(6)
        stars.mass = [1.0,2.0,3.0, 0.1, 0.2, 0.3]  | units.MSun
        
        binaries =  Particles(3)
        binaries.eccentricity = 0.0
        for i in range(3):
            binaries[i].child1 = stars[i]
            binaries[i].child2 = stars[i+3]
        orbital_period =   200.0 | units.day
        semi_major_axis = instance.orbital_period_to_semi_major_axis(
            orbital_period,  
            binaries.child1.as_set().mass , 
            binaries.child2.as_set().mass
        )
        binaries.semi_major_axis = semi_major_axis
        
        instance.particles.add_particles(stars)
        instance.binaries.add_particles(binaries)
        
        self.assertAlmostEqual(instance.binaries.age, [0.0, 0.0, 0.0] | units.yr)
        self.assertAlmostEqual(instance.binaries.time_step, [550.1565, 58.2081, 18.8768] | units.Myr, 3)
        
        print "evolve_model without arguments: use shared timestep = min(particles.time_step)"
        instance.evolve_model()
        self.assertAlmostEqual(instance.binaries.age, [18.8768, 18.8768, 18.8768] | units.Myr, 3)
        self.assertAlmostEqual(instance.binaries.time_step, [550.1565, 58.2081, 18.8768] | units.Myr, 3)
        self.assertAlmostEqual(instance.model_time, 18.8768 | units.Myr, 3)
        
        print "evolve_model with end_time: take timesteps, until end_time is reached exactly"
        instance.evolve_model(100 | units.Myr)
        self.assertAlmostEqual(instance.binaries.age, [100.0, 100.0, 100.0] | units.Myr, 3)
        self.assertAlmostEqual(instance.binaries.time_step, [550.1565, 58.2081, 18.8785] | units.Myr, 3)
        self.assertAlmostEqual(instance.model_time, 100.0 | units.Myr, 3)
        
        print "evolve_model with keep_synchronous: use non-shared timestep, particle ages will typically diverge"
        instance.evolve_model(keep_synchronous = False)
        self.assertAlmostEqual(instance.binaries.age, (100 | units.Myr) + ([550.1565, 58.2081, 18.8785] | units.Myr), 3)
        self.assertAlmostEqual(instance.binaries.time_step, [550.1565, 58.2081, 18.8785] | units.Myr, 3)
        self.assertAlmostEqual(instance.model_time, 100.0 | units.Myr, 3) # Unchanged!
        instance.stop()
        
    def test8(self):
        print "Testing adding and removing particles from stellar evolution code..."
        
        instance = MOBSE()
        instance.initialize_code()
        
        stars =  Particles(6)
        stars.mass = [1.0,1.0, 1.0, 0.2, 0.2, 0.2]  | units.MSun
        
        binaries =  Particles(3)
        binaries.eccentricity = 0.0
        for i in range(3):
            binaries[i].child1 = stars[i]
            binaries[i].child2 = stars[i+3]
        orbital_period =   200.0 | units.day
        semi_major_axis = instance.orbital_period_to_semi_major_axis(
            orbital_period,  
            binaries.child1.as_set().mass , 
            binaries.child2.as_set().mass
        )
        binaries.semi_major_axis = semi_major_axis

        instance.commit_parameters()
        self.assertEquals(len(instance.particles), 0)
        self.assertEquals(len(instance.binaries), 0) # before creation
        instance.particles.add_particles(stars)
        instance.binaries.add_particles(binaries[:-1])
        instance.commit_particles()
        instance.evolve_model(1.0 | units.Myr)
        self.assertEquals(len(instance.binaries), 2) # before remove
        self.assertAlmostEqual(instance.binaries.age, 1.0 | units.Myr)
        
        instance.binaries.remove_particle(binaries[0])
        self.assertEquals(len(instance.binaries), 1)
        instance.evolve_model(2.0 | units.Myr)
        self.assertAlmostEqual(instance.binaries[0].age, 2.0 | units.Myr)
        
        instance.binaries.add_particles(binaries[::2])
        self.assertEquals(len(instance.binaries), 3) # it's back...
        self.assertAlmostEqual(instance.binaries[0].age, 2.0 | units.Myr)
        self.assertAlmostEqual(instance.binaries[1].age, 0.0 | units.Myr)
        self.assertAlmostEqual(instance.binaries[2].age, 0.0 | units.Myr) # ... and rejuvenated.
        
        instance.evolve_model(3.0 | units.Myr) # The young stars keep their age offset from the old star
        self.assertAlmostEqual(instance.binaries.age, [3.0, 1.0, 1.0] | units.Myr)
        instance.evolve_model(4.0 | units.Myr)
        self.assertAlmostEqual(instance.binaries.age, [4.0, 2.0, 2.0] | units.Myr)
        instance.stop()
    
    def test9(self):
        print "Testing MOBSE states"
        instance = MOBSE()
        
        stars =  Particles(2)
        stars.mass = [1.0, 0.2]  | units.MSun
        
        binaries =  Particles(1)
        orbital_period =   200.0 | units.day
        semi_major_axis = instance.orbital_period_to_semi_major_axis(orbital_period,  stars[0].mass , stars[1].mass)
        binaries.semi_major_axis = semi_major_axis
        binaries.eccentricity = 0.0
        binaries[0].child1 = stars[0]
        binaries[0].child2 = stars[1]
        
        print "First do everything manually:",
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
        instance = MOBSE()
        self.assertEquals(instance.get_name_of_current_state(), 'UNINITIALIZED')
        instance.parameters.reimers_mass_loss_coefficient = 0.5
        self.assertEquals(instance.get_name_of_current_state(), 'INITIALIZED')
        instance.particles.add_particles(stars)
        instance.binaries.add_particles(binaries)
        self.assertEquals(instance.get_name_of_current_state(), 'RUN')
        instance.stop()
        self.assertEquals(instance.get_name_of_current_state(), 'STOPPED')
        print "ok"

