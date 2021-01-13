import numpy

from amuse.test.amusetest import TestWithMPI

from amuse.community.seba.interface import SeBaInterface, SeBa

from amuse.units import units
from amuse.units import constants
from amuse.datamodel import Particle
from amuse.datamodel import Particles

class TestSeBaInterface(TestWithMPI):

    def test1(self):
        instance = self.new_instance_of_an_optional_code(SeBaInterface)

        endtime, mass, radius, luminosity, temperature, time_step, stellar_type, error = instance.evolve_star(1, 4600, 0.02)
        self.assertEqual(error, 0)
        self.assertTrue( endtime <= 4600.0)
        self.assertAlmostRelativeEqual(endtime, 4600.0, 4)
        self.assertAlmostRelativeEqual(mass, 1.0, 6)
        self.assertAlmostRelativeEqual(radius, 0.9856, 4)
        self.assertAlmostRelativeEqual(luminosity, 0.9585, 4)
        self.assertAlmostRelativeEqual(temperature, 5751, 4)
        self.assertAlmostRelativeEqual(time_step, 1089.3, 4)
        self.assertEqual(stellar_type, 1)

        instance.stop()

    def test2(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEqual(error, 0)

        index,error = instance.new_particle(1.)
        self.assertEqual(error, 0)
        self.assertEqual(index, 1)
        mass, error = instance.get_mass(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 1.0, 6)
        value, error = instance.get_radius(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(value, 0.88824945029751212, 6)

        stellar_type, error = instance.get_stellar_type(index)
        self.assertEqual(error, 0)
        self.assertEqual(stellar_type, 1)

        instance.stop()

    def test3(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEqual(error, 0)

        index,error = instance.new_particle(1.)
        self.assertEqual(error, 0)
        self.assertEqual(index, 1)
        error = instance.evolve_model(4600)
        self.assertEqual(error, 0)

        mass, error = instance.get_mass(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 1.0, 6)
        value, error = instance.get_radius(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(value, 0.9856, 4)
        value, error = instance.get_temperature(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(value, 5751, 4)
        value, error = instance.get_time_step(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(value, 1089.3, 4)

        stellar_type, error = instance.get_stellar_type(index)
        self.assertEqual(error, 0)
        self.assertEqual(stellar_type, 1)

        instance.stop()

    def test4(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEqual(error, 0)

        index,error = instance.new_particle(1.)
        self.assertEqual(error, 0)
        self.assertEqual(index, 1)
        for t in range(46):
            error = instance.evolve_model((t+1) * 100)
            self.assertEqual(error, 0)

        mass, error = instance.get_mass(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 1.0, 6)
        value, error = instance.get_radius(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(value, 0.9856, 4)
        value, error = instance.get_temperature(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(value, 5751, 4)
        value, error = instance.get_time_step(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(value, 1089.3, 4)

        stellar_type, error = instance.get_stellar_type(index)
        self.assertEqual(error, 0)
        self.assertEqual(stellar_type, 1)

        instance.stop()

    def test5(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEqual(error, 0)

        index,error = instance.new_particle([1., 2., 3.])
        self.assertEqual(error, 0)
        self.assertEqual(index, [1,2,3])

        mass, error = instance.get_mass(2)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 2 , 6)

        mass, error = instance.get_mass(3)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 3, 6)


        error = instance.evolve_model(4600)
        self.assertEqual(error, 0)

        mass, error = instance.get_mass(index)
        print(mass)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass[0], 1.0, 6)
        self.assertAlmostRelativeEqual(mass[1], 0.62973, 4)
        self.assertAlmostRelativeEqual(mass[2], 0.75012, 4)


        instance.stop()

    def test6(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEqual(error, 0)

        index,error = instance.new_particle([1., 2., 3.])
        self.assertEqual(error, 0)
        self.assertEqual(index, [1,2,3])

        for t in range(46):
            error = instance.evolve_model((t+1) * 100)
            self.assertEqual(error, 0)

        mass, error = instance.get_mass(index)
        print(mass)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, [1.0, 0.62973, 0.75072], 4)

        instance.stop()

    def test7(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEqual(error, 0)

        index,error = instance.new_particle([1., 2., 3.])
        self.assertEqual(error, 0)
        self.assertEqual(index, [1,2,3])

        mass, error = instance.get_mass(2)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 2 , 6)

        mass, error = instance.get_mass(3)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 3, 6)

        mass, error = instance.get_mass(4)
        self.assertEqual(error, -1)

        error = instance.delete_star(2)
        self.assertEqual(error, 0)

        mass, error = instance.get_mass(2)
        self.assertEqual(error, -1)

        mass, error = instance.get_mass(3)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 3, 6)

        index, error = instance.new_particle(4.)
        self.assertEqual(error, 0)
        self.assertEqual(index, 4)

        instance.stop()

    def test8(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEqual(error, 0)

        index,error = instance.new_particle([3.0,1.0,2.0])
        self.assertEqual(error, 0)
        self.assertEqual(index, [1,2,3])

        error = instance.delete_star(1)
        self.assertEqual(error, 0)

        error = instance.evolve_model(4600);

        mass, error = instance.get_mass(2)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 1, 6)

        error = instance.delete_star(3)
        self.assertEqual(error, 0)


        index,error = instance.new_particle([5.0])
        self.assertEqual(error, 0)

        mass, error = instance.get_mass(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 5.0, 6)
        error = instance.evolve_model(5000);


        mass, error = instance.get_mass(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 0.99057, 4)

        error = instance.delete_star(2)
        self.assertEqual(error, 0)
        error = instance.delete_star(index)
        self.assertEqual(error, 0)

        for i in range(4):
            mass, error = instance.get_mass(index+1)
            self.assertEqual(error, -1)

        index,error = instance.new_particle([5.0])
        self.assertEqual(error, 0)

        error = instance.evolve_model(10000);

        mass, error = instance.get_mass(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 0.99057, 4)

        instance.stop()

    def test9(self):
        instance = SeBaInterface() #self.new_instance_of_an_optional_code(SeBaInterface)
        error = instance.initialize_code()
        self.assertEqual(error, 0)
        instance.set_metallicity(0.001)

        index,error = instance.new_particle([3.0,0.3])
        self.assertEqual(error, 0)
        self.assertEqual(index, [1,2])

        mu = (3.3 | units.MSun) * constants.G
        orbital_period = 200.0 | units.day
        semi_major_axis = (((orbital_period / 2.0 * numpy.pi)**2)*mu)**(1.0/3.0)
        print(semi_major_axis.value_in(units.RSun))


        eccentricity = 0.5
        index,error = instance.new_binary(
            semi_major_axis.value_in(units.RSun),
            eccentricity,
            index[0],
            index[1]
        )
        self.assertEqual(error, 0)
        self.assertEqual(index, 3)

        mass, error = instance.get_mass(index)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 3.3, 4)
        mass, error = instance.get_mass(2)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 0.3, 4)

        error = instance.evolve_model(300)
        self.assertEqual(error, 0)
        mass, error = instance.get_mass(1)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 2.98777, 4)
        mass, error = instance.get_mass(2)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 0.29999, 4)


        error = instance.evolve_model(400)
        self.assertEqual(error, 0)
        mass, error = instance.get_mass(1)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 0.86679, 4)
        mass, error = instance.get_mass(2)
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(mass, 0.3, 4)

        error = instance.delete_binary(index)
        self.assertEqual(error, 0)
        mass, error = instance.get_mass(index)
        self.assertEqual(error, -1)

        # check if singles are still in the mode and evolve
        value, error = instance.get_age([1,2])
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(value, 400, 4)
        error = instance.evolve_model(500)
        self.assertEqual(error, 0)
        value, error = instance.get_age([1,2])
        self.assertEqual(error, 0)
        self.assertAlmostRelativeEqual(value, 500, 4)

class TestSeBa(TestWithMPI):

    def test_evolve_star(self):

        instance = self.new_instance_of_an_optional_code(SeBa)

        endtime, mass, radius, luminosity, temperature, time_step, stellar_type = instance.evolve_star(1 | units.MSun, 4600 | units.Myr, 0.02)

        self.assertTrue( endtime <= 4600 | units.Myr)
        self.assertAlmostRelativeEqual(mass, 1.0 | units.MSun, 4)
        self.assertAlmostRelativeEqual(radius, 0.9856 | units.RSun, 4)
        self.assertAlmostRelativeEqual(luminosity, 0.9585 | units.LSun, 4)
        self.assertAlmostRelativeEqual(temperature, 5751 | units.K, 4)
        self.assertAlmostRelativeEqual(time_step, 1089.3 | units.Myr, 4)
        self.assertEqual(stellar_type, 1 | units.stellar_type)


    def test_add_particle(self):
        instance = self.new_instance_of_an_optional_code(SeBa)

        p = Particle()
        p.mass = 5 | units.MSun
        p.metallicity = 0.02

        p = instance.particles.add_particle(p)
        instance.evolve_model(130 | units.Myr)
        print(p)

        self.assertAlmostRelativeEqual(p.mass, 0.9906 | units.MSun, 4)

    def test_evolution_of_close_binary_system(self):
        # print("Testing evolution of a close binary system...")
        instance = self.new_instance_of_an_optional_code(SeBa)
        instance.commit_parameters()
        stars =  Particles(2)
        stars[0].mass = 3.0 | units.MSun
        stars[1].mass = 0.3 | units.MSun


        mu = (3.3 | units.MSun) * constants.G
        orbital_period = 200.0 | units.day
        semi_major_axis = (((orbital_period / (2.0 * numpy.pi))**2)*mu)**(1.0/3.0)

        instance.particles.add_particles(stars)

        binaries =  Particles(1)

        binary = binaries[0]
        binary.semi_major_axis = semi_major_axis
        binary.eccentricity = 0.5
        binary.child1 = stars[0]
        binary.child2 = stars[1]

        instance.binaries.add_particles(binaries)

        from_seba_to_model = instance.particles.new_channel_to(stars)
        from_seba_to_model.copy()

        from_seba_to_model_binaries = instance.binaries.new_channel_to(binaries)
        from_seba_to_model_binaries.copy()

        previous_type = binary.child1.stellar_type
        results = []
        current_time = 0 | units.Myr

        while current_time < (480 | units.Myr):
            instance.update_time_steps()
            # The next line appears a bit weird, but saves time for this simple test.
            deltat = max(1.0*instance.binaries[0].time_step, 0.1| units.Myr)
            current_time = current_time + deltat
            instance.evolve_model(current_time)
            from_seba_to_model.copy()
            from_seba_to_model_binaries.copy()
            if not binary.child1.stellar_type == previous_type:
                results.append((binary.age, binary.child1.mass, binary.child1.stellar_type))
                previous_type = binary.child1.stellar_type

        self.assertEqual(len(results), 6)
        for x in results:
            print(x)

        types = (
            "Hertzsprung Gap",
            "First Giant Branch",
            "Core Helium Burning",
            "First Asymptotic Giant Branch",
            "Giant Branch Naked Helium star",
            "Carbon/Oxygen White Dwarf",
        )


        for result, expected in zip(results, types):
            self.assertEqual(str(result[2]), expected)

        times = (
            377.6369 | units.Myr,
            379.8877 | units.Myr,
            382.3112 | units.Myr,
            473.4804 | units.Myr,
            475.4766 | units.Myr,
            476.6182 | units.Myr,
        )
        for result, expected in zip(results, times):
            self.assertAlmostEqual(result[0].value_in(units.Myr), expected.value_in(units.Myr), 0)

        masses = (
            3.0000 | units.MSun,
            3.0000 | units.MSun,
            2.9983 | units.MSun,
            2.9741 | units.MSun,
            0.6710 | units.MSun,
            0.6596 | units.MSun,
        )
        for result, expected in zip(results, masses):
            self.assertAlmostEqual(result[1].value_in(units.MSun), expected.value_in(units.MSun), 2)

        instance.stop()


    def test_set_parameter_metallicity(self):
        instance = self.new_instance_of_an_optional_code(SeBa)
        self.assertAlmostRelativeEquals(instance.parameters.metallicity , 0.02)
        instance.parameters.metallicity = 0.04
        self.assertAlmostRelativeEquals(instance.parameters.metallicity , 0.04)

    def test_set_parameter_logging(self):
        instance = self.new_instance_of_an_optional_code(SeBa)
        self.assertFalse(instance.parameters.is_logging_of_evolve_enabled)
        instance.parameters.is_logging_of_evolve_enabled = True
        self.assertTrue(instance.parameters.is_logging_of_evolve_enabled)

    def test_add_binary_particles(self):
        instance = self.new_instance_of_an_optional_code(SeBa)
        instance.commit_parameters()
        stars =  Particles(2)
        stars[0].mass = 3.0 | units.MSun
        stars[1].mass = 0.3 | units.MSun


        mu = (3.3 | units.MSun) * constants.G
        orbital_period = 200.0 | units.day
        semi_major_axis = (((orbital_period / (2.0 * numpy.pi))**2)*mu)**(1.0/3.0)

        instance.particles.add_particles(stars)

        binaries =  Particles(1)

        binary = binaries[0]
        binary.semi_major_axis = semi_major_axis
        binary.eccentricity = 0.5
        binary.child1 = stars[0]
        binary.child2 = stars[1]

        instance.binaries.add_particles(binaries)

        self.assertAlmostRelativeEquals(instance.binaries[0].child1.mass, 3.0 | units.MSun, 4)
        self.assertAlmostRelativeEquals(instance.binaries[0].child2.mass, 0.3 | units.MSun, 4)


    def xtest7(self):
        instance = self.new_instance_of_an_optional_code(SeBa)
        instance.parameters.metallicity = 0.03
        p = Particle()
        p.mass = 99.1605930967 | units.MSun

        p = instance.particles.add_particle(p)
        instance.evolve_model(614 | units.Myr)
        print(p.stellar_type)
        self.assertEqual(str(p.stellar_type),'Black Hole')
        self.assertAlmostRelativeEqual(p.mass, 0.9906 | units.MSun, 4)

    def test_set_semi_major_axis_eccentricity(self):
        instance = self.new_instance_of_an_optional_code(SeBa)
        instance.parameters.supernova_kick_velocity = 0 | units.kms
        instance.commit_parameters()
        print("v_kick=", instance.parameters.supernova_kick_velocity)
        stars = Particles(2)
        stars[0].mass = 10 | units.MSun
        stars[1].mass = 9 | units.MSun

        semi_major_axis = 10000 | units.AU

        instance.particles.add_particles(stars)

        binaries = Particles(1)

        binary = binaries[0]
        binary.semi_major_axis = semi_major_axis
        binary.eccentricity = 0
        binary.child1 = stars[0]
        binary.child2 = stars[1]

        instance.binaries.add_particles(binaries)
        instance.evolve_model(30 | units.Myr)
        print(instance.particles)
        print(instance.binaries)

        self.assertAlmostRelativeEquals(instance.binaries[0].eccentricity, 1, 4)

    def test_add_stars_at_different_times(self):
        instance = self.new_instance_of_an_optional_code(SeBa)
        stars = Particles(2)
        stars[0].mass = 10.0 | units.MSun
        stars[1].mass = 9 | units.MSun

        instance.particles.add_particles(stars)
        instance.evolve_model(30|units.Myr)

        self.assertAlmostRelativeEquals(instance.particles.age, [30,30] |units.Myr)
        self.assertAlmostRelativeEquals(instance.model_time, 30 | units.Myr)
        self.assertAlmostRelativeEquals(instance.particles[0].mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEquals(instance.particles[1].mass, 8.8682 | units.MSun, 4)
        stars = Particles(2)
        stars[0].mass = 10.0 | units.MSun
        stars[1].mass = 9 | units.MSun

        instance.particles.add_particles(stars)
        instance.evolve_model(60 | units.Myr)
        print(instance.particles.age)
        print(instance.particles.mass)
        self.assertAlmostRelativeEquals(instance.model_time, 60 | units.Myr)
        self.assertAlmostRelativeEquals(instance.particles.age, [60,60,30,30] | units.Myr)
        self.assertAlmostRelativeEquals(instance.particles[2].mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEquals(instance.particles[3].mass, 8.8682 | units.MSun, 4)

    def test_supernova_stopping_condition(self):
        """ Test supernova stopping condition """
        instance = self.new_instance_of_an_optional_code(SeBa)
        instance.stopping_conditions.supernova_detection.enable()

        p = Particle()
        p.mass = 10 | units.MSun
        p.metallicity = 0.02

        p = instance.particles.add_particle(p)
        instance.set_supernova_kick_velocity(0.0|units.kms)
        instance.evolve_model(30 | units.Myr)
        self.assertEqual(instance.stopping_conditions.supernova_detection.is_set(), True)
        self.assertEqual(instance.stopping_conditions.supernova_detection.particles(0)[0].key, p.key)

        self.assertAlmostRelativeEqual(p.age, 27.35866 | units.Myr, 4)

        self.assertAlmostRelativeEqual(p.mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEqual(p.natal_kick_velocity, [-54.070,12.818,42.583] | units.kms, 4)


    def test_supernova_stopping_condition_in_a_binary(self):
        """ Test supernova stopping condition in a binary """
        instance = self.new_instance_of_an_optional_code(SeBa)
        instance.stopping_conditions.supernova_detection.enable()

        stars =  Particles(2)
        stars[0].mass = 10.0 | units.MSun
        stars[1].mass = 5.0 | units.MSun
        instance.particles.add_particles(stars)

        binaries =  Particles(1)
        binary = binaries[0]
        binary.semi_major_axis = 1.e6|units.RSun
        binary.eccentricity = 0.1
        binary.child1 = stars[0]
        binary.child2 = stars[1]
        instance.binaries.add_particles(binaries)

        instance.evolve_model(30 | units.Myr)
        self.assertEqual(instance.stopping_conditions.supernova_detection.is_set(), True)
        self.assertEqual(instance.stopping_conditions.supernova_detection.particles(0)[0].key, instance.binaries[0].child1.key)
        self.assertEqual(instance.stopping_conditions.supernova_detection.particles(0)[0].key, instance.particles[0].key)

        print(instance.parameters)
        self.assertAlmostRelativeEqual(instance.particles[0].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[0].mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].mass, 5.0 | units.MSun, 4)
#        self.assertAlmostRelativeEqual(instance.particles[0].natal_kick_velocity, [0,0,0] | units.kms, 4)

        self.assertAlmostRelativeEqual(instance.binaries[0].child1.age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.binaries[0].child2.age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.binaries[0].child1.mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.binaries[0].child2.mass, 5.0 | units.MSun, 4)
#        self.assertAlmostRelativeEqual(instance.binaries[0].child1.natal_kick_velocity, [0,0,0] | units.kms, 4)

    def test_supernova_stopping_condition_with_multiple_stars(self):
        """ Test supernova stopping condition with multiple stars """
        instance = self.new_instance_of_an_optional_code(SeBa)
        instance.stopping_conditions.supernova_detection.enable()

        stars =  Particles(3)
        stars[0].mass = 10.0 | units.MSun
        stars[1].mass = 5.0 | units.MSun
        stars[2].mass = 0.5 | units.MSun
        instance.particles.add_particles(stars)


        instance.evolve_model(30 | units.Myr)
        self.assertEqual(instance.stopping_conditions.supernova_detection.is_set(), True)
        self.assertEqual(instance.stopping_conditions.supernova_detection.particles(0)[0].key, instance.particles[0].key)

        self.assertAlmostRelativeEqual(instance.particles[0].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[0].mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].mass, 5.0 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].mass, 0.5 | units.MSun, 4)
#        self.assertAlmostRelativeEqual(instance.particles[0].natal_kick_velocity, [0,0,0] | units.kms, 4)

    def test_supernova_stopping_condition_with_multiple_stars_multiple_supernovae(self):
        """ Test supernova stopping condition with multiple stars """
        instance = self.new_instance_of_an_optional_code(SeBa)
        instance.stopping_conditions.supernova_detection.enable()

        stars =  Particles(3)
        stars[0].mass = 10.0 | units.MSun
        stars[1].mass = 11.0 | units.MSun
        stars[2].mass = 0.5 | units.MSun
        instance.particles.add_particles(stars)


        instance.evolve_model(30 | units.Myr)
        self.assertEqual(instance.stopping_conditions.supernova_detection.is_set(), True)
        self.assertEqual(instance.stopping_conditions.supernova_detection.particles(0)[0].key, instance.particles[1].key)

        self.assertAlmostRelativeEqual(instance.particles[0].age, 23.08688 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].age, 23.08688 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].age, 23.08688 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[0].mass, 9.92275 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].mass, 1.28809 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].mass, 0.5 | units.MSun, 4)
#        self.assertAlmostRelativeEqual(instance.particles[0].natal_kick_velocity, [0,0,0] | units.kms, 4)

        instance.evolve_model(30 | units.Myr)
        self.assertEqual(instance.stopping_conditions.supernova_detection.is_set(), True)
        self.assertEqual(instance.stopping_conditions.supernova_detection.particles(0)[0].key, instance.particles[0].key)

        self.assertAlmostRelativeEqual(instance.particles[0].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[0].mass, 1.24766 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].mass, 1.28809 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].mass, 0.5 | units.MSun, 4)
#        self.assertAlmostRelativeEqual(instance.particles[0].natal_kick_velocity, [0,0,0] | units.kms, 4)

    def test_supernova_stopping_condition_with_multiple_stars_of_equal_mass(self):
        """ Test supernova stopping condition with multiple stars of equal mass """
        instance = self.new_instance_of_an_optional_code(SeBa)
        instance.stopping_conditions.supernova_detection.enable()

        stars =  Particles(3)
        stars[0].mass = 10.0 | units.MSun
        stars[1].mass = 10.0 | units.MSun
        stars[2].mass = 0.5 | units.MSun
        instance.particles.add_particles(stars)


        instance.evolve_model(30 | units.Myr)
        self.assertEqual(instance.stopping_conditions.supernova_detection.is_set(), True)
        self.assertEqual(len(instance.stopping_conditions.supernova_detection.particles(0)), 2)
        self.assertEqual(instance.stopping_conditions.supernova_detection.particles(0)[0].key, instance.particles[0].key)
        self.assertEqual(instance.stopping_conditions.supernova_detection.particles(0)[1].key, instance.particles[1].key)

        self.assertAlmostRelativeEqual(instance.particles[0].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].age, 27.35866 | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[0].mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].mass, 0.5 | units.MSun, 4)

#        self.assertAlmostRelativeEqual(instance.particles[0].natal_kick_velocity, [0,0,0] | units.kms, 4)
#        self.assertAlmostRelativeEqual(instance.particles[1].natal_kick_velocity, [0,0,0] | units.kms, 4)

        instance.evolve_model(30 | units.Myr)
        self.assertEqual(instance.stopping_conditions.supernova_detection.is_set(), False)

        self.assertAlmostRelativeEqual(instance.particles[0].age, 30. | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].age, 30. | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].age, 30. | units.Myr, 4)
        self.assertAlmostRelativeEqual(instance.particles[0].mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[1].mass, 1.2476 | units.MSun, 4)
        self.assertAlmostRelativeEqual(instance.particles[2].mass, 0.5 | units.MSun, 4)

