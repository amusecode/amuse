import numpy

from amuse.test import amusetest

from amuse.ext.orbital_elements import (
        generate_binaries,
        get_orbital_elements_from_binaries,
        new_binary_from_orbital_elements,
        # get_orbital_elements_from_binary,
        orbital_elements_for_rel_posvel_arrays,
        orbital_elements,
        rel_posvel_arrays_from_orbital_elements,
        )

from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse import datamodel

from numpy import random


class KeplerTests(amusetest.TestCase):

    def test1(self):
        mass1 = 1 | nbody_system.mass
        mass2 = 1 | nbody_system.mass

        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length
        )

        self.assertEqual(len(binary), 2)

        binary.position -= binary[0].position
        binary.velocity -= binary[0].velocity

        self.assertAlmostRelativeEquals(
                binary[0].position,
                [0, 0, 0] | nbody_system.length)
        self.assertAlmostRelativeEquals(
                binary[1].position,
                [1, 0, 0] | nbody_system.length)
        self.assertAlmostRelativeEquals(
                binary[0].velocity,
                [0, 0, 0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(
                binary[1].velocity,
                [0, numpy.sqrt(2), 0] | nbody_system.speed)

    def test2(self):
        # test going around in a circular orbit
        mass1 = 1 | nbody_system.mass
        mass2 = 1 | nbody_system.mass

        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity=0,
            true_anomaly=90,
        )

        self.assertEqual(len(binary), 2)

        binary.position -= binary[0].position
        binary.velocity -= binary[0].velocity

        self.assertAlmostRelativeEquals(
                binary[0].position,
                [0, 0, 0] | nbody_system.length)
        self.assertAlmostRelativeEquals(
                binary[1].position,
                [0, 1, 0] | nbody_system.length)
        self.assertAlmostRelativeEquals(
                binary[0].velocity,
                [0, 0, 0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(
                binary[1].velocity,
                [-numpy.sqrt(2), 0, 0] | nbody_system.speed)

        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity=0,
            true_anomaly=180,
        )

        self.assertEqual(len(binary), 2)

        binary.position -= binary[0].position
        binary.velocity -= binary[0].velocity

        self.assertAlmostRelativeEquals(
                binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(
                binary[1].position, [-1,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(
                binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(
                binary[1].velocity, [0,-numpy.sqrt(2),0] | nbody_system.speed)
        
        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 270,
        )
        
        self.assertEqual(len(binary), 2)

        binary.position-=binary[0].position
        binary.velocity-=binary[0].velocity

        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0,-1,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [numpy.sqrt(2),0,0] | nbody_system.speed)
        
        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1 | nbody_system.length,
            eccentricity = 0,
            true_anomaly = 45,
        )
        
        self.assertEqual(len(binary), 2)

        binary.position-=binary[0].position
        binary.velocity-=binary[0].velocity

        self.assertAlmostRelativeEquals(binary[0].position, [0,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0.5 * numpy.sqrt(2),0.5 * numpy.sqrt(2),0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,0,0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [-1,1,0] | nbody_system.speed)


    def test3(self):
        mass1 = 1. | nbody_system.mass 
        mass2 = 1. | nbody_system.mass
        
        binary = new_binary_from_orbital_elements(
            mass1,
            mass2,
            1. | nbody_system.length,
            eccentricity = 0.
        )
        
        self.assertEqual(len(binary), 2)
        self.assertAlmostRelativeEquals(binary[0].position, [-0.5,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[1].position, [0.5,0,0] | nbody_system.length)
        self.assertAlmostRelativeEquals(binary[0].velocity, [0,-1/numpy.sqrt(2),0] | nbody_system.speed)
        self.assertAlmostRelativeEquals(binary[1].velocity, [0,1/numpy.sqrt(2),0] | nbody_system.speed)
    
    def test4(self):
        numpy.random.seed(3456789)
        N = 100

        mass1 = random.random(N) | nbody_system.mass
        mass2 = random.random(N) | nbody_system.mass
        semi_major_axis = (-numpy.log(random.random(N))) | nbody_system.length
        eccentricity = random.random(N)
        true_anomaly = (360.*random.random(N)-180.) | units.deg
        inclination = (180*random.random(N)) | units.deg
        longitude_of_the_ascending_node = (
                360*random.random(N)-180) | units.deg
        argument_of_periapsis = (360*random.random(N)-180) | units.deg

        for arg in zip(
                mass1, mass2, semi_major_axis, eccentricity, true_anomaly,
                inclination, longitude_of_the_ascending_node,
                argument_of_periapsis):
            arg_ = orbital_elements(
                    new_binary_from_orbital_elements(*arg))
            for i, (copy, org) in enumerate(zip(arg_, arg)):
                self.assertAlmostEqual(copy, org)

    def test5(self):
        numpy.random.seed(4567893)
        N = 100

        mass1 = random.random(N) | units.MSun
        mass2 = random.random(N) | units.MSun
        semi_major_axis = (-numpy.log(random.random(N))) | units.AU
        eccentricity = random.random(N)
        true_anomaly = (360.*random.random(N)-180.) | units.deg
        inclination = (180*random.random(N)) | units.deg
        longitude_of_the_ascending_node = (
                360*random.random(N)-180
                ) | units.deg
        argument_of_periapsis = (360*random.random(N)-180) | units.deg

        for arg in zip(
                mass1, mass2, semi_major_axis, eccentricity, true_anomaly,
                inclination, longitude_of_the_ascending_node,
                argument_of_periapsis):
            arg_ = orbital_elements(
                new_binary_from_orbital_elements(*arg, G=constants.G),
                G=constants.G)
            for i, (copy, org) in enumerate(zip(arg_, arg)):
                self.assertAlmostEqual(copy, org)
    
    def test6(self):
        """
        testing orbital_elements_for_rel_posvel_arrays for N particles 
        with random orbital elements
        """
        numpy.random.seed(666)
        N = 100
        
        mass_sun = 1. | units.MSun
        mass1 = numpy.ones(N) * mass_sun
        mass2 = numpy.zeros(N) | units.MSun
        semi_major_axis=(-numpy.log(random.random(N))) | units.AU 
        eccentricity = random.random(N)
        true_anomaly = 360.*random.random(N)-180.
        inclination = 180*random.random(N)
        longitude_of_the_ascending_node = 360*random.random(N)-180
        argument_of_periapsis = 360*random.random(N)-180       
        
        comets = datamodel.Particles(N)
        for i,arg in enumerate(zip(mass1,mass2,semi_major_axis,eccentricity,true_anomaly,inclination, 
                                   longitude_of_the_ascending_node,argument_of_periapsis)):
            sun_and_comet = new_binary_from_orbital_elements(*arg,G=constants.G)
            comets[i].mass = sun_and_comet[1].mass
            comets[i].position = sun_and_comet[1].position
            comets[i].velocity = sun_and_comet[1].velocity
        
        semi_major_axis_ext, eccentricity_ext, ta_ext, inclination_ext, \
        longitude_of_the_ascending_node_ext, argument_of_periapsis_ext = \
        orbital_elements_for_rel_posvel_arrays(comets.position,
                                               comets.velocity,
                                               comets.mass + mass_sun,
                                               G=constants.G)        

        self.assertAlmostEqual(semi_major_axis,semi_major_axis_ext)
        self.assertAlmostEqual(eccentricity,eccentricity_ext)
        self.assertAlmostEqual(inclination ,inclination_ext)
        self.assertAlmostEqual(longitude_of_the_ascending_node ,longitude_of_the_ascending_node_ext)
        self.assertAlmostEqual(argument_of_periapsis ,argument_of_periapsis_ext)
        self.assertAlmostEqual(true_anomaly,ta_ext)

    def test7(self):
        """
        testing orbital_elements_for_rel_posvel_arrays for the case of one particle
        """
        numpy.random.seed(999)
        
        mass1 = 0.5 | units.MSun
        mass2 = 0.8 | units.MSun
        sem = 12. | units.AU
        ecc = 0.05
        inc = 20.
        lon = 10.
        arg = 0.4
        ta = 360.*random.random()-180.
        
        binary = new_binary_from_orbital_elements(mass1, 
                                                  mass2,
                                                  sem,
                                                  ecc,
                                                  ta,
                                                  inc,
                                                  lon,
                                                  arg,
                                                  G=constants.G)

        rel_pos = binary[1].position - binary[0].position
        rel_vel = binary[1].velocity - binary[0].velocity
        mass_12 = binary[1].mass + binary[0].mass
        sem_ext, ecc_ext, ta_ext, inc_ext, lon_ext, arg_ext = \
        orbital_elements_for_rel_posvel_arrays(rel_pos, rel_vel, mass_12, G=constants.G)
        
        self.assertAlmostEqual(sem, sem_ext)
        self.assertAlmostEqual(ecc, ecc_ext)
        self.assertAlmostEqual(inc, inc_ext)
        self.assertAlmostEqual(lon, lon_ext)
        self.assertAlmostEqual(arg, arg_ext)
        self.assertAlmostEqual(ta,ta_ext)
    
    def test8(self):
        """
        testing orbital_elements_for_rel_posvel_arrays for extreme cases
        """
        N = 3
        mass1 = (1.2*numpy.ones(N)) | units.MSun
        mass2 = (0.1, 0.05, 0.003) | units.MSun
        semi_major_axis = (1., 2., 3.) | units.AU
        eccentricity = (0., 0.5, 0.6)
        true_anomaly = (0., 0., 66.)
        inclination = (12., 0., 180.)
        longitude_of_the_ascending_node = (0., 0., 0.,)
        argument_of_periapsis = (0., 23., 90.)
        mass12 = []
        rel_position = []
        rel_velocity = []
        for i,arg in enumerate(zip(mass1,mass2,semi_major_axis,eccentricity,true_anomaly,inclination, 
                                   longitude_of_the_ascending_node,argument_of_periapsis)):
            sun_and_comet = new_binary_from_orbital_elements(*arg,G=constants.G)
            mass12.append(sun_and_comet[0].mass + sun_and_comet[1].mass)
            rel_position.append(sun_and_comet[1].position - sun_and_comet[0].position)
            rel_velocity.append(sun_and_comet[1].velocity - sun_and_comet[0].velocity)
            
        # to convert lists to vector quantities
        rel_pos = numpy.array([vec_i.value_in(units.AU) for vec_i in rel_position]) | units.AU
        rel_vel = numpy.array([vec_i.value_in(units.kms) for vec_i in rel_velocity]) | units.kms
        mass_12 = numpy.array([m_i.value_in(units.MSun) for m_i in mass12]) | units.MSun
        
        semi_major_axis_ext, eccentricity_ext, ta_ext, inclination_ext, \
        longitude_of_the_ascending_node_ext, argument_of_periapsis_ext = \
        orbital_elements_for_rel_posvel_arrays(rel_pos,
                                               rel_vel,
                                               mass_12,
                                               G=constants.G)        

        self.assertAlmostEqual(semi_major_axis,semi_major_axis_ext)
        self.assertAlmostEqual(eccentricity,eccentricity_ext)
        self.assertAlmostEqual(inclination,inclination_ext)
        self.assertAlmostEqual(longitude_of_the_ascending_node,longitude_of_the_ascending_node_ext)
        self.assertAlmostEqual(argument_of_periapsis,argument_of_periapsis_ext)
        self.assertAlmostEqual(true_anomaly,ta_ext)

    def test9(self):
        """
        testing orbital_elements_for_rel_posvel_arrays for N particles 
        with random orbital elements, nbody_system
        """
        numpy.random.seed(666)
        N = 100
        
        mass_sun = 1. | nbody_system.mass
        mass1 = numpy.ones(N) * mass_sun
        mass2 = numpy.zeros(N) | nbody_system.mass
        semi_major_axis=(-numpy.log(random.random(N))) | nbody_system.length 
        eccentricity = random.random(N)
        true_anomaly = 360.*random.random(N)-180.
        inclination = 180*random.random(N)
        longitude_of_the_ascending_node = 360*random.random(N)-180
        argument_of_periapsis = 360*random.random(N)-180       
        
        comets = datamodel.Particles(N)
        for i,arg in enumerate(zip(mass1,mass2,semi_major_axis,eccentricity,true_anomaly,inclination, 
                                   longitude_of_the_ascending_node,argument_of_periapsis)):
            sun_and_comet = new_binary_from_orbital_elements(*arg,G=nbody_system.G)
            comets[i].mass = sun_and_comet[1].mass
            comets[i].position = sun_and_comet[1].position
            comets[i].velocity = sun_and_comet[1].velocity
        
        semi_major_axis_ext, eccentricity_ext, ta_ext, inclination_ext, \
        longitude_of_the_ascending_node_ext, argument_of_periapsis_ext = \
        orbital_elements_for_rel_posvel_arrays(comets.position,
                                               comets.velocity,
                                               comets.mass + mass_sun,
                                               G=nbody_system.G)
        self.assertAlmostEqual(semi_major_axis,semi_major_axis_ext)
        self.assertAlmostEqual(eccentricity,eccentricity_ext)
        self.assertAlmostEqual(inclination,inclination_ext)
        self.assertAlmostEqual(longitude_of_the_ascending_node,longitude_of_the_ascending_node_ext)
        self.assertAlmostEqual(argument_of_periapsis,argument_of_periapsis_ext)
        self.assertAlmostEqual(true_anomaly,ta_ext)

    def xtest10(self):
        """
        testing orbital_elements_for_rel_posvel_arrays for N particles 
        with random orbital elements, unitless
        """
        numpy.random.seed(666)
        N = 100
        
        mass_sun = 1. 
        mass1 = numpy.ones(N) * mass_sun
        mass2 = numpy.zeros(N) 
        semi_major_axis=(-numpy.log(random.random(N)))  
        eccentricity = random.random(N)
        true_anomaly = 360.*random.random(N)-180.
        inclination = 180*random.random(N)
        longitude_of_the_ascending_node = 360*random.random(N)-180
        argument_of_periapsis = 360*random.random(N)-180       
        
        comets = datamodel.Particles(N)
        for i,arg in enumerate(zip(mass1,mass2,semi_major_axis,eccentricity,true_anomaly,inclination, 
                                   longitude_of_the_ascending_node,argument_of_periapsis)):
            sun_and_comet = new_binary_from_orbital_elements(*arg,G=1)
            comets[i].mass = sun_and_comet[1].mass
            comets[i].position = sun_and_comet[1].position
            comets[i].velocity = sun_and_comet[1].velocity
        
        semi_major_axis_ext, eccentricity_ext, ta_ext, inclination_ext, \
        longitude_of_the_ascending_node_ext, argument_of_periapsis_ext = \
        orbital_elements_for_rel_posvel_arrays(comets.position,
                                               comets.velocity,
                                               comets.mass + mass_sun,
                                               G=1)

        self.assertAlmostEqual(semi_major_axis,semi_major_axis_ext)
        self.assertAlmostEqual(eccentricity,eccentricity_ext)
        self.assertAlmostEqual(inclination,inclination_ext)
        self.assertAlmostEqual(longitude_of_the_ascending_node,longitude_of_the_ascending_node_ext)
        self.assertAlmostEqual(argument_of_periapsis,argument_of_periapsis_ext)
        self.assertAlmostEqual(true_anomaly,ta_ext)
            
            
    def test11(self):
        """
        testing orbital_elements_for_rel_posvel_arrays for unbound orbits
        """
        
        from amuse.community.kepler.interface import Kepler
        
        numpy.random.seed(66)
        N = 10
        
        mass_sun = 1. | units.MSun
        mass1 = numpy.ones(N) * mass_sun
        mass2 = numpy.zeros(N) | units.MSun
        semi_major_axis=-1000.*(random.random(N)) | units.AU 
        eccentricity = (1.+random.random(N))*10.-9.
        inclination = numpy.pi*random.random(N)
        longitude_of_the_ascending_node = 2.*numpy.pi*random.random(N)-numpy.pi
        argument_of_periapsis = 2.*numpy.pi*random.random(N)-numpy.pi      
        
        # kepler.initialize_from_elements initializes orbits with mean_anomaly=0 and true_anomaly=0
        true_anomaly = 0.*(360.*random.random(N)-180.)
        
        comets = datamodel.Particles(N)
        
        converter = nbody_system.nbody_to_si(1|units.MSun,1|units.AU)
        kepler = Kepler(converter)
        kepler.initialize_code()
        for i,arg in enumerate(zip(mass1,mass2,semi_major_axis,eccentricity,true_anomaly,inclination, 
                                   longitude_of_the_ascending_node,argument_of_periapsis)):
            kepler.initialize_from_elements(mass=(mass1[i]+mass2[i]),
                                            semi=semi_major_axis[i],
                                            ecc=eccentricity[i])
            ri = kepler.get_separation_vector()
            vi = kepler.get_velocity_vector()
            
            om = longitude_of_the_ascending_node[i]
            w = argument_of_periapsis[i]
            incl = inclination[i]
            a1 = ([numpy.cos(om), -numpy.sin(om), 0.0], [numpy.sin(om), numpy.cos(om), 0.0], [0.0, 0.0, 1.0])
            a2 = ([1.0, 0.0, 0.0], [0.0, numpy.cos(incl), -numpy.sin(incl)], [0.0, numpy.sin(incl), numpy.cos(incl)])
            a3 = ([numpy.cos(w), -numpy.sin(w), 0.0], [numpy.sin(w), numpy.cos(w), 0.0], [0.0, 0.0, 1.0])
            A = numpy.dot(numpy.dot(a1,a2),a3)

            r_vec = numpy.dot(A,numpy.reshape(ri,3,1))
            v_vec = numpy.dot(A,numpy.reshape(vi,3,1))
          
            r = (0.0, 0.0, 0.0) | units.AU
            v = (0.0, 0.0, 0.0) | (units.AU / units.day)
            r[0] = r_vec[0]
            r[1] = r_vec[1]
            r[2] = r_vec[2]
            v[0] = v_vec[0]
            v[1] = v_vec[1]
            v[2] = v_vec[2]
  
            comets[i].mass = mass2[i]
            comets[i].position = r_vec
            comets[i].velocity = v_vec
        
        kepler.stop()
        
        semi_major_axis_ext, eccentricity_ext, ta_ext, inclination_ext, \
        longitude_of_the_ascending_node_ext, argument_of_periapsis_ext = \
        orbital_elements(comets.position,
                                               comets.velocity,
                                               comets.mass + mass_sun,
                                               G=constants.G)
        
        self.assertAlmostEqual(semi_major_axis,semi_major_axis_ext.in_(units.AU))
        self.assertAlmostEqual(eccentricity,eccentricity_ext)
        self.assertAlmostEqual(inclination,inclination_ext)
        self.assertAlmostEqual(longitude_of_the_ascending_node,longitude_of_the_ascending_node_ext)
        self.assertAlmostEqual(argument_of_periapsis,argument_of_periapsis_ext)
        self.assertAlmostEqual(true_anomaly,ta_ext)

    def test12(self):
        """
        tests generating cartesian coordinates from orbital elements
        """
        numpy.random.seed(1701)

        mass1 = 1.0 | units.MSun
        mass2 = 0.1 | units.MEarth
        sem = 2. | units.AU
        ecc = 0.15
        inc = 11. | units.deg
        lon = 30. | units.deg
        arg = 0.3 | units.deg
        ta = (360.*random.random()-180.) | units.deg

        rel_pos, rel_vel = rel_posvel_arrays_from_orbital_elements(
                mass1,
                mass2,
                sem,
                ecc,
                ta,
                inc,
                lon,
                arg,
                G=constants.G)

        mass_12 = mass1 + mass2
        sem_ext, ecc_ext, ta_ext, inc_ext, lon_ext, arg_ext = \
            orbital_elements(
                    rel_pos, rel_vel, mass_12, G=constants.G)

        self.assertAlmostEqual(sem, sem_ext)
        self.assertAlmostEqual(ecc, ecc_ext)
        self.assertAlmostEqual(inc, inc_ext)
        self.assertAlmostEqual(lon, lon_ext)
        self.assertAlmostEqual(arg, arg_ext)
        self.assertAlmostEqual(ta, ta_ext)

    def test13(self):
        """
        tests generating cartesian coordinates from orbital elements
        """
        numpy.random.seed(17014)
        N = 5

        mass1 = 1.0 | units.MSun
        mass2 = numpy.ones(N) * 0.01 | units.MEarth
        sem = 2. | units.AU
        ecc = 0.15
        inc = 11. | units.deg
        lon = 30. | units.deg
        arg = 0.3 | units.deg
        ta = (360.*random.random()-180.) | units.deg

        rel_pos, rel_vel = rel_posvel_arrays_from_orbital_elements(
                mass1,
                mass2,
                sem,
                ecc,
                ta,
                inc,
                lon,
                arg,
                G=constants.G)

        mass_12 = mass1 + mass2
        sem_ext, ecc_ext, ta_ext, inc_ext, lon_ext, arg_ext = \
            orbital_elements(
                    rel_pos, rel_vel, mass_12, G=constants.G)

        self.assertAlmostEqual(
                sem, sem_ext)
        self.assertAlmostEqual(ecc, ecc_ext)
        self.assertAlmostEqual(inc, inc_ext)
        self.assertAlmostEqual(lon, lon_ext)
        self.assertAlmostEqual(arg, arg_ext)
        self.assertAlmostEqual(ta, ta_ext)

    def test14(self):
        """
        tests generating cartesian coordinates from orbital elements
        """
        numpy.random.seed(17018)
        N = 5

        mass1 = numpy.ones(N) * 1.0 | units.MSun
        mass2 = random.random(N) | units.MEarth
        sem = numpy.array([2., 1.0, 1.1, 1.2, 4.0]) | units.AU
        ecc = numpy.array([0.15, 0.01, 0.5, 0.9, 0.99])
        inc = numpy.array([11., 0.1, 20, 90, 180.]) | units.deg
        lon = numpy.array([31., 32., 33., 45., 30.]) | units.deg
        arg = numpy.array([0.3, 11., 15., 30., 95.]) | units.deg
        ta = (360.*random.random(N)-180.) | units.deg

        rel_pos, rel_vel = rel_posvel_arrays_from_orbital_elements(
                mass1,
                mass2,
                sem,
                ecc,
                ta,
                inc,
                lon,
                arg,
                G=constants.G)

        mass_12 = mass1 + mass2
        sem_ext, ecc_ext, ta_ext, inc_ext, lon_ext, arg_ext = \
            orbital_elements(
                    rel_pos, rel_vel, mass_12, G=constants.G)
        self.assertAlmostEqual(
                sem.value_in(units.AU), sem_ext.value_in(units.AU))
        self.assertAlmostEqual(ecc, ecc_ext)
        self.assertAlmostEqual(inc, inc_ext)
        self.assertAlmostEqual(lon, lon_ext)
        self.assertAlmostEqual(arg, arg_ext)
        self.assertAlmostEqual(ta, ta_ext)

    def test15(self):
        """
        testing orbital_elements_for_rel_posvel_arrays for N particles
        with random orbital elements
        """
        numpy.random.seed(666)
        N = 100

        mass_sun = 1. | units.MSun
        mass1 = numpy.ones(N) * mass_sun
        mass2 = numpy.zeros(N) | units.MSun
        semi_major_axis = (-numpy.log(random.random(N))) | units.AU
        eccentricity = random.random(N)
        true_anomaly = 360.*random.random(N)-180.
        inclination = 180*random.random(N)
        longitude_of_the_ascending_node = 360*random.random(N)-180
        argument_of_periapsis = 360*random.random(N)-180

        comets = datamodel.Particles(N)
        suns = datamodel.Particles(N)
        for i, arg in enumerate(
                zip(mass1, mass2, semi_major_axis, eccentricity, true_anomaly,
                    inclination, longitude_of_the_ascending_node,
                    argument_of_periapsis)):
            sun_and_comet = new_binary_from_orbital_elements(
                    *arg, G=constants.G)
            comets[i].mass = sun_and_comet[1].mass
            comets[i].position = sun_and_comet[1].position
            comets[i].velocity = sun_and_comet[1].velocity

        suns.mass = mass1
        suns.position = 0*comets.position
        suns.velocity = 0*comets.velocity

        mass1_ext, mass2_ext, semi_major_axis_ext, eccentricity_ext, ta_ext,\
            inclination_ext, longitude_of_the_ascending_node_ext,\
            argument_of_periapsis_ext = orbital_elements(
                    suns, comets, G=constants.G)
        rad_to_deg = 180./numpy.pi
        for i in range(N):
            self.assertAlmostEqual(
                    semi_major_axis[i].value_in(units.AU),
                    semi_major_axis_ext[i].value_in(units.AU))
            self.assertAlmostEqual(eccentricity[i], eccentricity_ext[i])
            self.assertAlmostEqual(
                    inclination[i], rad_to_deg*inclination_ext[i])
            self.assertAlmostEqual(
                    longitude_of_the_ascending_node[i],
                    rad_to_deg*longitude_of_the_ascending_node_ext[i])
            self.assertAlmostEqual(
                    argument_of_periapsis[i],
                    rad_to_deg*argument_of_periapsis_ext[i])
            self.assertAlmostEqual(true_anomaly[i], rad_to_deg*ta_ext[i])

    def test16(self):
        """ tests a mismatch in shape in generate_binaries """
        m1=[1]*5 | nbody_system.mass
        m2=[0]*5 | nbody_system.mass
        a=[1.]*5 | nbody_system.length
        ecc=numpy.array([0,0,.99999,0.1,0.5])
        ta=[180,180,20,30,0]| units.deg
        primaries,secondaries=generate_binaries(m1,m2,a,eccentricity=ecc,true_anomaly=ta)
        m1_,m2_,a_,ecc_,ta_,i_,lasc_,ap_= get_orbital_elements_from_binaries(primaries,secondaries)
        self.assertAlmostEqual(ecc,ecc_)
      
      
