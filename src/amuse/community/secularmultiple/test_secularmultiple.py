from amuse.test.amusetest import TestWithMPI
import os
import sys
import numpy
import math

from amuse.community.secularmultiple.interface import SecularMultipleInterface, SecularMultiple

from amuse.units import units,quantities,constants
from amuse.datamodel import Particles

HAS_MATPLOTLIB = False  # disable plotting in the test script

def parse_arguments():
    from amuse.units.optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--test",
                    dest="test",type="int",default=1,
                    help="")

    options, args = parser.parse_args()
    return options.__dict__

class TestSecularMultipleInterface(TestWithMPI):

    def test0(self):
        instance = SecularMultipleInterface()
        instance.stop()

def create_binary(m1,m2,a,e,i,g,h):
    particles = Particles(3)

    for index in range(2):
        particle = particles[index]
        particle.is_binary = False
        if index == 0:
            particle.mass = m1
        else:
            particle.mass = m2
        particle.child1 = None
        particle.child2 = None

    particles[2].is_binary = True
    particles[2].semimajor_axis = a
    particles[2].eccentricity = e
    particles[2].inclination = i
    particles[2].argument_of_pericenter = g
    particles[2].longitude_of_ascending_node = h
    particles[2].child1 = particles[0]
    particles[2].child2 = particles[1]

    return particles

def create_triple(m1,m2,m3,a_in,a_out,e_in,e_out,i_in,i_out,g_in,g_out,h_in,h_out):
    N_bodies = 3
    N_binaries = 2

    m_list = [m1,m2,m3]
    a_list = [a_in,a_out]
    e_list = [e_in,e_out]
    i_list = [i_in,i_out]
    g_list = [g_in,g_out]
    h_list = [h_in,h_out]

    particles = Particles(N_bodies + N_binaries)

    for index in range(N_bodies):
        particle = particles[index]
        particle.is_binary = False
        particle.mass = m_list[index]
        particle.child1 = None
        particle.child2 = None

    for index in range(N_binaries):
        particles[index+N_bodies].is_binary = True
        particles[index+N_bodies].semimajor_axis = a_list[index]
        particles[index+N_bodies].eccentricity = e_list[index]
        particles[index+N_bodies].inclination = i_list[index]
        particles[index+N_bodies].argument_of_pericenter = g_list[index]
        particles[index+N_bodies].longitude_of_ascending_node = h_list[index]

        if index==0:
            particles[index+N_bodies].child1 = particles[0]
            particles[index+N_bodies].child2 = particles[1]
        elif index==1:
            particles[index+N_bodies].child1 = particles[2]
            particles[index+N_bodies].child2 = particles[3]

    return particles

def create_quadruple_triple_single(m1,m2,m3,m4,a_A,a_B,a_C,e_A,e_B,e_C,i_A,i_B,i_C,g_A,g_B,g_C,h_A,h_B,h_C):
    N_bodies = 4
    N_binaries = 3

    m_list = [m1,m2,m3,m4]
    a_list = [a_A,a_B,a_C]
    e_list = [e_A,e_B,e_C]
    i_list = [i_A,i_B,i_C]
    g_list = [g_A,g_B,g_C]
    h_list = [h_A,h_B,h_C]

    particles = Particles(N_bodies + N_binaries)

    for index in range(N_bodies):
        particle = particles[index]
        particle.is_binary = False
        particle.mass = m_list[index]
        particle.child1 = None
        particle.child2 = None

    for index in range(N_binaries):
        particle = particles[index+N_bodies]
        particle.is_binary = True
        particle.semimajor_axis = a_list[index]
        particle.eccentricity = e_list[index]
        particle.inclination = i_list[index]
        particle.argument_of_pericenter = g_list[index]
        particle.longitude_of_ascending_node = h_list[index]

        if index==0:
            particle.child1 = particles[1]
            particle.child2 = particles[0]
        elif index==1:
            particle.child1 = particles[2]
            particle.child2 = particles[4]
        elif index==2:
            particle.child1 = particles[5]
            particle.child2 = particles[3]

    return particles


def create_nested_multiple(N,masses,semimajor_axes,eccentricities,inclinations,arguments_of_pericentre,longitudes_of_ascending_node):
    """
    N is number of bodies
    masses should be N-sized array
    the other arguments should be (N-1)-sized arrays
    """

    N_bodies = N
    N_binaries = N-1

    particles = Particles(N_bodies + N_binaries)

    for index in range(N_bodies):
        particle = particles[index]
        particle.is_binary = False
        particle.mass = masses[index]
        particle.child1 = None
        particle.child2 = None

    previous_binary = particles
    for index in range(N_binaries):
        particle = particles[index+N_bodies]
        particle.is_binary = True
        particle.semimajor_axis = semimajor_axes[index]
        particle.eccentricity = eccentricities[index]
        particle.inclination = inclinations[index]
        particle.argument_of_pericenter = arguments_of_pericentre[index]
        particle.longitude_of_ascending_node = longitudes_of_ascending_node[index]


        if index==0:
            particle.child1 = particles[0]
            particle.child2 = particles[1]
        else:
            particle.child1 = previous_binary
            particle.child2 = particles[index+1]

        previous_binary = particle

    return particles

class TestSecularMultiple(TestWithMPI):
    def test0(self):
        instance = SecularMultiple()
        print(instance.parameters)


    def test1(self):
        """
        test reference system of Naoz et al. (2009)
        """
        particles = create_nested_multiple(3,[1.0|units.MSun, 1.0|units.MJupiter, 40.0|units.MJupiter], [6.0|units.AU,100.0|units.AU], [0.001,0.6], [0.0001,65.0*numpy.pi/180.0],[45.0*numpy.pi/180.0,0.0001],[0.01,0.01])
        binaries = particles[particles.is_binary]

        binaries.include_pairwise_1PN_terms = False

        code = SecularMultiple(redirection='none')
        code.particles.add_particles(particles)

        channel_to_code = particles.new_channel_to(code.particles)
        channel_from_code = code.particles.new_channel_to(particles)


        channel_to_code.copy()
        self.assertEqual(0.6, code.particles[particles.is_binary][1].eccentricity)

        t = 0.0 | units.Myr
        dt = 1.0e-2 | units.Myr
        tend = 3.0e0 | units.Myr

        t_print_array = quantities.AdaptingVectorQuantity()
        a_print_array = quantities.AdaptingVectorQuantity()
        e_print_array = quantities.AdaptingVectorQuantity()
        INCL_print_array = quantities.AdaptingVectorQuantity()
        AP_print_array = quantities.AdaptingVectorQuantity()

        N = 0
        while (t<tend):
            t+=dt
            N+=1
            code.evolve_model(t)
            print('t/Myr = ',code.model_time.value_in(units.Myr))

            channel_from_code.copy()
            t_print_array.append(t)
            a_print_array.append(binaries[0].semimajor_axis)
            e_print_array.append(binaries[0].eccentricity | units.none)
            INCL_print_array.append(binaries[0].inclination_relative_to_parent | units.none)
            AP_print_array.append(binaries[0].argument_of_pericenter | units.none)


    def test2(self):
        """
        test 1PN precession in 2-body system
        """
        particles = create_nested_multiple(2,[1.0|units.MSun, 1.0|units.MJupiter], [1.0|units.AU], [0.99], [0.01], [0.01], [0.01])
        binaries = particles[particles.is_binary]

        binaries.include_pairwise_1PN_terms = True

        code = SecularMultiple(redirection='none')
        code.particles.add_particles(particles)

        channel_to_code = particles.new_channel_to(code.particles)
        channel_from_code = code.particles.new_channel_to(particles)

        channel_to_code.copy()

        t = 0.0 | units.Myr
        dt = 1.0e-1 | units.Myr
        tend = 1.0e0 | units.Myr

        t_print_array = quantities.AdaptingVectorQuantity()
        a_print_array = quantities.AdaptingVectorQuantity()
        e_print_array = quantities.AdaptingVectorQuantity()
        AP_print_array = quantities.AdaptingVectorQuantity()

        N = 0
        while (t<tend):
            t+=dt
            N+=1
            code.evolve_model(t)

            channel_from_code.copy()

            print('t/Myr',t.value_in(units.Myr),'omega',binaries[0].argument_of_pericenter | units.none)
            t_print_array.append(t)
            a_print_array.append(binaries[0].semimajor_axis)
            e_print_array.append(binaries[0].eccentricity | units.none)
            AP_print_array.append(binaries[0].argument_of_pericenter | units.none)


        a = binaries[0].semimajor_axis
        e = binaries[0].eccentricity
        M = binaries[0].mass
        rg = constants.G*M/(constants.c**2)
        P = 2.0*numpy.pi*numpy.sqrt(a**3/(constants.G*M))
        t_1PN = (1.0/3.0)*P*(1.0-e**2)*(a/rg)


    def test3(self):
        """
        test GW emission in 2-body system + collision detection
        """
        particles = create_nested_multiple(2,[1.0|units.MSun, 1.0|units.MJupiter], [0.1|units.AU], [0.994], [0.01], [0.01], [0.01])
        binaries = particles[particles.is_binary]
        stars = particles - binaries
        stars.radius = 0.0001 | units.AU

        binaries.check_for_physical_collision_or_orbit_crossing = True
        binaries.include_pairwise_25PN_terms = True

        code = SecularMultiple(redirection='none')
        code.particles.add_particles(particles)

        channel_to_code = particles.new_channel_to(code.particles)
        channel_from_code = code.particles.new_channel_to(particles)

        channel_to_code.copy()

        t_print_array = quantities.AdaptingVectorQuantity()
        a_print_array = quantities.AdaptingVectorQuantity()
        e_print_array = quantities.AdaptingVectorQuantity()
        AP_print_array = quantities.AdaptingVectorQuantity()

        tend = 1.0 | units.Gyr
        N = 0
        t = 0.0 | units.Myr
        dt = 100.0 | units.Myr
        while (t<tend):
            t+=dt
            N+=1
            code.evolve_model(t)
            flag = code.flag

            if flag == 2:
                print('root found')
                break

            channel_from_code.copy()
            t_print_array.append(t)
            a_print_array.append(binaries[0].semimajor_axis)
            e_print_array.append(binaries[0].eccentricity | units.none)
            AP_print_array.append(binaries[0].argument_of_pericenter | units.none)

            print('e',binaries.eccentricity,'a/AU',binaries.semimajor_axis.value_in(units.AU),'rp/AU',(binaries.semimajor_axis*(1.0-binaries.eccentricity)).value_in(units.AU))


    def test4(self):
        """
        test tidal friction in 2-body system
        """

        M = 1.0|units.MJupiter
        R = 40.0|units.RJupiter
        m_per = 1.0|units.MSun
        m = m_per
        mu = m*M/(m+M)
        a0 = 0.1 | units.AU
        e0 = 0.3
        P0 = 2.0*numpy.pi*numpy.sqrt(a0**3/(constants.G*(M+m_per)))
        n0 = 2.0*numpy.pi/P0

        aF = a0*(1.0-e0**2)
        nF = numpy.sqrt( constants.G*(M+m_per)/(aF**3) )

        particles = create_nested_multiple(2, [m_per, M], [a0], [e0], [0.01], [0.01], [0.01])
        binaries = particles[particles.is_binary]
        particles[0].radius = 1.0 | units.RSun
        particles[1].radius = R
        particles[1].spin_vec_x = 0.0 | 1.0/units.day
        particles[1].spin_vec_y = 0.0 | 1.0/units.day
        particles[1].spin_vec_z = 4.0e-2 | 1.0/units.day

        k_L = 0.38
        k_AM = k_L/2.0
        rg = 0.25
        tau = 0.66 | units.s

        I = rg*M*R**2
        alpha = I/(mu*a0**2)
        T = R**3/(constants.G*M*tau)
        t_V = 3.0*(1.0 + 1.0/k_L)*T

        particles[1].tides_method = 2
        particles[1].include_tidal_friction_terms = True
        particles[1].include_tidal_bulges_precession_terms = False
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-8

        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_viscous_time_scale = t_V
        particles[1].tides_gyration_radius = rg

        tD = M*aF**8/(3.0*k_L*tau*constants.G*m_per*(M+m_per)*R**5)
        particles[2].check_for_physical_collision_or_orbit_crossing = True

        code = SecularMultiple(redirection='none')
        code.particles.add_particles(particles)

        code.parameters.relative_tolerance = 1.0e-14

        channel_to_code = particles.new_channel_to(code.particles)
        channel_from_code = code.particles.new_channel_to(particles)

        channel_to_code.copy()

        t = 0.0 | units.Myr
        dt = 1.0e-5 | units.Myr
        tend = 1.0e-4 | units.Myr

        t_print_array = quantities.AdaptingVectorQuantity()
        a_print_array = quantities.AdaptingVectorQuantity()
        n_print_array = quantities.AdaptingVectorQuantity()
        e_print_array = quantities.AdaptingVectorQuantity()
        AP_print_array = quantities.AdaptingVectorQuantity()
        spin_print_array = quantities.AdaptingVectorQuantity()

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            print('flag',code.flag,t,'a/AU',binaries[0].semimajor_axis,'e',binaries[0].eccentricity)

            channel_from_code.copy()

            t_print_array.append(t)
            a_print_array.append(binaries[0].semimajor_axis)
            n_print_array.append(numpy.sqrt(constants.G*(M+m)/(binaries[0].semimajor_axis**3)))
            e_print_array.append(binaries[0].eccentricity | units.none)
            AP_print_array.append(binaries[0].argument_of_pericenter | units.none)
            spin_print_array.append( numpy.sqrt( particles[1].spin_vec_x**2 + particles[1].spin_vec_y**2 + particles[1].spin_vec_z**2) )

            binaries = particles[particles.is_binary]
            bodies = particles - binaries
            print('S_x',bodies.spin_vec_x.value_in(1.0/units.day))
            print('S_y',bodies.spin_vec_y.value_in(1.0/units.day))
            print('S_z',bodies.spin_vec_z.value_in(1.0/units.day))
            print('='*50)

    def test5(self):
        """
        test precession due to tidal bulges
        """

        M = 1.0|units.MJupiter
        R = 1.0|units.RJupiter
        m_per = 1.0|units.MSun
        a0 = 30.0 | units.AU
        e0 = 0.999
        P0 = 2.0*numpy.pi*numpy.sqrt(a0**3/(constants.G*(M+m_per)))
        n0 = 2.0*numpy.pi/P0

        particles = create_nested_multiple(2, [m_per, M], [a0], [e0], [0.01], [0.01], [0.01])
        binaries = particles[particles.is_binary]
        particles[0].radius = 1.0 | units.RSun
        particles[1].radius = R

        k_L = 0.41
        k_AM = k_L/2.0

        particles[1].tides_method = 0
        particles[1].include_tidal_friction_terms = False
        particles[1].include_tidal_bulges_precession_terms = True
        particles[1].include_rotation_precession_terms = False
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-5
        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_gyration_radius = 0.25
        
        code = SecularMultiple(redirection='none')
        code.particles.add_particles(particles)

        channel_to_code = particles.new_channel_to(code.particles)
        channel_from_code = code.particles.new_channel_to(particles)

        channel_to_code.copy()

        t = 0.0 | units.Myr
        dt = 1.0e1 | units.Myr
        tend = 1.0e2 | units.Myr

        t_print_array = quantities.AdaptingVectorQuantity()
        a_print_array = quantities.AdaptingVectorQuantity()
        e_print_array = quantities.AdaptingVectorQuantity()
        AP_print_array = quantities.AdaptingVectorQuantity()

        g_dot_TB = (15.0/8.0)*n0*(8.0+12.0*e0**2+e0**4)*(m_per/M)*k_AM*pow(R/a0,5.0)/pow(1.0-e0**2,5.0)
        t_TB = 2.0*numpy.pi/g_dot_TB

        N=0
        while (t<tend):
            t+=dt
            code.evolve_model(t)
            print('flag',code.flag,t,binaries[0].semimajor_axis,binaries[0].eccentricity)

            channel_from_code.copy()

            t_print_array.append(t)
            a_print_array.append(binaries[0].semimajor_axis)
            e_print_array.append(binaries[0].eccentricity | units.none)
            AP_print_array.append(binaries[0].argument_of_pericenter | units.none)

            N+=1

    def test6(self):
        """
        test precession due to rotation
        """

        M = 1.0|units.MJupiter
        R = 1.5|units.RJupiter
        m_per = 1.0|units.MSun
        a0 = 30.0 | units.AU
        e0 = 0.999
        P0 = 2.0*numpy.pi*numpy.sqrt(a0**3/(constants.G*(M+m_per)))
        n0 = 2.0*numpy.pi/P0

        aF = a0*(1.0-e0**2)
        nF = numpy.sqrt( constants.G*(M+m_per)/(aF**3) )

        particles = create_nested_multiple(2, [m_per, M], [a0], [e0], [1.0e-5], [1.0e-5], [1.0e-5])
        binaries = particles[particles.is_binary]
        particles[0].radius = 1.0 | units.RSun
        particles[1].radius = R
        particles[1].spin_vec_x = 0.0 | 1.0/units.day
        particles[1].spin_vec_y = 0.0 | 1.0/units.day
        Omega_PS0 = n0*(33.0/10.0)*pow(a0/aF,3.0/2.0)
        particles[1].spin_vec_z = Omega_PS0


        k_L = 0.51
        k_AM = k_L/2.0
        rg = 0.25
        particles[1].tides_method = 1
        particles[1].include_tidal_friction_terms = False
        particles[1].include_tidal_bulges_precession_terms = False
        particles[1].include_rotation_precession_terms = True
        particles[1].minimum_eccentricity_for_tidal_precession = 1.0e-5
        particles[1].tides_apsidal_motion_constant = k_AM
        particles[1].tides_gyration_radius = rg

        code = SecularMultiple(redirection='none')
        code.particles.add_particles(particles)

        channel_to_code = particles.new_channel_to(code.particles)
        channel_from_code = code.particles.new_channel_to(particles)

        channel_to_code.copy()

        t = 0.0 | units.Myr
        dt = 1.0e1 | units.Myr
        tend = 1.0e2 | units.Myr

        t_print_array = quantities.AdaptingVectorQuantity()
        a_print_array = quantities.AdaptingVectorQuantity()
        e_print_array = quantities.AdaptingVectorQuantity()
        AP_print_array = quantities.AdaptingVectorQuantity()

        Omega_vec = [particles[1].spin_vec_x,particles[1].spin_vec_y,particles[1].spin_vec_z]
        Omega = numpy.sqrt(Omega_vec[0]**2 + Omega_vec[1]**2 + Omega_vec[2]**2)
        print('Omega/n',Omega/n0)

        g_dot_rot = n0*(1.0 + m_per/M)*k_AM*pow(R/a0,5.0)*(Omega/n0)**2/((1.0-e0**2)**2)
        t_rot = 2.0*numpy.pi/g_dot_rot
        print('t_rot/Myr',t_rot.value_in(units.Myr))

        N=0
        while (t<tend):
            t+=dt
            code.evolve_model(t)
            print('flag',code.flag,t,binaries[0].semimajor_axis,binaries[0].eccentricity)

            channel_from_code.copy()

            t_print_array.append(t)
            a_print_array.append(binaries[0].semimajor_axis)
            e_print_array.append(binaries[0].eccentricity | units.none)
            AP_print_array.append(binaries[0].argument_of_pericenter | units.none)

            N+=1

    def test7(self):
        """
        test collision detection in 3-body system
        """

        particles = create_nested_multiple(3,[1.0|units.MSun, 1.2|units.MSun, 0.9|units.MSun], [1.0|units.AU, 100.0|units.AU], [0.1, 0.5], [0.01, 80.0*numpy.pi/180.0], [0.01, 0.01], [0.01, 0.01])
        binaries = particles[particles.is_binary]
        stars = particles - binaries

        binaries.check_for_physical_collision_or_orbit_crossing = True
        stars.radius = 0.03 | units.AU

        code = SecularMultiple(redirection='none')
        code.particles.add_particles(particles)

        channel_to_code = particles.new_channel_to(code.particles)
        channel_from_code = code.particles.new_channel_to(particles)

        channel_to_code.copy()

        t = 0.0 | units.Myr
        dt = 1.0e-2 | units.Myr
        tend = 1.0e-1 | units.Myr

        t_print_array = quantities.AdaptingVectorQuantity()
        a_print_array = quantities.AdaptingVectorQuantity()
        e_print_array = quantities.AdaptingVectorQuantity()

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.flag
            channel_from_code.copy()

            print('secular_breakdown_has_occurred',binaries.secular_breakdown_has_occurred)
            print('dynamical_instability_has_occurred',binaries.dynamical_instability_has_occurred)
            print('physical_collision_or_orbit_crossing_has_occurred',binaries.physical_collision_or_orbit_crossing_has_occurred)
            print('minimum_periapse_distance_has_occurred',binaries.minimum_periapse_distance_has_occurred)
            print('RLOF_at_pericentre_has_occurred',binaries.RLOF_at_pericentre_has_occurred)

            if flag == 2:
                print('root found')
                break
            print('t_end',code.model_time.value_in(units.Myr))

            t_print_array.append(t)
            a_print_array.append(binaries[0].semimajor_axis)
            e_print_array.append(binaries[0].eccentricity | units.none)

    def test8(self):
        """
        test minimum periapse occurrence
        """

        particles = create_nested_multiple(3,[1.0|units.MSun, 1.2|units.MSun, 0.9|units.MSun], [1.0|units.AU, 100.0|units.AU], [0.1, 0.5], [0.01, 80.0*numpy.pi/180.0], [0.01, 0.01], [0.01, 0.01])
        binaries = particles[particles.is_binary]
        stars = particles - binaries

        binaries.check_for_minimum_periapse_distance = True
        rp_min = 0.1 | units.AU
        binaries.check_for_minimum_periapse_distance_value = rp_min

        code = SecularMultiple(redirection='none')
        code.particles.add_particles(particles)

        channel_to_code = particles.new_channel_to(code.particles)
        channel_from_code = code.particles.new_channel_to(particles)

        channel_to_code.copy()

        t = 0.0 | units.Myr
        dt = 5.0e-3 | units.Myr
        tend = 1.0e-1 | units.Myr

        t_print_array = quantities.AdaptingVectorQuantity()
        a_print_array = quantities.AdaptingVectorQuantity()
        e_print_array = quantities.AdaptingVectorQuantity()

        while (t<tend):
            t+=dt
            code.evolve_model(t)
            flag = code.flag
            channel_from_code.copy()

            print('secular_breakdown_has_occurred',binaries.secular_breakdown_has_occurred)
            print('dynamical_instability_has_occurred',binaries.dynamical_instability_has_occurred)
            print('physical_collision_or_orbit_crossing_has_occurred',binaries.physical_collision_or_orbit_crossing_has_occurred)
            print('minimum_periapse_distance_has_occurred',binaries.minimum_periapse_distance_has_occurred)
            print('RLOF_at_pericentre_has_occurred',binaries.RLOF_at_pericentre_has_occurred)

            if flag == 2:
                print('root found')
                break
            print('t_end',code.model_time.value_in(units.Myr))

            t_print_array.append(t)
            a_print_array.append(binaries[0].semimajor_axis)
            e_print_array.append(binaries[0].eccentricity | units.none)


if __name__ in ('__main__','__plot__'):
    testSecularMultiple = TestSecularMultiple()

    cmd_options = parse_arguments()

    if cmd_options["test"] == 1:
        testSecularMultiple.test1()
    if cmd_options["test"] == 2:
        testSecularMultiple.test2()
    if cmd_options["test"] == 3:
        testSecularMultiple.test3()
    if cmd_options["test"] == 4:
        testSecularMultiple.test4()
    if cmd_options["test"] == 5:
        testSecularMultiple.test5()
    if cmd_options["test"] == 6:
        testSecularMultiple.test6()
    if cmd_options["test"] == 7:
        testSecularMultiple.test7()
    if cmd_options["test"] == 8:
        testSecularMultiple.test8()
