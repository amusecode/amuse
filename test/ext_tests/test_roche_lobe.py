import os.path
import math
from amuse.test.amusetest import get_path_to_results, TestWithMPI
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, semilogy, xlabel, ylabel, loglog
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.data.core import Particles, Particle, ParticlesSuperset
from amuse.support.units import units, generic_unit_system, nbody_system, constants
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.support.exceptions import AmuseException
from amuse.community.mesa.interface import MESA
from amuse.community.hermite0.interface import Hermite
from amuse.ext.roche_lobe import *

turn_tests_into_slowtests = False
# False True

class TestRocheLobeOverflow(TestWithMPI):
    
    def test1(self):
        print "Testing the stand-alone determine_RLOF_mass_excess and do_roche_lobe_overflow functions"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            self.skip("MESA was not built. Skipping test.")
            return
        stellar_evolution.initialize_code()
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.commit_parameters() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        star = stellar_evolution.particles[0]
        
        mass_profile_plot(star, os.path.join(get_path_to_results(), "lmx_binary_test_1_cumulative_mass.png"))
        
        current_radius = star.radius
        d_mass, (current_mass, r_profile, mass_transfer_flag) = \
            determine_RLOF_mass_excess(star, current_radius)
        self.assertEqual(current_mass, 1.0 | units.MSun)
        self.assertAlmostEqual(d_mass, 0.0 | units.MSun)
        self.assertAlmostEqual(r_profile[-1], current_radius, places=2)
        d_mass, (current_mass, r_profile, mass_transfer_flag) = \
            determine_RLOF_mass_excess(star, 2*current_radius)
        self.assertFalse(mass_transfer_flag)
        self.assertAlmostEqual(d_mass, 0.0 | units.MSun)
        roche_radii = [0.5, 0.3, 0.2] * current_radius
        expected_mass_excess = [0.14, 0.49, 0.77] | units.MSun
        for r_roche, expected in zip(roche_radii, expected_mass_excess):
            d_mass, (current_mass, r_profile, mass_transfer_flag) = \
                determine_RLOF_mass_excess(star, r_roche)
            self.assertTrue(mass_transfer_flag)
            self.assertAlmostRelativeEqual(d_mass, expected, places=2)
        
        do_roche_lobe_overflow(star, r_roche)
        self.assertEqual(star.mass, current_mass - d_mass)
        
        if turn_tests_into_slowtests:
            print "1 - original radius:  ", current_radius
            print "2 - radius after RLOF:", r_roche 
            stellar_evolution.evolve_model()
            print "3 - radius after recovery:", star.radius
            self.assertEqual(star.mass, current_mass - d_mass)
            print current_mass, "MESA star succesfully recovered after ripping off", expected
            print "Number of backup steps taken by MESA:", star.get_number_of_backups_in_a_row().number
    
    def test2(self):
        print "Testing the RocheLobeOverflow class"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            self.skip("MESA was not built. Skipping test.")
            return
        stellar_evolution.initialize_code()
        stars =  Particles(3)
        stars.mass = [1.0, 2.0, 3.0] | units.MSun
        stellar_evolution.commit_parameters() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        first  = stars[0:1].get_intersecting_subset_in(stellar_evolution.particles)[0]
        second = stars[1:2].get_intersecting_subset_in(stellar_evolution.particles)[0]
        third  = stars[2:].get_intersecting_subset_in(stellar_evolution.particles)[0]
        
        rlof = RocheLobeOverflow()
        self.assertEqual(rlof.accretion_efficiency, 1.0)
        self.assertEqual(len(rlof.particles), 0)
        self.assertEqual(len(rlof.companions), 0)
        
        rlof.add_particle(first, 10.0 | units.RSun)
        rlof.add_particle(second, 20.0 | units.RSun, companion = third)
        self.assertEqual(len(rlof.particles), 2)
        self.assertEqual(rlof.particles, [first, second])
        self.assertEqual(rlof.companions, [None, third])
        self.assertEqual(rlof.overflow_radii, [10.0, 20.0] | units.RSun)
        
        rlof.set_roche_radii([first, second], [100.0, 200.0] | units.RSun)
        self.assertEqual(rlof.overflow_radii, [100.0, 200.0] | units.RSun)
        self.assertRaises(AmuseException, rlof.set_roche_radii, [third], [30.0] | units.RSun, 
            expected_message="A particle (key: {0}) wasn't found.".format(third.key))
        
    
    def test3(self):
        print "Testing RocheLobeOverflow with primaries only"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            self.skip("MESA was not built. Skipping test.")
            return
        stellar_evolution.initialize_code()
        stars =  Particles(3)
        stars.mass = [1.0, 2.0, 3.0] | units.MSun
        stellar_evolution.commit_parameters() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        
        instance = RocheLobeOverflow()
        radii = stellar_evolution.particles.radius
        instance.add_particles(stellar_evolution.particles, 0.9*radii)
        mass_lost = instance.do_roche_lobe_overflow()
        self.assertAlmostEqual(stellar_evolution.particles.mass, stars.mass - mass_lost)
        if turn_tests_into_slowtests:
            stellar_evolution.evolve_model()
            self.assertAlmostEqual(stellar_evolution.particles.mass, stars.mass - mass_lost)
            print stars.mass, "MESA stars succesfully recovered after ripping off", mass_lost
            print "Number of backup steps taken by MESA:", stellar_evolution.particles.get_number_of_backups_in_a_row().number
    
    def test4(self):
        print "Testing RocheLobeOverflow with companions"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            self.skip("MESA was not built. Skipping test.")
            return
        stellar_evolution.initialize_code()
        stars =  Particles(4)
        primaries = stars[:2]
        companions = stars[2:]
        primaries.mass = [3.0, 4.0] | units.MSun
        companions.mass = [1.0, 1.0] | units.MSun
        stellar_evolution.commit_parameters() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        se_primaries = primaries.get_intersecting_subset_in(stellar_evolution.particles)
        se_companions = companions.get_intersecting_subset_in(stellar_evolution.particles)
        
        accretion_efficiency = 0.1
        instance = RocheLobeOverflow(accretion_efficiency = accretion_efficiency)
        radii = se_primaries.radius
        instance.add_particles(se_primaries, 0.8*radii, companions = se_companions)
        mass_lost = instance.do_roche_lobe_overflow()
        self.assertAlmostEqual(se_primaries.mass, primaries.mass - mass_lost)
        self.assertAlmostEqual(se_companions.mass, companions.mass + mass_lost * accretion_efficiency)
        if turn_tests_into_slowtests:
            stellar_evolution.evolve_model()
            self.assertAlmostEqual(se_primaries.mass, primaries.mass - mass_lost)
            self.assertAlmostEqual(se_companions.mass, companions.mass + mass_lost * accretion_efficiency)
            print primaries.mass, "MESA stars succesfully recovered after ripping off", mass_lost
            print companions.mass, "MESA stars succesfully recovered after accreting", mass_lost * accretion_efficiency
            print "Number of backup steps taken by MESA:", stellar_evolution.particles.get_number_of_backups_in_a_row().number
    
    def test5(self):
        print "Testing RocheLobeOverflow with variable roche-radii"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            self.skip("MESA was not built. Skipping test.")
            return
        stellar_evolution.initialize_code()
        stars =  Particles(3)
        stars.mass = [1.0, 2.0, 3.0] | units.MSun
        stellar_evolution.commit_parameters() 
        stellar_evolution.particles.add_particles(stars[::-1])
        stellar_evolution.commit_particles()
        se_stars = stars.get_intersecting_subset_in(stellar_evolution.particles)
        # we added the stars to MESA in reverse order, so the orders of stars and stellar_evolution.particles don't match:
        self.assertEqual(stellar_evolution.particles.mass, stars.mass[::-1])
        # 'se_stars' are MESA particles, stored in the same order as 'stars'
        self.assertEqual(se_stars.mass, stars.mass)
        # Concluding, we can use 'se_stars' as MESA particle set, without the need to worry whether they are in the right order
        first  = se_stars[0]
        second = se_stars[1]
        third  = se_stars[2]
        
        rlof = RocheLobeOverflow()
        radii = se_stars.radius
        rlof.add_particles(se_stars, 1.1*radii)
        mass_lost = rlof.do_roche_lobe_overflow()
        self.assertEqual(mass_lost, [0, 0, 0] | units.MSun)
        
        rlof.set_roche_radii([third, first], 0.9*radii[::-2])
        mass_lost = rlof.do_roche_lobe_overflow()
        self.assertAlmostEqual(mass_lost, [0.00218690, 0.0, 0.00004012] | units.MSun)
        
        rlof.set_roche_radii([second, third], 0.85*radii[1:])
        mass_lost_2 = rlof.do_roche_lobe_overflow()
        self.assertAlmostEqual(mass_lost_2, [0.0, 0.00014400, 0.00000335] | units.MSun)
        
        self.assertAlmostEqual(se_stars.mass, stars.mass - mass_lost - mass_lost_2)
        
        if turn_tests_into_slowtests:
            stellar_evolution.evolve_model()
            self.assertAlmostEqual(se_stars.mass, stars.mass - mass_lost - mass_lost_2)
            print stars.mass, "MESA stars succesfully recovered after ripping off", mass_lost + mass_lost_2
            print "Number of backup steps taken by MESA:", [star.get_number_of_backups_in_a_row().number for star in se_stars]
    
    def test6(self):
        print "Testing RocheLobeOverflow with companions in a dynamics code"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            self.skip("MESA was not built. Skipping test.")
            return
        stellar_evolution.initialize_code()
        stars =  Particles(4)
        stars.mass = [1.0, 3.0, 1.0, 3.0] | units.MSun
        stellar_evolution.commit_parameters() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        se_stars = stars.get_intersecting_subset_in(stellar_evolution.particles)
        
        # two binaries: 1 circular, 1 elliptical with periastron within roche radius
        apastron = 40.0 | units.RSun
        mass_binary = stars.mass[0] + stars.mass[1]
        q = (stars.mass[1] / stars.mass[0]).value_in(units.none) # mass ratio
        Eggleton_roche_estimate = 0.49 * q**(2/3.0) / (0.6 * q**(2/3.0) + math.log(1 + q**(1/3.0))) | units.none
        periastron = 0.9 * se_stars[1].radius / Eggleton_roche_estimate
        
        v_outer = 0.5 * (constants.G * 2 * mass_binary / (2000.0 | units.RSun)).sqrt()
        v_circle = 0.25 * (constants.G * mass_binary / apastron).sqrt()
        v_ellipse = 0.25 * (constants.G * mass_binary * (2/apastron - 2/(apastron + periastron))).sqrt()
        stars.position = [[1030.0, 0, 0], [990.0, 0, 0], [-970.0, 0, 0], [-1010.0, 0, 0]] | units.RSun
        stars.vx = [0, 0, 0, 0] | units.km / units.s
        stars.vy = [3 * v_ellipse, -v_ellipse, 3 * v_circle, -v_circle]
        stars.vz = [1, 1, -1, -1] * v_outer
        stars.radius = 0.0 | units.RSun
        
        converter = nbody_system.nbody_to_si(2 * mass_binary, 100.0 | units.RSun)
        dynamics = Hermite(converter)
        dynamics.initialize_code()
        dynamics.parameters.epsilon_squared = 0.0 | units.RSun ** 2
        dynamics.particles.add_particles(stars)
        dynamics.commit_particles()
    
        rlof = RocheLobeOverflow(dynamics_code = dynamics)
        rlof.add_particles([se_stars[1], se_stars[3]], [0.0, 0.0] | units.RSun, companions = [se_stars[0], se_stars[2]])
        self.assertEqual(rlof.overflow_radii, [0.0, 0.0] | units.RSun)
        rlof.update_roche_radii()
        self.assertAlmostEqual(rlof.overflow_radii, Eggleton_roche_estimate * [1, 1] * apastron)
        mass_lost = rlof.do_roche_lobe_overflow()
        self.assertEqual(mass_lost, [0.0, 0.0] | units.MSun)
        
        half_period = math.pi * (((apastron + periastron)/2)**3 / (constants.G * mass_binary)).sqrt()
        dynamics.evolve_model(half_period)
        
        separations = (dynamics.particles.position[::2]-dynamics.particles.position[1::2]).lengths().as_quantity_in(units.RSun)
        self.assertAlmostEqual(separations[0], periastron, 2)
        self.assertAlmostEqual(separations[1], apastron, 2)
        rlof.update_roche_radii()
        self.assertAlmostEqual(rlof.overflow_radii, Eggleton_roche_estimate * separations)
        mass_lost = rlof.do_roche_lobe_overflow()
        self.assertAlmostEqual(mass_lost, [0.00004012, 0.0] | units.MSun)


class TestLowMassXrayBinary(TestWithMPI):
    
    def test1(self):
        
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            self.skip("MESA was not built. Skipping test.")
            return
        stellar_evolution.stop()
        
        lmxb = LowMassXrayBinary(5 | units.MSun, 1.1 | units.MSun)
        print lmxb.secondary
        lmxb.initialize(1 | units.Gyr)
        print lmxb.secondary
        self.assertTrue(lmxb.secondary.age >= 1 | units.Gyr)
        next_age = lmxb.secondary.age + lmxb.secondary.time_step
        primary_mass, secondary_mass, age = lmxb.roche_lobe_overflow_evolve(1.5 * lmxb.secondary.radius)
        self.assertAlmostRelativeEqual(age, next_age)
        self.assertEqual(primary_mass, 5 | units.MSun)
        self.assertAlmostEqual(secondary_mass, 1.1 | units.MSun)
        
    def test2(self):
        
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            self.skip("MESA was not built. Skipping test.")
            return
        stellar_evolution.stop()
        
        lmxb = LowMassXrayBinary(5 | units.MSun, 1.1 | units.MSun)
        lmxb.initialize(1 | units.Gyr)
        next_age = lmxb.secondary.age + lmxb.secondary.time_step
        primary_mass, secondary_mass, age = lmxb.roche_lobe_overflow_evolve(1.5 * lmxb.secondary.radius)
        self.assertAlmostRelativeEqual(age, next_age)
        self.assertEqual(primary_mass, 5 | units.MSun)
        self.assertAlmostEqual(secondary_mass, 1.1 | units.MSun)
    


def mass_profile_plot(star, figname):
    if not HAS_MATPLOTLIB:
        return
    total_mass = star.mass
    radius_profile = star.get_radius_profile().as_quantity_in(units.RSun)
    mass_profile   = star.get_mass_profile() * total_mass
    
    # Reverse cumulative mass: m[i] = all mass OUTside of r[i-1]
    rev_cumul_mass = [] | units.MSun
    tmp = 0 | units.MSun
    for mass_i in mass_profile:
        rev_cumul_mass.append(total_mass - tmp)
        tmp += mass_i
        
    pyplot.figure(figsize = (5, 5))
    plot(radius_profile, rev_cumul_mass)
    xlabel('radius')
    ylabel('mass')
    pyplot.savefig(figname)
    print "\nPlot of mass profile was saved to: ", figname
    pyplot.close()

