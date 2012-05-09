import os.path
from amuse.test.amusetest import get_path_to_results, TestWithMPI
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, semilogy, xlabel, ylabel, loglog
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.exceptions import AmuseException
from amuse.community.mesa.interface import MESA
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.ext.star_to_sph import *
from amuse.units import units
from amuse.units import generic_unit_system
from amuse.units import nbody_system
from amuse.units import constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.datamodel import Particles
from amuse.datamodel import Particle
from amuse.datamodel import ParticlesSuperset
from amuse.datamodel import Grid
class TestStellarModel2SPH(TestWithMPI):
    
    class StarParticleWithStructure(Particle):
    
        def __init__(self, number_of_species = 3, **keyword_arguments):
            Particle.__init__(self, **keyword_arguments)
            self.particles_set._private.number_of_species = number_of_species
            self.mass = 4.0/3.0 * numpy.pi * (9.0 / 8.0) | units.MSun
            self.radius = 1.0 | units.RSun
        
        def get_number_of_zones(self):
            return 4
        
        def get_number_of_species(self):
            return self.particles_set._private.number_of_species
        
        def get_names_of_species(self, number_of_species = None):
            return (['h1', 'he3', 'he4', 'c12'])[:int(self.particles_set._private.number_of_species)]
        
        def get_masses_of_species(self, number_of_species = None):
            return ([1.0078250, 3.0160293, 4.0026032, 12.0] | units.amu)[:int(self.particles_set._private.number_of_species)]
        
        def get_mass_profile(self, number_of_zones = None):
            return numpy.asarray([2.0, 14.0, 112.0, 448.0]) / sum([2.0, 14.0, 112.0, 448.0])
        
        def get_density_profile(self, number_of_zones = None):
            return [2.0, 2.0, 2.0, 1.0] | units.MSun/units.RSun**3
        
        def get_radius_profile(self, number_of_zones = None):
            return ([1.0, 2.0, 4.0, 8.0] | units.RSun) / 8.0
        
        def get_temperature_profile(self, number_of_zones = None):
            return [1e7, 1e6, 1e5, 1e4] | units.K
        
        def get_luminosity_profile(self, number_of_zones = None):
            return [1.0, 1.0, 1.0, 1.0] | units.LSun
        
        def get_mu_profile(self, number_of_zones = None):
            return [0.8, 0.6, 0.6, 1.3] | units.amu
        
        def get_chemical_abundance_profiles(self, number_of_zones = None, number_of_species = None):
            return numpy.asarray([[0.0, 0.7, 0.7, 0.7], [0.05, 0.01, 0.01, 0.01], [0.95, 0.29, 0.29, 0.29], 
                [0.0, 0.0, 0.0, 0.0]] )[:int(self.particles_set._private.number_of_species)]
    
    def test1(self):
        star = self.StarParticleWithStructure()
        number_of_zones = star.get_number_of_zones()
        delta_mass =  star.mass * star.get_mass_profile()
        outer_radius = star.get_radius_profile()
        inner_radius = [0.0] | units.RSun
        inner_radius.extend(outer_radius[:-1])
        delta_radius_cubed = (outer_radius**3 - inner_radius**3)
        self.assertAlmostEquals(star.get_density_profile() / (delta_mass/(4./3.*numpy.pi*delta_radius_cubed)), 
                                [1]*number_of_zones)
    
    def test2(self):
        star = self.StarParticleWithStructure()
        number_of_sph_particles = 100 # only few particles for test speed-up
        converter = StellarModel2SPH(star, number_of_sph_particles, seed=12345)
        converter.retrieve_stellar_structure()
        self.assertAlmostEqual(converter.specific_internal_energy_profile, 
            [155896.35894, 20786.18119, 2078.61812, 95.93622] | (units.km/units.s)**2, places = 1)
    
    def test3(self):
        print "Test interpolate_hydro_quantities"
        star = self.StarParticleWithStructure()
        number_of_sph_particles = 100 # only few particles for test speed-up
        converter = StellarModel2SPH(star, number_of_sph_particles, seed=12345)
        converter.retrieve_stellar_structure()
        outer_radii = star.get_radius_profile()
        inner_radii = [0.0] | units.RSun
        inner_radii.extend(outer_radii[:-1])
        self.assertEqual(outer_radii, [0.125, 0.25, 0.5, 1.0] | units.RSun)
        self.assertEqual(inner_radii, [0.0, 0.125, 0.25, 0.5] | units.RSun)
        radial_positions = (outer_radii + inner_radii) / 2
        int_specific_internal_energy, int_composition, int_mu = converter.interpolate_internal_energy(radial_positions)
        self.assertEqual( converter.specific_internal_energy_profile, int_specific_internal_energy)
    
    def test4(self):
        print "Test convert_stellar_model_to_SPH"
        star = self.StarParticleWithStructure()
        number_of_sph_particles = 100 # only few particles for test speed-up
        sph_particles = convert_stellar_model_to_SPH(
            star, 
            number_of_sph_particles, 
            seed = 12345
        ).gas_particles
        self.assertEqual(len(sph_particles), number_of_sph_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(), star.mass)
        self.assertAlmostEqual(sph_particles.center_of_mass(), [0,0,0] | units.RSun, 1)
        self.assertIsOfOrder(max(sph_particles.x), star.radius)
        aa = sph_particles.composition.sum(axis=1) - numpy.asarray([1.0]*number_of_sph_particles)
        self.assertAlmostEqual(sph_particles.composition.sum(axis=1), 1.0)
        self.assertTrue(numpy.all( sph_particles.h1  <= 0.7001 ))
        self.assertTrue(numpy.all( sph_particles.he3 <= 0.0501 ))
        self.assertTrue(numpy.all( sph_particles.he4 >= 0.2899 ))
        self.assertTrue(numpy.all( sph_particles.h1[1:]  - sph_particles.h1[:-1]  >= -0.0001 ))
        self.assertTrue(numpy.all( sph_particles.he3[1:] - sph_particles.he3[:-1] <=  0.0001 ))
        self.assertTrue(numpy.all( sph_particles.he4[1:] - sph_particles.he4[:-1] <=  0.0001 ))
    
    def test5(self):
        print "Test evolving created SPH particles in Gadget"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        number_of_sph_particles = 200 # only few particles for test speed-up
        sph_particles = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345
        ).gas_particles
        self.assertEqual(len(sph_particles), number_of_sph_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(), stars.mass)
        self.assertAlmostEqual(sph_particles.center_of_mass(), [0,0,0] | units.RSun, 1)
        self.assertIsOfOrder(max(sph_particles.x), 1 | units.RSun)
        self.assertAlmostEqual(sph_particles.h1, 0.7, places=2)
        stellar_evolution.stop()
        
        time_end = 1.0 | units.s
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, time_end*1000)
        hydrodynamics = Gadget2(unit_converter)
        hydrodynamics.initialize_code()
        hydrodynamics.gas_particles.add_particles(sph_particles)
        hydrodynamics.evolve_model(time_end)
        hydrodynamics.stop()
    
    def slowtest6(self):
        print "Compare composition profile of stellar model to SPH model"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(5.0 | units.Gyr)
        number_of_sph_particles = 1000 # only few particles for test speed-up
        composition = stellar_evolution.particles[0].get_chemical_abundance_profiles()
        outer_radii = stellar_evolution.particles[0].get_radius_profile()
        outer_radii.prepend(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        sph_particles = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345
        ).gas_particles
        temperature = stellar_evolution.particles[0].get_temperature_profile()
        mu          = stellar_evolution.particles[0].get_mu_profile()
        specific_internal_energy = (1.5 * constants.kB * temperature / mu).as_quantity_in(units.J/units.kg)
        stellar_evolution.stop()
        sph_midpoints = sph_particles.position.lengths()
        
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, sph_particles.h1, 
            os.path.join(get_path_to_results(), "star2sph_test_6_composition_h1.png")
        )
        composition_comparison_plot(
            midpoints, composition[2], 
            sph_midpoints, sph_particles.he4, 
            os.path.join(get_path_to_results(), "star2sph_test_6_composition_he4.png")
        )
        
        time_end = 4000.0 | units.s
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, time_end)
        hydrodynamics = Gadget2(unit_converter)
        hydrodynamics.initialize_code()
        hydrodynamics.gas_particles.add_particles(sph_particles)
        hydrodynamics.evolve_model(time_end)
        sph_midpoints = hydrodynamics.gas_particles.position.lengths()
        hydrodynamics.stop()
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, sph_particles.h1, 
            os.path.join(get_path_to_results(), "star2sph_test_6_relaxed_composition_h1.png")
        )
        composition_comparison_plot(
            midpoints, composition[2], 
            sph_midpoints, sph_particles.he4, 
            os.path.join(get_path_to_results(), "star2sph_test_6_relaxed_composition_he4.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, sph_particles.u, 
            os.path.join(get_path_to_results(), "star2sph_test_6_relaxed_composition_u.png")
        )
        
    def slowtest7(self):
        print "Relaxation of stellar evolution model (Gadget2)"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(10.0 | units.Gyr)
        
        composition = stellar_evolution.particles[0].get_chemical_abundance_profiles()
        outer_radii = stellar_evolution.particles[0].get_radius_profile()
        outer_radii.prepend(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        temperature = stellar_evolution.particles[0].get_temperature_profile()
        mu          = stellar_evolution.particles[0].get_mu_profile()
        specific_internal_energy = (1.5 * constants.kB * temperature / mu).as_quantity_in(units.J/units.kg) # units.m**2/units.s**2)
        
        number_of_sph_particles = 10000
        print "Creating initial conditions from a MESA stellar evolution model:"
        print stars.mass[0], "star consisting of", number_of_sph_particles, "particles."
        gas = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345
        ).gas_particles
        stellar_evolution.stop()
        
        t_end = 1.0e4 | units.s
        print "Evolving to:", t_end
        n_steps = 100
        
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, t_end)
        hydro_code = Gadget2(unit_converter)
        hydro_code.gas_particles.add_particles(gas)
        
        times = [] | units.Myr
        kinetic_energies =   [] | units.J
        potential_energies = [] | units.J
        thermal_energies =   [] | units.J
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_code.evolve_model(time)
            times.append(time)
            kinetic_energies.append(   hydro_code.kinetic_energy)
            potential_energies.append( hydro_code.potential_energy)
            thermal_energies.append(   hydro_code.thermal_energy)
        
        sph_midpoints = hydro_code.gas_particles.position.lengths()
        energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_7_after_t1e4_gadget_energy_evolution.png"))
        thermal_energy_plot(times, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_7_after_t1e4_gadget_thermal_energy_evolution.png"))
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, gas.h1, 
            os.path.join(get_path_to_results(), "star2sph_test_7_after_t1e4_gadget_composition_h1.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, gas.u, 
            os.path.join(get_path_to_results(), "star2sph_test_7_after_t1e4_gadget_original_u.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, hydro_code.gas_particles.u, 
            os.path.join(get_path_to_results(), "star2sph_test_7_after_t1e4_gadget_new_u.png")
        )
        hydro_code.stop()
        print "All done!\n"
         
    def slowtest8(self):
        print "Isothermal relaxation of stellar evolution model (Fi)"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        number_of_sph_particles = 10000
        print "Creating initial conditions from a MESA stellar evolution model:"
        print stars.mass[0], "star consisting of", number_of_sph_particles, "particles."
        gas = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345
        ).gas_particles
        stellar_evolution.stop()
        
        t_end = 1000.0 | units.s
        n_steps = 100
        
        print "Evolving to:", t_end
        
        gas.h_smooth = 0.01 | units.RSun
        
        unit_converter = nbody_system.nbody_to_si(1000.0 | units.s, 1.0 | units.RSun)
        hydro_code = Fi(unit_converter)
        hydro_code.parameters.timestep = t_end / n_steps
        hydro_code.parameters.isothermal_flag = True
        hydro_code.parameters.gamma = 1.0
        hydro_code.gas_particles.add_particles(gas)
        
        times = [] | units.Myr
        kinetic_energies =   [] | units.J
        potential_energies = [] | units.J
        thermal_energies =   [] | units.J
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_code.evolve_model(time)
            print "Evolved model to:", time
            times.append(time)
            kinetic_energies.append(   hydro_code.kinetic_energy)
            potential_energies.append( hydro_code.potential_energy)
            thermal_energies.append(   hydro_code.thermal_energy)
        hydro_code.stop()
        energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_8_fi_star2sph.png"))
        print "All done!\n"
    
    def slowtest9(self):
        print "Test convert_stellar_model_to_SPH and relaxation"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(10.0 | units.Gyr)
        composition = stellar_evolution.particles[0].get_chemical_abundance_profiles()
        outer_radii = stellar_evolution.particles[0].get_radius_profile()
        outer_radii.prepend(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        temperature = stellar_evolution.particles[0].get_temperature_profile()
        mu          = stellar_evolution.particles[0].get_mu_profile()
        specific_internal_energy = (1.5 * constants.kB * temperature / mu).as_quantity_in(units.J/units.kg) # units.m**2/units.s**2)
        
        number_of_sph_particles = 1000 # only few particles for test speed-up
        print "Creating initial conditions from a MESA stellar evolution model:"
        print stars.mass[0], "star consisting of", number_of_sph_particles, "particles, with relaxation turned ON."
        sph_particles = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed = 12345,
            do_relax = True
        ).gas_particles
        stellar_evolution.stop()
        sph_midpoints = sph_particles.position.lengths()
        
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, sph_particles.h1, 
            os.path.join(get_path_to_results(), "star2sph_test_9_relax_h1.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, sph_particles.u, 
            os.path.join(get_path_to_results(), "star2sph_test_9_relax_u.png")
        )
    
    def slowtest10(self):
        print "Test convert_stellar_model_to_SPH with relaxation, and subsequently relax with Gadget2"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(10.0 | units.Gyr)
        composition = stellar_evolution.particles[0].get_chemical_abundance_profiles()
        outer_radii = stellar_evolution.particles[0].get_radius_profile()
        outer_radii.prepend(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        temperature = stellar_evolution.particles[0].get_temperature_profile()
        mu          = stellar_evolution.particles[0].get_mu_profile()
        specific_internal_energy = (1.5 * constants.kB * temperature / mu).as_quantity_in(units.J/units.kg) # units.m**2/units.s**2)
        
        number_of_sph_particles = 10000 # only few particles for test speed-up
        print "Creating initial conditions from a MESA stellar evolution model:"
        print stars.mass[0], "star consisting of", number_of_sph_particles, "particles, with relaxation turned ON."
        sph_particles = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed = 12345,
            do_relax = True
        ).gas_particles
        stellar_evolution.stop()
        sph_midpoints = sph_particles.position.lengths()
        
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, sph_particles.h1, 
            os.path.join(get_path_to_results(), "star2sph_test_10_before_h1_new.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, sph_particles.u, 
            os.path.join(get_path_to_results(), "star2sph_test_10_before_u_new.png")
        )
        t_end = 1.0e4 | units.s
        print "Evolving to:", t_end
        n_steps = 100
        
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, t_end)
        hydro_code = Gadget2(unit_converter)
        hydro_code.gas_particles.add_particles(sph_particles)
        
        times = [] | units.Myr
        kinetic_energies =   [] | units.J
        potential_energies = [] | units.J
        thermal_energies =   [] | units.J
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_code.evolve_model(time)
            times.append(time)
            kinetic_energies.append(   hydro_code.kinetic_energy)
            potential_energies.append( hydro_code.potential_energy)
            thermal_energies.append(   hydro_code.thermal_energy)
        
        sph_midpoints = hydro_code.gas_particles.position.lengths()
        energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_10_after_t1e4_gadget_energy_evolution.png"))
        thermal_energy_plot(times, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_10_after_t1e4_gadget_thermal_energy_evolution.png"))
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, sph_particles.h1, 
            os.path.join(get_path_to_results(), "star2sph_test_10_after_t1e4_gadget_h1.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, sph_particles.u, 
            os.path.join(get_path_to_results(), "star2sph_test_10_after_t1e4_gadget_original_u.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, hydro_code.gas_particles.u, 
            os.path.join(get_path_to_results(), "star2sph_test_10_after_t1e4_gadget_new_u.png")
        )
        hydro_code.stop()
        print "All done!\n"
     
    def test11(self):
        print "Test convert_stellar_model_to_SPH with two stars"
        star = self.StarParticleWithStructure()
        number_of_sph_particles = 100 # only few particles for test speed-up
        some_sph_particles = convert_stellar_model_to_SPH(
            star, 
            number_of_sph_particles, 
            seed = 12345
        ).gas_particles
        
        another_star = self.StarParticleWithStructure(number_of_species = 4)
        more_sph_particles = convert_stellar_model_to_SPH(
            another_star, 
            number_of_sph_particles, 
            seed = 12345
        ).gas_particles
        more_sph_particles.x += 100.0 | units.RSun
        
        sph_particles = ParticlesSuperset([some_sph_particles, more_sph_particles])
        string_produced_by_print = sph_particles.__str__()
        self.assertTrue("he3" in string_produced_by_print)
        self.assertFalse("c12" in string_produced_by_print)
        
        self.assertEqual(len(sph_particles), 2 * number_of_sph_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(), 2 * star.mass)
        self.assertIsOfOrder(max(sph_particles.x), (100.0 | units.RSun) + another_star.radius)
        self.assertIsOfOrder(min(sph_particles.x), -star.radius)
        
        self.assertEqual(len(some_sph_particles.composition), number_of_sph_particles)
        self.assertEqual(len(some_sph_particles[0].composition), 3)
        self.assertEqual(len(more_sph_particles[0].composition), 4)
        self.assertRaises(AttributeError, getattr, sph_particles, "composition", 
            expected_message = "Subsets return incompatible quantities for attribute 'composition', attribute cannot be queried from the superset")
        
        
        self.assertAlmostEqual(some_sph_particles.composition.sum(axis=1), [1.0]*number_of_sph_particles)
        self.assertAlmostEqual(more_sph_particles.composition.sum(axis=1), [1.0]*number_of_sph_particles)
        self.assertTrue(numpy.all( sph_particles.h1  <= 0.7001 ))
        self.assertTrue(numpy.all( sph_particles.he3 <= 0.0501 ))
        self.assertTrue(numpy.all( sph_particles.he4 >= 0.2899 ))
        
    def slowtest12(self):
        print "Test merge two stars"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(2)
        stars.mass = [1.0, 1.0] | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(10.0 | units.Gyr)
        
        composition = stellar_evolution.particles[0].get_chemical_abundance_profiles()
        outer_radii = stellar_evolution.particles[0].get_radius_profile()
        outer_radii.prepend(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        temperature = stellar_evolution.particles[0].get_temperature_profile()
        mu          = stellar_evolution.particles[0].get_mu_profile()
        specific_internal_energy = (1.5 * constants.kB * temperature / mu).as_quantity_in(units.J/units.kg)
        
        number_of_sph_particles = 4000
        n_string = "n4e3"
        print "Creating initial conditions from a MESA stellar evolution model:"
        print stars.mass[0], "star consisting of", number_of_sph_particles, "particles."
        sph_particles_1 = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345
        ).gas_particles
        print stars.mass[1], "star consisting of", number_of_sph_particles, "particles."
        sph_particles_2 = convert_stellar_model_to_SPH(
            stellar_evolution.particles[1], 
            number_of_sph_particles, 
            seed=12345
        ).gas_particles
        stellar_evolution.stop()
        initial_separation = 4.0 | units.RSun
        initial_speed = 100.0 | units.km / units.s
        sph_particles_2.x  += initial_separation
        sph_particles_1.vx += initial_speed
        all_sph_particles = ParticlesSuperset([sph_particles_1, sph_particles_2])
        
        t_end = 4.0e4 | units.s
        t_end_string = "t4e4"
        print "Evolving to:", t_end
        n_steps = 100
        
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, t_end)
        hydro_code = Gadget2(unit_converter)
        hydro_code.gas_particles.add_particles(all_sph_particles)
        
        times = [] | units.Myr
        kinetic_energies =   [] | units.J
        potential_energies = [] | units.J
        thermal_energies =   [] | units.J
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_code.evolve_model(time)
            times.append(time)
            kinetic_energies.append(   hydro_code.kinetic_energy)
            potential_energies.append( hydro_code.potential_energy)
            thermal_energies.append(   hydro_code.thermal_energy)
        
        energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_12_merger_"+n_string+"_"+t_end_string+"_energy_evolution.png"))
        thermal_energy_plot(times, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_12_merger_"+n_string+"_"+t_end_string+"_thermal_energy_evolution.png"))
        
        channel = hydro_code.gas_particles.new_channel_to(all_sph_particles)
        channel.copy_attributes(['mass', 'rho', 'x','y','z', 'vx','vy','vz', 'u'])   
        center_of_mass = all_sph_particles.center_of_mass().as_quantity_in(units.RSun)
        center_of_mass_velocity = all_sph_particles.center_of_mass_velocity().as_quantity_in(units.km / units.s)
        print "center_of_mass:", center_of_mass
        print "center_of_mass_velocity:", center_of_mass_velocity
        self.assertIsOfOrder(center_of_mass[0], 0.5 * (initial_separation + t_end * initial_speed))
        self.assertIsOfOrder(center_of_mass_velocity[0], 0.5 * initial_speed)
        all_sph_particles.position -= center_of_mass
        sph_midpoints = all_sph_particles.position.lengths()
        
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, all_sph_particles.h1, 
            os.path.join(get_path_to_results(), "star2sph_test_12_merger_"+n_string+"_"+t_end_string+"_composition_h1.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, all_sph_particles.u, 
            os.path.join(get_path_to_results(), "star2sph_test_12_merger_"+n_string+"_"+t_end_string+"_new_u.png")
        )
        hydro_code.stop()
        print "All done!\n"
        
    def slowtest13(self):
        print "Super giant model in SPH"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        star =  Particle()
        star.mass = 10.0 | units.MSun
        stellar_evolution.initialize_code() 
        se_star = stellar_evolution.particles.add_particle(star)
        stellar_evolution.commit_particles()
        original_outer_radii = se_star.get_radius_profile().as_quantity_in(units.RSun)
        original_density     = se_star.get_density_profile()
        try:
            while True:
                stellar_evolution.evolve_model()
        except AmuseException as ex:
            self.assertEqual(str(ex), "Error when calling 'evolve' of a 'MESA', errorcode is -14, error "
            "is 'Evolve terminated: Maximum number of backups reached.'")
        
        composition = se_star.get_chemical_abundance_profiles()
        density     = se_star.get_density_profile()
        outer_radii = se_star.get_radius_profile()
        outer_radii.prepend(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        temperature = se_star.get_temperature_profile()
        mu          = se_star.get_mu_profile()
        specific_internal_energy = (1.5 * constants.kB * temperature / mu).as_quantity_in(units.J/units.kg)
        
        pyplot.figure(figsize = (5, 5))
        loglog(original_outer_radii, original_density, label = "t = "+str(0|units.Myr))
        loglog(outer_radii[1:], density, label = "t = "+str(se_star.age.as_quantity_in(units.Myr)))
        xlabel('radius')
        ylabel('density')
        pyplot.legend(loc=3)
        figname = os.path.join(get_path_to_results(), "star2sph_test_13_density.png")
        pyplot.savefig(figname)
        print "\nPlot of density profile was saved to: ", figname
        pyplot.close()
        
        number_of_sph_particles = 1000
        print "Creating initial conditions from a MESA stellar evolution model:"
        print star.mass[0], "star consisting of", number_of_sph_particles, "particles."
        gas = convert_stellar_model_to_SPH(
            se_star, 
            number_of_sph_particles
        ).gas_particles
        stellar_evolution.stop()
        
        t_end = 1.0e3 | units.s
        print "Evolving to:", t_end
        n_steps = 100
        
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, t_end)
        hydro_code = Gadget2(unit_converter)
        hydro_code.gas_particles.add_particles(gas)
        
        times = [] | units.Myr
        kinetic_energies =   [] | units.J
        potential_energies = [] | units.J
        thermal_energies =   [] | units.J
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_code.evolve_model(time)
            times.append(time)
            kinetic_energies.append(   hydro_code.kinetic_energy)
            potential_energies.append( hydro_code.potential_energy)
            thermal_energies.append(   hydro_code.thermal_energy)
        
        sph_midpoints = hydro_code.gas_particles.position.lengths()
        energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_13_n1e3_after_t1e3_gadget_energy_evolution.png"))
        thermal_energy_plot(times, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_13_n1e3_after_t1e3_gadget_thermal_energy_evolution.png"))
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, gas.h1, 
            os.path.join(get_path_to_results(), "star2sph_test_13_n1e3_after_t1e3_gadget_composition_h1.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, gas.u, 
            os.path.join(get_path_to_results(), "star2sph_test_13_n1e3_after_t1e3_gadget_original_u.png")
        )
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, hydro_code.gas_particles.u, 
            os.path.join(get_path_to_results(), "star2sph_test_13_n1e3_after_t1e3_gadget_new_u.png")
        )
        hydro_code.stop()
        print "All done!\n"
    
    def slowtest14(self):
        print "SPH model with core"
        # options:
        with_core = True # set to False to do a comparison run without a core (True)
        hydro_code = Gadget2 # Fi -or- Gadget2
        
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        try:
            while True:
                stellar_evolution.evolve_model()
        except AmuseException as ex:
            self.assertEqual(str(ex), "Error when calling 'evolve' of a 'MESA', errorcode is -14, error "
            "is 'Evolve terminated: Maximum number of backups reached.'")
        
        outer_radii = stellar_evolution.particles[0].get_radius_profile()
        outer_radii.prepend(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        temperature = stellar_evolution.particles[0].get_temperature_profile()
        mu          = stellar_evolution.particles[0].get_mu_profile()
        specific_internal_energy = (1.5 * constants.kB * temperature / mu).as_quantity_in(units.J/units.kg) # units.m**2/units.s**2)
        
        number_of_sph_particles = 3000
        print "Creating initial conditions from a MESA stellar evolution model:"
        print stars.mass[0], "star consisting of", number_of_sph_particles, "particles."
        
        
        stellar_model_in_SPH = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed = 12345,
            with_core_particle = with_core
        )
        if len(stellar_model_in_SPH.core_particle):
            print "Created", len(stellar_model_in_SPH.gas_particles), 
            print "SPH particles and one 'core-particle':\n", stellar_model_in_SPH.core_particle
            core_radius = stellar_model_in_SPH.core_radius
        else:
            print "Only SPH particles created."
            core_radius = 1.0 | units.RSun
        print "Setting gravitational smoothing to:", core_radius
        
        stellar_evolution.stop()
        
        t_end = 1.0e2 | units.s
        print "Evolving to:", t_end
        n_steps = 100
        
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, constants.G, t_end)
        hydro_code = hydro_code(unit_converter)
        
        try:
            hydro_code.parameters.timestep = t_end / n_steps
        except Exception as exc:
            if not "parameter is read-only" in str(exc): raise
        
        hydro_code.parameters.epsilon_squared = core_radius**2
        hydro_code.gas_particles.add_particles(stellar_model_in_SPH.gas_particles)
        if len(stellar_model_in_SPH.core_particle):
            hydro_code.dm_particles.add_particles(stellar_model_in_SPH.core_particle)
        
        self.assertAlmostRelativeEqual(stars.mass, hydro_code.particles.total_mass(), places=7)
        
        times = [] | units.s
        kinetic_energies =   [] | units.J
        potential_energies = [] | units.J
        thermal_energies =   [] | units.J
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_code.evolve_model(time)
            times.append(time)
            kinetic_energies.append(   hydro_code.kinetic_energy)
            potential_energies.append( hydro_code.potential_energy)
            thermal_energies.append(   hydro_code.thermal_energy)
        
        sph_midpoints = hydro_code.gas_particles.position.lengths()
        energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_14_after_t1e2_gadget_energy_evolution.png"))
        thermal_energy_plot(times, thermal_energies, 
            os.path.join(get_path_to_results(), "star2sph_test_14_after_t1e2_gadget_thermal_energy_evolution.png"))
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, hydro_code.gas_particles.u, 
            os.path.join(get_path_to_results(), "star2sph_test_14_after_t1e2_gadget_internal_energy.png")
        )
        hydro_code.stop()
        print "All done!\n"
    
    def test15(self):
        print "Test pickling of stellar structure"
        star = self.StarParticleWithStructure()
        test_pickle_file = os.path.join(get_path_to_results(), "test_star_structure.pkl")
        if os.path.exists(test_pickle_file):
            os.remove(test_pickle_file)
        pickle_stellar_model(star, test_pickle_file)
        converter = StellarModel2SPH(None, 100, seed=12345, 
            pickle_file = test_pickle_file)
        converter.unpickle_stellar_structure()
        self.assertEqual(converter.mass, numpy.pi * 1.5 | units.MSun)
        self.assertEqual(converter.radius, 1.0 | units.RSun)
        self.assertEqual(converter.number_of_zones, 4)
        self.assertEqual(converter.number_of_species, 3)
        self.assertEqual(converter.species_names, ['h1', 'he3', 'he4'])
        self.assertEqual(converter.density_profile, [2.0, 2.0, 2.0, 1.0] | units.MSun/units.RSun**3)
        self.assertEqual(converter.radius_profile, [1.0/8.0, 1.0/4.0, 1.0/2.0, 1.0] | units.RSun)
        self.assertEqual(converter.temperature_profile, [1e7, 1e6, 1e5, 1e4] | units.K)
        self.assertEqual(converter.mu_profile, [0.8, 0.6, 0.6, 1.3] | units.amu)
        self.assertEqual(converter.composition_profile, [[0.0, 0.7, 0.7, 0.7], 
            [0.05, 0.01, 0.01, 0.01], [0.95, 0.29, 0.29, 0.29]])
        self.assertAlmostEqual(converter.specific_internal_energy_profile, 
            [155896.35894, 20786.18119, 2078.61812, 95.93622] | (units.km/units.s)**2, places = 1)
        
        self.assertRaises(AmuseWarning, pickle_stellar_model, star, test_pickle_file, expected_message = 
            "Incorrect file name '{0}'; directory must exist and file may not exist".format(test_pickle_file))
        bogus_pickle_file = os.path.join(get_path_to_results(), "bogus.pkl")
        converter = StellarModel2SPH(None, 100, seed=12345, 
            pickle_file = bogus_pickle_file)
        self.assertRaises(AmuseException, converter.unpickle_stellar_structure, expected_message = 
            "Input pickle file '{0}' does not exist".format(bogus_pickle_file))
    
    def test16(self):
        print "Test convert_stellar_model_to_SPH with pickled stellar structure"
        star = self.StarParticleWithStructure()
        test_pickle_file = os.path.join(get_path_to_results(), "test_star_structure.pkl")
        if os.path.exists(test_pickle_file):
            os.remove(test_pickle_file)
        pickle_stellar_model(star, test_pickle_file)
        
        number_of_sph_particles = 100 # only few particles for test speed-up
        sph_particles = convert_stellar_model_to_SPH(
            None, 
            number_of_sph_particles, 
            seed = 12345,
            pickle_file = test_pickle_file
        ).gas_particles
        self.assertEqual(len(sph_particles), number_of_sph_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(), star.mass)
        self.assertAlmostEqual(sph_particles.center_of_mass(), [0,0,0] | units.RSun, 1)
        self.assertIsOfOrder(max(sph_particles.x), star.radius)
        self.assertAlmostEqual(sph_particles.composition.sum(axis=1), [1.0]*number_of_sph_particles)
        self.assertTrue(numpy.all( sph_particles.h1  <= 0.7001 ))
        self.assertTrue(numpy.all( sph_particles.he3 <= 0.0501 ))
        self.assertTrue(numpy.all( sph_particles.he4 >= 0.2899 ))
        self.assertTrue(numpy.all( sph_particles.h1[1:]  - sph_particles.h1[:-1]  >= -0.0001 ))
        self.assertTrue(numpy.all( sph_particles.he3[1:] - sph_particles.he3[:-1] <=  0.0001 ))
        self.assertTrue(numpy.all( sph_particles.he4[1:] - sph_particles.he4[:-1] <=  0.0001 ))
    
    def slowtest17(self):
        print "SPH red super giant model with core"
        # options:
        with_core = True # set to False to do a comparison run without a core (True)
        use_hydro_code = Gadget2 # Fi -or- Gadget2
        use_stellar_evolution_code = MESA # to be implemented as option...
        number_of_sph_particles = 3000
        t_end = 3.0e6 | units.s
        
        # Convert some of the parameters to string, for use in output file names:    
        hydro_code_string = "_" + str(use_hydro_code.__name__)
        n_string = "_n" + ("%1.0e"%(number_of_sph_particles)).replace("+0","").replace("+","")
        t_end_string = "_t" + ("%1.0e"%(t_end.value_in(units.s))).replace("+0","").replace("+","")
        base_plotfile_string = "star2sph_test_17" + n_string + hydro_code_string + t_end_string
        
        stellar_evolution = self.new_instance(MESA, redirection = "none")
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 10.0 | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()

        star_age = [] | units.yr
        star_radius = [] | units.RSun
        try:
            while True:
                stellar_evolution.evolve_model()
                star_age.append(stellar_evolution.particles[0].age)
                star_radius.append(stellar_evolution.particles[0].radius)
        except AmuseException as ex:
            self.assertEqual(str(ex), "Error when calling 'evolve' of a 'MESA', errorcode is -14, error "
            "is 'Evolve terminated: Maximum number of backups reached.'")
        
        radius_evolution_plot(star_age, star_radius, 
            os.path.join(get_path_to_results(), base_plotfile_string + "_radius_evolution.png"))
        self.assertIsOfOrder(stellar_evolution.particles[0].age, # MS lifetime:
            (1.0e10 | units.yr) * (stars.mass.value_in(units.MSun)) ** -2.5)
        
        outer_radii = stellar_evolution.particles[0].get_radius_profile()
        outer_radii.prepend(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        temperature = stellar_evolution.particles[0].get_temperature_profile()
        mu          = stellar_evolution.particles[0].get_mu_profile()
        specific_internal_energy = (1.5 * constants.kB * temperature / mu).as_quantity_in(units.J/units.kg) # units.m**2/units.s**2)
        
        print "Creating initial conditions from a MESA stellar evolution model:"
        print stars.mass[0], "star consisting of", number_of_sph_particles, "particles."
        
        stellar_model_in_SPH = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed = 12345,
            base_grid_options = dict(type = "glass", target_rms = 0.04),
            with_core_particle = with_core
        )
        if len(stellar_model_in_SPH.core_particle):
            print "Created", len(stellar_model_in_SPH.gas_particles), 
            print "SPH particles and one 'core-particle':\n", stellar_model_in_SPH.core_particle
            core_radius = stellar_model_in_SPH.core_radius
        else:
            print "Only SPH particles created."
            core_radius = 1.0 | units.RSun
        print "Setting gravitational smoothing to:", core_radius
        
        t_dyn = (stellar_evolution.particles[0].radius**3 / (2*constants.G*stars.mass[0])).sqrt()
        print "Dynamical timescale:", t_dyn.as_quantity_in(units.yr)
        stellar_evolution.stop()
        
        print "Evolving to:", t_end, "("+str((t_end/t_dyn)), "dynamical timescales)"
        n_steps = 100
        
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, constants.G, t_end)
        hydro_code = use_hydro_code(unit_converter, redirection = "none")
        
        try:
            hydro_code.parameters.timestep = t_end / n_steps
        except Exception as exc:
            if not "parameter is read-only" in str(exc): raise
        
        hydro_code.parameters.epsilon_squared = core_radius**2
        hydro_code.gas_particles.add_particles(stellar_model_in_SPH.gas_particles)
        if len(stellar_model_in_SPH.core_particle):
            hydro_code.dm_particles.add_particles(stellar_model_in_SPH.core_particle)
        
        self.assertAlmostRelativeEqual(stars.mass, hydro_code.particles.total_mass(), places=7)
        
        times = [] | units.s
        kinetic_energies =   [] | units.J
        potential_energies = [] | units.J
        thermal_energies =   [] | units.J
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_code.evolve_model(time)
            times.append(time)
            kinetic_energies.append(   hydro_code.kinetic_energy)
            potential_energies.append( hydro_code.potential_energy)
            thermal_energies.append(   hydro_code.thermal_energy)
        
        sph_midpoints = hydro_code.gas_particles.position.lengths()
        energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
            os.path.join(get_path_to_results(), base_plotfile_string + "_energy_evolution.png"))
        thermal_energy_plot(times, thermal_energies, 
            os.path.join(get_path_to_results(), base_plotfile_string + "_thermal_energy_evolution.png"))
        internal_energy_comparison_plot(
            midpoints, specific_internal_energy, 
            sph_midpoints, hydro_code.gas_particles.u, 
            os.path.join(get_path_to_results(), base_plotfile_string + "_internal_energy.png")
        )
        hydro_code.stop()
        print "All done!\n"
    
    def slowtest18(self):
        print "SPH red super giant model with core (fixed core mass)"
        number_of_sph_particles = 300
        
        stellar_evolution = self.new_instance(MESA, redirection = "none")
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 50.0 | units.MSun
        stellar_evolution.initialize_code() 
        stellar_evolution.particles.add_particles(stars)
        stellar_evolution.commit_particles()
        stellar_evolution.evolve_model(3.927 | units.Myr)
        
        expected_core_radii = [0.37648393, 0.58139942, 3.26189210, 31.89893263] | units.RSun
        for i, tgt_core_mass in enumerate([1.0, 5.0, 25.0, 49.0] | units.MSun):
            stellar_model_in_SPH = convert_stellar_model_to_SPH(
                stellar_evolution.particles[0], 
                number_of_sph_particles, 
                seed = 12345,
                with_core_particle = True,
                target_core_mass = tgt_core_mass
            )
            self.assertAlmostRelativeEqual(stellar_model_in_SPH.core_particle[0].mass, tgt_core_mass, 1)
            self.assertAlmostEqual(stellar_model_in_SPH.core_radius, expected_core_radii[i])
        stellar_evolution.stop()
    
    def test19(self):
        print "Test convert_stellar_model_to_SPH with do_store_composition"
        star = self.StarParticleWithStructure(number_of_species = 4)
        number_of_sph_particles = 100 # only few particles for test speed-up
        
        sph_particles = convert_stellar_model_to_SPH(
            star, 
            number_of_sph_particles, 
            seed = 12345,
            do_store_composition = False
        ).gas_particles
        self.assertEqual(len(sph_particles), number_of_sph_particles)
        self.assertTrue(set(sph_particles.get_attribute_names_defined_in_store()) ==
            set(["mass", "x", "y", "z", "vx", "vy", "vz", "u", "h_smooth"]))
        self.assertTrue(set(sph_particles.get_attribute_names_defined_in_store()).isdisjoint(
            set(["h1", "he3", "he4", "c12"])))
        
        sph_particles = convert_stellar_model_to_SPH(
            star, 
            number_of_sph_particles, 
            seed = 12345,
        ).gas_particles
        self.assertEqual(len(sph_particles), number_of_sph_particles)
        self.assertTrue(set(sph_particles.get_attribute_names_defined_in_store()) >
            set(["mass", "x", "y", "z", "vx", "vy", "vz", "u", "h_smooth"]))
        self.assertTrue(set(sph_particles.get_attribute_names_defined_in_store()) >
            set(["h1", "he3", "he4", "c12"]))
    

def composition_comparison_plot(radii_SE, comp_SE, radii_SPH, comp_SPH, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (7, 5))
    plot(radii_SE.as_quantity_in(units.RSun), comp_SE, 
        label='stellar evolution model')
    plot(radii_SPH, comp_SPH, 'go', label='SPH model')
    xlabel('radius')
    ylabel('mass fraction')
    pyplot.legend()
    pyplot.savefig(figname)
    print "\nPlot of composition profiles was saved to: ", figname
    pyplot.close()

def internal_energy_comparison_plot(radii_SE, u_SE, radii_SPH, u_SPH, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (7, 5))
    semilogy(radii_SE.as_quantity_in(units.RSun), u_SE, 
        label='stellar evolution model')
    semilogy(radii_SPH, u_SPH, 'go', label='SPH model')
    xlabel('radius')
    ylabel('internal energy')
    pyplot.legend()
    pyplot.savefig(figname)
    print "\nPlot of internal energy profiles was saved to: ", figname
    pyplot.close()

def energy_plot(time, E_kin, E_pot, E_therm, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (5, 5))
    plot(time, E_kin.as_quantity_in(units.erg), label='E_kin')
    plot(time, E_pot, label='E_pot')
    plot(time, E_therm, label='E_therm')
    plot(time, E_kin+E_pot+E_therm, label='E_total')
    xlabel('Time')
    ylabel('Energy')
    pyplot.legend(loc=3)
    pyplot.savefig(figname)
    print "\nPlot of energy evolution was saved to: ", figname
    pyplot.close()

def thermal_energy_plot(time, E_therm, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (5, 5))
    plot(time, E_therm.as_quantity_in(units.erg), label='E_therm')
    xlabel('Time')
    ylabel('Energy')
    pyplot.legend(loc=3)
    pyplot.savefig(figname)
    print "\nPlot of thermal energy evolution was saved to: ", figname
    pyplot.close()
    
def radius_evolution_plot(star_age, star_radius, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (5, 5))
    plot(star_age, star_radius)
    xlabel('Time')
    ylabel('Radius')
    pyplot.savefig(figname)
    print "\nPlot of radius evolution was saved to: ", figname
    pyplot.close()
