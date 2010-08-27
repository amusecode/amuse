import os.path
from amuse.test.amusetest import get_path_to_results, TestWithMPI
try:
    from matplotlib import pyplot
    HAS_MATPLOTLIB = True
    from amuse.plot import plot, xlabel, ylabel
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.data.core import Particles, Particle
from amuse.support.units import units, generic_unit_system, nbody_system
from amuse.support.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.legacy.mesa.interface import MESA
from amuse.legacy.gadget2.interface import Gadget2
from amuse.legacy.fi.interface import Fi
from amuse.ext.star_to_sph import *


class TestStellarModel2SPH(TestWithMPI):
    
    class StarParticleWithStructure(Particle):
    
        def __init__(self, **keyword_arguments):
            Particle.__init__(self, **keyword_arguments)
            self.mass = 4.0/3.0 * numpy.pi * (9.0 / 8.0) | units.MSun
            self.radius = 1.0 | units.RSun
        
        def get_number_of_zones(self):
            return 4 | units.none
        
        def get_number_of_species(self):
            return 3 | units.none
        
        def get_names_of_species(self, number_of_species = None):
            return ['h1', 'he3', 'he4']
        
        def get_IDs_of_species(self, number_of_species = None):
            return [2,    5,     6]
        
        def get_masses_of_species(self, number_of_species = None):
            return [1.0078250, 3.0160293, 4.0026032] | units.amu
        
        def get_mass_profile(self, number_of_zones = None):
            return ([448.0, 112.0, 14.0, 2.0] | units.none) / sum([448.0, 112.0, 14.0, 2.0])
        
        def get_density_profile(self, number_of_zones = None):
            return [1.0, 2.0, 2.0, 2.0] | units.MSun/units.RSun**3
        
        def get_radius_profile(self, number_of_zones = None):
            return ([8.0, 4.0, 2.0, 1.0] | units.RSun) / 8.0
        
        def get_temperature_profile(self, number_of_zones = None):
            return [1e4, 1e5, 1e6, 1e7] | units.K
        
        def get_luminosity_profile(self, number_of_zones = None):
            return [1.0, 1.0, 1.0, 1.0] | units.LSun
        
        def get_mu_profile(self, number_of_zones = None):
            return [1.3, 0.6, 0.6, 0.8] | units.amu
        
        def get_chemical_abundance_profiles(self, number_of_zones = None, number_of_species = None):
            return [[0.7, 0.7, 0.7, 0.0], [0.01, 0.01, 0.01, 0.05], [0.29, 0.29, 0.29, 0.95]] | units.none
    
    def test1(self):
        star = self.StarParticleWithStructure()
        number_of_zones = star.get_number_of_zones().number
        delta_mass = star.get_mass_profile() * star.mass
        outer_radius = star.get_radius_profile()
        inner_radius = outer_radius[1:]
        inner_radius.append(0|units.m)
        delta_radius_cubed = (outer_radius**3 - inner_radius**3)
        self.assertAlmostEquals(star.get_density_profile() / (delta_mass/(4./3.*numpy.pi*delta_radius_cubed)), 
                                [1]*number_of_zones|units.none)
    
    def test2(self):
        star = self.StarParticleWithStructure()
        number_of_sph_particles = 100 # only few particles for test speed-up
        converter = StellarModel2SPH(star, number_of_sph_particles, seed=12345, mode = "scaling method")
        converter.retrieve_stellar_structure()
        self.assertAlmostEqual(converter.specific_internal_energy, 
            [95.93622, 2078.61812, 20786.18119, 155896.35894] | (units.km/units.s)**2, places = 1)
    
    def test3(self):
        print "Test interpolate_hydro_quantities with scaling method"
        star = self.StarParticleWithStructure()
        number_of_sph_particles = 100 # only few particles for test speed-up
        converter = StellarModel2SPH(star, number_of_sph_particles, seed=12345, mode = "scaling method")
        converter.retrieve_stellar_structure()
        outer_radii = star.get_radius_profile()
        inner_radii = outer_radii[1:]
        inner_radii.append(0 | units.RSun)
        self.assertEqual(outer_radii, [1.0, 0.5, 0.25, 0.125] | units.RSun)
        self.assertEqual(inner_radii, [0.5, 0.25, 0.125, 0.0] | units.RSun)
        radial_positions = (outer_radii + inner_radii) / 2
        (int_specific_internal_energy,), int_composition = converter.interpolate_hydro_quantities(radial_positions)
        self.assertEqual( converter.specific_internal_energy, int_specific_internal_energy)
    
    def test4(self):
        print "Test interpolate_hydro_quantities with random sampling"
        star = self.StarParticleWithStructure()
        number_of_sph_particles = 100 # only few particles for test speed-up
        converter = StellarModel2SPH(star, number_of_sph_particles, seed=12345, mode = "random sampling")
        converter.retrieve_stellar_structure()
        converter.set_zone_indices_and_interpolation_coeffs()
        radius_profile = star.get_radius_profile()
        radius_profile.append(0|units.m)
        radial_positions = (    converter.delta    * radius_profile[converter.zone_index] + 
                             (1 - converter.delta) * radius_profile[converter.zone_index+1] ).value_in(units.RSun)
        (int_specific_internal_energy,), int_composition = converter.interpolate_hydro_quantities(radial_positions)
        self.assertEqual(len(int_specific_internal_energy), number_of_sph_particles)
        self.assertTrue(numpy.all( int_specific_internal_energy >= min(converter.specific_internal_energy) ))
        self.assertTrue(numpy.all( int_specific_internal_energy <= max(converter.specific_internal_energy) ))
        dtype = [('r', 'float64'), ('index', 'int32')]
        sorted = numpy.sort(numpy.array(zip(radial_positions, range(number_of_sph_particles)), dtype=dtype), order='r')
        sorted_energies = int_specific_internal_energy[sorted['index']]
        self.assertTrue(numpy.all( sorted_energies[1:]  - sorted_energies[:-1]  <= 0.0 | (units.m/units.s)**2 ))
    
    def test5(self):
        print "Test convert_stellar_model_to_SPH with scaling method"
        star = self.StarParticleWithStructure()
        number_of_sph_particles = 100 # only few particles for test speed-up
        sph_particles = convert_stellar_model_to_SPH(
            star, 
            number_of_sph_particles, 
            seed = 12345,
            mode = "scaling method"
        )
        self.assertEqual(len(sph_particles), number_of_sph_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(), star.mass)
        self.assertAlmostEqual(sph_particles.center_of_mass(), [0,0,0] | units.RSun, 1)
        self.assertIsOfOrder(max(sph_particles.x), star.radius)
        self.assertAlmostEqual(sph_particles.composition.sum(axis=1), [1.0]*number_of_sph_particles | units.none)
        self.assertTrue(numpy.all( sph_particles.h1  <= 0.7001 | units.none ))
        self.assertTrue(numpy.all( sph_particles.he3 <= 0.0501 | units.none ))
        self.assertTrue(numpy.all( sph_particles.he4 >= 0.2899 | units.none ))
        self.assertTrue(numpy.all( sph_particles.h1[1:]  - sph_particles.h1[:-1]  >= -0.0001 | units.none ))
        self.assertTrue(numpy.all( sph_particles.he3[1:] - sph_particles.he3[:-1] <=  0.0001 | units.none ))
        self.assertTrue(numpy.all( sph_particles.he4[1:] - sph_particles.he4[:-1] <=  0.0001 | units.none ))
    
    def test6(self):
        print "Test convert_stellar_model_to_SPH with random sampling"
        star = self.StarParticleWithStructure()
        number_of_sph_particles = 100 # only few particles for test speed-up
        sph_particles = convert_stellar_model_to_SPH(
            star, 
            number_of_sph_particles, 
            seed = 12345,
            mode = "random sampling"
        )
        
        self.assertEqual(len(sph_particles), number_of_sph_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(), star.mass)
        self.assertIsOfOrder(max(sph_particles.x), star.radius)
        self.assertAlmostEqual
        (sph_particles.composition.sum(axis=1), [1.0]*number_of_sph_particles | units.none)
        self.assertTrue(numpy.all( sph_particles.h1  <= 0.7001 | units.none ))
        self.assertTrue(numpy.all( sph_particles.he3 <= 0.0501 | units.none ))
        self.assertTrue(numpy.all( sph_particles.he4 >= 0.2899 | units.none ))
        
        r_squared = sph_particles.position.lengths_squared().value_in(units.RSun**2)
        dtype = [('r_squared', 'float64'), ('key', 'uint64')]
        sorted = numpy.sort(numpy.array(zip(r_squared, sph_particles._get_keys()), dtype=dtype), order='r_squared')
        sorted_particles = sph_particles._subset(sorted['key'])
        sorted_particles.r = sorted_particles.position.lengths()
        self.assertTrue(numpy.all( sorted_particles.h1[1:]  - sorted_particles.h1[:-1]  >= -0.0001 | units.none ))
        self.assertTrue(numpy.all( sorted_particles.he3[1:] - sorted_particles.he3[:-1] <=  0.0001 | units.none ))
        self.assertTrue(numpy.all( sorted_particles.he4[1:] - sorted_particles.he4[:-1] <=  0.0001 | units.none ))
    
    def slowtest7(self):
        print "Test relaxation with scaling method"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_module_with_default_parameters() 
        stellar_evolution.setup_particles(stars)
        stellar_evolution.initialize_stars()
        number_of_sph_particles = 1000 # only few particles for test speed-up
        sph_particles = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345
        )
        self.assertEqual(len(sph_particles), number_of_sph_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(), stars.mass)
        self.assertAlmostEqual(sph_particles.center_of_mass(), [0,0,0] | units.RSun, 1)
        self.assertIsOfOrder(max(sph_particles.x), 1 | units.RSun)
        self.assertAlmostEqual(sph_particles.h1, 0.7 | units.none, places=2)
        stellar_evolution.stop()
        
        time_end = 10.0 | units.s
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, time_end*100)
        hydrodynamics = Gadget2(unit_converter)
        hydrodynamics.initialize_code()
        hydrodynamics.gas_particles.add_particles(sph_particles)
        hydrodynamics.evolve_model(time_end)
        hydrodynamics.stop()
    
    def slowtest8(self):
        print "Test relaxation with random sampling"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_module_with_default_parameters() 
        stellar_evolution.setup_particles(stars)
        stellar_evolution.initialize_stars()
        number_of_sph_particles = 1000 # only few particles for test speed-up
        sph_particles = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed = 12345,
            mode = "random sampling"
        )
        self.assertEqual(len(sph_particles), number_of_sph_particles)
        self.assertAlmostEqual(sph_particles.mass.sum(), stars.mass)
        self.assertAlmostEqual(sph_particles.center_of_mass(), [0,0,0] | units.RSun, 1)
        self.assertIsOfOrder(max(sph_particles.x), 1 | units.RSun)
        self.assertAlmostEqual(sph_particles.h1, 0.7 | units.none, places=2)
        stellar_evolution.stop()
        
        time_end = 10.0 | units.s
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, time_end*100)
        hydrodynamics = Gadget2(unit_converter)
        hydrodynamics.initialize_code()
        hydrodynamics.gas_particles.add_particles(sph_particles)
        hydrodynamics.evolve_model(time_end)
        hydrodynamics.stop()
    
    def slowtest9(self):
        print "Compare composition profile of stellar model to SPH model"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_module_with_default_parameters() 
        stellar_evolution.setup_particles(stars)
        stellar_evolution.initialize_stars()
        stellar_evolution.evolve_model(5.0 | units.Gyr)
        number_of_sph_particles = 10000 # only few particles for test speed-up
        composition = stellar_evolution.particles[0].get_chemical_abundance_profiles()
        outer_radii = stellar_evolution.particles[0].get_radius_profile()
        outer_radii.append(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        sph_particles = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345
        )
        stellar_evolution.stop()
        sph_midpoints = sph_particles.position.lengths()
        
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, sph_particles.h1, 
            os.path.join(get_path_to_results(), "test_composition_h1.png")
        )
        composition_comparison_plot(
            midpoints, composition[2], 
            sph_midpoints, sph_particles.he4, 
            os.path.join(get_path_to_results(), "test_composition_he4.png")
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
            os.path.join(get_path_to_results(), "test_relaxed_composition_h1.png")
        )
        composition_comparison_plot(
            midpoints, composition[2], 
            sph_midpoints, sph_particles.he4, 
            os.path.join(get_path_to_results(), "test_relaxed_composition_he4.png")
        )
    
    def slowtest10(self):
        print "Compare composition profile of stellar model to SPH model (random sampling)"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_module_with_default_parameters() 
        stellar_evolution.setup_particles(stars)
        stellar_evolution.initialize_stars()
        stellar_evolution.evolve_model(5.0 | units.Gyr)
        number_of_sph_particles = 10000 # only few particles for test speed-up
        composition = stellar_evolution.particles[0].get_chemical_abundance_profiles()
        outer_radii = stellar_evolution.particles[0].get_radius_profile()
        outer_radii.append(0.0 | units.m)
        midpoints = (outer_radii[:-1] + outer_radii[1:]) / 2
        sph_particles = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345, 
            mode = "random sampling"
        )
        stellar_evolution.stop()
        sph_midpoints = sph_particles.position.lengths()
        
        composition_comparison_plot(
            midpoints, composition[0], 
            sph_midpoints, sph_particles.h1, 
            os.path.join(get_path_to_results(), "test_rs_composition_h1.png")
        )
        composition_comparison_plot(
            midpoints, composition[2], 
            sph_midpoints, sph_particles.he4, 
            os.path.join(get_path_to_results(), "test_rs_composition_he4.png")
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
            os.path.join(get_path_to_results(), "test_rs_relaxed_composition_h1.png")
        )
        composition_comparison_plot(
            midpoints, composition[2], 
            sph_midpoints, sph_particles.he4, 
            os.path.join(get_path_to_results(), "test_rs_relaxed_composition_he4.png")
        )
        
        
    def slowtest11(self):
        print "Relaxation of stellar evolution model (Gadget2)"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_module_with_default_parameters() 
        stellar_evolution.setup_particles(stars)
        stellar_evolution.initialize_stars()
        number_of_sph_particles = 1000
        print "Creating initial conditions from a MESA stellar evolution model:"
        print stars.mass[0], "star consisting of", number_of_sph_particles, "particles."
        gas = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345
        )
        stellar_evolution.stop()
        
        t_end = 1000.0 | units.s
        print "Evolving to:", t_end
        n_steps = 100
        
        gas.radius = 0.01 | units.RSun # Fi uses radius as smoothing length. Like Gadget2, we use 0.01 | generic length
        
        unit_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, t_end)
        hydro_legacy_code = Gadget2(unit_converter)
        try:
            hydro_legacy_code.parameters.timestep = t_end / n_steps
        except Exception as exc:
            if not "parameter is read-only" in str(exc): raise
        hydro_legacy_code.gas_particles.add_particles(gas)
        
        times = [] | units.Myr
        kinetic_energies =   [] | units.J
        potential_energies = [] | units.J
        thermal_energies =   [] | units.J
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_legacy_code.evolve_model(time)
            times.append(time)
            kinetic_energies.append(   hydro_legacy_code.kinetic_energy)
            potential_energies.append( hydro_legacy_code.potential_energy)
            thermal_energies.append(   hydro_legacy_code.thermal_energy)
        hydro_legacy_code.stop()
        energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
            os.path.join(get_path_to_results(), "test_gadget_star2sph.png"))
        print "All done!\n"
         
    def slowtest12(self):
        print "Isothermal relaxation of stellar evolution model (Fi)"
        stellar_evolution = self.new_instance(MESA)
        if stellar_evolution is None:
            print "MESA was not built. Skipping test."
            return
        stars =  Particles(1)
        stars.mass = 1.0 | units.MSun
        stellar_evolution.initialize_module_with_default_parameters() 
        stellar_evolution.setup_particles(stars)
        stellar_evolution.initialize_stars()
        number_of_sph_particles = 10000
        print "Creating initial conditions from a MESA stellar evolution model:"
        print stars.mass[0], "star consisting of", number_of_sph_particles, "particles."
        gas = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            number_of_sph_particles, 
            seed=12345
        )
        stellar_evolution.stop()
        
        t_end = 10.0 | units.s
        n_steps = 100
        
        print "Evolving to:", t_end
        
        gas.radius = 0.01 | units.RSun # Fi uses radius as smoothing length. Like Gadget2, we use 0.01 | generic length
        
        unit_converter = nbody_system.nbody_to_si(1000.0 | units.s, 1.0 | units.RSun)
        hydro_legacy_code = Fi(unit_converter)
        hydro_legacy_code.parameters.timestep = t_end / n_steps
        hydro_legacy_code.parameters.isothermal_flag = True
        hydro_legacy_code.parameters.gamma = 1.0 | units.none
        hydro_legacy_code.gas_particles.add_particles(gas)
        
        times = [] | units.Myr
        kinetic_energies =   [] | units.J
        potential_energies = [] | units.J
        thermal_energies =   [] | units.J
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_legacy_code.evolve_model(time)
            print "Evolved model to:", time
            times.append(time)
            kinetic_energies.append(   hydro_legacy_code.kinetic_energy)
            potential_energies.append( hydro_legacy_code.potential_energy)
            thermal_energies.append(   hydro_legacy_code.thermal_energy)
        hydro_legacy_code.stop()
        energy_plot(times, kinetic_energies, potential_energies, thermal_energies, 
            os.path.join(get_path_to_results(), "test_fi_star2sph.png"))
        print "All done!\n"
    

def composition_comparison_plot(radii_SE, comp_SE, radii_SPH, comp_SPH, figname):
    if not HAS_MATPLOTLIB:
        return
    pyplot.figure(figsize = (7, 5))
    plot(radii_SE.as_quantity_in(units.RSun), comp_SE.value_in(units.none), 
        label='stellar evolution model')
    plot(radii_SPH, comp_SPH.value_in(units.none), 'go', label='SPH model')
    xlabel('radius')
    ylabel('mass fraction')
    pyplot.legend()
    pyplot.savefig(figname)
    print "\nPlot of composition profiles was saved to: ", figname
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
    
