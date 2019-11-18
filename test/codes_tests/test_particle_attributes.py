import random
import numpy.random
import sys

from amuse.test import amusetest
from amuse.units import units
from amuse.units import constants
from amuse.units import nbody_system
from amuse.support.interface import ConvertArgumentsException

from amuse.ic.plummer import new_plummer_sphere
from amuse.ic.salpeter import new_salpeter_mass_distribution_nbody
from amuse.datamodel import Particle, Particles, ParticlesWithUnitsConverted
from amuse.datamodel import particle_attributes

class TestParticlesAttributes(amusetest.TestCase):
   
    def test3(self):
        print("Test new_particle_from_cluster_core - nbody units")
        numpy.random.seed(123)
        plummer = new_plummer_sphere(10000)
        result = plummer.new_particle_from_cluster_core(density_weighting_power=1, reuse_hop=True)
        self.assertTrue(isinstance(result, Particle))
        
        # Casertano & Hut (1985, ApJ, 298, 80):  density weighted core radius = 0.6791 * r_plummer
        plummer_radius = 3 * constants.pi / 16.0 | nbody_system.length
        self.assertAlmostRelativeEqual(result.radius, 0.6791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [0.0, 0.0, 0.0] | nbody_system.length, 1)
        self.assertAlmostEqual(result.velocity, [0.0, 0.0, 0.0] | nbody_system.speed, 1)
        self.assertAlmostEqual(result.density, 3.55015420914 | nbody_system.density)
        
        plummer.vx += 42 | nbody_system.speed
        plummer.vy += (1 - plummer.y / abs(plummer.y).amax()) * (13|nbody_system.speed)
        plummer.position *= 0.1
        plummer.position += [1.0, 2.0, 3.0] | nbody_system.length
        result = plummer.new_particle_from_cluster_core(density_weighting_power=1, reuse_hop=False)
        self.assertTrue(isinstance(result, Particle))
        
        self.assertAlmostRelativeEqual(result.radius, 0.06791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [1.0, 2.0, 3.0] | nbody_system.length, 1)
        self.assertAlmostEqual(result.velocity, [42.0, 13.0, 0.0] | nbody_system.speed, 1)
        self.assertAlmostRelativeEqual(result.density, 3.55015420914e3 | nbody_system.density, 4)
    
    def test4(self):
        print("Test new_particle_from_cluster_core - SI units")
        numpy.random.seed(123)
        converter = nbody_system.nbody_to_si(1000.0|units.MSun, 1.0 | units.parsec)
        plummer = new_plummer_sphere(10000, convert_nbody=converter)
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, density_weighting_power=1, reuse_hop=True)
        self.assertTrue(isinstance(result, Particle))
        
        # Casertano & Hut (1985, ApJ, 298, 80):  density weighted core radius = 0.6791 * r_plummer
        plummer_radius = 3 * constants.pi / 16.0 | units.parsec
        self.assertAlmostRelativeEqual(result.radius, 0.6791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [0.0, 0.0, 0.0] | units.parsec, 1)
        self.assertAlmostEqual(result.velocity, [0.0, 0.0, 0.0] | units.km / units.s, 1)
        self.assertAlmostEqual(result.density, 3.55015420914e3 | units.MSun * units.parsec**-3)
        
        plummer.vx += 42 | units.km / units.s
        plummer.vy += (1 - plummer.y / abs(plummer.y).amax()) * (13|units.km / units.s)
        plummer.position *= 0.1
        plummer.position += [1.0, 2.0, 3.0] | units.parsec
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, density_weighting_power=1, reuse_hop=False)
        self.assertTrue(isinstance(result, Particle))
        
        self.assertAlmostRelativeEqual(result.radius, 0.06791 * plummer_radius, 2)
        self.assertAlmostEqual(result.position, [1.0, 2.0, 3.0] | units.parsec, 1)
        self.assertAlmostEqual(result.velocity, [42.0, 13.0, 0.0] | units.km / units.s, 1)
        self.assertAlmostRelativeEqual(result.density, 3.55015420914e6 | units.MSun * units.parsec**-3, 4)
    
    def test5(self):
        print("Test new_particle_from_cluster_core - reuse_hop or not")
        converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0 | units.parsec)
        plummer = new_plummer_sphere(100, convert_nbody=converter)
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, reuse_hop=True)
        
        # Hop wasn't stopped, will use same Hop instance:
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, reuse_hop=True)
        
        nbody_plummer = new_plummer_sphere(100)
        # Hop wasn't stopped, unit_converters don't match:
        self.assertRaises(AttributeError, nbody_plummer.new_particle_from_cluster_core, 
            expected_message="Cannot combine units from different systems: m and length")
        
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, reuse_hop=False)
        
        # Hop was stopped, new instance will be made with supplied unit_converter (None in this case):
        result = nbody_plummer.new_particle_from_cluster_core(reuse_hop=True)
        
        self.assertRaises(ConvertArgumentsException, plummer.new_particle_from_cluster_core, unit_converter=converter,#,
            expected_message="error while converting parameter 'mass', error: Cannot express kg in mass, the units do not have the same bases")
        
        result = nbody_plummer.new_particle_from_cluster_core(reuse_hop=False)
        result = plummer.new_particle_from_cluster_core(unit_converter=converter, reuse_hop=False)
    
    def test6(self):
        print("Test all particle attributes using Hop - each different function creates its own instance of Hop")
        numpy.random.seed(123)
        nbody_plummer = new_plummer_sphere(100)
        nbody_plummer.mass = new_salpeter_mass_distribution_nbody(100)
        
        # Each different function creates its own instance of Hop
        result = nbody_plummer.new_particle_from_cluster_core(reuse_hop=True)
        result = nbody_plummer.bound_subset(G=nbody_system.G, reuse_hop=True)
        result = nbody_plummer.mass_segregation_Gini_coefficient(reuse_hop=True)
        result = nbody_plummer.LagrangianRadii(reuse_hop=True)
        result = nbody_plummer.densitycentre_coreradius_coredens(reuse_hop=True)
        
        converter = nbody_system.nbody_to_si(1.0|units.MSun, 1.0 | units.parsec)
        si_plummer = ParticlesWithUnitsConverted(nbody_plummer, converter.as_converter_from_si_to_nbody())
        functions_using_hop = [
            si_plummer.new_particle_from_cluster_core, 
            si_plummer.bound_subset, 
            si_plummer.mass_segregation_Gini_coefficient, 
            si_plummer.LagrangianRadii, 
            si_plummer.densitycentre_coreradius_coredens
        ]
        
        # Each fails since the Hop instance it tries to reuse has a different unit_converter
        for function_using_hop in functions_using_hop:
            self.assertRaises(Exception, function_using_hop, unit_converter=converter)
        #        expected_message="error while converting parameter 'mass', error: Cannot express kg in mass, the units do not have the same bases (note: check whether Hop needs a converter here)")
        
        # Close all Hop instances:
        nbody_results = []
        nbody_results.append(nbody_plummer.new_particle_from_cluster_core(reuse_hop=False))
        nbody_results.append(nbody_plummer.bound_subset(G=nbody_system.G, reuse_hop=False))
        nbody_results.append(nbody_plummer.mass_segregation_Gini_coefficient(reuse_hop=False))
        nbody_results.append(nbody_plummer.LagrangianRadii(reuse_hop=False))
        nbody_results.append(nbody_plummer.densitycentre_coreradius_coredens(reuse_hop=False))
        
        # Now it works, because the Hop instances were closed, and new ones will be instantiated
        si_results = []
        for function_using_hop in functions_using_hop:
            si_results.append(function_using_hop(unit_converter=converter))
        
        convert = converter.as_converter_from_si_to_nbody()
        self.assertAlmostRelativeEqual(si_results[0].position, 
            ParticlesWithUnitsConverted(nbody_results[0].as_set(), convert)[0].position)
        self.assertAlmostRelativeEqual(si_results[1].position, 
            ParticlesWithUnitsConverted(nbody_results[1], convert).position)
        self.assertAlmostRelativeEqual(si_results[2], nbody_results[2], places=10)
        self.assertAlmostRelativeEqual(si_results[3][0], convert.from_target_to_source(nbody_results[3][0]), places=10)
        for in_si, in_nbody in zip(si_results[4], nbody_results[4]):
            self.assertAlmostRelativeEqual(in_si, convert.from_target_to_source(in_nbody), places=10)
    
    
