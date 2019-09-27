import os.path
import numpy
from time import sleep

from amuse.test.amusetest import get_path_to_results, TestWithMPI
try:
    from matplotlib import pyplot
    from amuse.plot import scatter, xlabel, ylabel, plot,loglog,semilogx,semilogy, sph_particles_plot
    from amuse.plot import pynbody_column_density_plot, HAS_PYNBODY
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

from amuse.support.exceptions import AmuseException
from amuse.io import write_set_to_file, read_set_from_file
from amuse.community.evtwin.interface import EVtwin
from amuse.community.mesa.interface import MESA
from amuse.community.gadget2.interface import Gadget2
from amuse.community.fi.interface import Fi
from amuse.units import units, generic_unit_system, nbody_system, constants
from amuse.units.generic_unit_converter import ConvertBetweenGenericAndSiUnits
from amuse.datamodel import Particle, Particles, ParticlesSuperset


from amuse.ext.evrard_test import new_evrard_gas_sphere
from amuse.ext.star_to_sph import convert_stellar_model_to_SPH

from amuse.ext.sph_to_star import SPH2StellarModel, convert_SPH_to_stellar_model


class TestSPH2StellarModel(TestWithMPI):
    
    def new_particles(self):
        input_file = os.path.join(get_path_to_results(), "test_sph_to_star_input.hdf5")
        if os.path.exists(input_file):
            return read_set_from_file(input_file, "hdf5")
        
        stellar_evolution = EVtwin()
        stellar_evolution.particles.add_particle(Particle(mass=1.0|units.MSun))
        stellar_evolution.evolve_model(100.0|units.Myr)
        particles = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            500, 
            seed=12345
        ).gas_particles
        stellar_evolution.stop()
        
        hydrodynamics = Gadget2(ConvertBetweenGenericAndSiUnits(1.0|units.MSun, 1.0|units.RSun, 1.0e3|units.s))
        hydrodynamics.gas_particles.add_particles(particles)
        hydrodynamics.evolve_model(1.0|units.s)
        hydrodynamics.gas_particles.copy_values_of_attributes_to(["density", "u", "pressure"], particles)
        hydrodynamics.stop()
        write_set_to_file(particles, input_file, "hdf5")
        return particles
    
    def test1(self):
        print("Test SPH2StellarModel")
        converter = SPH2StellarModel(self.new_particles())
        model = converter.derive_stellar_structure()
        for variable in ['dmass', 'radius', 'rho', 'temperature', 'luminosity', 'X_H', 
                'X_He', 'X_C', 'X_N', 'X_O', 'X_Ne', 'X_Mg', 'X_Si', 'X_Fe']:
            self.assertTrue(hasattr(model, variable))
            self.assertEqual(len(getattr(model, variable)), 500)
    
    def test2(self):
        print("Test SPH2StellarModel result properties")
        converter = SPH2StellarModel(self.new_particles())
        model = converter.derive_stellar_structure() # model is from center to surface
        self.assertTrue(numpy.all(model.radius[:-1] <= model.radius[1:])) # monotonically increasing
        self.assertTrue(numpy.all(model.mass[:-1] <= model.mass[1:])) # monotonically increasing
        n=50 # following properties are not strictly monotonic, because of randomness in the particle distribution
        self.assertTrue(numpy.all(model.temperature[:-n:n] >= model.temperature[n::n]))
        self.assertTrue(numpy.all(model.rho[:-n:n] >= model.rho[n::n]))
        self.assertTrue(numpy.all(model.X_H[:-n:n] <= model.X_H[n::n]))
    
    def test3(self):
        print("Test convert_SPH_to_stellar_model with particles_per_zone")
        model = convert_SPH_to_stellar_model(self.new_particles(), particles_per_zone=50)
        for variable in ['dmass', 'radius', 'rho', 'temperature', 'luminosity', 'X_H', 
                'X_He', 'X_C', 'X_N', 'X_O', 'X_Ne', 'X_Mg', 'X_Si', 'X_Fe']:
            self.assertTrue(hasattr(model, variable))
            self.assertEqual(len(getattr(model, variable)), 10)
        
        self.assertTrue(numpy.all(model.radius[:-1] <= model.radius[1:])) # monotonically increasing
        self.assertTrue(numpy.all(model.mass[:-1] <= model.mass[1:])) # monotonically increasing
        self.assertTrue(numpy.all(model.temperature[:-1] >= model.temperature[1:]))
        self.assertTrue(numpy.all(model.rho[:-1] >= model.rho[1:]))
        self.assertTrue(numpy.all(model.X_H[:-1] <= model.X_H[1:]))
        
        lowres_model = convert_SPH_to_stellar_model(self.new_particles(), particles_per_zone=100)
        self.assertAlmostRelativeEqual(model.dmass.sum(), 1|units.MSun, 3)
        self.assertAlmostRelativeEqual(lowres_model.dmass.sum(), 1|units.MSun, 3)
        self.assertAlmostRelativeEqual(lowres_model.radius, model.radius[1::2], 7)
        self.assertAlmostRelativeEqual(lowres_model.X_H, (model.X_H[0::2]+model.X_H[1::2])/2.0, 7)
    
    def slowtest4(self):
        print("Test convert_SPH_to_stellar_model result in MESA")
        stellar_evolution = self.new_instance(MESA)
        stellar_evolution.particles.add_particle(Particle(mass=1.0|units.MSun)) # reference particle
        stellar_evolution.evolve_model(100.0|units.Myr)
        
        model = convert_SPH_to_stellar_model(self.new_particles()) # model is from center to surface
        stellar_evolution.new_particle_from_model(model, 0.0|units.Myr)
        print(stellar_evolution.particles)
        self.assertAlmostEqual(stellar_evolution.particles.age, [118.18, 0.0] | units.Myr, 1)
        stellar_evolution.evolve_model(200.0|units.Myr)
        print(stellar_evolution.particles)
        self.assertAlmostEqual(stellar_evolution.particles.age, [204.59, 103.02] | units.Myr, 1)
        self.assertAlmostRelativeEqual(stellar_evolution.particles[0].temperature, 
            stellar_evolution.particles[1].temperature, 2)
        self.assertAlmostRelativeEqual(stellar_evolution.particles[0].luminosity, 
            stellar_evolution.particles[1].luminosity, 2)
        stellar_evolution.stop()
    
    def slowtest5(self):
        print("Test convert_SPH_to_stellar_model result in EVtwin")
        stellar_evolution = EVtwin()
        stellar_evolution.parameters.verbosity = True
        stellar_evolution.particles.add_particle(Particle(mass=1.0|units.MSun)) # reference particle
        stellar_evolution.evolve_model(100.0|units.Myr)
        
        model = convert_SPH_to_stellar_model(self.new_particles()) # model is from center to surface
        stellar_evolution.new_particle_from_model(model, 0.0|units.Myr)
        print(stellar_evolution.particles)
        self.assertAlmostEqual(stellar_evolution.particles.age, [100.0, 0.0] | units.Myr, 1)
        stellar_evolution.evolve_model(200.0|units.Myr)
        print(stellar_evolution.particles)
        self.assertAlmostEqual(stellar_evolution.particles.age, [200.0, 100.0] | units.Myr, 1)
        self.assertAlmostRelativeEqual(stellar_evolution.particles[0].temperature, 
            stellar_evolution.particles[1].temperature, 2)
        self.assertAlmostRelativeEqual(stellar_evolution.particles[0].luminosity, 
            stellar_evolution.particles[1].luminosity, 2)
        stellar_evolution.stop()
    


class TestMergerProductToStar(TestWithMPI):
    
    def slowtest1(self):
        stellar_evolution = self.new_instance(MESA)
        stellar_evolution.particles.add_particles(Particles(2, mass=[1.0, 5.0]|units.MSun))
        stellar_evolution.evolve_model(10.0 | units.Myr)
        initial_separation = stellar_evolution.particles.radius.sum()
        
        sph_particles_1 = convert_stellar_model_to_SPH(
            stellar_evolution.particles[0], 
            200, 
            seed=12345
        ).gas_particles
        sph_particles_2 = convert_stellar_model_to_SPH(
            stellar_evolution.particles[1], 
            1000, 
            seed=12345
        ).gas_particles
        stellar_evolution.stop()
        
        initial_speed = 10.0 | units.km / units.s
        sph_particles_2.x  += initial_separation
        sph_particles_1.vx += initial_speed
        all_sph_particles = ParticlesSuperset([sph_particles_1, sph_particles_2])
        all_sph_particles.move_to_center()
        
        t_end = 4.0e3 | units.s
        print("Evolving to:", t_end)
        n_steps = 4
        
        unit_system_converter = ConvertBetweenGenericAndSiUnits(1.0 | units.RSun, 1.0 | units.MSun, t_end)
        hydro_code = Gadget2(unit_system_converter)
        hydro_code.gas_particles.add_particles(all_sph_particles)
        
        pyplot.ion()
        for time in [i*t_end/n_steps for i in range(1, n_steps+1)]:
            hydro_code.evolve_model(time)
            pyplot.close('all')
            pyplot.figure()
            pynbody_column_density_plot(hydro_code.gas_particles, width=10|units.RSun)
            pyplot.draw()
            pyplot.figure()
            loglog(hydro_code.gas_particles.position.lengths_squared(), hydro_code.gas_particles.pressure, 'bo')
            pyplot.draw()
        
        hydro_code.stop()
        sleep(3)
        pyplot.ioff()
