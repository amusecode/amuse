import os
import sys
import numpy
import math

from amuse.test.amusetest import TestWithMPI
from amuse.community.mocassin.interface import MocassinInterface, Mocassin



from mpi4py import MPI
from amuse.units import generic_unit_system
from amuse.units import units
from amuse import datamodel

class TestMocassinInterface(TestWithMPI):
    
    def test0(self):
        instance=self.new_instance(MocassinInterface)
        instance.initialize_code()
        instance.stop()
        
    def test1(self):
        instance=self.new_instance(MocassinInterface)
        instance.initialize_code()
        #instance.redirect_outputs_to("moc-out.txt", "moc-err.txt")
        
        instance.setup_mesh(11,11,11, 100,100,100)
        instance.setup_abundancies()
        print instance.get_default_input_directory()
        instance.set_input_directory(instance.get_default_input_directory())
        instance.set_constant_hydrogen_density(900.0)
        instance.commit_parameters()
        indices_x = list(range(1,12,1))
        x,y,z,error = instance.get_position_of_index(indices_x,[1]*len(indices_x), [1]*len(indices_x))
    
        self.assertEquals(error, 0)
        for index, expected_x in enumerate(list(range(-100, 120, 20))):
            self.assertAlmostRelativeEqual(x[index], expected_x, 6)
            self.assertAlmostRelativeEqual(y[index], -100)
            self.assertAlmostRelativeEqual(z[index], -100)
        
        
        
        instance.stop()
    
    def test2(self):
        instance=self.new_instance(MocassinInterface) #, debugger = "ddd")
        instance.initialize_code()
        instance.set_symmetricXYZ(True)
        
      
        instance.setup_mesh(13,13,13,0.95E+19,0.95E+19,0.95E+19)
        instance.setup_abundancies()
        instance.set_input_directory(instance.get_default_input_directory())
        instance.set_constant_hydrogen_density(100.0)
        instance.set_initial_nebular_temperature(6000.0)
        instance.set_maximum_number_of_monte_carlo_iterations(20)
        instance.set_minimum_convergence_level(100)
        instance.set_total_number_of_photons(10000000)
        instance.set_total_number_of_points_in_frequency_mesh(600)
        instance.set_high_limit_of_the_frequency_mesh(15)
        instance.set_low_limit_of_the_frequency_mesh(1.001e-5)
        instance.set_convergence_limit(0.09)
        instance.set_number_of_ionisation_stages(6)
        instance.setup_auto_convergence(0.2, 2.0, 1000000000)
        instance.commit_parameters()
        instance.define_stars(0.0, 0.0, 0.0, 20000, 6003.6396)
        instance.commit_particles()
        instance.commit_grid()
        
        x, error = instance.get_number_of_elements_used()
        self.assertEquals(0, error)
        self.assertEquals(x, 7)
        indices_x = list(range(1,12,1))
        is_active,error = instance.get_grid_active(indices_x,[1]*len(indices_x), [1]*len(indices_x), 1)
        self.assertEquals(0, error)
        self.assertEquals([True,True,True,True,True,True,True,True,True,True,True] , is_active)
        
        
        indices_x = list(range(5,12,1))
        temperatures,error = instance.get_grid_electron_temperature(indices_x,[1]*len(indices_x), [1]*len(indices_x), 1)
        self.assertEquals(0, error)
        self.assertEquals([6000.0] * len(indices_x), temperatures)
        indices_x = list(range(1,5,1))
        temperatures,error = instance.get_grid_electron_temperature(indices_x,[1]*len(indices_x), [1]*len(indices_x), 1)
        self.assertEquals(0, error)
        self.assertEquals([6000.0] * len(indices_x), temperatures)
        
        ni, nj, nk, error = instance.get_max_indices(1)
        self.assertEquals(0, error)
        
        self.assertEquals(ni, 13)
        self.assertEquals(nj, 13)
        self.assertEquals(nj, 13)
        instance.stop()
        
        
    def xtest3(self):
        instance=self.new_instance(MocassinInterface, debugger = "ddd")
        #instance.redirect_outputs_to("moc3-out.txt", "moc3-err.txt")
        instance.initialize_code()
        instance.set_symmetricXYZ(True)
        instance.setup_mesh(13,13,13,0.95E+19,0.95E+19,0.95E+19)
        instance.setup_abundancies()
        #instance.set_abundancies_filename('abunHII20.in')
        instance.set_input_directory(instance.get_default_input_directory())
        instance.set_constant_hydrogen_density(100.0)
        instance.set_initial_nebular_temperature(6000.0)
        instance.set_maximum_number_of_monte_carlo_iterations(20)
        instance.set_minimum_convergence_level(100)
        instance.set_total_number_of_photons(10000000)
        instance.set_total_number_of_points_in_frequency_mesh(600)
        instance.set_high_limit_of_the_frequency_mesh(15.)
        instance.set_low_limit_of_the_frequency_mesh(1.001e-5)
        instance.set_convergence_limit(0.09)
        instance.set_number_of_ionisation_stages(6)
        instance.setup_auto_convergence(0.8, 2.0, 1000000000)
        #instance.set_emit_rate_of_photons(1.006e13)
        instance.commit_parameters()
        instance.define_stars(0.0, 0.0, 0.0, 20000.0, 6003.6396)
        instance.commit_particles()
        instance.commit_grid()
        
        x, error = instance.get_number_of_elements_used()
        self.assertEquals(0, error)
        self.assertEquals(x, 7)
        
        instance.iterate()
        
        
        indices_x = list(range(1,12,1))
        temperatures,error = instance.get_grid_electron_temperature(indices_x,[1]*len(indices_x), [1]*len(indices_x), 1)
        self.assertEquals(0, error)
        print temperatures
        instance.stop()



    def test4(self):
        instance=self.new_instance(MocassinInterface)
        instance.initialize_code()
        instance.setup_mesh(3,3,3,0.95E+19,0.95E+19,0.95E+19)
        instance.setup_abundancies()
        instance.set_input_directory(instance.get_default_input_directory())
        instance.set_initial_nebular_temperature(6000.0)
        instance.set_maximum_number_of_monte_carlo_iterations(20)
        instance.set_minimum_convergence_level(100)
        instance.set_total_number_of_photons(10000000)
        instance.set_total_number_of_points_in_frequency_mesh(600)
        instance.set_high_limit_of_the_frequency_mesh(15)
        instance.set_low_limit_of_the_frequency_mesh(1.001e-5)
    
        instance.set_convergence_limit(0.09)
        instance.set_number_of_ionisation_stages(6)
        instance.setup_auto_convergence(0.2, 2.0, 1000000000)
        instance.commit_parameters()
        instance.set_grid_hydrogen_density(1,1,1, 100)
    
        value,error = instance.get_grid_hydrogen_density(1,1,1)
    
        self.assertEquals(error, 0)
        self.assertEquals(value, 100)
    
        is_active,error = instance.get_grid_active(1,1,1)
        self.assertEquals(error, 0)
        self.assertFalse(is_active)
        
        is_active,error = instance.get_grid_active(1,2,1)
    
        self.assertEquals(error, 0)
        self.assertFalse(is_active)
        
        value,error = instance.get_grid_hydrogen_density(1,2,1,1)
    
        self.assertEquals(error, 0)
        self.assertEquals(value, 0)
    
        instance.commit_particles()
        instance.commit_grid()
        
        is_active,error = instance.get_grid_active(1,1,1)
    
        self.assertEquals(error, 0)
        self.assertTrue(is_active)
    
        value,error = instance.get_grid_hydrogen_density(1,1,1)
        self.assertEquals(error, 0)
        self.assertEquals(value, 100)
        
        value,error = instance.get_grid_ion_density(1,1,1,1,1)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEquals(value, 1e-5, 6)
    
        is_active,error = instance.get_grid_active(1,2,1)
    
        self.assertEquals(error, 0)
        self.assertFalse(is_active)
        
        value,error = instance.get_grid_hydrogen_density(1,2,1)
        self.assertEquals(error, 0)
        self.assertEquals(value, 0)
    
    
        instance.stop()

class TestMocassin(TestWithMPI):
    
    def test1(self):
        instance=self.new_instance(Mocassin)
        instance.initialize_code()
        
        self.assertEquals(0.0 | units.cm**-3, instance.parameters.constant_hydrogen_density)
        print instance.parameters.abundancies_filename
        instance.parameters.constant_hydrogen_density = 100.0 | units.cm**-3
        self.assertEquals(100.0 | units.cm**-3, instance.parameters.constant_hydrogen_density)
        
        self.assertEquals(10000 | units.K, instance.parameters.initial_nebular_temperature)
        
        self.assertEquals("" | units.string, instance.parameters.abundancies_filename)
        instance.stop()
        
        
    def test2(self):
        instance=self.new_instance(Mocassin) #, redirection = "none")
        instance.initialize_code()
        instance.set_input_directory(instance.get_default_input_directory())
        
        
        instance.set_symmetricXYZ(True)
        
        instance.parameters.nx = 11
        instance.parameters.ny = 12
        instance.parameters.nz = 13
        
        instance.parameters.length_x = 1 | units.km
        instance.parameters.length_y = 1 | units.km
        instance.parameters.length_z = 1 | units.km
        
        instance.commit_parameters()
        self.assertEquals(instance.grid.shape[0], 11)
        self.assertEquals(instance.grid.shape[1], 12)
        self.assertEquals(instance.grid.shape[2], 13)
        instance.stop()
        
    def test3(self):
        instance=self.new_instance(Mocassin) #, debugger = "xterm")
        instance.initialize_code()
        instance.set_random_seed(1)
        instance.set_input_directory(instance.get_default_input_directory())
        instance.set_output_directory(instance.get_default_output_directory())
        
        instance.set_initial_nebular_temperature(200.0 | units.K)
        
        instance.parameters.nx = 7
        instance.parameters.ny = 7
        instance.parameters.nz = 7
        
        print (0.95E+19 | units.cm).value_in(units.parsec)
        instance.parameters.length_x = 0.95E+19 | units.cm
        instance.parameters.length_y = 0.95E+19 | units.cm
        instance.parameters.length_z = 0.95E+19 | units.cm
        
        instance.set_high_limit_of_the_frequency_mesh(15 | units.ryd)
        instance.set_low_limit_of_the_frequency_mesh(1.001e-5| units.ryd)
        instance.set_maximum_number_of_monte_carlo_iterations(1)
        instance.set_total_number_of_photons(100)
        #instance.set_constant_hydrogen_density(100 | units.cm**-3)
        instance.commit_parameters()
        instance.grid.hydrogen_density = 100 | units.cm**-3
        instance.commit_grid()
        
        p = datamodel.Particle()
        p.x = 0 | units.cm
        p.y = 0 | units.cm
        p.z = 0 | units.cm
        p.temperature = 20000 | units.K
        p.luminosity = 1.0  | units.LSun
        
        instance.particles.add_particle(p)
        
        instance.commit_particles()
        
        #print instance.grid.electron_density[3]
        #print instance.grid.electron_temperature[3]
        print instance.ion_density_grid.density[3][1][2][0][0]
        self.assertAlmostRelativeEquals(1e-5 | units.none, instance.ion_density_grid.density[3][1][2][0][0], 7)
        self.assertAlmostRelativeEquals(1e-5 | units.none , instance.ion_density_grid.density[3][1][3][0][0], 7)
        
        instance.step()
        
        print instance.ion_density_grid.density[3][1][2][0][0]
        
        self.assertAlmostRelativeEquals(0.0,  instance.get_percentage_converged())
        self.assertAlmostRelativeEquals(0.00297847623006 | units.cm**-3 , instance.grid.electron_density[3][1][2], 5)
        self.assertAlmostRelativeEquals(0.00351035199128 | units.cm**-3 , instance.grid.electron_density[3][1][3], 5)
        self.assertAlmostRelativeEquals(0.99998998642 | units.none, instance.ion_density_grid.density[3][1][2][0][0], 7)
        self.assertAlmostRelativeEquals(0.99998998642 | units.none , instance.ion_density_grid.density[3][1][3][0][0], 7)
        instance.stop()
