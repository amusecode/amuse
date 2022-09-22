from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from amuse import datamodel

from amuse.community.mpiamrvac.interface import MpiAmrVacInterface
from amuse.community.mpiamrvac.interface import MpiAmrVac
from amuse import io
from amuse.units import generic_unit_system

import os
import numpy
class TestMpiAmrVacInterface(TestWithMPI):
    
    def test1(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVacInterface)
        instance.initialize_code()
        instance.stop()
    
    def test2(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVacInterface)
        filename, error =  instance.get_parameters_filename()
        self.assertEqual(error, 0)
        self.assertEqual(filename, "amrvac.par")
        error =  instance.set_parameters_filename("dontexists.file")
        self.assertEqual(error, -1)
        
        filename, error =  instance.get_parameters_filename()
        self.assertEqual(error, 0)
        self.assertEqual(filename, "amrvac.par")
        
        name_of_parametersfile = 'amrvac.tst.par'
        with open(name_of_parametersfile, 'w') as f:
            f.write('test param')
        error =  instance.set_parameters_filename(name_of_parametersfile)
        self.assertEqual(error, 0)
        
        filename, error =  instance.get_parameters_filename()
        self.assertEqual(error, 0)
        self.assertEqual(filename, name_of_parametersfile)
        
        os.remove(name_of_parametersfile)
        
        instance.initialize_code()
        instance.stop()
        

    def test3(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVacInterface)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEqual(error, 0)
        instance.stop()
        
    def test4(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVacInterface)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEqual(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.commit_parameters()
        self.assertEqual(error, 0)
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEqual(error, 0)
        self.assertEqual(number_of_grids, 8)
        
        lx, ly, lz, error = instance.get_mesh_size(1)
        self.assertEqual(lx, 10)
        self.assertEqual(ly, 10)
        self.assertEqual(lz, 10)
        x,y,z, error = instance.get_position_of_index(0,0,0,1)
        self.assertEqual(error, 0)
        self.assertEqual(x % 0.5, 0)
        self.assertEqual(y % 0.5, 0)
        self.assertEqual(z % 0.5, 0)
        instance.stop()
        #self.assertTrue(False) 
    
    
    def test5(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVacInterface)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEqual(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.commit_parameters()
        self.assertEqual(error, 0)
        
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEqual(error, 0)
        self.assertEqual(number_of_grids, 8)
        
        rho, error = instance.get_grid_density(1,1,1, 1)
        self.assertEqual(error, 0)
        self.assertEqual(rho, 0)
        error = instance.set_grid_density(1,1,1, 0.1, 1)
        self.assertEqual(error, 0)
        rho, error = instance.get_grid_density(1,1,1, 1)
        self.assertEqual(error, 0)
        self.assertEqual(rho, 0.1)
        rho, error = instance.get_grid_density(1,1,1, 2)
        self.assertEqual(error, 0)
        self.assertEqual(rho, 0.0)
        
        instance.stop()
        
    
    
    def test6(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVacInterface)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEqual(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.commit_parameters()
        self.assertEqual(error, 0)
        
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEqual(error, 0)
        self.assertEqual(number_of_grids, 8)
        
        rho, error = instance.get_grid_energy_density(1,1,1, 1)
        self.assertEqual(error, 0)
        self.assertEqual(rho, 0)
        error = instance.set_grid_energy_density(1,1,1, 0.1, 1)
        self.assertEqual(error, 0)
        rho, error = instance.get_grid_energy_density(1,1,1, 1)
        self.assertEqual(error, 0)
        self.assertEqual(rho, 0.1)
        rho, error = instance.get_grid_density(1,1,1, 1)
        self.assertEqual(error, 0)
        self.assertEqual(rho, 0.0)
        rho, error = instance.get_grid_energy_density(1,1,1, 2)
        self.assertEqual(error, 0)
        self.assertEqual(rho, 0.0)
        instance.stop()

    def test7(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVacInterface)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEqual(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.commit_parameters()
        self.assertEqual(error, 0)
        
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEqual(error, 0)
        self.assertEqual(number_of_grids, 8)
        
        rhovx, rhovy, rhovz, error = instance.get_grid_momentum_density(1,1,1, 1)
        self.assertEqual(error, 0)
        self.assertEqual(rhovx, 0.0)
        self.assertEqual(rhovy, 0.0)
        self.assertEqual(rhovz, 0.0)
        error = instance.set_grid_momentum_density(1,1,1, 0.1, 0.2, 0.3,  1)
        self.assertEqual(error, 0)
        rhovx, rhovy, rhovz, error = instance.get_grid_momentum_density(1,1,1, 1)
        self.assertEqual(error, 0)
        self.assertEqual(rhovx, 0.1)
        self.assertEqual(rhovy, 0.2)
        self.assertEqual(rhovz, 0.3)
        rhovx, rhovy, rhovz, error = instance.get_grid_momentum_density(1,1,1, 2)
        self.assertEqual(error, 0)
        self.assertEqual(rhovx, 0.0)
        self.assertEqual(rhovy, 0.0)
        self.assertEqual(rhovz, 0.0)
        rho, error = instance.get_grid_energy_density(1,1,1, 1)
        self.assertEqual(error, 0)
        self.assertEqual(rho, 0.0)
                        
        instance.stop()
        
    def test8(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVacInterface)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEqual(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.set_boundary("periodic", "periodic", "periodic", "periodic", "periodic", "periodic")
        
        self.assertEqual(error, 0)
        error = instance.commit_parameters()
        self.assertEqual(error, 0)
        
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEqual(error, 0)
        self.assertEqual(number_of_grids, 8)
        error = instance.initialize_grid()
        self.assertEqual(error, 0)
                                
        instance.stop()
        
    
        
    def test9(self): 

        instance = self.new_instance_of_an_optional_code(MpiAmrVacInterface)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEqual(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        instance.commit_parameters()
        n1, n2, n3, error = instance.get_acceleration_grid_size()
        self.assertEqual(error, 0)
        self.assertEqual(n1, 50)
        self.assertEqual(n2, 50)
        self.assertEqual(n3, 50)
        
        a1, a2, a3, error = instance.get_acceleration_grid_acceleration(1,1,1)
        self.assertEqual(error, 0)
        self.assertEqual(a1, 0.0)
        self.assertEqual(a2, 0)
        self.assertEqual(a3, 0)
        
        error = instance.set_acceleration_grid_acceleration(1,1,1, 10.0, 20.0, 30.0)
        self.assertEqual(error, 0)
        
        a1, a2, a3, error = instance.get_acceleration_grid_acceleration(1,1,1)
        self.assertEqual(error, 0)
        self.assertEqual(a1, 10.0)
        self.assertEqual(a2, 20.0)
        self.assertEqual(a3, 30.0)
        
        x,y,z, error = instance.get_acceleration_grid_position_of_index(1,1,1)
        self.assertEqual(error, 0)
        
        self.assertEqual(x, -1)
        self.assertEqual(y, -1)
        self.assertEqual(z, -1)
        instance.stop()
        
    


class TestMpiAmrVac(TestWithMPI):
    
    def test1(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVac)
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.setup_mesh(20,20,20, 20.0 | generic_unit_system.length, 20.0 | generic_unit_system.length, 20.0 | generic_unit_system.length)
        error = instance.set_boundary("periodic", "periodic", "periodic", "periodic", "periodic", "periodic")
        
        rhovx, rhovy, rhovz = instance.get_grid_momentum_density(1,1,1, 1)

        self.assertEqual(rhovx, 0.0 | generic_unit_system.momentum_density)
        self.assertEqual(rhovy, 0.0 | generic_unit_system.momentum_density)
        self.assertEqual(rhovz, 0.0 | generic_unit_system.momentum_density)

        instance.stop()
        
    def test2(self):
    
        instance = self.new_instance_of_an_optional_code(MpiAmrVac)
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.parameters.mesh_length = (20.0, 20.0, 20.0) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 20, 20) 
        instance.set_boundary("periodic", "periodic", "periodic", "periodic", "periodic", "periodic")
        
        grids = list(instance.itergrids())
        
        self.assertEqual(len(grids), 8)
        
        for grid in grids:
            position = grids[0].position
            for i in range(3):
                max_x = position[...,i].amax()
                min_x = position[...,i].amin()
                
                self.assertTrue(min_x >= 0.5 | generic_unit_system.length)
                self.assertTrue(min_x <= 11.5 | generic_unit_system.length)
                self.assertTrue(max_x >= 9.5 | generic_unit_system.length)
                self.assertTrue(max_x <= 19.5 | generic_unit_system.length)
                            
        
        
        
        self.assertEqual(grids[0][0][0][0].rho,  0.0 | generic_unit_system.density)
        
        grids[0].rho = 0.2 | generic_unit_system.density
        
        rho1 = grids[1].rho
        rho0 = grids[0].rho
        for i in range(10):
            for j in range(10):
                for k in range(10):
                    self.assertEqual(rho1[0][0][0],  0.0 | generic_unit_system.density)
                    self.assertEqual(rho0[0][0][0],  0.2 | generic_unit_system.density)
                    
        instance.stop()
    
    
    def test3(self):
    
        for number_of_workers in range(2,6):
            instance = self.new_instance_of_an_optional_code(MpiAmrVac, number_of_workers = number_of_workers)
            instance.set_parameters_filename(instance.default_parameters_filename)
            instance.parameters.mesh_length = (20.0, 20.0, 20.0) | generic_unit_system.length
            instance.parameters.mesh_size = (20, 20, 20) 
            instance.parameters.x_boundary_conditions = ("periodic","periodic")
            instance.parameters.y_boundary_conditions = ("periodic","periodic")
            instance.parameters.z_boundary_conditions = ("periodic","periodic")

            grids = list(instance.itergrids())
            
            self.assertEqual(len(grids), 8)
            
            for index, grid in enumerate(grids):
                position = grid.position
                #print instance.get_level_of_grid(index + 1)
                level = instance.get_level_of_grid(index + 1)
                self.assertEqual(level, 1)
                for i in range(3):
                    max_x = position[...,i].amax()
                    min_x = position[...,i].amin()
                    
                    self.assertTrue(min_x >= 0.5 | generic_unit_system.length)
                    self.assertTrue(min_x <= 10.5 | generic_unit_system.length)
                    self.assertTrue(max_x >= 9.5 | generic_unit_system.length)
                    self.assertTrue(max_x <= 19.5 | generic_unit_system.length)
                    self.assertEqual(max_x - min_x , 9.0 | generic_unit_system.length)
                    
            self.assertEqual(grids[0][0][0][0].rho,  0.0 | generic_unit_system.density)
            
            grids[0].rho = 0.2 | generic_unit_system.density
            
            rho1 = grids[4].rho
            rho0 = grids[0].rho
            for i in range(10):
                for j in range(10):
                    for k in range(10):
                        self.assertEqual(rho1[0][0][0],  0.0 | generic_unit_system.density)
                        self.assertEqual(rho0[0][0][0],  0.2 | generic_unit_system.density)
                        
            instance.stop()
        
    def test4(self):
    
        instance = self.new_instance_of_an_optional_code(MpiAmrVac, number_of_workers = 1)
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.parameters.mesh_length = (10.0,10.0, 10.0) | generic_unit_system.length
        instance.parameters.mesh_size = (10, 10, 10)
        instance.parameters.maximum_number_of_grid_levels = 5
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
     
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        
        levels = [instance.get_level_of_grid(i+1) for i in range(len(grids))]
        self.assertEqual(levels, [1] * 1)
        
        gamma = 5.0 / 3.0
        energy =  generic_unit_system.mass / (generic_unit_system.time**2 * generic_unit_system.length)
        grids_in_memory = []
        for grid in grids:
            grid.rho = 1.0  | generic_unit_system.density
            grid.energy = (0.1795 | energy)/ (gamma - 1)
    
        
        has_advanced = instance.refine_grid()
        grids = list(instance.itergrids())
        self.assertFalse(has_advanced)
        instance.stop()

    def test5(self):
        def fill_grids(grids):
            for grid in grids:
                firsthalf = grid.x > 5.0 | generic_unit_system.length
                secondhalf = numpy.logical_not(firsthalf)
                if(numpy.any(firsthalf)):
                    grid[firsthalf].rho = 1.0  | generic_unit_system.density
        
        instance = self.new_instance_of_an_optional_code(MpiAmrVac, number_of_workers = 1)
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.parameters.mesh_length = (10.0,10.0, 10.0) | generic_unit_system.length
        instance.parameters.mesh_size = (10, 10, 10)
        instance.parameters.maximum_number_of_grid_levels = 4
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
     

        n = 1
        for i in range(3):
            expected_length = n ** 3
            grids = list(instance.itergrids())
            levels = [instance.get_level_of_grid(grid+1) for grid in range(len(grids))]
            self.assertEqual(len(grids), expected_length)
            
            self.assertEqual(levels, [i+1] * expected_length)
            
            n *= 2
            
            fill_grids(grids)
            has_advanced = instance.refine_grid()
            self.assertTrue(has_advanced)
        
        has_advanced = instance.refine_grid()
        self.assertFalse(has_advanced)
        
        
        instance.stop()
        
    
    def test6(self):
        instance = self.new_instance_of_an_optional_code(MpiAmrVac, number_of_workers = 1)
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.parameters.mesh_length = (10.0,10.0, 10.0) | generic_unit_system.length
        instance.parameters.mesh_size = (10, 10, 10)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        
        self.assertEqual(instance.parameters.maximum_number_of_grid_levels, 3)
        instance.parameters.maximum_number_of_grid_levels = 2
        self.assertEqual(instance.parameters.maximum_number_of_grid_levels, 2)
        
        self.assertEqual(instance.parameters.entropy_type , 'nul')
        instance.parameters.entropy_type = 'powell'
        self.assertEqual(instance.parameters.entropy_type , 'powell')
        
        self.assertEqual(instance.parameters.time_integration_procedure , 'twostep')
        instance.parameters.time_integration_procedure = 'onestep'
        self.assertEqual(instance.parameters.time_integration_procedure , 'onestep')
        
        self.assertEqual(instance.parameters.spatial_discretization_method , 'tvdmu')
        instance.parameters.spatial_discretization_method = 'tvdlf'
        self.assertEqual(instance.parameters.spatial_discretization_method , 'tvdlf')
        
        
        
        self.assertEqual(instance.parameters.predictor_step_discretization_method , 'tvdmu')
        instance.parameters.predictor_step_discretization_method = 'hancock'
        self.assertEqual(instance.parameters.predictor_step_discretization_method , 'hancock')
        
        instance.commit_parameters()
        
        self.assertEqual(instance.parameters.entropy_type , 'powell')
        self.assertEqual(instance.parameters.maximum_number_of_grid_levels, 2)
        self.assertEqual(instance.parameters.time_integration_procedure , 'onestep')
        self.assertEqual(instance.parameters.spatial_discretization_method , 'tvdlf')
        self.assertEqual(instance.parameters.predictor_step_discretization_method , 'hancock')
        
    def test7(self):
        instance=self.new_instance_of_an_optional_code(MpiAmrVac)
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.parameters.mesh_length = (10.0,10.0, 10.0) | generic_unit_system.length
        instance.parameters.mesh_size = (10,10,10)
        instance.parameters.maximum_number_of_grid_levels = 1
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        
    
        grid = datamodel.Grid(10,10,10)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy =  1.0 | generic_unit_system.energy_density
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        channel = grid.new_channel_to(igrid)
        channel.copy()
        
                   
        instance.evolve_model(0.1 | generic_unit_system.time)
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.evolve_model(0.3 | generic_unit_system.time)
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.stop()
        
    def test8(self):
        instance=self.new_instance_of_an_optional_code(MpiAmrVac)
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.parameters.mesh_size = (10,10,10)
        instance.parameters.maximum_number_of_grid_levels = 1
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 1.0 | generic_unit_system.length
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        
        name = instance.get_typeghostfill()
        self.assertEqual(name, 'linear')
        instance.commit_parameters()
        
        name = instance.get_typeghostfill()
        self.assertEqual(name, 'linear')
        grid = datamodel.Grid(10,10,10)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy =  1.0 | generic_unit_system.energy_density
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        channel = grid.new_channel_to(igrid)
        channel.copy()
        
        name = instance.get_typeghostfill()
        self.assertEqual(name, 'linear')
    
        self.assertEqual((50,50,50), instance.acceleration_grid.shape)
        
        acc_grid = datamodel.Grid(50,50,50)
        acceleration = 1 | generic_unit_system.acceleration
        acc_grid.ax = acceleration
        acc_grid.ay = acceleration
        acc_grid.az = acceleration
        #self.assertEquals(acc_grid.acceleration[0][0][0], ( 1,1,1) | generic_unit_system.acceleration)
        channel = acc_grid.new_channel_to(instance.acceleration_grid)
        channel.copy()
        
        name = instance.get_typeghostfill()
        self.assertEqual(name, 'linear')
        
        result = instance.initialize_grid()
                   
        name = instance.get_typeghostfill()
        self.assertEqual(name, 'linear')
    
        instance.evolve_model(0.1 | generic_unit_system.time)
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho);
        #self.assertAlmostRelativeEquals(igrid.rhovx, 0.1 * 1.0 * 0.1 | generic_unit_system.momentum_density);
        #self.assertAlmostRelativeEquals(igrid.rhovy, 0.1 * 1.0 * 0.1 | generic_unit_system.momentum_density);
        #self.assertAlmostRelativeEquals(igrid.rhovz, 0.1 * 1.0 * 0.1 | generic_unit_system.momentum_density);

        instance.evolve_model(0.3 | generic_unit_system.time)
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho);
        #self.assertAlmostRelativeEquals(igrid.rhovx, grid.rho *  instance.model_time * acceleration,2);
        #self.assertAlmostRelativeEquals(igrid.rhovy, grid.rho *  instance.model_time * acceleration,2);
        #self.assertAlmostRelativeEquals(igrid.rhovz, grid.rho *  instance.model_time * acceleration,2);
        instance.stop()
        
    
    def test9(self):
        instance=self.new_instance_of_an_optional_code(MpiAmrVac, mode="2d")
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.parameters.mesh_length = (10.0,10.0, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (10,10,1)
        instance.parameters.maximum_number_of_grid_levels = 1
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        
    
        grid = datamodel.Grid(10,10,1)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy =  1.0 | generic_unit_system.energy_density
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        channel = grid.new_channel_to(igrid)
        channel.copy()
        
                   
        instance.evolve_model(0.1 | generic_unit_system.time)
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.evolve_model(0.3 | generic_unit_system.time)
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.stop()
        
    
    def test10(self): 
        instance=self.new_instance_of_an_optional_code(MpiAmrVac)
        instance.parameters.mesh_size = (10,10,10)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.maximum_number_of_grid_levels = 1
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("periodic", "periodic")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        gamma = 5.0 / 3.0
        
        grid = datamodel.new_regular_grid((10,10,10), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.1 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid[0:5].rho = 0.015 | density
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        channel = grid.new_channel_to(igrid)
        channel.copy()
        
        instance.stopping_conditions.number_of_steps_detection.enable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        rho = next(instance.itergrids()).rho[...,0,0]
        self.assertAlmostRelativeEquals(rho[7], 0.01 | density)
        self.assertTrue(rho[0] < 0.015 | density)
        self.assertTrue( instance.model_time < 1.0 | generic_unit_system.time)
        
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(0.1 | generic_unit_system.time)
        rho = next(instance.itergrids()).rho[...,0,0]
        self.assertAlmostRelativeEquals( instance.model_time, 0.1 | generic_unit_system.time)
        self.assertAlmostRelativeEquals(rho[7], 0.012812 | density, 3)
        self.assertTrue(rho[0] < 0.015 | density)
        
        instance.stop()

    def test11(self): 
        instance=self.new_instance_of_an_optional_code(MpiAmrVac)
        instance.parameters.mesh_size = (10,10,10)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.maximum_number_of_grid_levels = 1
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("periodic", "periodic")
        instance.parameters.stopping_conditions_timeout = 0.1 | units.s
        
        gamma = 5.0 / 3.0
        
        grid = datamodel.new_regular_grid((10,10,10), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.1 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid[0:5].rho = 0.015 | density
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        channel = grid.new_channel_to(igrid)
        channel.copy()
        
        instance.stopping_conditions.timeout_detection.enable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        
        self.assertTrue(instance.stopping_conditions.timeout_detection.is_set())
        rho = next(instance.itergrids()).rho[...,0,0]
        self.assertAlmostRelativeEquals(rho.mean(), 0.0125 | density)
        self.assertTrue( instance.model_time < 1.0 | generic_unit_system.time)
        
        
        instance.stopping_conditions.timeout_detection.disable()
        tnext =  instance.model_time.round(2) + (0.2 | generic_unit_system.time)
        instance.evolve_model(tnext)
        rho = next(instance.itergrids()).rho[...,0,0]
        self.assertAlmostRelativeEquals( instance.model_time.round(2) ,tnext)
        self.assertAlmostRelativeEquals(rho.mean(), 0.0125 | density)
        self.assertTrue(rho[0] < 0.015 | density)
        
        instance.stop()
        
    
    def test12(self):
        instance=self.new_instance_of_an_optional_code(MpiAmrVac, mode="1d")
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.parameters.mesh_length = (10.0, 1, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (10, 1, 1)
        instance.parameters.maximum_number_of_grid_levels = 1
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        
    
        grid = datamodel.Grid(10,1,1)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy =  1.0 | generic_unit_system.energy_density
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        channel = grid.new_channel_to(igrid)
        channel.copy()
        
                   
        instance.evolve_model(0.1 | generic_unit_system.time)
        
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.evolve_model(0.3 | generic_unit_system.time)
        grids = list(instance.itergrids())
        self.assertEqual(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.stop()

    
    def test13(self):
        instance=self.new_instance_of_an_optional_code(MpiAmrVac, mode="2d-acc")
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.parameters.mesh_length = (10.0, 10.0, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 20, 1)
        instance.parameters.maximum_number_of_grid_levels = 2
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        
        must_refine = True      
        middle = 5.0 | generic_unit_system.length    
        while must_refine:
            must_refine = instance.refine_grid()
            
            for x in instance.itergrids():
                inmem = x.copy()
                inmem[inmem.x <  middle].rho= 0.3 | generic_unit_system.density
                inmem[inmem.x >= middle].rho = 0.1 | generic_unit_system.density
                inmem.rhovx = 0.0 | generic_unit_system.momentum_density
                inmem.rhovy = 0.0 |  generic_unit_system.momentum_density
                inmem.rhovz = 0.0 |  generic_unit_system.momentum_density
                inmem.ax = 0.2 | generic_unit_system.acceleration
                inmem.ay = 0.0 |  generic_unit_system.acceleration
                inmem.az = 0.0 |  generic_unit_system.acceleration
                
                inmem.energy =  1.0 | generic_unit_system.energy_density
                from_model_to_code = inmem.new_channel_to(x)
                from_model_to_code.copy()
                self.assertAlmostRelativeEquals(x.ax, 0.2 | generic_unit_system.acceleration)
                self.assertAlmostRelativeEquals(x.ay, 0.0 | generic_unit_system.acceleration)
                self.assertAlmostRelativeEquals(x.az, 0.0 | generic_unit_system.acceleration)
        
        self.assertEqual(len(list(instance.itergrids())), 4)
        for igrid in instance.itergrids():
            self.assertAlmostRelativeEquals(igrid.ax, 0.2 | generic_unit_system.acceleration)
            self.assertAlmostRelativeEquals(igrid.ay, 0.0 | generic_unit_system.acceleration)
            self.assertAlmostRelativeEquals(igrid.az, 0.0 | generic_unit_system.acceleration)
                   
        instance.evolve_model(0.1 | generic_unit_system.time)
        
        self.assertEqual(len(list(instance.itergrids())), 16)
        for igrid in instance.itergrids():
            self.assertAlmostRelativeEquals(igrid.ax, 0.2 | generic_unit_system.acceleration)
            self.assertAlmostRelativeEquals(igrid.ay, 0.0 | generic_unit_system.acceleration)
            self.assertAlmostRelativeEquals(igrid.az, 0.0 | generic_unit_system.acceleration)
                   
      
        instance.stop()
        
    
    def test14(self):
        for ax, ay in ((0.2,0.0), (0.0, 0.2), (0.2,0.2)) | generic_unit_system.acceleration:
            instance=self.new_instance_of_an_optional_code(MpiAmrVac, mode="2d-acc")
            instance.set_parameters_filename(instance.default_parameters_filename)
            instance.parameters.mesh_length = (10.0, 10.0, 1) | generic_unit_system.length
            instance.parameters.mesh_size = (10, 10, 1)
            instance.parameters.maximum_number_of_grid_levels = 1
            instance.parameters.x_boundary_conditions = ("periodic","periodic")
            instance.parameters.y_boundary_conditions = ("periodic","periodic")
            rho = 0.1 | generic_unit_system.density
            middle = 5.0 | generic_unit_system.length    
            for x in instance.itergrids():
                inmem = x.copy()
                inmem.rho = 0.1 | generic_unit_system.density
                inmem.rhovx = 0.0 | generic_unit_system.momentum_density
                inmem.rhovy = 0.0 |  generic_unit_system.momentum_density
                inmem.ax = ax
                inmem.ay = ay
                
                inmem.energy =  1.0 | generic_unit_system.energy_density
                from_model_to_code = inmem.new_channel_to(x)
                from_model_to_code.copy()
    
            self.assertEqual(len(list(instance.itergrids())), 1)
            
            dt = 0.1 | generic_unit_system.time
            instance.evolve_model(dt)
            
            self.assertEqual(len(list(instance.itergrids())), 1)
            igrid = list(instance.itergrids())[0]
            
            self.assertAlmostRelativeEquals(igrid.rho , 0.1 | generic_unit_system.density)
            self.assertAlmostRelativeEquals(igrid.rhovx , ax * dt* rho)
            self.assertAlmostRelativeEquals(igrid.rhovy , ay * dt* rho)
            instance.stop()
            
    
    
    def test15(self):
        for ax, ay in ((0.2,0.0), (0.0, 0.2), (0.2,0.2)) | generic_unit_system.acceleration: #(0.0,0.0), 
            instance=self.new_instance_of_an_optional_code(MpiAmrVac, mode="2d-acc")
            instance.set_parameters_filename(instance.default_parameters_filename)
            instance.parameters.mesh_length = (10.0, 10.0, 1) | generic_unit_system.length
            instance.parameters.mesh_size = (20, 20, 1)
            instance.parameters.maximum_number_of_grid_levels = 1
            instance.parameters.time_accurate = False
            instance.parameters.x_boundary_conditions = ("periodic","periodic")
            instance.parameters.y_boundary_conditions = ("periodic","periodic")
            rho = 0.1 | generic_unit_system.density
            middle = 5.0 | generic_unit_system.length    
            for x in instance.itergrids():
                inmem = x.copy()
                inmem.rho = 0.1 | generic_unit_system.density
                inmem.rhovx = 0.0 | generic_unit_system.momentum_density
                inmem.rhovy = 0.0 |  generic_unit_system.momentum_density
                inmem.ax = ax
                inmem.ay = ay
                
                inmem.energy =  1.0 | generic_unit_system.energy_density
                from_model_to_code = inmem.new_channel_to(x)
                from_model_to_code.copy()
    
            self.assertEqual(len(list(instance.itergrids())), 4)
            
            dt = 0.1 | generic_unit_system.time
            instance.evolve_model(dt)
            dt = instance.model_time
            #~ for igrid in instance.itergrids():
                #~ print(ax * dt* rho, ax, instance.model_time)
            self.assertEqual(len(list(instance.itergrids())), 4)
              
            for igrid in instance.itergrids():
                #print igrid.rho[2]
                self.assertAlmostRelativeEquals(igrid.rho  , rho)
                self.assertAlmostRelativeEquals(igrid.rhovx, ax * dt* rho)
                self.assertAlmostRelativeEquals(igrid.rhovy, ay * dt* rho)
            instance.stop()
            
    
    def test16(self):
        
        instance=self.new_instance_of_an_optional_code(MpiAmrVac, mode="1d")
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (20.0, 1, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 1, 1)
        instance.parameters.maximum_number_of_grid_levels = 1
        
        for x in instance.itergrids():
            inmem = x.copy()
            inmem.rho = inmem.x/(1| generic_unit_system.length) | generic_unit_system.density
            inmem.rhovx = 0.0 | generic_unit_system.momentum_density
            inmem.energy =  1.0 | generic_unit_system.energy_density
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(0.5| generic_unit_system.length,0.0| generic_unit_system.length,0.0| generic_unit_system.length)
        
        self.assertEqual(rho , 0.5 | generic_unit_system.density)
        
        for value in numpy.arange(0.5, 19.6, 0.1):
            
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value | generic_unit_system.length,
                0.0 | generic_unit_system.length,
                0.0 | generic_unit_system.length
            )
        
            self.assertAlmostRelativeEquals(rho , value | generic_unit_system.density)
        
        for value in numpy.arange(0.0, 0.6, 0.1):
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value | generic_unit_system.length,
                0.0 | generic_unit_system.length,
                0.0 | generic_unit_system.length
            )
            self.assertAlmostRelativeEquals(rho , ((0.5 + value) * 0.5 + (0.5-value) * 19.5) | generic_unit_system.density)
        
        
        for value in numpy.arange(0.0, 0.5, 0.1):
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value + 19.5| generic_unit_system.length,
                0.0 | generic_unit_system.length,
                0.0 | generic_unit_system.length
            )
            self.assertAlmostRelativeEquals(rho , (19.5 - (value * 19))  | generic_unit_system.density, 9)
        
        # out of range
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
            20.0| generic_unit_system.length,
            0.0 | generic_unit_system.length,
            0.0 | generic_unit_system.length
        )
        self.assertAlmostRelativeEquals(rho , 0.0 | generic_unit_system.density, 9)
        
    def test17(self):
        
        instance=self.new_instance_of_an_optional_code(MpiAmrVac, mode="2d", number_of_workers=2)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (20.0, 20.0, 1) | generic_unit_system.length
        instance.parameters.mesh_length = (20.0, 20.0, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 20, 1)
        instance.parameters.maximum_number_of_grid_levels = 1
        
        for x in instance.itergrids():
            inmem = x.copy()
            inmem.rho = (inmem.x + ((inmem.y - (0.5| generic_unit_system.length))* 20.0))/(1| generic_unit_system.length) | generic_unit_system.density
            inmem.rhovx = 0.0 | generic_unit_system.momentum_density
            inmem.energy =  1.0 | generic_unit_system.energy_density
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(0.5| generic_unit_system.length,0.5| generic_unit_system.length,0.0| generic_unit_system.length)
        
        self.assertEqual(rho , 0.5 | generic_unit_system.density)
        
        for value in numpy.arange(0.5, 19.6, 0.1):
            
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value | generic_unit_system.length,
                0.5 | generic_unit_system.length,
                0.0 | generic_unit_system.length
            )
        
            self.assertAlmostRelativeEquals(rho , value | generic_unit_system.density)
        
        for x in numpy.arange(8.5, 11.5, 0.25):
            for y in numpy.arange(0.5, 19.6, 0.25):
                rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                    x | generic_unit_system.length,
                    y | generic_unit_system.length,
                    0.0 | generic_unit_system.length
                )
            
                self.assertAlmostRelativeEquals(rho , x + (20 * (y-0.5))  | generic_unit_system.density)
            
    
    def test18(self):
        
        instance=self.new_instance_of_an_optional_code(MpiAmrVac, number_of_workers=3)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (20.0, 20.0, 20.0) | generic_unit_system.length
        instance.parameters.mesh_length = (20.0, 20.0, 20.0) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 20, 20)
        instance.parameters.maximum_number_of_grid_levels = 1
        
        for x in instance.itergrids():
            inmem = x.copy()
            inmem.rho = (
                (
                    inmem.x + 
                    ((inmem.y - (0.5| generic_unit_system.length))* 20.0) +
                    ((inmem.z - (0.5| generic_unit_system.length))* 400.0)
                )
                /(1| generic_unit_system.length) | generic_unit_system.density
            )
            inmem.rhovx = 0.0 | generic_unit_system.momentum_density
            inmem.energy =  1.0 | generic_unit_system.energy_density
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(0.5| generic_unit_system.length,0.5| generic_unit_system.length,0.5| generic_unit_system.length)
        
        self.assertEqual(rho , 0.5 | generic_unit_system.density)
        
        for value in numpy.arange(0.5, 19.6, 0.1):
            
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value | generic_unit_system.length,
                0.5 | generic_unit_system.length,
                0.5 | generic_unit_system.length
            )
        
            self.assertAlmostRelativeEquals(rho , value | generic_unit_system.density)
        
        sample = sample = datamodel.new_regular_grid(
            (4, 4, 76),
            (2, 2, 19) | generic_unit_system.length
        )
        sample.x += 9.5 | generic_unit_system.length
        sample.y += 9.5 | generic_unit_system.length
        sample.z += 0.5 | generic_unit_system.length
        x = sample.x.flatten()
        y = sample.y.flatten()
        z = sample.z.flatten()
        
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
            x,
            y,
            z
        )
        half = 0.5 | generic_unit_system.length
        
        self.assertAlmostRelativeEquals(rho , (x + (20 * (y-half)) + (400 * (z-half)))/(1| generic_unit_system.length) | generic_unit_system.density )
            




