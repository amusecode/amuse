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
        instance = self.new_instance(MpiAmrVacInterface)
        instance.initialize_code()
        instance.stop()
    
    def test2(self):
        instance = self.new_instance(MpiAmrVacInterface)
        filename, error =  instance.get_parameters_filename()
        self.assertEquals(error, 0)
        self.assertEquals(filename, "amrvac.par")
        error =  instance.set_parameters_filename("dontexists.file")
        self.assertEquals(error, -1)
        
        filename, error =  instance.get_parameters_filename()
        self.assertEquals(error, 0)
        self.assertEquals(filename, "amrvac.par")
        
        name_of_parametersfile = 'amrvac.tst.par'
        with open(name_of_parametersfile, 'w') as f:
            f.write('test param')
        error =  instance.set_parameters_filename(name_of_parametersfile)
        self.assertEquals(error, 0)
        
        filename, error =  instance.get_parameters_filename()
        self.assertEquals(error, 0)
        self.assertEquals(filename, name_of_parametersfile)
        
        os.remove(name_of_parametersfile)
        
        instance.initialize_code()
        instance.stop()
        

    def test3(self):
        instance = self.new_instance(MpiAmrVacInterface)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.stop()
        
    def test4(self):
        instance = self.new_instance(MpiAmrVacInterface)
        print instance.default_parameters_filename
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.commit_parameters()
        self.assertEquals(error, 0)
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEquals(error, 0)
        self.assertEquals(number_of_grids, 8)
        
        lx, ly, lz, error = instance.get_mesh_size(1)
        self.assertEquals(lx, 10)
        self.assertEquals(ly, 10)
        self.assertEquals(lz, 10)
        x,y,z, error = instance.get_position_of_index(0,0,0,1)
        self.assertEquals(error, 0)
        self.assertEquals(x % 0.5, 0)
        self.assertEquals(y % 0.5, 0)
        self.assertEquals(z % 0.5, 0)
        instance.stop()
        #self.assertTrue(False) 
    
    
    def test5(self):
        instance = self.new_instance(MpiAmrVacInterface)
        print instance.default_parameters_filename
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.commit_parameters()
        self.assertEquals(error, 0)
        
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEquals(error, 0)
        self.assertEquals(number_of_grids, 8)
        
        rho, error = instance.get_grid_density(1,1,1, 1)
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0)
        error = instance.set_grid_density(1,1,1, 0.1, 1)
        self.assertEquals(error, 0)
        rho, error = instance.get_grid_density(1,1,1, 1)
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0.1)
        rho, error = instance.get_grid_density(1,1,1, 2)
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0.0)
        
        instance.stop()
        
    
    
    def test6(self):
        instance = self.new_instance(MpiAmrVacInterface)
        print instance.default_parameters_filename
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.commit_parameters()
        self.assertEquals(error, 0)
        
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEquals(error, 0)
        self.assertEquals(number_of_grids, 8)
        
        rho, error = instance.get_grid_energy_density(1,1,1, 1)
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0)
        error = instance.set_grid_energy_density(1,1,1, 0.1, 1)
        self.assertEquals(error, 0)
        rho, error = instance.get_grid_energy_density(1,1,1, 1)
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0.1)
        rho, error = instance.get_grid_density(1,1,1, 1)
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0.0)
        rho, error = instance.get_grid_energy_density(1,1,1, 2)
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0.0)
        instance.stop()

    def test7(self):
        instance = self.new_instance(MpiAmrVacInterface)
        print instance.default_parameters_filename
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.commit_parameters()
        self.assertEquals(error, 0)
        
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEquals(error, 0)
        self.assertEquals(number_of_grids, 8)
        
        rhovx, rhovy, rhovz, error = instance.get_grid_momentum_density(1,1,1, 1)
        self.assertEquals(error, 0)
        self.assertEquals(rhovx, 0.0)
        self.assertEquals(rhovy, 0.0)
        self.assertEquals(rhovz, 0.0)
        error = instance.set_grid_momentum_density(1,1,1, 0.1, 0.2, 0.3,  1)
        self.assertEquals(error, 0)
        rhovx, rhovy, rhovz, error = instance.get_grid_momentum_density(1,1,1, 1)
        self.assertEquals(error, 0)
        self.assertEquals(rhovx, 0.1)
        self.assertEquals(rhovy, 0.2)
        self.assertEquals(rhovz, 0.3)
        rhovx, rhovy, rhovz, error = instance.get_grid_momentum_density(1,1,1, 2)
        self.assertEquals(error, 0)
        self.assertEquals(rhovx, 0.0)
        self.assertEquals(rhovy, 0.0)
        self.assertEquals(rhovz, 0.0)
        rho, error = instance.get_grid_energy_density(1,1,1, 1)
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0.0)
                        
        instance.stop()
        
    def test8(self):
        instance = self.new_instance(MpiAmrVacInterface)
        print instance.default_parameters_filename
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.set_boundary("periodic", "periodic", "periodic", "periodic", "periodic", "periodic")
        
        self.assertEquals(error, 0)
        error = instance.commit_parameters()
        self.assertEquals(error, 0)
        
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEquals(error, 0)
        self.assertEquals(number_of_grids, 8)
        error = instance.initialize_grid()
        self.assertEquals(error, 0)
                                
        instance.stop()
        
    
        
    def test9(self): 

        instance = self.new_instance(MpiAmrVacInterface)
        print instance.default_parameters_filename
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        instance.commit_parameters()
        n1, n2, n3, error = instance.get_acceleration_grid_size()
        self.assertEquals(error, 0)
        self.assertEquals(n1, 50)
        self.assertEquals(n2, 50)
        self.assertEquals(n3, 50)
        
        a1, a2, a3, error = instance.get_acceleration_grid_acceleration(1,1,1)
        self.assertEquals(error, 0)
        self.assertEquals(a1, 0.0)
        self.assertEquals(a2, 0)
        self.assertEquals(a3, 0)
        
        error = instance.set_acceleration_grid_acceleration(1,1,1, 10.0, 20.0, 30.0)
        self.assertEquals(error, 0)
        
        a1, a2, a3, error = instance.get_acceleration_grid_acceleration(1,1,1)
        self.assertEquals(error, 0)
        self.assertEquals(a1, 10.0)
        self.assertEquals(a2, 20.0)
        self.assertEquals(a3, 30.0)
        
        x,y,z, error = instance.get_acceleration_grid_position_of_index(1,1,1)
        self.assertEquals(error, 0)
        
        self.assertEquals(x, -1)
        self.assertEquals(y, -1)
        self.assertEquals(z, -1)
        instance.stop()
        
    


class TestMpiAmrVac(TestWithMPI):
    
    def test1(self):
        instance = self.new_instance(MpiAmrVac)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.setup_mesh(20,20,20, 20.0 | generic_unit_system.length, 20.0 | generic_unit_system.length, 20.0 | generic_unit_system.length)
        error = instance.set_boundary("periodic", "periodic", "periodic", "periodic", "periodic", "periodic")
        error = instance.commit_parameters()
        
        rhovx, rhovy, rhovz = instance.get_grid_momentum_density(1,1,1, 1)

        self.assertEquals(rhovx, 0.0 | generic_unit_system.momentum_density)
        self.assertEquals(rhovy, 0.0 | generic_unit_system.momentum_density)
        self.assertEquals(rhovz, 0.0 | generic_unit_system.momentum_density)

        instance.stop()
        
    def test2(self):
    
        instance = self.new_instance(MpiAmrVac)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.parameters.mesh_length = (20.0, 20.0, 20.0) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 20, 20) 
        
        #instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
        error = instance.set_boundary("periodic", "periodic", "periodic", "periodic", "periodic", "periodic")
        error = instance.commit_parameters()

        grids = list(instance.itergrids())
        
        self.assertEquals(len(grids), 8)
        
        for grid in grids:
            position = grids[0].position
            print position.shape
            for i in range(3):
                max_x = position[...,...,...,i].amax()
                min_x = position[...,...,...,i].amin()
                print max_x, min_x
                
                self.assertTrue(min_x >= 0.5 | generic_unit_system.length)
                self.assertTrue(min_x <= 11.5 | generic_unit_system.length)
                self.assertTrue(max_x >= 9.5 | generic_unit_system.length)
                self.assertTrue(max_x <= 19.5 | generic_unit_system.length)
                            
        
        
        
        self.assertEquals(grids[0][0][0][0].rho,  0.0 | generic_unit_system.density)
        
        grids[0].rho = 0.2 | generic_unit_system.density
        
        rho1 = grids[1].rho
        rho0 = grids[0].rho
        for i in range(10):
            for j in range(10):
                for k in range(10):
                    self.assertEquals(rho1[0][0][0],  0.0 | generic_unit_system.density)
                    self.assertEquals(rho0[0][0][0],  0.2 | generic_unit_system.density)
                    
        instance.stop()
    
    
    def test3(self):
    
        for number_of_workers in range(2,6):
            instance = self.new_instance(MpiAmrVac, number_of_workers = number_of_workers)
            instance.set_parameters_filename(instance.default_parameters_filename)
            error = instance.initialize_code()
            self.assertEquals(error, 0)
            instance.parameters.mesh_length = (20.0, 20.0, 20.0) | generic_unit_system.length
            instance.parameters.mesh_size = (20, 20, 20) 
            instance.parameters.x_boundary_conditions = ("periodic","periodic")
            instance.parameters.y_boundary_conditions = ("periodic","periodic")
            instance.parameters.z_boundary_conditions = ("periodic","periodic")
            error = instance.commit_parameters()

            grids = list(instance.itergrids())
            
            self.assertEquals(len(grids), 8)
            
            
            for index, grid in enumerate(grids):
                position = grid.position
                #print instance.get_level_of_grid(index + 1)
                level = instance.get_level_of_grid(index + 1)
                self.assertEquals(level, 1)
                for i in range(3):
                    max_x = position[...,...,...,i].amax()
                    min_x = position[...,...,...,i].amin()
                    
                    self.assertTrue(min_x >= 0.5 | generic_unit_system.length)
                    self.assertTrue(min_x <= 10.5 | generic_unit_system.length)
                    self.assertTrue(max_x >= 9.5 | generic_unit_system.length)
                    self.assertTrue(max_x <= 19.5 | generic_unit_system.length)
                    self.assertEquals(max_x - min_x , 9.0 | generic_unit_system.length)
                    
            self.assertEquals(grids[0][0][0][0].rho,  0.0 | generic_unit_system.density)
            
            grids[0].rho = 0.2 | generic_unit_system.density
            
            rho1 = grids[4].rho
            rho0 = grids[0].rho
            for i in range(10):
                for j in range(10):
                    for k in range(10):
                        self.assertEquals(rho1[0][0][0],  0.0 | generic_unit_system.density)
                        
                        self.assertEquals(rho0[0][0][0],  0.2 | generic_unit_system.density)
                        
            instance.stop()
        
    def test4(self):
    
        instance = self.new_instance(MpiAmrVac, number_of_workers = 1)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.parameters.mesh_length = (10.0,10.0, 10.0) | generic_unit_system.length
        instance.parameters.mesh_size = (10, 10, 10)
        instance.parameters.maximum_number_of_grid_levels = 5
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
     
        instance.commit_parameters()

        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        
        levels = [instance.get_level_of_grid(i+1) for i in range(len(grids))]
        self.assertEquals(levels, [1] * 1)
        
        gamma = 5.0 / 3.0
        energy =  generic_unit_system.mass / (generic_unit_system.time**2 * generic_unit_system.length)
        grids_in_memory = []
        for grid in grids:
            grid.rho = 1.0  | generic_unit_system.density
            grid.energy = (0.1795 | energy)/ (gamma - 1)
    
        
        has_advanced = instance.refine_grid()
        grids = list(instance.itergrids())
        print len(grids)
        self.assertFalse(has_advanced)
        instance.stop()

    def test5(self):
        def fill_grids(grids):
            for grid in grids:
                firsthalf = grid.x > 5.0 | generic_unit_system.length
                secondhalf = numpy.logical_not(firsthalf)
                if(numpy.any(firsthalf)):
                    grid[firsthalf].rho = 1.0  | generic_unit_system.density
        
        instance = self.new_instance(MpiAmrVac, number_of_workers = 1)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.parameters.mesh_length = (10.0,10.0, 10.0) | generic_unit_system.length
        instance.parameters.mesh_size = (10, 10, 10)
        instance.parameters.maximum_number_of_grid_levels = 4
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
     
        instance.commit_parameters()

        n = 1
        for i in range(3):
            expected_length = n ** 3
            grids = list(instance.itergrids())
            levels = [instance.get_level_of_grid(grid+1) for grid in range(len(grids))]
            self.assertEquals(len(grids), expected_length)
            
            self.assertEquals(levels, [i+1] * expected_length)
            
            n *= 2
            
            fill_grids(grids)
            has_advanced = instance.refine_grid()
            self.assertTrue(has_advanced)
        
        has_advanced = instance.refine_grid()
        self.assertFalse(has_advanced)
        
        
        instance.stop()
        
    
    def test6(self):
        instance = self.new_instance(MpiAmrVac, number_of_workers = 1)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        instance.parameters.mesh_length = (10.0,10.0, 10.0) | generic_unit_system.length
        instance.parameters.mesh_size = (10, 10, 10)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        
        self.assertEquals(instance.parameters.maximum_number_of_grid_levels, 3)
        instance.parameters.maximum_number_of_grid_levels = 2
        self.assertEquals(instance.parameters.maximum_number_of_grid_levels, 2)
        
        self.assertEquals(instance.parameters.entropy_type , 'nul')
        instance.parameters.entropy_type = 'powell'
        self.assertEquals(instance.parameters.entropy_type , 'powell')
        
        self.assertEquals(instance.parameters.time_integration_procedure , 'twostep')
        instance.parameters.time_integration_procedure = 'onestep'
        self.assertEquals(instance.parameters.time_integration_procedure , 'onestep')
        
        self.assertEquals(instance.parameters.spatial_discretization_method , 'tvdmu')
        instance.parameters.spatial_discretization_method = 'tvdlf'
        self.assertEquals(instance.parameters.spatial_discretization_method , 'tvdlf')
        
        
        
        self.assertEquals(instance.parameters.predictor_step_discretization_method , 'tvdmu')
        instance.parameters.predictor_step_discretization_method = 'hancock'
        self.assertEquals(instance.parameters.predictor_step_discretization_method , 'hancock')
        
        instance.commit_parameters()
        
        self.assertEquals(instance.parameters.entropy_type , 'powell')
        self.assertEquals(instance.parameters.maximum_number_of_grid_levels, 2)
        self.assertEquals(instance.parameters.time_integration_procedure , 'onestep')
        self.assertEquals(instance.parameters.spatial_discretization_method , 'tvdlf')
        self.assertEquals(instance.parameters.predictor_step_discretization_method , 'hancock')
        
    def test7(self):
        instance=self.new_instance(MpiAmrVac)
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.initialize_code()
        instance.parameters.mesh_length = (10.0,10.0, 10.0) | generic_unit_system.length
        instance.parameters.mesh_size = (10,10,10)
        instance.parameters.maximum_number_of_grid_levels = 1
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        
        instance.commit_parameters()
    
        grid = datamodel.Grid(10,10,10)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy =  1.0 | generic_unit_system.energy_density
        
        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        igrid = grids[0]
        channel = grid.new_channel_to(igrid)
        channel.copy()
        
                   
        instance.evolve_model(0.1 | generic_unit_system.time)
        
        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.evolve_model(0.3 | generic_unit_system.time)
        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.stop()
        
    def test8(self):
        instance=self.new_instance(MpiAmrVac)
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.initialize_code()
        instance.parameters.mesh_size = (10,10,10)
        instance.parameters.maximum_number_of_grid_levels = 1
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 1.0 | generic_unit_system.length
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        
        name, error = instance.get_typeghostfill()
        self.assertEquals(name, 'linear')
        instance.commit_parameters()
        
        name, error = instance.get_typeghostfill()
        self.assertEquals(name, 'linear')
        grid = datamodel.Grid(10,10,10)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy =  1.0 | generic_unit_system.energy_density
        
        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        igrid = grids[0]
        channel = grid.new_channel_to(igrid)
        channel.copy()
        
        name, error = instance.get_typeghostfill()
        self.assertEquals(name, 'linear')
    
        self.assertEquals((50,50,50), instance.acceleration_grid.shape)
        
        acc_grid = datamodel.Grid(50,50,50)
        acceleration = 1 | generic_unit_system.acceleration
        acc_grid.ax = acceleration
        acc_grid.ay = acceleration
        acc_grid.az = acceleration
        #self.assertEquals(acc_grid.acceleration[0][0][0], ( 1,1,1) | generic_unit_system.acceleration)
        channel = acc_grid.new_channel_to(instance.acceleration_grid)
        channel.copy()
        
        name, error = instance.get_typeghostfill()
        self.assertEquals(name, 'linear')
        
        result = instance.initialize_grid()
                   
        name, error = instance.get_typeghostfill()
        self.assertEquals(name, 'linear')
    
        instance.evolve_model(0.1 | generic_unit_system.time)
        
        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho);
        self.assertAlmostRelativeEquals(igrid.rhovx, 0.1 * 1.0 * 0.1 | generic_unit_system.momentum_density);
        self.assertAlmostRelativeEquals(igrid.rhovy, 0.1 * 1.0 * 0.1 | generic_unit_system.momentum_density);
        self.assertAlmostRelativeEquals(igrid.rhovz, 0.1 * 1.0 * 0.1 | generic_unit_system.momentum_density);

        instance.evolve_model(0.3 | generic_unit_system.time)
        print instance.model_time
        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho);
        self.assertAlmostRelativeEquals(igrid.rhovx, grid.rho *  instance.model_time * acceleration,2);
        self.assertAlmostRelativeEquals(igrid.rhovy, grid.rho *  instance.model_time * acceleration,2);
        self.assertAlmostRelativeEquals(igrid.rhovz, grid.rho *  instance.model_time * acceleration,2);
        instance.stop()
        
    
    def test9(self):
        instance=self.new_instance(MpiAmrVac, mode="2d")
        instance.set_parameters_filename(instance.default_parameters_filename)
        instance.initialize_code()
        instance.parameters.mesh_length = (10.0,10.0, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (10,10,1)
        instance.parameters.maximum_number_of_grid_levels = 1
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        
        instance.commit_parameters()
    
        grid = datamodel.Grid(10,10,1)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy =  1.0 | generic_unit_system.energy_density
        
        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        igrid = grids[0]
        channel = grid.new_channel_to(igrid)
        channel.copy()
        
                   
        instance.evolve_model(0.1 | generic_unit_system.time)
        
        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.evolve_model(0.3 | generic_unit_system.time)
        grids = list(instance.itergrids())
        self.assertEquals(len(grids), 1)
        igrid = grids[0]
        self.assertAlmostRelativeEquals(igrid.rho, grid.rho)
        
        instance.stop()
