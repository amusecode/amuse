from amuse.legacy import *
from amuse.test.amusetest import TestWithMPI

from amuse.legacy.mpiamrvac.interface import MpiAmrVacInterface
from amuse.legacy.mpiamrvac.interface import MpiAmrVac

from amuse.support.units import generic_unit_system
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

class TestMpiAmrVac(TestWithMPI):
    
    def test1(self):
        instance = self.new_instance(MpiAmrVac)
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        instance.setup_mesh(20,20,20, 20.0, 20.0, 20.0)
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
            
            error = instance.set_boundary("periodic", "periodic", "periodic", "periodic", "periodic", "periodic")
            error = instance.commit_parameters()

            grids = list(instance.itergrids())
            
            self.assertEquals(len(grids), 8)
            
            
            for index, grid in enumerate(grids):
                position = grid.position
                print instance.get_level_of_grid(index)
                level, error = instance.get_level_of_grid(index)
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
        
