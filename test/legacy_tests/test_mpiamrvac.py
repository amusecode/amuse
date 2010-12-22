from amuse.legacy import *
from amuse.test.amusetest import TestWithMPI

from amuse.legacy.mpiamrvac.interface import MpiAmrVacInterface
from amuse.legacy.mpiamrvac.interface import MpiAmrVac

import os

class MpiAmrVacInterfaceTests(TestWithMPI):
    
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

        

        
