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
        instance = self.new_instance(MpiAmrVacInterface) #, redirection="none")
        print instance.default_parameters_filename
        instance.set_parameters_filename(instance.default_parameters_filename)
        error = instance.initialize_code()
        self.assertEquals(error, 0)
        
        error = instance.commit_parameters()
        self.assertEquals(error, 0)
        
        print instance.get_mesh_size(1)
        if True:
            return
        error = instance.initialize_grid()
        self.assertEquals(error, 0)
        print instance.get_mesh_size(1)
        print instance.get_position_of_index(0,0,0,1)
        print instance.get_position_of_index(11,11,11,1)
        print instance.get_position_of_index(0,0,0,2)
        instance.stop()
        self.assertTrue(False)        
        
