from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import mmcInterface
from .interface import mmc

class mmcInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = mmcInterface(redirection="file", redirect_file = "junk.txt")
        instance.set_mmc_data_directory(instance.data_directory)
        instance.nonstandard_init()
        instance.set_nt(1010)
        while (1):
            print instance.get_time()
            print instance.get_kinetic_energy()
            res, err =  instance.initial_run()
            print res
            if (res<0): break
            
        instance.stop()
    
