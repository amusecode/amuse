from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import mmcInterface
from .interface import mmc

class mmcInterfaceTests(TestWithMPI):
    
    def test1(self):
        #instance = mmcInterface(redirection="file", redirect_file = "junk.txt")
        instance = mmcInterface(redirection="null")
        instance.set_mmc_data_directory(instance.data_directory)
        instance.nonstandard_init()
        instance.set_irun(10)
        instance.set_nt(1000)
        instance.set_istart(1)
        instance.set_imodel(2)#plummer
        while (1):
            res, err =  instance.initial_run()
            if (res<0): break

        instance.set_istart(2)
        print "restart"
        s = raw_input()
        instance.set_time(1.9)
        instance.nonstandard_init()
        while (1):
            res, err =  instance.initial_run()
            print instance.get_time()
            print instance.get_state(1)
            if (res<0): break
            
        instance.stop()
    
