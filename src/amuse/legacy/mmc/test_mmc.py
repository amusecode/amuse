from amuse.legacy import *
from amuse.test.amusetest import TestWithMPI

from .interface import mmcInterface
from .interface import mmc

class mmcInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = mmcInterface(redirection="null")
        print "enter"
        instance.nonstandard_init()
        print "leave"
        print "leave"
        print "leave"
        print "leave"
        print instance.get_time()
        instance.stop()
    
