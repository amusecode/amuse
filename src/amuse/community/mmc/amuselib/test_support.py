from amuse.community import *
from amuse.test.amusetest import TestWithMPI

from .interface import supportInterface
import numpy as np

class supportInterfaceTests(TestWithMPI):
    
    def test1(self):
        instance = supportInterface()
        print instance.add([1,1,1,1,1],[1,1,1,1,1]).sum
        instance.stop()

    def test2(self):
        instance = supportInterface(redirection='none')
        x = np.zeros(10)
        y = np.zeros(10)
        z = np.zeros(10)
        print instance.many_points_on_sphere(x,y,z)
        instance.stop()


    
    
