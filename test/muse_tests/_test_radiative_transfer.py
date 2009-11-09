import unittest
import os

from radiative_transfer.mpirad.muse_radiative_transfer import MPIRad
from radiative_transfer.mpirad.mpirad_interface import *


class TestRadiativeTransfer(unittest.TestCase):
   pass

class TestMPIRad(TestRadiativeTransfer):
    def test1(self):
        grid_x = 1000                             
        grid_y = 1000
        grid_z = 100
        grid = new_doubleArray(1);
        interface = MPIRad()
        interface.SetDimensions(3.0, 3.0, 3.0, 0.5, 0.0, -1.0)
        interface.SetObserver(0.0, 0.0, 1.5)
        interface.SetNumberOfPhotons(1000)
        interface.SetRandomSeed(31415)
        interface.SetAlbedo(1.0)
        interface.SetOutputFiles("rtshot.dat", "luminosities.dat")
        #scattering type: 0 = isotropic, 1 = gas (electron) scattering, 2 = dust scattering
        interface.SetScatteringType(1)
        n = 2
        data = new_doubleArray(5 * n)
        for x in range(0,n):
            doubleArray_setitem(data, 5 * x, 3000)
            doubleArray_setitem(data, 5 * x + 1, 4)
            doubleArray_setitem(data, 5 * x + 2, 0.6+0.1*n)
            doubleArray_setitem(data, 5 * x + 3, 0.2)
            doubleArray_setitem(data, 5 * x + 4, -0.4)
        interface.ComputeRT(data, 2, grid, grid_x, grid_y, grid_z)
        os.remove('rtshot0.dat')
        os.remove('luminosities0.dat')

        
        