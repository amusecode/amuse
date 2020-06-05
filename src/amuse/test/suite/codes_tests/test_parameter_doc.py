import os
import sys

from amuse.test import amusetest
from amuse.community.bhtree.interface import BHTreeInterface, BHTree
from amuse.community.hermite.interface import HermiteInterface, Hermite
from amuse.community.phigrape.interface import PhiGRAPEInterface, PhiGRAPE
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel
class TestParameterDoc(amusetest.TestCase):

    def test_bhtree(self):

        convert_nbody = nbody_system.nbody_to_si(1.0 | units.kg, 1.0 | units.km)

        bhtree = BHTree(convert_nbody)

        bhtree.parameters.epsilon_squared = 10 | units.km**2
        bhtree.parameters.timestep = 1.0 | units.s
        bhtree.parameters.opening_angle = 0.1

        docstring =  bhtree.parameters.__doc__

        self.assertTrue("smoothing parameter for gravity calculations (default value:125000.0 m**2)" in docstring)

        parameter_str_method_output = str(bhtree.parameters)

        self.assertTrue("epsilon_squared: 10000000.0 m**2" in parameter_str_method_output)
        self.assertTrue("timestep: 1.0 s" in parameter_str_method_output)
        self.assertTrue("opening_angle: 0.1" in parameter_str_method_output)

    def test_hermite(self):

        convert_nbody = nbody_system.nbody_to_si(1.0 | units.kg, 1.0 | units.km)

        hermite = Hermite(convert_nbody)

        hermite.parameters.epsilon_squared = 10 | units.km**2

        docstring =  hermite.parameters.__doc__
        self.assertTrue("smoothing parameter for gravity calculations (default value:0.0 m**2)" in docstring)

        parameter_str_method_output = str(hermite.parameters)
        self.assertTrue("epsilon_squared: 10000000.0 m**2" in parameter_str_method_output)


    def test_phigrape(self):

        convert_nbody = nbody_system.nbody_to_si(1.0 | units.kg, 1.0 | units.km)

        phigrape = PhiGRAPE(convert_nbody)
        phigrape.parameters.epsilon_squared = 10 | units.km**2

        docstring =  phigrape.parameters.__doc__
        print(docstring)
        self.assertTrue("smoothing parameter for gravity calculations (default value:0.0 m**2)" in docstring)
        self.assertTrue("timestep parameter (default value:0.02" in docstring)
        self.assertTrue("parameter to determine the initial timestep (default value:0.01" in docstring)

        parameter_str_method_output = str(phigrape.parameters)
        self.assertTrue("epsilon_squared: 10000000.0 m**2" in parameter_str_method_output)
        self.assertTrue("timestep_parameter: 0.0" in parameter_str_method_output)
        self.assertTrue("initial_timestep_parameter: 0.0" in parameter_str_method_output)
