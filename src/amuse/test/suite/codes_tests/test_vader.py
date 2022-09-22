import numpy as np

from amuse.community import *
from amuse.test.amusetest import TestWithMPI
from amuse.units import units, constants

from amuse.community.vader.interface import VaderInterface, Vader

default_options = dict()

class TestVaderInterface(TestWithMPI):

    def test1(self):
        instance = self.new_instance_of_an_optional_code(
            VaderInterface, **default_options)
        instance.initialize_code()
        instance.set_nUserOut(nUserOut=1)
        instance.initialize_flat_grid(n=99, linear=True, rmin=1., rmax=100.,
            vphi=1.23)

        index = 25

        self.assertEqual(1.23, 
            instance.get_rotational_velocity_of_index([index])['vphi'][0])
        self.assertEqual(float(index)+1.5, 
            instance.get_position_of_index([index])['r'][0])

        instance.cleanup_code()
        instance.stop()

    def test2(self):
        instance = self.new_instance_of_an_optional_code(
            VaderInterface, **default_options)
        instance.initialize_code()
        instance.set_nUserOut(nUserOut=1)
        instance.initialize_keplerian_grid(n=100, linear=False, rmin=1., rmax=100.,
            m=1.)

        G = 6.67384e-8 # From VADER source
        index = 50

        r = instance.get_position_of_index([index])['r'][0]

        v_kepler = (G/r)**0.5

        # Should be exactly equal, but different roundoff in python/c could ruin this
        #self.assertEqual(v_kepler, 
        #    instance.get_rotational_velocity_of_index([index])['vphi'][0])
        self.assertAlmostEqual(v_kepler, 
            instance.get_rotational_velocity_of_index([index])['vphi'][0])
        self.assertAlmostEqual(v_kepler*v_kepler,
            instance.get_gravitational_potential_of_index([index])['psi_grav'][0])

        instance.cleanup_code()
        instance.stop()

    def test3(self):
        instance = self.new_instance_of_an_optional_code(
            VaderInterface, **default_options)
        instance.initialize_code()
        instance.set_nUserOut(nUserOut=1)
        instance.initialize_flat_grid(n=100, linear=False, rmin=1., rmax=100.,
            vphi=1.23)

        instance.set_grid_column_density(np.arange(100), np.logspace(0, -2, num=100))
        instance.set_grid_pressure(np.arange(100), 10.*np.logspace(0, -2, num=100))

        index = 75

        self.assertEqual(10., 
            instance.get_grid_pressure([index])['pressure'][0] / \
            instance.get_grid_column_density([index])['sigma'][0])

        instance.cleanup_code()
        instance.stop()

    def test4(self):
        instance = self.new_instance_of_an_optional_code(
            VaderInterface, **default_options)
        instance.initialize_code()
        instance.set_nUserOut(nUserOut=1)
        instance.initialize_flat_grid(n=100, linear=False, rmin=1., rmax=100.,
            vphi=1.23)

        instance.set_number_of_user_parameters(5)
        instance.set_parameter(0, 1.)
        instance.set_parameter(4, 2.)

        self.assertEqual(1., instance.get_parameter(0)['param'])
        self.assertEqual(2., instance.get_parameter(4)['param'])

        instance.cleanup_code()
        instance.stop()

class TestVader(TestWithMPI):

    def test1(self):
        instance = self.new_instance_of_an_optional_code(
            Vader, **default_options)
        instance.initialize_flat_grid(100, False, 1.|units.AU, 100.|units.AU,
            1.|units.kms)

        self.assertEqual(1.|units.kms, instance.grid.rotational_velocity[0])

        instance.update_flat_grid(2.|units.kms)

        self.assertEqual(2.|units.kms, instance.grid.rotational_velocity[0])

        instance.stop()

    def test2(self):
        instance = self.new_instance_of_an_optional_code(
            Vader, **default_options)
        instance.initialize_keplerian_grid(100, False, 1.|units.AU, 100.|units.AU,
            1.|units.MSun)

        G = 6.67384e-8 | units.cm**3/units.g/units.s**2 # From VADER source
        index = 25

        v_kepler = (G*(1.|units.MSun)/instance.grid.r[index])**0.5

        self.assertEqual(v_kepler, instance.grid.rotational_velocity[index])

        instance.update_keplerian_grid(2.|units.MSun)

        v_kepler = (G*(2.|units.MSun)/instance.grid.r[index])**0.5

        self.assertEqual(v_kepler, instance.grid.rotational_velocity[index])
        self.assertEqual((v_kepler*v_kepler).value_in((units.cm/units.s)**2),
            -instance.grid.gravitational_potential[index].value_in(
            (units.cm/units.s)**2))

        instance.stop()

    def test3(self):
        instance = self.new_instance_of_an_optional_code(
            Vader, **default_options)
        instance.initialize_keplerian_grid(100, False, 1.|units.AU, 100.|units.AU,
            1.|units.MSun)

        instance.parameters.alpha = 1e-3
        instance.parameters.gamma = 5./3.
        instance.parameters.delta = 1. | (units.cm/units.s)**2

        self.assertEqual(0.001, instance.parameters.alpha)
        self.assertEqual(5./3., instance.parameters.gamma)
        self.assertEqual(1., 
            instance.parameters.delta.value_in((units.cm/units.s)**2))

        col = np.logspace(0, -2, num=100) | units.g/units.cm**2
        pres = col * constants.kB*(100.|units.K)/(constants.u)

        instance.grid.column_density = col            
        instance.grid.pressure = pres

        mass0 = instance.mass

        instance.evolve_model( 1.|units.kyr )

        mass1 = instance.mass

        self.assertAlmostEqual(mass0.value_in(units.MSun),
                               mass1.value_in(units.MSun))

        instance.stop()

    def test4(self):
        instance = self.new_instance_of_an_optional_code(
            Vader, **default_options)
        instance.initialize_keplerian_grid(100, False, 1.|units.AU, 100.|units.AU,
            1.|units.MSun)

        instance.parameters.alpha = 1e-3
        instance.parameters.inner_pressure_boundary_mass_flux = \
            -1e-10 | units.MSun/units.yr
        instance.parameters.outer_pressure_boundary_mass_flux = \
             1e-12 | units.MSun/units.yr

        col = np.logspace(0, -2, num=100) | units.g/units.cm**2
        pres = col * constants.kB*(100.|units.K)/(constants.u)

        instance.grid.column_density = col            
        instance.grid.pressure = pres

        mass0 = instance.mass

        instance.evolve_model( 1.|units.kyr )

        mass1 = instance.mass

        self.assertAlmostEqual(mass0.value_in(units.MSun),
                               (mass1 - \
                                instance.inner_boundary_mass_out - \
                                instance.outer_boundary_mass_out
                               ).value_in(units.MSun))

        instance.stop()
