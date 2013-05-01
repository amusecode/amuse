import os
import sys
import numpy
import math

from amuse.test.amusetest import TestWithMPI
from amuse.community.athena.interface import AthenaInterface, Athena
from amuse.units.quantities import VectorQuantity
from amuse.units import generic_unit_system
from amuse.units import units
from amuse.units import generic_unit_converter
from amuse import datamodel

class TestAthenaInterface(TestWithMPI):
    
    def test0(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.stop()
        
    def test1(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.par_seti("test", "testname", "%d", 10, "a test parameter")
        x = instance.par_geti("test", "testname")
        
        self.assertEquals(x, 10)
        
        instance.stop()
        
    def test2(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.par_setd("test", "test2", "%.15e", 1.123, "a test parameter")
        x = instance.par_getd("test", "test2")
        
        self.assertEquals(x, 1.123)
        instance.stop()
        
        
    def test3(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.setup_mesh(5, 1, 1, 1.0, 0.0, 0.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        x = instance.par_geti("domain1", "Nx1")
        self.assertEquals(x, 5)
        
        error = instance.commit_parameters()
        self.assertEquals(error, 0)
        
        number_of_grids, error = instance.get_number_of_grids()
        self.assertEquals(error, 0)
        self.assertEquals(number_of_grids, 1)

        x,y,z,error = instance.get_position_of_index(0,0,0,1)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEquals(x, 0.1)
        
        x,y,z,error = instance.get_position_of_index(1,0,0,1)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEquals(x, 0.3)
        
        x,y,z,error = instance.get_position_of_index(2,0,0,1)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEquals(x, 0.5)
        
        x,y,z,error = instance.get_position_of_index(3,0,0,1)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEquals(x, 0.7)
        
        x,y,z,error = instance.get_position_of_index(4,0,0,1)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEquals(x, 0.9)
        
        x,y,z,error = instance.get_position_of_index(5,0,0,1)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEquals(x, 1.1)
        
        instance.stop()
        
        
    def test4(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.setup_mesh(10, 20, 40, 1.0, 1.0, 1.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        imin, imax, jmin, jmax, kmin, kmax = instance.get_index_range_inclusive()
        
        x,y,z, error= instance.get_position_of_index(2,2,2)
        self.assertEquals(error, -1)
        
        result = instance.commit_parameters()
        self.assertEquals(result, 0)
        
        x,y,z, error= instance.get_position_of_index(0,0,0)
        self.assertEquals(error, 0)
        print x,y,z
        self.assertAlmostRelativeEquals(0.05, x)
        self.assertAlmostRelativeEquals(0.025, y)
        self.assertAlmostRelativeEquals(0.0125, z)
        
        x,y,z, error= instance.get_position_of_index(10,20,40)
        self.assertEquals(error, 0)
        print x,y,z
        self.assertAlmostRelativeEquals(1.05, x)
        self.assertAlmostRelativeEquals(1.025, y)
        self.assertAlmostRelativeEquals(1.0125, z)

        instance.stop()
        
    



    def test5(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(5, 5, 5, 1.0, 1.0, 1.0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.4)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        result = instance.commit_parameters()
        self.assertEquals(result, 0)
        
        time, error = instance.get_time()
        self.assertEquals(error,0)
        self.assertEquals(time, 0.0)
        
        error = instance.set_grid_state(1,1,1,0.1, 0.2, 0.3, 0.4, 0.5)
        self.assertEquals(error, 0)
        
        rho, rhovx, rhovy, rhovz, energy, error = instance.get_grid_state(1,1,1)
        
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0.1)
        self.assertEquals(rhovx, 0.2)
        self.assertEquals(rhovy, 0.3)
        self.assertEquals(rhovz, 0.4)
        self.assertEquals(energy, 0.5)
        
        rhovx, rhovy, rhovz, error = instance.get_grid_momentum_density(1,1,1)
        
        self.assertEquals(error, 0)
        self.assertEquals(rhovx, 0.2)
        self.assertEquals(rhovy, 0.3)
        self.assertEquals(rhovz, 0.4)
        
        rho, error = instance.get_grid_density(1,1,1)
        
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0.1)
        
        energy, error = instance.get_grid_energy_density(1,1,1)
        
        self.assertEquals(error, 0)
        self.assertEquals(energy, 0.5)
        
        
        rho, rhovx, rhovy, rhovz, energy, error = instance.get_grid_state([1],[1],[1])
        self.assertEquals(error[0], 0)
        self.assertEquals(rho[0], 0.1)
        error = instance.initialize_grid()
        self.assertEquals(error, 0)
        
        timestep, error =  instance.get_timestep()
        self.assertEquals(error, 0)
        
        
        instance.stop()
    
    def test5a(self):
        instance=self.new_instance(AthenaInterface, mode="mhd")
        instance.initialize_code()
        instance.setup_mesh(5, 5, 5, 1.0, 1.0, 1.0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.4)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        result = instance.commit_parameters()
        self.assertEquals(result, 0)
        
        time, error = instance.get_time()
        self.assertEquals(error,0)
        self.assertEquals(time, 0.0)
        
        error = instance.set_grid_magnetic_field(1,1,1,0.1, 0.2, 0.3)
        self.assertEquals(error, 0)
        
        B1i, B2i, B3i, error = instance.get_grid_magnetic_field(1,1,1)
        
        self.assertEquals(error, 0)
        self.assertEquals(B1i, 0.1)
        self.assertEquals(B2i, 0.2)
        self.assertEquals(B3i, 0.3)
        
        instance.stop()
    
    

    def test7(self):
        results = []
        for x in range(1,5):
            instance=self.new_instance(AthenaInterface, number_of_workers=x)
            instance.initialize_code()
            instance.setup_mesh(128,1,1,1.0,0,0)
            instance.set_gamma(1.6666666666666667)
            instance.set_courant_friedrichs_lewy_number(0.8)
            instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
            result = instance.commit_parameters()
            self.assertEquals(result, 0)
            
            nghost, error = instance.get_nghost()
            self.assertEquals(4, nghost)
            instance.set_grid_state(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128),0.1, 0.2, 0.3, 0.4, 0.5)
            error = instance.initialize_grid()
            self.assertEquals(error, 0)
            
            result = instance.get_grid_state(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128))
            results.append(list(result))
            
            
            instance.stop()
        
        for x in range(128):
            for y in range(6):
                self.assertEquals(results[1][y][x], results[0][y][x])
                self.assertEquals(results[2][y][x], results[0][y][x])
                self.assertEquals(results[3][y][x], results[0][y][x])
    
    

    def test8(self):
        instance=self.new_instance(AthenaInterface, number_of_workers=1)
        instance.initialize_code()
        instance.setup_mesh(128,1,1,1.0,0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        result = instance.commit_parameters()
        self.assertEquals(result, 0)
        
        instance.fill_grid_linearwave_1d(0, 1e-06, 0.0, 1)
        error = instance.initialize_grid()
        self.assertEquals(error, 0)
        
        print instance.get_grid_state(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128))
        
        timestep, error =  instance.get_timestep()
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(timestep,  0.006249991, 5)
        
        rho0, rhovx0, rhovy0, rhovz0, energy0, error0  = instance.get_grid_state(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128))
        instance.evolve_model(5.0)
        rho, rhovx, rhovy, rhovz, energy, error  = instance.get_grid_state(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128))
        
        error_rho = numpy.sum(numpy.abs(rho - rho0))
        error_rhovx =  numpy.sum(numpy.abs(rhovx - rhovx0))
        error_rhovy =  numpy.sum(numpy.abs(rhovy - rhovy0))
        error_rhovz =  numpy.sum(numpy.abs(rhovz - rhovz0))
        error_energy =  numpy.sum(numpy.abs(energy - energy0))
        self.assertAlmostRelativeEquals(error_rho / 128.0, 1.877334e-09, 6)
        self.assertAlmostRelativeEquals(error_rhovx / 128.0, 1.877334e-09, 6)
        self.assertAlmostRelativeEquals(error_energy / 128.0, 2.816001e-09, 6)
        
        instance.stop()
    
    

    def test9(self):
        results = []
        for x in range(1,5):
            instance=self.new_instance(AthenaInterface, number_of_workers=x)
            instance.initialize_code()
            instance.setup_mesh(128,1,1,1.0,0,0)
            instance.set_gamma(1.6666666666666667)
            instance.set_courant_friedrichs_lewy_number(0.8)
            instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
            result = instance.commit_parameters()
            self.assertEquals(result, 0)
            
            nghost, error = instance.get_nghost()
            self.assertEquals(4, nghost)
            instance.fill_grid_linearwave_1d(0, 1e-06, 0.0, 1)
            error = instance.initialize_grid()
            self.assertEquals(error, 0)
            instance.evolve_model(5.0)
        
            result = instance.get_grid_state(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128))
            results.append(list(result))
            
            
            instance.stop()
        
        for x in range(128):
            for y in range(6):
                self.assertEquals(results[1][y][x], results[0][y][x])
                self.assertEquals(results[2][y][x], results[0][y][x])
                self.assertEquals(results[3][y][x], results[0][y][x])
    
    

    def test10(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.setup_mesh(5, 5, 5, 1.0, 1.0, 1.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        instance.initialize_grid()
        x,y,z,error = instance.get_position_of_index([0,1,2,3,4],[0,0,0,0,0],[0,0,0,0,0])
        
        for x0, x1 in zip(x, [0.1, 0.3, 0.5, 0.7, 0.9]):
            self.assertAlmostRelativeEqual(x0, x1)
        
        for y0, y1 in zip(y, [0.1, 0.1, 0.1, 0.1, 0.1]):
            self.assertAlmostRelativeEqual(y0, y1)
        
        i,j,k,error = instance.get_index_of_position(0.3, 0.1, 0.1)
        
        print i,j,k
        
        self.assertAlmostRelativeEqual(i, 1)
        self.assertAlmostRelativeEqual(j, 0)
        self.assertAlmostRelativeEqual(k, 0)
        
        i,j,k,error = instance.get_index_of_position(0.39, 0.1, 0.1)
        self.assertAlmostRelativeEqual(i, 1.0)
        self.assertAlmostRelativeEqual(j, 0)
        self.assertAlmostRelativeEqual(k, 0)
        

        i,j,k,error = instance.get_index_of_position(0.4, 0.1, 0.1)
        self.assertAlmostRelativeEqual(i, 2.0)
        self.assertAlmostRelativeEqual(j, 0)
        self.assertAlmostRelativeEqual(k, 0)
        
        x,y,z,error = instance.get_position_of_index(-1,0,0)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(x, -0.1)
        x,y,z,error = instance.get_position_of_index(5,0,0)
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(x, 1.1)
        instance.stop()
    
    

    def test6(self):
        instance=self.new_instance(AthenaInterface, number_of_workers = 5)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.setup_mesh(100, 200, 400, 10.0, 10.0, 10.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        imin, imax, jmin, jmax, kmin, kmax = instance.get_index_range_inclusive()
        
        x,y,z, error= instance.get_position_of_index(2,2,2)
        self.assertEquals(error, -1)
        
        result = instance.commit_parameters()
        self.assertEquals(result, 0)
        
        x,y,z, error= instance.get_position_of_index(0,0,0)
        self.assertEquals(error, 0)
        print x,y,z
        self.assertAlmostRelativeEquals(0.05, x)
        self.assertAlmostRelativeEquals(0.025, y)
        self.assertAlmostRelativeEquals(0.0125, z)
        
        
        
        x,y,z, error= instance.get_position_of_index(100,200,400)
        self.assertEquals(error, 0)
        print x,y,z
        self.assertAlmostRelativeEquals(10.05, x)
        self.assertAlmostRelativeEquals(10.025, y)
        self.assertAlmostRelativeEquals(10.0125, z)
        
        
        instance.stop()
    
    

    def test11(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.setup_mesh(10, 20, 40, 1.0, 1.0, 1.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        imin, imax, jmin, jmax, kmin, kmax = instance.get_index_range_inclusive()
        
        self.assertEquals(imin, 0)
        self.assertEquals(jmin, 0)
        self.assertEquals(kmin, 0)
        self.assertEquals(imax, 9)
        self.assertEquals(jmax, 19)
        self.assertEquals(kmax, 39)
        
        imin, imax, jmin, jmax, kmin, kmax = instance.get_index_range_for_potential()
        
        self.assertEquals(imin, -1)
        self.assertEquals(jmin, -1)
        self.assertEquals(kmin, -1)
        self.assertEquals(imax, 10)
        self.assertEquals(jmax, 20)
        self.assertEquals(kmax, 40)
        instance.stop()
    
    

    def test12(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.setup_mesh(10, 20, 40, 1.0, 1.0, 1.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        instance.commit_parameters()
        instance.initialize_grid()
        x,y,z, error= instance.get_position_of_index(-1,-1,-1)
        self.assertEquals(error, 0)
        print x,y,z
        self.assertAlmostRelativeEquals(-0.05, x)
        self.assertAlmostRelativeEquals(-0.025, y)
        self.assertAlmostRelativeEquals(-0.0125, z)
        instance.stop()
    
    

    def test13(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.setup_mesh(2, 2, 2, 1.0, 1.0, 1.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        instance.commit_parameters()
        instance.initialize_grid()
        
        potential_along_one_axis = [-1.0,0.0,1.0,2.0]
        
        instance.set_potential(
            [-1,0,1,2],
            [0,0,0,0],
            [0,0,0,0], 
            potential_along_one_axis
        )
        got_potential,error = instance.get_potential(
            [-1,0,1,2],
            [0,0,0,0],
            [0,0,0,0])
            
        print got_potential, error
        for expected, actual in zip(potential_along_one_axis, got_potential):
            self.assertEquals(expected, actual)
        
        x,y,z,error = instance.get_position_of_index(
            [-1,0,1,2],
            [0,0,0,0],
            [0,0,0,0])
        print x,y,z, error
        for expected, actual in zip([-0.25,0.25,0.75,1.25], x):
            self.assertEquals(expected, actual)
        for expected, actual in zip([0.25,0.25,0.25,0.25], y):
            self.assertEquals(expected, actual)
        for expected, actual in zip([0.25,0.25,0.25,0.25], z):
            self.assertEquals(expected, actual)
            
        potential, error = instance.get_interpolated_gravitational_potential(0, 0.25, 0.25)
        print potential, error
        self.assertEquals(error, 0)
        self.assertEquals(potential, -0.5)
        potential, error = instance.get_interpolated_gravitational_potential(0.75, 0.5, 0.25)
        print potential, error
        self.assertEquals(error, 0)
        self.assertEquals(potential, 0.5)
        potential, error = instance.get_interpolated_gravitational_potential(0.75, 0.25, 0.5)
        print potential, error
        self.assertEquals(error, 0)
        self.assertEquals(potential, 0.5)
        potential, error = instance.get_interpolated_gravitational_potential(0.75, 0.25, 0.0)
        print potential, error
        self.assertEquals(error, 0)
        self.assertEquals(potential, 0.5)
        instance.stop()
        
    def test14(self):
        instance=self.new_instance(AthenaInterface, mode="scalar")
        instance.initialize_code()
        instance.setup_mesh(5, 5, 5, 1.0, 1.0, 1.0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.4)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        result = instance.commit_parameters()
        self.assertEquals(result, 0)
        
        error = instance.set_grid_scalar(1,1,1,0.45)
        self.assertEquals(error, 0)
        
        scalar, error = instance.get_grid_scalar(1,1,1)
        
        self.assertEquals(error, 0)
        self.assertEquals(scalar, 0.45)
        
        scalar, error = instance.get_grid_scalar(1,1,2)
        
        self.assertEquals(error, 0)
        self.assertEquals(scalar, 0)
        
        instance.stop()
        
    def test15(self):
        results = []
        for x in range(1,6):
            instance=self.new_instance(AthenaInterface, number_of_workers=x)
            instance.initialize_code()
            instance.setup_mesh(100,1,1,100.0,0,0)
            instance.set_gamma(1.6666666666666667)
            instance.set_courant_friedrichs_lewy_number(0.8)
            instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
            result = instance.commit_parameters()
            
            for index in range(100):
                x,y,z,error = instance.get_position_of_index(index,0,0)
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(x, index + 0.5)
                
                i,j,k,error = instance.get_index_of_position(x,y,z)
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(i, index)
                
            
            instance.stop()
        
    def test16(self):
        for x in range(1,6):
            instance=self.new_instance(AthenaInterface, number_of_workers=x)
            instance.initialize_code()
            instance.setup_mesh(10,100,1,100.0,100.0,0)
            instance.set_gamma(1.6666666666666667)
            instance.set_courant_friedrichs_lewy_number(0.8)
            instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
            result = instance.commit_parameters()
            
            for index in range(100):
                x,y,z,error = instance.get_position_of_index(0,index,0)
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(y, index + 0.5)
                
                i,j,k,error = instance.get_index_of_position(x,y,z)
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(j, index)
                
            
            instance.stop()
    
    def test17(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(100,1,1,100.0,0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","periodic","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        minx, maxx, miny, maxy, minz, maxz, error = instance.get_boundary_index_range_inclusive(1,1)
        self.assertEquals(error, 0)
        self.assertEquals(minx, 0)
        self.assertEquals(maxx, 3)
        self.assertEquals(miny, 0)
        self.assertEquals(maxy, 0)
        self.assertEquals(minz, 0)
        self.assertEquals(maxz, 0)
        
        for i in range(2,7):
            minx, maxx, miny, maxy, minz, maxz, error = instance.get_boundary_index_range_inclusive(i,1)
            self.assertEquals(error, 0)
            self.assertEquals(minx, 0)
            self.assertEquals(maxx, 0)
            self.assertEquals(miny, 0)
            self.assertEquals(maxy, 0)
            self.assertEquals(minz, 0)
            self.assertEquals(maxz, 0)
    
    def test18(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(100,5,6,100.0,0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","periodic","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        minx, maxx, miny, maxy, minz, maxz, error = instance.get_boundary_index_range_inclusive(1,1)
        self.assertEquals(error, 0)
        self.assertEquals(minx, 0)
        self.assertEquals(maxx, 3)
        self.assertEquals(miny, 0)
        self.assertEquals(maxy, 4)
        self.assertEquals(minz, 0)
        self.assertEquals(maxz, 5)
        
        for i in range(2,7):
            minx, maxx, miny, maxy, minz, maxz, error = instance.get_boundary_index_range_inclusive(i,1)
            self.assertEquals(error, 0)
            self.assertEquals(minx, 0)
            self.assertEquals(maxx, 0)
            self.assertEquals(miny, 0)
            self.assertEquals(maxy, 0)
            self.assertEquals(minz, 0)
            self.assertEquals(maxz, 0)
    
    def test19(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(100,5,6,100.0,0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
        
        for i in range(1,7):
            minx, maxx, miny, maxy, minz, maxz, error = instance.get_boundary_index_range_inclusive(i,1)
            self.assertEquals(error, 0)
            self.assertEquals(minx, 0)
            self.assertEquals(miny, 0)
            self.assertEquals(minz, 0)
            if i == 1 or i == 2:
                self.assertEquals(maxx, 3)
                self.assertEquals(maxy, 4)
                self.assertEquals(maxz, 5)
            elif i == 3 or i == 4:
                self.assertEquals(maxx, 99+8)
                self.assertEquals(maxy, 3)
                self.assertEquals(maxz, 5)
            elif i == 5 or i == 6:
                self.assertEquals(maxx, 99+8)
                self.assertEquals(maxy, 4+8)
                self.assertEquals(maxz, 3)
    
    def test20(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(100,1,1,100.0,100.0,100.0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","periodic","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        for i in range(4):
            error = instance.set_boundary_state(
                i,0,0,       #  index
                1.0 * (i+1),         #  density
                2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                5.0 * (i+1),         #  energy
                1.0, 1.0     #  boundary + grid
            )
            self.assertEquals(error, 0)
            rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                i, 0, 0,
                1.0, 1.0
            )
            print rho, rhovx, rhovy, rhovz, rhoen, error 
            self.assertEquals(error, 0)
            self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
            self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
            self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
            self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
            self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    def test21(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(100,1,1,100.0,0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        for i in range(4):
            for j in [1,2]:
                error = instance.set_boundary_state(
                    i,0,0,       #  index
                    1.0 * (i+1),         #  density
                    2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                    5.0 * (i+1),         #  energy
                    j, 1.0     #  boundary + grid
                )
                self.assertEquals(error, 0)
                rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                    i, 0, 0,
                    j, 1.0
                )
                print j
                self.assertEquals(error, 0)
                
                self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
                
    
    
    def test22(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(5,6,7,100.0,100.0,100.0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
        x1range = (4,6,7)
        x2range = (5,4,7)
        x3range = (5,6,4)
    
        for xrange, j in zip([x1range, x1range, x2range, x2range, x3range, x3range], [1,2,3,4,5,6]):
            for i0 in range(xrange[0]):
                for j0 in range(xrange[1]):
                    for k0 in range(xrange[2]):
                        i = (i0 * (xrange[2] * xrange[1])) + (j0 * xrange[2]) + k0
                        
                        error = instance.set_boundary_state(
                            i0, j0, k0,       #  index
                            1.0 * (i+1),         #  density
                            2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                            5.0 * (i+1),         #  energy
                            j, 1.0     #  boundary + grid
                        )
                        self.assertEquals(error, 0)
                        rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                            i0, j0, k0,       #  index
                            j, 1.0
                        )
                        self.assertEquals(error, 0)
                        
                        self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
                
    def test24(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(100,5,1,100.0,0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
        
        for i in range(1,7):
            minx, maxx, miny, maxy, minz, maxz, error = instance.get_boundary_index_range_inclusive(i,1)
            self.assertEquals(error, 0)
            self.assertEquals(minx, 0)
            self.assertEquals(miny, 0)
            self.assertEquals(minz, 0)
            if i == 1 or i == 2:
                self.assertEquals(maxx, 3)
                self.assertEquals(maxy, 4)
                self.assertEquals(maxz, 0)
            elif i == 3 or i == 4:
                self.assertEquals(maxx, 99+8)
                self.assertEquals(maxy, 3)
                self.assertEquals(maxz, 0)
            elif i == 5 or i == 6:
                self.assertEquals(maxx, 99+8)
                self.assertEquals(maxy, 4 +8)
                self.assertEquals(maxz, 3)
        
    def test25(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(100,1,1,100.0,0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        dx = 1.0
        for i in range(4):
            x,y,z,error = instance.get_boundary_position_of_index(
                i,0,0, 
                1, 1 
            )
            self.assertEquals(error, 0)
            self.assertAlmostRelativeEquals(x, (0.5 * dx) - ((4 -i)*dx))
            self.assertAlmostRelativeEquals(y, 0.0)
            self.assertAlmostRelativeEquals(z, 0.0)
                
    
    def test26(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(100,1,1,100.0,0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        dx = 1.0
        for i in range(4):
            x,y,z,error = instance.get_boundary_position_of_index(
                i,0,0, 
                2, 1 
            )
            self.assertEquals(error, 0)
            self.assertAlmostRelativeEquals(x, 100.0 + (0.5 * dx) + (i * dx))
            self.assertAlmostRelativeEquals(y, 0.0)
            self.assertAlmostRelativeEquals(z, 0.0)
    
    def test27(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(100,5,1,100.0,100.0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","periodic","periodic")
        instance.commit_parameters()
        
        dx = 1.0
        dy = 100.0 / 5.0
        for i in range(4):
            for j in range(5):
                x,y,z,error = instance.get_boundary_position_of_index(
                    i, j, 1, 
                    2, 1 
                )
                print y, j, (0.5 * dy) - ((4 - j) * dy)
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(x, 100.0 + (0.5 * dx) + (i * dx))
                self.assertAlmostRelativeEquals(y, (0.5 * dy) + (j * dy))
                self.assertAlmostRelativeEquals(z, 0.0)
        
        for i in range(100 + 8):
            for j in range(4):
                x,y,z,error = instance.get_boundary_position_of_index(
                    i, j, 1, 
                    3, 1 
                )
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-4) * dx))
                self.assertAlmostRelativeEquals(y, ((0.5 * dy) - ((4-j) * dy)))
                self.assertAlmostRelativeEquals(z, 0.0)
                
                
                x,y,z,error = instance.get_boundary_position_of_index(
                    i, j, 1, 
                    4, 1 
                )
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-4) * dx))
                self.assertAlmostRelativeEquals(y, 100.0 + (0.5 * dy) + (j * dy))
                self.assertAlmostRelativeEquals(z, 0.0)
    
    
    def test28(self):
        results = []
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(3, 3, 3, 6,12,18)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
        dx = 6.0 / 3.0
        dy = 12.0 / 3.0
        dz = 18.0 / 3.0
        for i in range(4):
            for j in range(3):
                for k in range(3):
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        2, 1 
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, 6.0 + (0.5 * dx) + (i * dx))
                    self.assertAlmostRelativeEquals(y, (0.5 * dy) + (j * dy))
                    self.assertAlmostRelativeEquals(z, (0.5 * dz) + (k * dz))
        
        for i in range(3 + 8):
            for j in range(4):
                for k in range(3):
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        3, 1 
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-4) * dx))
                    self.assertAlmostRelativeEquals(y, ((0.5 * dy) - ((4-j) * dy)))
                    self.assertAlmostRelativeEquals(z, (0.5 * dz) + (k * dz))
                    
                    
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        4, 1 
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-4) * dx))
                    self.assertAlmostRelativeEquals(y, 12.0 + (0.5 * dy) + (j * dy))
                    self.assertAlmostRelativeEquals(z, (0.5 * dz) + (k * dz))
        
        for i in range(3 + 8):
            for j in range(3 + 8):
                for k in range(4):
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        5, 1 
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-4) * dx))
                    self.assertAlmostRelativeEquals(y, (0.5 * dy) + ((j-4) * dy))
                    self.assertAlmostRelativeEquals(z, ((0.5 * dz) - ((4-k) * dz)))
                    
                    
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        6, 1 
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-4) * dx))
                    self.assertAlmostRelativeEquals(y, (0.5 * dy) + ((j-4) * dy))
                    self.assertAlmostRelativeEquals(z, 18.0 + (0.5 * dz) + (k * dz))
    
    def test29(self):
        results = []
        instance=self.new_instance(AthenaInterface, number_of_workers = 3)
        instance.initialize_code()
        instance.setup_mesh(300,1,1,300.0,0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
           
        for j in [1,2]:
            print j
            for i in range(4):
                error = instance.set_boundary_state(
                    i,0,0,       #  index
                    1.0 * (i+1),         #  density
                    2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                    5.0 * (i+1),         #  energy
                    j, 1.0     #  boundary + grid
                )
                self.assertEquals(error, 0)
                rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                    i, 0, 0,
                    j, 1.0
                )
                self.assertEquals(error, 0)
                
                self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
                
    def test30(self):
        results = []
        instance=self.new_instance(AthenaInterface, number_of_workers = 3)
        instance.initialize_code()
        instance.setup_mesh(30,10,1,30.0,10.0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","periodic","periodic")
        instance.commit_parameters()
        
           
        for boundaryindex in [3,4]:
            for i0 in range(38):
                for j0 in range(4):
                    i = (i0 * (4*38)) + j0
                    error = instance.set_boundary_state(
                        i0,j0,0,       #  index
                        1.0 * (i+1),         #  density
                        2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                        5.0 * (i+1),         #  energy
                        boundaryindex, 1     #  boundary + grid
                    )
                    self.assertEquals(error, 0)
                    rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                        i0, j0, 0,
                        boundaryindex, 1.0
                    )
                    self.assertEquals(error, 0)
                    
                    self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
                    
    def test31(self):
        results = []
        instance=self.new_instance(AthenaInterface, number_of_workers = 3)
        instance.initialize_code()
        instance.set_auto_decomposition(0)
        instance.set_parallel_decomposition(1,3,1)
        instance.setup_mesh(5,6,1,5.0,6.0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","periodic","periodic")
        instance.commit_parameters()
        
           
        for boundaryindex in [1,2]:
            for i0 in range(4):
                for j0 in range(6):
                    i = (i0 * (4*6)) + j0
                    error = instance.set_boundary_state(
                        i0,j0,0,       #  index
                        1.0 * (i+1),         #  density
                        2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                        5.0 * (i+1),         #  energy
                        boundaryindex, 1     #  boundary + grid
                    )
                    self.assertEquals(error, 0)
                    rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                        i0, j0, 0,
                        boundaryindex, 1.0
                    )
                    self.assertEquals(error, 0)
                    
                    self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
        
    def test32(self):
        results = []
        instance=self.new_instance(AthenaInterface, number_of_workers = 3*3)
        instance.initialize_code()
        instance.set_auto_decomposition(0)
        instance.set_parallel_decomposition(3,3,1)
        instance.setup_mesh(6,6,1,6.0,6.0,0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","periodic","periodic")
        instance.commit_parameters()
        
           
        for boundaryindex in [1,2]:
            for i0 in range(4):
                for j0 in range(6):
                    i = (i0 * (4*6)) + j0
                    error = instance.set_boundary_state(
                        i0,j0,0,       #  index
                        1.0 * (i+1),         #  density
                        2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                        5.0 * (i+1),         #  energy
                        boundaryindex, 1     #  boundary + grid
                    )
                    self.assertEquals(error, 0)
                    rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                        i0, j0, 0,
                        boundaryindex, 1.0
                    )
                    self.assertEquals(error, 0)
                    
                    self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
        for boundaryindex in [3,4]:
            for i0 in range(6+8):
                for j0 in range(4):
                    i = (i0 * (4*(6+8))) + j0
                    error = instance.set_boundary_state(
                        i0,j0,0,       #  index
                        1.0 * (i+1),         #  density
                        2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                        5.0 * (i+1),         #  energy
                        boundaryindex, 1     #  boundary + grid
                    )
                    self.assertEquals(error, 0)
                    rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                        i0, j0, 0,
                        boundaryindex, 1.0
                    )
                    self.assertEquals(error, 0)
                    
                    self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    def test33(self):
        results = []
        instance=self.new_instance(AthenaInterface, number_of_workers = 3)
        instance.initialize_code()
        instance.set_auto_decomposition(0)
        instance.set_parallel_decomposition(1,1,3)
        instance.setup_mesh(5,5,6,5.0,5.0,30.0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","periodic","periodic")
        instance.commit_parameters()
        
           
        for boundaryindex in [1,2]:
            for i0 in range(4):
                for j0 in range(5):
                    for z0 in range(6):
                        i = (i0 * (5*6)) + (j0 * 6) + z0
                        error = instance.set_boundary_state(
                            i0,j0,z0,       #  index
                            1.0 * (i+1),         #  density
                            2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                            5.0 * (i+1),         #  energy
                            boundaryindex, 1     #  boundary + grid
                        )
                        self.assertEquals(error, 0)
                        rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                            i0, j0, z0,
                            boundaryindex, 1.0
                        )
                        self.assertEquals(error, 0)
                        
                        self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    
    def test34(self):
        results = []
        instance=self.new_instance(AthenaInterface, number_of_workers = 3)
        instance.initialize_code()
        instance.set_auto_decomposition(0)
        instance.set_parallel_decomposition(3,1,1)
        instance.setup_mesh(6,5,5,6.0,5.0,5.0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
           
        for boundaryindex in [5,6]:
            for i0 in range(6+8):
                for j0 in range(5+8):
                    for z0 in range(4):
                        i = (i0 * (5*4)) + (j0 * 4) + z0
                        error = instance.set_boundary_state(
                            i0,j0,z0,       #  index
                            1.0 * (i+1),         #  density
                            2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                            5.0 * (i+1),         #  energy
                            boundaryindex, 1     #  boundary + grid
                        )
                        self.assertEquals(error, 0)
                        rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                            i0, j0, z0,
                            boundaryindex, 1.0
                        )
                        self.assertEquals(error, 0)
                        
                        self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
             
    
    def test35(self):
        results = []
        instance=self.new_instance(AthenaInterface, number_of_workers = 9)
        instance.initialize_code()
        instance.set_auto_decomposition(0)
        instance.set_parallel_decomposition(3,3,1)
        instance.setup_mesh(6,6,5,6.0,6.0,5.0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
        boundary_indices = []
        all_i0 = []
        all_j0 = []
        all_z0 = []
        all_i = []
        
        for boundaryindex in [5,6]:
            for i0 in range(6+8):
                for j0 in range(6+8):
                    for z0 in range(4):
                        boundary_indices.append(boundaryindex)
                        all_i0.append(i0)
                        all_j0.append(j0)
                        all_z0.append(z0)
                        
                        
                        i = (i0 * (5*4)) + (j0 * 4) + z0
                        
                        all_i.append(i)
        all_i = numpy.asarray(all_i)
        error = instance.set_boundary_state(
            all_i0,all_j0,all_z0,       #  index
            1.0 * (all_i+1),         #  density
            2.0 * (all_i+1), 3.0 * (all_i+1), 4.0 * (all_i+1), #  momentum
            5.0 * (all_i+1),         #  energy
            boundary_indices, 1     #  boundary + grid
        )
        self.assertEquals(error, 0)
        rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
            all_i0, all_j0, all_z0,
            boundaryindex, 1.0
        )
        self.assertEquals(error, 0)
        
        self.assertAlmostRelativeEquals(rho, 1.0 * (all_i+1))
        self.assertAlmostRelativeEquals(rhovx, 2.0 * (all_i+1))
        self.assertAlmostRelativeEquals(rhovy, 3.0 * (all_i+1))
        self.assertAlmostRelativeEquals(rhovz, 4.0 * (all_i+1))
        self.assertAlmostRelativeEquals(rhoen, 5.0 * (all_i+1))
                        
class TestAthena(TestWithMPI):
    
        
    def test0(self):
        instance=self.new_instance(Athena)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.parameters.nx = 10
        instance.parameters.ny = 20
        instance.parameters.nz = 40
        instance.parameters.length_x = 1 | generic_unit_system.length
        instance.parameters.length_y = 2 | generic_unit_system.length
        instance.parameters.length_z = 3 | generic_unit_system.length
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        
        #print instance.grid[0].y
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        firstx = instance.grid[0][0][0].x
        allx = instance.grid[0].x
        for j in range(20):
            for k in range(40):
                self.assertEquals(allx[j][k], firstx)
        
        print instance.grid[0][0].rho
        self.assertEquals(instance.grid[0][0][0].rho, 0.0 |generic_unit_system.mass / generic_unit_system.length ** 3)
        
        
        potential_grid = datamodel.Grid(12,22,42)
        potential_grid.potential = 2.0 | generic_unit_system.potential
        channel = potential_grid.new_channel_to(instance.potential_grid)
        channel.copy()
        self.assertEquals(instance.potential_grid[0][0][0].potential, 2.0 | generic_unit_system.potential)
        self.assertEquals(instance.potential_grid[0][2][20].potential, 2.0 | generic_unit_system.potential)
        
        instance.stop()
    
    def test1(self):
        instance=self.new_instance(Athena)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.parameters.nx = 10
        instance.parameters.ny = 20
        instance.parameters.nz = 40
        instance.parameters.length_x = 1 | generic_unit_system.length
        instance.parameters.length_y = 2 | generic_unit_system.length
        instance.parameters.length_z = 3 | generic_unit_system.length
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        result = instance.commit_parameters()
                
        firstx = instance.potential_grid[0][0][0].x
        print firstx
        self.assertEquals(firstx, -0.05 | generic_unit_system.length)
        allx = instance.potential_grid[0].x
        for j in range(20):
            for k in range(40):
                self.assertEquals(allx[j][k], firstx)
        instance.stop()
    
    

    def test2(self):
        instance=self.new_instance(Athena)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.setup_mesh(10, 10, 1, 1.0 | generic_unit_system.length, 1.0 | generic_unit_system.length , 0.0 | generic_unit_system.length)
        instance.parameters.mesh_size = (10,10,1)
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 0.0 | generic_unit_system.length
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
    
        grid = datamodel.Grid(10,10,1)
        grid.rho = 0.1 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
            
        
        print instance.grid[1].rho
        self.assertEquals(instance.grid[1][1][0].rho, 0.1 | density)
        for x in instance.grid[1].rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
    
        instance.evolve_model(10.0 | generic_unit_system.time)
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
        instance.stop()
    
    

    def test3(self):
        instance=self.new_instance(Athena)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.parameters.nx = 10
        instance.parameters.ny = 10
        instance.parameters.nz = 1
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 0.0 | generic_unit_system.length
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        instance.set_has_external_gravitational_potential(1)
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
    
        grid = datamodel.Grid(10,10,1)
        grid.rho = 0.1 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
        
        self.assertEquals(grid._get_writeable_attribute_names(), ['energy', 'rho', 'rhovx', 'rhovy', 'rhovz', ] )
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        potential_grid = datamodel.Grid(12,12,1)
        potential_grid.potential = 0.0 | generic_unit_system.potential
        channel = potential_grid.new_channel_to(instance.potential_grid)
        channel.copy()
        result = instance.initialize_grid()
        
        self.assertEquals(instance.grid[1][1][0].rho, 0.1 | density)
        for x in instance.grid[1].rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
    
        instance.evolve_model(10.0 | generic_unit_system.time)
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
        instance.stop()
    
    

    def test4(self):
        instance=self.new_instance(Athena)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.parameters.nx = 10
        instance.parameters.ny = 10
        instance.parameters.nz = 1
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 0.0 | generic_unit_system.length
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        instance.set_has_external_gravitational_potential(1)
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
    
        grid = datamodel.Grid(10,10,1)
        grid.rho = 0.1 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        potential_grid = datamodel.Grid(12,12,1)
        potential_grid.potential = 0.0 | generic_unit_system.potential
        x = instance.potential_grid.x
        y = instance.potential_grid.y
        
        
        for i in range(12):
            for j in range(12):
                px = x[i][j][0].value_in(generic_unit_system.length)
                py = y[i][j][0].value_in(generic_unit_system.length)
                potential =  (math.sin(py * math.pi)+math.sin(px *math.pi)) / 200.0
                if px < 0 or px > 1.0:
                    potential = 0.0
                if py < 0 or py > 1.0:
                    potential = 0.0
                instance.potential_grid[i][j][0].potential = generic_unit_system.potential.new_quantity([potential])
        #channel = potential_grid.new_channel_to(instance.potential_grid)
        #channel.copy()
        result = instance.initialize_grid()
        
        self.assertEquals(instance.grid[1][1][0].rho, 0.1 | density)
        for x in instance.grid[1].rho.value_in(density).flatten():
            self.assertAlmostRelativeEquals(x, 0.1)
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        #print instance.grid.rhox
        z = instance.grid.rho[...,...,0]
        z = instance.potential_grid.potential[...,...,0]
        z = z.value_in(generic_unit_system.potential)
        #from matplotlib import pyplot
        #pyplot.imshow(z)
        #pyplot.savefig("bla.png")
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertNotEquals(x, 0.1)
        instance.stop()
    
    

    def test5(self):
        instance=self.new_instance(Athena)
        self.assertAlmostRelativeEquals(instance.parameters.isothermal_sound_speed, 0.0 | generic_unit_system.speed)
        instance.parameters.isothermal_sound_speed = 0.1 | generic_unit_system.speed
        self.assertAlmostRelativeEquals(instance.parameters.isothermal_sound_speed, 0.1 | generic_unit_system.speed)
        self.assertAlmostRelativeEquals(instance.parameters.gamma, 1.66666666666666667)
        instance.parameters.gamma = 0.1
        self.assertAlmostRelativeEquals(instance.parameters.gamma, 0.1)
        self.assertAlmostRelativeEquals(instance.parameters.courant_number, 0.3)
        instance.parameters.courant_number = 0.1
        self.assertAlmostRelativeEquals(instance.parameters.courant_number, 0.1)
        
        print instance.parameters
        instance.stop()
    
    

    def test6(self):
        instance=self.new_instance(Athena)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.setup_mesh(10 , 20, 40, 1.0 | generic_unit_system.length, 1.0 | generic_unit_system.length, 1.0 | generic_unit_system.length)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        grid = datamodel.Grid.create((10,10,10), [10.0, 10.0, 10.0] | units.m)
        
        grid.rho = 0.4 | density
        grid.rhovx = 0.1 | momentum
        grid.rhovy = 0.2 | momentum
        grid.rhovz = 0.3 | momentum
        grid.energy = 0.5 | energy
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        self.assertEquals(instance.grid[0][0][0].rho, 0.4 | density)
        self.assertEquals(instance.grid.rho.number.ndim, 3)
        
        self.assertEquals(len(list(instance.itergrids())), 1)
        
        instance.stop()

    def test7(self):
        instance=self.new_instance(Athena)
        instance.parameters.isothermal_sound_speed = 0.1 | generic_unit_system.speed
        instance.parameters.gamma = 0.1
        instance.parameters.courant_number = 0.1
       
        instance.parameters.nx = 10
        instance.parameters.ny = 20
        instance.parameters.nz = 40
        
        instance.parameters.length_x = 10 | generic_unit_system.length
        instance.parameters.length_y = 20 | generic_unit_system.length
        instance.parameters.length_z = 30 | generic_unit_system.length
        
        print instance.parameters
        instance.commit_parameters()
        
        mini,maxi, minj,maxj, mink,maxk = instance.get_index_range_inclusive()
        
        self.assertEquals(mini, 0)
        self.assertEquals(maxi, 9)
        self.assertEquals(minj, 0)
        self.assertEquals(maxj, 19)
        self.assertEquals(mink, 0)
        self.assertEquals(maxk, 39)
        self.assertEquals(instance.parameters.mesh_size, (10,20,40))
        print instance.parameters
        instance.stop()
    
    def test8(self):
        instance=self.new_instance(Athena)
        instance.parameters.stopping_conditions_number_of_steps = 10
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps, 10)
        instance.stop()

    def test8a(self):
        instance=self.new_instance(Athena)
        instance.parameters.stopping_conditions_timeout = 10 | units.s
        self.assertEquals(instance.parameters.stopping_conditions_timeout, 10|units.s)
        instance.stop()

    def test9(self):
        instance=self.new_instance(Athena)
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        self.assertEquals(instance.parameters.xbound1, "periodic")
        instance.stop()

    def xtest10(self):
        instance=self.new_instance(Athena)
        instance.parameters.gamma = 5/3.0
        instance.parameters.courant_number=0.3
        
        n = 100
        
        instance.parameters.nx = n
        instance.parameters.ny = n
        instance.parameters.nz = n
        
        instance.parameters.length_x = 1 | generic_unit_system.length
        instance.parameters.length_y = 1 | generic_unit_system.length
        instance.parameters.length_z = 1 | generic_unit_system.length
        
        instance.x_boundary_conditions = ("periodic","periodic")
        instance.y_boundary_conditions = ("periodic","periodic")
        instance.z_boundary_conditions = ("periodic","periodic")
        
        result = instance.commit_parameters()
         
        density = generic_unit_system.mass / (generic_unit_system.length**3)
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        grid = datamodel.Grid(n,n,n)
        grid.rho =  0.0 | generic_unit_system.density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
        
        halfway = n/2 - 1
        grid[:halfway].rho = 4.0  | generic_unit_system.density
        grid[:halfway].energy = (1.0 | energy)/ (instance.parameters.gamma - 1)
        grid[halfway:].rho = 1.0  | generic_unit_system.density
        grid[halfway:].energy = (0.1795 | energy)/ (instance.parameters.gamma - 1)
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
    
        #from amuse import plot
        #from matplotlib import pyplot
        #print grid.rho[...,0,0]
        #plot.plot(instance.grid.x[...,0,0], grid.rho[...,0,0])
        #pyplot.savefig("bla1.png")
        
        error = instance.initialize_grid()
        
        instance.evolve_model(0.12 | generic_unit_system.time)
        
        
        channel = instance.grid.new_channel_to(grid)
        channel.copy()
        
        #print grid.rho[...,0,0]
        #plot.plot(instance.grid.x[...,0,0], grid.rho[...,0,0])
        #pyplot.savefig("bla2.png")
        instance.stop()
        
    
    def test11(self):
        instance=self.new_instance(Athena, mode=AthenaInterface.MODE_SELF_GRAVITY) #, redirection = "none") #, debugger="gdb")
        instance.parameters.gamma = 5/3.0
        instance.parameters.courant_number=0.3
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        instance.parameters.four_pi_G =  4 * numpy.pi * (1|(generic_unit_system.length**3) / (generic_unit_system.mass * (generic_unit_system.time**2))) # G = 1, like nbody
        instance.parameters.gravity_mean_rho = 0.0 | density
        

        datamodel.Grid.add_global_vector_attribute("position", ["x","y","z"])
        n = 10
        
        instance.parameters.nx = n
        instance.parameters.ny = n
        instance.parameters.nz = n
        
        instance.parameters.length_x = 4.0 | generic_unit_system.length
        instance.parameters.length_y = 4.0 | generic_unit_system.length
        instance.parameters.length_z = 4.0 | generic_unit_system.length
        
        instance.x_boundary_conditions = ("periodic","periodic")
        instance.y_boundary_conditions = ("periodic","periodic")
        instance.z_boundary_conditions = ("periodic","periodic")
        
        result = instance.commit_parameters()
    
        grid = datamodel.Grid.create((n,n,n), [4.0 , 4.0, 4.0] | generic_unit_system.length)
        
        grid.rho = 0.0 | generic_unit_system.density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.001 | energy
        
        scaled_radius = 1.0 / 1.695 | generic_unit_system.length
        total_mass = 1.0 | generic_unit_system.mass
        
        radii = (grid.position - ([2.0, 2.0, 2.0] | generic_unit_system.length)).lengths()
        
        rho_sphere = ((0.75 * total_mass /  (numpy.pi * (scaled_radius ** 3))))
        grid.rho =  (rho_sphere * ((1 + (radii ** 2) / (scaled_radius ** 2))**(-5.0/2.0)))
        
        internal_energy = (0.25 | generic_unit_system.time ** -2  * generic_unit_system.mass ** -1 * generic_unit_system.length **3) * total_mass / scaled_radius 
        grid.energy = grid.rho * internal_energy/(1+(radii/scaled_radius)**2)**(1.0/2.0)
        
      
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        
        instance.initialize_grid()
        
        instance.evolve_model(0.01 | generic_unit_system.time)
        G = 1.0 | generic_unit_system.length **3 * generic_unit_system.mass**-1 * generic_unit_system.time**-2
        a = instance.grid[5][5].gravitational_potential
        b = (-1 * G * total_mass / (radii**2+scaled_radius**2).sqrt()) [5][5]
        
        for x in a:
            self.assertTrue(x < 0 | generic_unit_system.potential)
            
        
        a = instance.grid[5][5].gravitational_acceleration_z
        for index, x in enumerate(a):
            if index < 5:
                self.assertTrue(x > 0 | generic_unit_system.acceleration)
            else:
                self.assertTrue(x < 0 | generic_unit_system.acceleration)
        instance.stop()
    
    def test12(self):
        print "Testing Athena grid setters"
        instance=self.new_instance(Athena)
        instance.parameters.isothermal_sound_speed = 0.1 | generic_unit_system.speed
        instance.parameters.gamma = 5/3.0
        instance.parameters.courant_number = 0.3
        instance.parameters.mesh_size = (2, 2, 2)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("periodic", "periodic")
        
        instance.grid.rho =    1.0 | generic_unit_system.density
        self.assertAlmostEquals(instance.grid.rho, 
            numpy.ones((2,2,2)) | generic_unit_system.density)
        
        instance.grid.momentum = numpy.reshape(
            numpy.arange(0.0, 3.0, 0.125), (2,2,2,3))  | generic_unit_system.momentum_density
        self.assertAlmostEquals(instance.grid.rhovx, 
            numpy.reshape(numpy.arange(0.000, 3.0, 0.375), (2,2,2)) | generic_unit_system.momentum_density)
        self.assertAlmostEquals(instance.grid.rhovy, 
            numpy.reshape(numpy.arange(0.125, 3.0, 0.375), (2,2,2)) | generic_unit_system.momentum_density)
        self.assertAlmostEquals(instance.grid.rhovz, 
            numpy.reshape(numpy.arange(0.250, 3.0, 0.375), (2,2,2)) | generic_unit_system.momentum_density)
        
        momentum = instance.grid.momentum
        rhovx = -momentum.x
        rhovy = 2 * momentum.z
        rhovz = -0.5 * momentum.y
        instance.grid.momentum = VectorQuantity.new_from_scalar_quantities(rhovx,rhovy,rhovz).transpose(axes=(1,2,3,0))
        self.assertAlmostEquals(instance.grid.rhovx, 
            numpy.reshape(numpy.arange(0.000, -3.0, -0.375), (2,2,2)) | generic_unit_system.momentum_density)
        self.assertAlmostEquals(instance.grid.rhovy, 
            numpy.reshape(numpy.arange(0.5, 6.0, 0.75), (2,2,2)) | generic_unit_system.momentum_density)
        self.assertAlmostEquals(instance.grid.rhovz, 
            numpy.reshape(numpy.arange(-0.0625, -1.5, -0.1875), (2,2,2)) | generic_unit_system.momentum_density)
        
        instance.grid[...,0,...].momentum =  [12.0, 13.0, 14.0] | generic_unit_system.momentum_density
        self.assertAlmostEquals(instance.grid[0,...,...].rhovx, 
            [[12.0, 12.0], [-0.75, -1.125]] | generic_unit_system.momentum_density)
        self.assertAlmostEquals(instance.grid[0,...,...].rhovy, 
            [[13.0, 13.0], [2.0, 2.75]] | generic_unit_system.momentum_density)
        self.assertAlmostEquals(instance.grid[...,...,0].rhovz, 
            [[14.0, -0.4375], [14.0, -1.1875]] | generic_unit_system.momentum_density)
        
        instance.grid.energy = numpy.reshape(numpy.arange(0.0, 1.0, 0.125), (2,2,2)) | generic_unit_system.energy_density
        self.assertAlmostEquals(instance.grid[...,0,0].energy, 
            [0.0, 0.5] | generic_unit_system.energy_density)
        self.assertAlmostEquals(instance.grid[0,...,0].energy, 
            [0.0, 0.25] | generic_unit_system.energy_density)
        self.assertAlmostEquals(instance.grid[0,0,...].energy, 
            [0.0, 0.125] | generic_unit_system.energy_density)
        instance.initialize_grid()
        instance.stop()
    
    
    def test13(self): 
        converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
            1 | units.parsec,
            1 | units.Myr,
            1 | units.MSun
        )
        instance=self.new_instance(Athena, unit_converter = converter)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        
        instance.parameters.mesh_size = (10 , 20, 40)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | units.parsec
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("periodic", "periodic")
        
        
        density = units.MSun / (units.parsec ** 3)
        momentum =  units.MSun / (units.Myr * units.parsec ** 2 )
        energy =  units.MSun / (units.parsec * units.Myr ** 2)
        
        grid = datamodel.Grid.create((10,20,40), [1.0, 1.0, 1.0] | units.parsec )
        
        grid.rho = 0.4 | density
        grid.rhovx = 0.1 | momentum
        grid.rhovy = 0.2 | momentum
        grid.rhovz = 0.3 | momentum
        grid.energy = 0.5 | energy
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        print instance.grid[0].rho
        self.assertAlmostRelativeEquals(instance.grid[0].rho, 0.4 | density)
        self.assertAlmostRelativeEquals(instance.grid[0].rhovx,  0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid[0].rhovy,  0.2 | momentum)
        self.assertAlmostRelativeEquals(instance.grid[0].rhovz,  0.3 | momentum)
        self.assertAlmostRelativeEquals(instance.grid[0].energy,  0.5 | energy)
        self.assertEquals(instance.grid.rho.number.ndim, 3)
        
        self.assertEquals(len(list(instance.itergrids())), 1)
        
        instance.stop()

        
    def test14(self): 
        instance=self.new_instance(Athena)
        instance.parameters.mesh_size = (10 , 1, 1)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("interface", "outflow")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        grid = datamodel.Grid.create((10,1,1), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.1 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (instance.parameters.gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        #instance.grid.boundaries.left.
        xbound1 = instance.get_boundary_grid('xbound1')
        self.assertEquals(xbound1.shape, (4,1,1))
        memxbound1 = xbound1.copy()
        memxbound1.rho = 0.02 | density
        memxbound1.rhovx = 0.2 | momentum
        memxbound1.rhovy = 0.0 | momentum
        memxbound1.rhovz = 0.0 | momentum
        memxbound1.energy =  p / (instance.parameters.gamma - 1)
        memxbound1.energy += 0.5 * (memxbound1.rhovx ** 2  + memxbound1.rhovy ** 2 + memxbound1.rhovz ** 2) / memxbound1.rho
        channel = memxbound1.new_channel_to(xbound1)
        channel.copy()
        
        instance.evolve_model(1.0 | generic_unit_system.time)
        print instance.stopping_conditions.number_of_steps_detection.is_set()
        rho = instance.grid.rho[...,0,0]
        self.assertAlmostRelativeEquals(rho[-1], 0.01 | density)
        self.assertTrue(rho[0] > 0.01 | density)
        self.assertTrue(instance.grid.rhovx[0,0,0] > 0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[-1,0,0] , 0.1 | momentum)
        print instance.model_time
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[...,0,0]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[...,0,0], 0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()
        
    def test15(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.stop()
        
    def test16(self): 
        instance=self.new_instance(Athena)
        instance.parameters.mesh_size = (10 , 1, 1)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("outflow", "interface")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        grid = datamodel.Grid.create((10,1,1), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = -0.1 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (instance.parameters.gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        #instance.grid.boundaries.left.
        xbound = instance.get_boundary_grid('xbound2')
        self.assertEquals(xbound.shape, (4,1,1))
        memxbound = xbound.copy()
        memxbound.rho = 0.02 | density
        memxbound.rhovx = -0.2 | momentum
        memxbound.rhovy = 0.0 | momentum
        memxbound.rhovz = 0.0 | momentum
        memxbound.energy =  p / (instance.parameters.gamma - 1)
        memxbound.energy += 0.5 * (memxbound.rhovx ** 2  + memxbound.rhovy ** 2 + memxbound.rhovz ** 2) / memxbound.rho
        channel = memxbound.new_channel_to(xbound)
        channel.copy()
        
        instance.evolve_model(1.0 | generic_unit_system.time)
        print instance.stopping_conditions.number_of_steps_detection.is_set()
        rho = instance.grid.rho[...,0,0]
        self.assertAlmostRelativeEquals(rho[0], 0.01 | density)
        self.assertTrue(rho[-1] > 0.01 | density)
        self.assertTrue(instance.grid.rhovx[-1,0,0] < -0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[0,0,0] , -0.1 | momentum)
        print instance.model_time
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[...,0,0]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[...,0,0], -0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()
    
    
    def test17(self): 
        instance=self.new_instance(Athena, number_of_workers = 2)
        instance.set_parallel_decomposition(1,2,1)
        instance.parameters.mesh_size = (10,4,1)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("interface", "outflow")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        grid = datamodel.Grid.create((10,4,1), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.1 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (instance.parameters.gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        xbound = instance.get_boundary_grid('xbound1')
        self.assertEquals(xbound.shape, (4,4,1))
        memxbound = xbound.copy()
        memxbound.rho = 0.02 | density
        memxbound.rhovx = 0.2 | momentum
        memxbound.rhovy = 0.0 | momentum
        memxbound.rhovz = 0.0 | momentum
        memxbound.energy =  p / (instance.parameters.gamma - 1)
        memxbound.energy += 0.5 * (memxbound.rhovx ** 2  + memxbound.rhovy ** 2 + memxbound.rhovz ** 2) / memxbound.rho
        channel = memxbound.new_channel_to(xbound)
        channel.copy()
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        print instance.stopping_conditions.number_of_steps_detection.is_set()
        print instance.grid.rho
        rho = instance.grid.rho[...,0,0]
        self.assertAlmostRelativeEquals(rho[-1], 0.01 | density)
        self.assertTrue(rho[0] > 0.01 | density)
        self.assertTrue(instance.grid.rhovx[0,0,0] > 0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[-1,0,0] , 0.1 | momentum)
        print instance.model_time
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[...,0,0]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[...,0,0], 0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()

    
    def test18(self): 
        instance=self.new_instance(Athena, number_of_workers = 2)
        instance.set_parallel_decomposition(2,1,1)
        instance.parameters.mesh_size = (4,10,1)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("interface", "outflow")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        grid = datamodel.Grid.create((4,10,1), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.1 | momentum
        grid.rhovz = 0.0 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (instance.parameters.gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        ybound = instance.get_boundary_grid('ybound1')
        self.assertEquals(ybound.shape, (4+8,4,1))
        memybound = ybound.copy()
        memybound.rho = 0.02 | density
        memybound.rhovx = 0.0 | momentum
        memybound.rhovy = 0.2 | momentum
        memybound.rhovz = 0.0 | momentum
        memybound.energy =  p / (instance.parameters.gamma - 1)
        memybound.energy += 0.5 * (memybound.rhovx ** 2  + memybound.rhovy ** 2 + memybound.rhovz ** 2) / memybound.rho
        
        channel = memybound.new_channel_to(ybound)
        channel.copy()
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        print instance.stopping_conditions.number_of_steps_detection.is_set()
        print instance.grid.rho
        rho = instance.grid.rho[0,...,0]
        self.assertAlmostRelativeEquals(rho[-1], 0.01 | density)
        self.assertTrue(rho[0] > 0.01 | density)
        self.assertTrue(instance.grid.rhovy[0,0,0] > 0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovy[0,-1,0] , 0.1 | momentum)
        print instance.model_time
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[0,...,0]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovy[0,...,0], 0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()
    
    
    def test19(self): 
        instance=self.new_instance(Athena, number_of_workers = 1)
        instance.set_parallel_decomposition(1,1,1)
        instance.parameters.mesh_size = (4,5,6)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("interface", "outflow")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        grid = datamodel.Grid.create((4,5,6), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.1 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (instance.parameters.gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        zbound = instance.get_boundary_grid('zbound1')
        self.assertEquals(zbound.shape, (4+8,5+8,4))
        memzbound = zbound.copy()
        memzbound.rho = 0.02 | density
        memzbound.rhovx = 0.0 | momentum
        memzbound.rhovy = 0.0 | momentum
        memzbound.rhovz = 0.2 | momentum
        memzbound.energy =  p / (instance.parameters.gamma - 1)
        memzbound.energy += 0.5 * (memzbound.rhovx ** 2  + memzbound.rhovy ** 2 + memzbound.rhovz ** 2) / memzbound.rho
       
        channel = memzbound.new_channel_to(zbound)
        channel.copy()
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[0,0,...]
        self.assertAlmostRelativeEquals(rho[-1], 0.01 | density)
        self.assertTrue(rho[0] > 0.01 | density)
        self.assertTrue(instance.grid.rhovz[0,0,0] > 0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovz[0,0,-1] , 0.1 | momentum)
        print instance.model_time
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[0,...,0]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovz[0,0,...], 0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()
    
    
    def test20(self): 
        instance=self.new_instance(Athena, number_of_workers = 4)
        instance.parameters.parallel_decomposition = (2,2,1)
        instance.parameters.mesh_size = (4,5,6)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("outflow", "interface")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        grid = datamodel.Grid.create((4,5,6), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = -0.1 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (instance.parameters.gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        zbound = instance.get_boundary_grid('zbound2')
        self.assertEquals(zbound.shape, (4+8,5+8,4))
        memzbound = zbound.copy()
        memzbound.rho = 0.02 | density
        memzbound.rhovx = 0.0 | momentum
        memzbound.rhovy = 0.0 | momentum
        memzbound.rhovz = -0.2 | momentum
        memzbound.energy =  p / (instance.parameters.gamma - 1)
        memzbound.energy += 0.5 * (memzbound.rhovx ** 2  + memzbound.rhovy ** 2 + memzbound.rhovz ** 2) / memzbound.rho
       
        channel = memzbound.new_channel_to(zbound)
        channel.copy()
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[0,0,...]
        self.assertAlmostRelativeEquals(rho[0], 0.01 | density)
        self.assertTrue(rho[-1] > 0.01 | density)
        self.assertTrue(instance.grid.rhovz[0,0,-1] < -0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovz[0,0,0] , -0.1 | momentum)
        print instance.model_time
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[0,...,0]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovz[0,0,...], -0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()
        
    
    def test21(self):
        
        instance=self.new_instance(Athena)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (20.0, 1, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 1, 1)
        
        for x in instance.itergrids():
            inmem = x.copy()
            inmem.rho = inmem.x/(1| generic_unit_system.length) | generic_unit_system.density
            inmem.rhovx = 0.0 | generic_unit_system.momentum_density
            inmem.energy =  1.0 | generic_unit_system.energy_density
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()
            print inmem.rho
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(0.5| generic_unit_system.length,0.0| generic_unit_system.length,0.0| generic_unit_system.length)
        
        self.assertEquals(rho , 0.5 | generic_unit_system.density)
        
        for value in numpy.arange(0.5, 19.6, 0.1):
            
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value | generic_unit_system.length,
                0.0 | generic_unit_system.length,
                0.0 | generic_unit_system.length
            )
        
            self.assertAlmostRelativeEquals(rho , value | generic_unit_system.density)
        
        for value in numpy.arange(0.0, 0.6, 0.1):
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value | generic_unit_system.length,
                0.0 | generic_unit_system.length,
                0.0 | generic_unit_system.length
            )
            self.assertAlmostRelativeEquals(rho , ((0.5 + value) * 0.5 + (0.5-value) * 19.5) | generic_unit_system.density)
        
        
        for value in numpy.arange(0.0, 0.5, 0.1):
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value + 19.5| generic_unit_system.length,
                0.0 | generic_unit_system.length,
                0.0 | generic_unit_system.length
            )
            self.assertAlmostRelativeEquals(rho , (19.5 - (value * 19))  | generic_unit_system.density, 9)
        
        # out of range
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
            20.0| generic_unit_system.length,
            0.0 | generic_unit_system.length,
            0.0 | generic_unit_system.length
        )
        self.assertAlmostRelativeEquals(rho , 0.0 | generic_unit_system.density, 9)

    
    def test22(self):
        
        instance=self.new_instance(Athena, number_of_workers=2)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (20.0, 20.0, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 20, 1)
        
        for x in instance.itergrids():
            inmem = x.copy()
            inmem.rho = (inmem.x + ((inmem.y - (0.5| generic_unit_system.length))* 20.0))/(1| generic_unit_system.length) | generic_unit_system.density
            inmem.rhovx = 0.0 | generic_unit_system.momentum_density
            inmem.energy =  1.0 | generic_unit_system.energy_density
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()
            print inmem.rho[0], inmem.y[0], inmem.x[0]
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(0.5| generic_unit_system.length,0.5| generic_unit_system.length,0.0| generic_unit_system.length)
        
        self.assertEquals(rho , 0.5 | generic_unit_system.density)
        
        for value in numpy.arange(0.5, 19.6, 0.1):
            
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value | generic_unit_system.length,
                0.5 | generic_unit_system.length,
                0.0 | generic_unit_system.length
            )
        
            self.assertAlmostRelativeEquals(rho , value | generic_unit_system.density)
        
        for x in numpy.arange(8.5, 11.5, 0.25):
            for y in numpy.arange(0.5, 19.6, 0.25):
                rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                    x | generic_unit_system.length,
                    y | generic_unit_system.length,
                    0.0 | generic_unit_system.length
                )
            
                self.assertAlmostRelativeEquals(rho , x + (20 * (y-0.5))  | generic_unit_system.density)
            
    def test23(self):
        
        instance=self.new_instance(Athena, number_of_workers=3)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (20.0, 20.0, 20.0) | generic_unit_system.length
        instance.parameters.mesh_length = (20.0, 20.0, 20.0) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 20, 20)
        
        for x in instance.itergrids():
            inmem = x.copy()
            inmem.rho = (
                (
                    inmem.x + 
                    ((inmem.y - (0.5| generic_unit_system.length))* 20.0) +
                    ((inmem.z - (0.5| generic_unit_system.length))* 400.0)
                )
                /(1| generic_unit_system.length) | generic_unit_system.density
            )
            inmem.rhovx = 0.0 | generic_unit_system.momentum_density
            inmem.energy =  1.0 | generic_unit_system.energy_density
            from_model_to_code = inmem.new_channel_to(x)
            from_model_to_code.copy()
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(0.5| generic_unit_system.length,0.5| generic_unit_system.length,0.5| generic_unit_system.length)
        
        self.assertEquals(rho , 0.5 | generic_unit_system.density)
        
        for value in numpy.arange(0.5, 19.6, 0.1):
            
            rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
                value | generic_unit_system.length,
                0.5 | generic_unit_system.length,
                0.5 | generic_unit_system.length
            )
        
            self.assertAlmostRelativeEquals(rho , value | generic_unit_system.density)
        
        sample = sample = datamodel.Grid.create(
            (4, 4, 76),
            (2, 2, 19) | generic_unit_system.length
        )
        sample.x += 9.5 | generic_unit_system.length
        sample.y += 9.5 | generic_unit_system.length
        sample.z += 0.5 | generic_unit_system.length
        x = sample.x.flatten()
        y = sample.y.flatten()
        z = sample.z.flatten()
        
        rho, rhovx, rhovy, rhovx, rhoenergy = instance.get_hydro_state_at_point(
            x,
            y,
            z
        )
        half = 0.5 | generic_unit_system.length
        
        self.assertAlmostRelativeEquals(rho , (x + (20 * (y-half)) + (400 * (z-half)))/(1| generic_unit_system.length) | generic_unit_system.density )
            

    def test24(self):
        
        instance=self.new_instance(Athena, number_of_workers=1)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (10.0, 1, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 1, 1)
        instance.set_has_external_gravitational_potential(1)
        instance.commit_parameters()
        potential_grid = instance.potential_grid
        factor =  (2 | generic_unit_system.length / generic_unit_system.time**2)
        potential_grid.potential = potential_grid.x * factor
        
        x = numpy.arange(0,10.25, 0.1) | generic_unit_system.length
        y = 0.5 |generic_unit_system.length
        z = 0.5 |generic_unit_system.length
        interpolated = instance.get_interpolated_gravitational_potential(x,y,z)
        self.assertAlmostRelativeEquals(interpolated, x * factor)
        
    
    def test25(self):
        
        instance=self.new_instance(Athena, number_of_workers=1)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (5.0, 10.0, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 20, 1)
        instance.set_has_external_gravitational_potential(1)
        instance.commit_parameters()
        potential_grid = instance.potential_grid
        factor =  (2 | generic_unit_system.length / generic_unit_system.time**2)
        potential_grid.potential = potential_grid.y * factor
        print  potential_grid.y * factor
        y = numpy.arange(0,10.25, 0.1) | generic_unit_system.length
        x = (y * 0) + (2 |generic_unit_system.length)
        z = 0.5 |generic_unit_system.length
        interpolated = instance.get_interpolated_gravitational_potential(x,y,z)
        print y*factor
        self.assertAlmostRelativeEquals(interpolated, y * factor)
