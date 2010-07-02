import os
import sys
import numpy
import math

from amuse.test.amusetest import TestWithMPI
from amuse.legacy.athena.interface import AthenaInterface, Athena

from amuse.support.units import generic_unit_system
from amuse.support.data import core

from mpi4py import MPI


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
        instance=self.new_instance(AthenaInterface, debugger="none")
        instance.initialize_code()
        instance.par_setd("test", "test2", "%.15e", 1.123, "a test parameter")
        x = instance.par_getd("test", "test2")
        
        self.assertEquals(x, 1.123)
        instance.stop()
        
        
    def test3(self):
        instance=self.new_instance(AthenaInterface, debugger="none")
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.setup_mesh(128, 1, 1, 1.0, 0.0, 0.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        
        x = instance.par_geti("grid", "Nx1")
        self.assertEquals(x, 128)
        
        result = instance.commit_parameters()
        self.assertEquals(result, 0)
        
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
        
        
        x,y,z, error= instance.get_position_of_index(15,5,5)
        self.assertEquals(error, -1)
        instance.stop()
        
    
    def test5(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(5, 5, 5, 1.0, 1.0, 1.0)
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        result = instance.commit_parameters()
        self.assertEquals(result, 0)
        
        error = instance.fill_grid_state_mpi(1,1,1, 0.1, 0.2, 0.3, 0.4, 0.5)
        self.assertEquals(error, 0)
        
        rho, rhovx, rhovy, rhovz, energy, error = instance.get_grid_state_mpi(1,1,1)
        
        self.assertEquals(error, 0)
        self.assertEquals(rho, 0.1)
        self.assertEquals(rhovx, 0.2)
        self.assertEquals(rhovy, 0.3)
        self.assertEquals(rhovz, 0.4)
        self.assertEquals(energy, 0.5)
        
        error = instance.initialize_grid()
        self.assertEquals(error, 0)
        
        timestep, error =  instance.get_timestep()
        self.assertEquals(error, 0)
        
        
        instance.stop()
        
    def test6(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.setup_mesh(5, 5, 5, 1.0, 1.0, 1.0)
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
        
        timestep, error =  instance.get_timestep()
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(timestep, 0.159999, 5)
        
        instance.stop()
        
    def test7(self):
        results = []
        for x in range(1,4):
            instance=self.new_instance(AthenaInterface, number_of_workers=x, debugger="none")
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
            
            result = instance.get_grid_state_mpi(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128))
            results.append(list(result))
            
            
            instance.stop()
        
        for x in range(128):
            for y in range(6):
                self.assertEquals(results[1][y][x], results[0][y][x])
                self.assertEquals(results[2][y][x], results[0][y][x])
        
    
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
        
        timestep, error =  instance.get_timestep()
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEqual(timestep,  0.006249991, 5)
        
        rho0, rhovx0, rhovy0, rhovz0, energy0, error0  = instance.get_grid_state_mpi(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128))
        instance.evolve(5.0)
        rho, rhovx, rhovy, rhovz, energy, error  = instance.get_grid_state_mpi(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128))
        
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
            instance=self.new_instance(AthenaInterface, number_of_workers=x, debugger="none")
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
            instance.evolve(5.0)
        
            result = instance.get_grid_state_mpi(numpy.arange(0,128), numpy.zeros(128), numpy.zeros(128))
            results.append(list(result))
            
            
            instance.stop()
        
        for x in range(128):
            for y in range(6):
                self.assertEquals(results[1][y][x], results[0][y][x])
                self.assertEquals(results[2][y][x], results[0][y][x])
                self.assertEquals(results[3][y][x], results[0][y][x])
    
    def test10(self):
        instance=self.new_instance(AthenaInterface, debugger="none")
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
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
        
        i,j,k,error = instance.get_index_of_position(0.4, 0.1, 0.1)
        self.assertAlmostRelativeEqual(i, 1.5)
        self.assertAlmostRelativeEqual(j, 0)
        self.assertAlmostRelativeEqual(k, 0)
        
        x,y,z,error = instance.get_position_of_index(-1,0,0)
        print error
        self.assertAlmostRelativeEqual(x, -0.1)
        instance.stop()
        
    def test11(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
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

    def test12(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
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
        
    
    def test13(self):
        instance=self.new_instance(AthenaInterface)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
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


class TestAthena(TestWithMPI):
    def test0(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.setup_mesh(10, 20, 40, 1.0, 1.0, 1.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        result = instance.commit_parameters()
        
        #print instance.grid[0].y
        
        firstx = instance.grid[0][0][0].x
        allx = instance.grid[0].x
        for j in range(20):
            for k in range(40):
                self.assertEquals(allx[j][k], firstx)
        
        print instance.grid[0][0].rho
        self.assertEquals(instance.grid[0][0][0].rho , 0.0 |generic_unit_system.mass / generic_unit_system.length ** 3)
        
        instance.stop()
        
    
    def test1(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.setup_mesh(10, 20, 40, 1.0, 1.0, 1.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
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
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.setup_mesh(10, 10, 1, 1.0, 1.0, 0.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        result = instance.commit_parameters()
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass * generic_unit_system.length / generic_unit_system.time

        grid = core.Grid(10,10,1)
        grid.rho = 0.1 | density
        grid.rhox = 0.0 | momentum
        grid.rhoy = 0.0 | momentum
        grid.rhoz = 0.0 | momentum
        grid.energy = 0.0 | generic_unit_system.energy
        
        self.assertEquals(grid._get_writeable_attribute_names(), set(['rhoz', 'rhoy', 'rhox', 'energy', 'rho']) )
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
            
        result = instance.initialize_grid()   
        
        print instance.grid[1].rho
        self.assertEquals(instance.grid[1][1][0].rho, 0.1 | density)
        for x in instance.grid[1].rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
            
        instance.evolve(1.0 | generic_unit_system.time)
        
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)

        instance.evolve(10.0 | generic_unit_system.time)
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
        instance.stop()
            
   
    def test3(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.setup_mesh(10, 10, 1, 1.0, 1.0, 0.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        instance.set_has_external_gravitational_potential(1)
        
        result = instance.commit_parameters()
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass * generic_unit_system.length / generic_unit_system.time

        grid = core.Grid(10,10,1)
        grid.rho = 0.1 | density
        grid.rhox = 0.0 | momentum
        grid.rhoy = 0.0 | momentum
        grid.rhoz = 0.0 | momentum
        grid.energy = 0.0 | generic_unit_system.energy
        
        self.assertEquals(grid._get_writeable_attribute_names(), set(['rhoz', 'rhoy', 'rhox', 'energy', 'rho']) )
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        potential_grid = core.Grid(12,12,3)
        potential_grid.potential = 0.0 | generic_unit_system.energy
        channel = potential_grid.new_channel_to(instance.potential_grid)
        channel.copy()
        result = instance.initialize_grid() 
        
        self.assertEquals(instance.grid[1][1][0].rho, 0.1 | density)
        for x in instance.grid[1].rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
            
        instance.evolve(1.0 | generic_unit_system.time)
        
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)

        instance.evolve(10.0 | generic_unit_system.time)
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
        instance.stop()
    
    def test4(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.8)
        instance.setup_mesh(10, 10, 1, 1.0, 1.0, 0.0)
        instance.set_boundary("periodic","periodic","periodic","periodic","periodic","periodic")
        instance.set_has_external_gravitational_potential(1)
        
        result = instance.commit_parameters()
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass * generic_unit_system.length / generic_unit_system.time

        grid = core.Grid(10,10,1)
        grid.rho = 0.1 | density
        grid.rhox = 0.0 | momentum
        grid.rhoy = 0.0 | momentum
        grid.rhoz = 0.0 | momentum
        grid.energy = 0.0 | generic_unit_system.energy
        
        self.assertEquals(grid._get_writeable_attribute_names(), set(['rhoz', 'rhoy', 'rhox', 'energy', 'rho']) )
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        potential_grid = core.Grid(12,12,1)
        potential_grid.potential = 0.0 | generic_unit_system.energy
        print instance.potential_grid.shape
        print instance.potential_grid.x[0].number.shape
        print instance.grid.shape
        print instance.grid.x[0].number.shape
        x = instance.potential_grid.x[0]
        y = instance.potential_grid.y[0]
        print x
        print y
        for i in range(12):
            for j in range(12):
                px = x[i][j][0].value_in(generic_unit_system.length)
                py = y[i][j][0].value_in(generic_unit_system.length)
                potential =  (math.sin(py * math.pi)+math.sin(px *math.pi)) / -20.0
                if px < 0 or px > 1.0:
                    potential = 0.0
                if py < 0 or py > 1.0:
                    potential = 0.0
                potential_grid.potential[i][j][0] = generic_unit_system.energy.new_quantity(potential)
        print potential_grid.potential
        channel = potential_grid.new_channel_to(instance.potential_grid)
        channel.copy()
        result = instance.initialize_grid() 
        
        self.assertEquals(instance.grid[1][1][0].rho, 0.1 | density)
        for x in instance.grid[1].rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
            
        instance.evolve(1.0 | generic_unit_system.time)
        print instance.grid.rhox
        print instance.grid.rho
        
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertNotEquals(x, 0.1)
        
            
            