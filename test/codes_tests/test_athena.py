import os
import sys
import numpy
import math

from amuse.test.amusetest import TestWithMPI
from amuse.community.athena.interface import AthenaInterface, Athena



from mpi4py import MPI
from amuse.units import generic_unit_system
from amuse.units import units
from amuse.support import data
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
    
    
class TestAthena(TestWithMPI):
    
    def test0(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.parameters.nx = 10
        instance.parameters.ny = 20
        instance.parameters.nz = 40
        instance.parameters.length_x = 1 | generic_unit_system.length
        instance.parameters.length_y = 2 | generic_unit_system.length
        instance.parameters.length_z = 3 | generic_unit_system.length
        instance.set_boundary("periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string)
        
        result = instance.commit_parameters()
        
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
        
        
        potential_grid = data.Grid(12,22,42)
        potential_grid.potential = 2.0 | generic_unit_system.potential
        channel = potential_grid.new_channel_to(instance.potential_grid)
        channel.copy()
        self.assertEquals(instance.potential_grid[0][0][0].potential, 2.0 | generic_unit_system.potential)
        self.assertEquals(instance.potential_grid[0][2][20].potential, 2.0 | generic_unit_system.potential)
        
        instance.stop()
    
    def test1(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
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
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.setup_mesh(10 | units.none, 10 | units.none, 1 | units.none, 1.0 | generic_unit_system.length, 1.0 | generic_unit_system.length , 0.0 | generic_unit_system.length)
        instance.parameters.mesh_size = (10,10,1)
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 0.0 | generic_unit_system.length
        instance.set_boundary("periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string)
        
        result = instance.commit_parameters()
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
    
        grid = data.Grid(10,10,1)
        grid.rho = 0.1 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
            
        result = instance.initialize_grid()
        
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
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.parameters.nx = 10
        instance.parameters.ny = 10
        instance.parameters.nz = 1
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 0.0 | generic_unit_system.length
        instance.set_boundary("periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string)
        instance.set_has_external_gravitational_potential(1)
        
        result = instance.commit_parameters()
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
    
        grid = data.Grid(10,10,1)
        grid.rho = 0.1 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
        
        self.assertEquals(grid._get_writeable_attribute_names(), ['energy', 'rho', 'rhovx', 'rhovy', 'rhovz', ] )
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        potential_grid = data.Grid(12,12,1)
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
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.parameters.nx = 10
        instance.parameters.ny = 10
        instance.parameters.nz = 1
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 0.0 | generic_unit_system.length
        instance.set_boundary("periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string)
        instance.set_has_external_gravitational_potential(1)
        
        result = instance.commit_parameters()
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
    
        grid = data.Grid(10,10,1)
        grid.rho = 0.1 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        grid.energy = 0.0 | energy
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        potential_grid = data.Grid(12,12,1)
        potential_grid.potential = 0.0 | generic_unit_system.potential
        print instance.potential_grid.shape
        print instance.potential_grid.x[0].number.shape
        print instance.grid.shape
        print instance.grid.x[0].number.shape
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
        print potential_grid.potential
        #channel = potential_grid.new_channel_to(instance.potential_grid)
        #channel.copy()
        result = instance.initialize_grid()
        
        self.assertEquals(instance.grid[1][1][0].rho, 0.1 | density)
        for x in instance.grid[1].rho.value_in(density).flatten():
            self.assertEquals(x, 0.1)
            
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
    
    

    def test5(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        self.assertAlmostRelativeEquals(instance.parameters.isothermal_sound_speed, 0.0 | generic_unit_system.speed)
        instance.parameters.isothermal_sound_speed = 0.1 | generic_unit_system.speed
        self.assertAlmostRelativeEquals(instance.parameters.isothermal_sound_speed, 0.1 | generic_unit_system.speed)
        self.assertAlmostRelativeEquals(instance.parameters.gamma, 1.66666666666666667 | units.none)
        instance.parameters.gamma = 0.1 | units.none
        self.assertAlmostRelativeEquals(instance.parameters.gamma, 0.1 | units.none)
        self.assertAlmostRelativeEquals(instance.parameters.courant_number, 0.3 | units.none)
        instance.parameters.courant_number = 0.1 | units.none
        self.assertAlmostRelativeEquals(instance.parameters.courant_number, 0.1 | units.none)
        
        print instance.parameters
        instance.stop()
    
    

    def test6(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.set_gamma(1.6666666666666667)
        instance.set_courant_friedrichs_lewy_number(0.3)
        instance.setup_mesh(10 | units.none , 20 | units.none, 40 | units.none, 1.0 | generic_unit_system.length, 1.0 | generic_unit_system.length, 1.0 | generic_unit_system.length)
        instance.set_boundary("periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string,"periodic" | units.string)
        result = instance.commit_parameters()
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        grid = data.Grid.create((10,10,10), [10.0, 10.0, 10.0] | units.m)
        
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
        instance.initialize_code()
        instance.parameters.isothermal_sound_speed = 0.1 | generic_unit_system.speed
        instance.parameters.gamma = 0.1 | units.none
        instance.parameters.courant_number = 0.1 | units.none
       
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
        self.assertEquals(instance.parameters.mesh_size, (10,20,40) | units.none)
        print instance.parameters
        instance.stop()
    
    def test8(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.parameters.stopping_conditions_number_of_steps = 10
        self.assertEquals(instance.parameters.stopping_conditions_number_of_steps, 10|units.none)
        instance.stop()

    def test8a(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.parameters.stopping_conditions_timeout = 10 | units.s
        self.assertEquals(instance.parameters.stopping_conditions_timeout, 10|units.s)
        instance.stop()

    def test9(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        self.assertEquals(instance.parameters.xbound1, "periodic" | units.string)
        instance.stop()

    def xtest10(self):
        instance=self.new_instance(Athena)
        instance.initialize_code()
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
        
        grid = data.Grid(n,n,n)
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
        
    
    def test11(self):
        instance=self.new_instance(Athena, mode=AthenaInterface.MODE_SELF_GRAVITY) #, redirection = "none") #, debugger="gdb")
        instance.initialize_code()
        instance.parameters.gamma = 5/3.0
        instance.parameters.courant_number=0.3
        
        instance.set_four_pi_G( 4 * numpy.pi ) # G = 1, like nbody
        instance.set_grav_mean_rho(0.0)
        

        data.Grid.add_global_vector_attribute("position", ["x","y","z"])
        n = 10
        
        density = generic_unit_system.mass / (generic_unit_system.length ** 3)
        momentum =  generic_unit_system.mass / (generic_unit_system.time * (generic_unit_system.length**2))
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
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
    
        grid = data.Grid.create((n,n,n), [4.0 , 4.0, 4.0] | generic_unit_system.length)
        
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
                
        
        
