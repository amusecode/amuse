import os
import sys
import numpy


from amuse.test.amusetest import TestWithMPI
from amuse.community.capreole.interface import CapreoleInterface
from amuse.community.capreole.interface import Capreole
from amuse.units import generic_unit_system
from amuse.units import generic_unit_converter
from amuse.units import units
from amuse import datamodel

class TestCapreoleInterface(TestWithMPI):
    
    def test0(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        instance.stop()
        
    
    def test1(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        instance.setup_mesh(50,40,30,1.,1.,1.)
        nx,ny,nz, error = instance.get_mesh_size()
        self.assertEquals(nx, 50)
        self.assertEquals(ny, 40)
        self.assertEquals(nz, 30)
        instance.stop()
    
    def test2(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        instance.setup_mesh(50,50,50,1.,1.,1.)
        instance.commit_parameters()
        instance.set_grid_state(1,1,1,1.,0.1,0.1,0.1,1.)
        rho, err=instance.get_grid_density(1,1,1)
        self.assertEqual(rho,1.)
        rhovx,rhovy,rhovz,err=instance.get_grid_momentum_density(1,1,1)
        self.assertEqual(rhovx,0.1)
        self.assertEqual(rhovy,0.1)
        self.assertEqual(rhovz,0.1)
        en,err=instance.get_grid_energy_density(1,1,1)
        self.assertEqual(en,1.)
        rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(1,1,1)
        self.assertEqual(rho,1.)
        self.assertEqual(rhovx,0.1)
        self.assertEqual(rhovy,0.1)
        self.assertEqual(rhovz,0.1)
        self.assertEqual(en,1.)
        instance.stop()
    
    def test3(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        instance.setup_mesh(50,50,50,1.,1.,1.)
        instance.commit_parameters()
        x,y,z=numpy.indices( (50,50,50) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=0.25*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=2.*numpy.ones_like(x)
        instance.set_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(1,1,1)
        self.assertEqual(rho,0.25)
        self.assertEqual(rhovx,0.)
        self.assertEqual(rhovy,0.)
        self.assertEqual(rhovz,0.)
        self.assertEqual(en,2.)
        instance.stop()
    
    def test4(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        instance.setup_mesh(50,50,50,1.,1.,1.)
        instance.commit_parameters()
        x,y,z=numpy.indices( (50,50,50) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=0.1*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=0.1*numpy.ones_like(x)
        instance.set_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        instance.initialize_grid()
        instance.stop()
    
    def test5(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        instance.setup_mesh(40,40,40,1.,1.,1.)
        instance.commit_parameters()
        x,y,z=numpy.indices( (40,40,40) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=0.1*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=0.1*numpy.ones_like(x)
        instance.set_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        instance.initialize_grid()
        instance.evolve_model(0.01)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.01,15)
        instance.evolve_model(0.025)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.025,15)
        instance.evolve_model(0.025001)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.025001,15)
        instance.evolve_model(0.0321)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.0321,15)
        instance.evolve_model(0.0321)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.0321,15)
        instance.evolve_model(0.07)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.07,15)
        instance.stop()
    
    def test6(self):
        instance=CapreoleInterface(number_of_workers=1)
        instance.initialize_code()
        instance.setup_mesh(150,10,30,1.,1.,1.)
        instance.commit_parameters()
        x,y,z=numpy.indices( (150,10,30) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=0.1*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=0.1*numpy.ones_like(x)
        instance.set_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        instance.initialize_grid()
        instance.evolve_model(0.01)
        x,y,z,err=instance.get_position_of_index(15,5,20)
        self.assertAlmostEqual(x,15/150.-1/300.,15)
        self.assertAlmostEqual(y,5/10.-1/20.,15)
        self.assertAlmostEqual(z,20/30.-1/60.,15)
        i,j,k,err=instance.get_index_of_position(x,y,z)
        self.assertEqual([i,j,k],[15,5,20])
        instance.stop()

    def test7(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        instance.setup_mesh(50,50,50,1.,1.,1.)
        instance.commit_parameters()
        err=instance.set_gravity_field(1,2,3,1.,0.5,0.25)
        self.assertEqual(err,0)
        fx,fy,fz,err=instance.get_gravity_field(1,2,3)
        self.assertEqual(fx,1.)
        self.assertEqual(fy,0.5)
        self.assertEqual(fz,0.25)
        instance.stop()

    def test8(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        err=instance.set_boundary("periodic","reflective",
        "periodic","reflective",
        "periodic","reflective")
        self.assertEqual(err,-1)
        instance.stop()
    
        instance=CapreoleInterface()
        instance.initialize_code()
        err=instance.set_boundary("reflective","periodic",
        "periodic","reflective",
        "periodic","reflective")
        self.assertEqual(err,-2)
        instance.stop()
    
        instance=CapreoleInterface()
        instance.initialize_code()
        err=instance.set_boundary("periodic","periodic",
        "periodic","periodic",
        "periodic","periodic")
        self.assertEqual(err,0)
        instance.stop()

    def test9(self):
        instance=CapreoleInterface(number_of_workers=2)
        instance.initialize_code()
        instance.setup_mesh(50,50,50,1.,1.,1.)
        instance.commit_parameters()
        instance.set_grid_state(1,1,1,1.,0.1,0.1,0.1,1.)
        instance.set_grid_state(50,50,50,2.,0.2,0.2,0.2,2.)

        rho, err=instance.get_grid_density(1,1,1)
        self.assertEqual(rho,1.)
        rhovx,rhovy,rhovz,err=instance.get_grid_momentum_density(1,1,1)
        self.assertEqual(rhovx,0.1)
        self.assertEqual(rhovy,0.1)
        self.assertEqual(rhovz,0.1)
        en,err=instance.get_grid_energy_density(1,1,1)
        self.assertEqual(en,1.)
        rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(1,1,1)
        self.assertEqual(rho,1.)
        self.assertEqual(rhovx,0.1)
        self.assertEqual(rhovy,0.1)
        self.assertEqual(rhovz,0.1)
        self.assertEqual(en,1.)

        rho, err=instance.get_grid_density(50,50,50)
        self.assertEqual(err,0)
        self.assertEqual(rho,2.)
        rhovx,rhovy,rhovz,err=instance.get_grid_momentum_density(50,50,50)
        self.assertEqual(rhovx,0.2)
        self.assertEqual(rhovy,0.2)
        self.assertEqual(rhovz,0.2)
        en,err=instance.get_grid_energy_density(50,50,50)
        self.assertEqual(en,2.)
        rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(50,50,50)
        self.assertEqual(rho,2.)
        self.assertEqual(rhovx,0.2)
        self.assertEqual(rhovy,0.2)
        self.assertEqual(rhovz,0.2)
        self.assertEqual(en,2.)


        instance.stop()
        
    def test10(self):
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.setup_mesh(100,5,5,100.0,0,0)
        instance.set_boundary("interface","periodic","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        minx, maxx, miny, maxy, minz, maxz, error = instance.get_boundary_index_range_inclusive(1)
        self.assertEquals(error, 0)
        self.assertEquals(minx, 1)
        self.assertEquals(maxx, 2)
        self.assertEquals(miny, 1)
        self.assertEquals(maxy, 5)
        self.assertEquals(minz, 1)
        self.assertEquals(maxz, 5)
        
        for i in range(2,7):
            minx, maxx, miny, maxy, minz, maxz, error = instance.get_boundary_index_range_inclusive(i)
            self.assertEquals(error, 0)
            self.assertEquals(minx, 1)
            self.assertEquals(maxx, 1)
            self.assertEquals(miny, 1)
            self.assertEquals(maxy, 1)
            self.assertEquals(minz, 1)
            self.assertEquals(maxz, 1)
    
    def test11(self):
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.setup_mesh(100,5,6,100.0,0,0)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
        
        for i in range(1,7):
            minx, maxx, miny, maxy, minz, maxz, error = instance.get_boundary_index_range_inclusive(i)
            self.assertEquals(error, 0),
            self.assertEquals(minx, 1)
            self.assertEquals(miny, 1)
            self.assertEquals(minz, 1)
            if i == 1 or i == 2:
                self.assertEquals(maxx, 2)
                self.assertEquals(maxy, 5)
                self.assertEquals(maxz, 6)
            elif i == 3 or i == 4:
                self.assertEquals(maxx, 100+4)
                self.assertEquals(maxy, 2)
                self.assertEquals(maxz, 6)
            elif i == 5 or i == 6:
                self.assertEquals(maxx, 100+4)
                self.assertEquals(maxy, 5+4)
                self.assertEquals(maxz, 2)
    
    def test12(self):
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.setup_mesh(100,2,2,100.0,100.0,100.0)
        instance.set_boundary("interface","periodic","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        for i in [1,2]:
            error = instance.set_boundary_state(
                i,1,1,       #  index
                1.0 * (i+1),         #  density
                2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                5.0 * (i+1),         #  energy
                1                    #  boundary
            )
            self.assertEquals(error, 0)
            rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                i, 1, 1,
                1
            )
            print rho, rhovx, rhovy, rhovz, rhoen, error 
            self.assertEquals(error, 0)
            self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
            self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
            self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
            self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
            self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    def test13(self):
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.setup_mesh(100,2,2,100.0,0,0)
        instance.set_boundary("interface","interface","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        for i in [1,2]:
            for j in [1,2]:
                error = instance.set_boundary_state(
                    i,1,1,       #  index
                    1.0 * (i+1),         #  density
                    2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                    5.0 * (i+1),         #  energy
                    j    #  boundary 
                )
                self.assertEquals(error, 0)
                rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                    i, 1,1,
                    j
                )
                print j
                self.assertEquals(error, 0)
                
                self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    def test14(self):
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.setup_mesh(5,6,7,100.0,100.0,100.0)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
        x1range = (2,6,7)
        x2range = (5+4,2,7)
        x3range = (5+4,6+4,2)
    
        for xrange, j in zip([x1range, x1range, x2range, x2range, x3range, x3range], [1,2,3,4,5,6]):
            for i0 in range(xrange[0]):
                for j0 in range(xrange[1]):
                    for k0 in range(xrange[2]):
                        i = (i0 * (xrange[2] * xrange[1])) + (j0 * xrange[2]) + k0
                        print "boundary:", j, i0+1, j0+1, k0+1
                        error = instance.set_boundary_state(
                            i0+1, j0+1, k0+1,       #  index
                            1.0 * (i+1),         #  density
                            2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                            5.0 * (i+1),         #  energy
                            j
                        )
                        self.assertEquals(error, 0)
                        rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                            i0+1, j0+1, k0+1,       #  index
                            j
                        )
                        self.assertEquals(error, 0)
                        
                        self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    def test15(self):
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.setup_mesh(100,5,4,100.0,100.0, 100.0)
        instance.set_boundary("interface","interface","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        dx = 100.0 / 100.0
        dy = 100.0 / 5.0
        dz = 100.0 / 4.0
        
        for i in [1,2]:
            x,y,z,error = instance.get_boundary_position_of_index(
                i, 1, 1,
                1
            )
            self.assertEquals(error, 0)
            self.assertAlmostRelativeEquals(x, (0.5 * dx) - (i * dx))
            self.assertAlmostRelativeEquals(y, (0.5 * dy))
            self.assertAlmostRelativeEquals(z, (0.5 * dz))
    
    def test16(self):
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.setup_mesh(100,5,4,100.0,100.0, 100.0)
        instance.set_boundary("interface","interface","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
        dx = 100.0 / 100.0
        dy = 100.0 / 5.0
        dz = 100.0 / 4.0
        
        for i in [1,2]:
            x,y,z,error = instance.get_boundary_position_of_index(
                i,1,1,
                2
            )
            self.assertEquals(error, 0)
            self.assertAlmostRelativeEquals(x, 100.0 + (0.5 * dx) + ((i-1) * dx))
            self.assertAlmostRelativeEquals(y, (0.5 * dy))
            self.assertAlmostRelativeEquals(z, (0.5 * dz))
    
    def test17(self):
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.setup_mesh(100,5,4,100.0,100.0, 100.0)
        instance.set_boundary("interface","interface","interface","interface","periodic","periodic")
        instance.commit_parameters()
        
        dx = 100.0 / 100.0
        dy = 100.0 / 5.0
        dz = 100.0 / 4.0
        
        for i in [1,2]:
            for j in range(1,6):
                x,y,z,error = instance.get_boundary_position_of_index(
                    i, j, 1, 
                    2
                )
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(x, 100.0 + (0.5 * dx) + ((i-1) * dx))
                self.assertAlmostRelativeEquals(y, (0.5 * dy) + ((j-1) * dy))
                self.assertAlmostRelativeEquals(z, (0.5 * dz))
        
        for i in range(1, 100 + 4 + 1):
            for j in [1,2]:
                x,y,z,error = instance.get_boundary_position_of_index(
                    i, j, 1, 
                    3
                )
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-2-1) * dx))
                self.assertAlmostRelativeEquals(y, 0.0 - ((0.5 * dy) + ((j-1) * dy)))
                self.assertAlmostRelativeEquals(z, (0.5 * dz))
                
                
                x,y,z,error = instance.get_boundary_position_of_index(
                    i, j, 1, 
                    4
                )
                self.assertEquals(error, 0)
                self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-2-1) * dx))
                self.assertAlmostRelativeEquals(y, 100.0 + (0.5 * dy) + ((j-1) * dy))
                self.assertAlmostRelativeEquals(z, (0.5 * dz))
    
    def test18(self):
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.setup_mesh(3, 3, 3, 6,12,18)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
        dx = 6.0 / 3.0
        dy = 12.0 / 3.0
        dz = 18.0 / 3.0
        for i in [1,2]:
            for j in range(1,3+1):
                for k in range(1,3+1):
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        2
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, 6.0 + (0.5 * dx) + ((i-1) * dx))
                    self.assertAlmostRelativeEquals(y, (0.5 * dy) + ((j-1) * dy))
                    self.assertAlmostRelativeEquals(z, (0.5 * dz) + ((k-1) * dz))
        
        for i in range(1,3 + 4 +1):
            for j in [1,2]:
                for k in range(1,3+1):
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        3
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-2-1) * dx))
                    self.assertAlmostRelativeEquals(y, 0.0 - ((0.5 * dy) + ((j-1) * dy)))
                    self.assertAlmostRelativeEquals(z, (0.5 * dz) + ((k-1) * dz))
                    
                    
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        4
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-2-1) * dx))
                    self.assertAlmostRelativeEquals(y, 12.0 + (0.5 * dy) + ((j-1) * dy))
                    self.assertAlmostRelativeEquals(z, (0.5 * dz) + ((k-1) * dz))
        
        for i in range(1,3 + 4 +1):
            for j in range(1,3 + 4 +1):
                for k in [1,2]:
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        5 
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-2-1) * dx))
                    self.assertAlmostRelativeEquals(y, (0.5 * dy) + ((j-2-1) * dy))
                    self.assertAlmostRelativeEquals(z,  0.0 - ((0.5 * dz) + ((k-1) * dz)))
                    
                    
                    x,y,z,error = instance.get_boundary_position_of_index(
                        i, j, k, 
                        6
                    )
                    self.assertEquals(error, 0)
                    self.assertAlmostRelativeEquals(x, (0.5 * dx) + ((i-2-1) * dx))
                    self.assertAlmostRelativeEquals(y, (0.5 * dy) + ((j-2-1) * dy))
                    self.assertAlmostRelativeEquals(z, 18.0 + (0.5 * dz) + ((k-1) * dz))
        
    def test19(self):
        results = []
        instance=self.new_instance(CapreoleInterface)
        instance.initialize_code()
        instance.commit_parameters()
        nx, ny, nz, error = instance.get_parallel_decomposition()
        self.assertEquals(error, 0)
        self.assertEquals(nx, 1)
        self.assertEquals(ny, 1)
        self.assertEquals(nz, 1)
        error = instance.set_parallel_decomposition(2,1,1)
        self.assertEquals(error, -1)
        
   
    def test20(self):
        results = []
        instance=self.new_instance(CapreoleInterface, number_of_workers = 4)
        instance.initialize_code()
        nx, ny, nz, error = instance.get_parallel_decomposition()
        self.assertEquals(error, 0)
        self.assertEquals(nx, 0)
        self.assertEquals(ny, 0)
        self.assertEquals(nz, 0)
        error = instance.set_parallel_decomposition(2,1,2)
        self.assertEquals(error, 0)
        nx, ny, nz, error = instance.get_parallel_decomposition()
        self.assertEquals(error, 0)
        self.assertEquals(nx, 2)
        self.assertEquals(ny, 1)
        self.assertEquals(nz, 2)
        error = instance.set_parallel_decomposition(10,3,2)
        self.assertEquals(error, -1)
        
    def test21(self):
        results = []
        instance=self.new_instance(CapreoleInterface, number_of_workers = 2)
        instance.initialize_code()
        error = instance.set_parallel_decomposition(1,2,1)
        self.assertEquals(error, 0)
        instance.setup_mesh(10,30,10,100.0, 300.0, 100.0)
        instance.set_boundary("interface","interface","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
           
        for boundary_index in [1,2]:
            for i0 in range(1,2):
                for j0 in range(1, 30+1):
                    i = j0 * 30 + i0
                    error = instance.set_boundary_state(
                        i0, j0, 1,       #  index
                        1.0 * (i+1),         #  density
                        2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                        5.0 * (i+1),         #  energy
                        boundary_index     #  boundary 
                    )
                    self.assertEquals(error, 0)
                    rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                        i0, j0, 1,
                        boundary_index
                    )
                    self.assertEquals(error, 0)
                    
                    self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    def test22(self):
        results = []
        instance=self.new_instance(CapreoleInterface, number_of_workers = 2)
        instance.initialize_code()
        error = instance.set_parallel_decomposition(2,1,1)
        self.assertEquals(error, 0)
        instance.setup_mesh(10,30,10,100.0, 300.0, 100.0)
        instance.set_boundary("interface","interface","periodic","periodic","periodic","periodic")
        instance.commit_parameters()
        
           
        for boundary_index in [1,2]:
            for i0 in range(1,2):
                for j0 in range(1, 30+1):
                    i = j0 * 30 + i0
                    error = instance.set_boundary_state(
                        i0, j0, 1,       #  index
                        1.0 * (i+1),         #  density
                        2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                        5.0 * (i+1),         #  energy
                        boundary_index     #  boundary 
                    )
                    self.assertEquals(error, 0)
                    rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                        i0, j0, 1,
                        boundary_index
                    )
                    self.assertEquals(error, 0)
                    
                    self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    def test23(self):
        results = []
        instance=self.new_instance(CapreoleInterface, number_of_workers = 3)
        instance.initialize_code()
        error = instance.set_parallel_decomposition(3,1,1)
        self.assertEquals(error, 0)
        instance.setup_mesh(12,20,10,100.0, 300.0, 100.0)
        
        instance.set_boundary("interface","interface","interface","interface","periodic","periodic")
        instance.commit_parameters()
        
           
        for boundaryindex in [3,4]:
            for i0 in range(1,12+4+1):
                for j0 in [1,2]:
                    i = (i0 * 15) + j0
                    error = instance.set_boundary_state(
                        i0,j0,1,       #  index
                        1.0 * (i+1),         #  density
                        2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                        5.0 * (i+1),         #  energy
                        boundaryindex    #  boundary
                    )
                    print i0, j0
                    self.assertEquals(error, 0)
                    rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                        i0, j0, 1,
                        boundaryindex
                    )
                    self.assertEquals(error, 0)
                    
                    self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    def test24(self):
        results = []
        instance=self.new_instance(CapreoleInterface, number_of_workers = 3)
        instance.initialize_code()
        error = instance.set_parallel_decomposition(1,3,1)
        self.assertEquals(error, 0)
        instance.setup_mesh(12,30,10,100.0, 300.0, 100.0)
        
        instance.set_boundary("interface","interface","interface","interface","periodic","periodic")
        instance.commit_parameters()
        
           
        for boundaryindex in [3,4]:
            for i0 in range(1,12+4+1):
                for j0 in [1,2]:
                    i = (i0 * 15) + j0
                    error = instance.set_boundary_state(
                        i0,j0,1,       #  index
                        1.0 * (i+1),         #  density
                        2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                        5.0 * (i+1),         #  energy
                        boundaryindex    #  boundary
                    )
                    print i0, j0
                    self.assertEquals(error, 0)
                    rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                        i0, j0, 1,
                        boundaryindex
                    )
                    self.assertEquals(error, 0)
                    
                    self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                    self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
                    
    def test25(self):
        results = []
        instance=self.new_instance(CapreoleInterface, number_of_workers = 3)
        instance.initialize_code()
        error = instance.set_parallel_decomposition(1,3,1)
        self.assertEquals(error, 0)
        instance.setup_mesh(6,5,5,6.0,5.0,5.0)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
           
        for boundaryindex in [5,6]:
            for i0 in range(1, 6+4+1):
                for j0 in range(1, 5+4+1):
                    for z0 in[1,2]:
                        i = (i0 * (5*4)) + (j0 * 4) + z0
                        error = instance.set_boundary_state(
                            i0,j0,z0,       #  index
                            1.0 * (i+1),         #  density
                            2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                            5.0 * (i+1),         #  energy
                            boundaryindex     #  boundary 
                        )
                        self.assertEquals(error, 0)
                        rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                            i0, j0, z0,
                            boundaryindex
                        
                        )
                        self.assertEquals(error, 0)
                        
                        self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
    
    def test26(self):
        results = []
        instance=self.new_instance(CapreoleInterface, number_of_workers = 9)
        instance.initialize_code()
        error = instance.set_parallel_decomposition(3,3,1)
        instance.setup_mesh(6,6,5,6.0,6.0,5.0)
        self.assertEquals(error, 0)
        instance.set_boundary("interface","interface","interface","interface","interface","interface")
        instance.commit_parameters()
        
           
        for boundaryindex in [5,6]:
            for i0 in range(1,6+4+1):
                for j0 in range(1,6+4+1):
                    for z0 in [1,2]:
                        i = (i0 * (5*4)) + (j0 * 4) + z0
                        error = instance.set_boundary_state(
                            i0,j0,z0,            #  index
                            1.0 * (i+1),         #  density
                            2.0 * (i+1), 3.0 * (i+1), 4.0 * (i+1), #  momentum
                            5.0 * (i+1),         #  energy
                            boundaryindex        #  boundary 
                        )
                        self.assertEquals(error, 0)
                        rho, rhovx, rhovy, rhovz, rhoen, error = instance.get_boundary_state(
                            i0, j0, z0,
                            boundaryindex
                        )
                        self.assertEquals(error, 0)
                        
                        self.assertAlmostRelativeEquals(rho, 1.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovx, 2.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovy, 3.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhovz, 4.0 * (i+1))
                        self.assertAlmostRelativeEquals(rhoen, 5.0 * (i+1))
                    
    def test27(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        gamma, error = instance.get_gamma()
        self.assertAlmostRelativeEquals(gamma, 5.0 / 3.0)
        instance.set_gamma(1.3)
        gamma, error = instance.get_gamma()
        self.assertEquals(error, 0)
        self.assertAlmostRelativeEquals(gamma, 1.3)
        instance.stop()
        
    
    def test28(self):
        instance=CapreoleInterface()
        instance.initialize_code()
        timestep, error = instance.get_timestep()
        self.assertEquals(timestep, 0.0)
        self.assertEquals(error, 0)
        instance.setup_mesh(10,10,10,1.,1.,1.)
        instance.commit_parameters()
        x,y,z=numpy.indices( (10,10,10) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=0.1*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=0.1*numpy.ones_like(x)
        instance.set_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        timestep, error = instance.get_timestep()
        self.assertEquals(timestep, 0.0)
        self.assertEquals(error, 0)
        instance.initialize_grid()
        timestep, error = instance.get_timestep()
        self.assertAlmostRelativeEquals(timestep, 0.03795, 4)
        self.assertEquals(error, 0)
        error = instance.set_timestep(0.25)
        self.assertEquals(error, 0)
        timestep, error = instance.get_timestep()
        self.assertEquals(timestep, 0.25)
        self.assertEquals(error, 0)
        
        instance.stop()
        
class TestSodShocktube(TestWithMPI):
    
    def test0(self):
        N=100
        gamma=5/3.
        g=(gamma-1)/(gamma+1)
        b=(gamma-1)/2/gamma
        
        instance=CapreoleInterface()
        instance.initialize_code()
        instance.setup_mesh(N,N/10,N/10,1.,0.1,0.1)
        instance.commit_parameters()
        x,y,z=numpy.indices( (N,N/10,N/10) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        gamma=5./3.
        rho=0.125*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=(0.1/(gamma-1))*numpy.ones_like(x)
        instance.set_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
    
        x,y,z=numpy.indices( (N/2,N/10,N/10) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=1.*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=(1./(gamma-1))*numpy.ones_like(x)
        instance.set_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        instance.initialize_grid()
        instance.evolve_model(0.2)
    
        x=numpy.array([0.1,0.9,0.6,0.8])
        y=0.05*numpy.ones_like(x)
        z=0.05*numpy.ones_like(x)
        i,j,k,err=instance.get_index_of_position(x,y,z)
        rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(i,j,k)
        vel=numpy.sqrt((rhovx**2+rhovy**2+rhovz**2))/rho
        pres=(gamma-1)*(en-0.5*rho*vel**2)
        u=pres/(gamma-1)/rho
        rhoexp=numpy.zeros_like(x)
        rhoexp[0]=1.
        rhoexp[1]=0.125
        rhoexp[2]=rho[0]*(pres[2]/pres[0])**(1/gamma)
        rhoexp[3]=rho[1]*(pres[3]+g*pres[1])/(pres[1]+g*pres[3])
        
        for i in range(len(rho)):
            self.assertAlmostEqual(rhoexp[i],rho[i],2)
            
    
    

class TestCapreole(TestWithMPI):
    
    def test0(self):
        instance=self.new_instance(Capreole)
        instance.initialize_code()
        instance.stop()
        
    def test1(self):
        instance=self.new_instance(Capreole)
        instance.parameters.mesh_size = (10,10,5)
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 1.0 | generic_unit_system.length
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        
    
        self.assertEquals(len(list(instance.itergrids())),1)
        grid = datamodel.Grid(10,10,10)
        grid.rho = 0.4 | generic_unit_system.density
        grid.rhovx = 0.1 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy = 0.0 | generic_unit_system.energy_density
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
            
        instance.initialize_grid()
        
        channel = instance.grid.new_channel_to(grid)
        channel.copy()
        
        self.assertEquals(grid[1][1][0].rho, 0.4 | generic_unit_system.density)
        for x in grid[1].rho.value_in(generic_unit_system.density).flatten():
            self.assertEquals(x, 0.4)
            
        #instance.evolve_model(0.12 | generic_unit_system.time)
        
        #for x in instance.grid.rho.value_in(generic_unit_system.density).flatten():
        #    self.assertEquals(x, 0.1)
    
        #instance.evolve_model(10.0 | generic_unit_system.time)
        #for x in instance.grid.rho.value_in(generic_unit_system.density).flatten():
        #    self.assertEquals(x, 0.1)
        instance.stop()


    def test2(self):
        instance=self.new_instance(Capreole)
        instance.parameters.mesh_size = (3,3,3)
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 1.0 | generic_unit_system.length
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        
        grid = datamodel.Grid(3,3,3)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy = 1.0 | generic_unit_system.energy_density
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
            
        
        print instance.grid[1].rho
        self.assertEquals(instance.grid[1][1][0].rho, 0.1 | generic_unit_system.density)
        for x in instance.grid[1].rho.value_in(generic_unit_system.density).flatten():
            self.assertEqual(x, 0.1)
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        
        for x in instance.grid.rho.value_in(generic_unit_system.density).flatten():
            self.assertEquals(x, 0.1)
    
        instance.evolve_model(10.0 | generic_unit_system.time)
        for x in instance.grid.rho.value_in(generic_unit_system.density).flatten():
            self.assertEquals(x, 0.1)
        instance.stop()
    
    def test3(self):
        instance=self.new_instance(Capreole)
        instance.parameters.mesh_size = (5,5,5)
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 1.0 | generic_unit_system.length
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        
    
        grid = datamodel.Grid(5,5,5)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy =  1.0 | generic_unit_system.energy_density
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        
        self.assertEquals((5,5,5), instance.acceleration_grid.shape)
        
        acc_grid = datamodel.Grid(5,5,5)
        acc_grid.ax = 1 | generic_unit_system.acceleration
        acc_grid.ay = 1 | generic_unit_system.acceleration
        acc_grid.az = 1 | generic_unit_system.acceleration
        #self.assertEquals(acc_grid.acceleration[0][0][0], ( 1,1,1) | generic_unit_system.acceleration)
        channel = acc_grid.new_channel_to(instance.acceleration_grid)
        channel.copy()
        
                   
        instance.evolve_model(0.1 | generic_unit_system.time)
        
        self.assertAlmostRelativeEquals(instance.grid.rho, grid.rho);
        self.assertAlmostRelativeEquals(instance.grid.rhovx, 0.1 * 1.0 * 0.1 | generic_unit_system.momentum_density);
        self.assertAlmostRelativeEquals(instance.grid.rhovy, 0.1 * 1.0 * 0.1 | generic_unit_system.momentum_density);
        self.assertAlmostRelativeEquals(instance.grid.rhovz, 0.1 * 1.0 * 0.1 | generic_unit_system.momentum_density);

        instance.evolve_model(0.3 | generic_unit_system.time)
        print instance.model_time
        self.assertAlmostRelativeEquals(instance.grid.rho, grid.rho);
        self.assertAlmostRelativeEquals(instance.grid.rhovx, grid.rho *  instance.model_time * acc_grid.ax,2);
        self.assertAlmostRelativeEquals(instance.grid.rhovy, grid.rho *  instance.model_time * acc_grid.ay,2);
        self.assertAlmostRelativeEquals(instance.grid.rhovz, grid.rho *  instance.model_time * acc_grid.az,2);
        instance.stop()
    
    def test4(self):
        converter = generic_unit_converter.ConvertBetweenGenericAndSiUnits(
            1 | units.parsec,
            1 | units.Myr,
            1 | units.MSun
        )
        instance=self.new_instance(Capreole, unit_converter = converter)
        instance.parameters.mesh_size = (3,3,3)
        instance.parameters.length_x = 1.0 | units.parsec
        instance.parameters.length_y = 1.0 | units.parsec
        instance.parameters.length_z = 1.0 | units.parsec
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        
        instance.commit_parameters()
        density = units.MSun / (units.parsec ** 3)
        grid = datamodel.Grid(3,3,3)
        grid.rho = 0.1 | density
        grid.rhovx = 0.0 | units.MSun / (units.Myr * units.parsec ** 2 )
        grid.rhovy = 0.0 | units.MSun / (units.Myr * units.parsec ** 2 )
        grid.rhovz = 0.0 | units.MSun / (units.Myr * units.parsec ** 2 )
        grid.energy = 1.0 | units.MSun / (units.parsec * units.Myr ** 2)
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
            
        
        print instance.grid[1].rho
        self.assertAlmostRelativeEquals(instance.grid[1][1][0].rho, 0.1 | density)
        for x in instance.grid[1].rho.value_in(density).flatten():
            self.assertAlmostRelativeEquals(x, 0.1)
            
        instance.evolve_model(1.0 | units.Myr)
        
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertAlmostRelativeEquals(x, 0.1)
    
        instance.evolve_model(10.0 | units.Myr)
        for x in instance.grid.rho.value_in(density).flatten():
            self.assertAlmostRelativeEquals(x, 0.1)
        instance.stop()
    
    
         
    def test5(self): 
        instance=self.new_instance(Capreole)
        instance.parameters.mesh_size = (10 , 4, 4)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("interface", "outflow")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("periodic", "periodic")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        gamma = 5.0 / 3.0
        
        grid = datamodel.Grid.create((10,4,4), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.1 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        #instance.grid.boundaries.left.
        
        xbound1 = instance.get_boundary_grid('xbound1')
        self.assertEquals(xbound1.shape, (2,4,4))
        memxbound1 = xbound1.copy()
        memxbound1.rho = 0.02 | density
        memxbound1.rhovx = 0.2 | momentum
        memxbound1.rhovy = 0.0 | momentum
        memxbound1.rhovz = 0.0 | momentum
        memxbound1.energy =  p / (gamma - 1)
        memxbound1.energy += 0.5 * (memxbound1.rhovx ** 2  + memxbound1.rhovy ** 2 + memxbound1.rhovz ** 2) / memxbound1.rho
        channel = memxbound1.new_channel_to(xbound1)
        channel.copy()
        
        instance.evolve_model(1.0 | generic_unit_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        
        rho = instance.grid.rho[...,0,0]
        print rho
        print instance.model_time
        self.assertAlmostRelativeEquals(rho[-1], 0.01 | density)
        self.assertTrue(rho[0] > 0.01 | density)
        self.assertTrue(instance.grid.rhovx[0,0,0] > 0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[-1,0,0] , 0.1 | momentum)
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        print instance.model_time
        rho = instance.grid.rho[...,0,0]
        print rho
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[...,0,0], 0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()
    
    def test6(self): 
        instance=self.new_instance(Capreole)
        instance.parameters.mesh_size = (10 , 4, 4)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("outflow", "interface")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("periodic", "periodic")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        gamma = 5.0 / 3.0
        
        grid = datamodel.Grid.create((10,4,4), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = -0.1 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        #instance.grid.boundaries.left.
        xbound = instance.get_boundary_grid('xbound2')
        self.assertEquals(xbound.shape, (2,4,4))
        memxbound = xbound.copy()
        memxbound.rho = 0.02 | density
        memxbound.rhovx = -0.2 | momentum
        memxbound.rhovy = 0.0 | momentum
        memxbound.rhovz = 0.0 | momentum
        memxbound.energy =  p / (gamma - 1)
        memxbound.energy += 0.5 * (memxbound.rhovx ** 2  + memxbound.rhovy ** 2 + memxbound.rhovz ** 2) / memxbound.rho
        channel = memxbound.new_channel_to(xbound)
        channel.copy()
        
        instance.evolve_model(1.0 | generic_unit_system.time)
        
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        rho = instance.grid.rho[...,0,0]
        print rho
        print instance.model_time
        self.assertAlmostRelativeEquals(rho[0], 0.01 | density)
        self.assertTrue(rho[-1] > 0.01 | density)
        self.assertTrue(instance.grid.rhovx[-1,0,0] < -0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[0,0,0] , -0.1 | momentum)
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[...,0,0]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[...,0,0], -0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()
    
    def test7(self): 
        instance=self.new_instance(Capreole, number_of_workers = 2)
        instance.set_parallel_decomposition(1,2,1)
        instance.parameters.mesh_size = (10,4,4)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("interface", "outflow")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        gamma = 5.0 / 3.0
        
        grid = datamodel.Grid.create((10,4,4), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.1 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.0 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        #instance.grid.boundaries.left.
        
        xbound1 = instance.get_boundary_grid('xbound1')
        self.assertEquals(xbound1.shape, (2,4,4))
        memxbound1 = xbound1.copy()
        memxbound1.rho = 0.02 | density
        memxbound1.rhovx = 0.2 | momentum
        memxbound1.rhovy = 0.0 | momentum
        memxbound1.rhovz = 0.0 | momentum
        memxbound1.energy =  p / (gamma - 1)
        memxbound1.energy += 0.5 * (memxbound1.rhovx ** 2  + memxbound1.rhovy ** 2 + memxbound1.rhovz ** 2) / memxbound1.rho
        channel = memxbound1.new_channel_to(xbound1)
        channel.copy()
        
        instance.evolve_model(1.0 | generic_unit_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        
        rho = instance.grid.rho[...,0,0]
        print rho
        print instance.model_time
        self.assertAlmostRelativeEquals(rho[-1], 0.01 | density)
        self.assertTrue(rho[0] > 0.01 | density)
        self.assertTrue(instance.grid.rhovx[0,0,0] > 0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[-1,0,0] , 0.1 | momentum)
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        print instance.model_time
        rho = instance.grid.rho[...,0,0]
        print rho
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovx[...,0,0], 0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()

    def test8(self): 
        instance=self.new_instance(Capreole, number_of_workers = 1)
        #instance.set_parallel_decomposition(2,1,1)
        instance.parameters.mesh_size = (4,10,4)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("interface", "outflow")
        instance.parameters.z_boundary_conditions = ("periodic", "periodic")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        gamma = 5.0 / 3.0
        
        grid = datamodel.Grid.create((4,10,4), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.1 | momentum
        grid.rhovz = 0.0 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        ybound = instance.get_boundary_grid('ybound1')
        self.assertEquals(ybound.shape, (4+4,2,4))
        memybound = ybound.copy()
        memybound.rho = 0.02 | density
        memybound.rhovx = 0.0 | momentum
        memybound.rhovy = 0.2 | momentum
        memybound.rhovz = 0.0 | momentum
        memybound.energy =  p / (gamma - 1)
        memybound.energy += 0.5 * (memybound.rhovx ** 2  + memybound.rhovy ** 2 + memybound.rhovz ** 2) / memybound.rho
        
        channel = memybound.new_channel_to(ybound)
        channel.copy()
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        print instance.stopping_conditions.number_of_steps_detection.is_set()
        
        print instance.grid.rho[0,...,0]
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

    
    def test9(self): 
        instance=self.new_instance(Capreole, number_of_workers = 1)
        #instance.set_parallel_decomposition(2,1,1)
        instance.parameters.mesh_size = (4,10,4)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("outflow", "interface")
        instance.parameters.z_boundary_conditions = ("periodic", "periodic")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        gamma = 5.0 / 3.0
        
        grid = datamodel.Grid.create((4,10,4), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = -0.1 | momentum
        grid.rhovz = 0.0 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        ybound = instance.get_boundary_grid('ybound2')
        self.assertEquals(ybound.shape, (4+4,2,4))
        memybound = ybound.copy()
        memybound.rho = 0.02 | density
        memybound.rhovx = 0.0 | momentum
        memybound.rhovy = -0.2 | momentum
        memybound.rhovz = 0.0 | momentum
        memybound.energy =  p / (gamma - 1)
        memybound.energy += 0.5 * (memybound.rhovx ** 2  + memybound.rhovy ** 2 + memybound.rhovz ** 2) / memybound.rho
        
        channel = memybound.new_channel_to(ybound)
        channel.copy()
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        print instance.stopping_conditions.number_of_steps_detection.is_set()
        print instance.grid.rho[0,...,0]
        rho = instance.grid.rho[0,...,0]
        self.assertAlmostRelativeEquals(rho[0], 0.01 | density)
        self.assertTrue(rho[-1] > 0.01 | density)
        self.assertTrue(instance.grid.rhovy[0,-1,0] < 0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovy[0,0,0] , -0.1 | momentum)
        print instance.model_time
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[0,...,0]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovy[0,...,0], -0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()
    
    def test10(self): 
        instance=self.new_instance(Capreole, number_of_workers = 1)
        #instance.set_parallel_decomposition(2,1,1)
        instance.parameters.mesh_size = (4,4,10)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("interface", "outflow")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        gamma = 5.0 / 3.0
        
        grid = datamodel.Grid.create((4,4,10), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = 0.1 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        ybound = instance.get_boundary_grid('zbound1')
        self.assertEquals(ybound.shape, (4+4,4+4,2))
        memybound = ybound.copy()
        memybound.rho = 0.02 | density
        memybound.rhovx = 0.0 | momentum
        memybound.rhovy = 0.0 | momentum
        memybound.rhovz = 0.2 | momentum
        memybound.energy =  p / (gamma - 1)
        memybound.energy += 0.5 * (memybound.rhovx ** 2  + memybound.rhovy ** 2 + memybound.rhovz ** 2) / memybound.rho
        
        channel = memybound.new_channel_to(ybound)
        channel.copy()
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        
        rho = instance.grid.rho[0,0,...]
        print rho
        self.assertAlmostRelativeEquals(rho[-1], 0.01 | density)
        self.assertTrue(rho[0] > 0.01 | density)
        self.assertTrue(instance.grid.rhovz[0,0,0] > 0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovz[0,0,-1] , 0.1 | momentum)
        print instance.model_time
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[0,0,...]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovz[0,...,0], 0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()

    
    def test11(self): 
        instance=self.new_instance(Capreole, number_of_workers = 2)
        instance.set_parallel_decomposition(2,1,1)
        instance.parameters.mesh_size = (4,4,10)
        instance.parameters.mesh_length = [1.0, 1.0, 1.0] | generic_unit_system.length
        instance.parameters.x_boundary_conditions = ("periodic", "periodic")
        instance.parameters.y_boundary_conditions = ("periodic", "periodic")
        instance.parameters.z_boundary_conditions = ("outflow", "interface")
        instance.parameters.stopping_conditions_number_of_steps = 1
        
        
        gamma = 5.0 / 3.0
        
        grid = datamodel.Grid.create((4,4,10), [1.0, 1.0, 1.0] | generic_unit_system.length )
        
        density = generic_unit_system.density
        momentum =  generic_unit_system.speed * generic_unit_system.density
        energy =  generic_unit_system.mass / ((generic_unit_system.time**2) * generic_unit_system.length)
        
        
        grid.rho = 0.01 | density
        grid.rhovx = 0.0 | momentum
        grid.rhovy = 0.0 | momentum
        grid.rhovz = -0.10 | momentum
        
        p = 1.0 | (generic_unit_system.mass / (generic_unit_system.length * generic_unit_system.time**2))
        
        grid.energy =  p / (gamma - 1)
        grid.energy += 0.5 * (grid.rhovx ** 2  + grid.rhovy ** 2 + grid.rhovz ** 2) / grid.rho
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        instance.stopping_conditions.number_of_steps_detection.enable()
        
        ybound = instance.get_boundary_grid('zbound2')
        self.assertEquals(ybound.shape, (4+4,4+4,2))
        memybound = ybound.copy()
        memybound.rho = 0.02 | density
        memybound.rhovx = 0.0 | momentum
        memybound.rhovy = 0.0 | momentum
        memybound.rhovz = -0.2 | momentum
        memybound.energy =  p / (gamma - 1)
        memybound.energy += 0.5 * (memybound.rhovx ** 2  + memybound.rhovy ** 2 + memybound.rhovz ** 2) / memybound.rho
        
        channel = memybound.new_channel_to(ybound)
        channel.copy()
            
        instance.evolve_model(1.0 | generic_unit_system.time)
        self.assertTrue(instance.stopping_conditions.number_of_steps_detection.is_set())
        rho = instance.grid.rho[0,0,...]
        self.assertAlmostRelativeEquals(rho[0], 0.01 | density)
        self.assertTrue(rho[-1] > 0.01 | density)
        self.assertTrue(instance.grid.rhovz[0,0,-1] < 0.1 | momentum)
        self.assertAlmostRelativeEquals(instance.grid.rhovz[0,0,0] , -0.1 | momentum)
        print instance.model_time
        
        instance.stopping_conditions.number_of_steps_detection.disable()
        instance.evolve_model(1.0 | generic_unit_system.time)
        rho = instance.grid.rho[0,0,...]
        self.assertAlmostRelativeEquals(rho, 0.02 | density, 8)
        self.assertAlmostRelativeEquals(instance.grid.rhovz[0,0,...], -0.2 | momentum, 8)
        print instance.model_time
        
        instance.stop()
        
    
    def test12(self):
        
        instance=self.new_instance(Capreole)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.z_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (20.0, 1, 1) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 2, 2)
        
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
                0.5 | generic_unit_system.length,
                0.5 | generic_unit_system.length
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

    def test13(self):
        
        instance=self.new_instance(Capreole,  number_of_workers=2)
        instance.parameters.x_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.y_boundary_conditions = ("periodic","periodic")
        instance.parameters.mesh_length = (20.0, 20.0, 4) | generic_unit_system.length
        instance.parameters.mesh_length = (20.0, 20.0, 4) | generic_unit_system.length
        instance.parameters.mesh_size = (20, 20, 2)
        
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
            
    def test14(self):
        
        instance=self.new_instance(Capreole, number_of_workers=3)
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
        
        sample = datamodel.Grid.create(
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
            

    
    def test15(self): 
        instance=self.new_instance(Capreole, number_of_workers = 1)
        self.assertAlmostRelativeEquals(instance.parameters.gamma, 5.0 / 3.0)
        instance.parameters.gamma = 1.2
        self.assertAlmostRelativeEquals(instance.parameters.gamma, 1.2)
        #self.assertAlmostRelativeEquals(instance.parameters.timestep, 0.1 | generic_unit_system.time)
        #instance.parameters.timestep = 0.2 | generic_unit_system.time
        #self.assertAlmostRelativeEquals(instance.parameters.timestep, 0.2 | generic_unit_system.time)
