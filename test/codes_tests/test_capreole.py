import os
import sys
import numpy

from amuse.support.units import generic_unit_system
from amuse.support.data import core

from amuse.test.amusetest import TestWithMPI
from amuse.community.capreole.interface import CapreoleInterface
from amuse.community.capreole.interface import Capreole

class TestMPIInterface(TestWithMPI):
    
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
        rho, err=instance.get_density(1,1,1)
        self.assertEqual(rho,1.)
        rhovx,rhovy,rhovz,err=instance.get_momentum_density(1,1,1)
        self.assertEqual(rhovx,0.1)
        self.assertEqual(rhovy,0.1)
        self.assertEqual(rhovz,0.1)
        en,err=instance.get_energy_density(1,1,1)
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
        instance.evolve(0.01)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.01,15)
        instance.evolve(0.025)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.025,15)
        instance.evolve(0.025001)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.025001,15)
        instance.evolve(0.0321)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.0321,15)
        instance.evolve(0.0321)
        tnow,err=instance.get_time()
        self.assertAlmostEqual(tnow,0.0321,15)
        instance.evolve(0.07)
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
        instance.evolve(0.01)
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
        instance.evolve(0.2)
    
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
        instance.initialize_code()
        instance.parameters.mesh_size = (10,10,5)
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 1.0 | generic_unit_system.length
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        
        instance.commit_parameters()
    
        self.assertEquals(len(list(instance.itergrids())),1)
        grid = core.Grid(10,10,10)
        grid.rho = 0.4 | generic_unit_system.density
        grid.rhovx = 0.1 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy = 0.0 | generic_unit_system.energy_density
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
            
        result = instance.initialize_grid()
        
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
        instance.initialize_code()
        instance.parameters.mesh_size = (3,3,3)
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 1.0 | generic_unit_system.length
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        
        instance.commit_parameters()
    
        grid = core.Grid(3,3,3)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy = 1.0 | generic_unit_system.energy_density
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
            
        result = instance.initialize_grid()
        
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
        instance.initialize_code()
        instance.parameters.mesh_size = (5,5,5)
        instance.parameters.length_x = 1.0 | generic_unit_system.length
        instance.parameters.length_y = 1.0 | generic_unit_system.length
        instance.parameters.length_z = 1.0 | generic_unit_system.length
        instance.parameters.x_boundary_conditions = "periodic","periodic"
        instance.parameters.y_boundary_conditions = "periodic","periodic"
        instance.parameters.z_boundary_conditions = "periodic","periodic"
        
        instance.commit_parameters()
    
        grid = core.Grid(5,5,5)
        grid.rho = 0.1 | generic_unit_system.density
        grid.rhovx = 0.0 | generic_unit_system.momentum_density
        grid.rhovy = 0.0 |  generic_unit_system.momentum_density
        grid.rhovz = 0.0 |  generic_unit_system.momentum_density
        grid.energy =  1.0 | generic_unit_system.energy_density
        
        channel = grid.new_channel_to(instance.grid)
        channel.copy()
        
        self.assertEquals((5,5,5), instance.acceleration_grid.shape)
        
        acc_grid = core.Grid(5,5,5)
        acc_grid.ax = 1 | generic_unit_system.acceleration
        acc_grid.ay = 1 | generic_unit_system.acceleration
        acc_grid.az = 1 | generic_unit_system.acceleration
        #self.assertEquals(acc_grid.acceleration[0][0][0], ( 1,1,1) | generic_unit_system.acceleration)
        channel = acc_grid.new_channel_to(instance.acceleration_grid)
        channel.copy()
        
        result = instance.initialize_grid()
                   
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
    
    

    
