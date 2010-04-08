import os
import sys
import numpy

from legacy_support import TestWithMPI
from amuse.legacy.capreole.interface import capreole
from amuse.legacy.support.channel import MessageChannel


#MessageChannel.no_redirection()

class TestMPIInterface(TestWithMPI):
    
    def test0(self):
        instance=capreole()
        instance.setup_module()
        instance.stop()
    
    def test1(self):
        instance=capreole()
        instance.setup_module()
        instance.setup_mesh(50,50,50,1.,1.,1.)
        instance.stop()
    
    def xtest2(self):
        instance=capreole()
        instance.setup_module()
        instance.setup_mesh(50,50,50,1.,1.,1.)
        instance.fill_grid_state(1,1,1,1.,0.,0.,0.,1.)
        rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(1,1,1)
        self.assertEqual(rho,1.)
        self.assertEqual(rhovx,0.)
        self.assertEqual(rhovy,0.)
        self.assertEqual(rhovz,0.)
        self.assertEqual(en,1.)
        instance.stop()
    
    def test3(self):
        instance=capreole()
        instance.setup_module()
        instance.setup_mesh(50,50,50,1.,1.,1.)
        x,y,z=numpy.indices( (50,50,50) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=0.25*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=2.*numpy.ones_like(x)
        instance.fill_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        rho,rhovx,rhovy,rhovz,en,err=instance.get_grid_state(1,1,1)
        self.assertEqual(rho,0.25)
        self.assertEqual(rhovx,0.)
        self.assertEqual(rhovy,0.)
        self.assertEqual(rhovz,0.)
        self.assertEqual(en,2.)
        instance.stop()
    
    def test4(self):
        instance=capreole()
        instance.setup_module()
        instance.setup_mesh(50,50,50,1.,1.,1.)
        x,y,z=numpy.indices( (50,50,50) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=0.1*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=0.1*numpy.ones_like(x)
        instance.fill_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        instance.initialize_grid(0.0)
        instance.stop()
    
    def test5(self):
        instance=capreole()
        instance.setup_module()
        instance.setup_mesh(40,40,40,1.,1.,1.)
        x,y,z=numpy.indices( (40,40,40) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=0.1*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=0.1*numpy.ones_like(x)
        instance.fill_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        instance.initialize_grid(0.0)
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
        instance=capreole(number_of_workers=1)
        instance.setup_module()
        instance.setup_mesh(150,10,30,1.,1.,1.)
        x,y,z=numpy.indices( (150,10,30) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=0.1*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=0.1*numpy.ones_like(x)
        instance.fill_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        instance.initialize_grid(0.0)
        instance.evolve(0.01)
        x,y,z,err=instance.get_position_of_index(15,5,20)
        self.assertAlmostEqual(x,15/150.-1/300.,15)
        self.assertAlmostEqual(y,5/10.-1/20.,15)
        self.assertAlmostEqual(z,20/30.-1/60.,15)
        i,j,k,err=instance.get_index_of_position(x,y,z)
        self.assertEqual([i,j,k],[15,5,20])
        instance.stop()

class TestSodShocktube(TestWithMPI):
    
    def test0(self):
        N=100
        gamma=5/3.
        g=(gamma-1)/(gamma+1)
        b=(gamma-1)/2/gamma
        
        instance=capreole()
        instance.setup_module()
        instance.setup_mesh(N,N/10,N/10,1.,0.1,0.1)
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
        instance.fill_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)

        x,y,z=numpy.indices( (N/2,N/10,N/10) )
        x=x.flatten()+1
        y=y.flatten()+1
        z=z.flatten()+1
        rho=1.*numpy.ones_like(x)
        rhvx=0.*numpy.ones_like(x)
        rhvy=0.*numpy.ones_like(x)
        rhvz=0.*numpy.ones_like(x)
        en=(1./(gamma-1))*numpy.ones_like(x)
        instance.fill_grid_state(x,y,z,rho,rhvx,rhvy,rhvz,en)
        instance.initialize_grid(0.0)
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
    
