"""
initial conditions for the SPH evrard collapse test
"""
from math import *
import numpy

class MakeGridSphere(object):
    def __init__(self, number_of_particles,size=1.):
        self.target_number_of_particles = number_of_particles
        self.size=size

    def new_model(self):
        cube_sphere_ratio=4/3.*numpy.pi*0.5**3
        nf=(self.target_number_of_particles/cube_sphere_ratio)**(1/3.)
        nf=long(nf+0.5)
        
        x,y,z=numpy.mgrid[-1.:1.:nf*1j,-1.:1.:nf*1j,-1.:1.:nf*1j] 
        x=x.flatten()
        y=y.flatten()
        z=z.flatten()
        r=x**2+y**2+z**2
        selection=r < numpy.ones_like(r)
        x=self.size*x.compress(selection)
        y=self.size*y.compress(selection)
        z=self.size*z.compress(selection)
        self.actual_number_of_particles=len(x)
        return x,y,z

class MakeRandomSphere(object):
    def __init__(self, number_of_particles,size=1.,seed=634567):
        self.target_number_of_particles = number_of_particles
        self.size=size
        numpy.random.seed(seed)

    def new_model(self):
        cube_sphere_ratio=4/3.*numpy.pi*0.5**3
        
        x=numpy.random.uniform(-1.,1.,long(self.target_number_of_particles/cube_sphere_ratio))
        y=numpy.random.uniform(-1.,1.,long(self.target_number_of_particles/cube_sphere_ratio))
        z=numpy.random.uniform(-1.,1.,long(self.target_number_of_particles/cube_sphere_ratio))
        r=x**2+y**2+z**2
        selection=r < numpy.ones_like(r)
        x=self.size*x.compress(selection)
        y=self.size*y.compress(selection)
        z=self.size*z.compress(selection)
        self.actual_number_of_particles=len(x)
        return x,y,z
        
class MakeEvrardTest(object):
    def __init__(self, number_of_particles, grid=True, size=1.,
                   mass=1.,internal_energy=0.05,seed=345672):
        self.target_number_of_particles = number_of_particles
        self.size=size
        self.mass=mass
        self.internal_energy=internal_energy
        if(grid):
          self.sphere=MakeGridSphere(self.target_number_of_particles,
                                       size=self.size)
        else:
          self.sphere=MakeRandomSphere(self.target_number_of_particles,
                                         size=self.size,seed=seed)
           
    def new_model(self):
        x,y,z=self.sphere.new_model()
        r=numpy.sqrt(x**2+y**2+z**2)
        rtarget=self.size*r**1.5
        self.actual_number_of_particles=len(x)
        mass=numpy.ones_like(x)/self.actual_number_of_particles
        internal_energy=numpy.ones_like(x)*self.internal_energy
        r=r.clip(1.e-8,2*self.size)
        x=rtarget*x/r
        y=rtarget*y/r
        z=rtarget*z/r
        vx=numpy.zeros_like(x)
        vy=numpy.zeros_like(x)
        vz=numpy.zeros_like(x)
        return (mass,x,y,z,vx,vy,vz,internal_energy)

if __name__=="__main__":
  sphere=MakeEvrardTest(10000,grid=False)
  mass,x,y,z,vx,vy,vz,u=sphere.new_model()
  print mass[0],x[0],vx[0],u[0]
  print len(mass),len(x),len(vx),len(u)
  
