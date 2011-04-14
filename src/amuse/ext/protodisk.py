import numpy

from matplotlib import pyplot

from amuse.support.units import units
from amuse.support.units import nbody_system
from amuse.support.data.core import Particles, ParticlesWithUnitsConverted

from amuse.ext.evrard_test import uniform_random_unit_cube,regular_grid_unit_cube,body_centered_grid_unit_cube

def approximate_inverse_error_function(x):
  a=8*(numpy.pi-3)/3*numpy.pi*(4-numpy.pi)
  return numpy.sign(x)*numpy.sqrt(
    numpy.sqrt((2/numpy.pi/a+numpy.log(1-x**2)/2)**2-numpy.log(1-x**2)/a)-(2/numpy.pi/a+numpy.log(1-x**2)/2)
   )

class uniform_unit_cylinder(object):
    def __init__(self,targetN, base_grid=None):
        cube_cylinder_ratio=numpy.pi*0.5**2
        self.targetN=targetN
        self.estimatedN=targetN/cube_cylinder_ratio
        if base_grid is None:
            self.base_grid=uniform_random_unit_cube
        else:
            self.base_grid=base_grid
   
    def cutout_cylinder(self,x,y,z):
        r=x**2+y**2
        selection=r < numpy.ones_like(r)        
        x=x.compress(selection)
        y=y.compress(selection)
        z=z.compress(selection)
        return x,y,z

    def make_xyz(self):
        if(self.base_grid==uniform_random_unit_cube):
            estimatedN=self.estimatedN
            x=[]
            while len(x) < self.targetN:
                estimadedN=estimatedN*1.1+1
                x,y,z=self.cutout_cylinder(*(self.base_grid(estimatedN)).make_xyz())
            return x[0:self.targetN],y[0:self.targetN],z[0:self.targetN]  
        else:
            return self.cutout_cylinder(*(self.base_grid(self.estimatedN)).make_xyz())


class ProtoPlanetaryDisk(object):

    def __init__(self, targetN, convert_nbody = None, discfraction=0.1,
                   densitypower=1., thermalpower=0.5, Rmin=1,Rmax=100,
                   gamma=1.,q_out=2.,base_grid=None):
        self.targetN=targetN
        self.convert_nbody=convert_nbody
        self.densitypower=densitypower
        self.thermalpower=thermalpower
        self.Rmin=Rmin
        self.Rmax=Rmax
        self.gamma=gamma
        self.q_out=q_out
        self.discfraction=discfraction

        self.a=self.thermalpower
        self.a2=self.thermalpower/2
        self.g=densitypower
        self.g2=2-densitypower
        self.k_out=((1+discfraction)/Rmax**3)**0.5
        self.sigma_out=self.g2*discfraction/(2*numpy.pi*Rmax**self.g*(Rmax**self.g2-Rmin**self.g2))
        self.cs_out=self.q_out*numpy.pi*self.sigma_out/self.k_out
        
        self.base_cylinder=uniform_unit_cylinder(targetN,base_grid)   


    def sigma(self,r):
        return self.sigma_out*(self.Rmax/r)**self.g

    def csound(self,r):
        return self.cs_out*(self.Rmax/r)**self.a2

    def cmass(self,r):
        return self.discfraction*(r**self.g2-self.Rmin**self.g2)/(self.Rmax**self.g2-self.Rmin**self.g2)

    def mass_encl(self,r):
        return 1+self.cmass(r)

    def kappa(self,r):
        return (self.mass_encl(r)/r**3)**0.5
  
    def toomreQ(self,r):
        return self.csound(r)*self.kappa(r)/numpy.pi/self.sigma(r)

    def getradius(self,f):
        return ((self.Rmax**self.g2-self.Rmin**self.g2)*f+self.Rmin**self.g2)**(1./self.g2)

    def zscale(self,r):
        return self.csound(r)/self.kappa(r) 

    def u(self,r):
        if self.gamma ==1.:
            return self.csound(r)**2
        else:
            return self.csound(r)**2/(self.gamma-1)

    def vcirc(self,r):
        return (self.mass_encl(r)/r)**0.5
    
    def new_model(self):
        x,y,z=self.base_cylinder.make_xyz()
        self.actualN=len(x)
        f=x**2+y**2
        r=f**0.5
        rtarget=self.getradius(f)
        
        mass=self.discfraction*numpy.ones_like(x)/self.actualN
        internal_energy=self.u(rtarget)
        zscale=self.zscale(rtarget)
        r=r.clip(1.e-8,2.)
        x=x/r
        y=y/r

        vx=-y*self.vcirc(rtarget)
        vy=x*self.vcirc(rtarget)
        vz=numpy.zeros_like(x)
         
        x=rtarget*x 
        y=rtarget*y                  
        z=approximate_inverse_error_function(z)*zscale*2.**0.5

        return (mass, x, y, z, vx, vy, vz, internal_energy)
  
    @property
    def result(self):
        masses, x,y,z, vx,vy,vz, internal_energies = self.new_model()
        result = Particles(self.actualN)
        result.mass = nbody_system.mass.new_quantity(masses)
        result.x = nbody_system.length.new_quantity(x)
        result.y = nbody_system.length.new_quantity(y)
        result.z = nbody_system.length.new_quantity(z)
        result.vx = nbody_system.speed.new_quantity(vx)
        result.vy = nbody_system.speed.new_quantity(vy)
        result.vz = nbody_system.speed.new_quantity(vz)
        result.u = nbody_system.specific_energy.new_quantity(internal_energies)
                
        if not self.convert_nbody is None:
            result = ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_nbody())
            result = result.copy_to_memory()
            
        return result

