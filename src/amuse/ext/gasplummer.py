import numpy

from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.units import nbody_system
from amuse.units import units

from amuse import datamodel
class MakePlummerGasModel(object):
    def __init__(self, targetN, convert_nbody = None, base_grid=None, rscale=1/1.695,
                   mass=1.,seed=345672,mass_frac=.999):
        numpy.random.seed(seed)
        self.targetN = targetN
        self.convert_nbody = convert_nbody
        self.rscale=rscale
        self.mass=mass
        self.mass_frac=mass_frac
        self.internal_energy=0.25*self.mass/self.rscale
        self.base_sphere=uniform_unit_sphere(targetN,base_grid)   
           
    def new_model(self):
        x,y,z=self.base_sphere.make_xyz()
        self.actualN=len(x)
        r=numpy.sqrt(x**2+y**2+z**2)*self.mass_frac**(1/3.)
        rtarget=self.rscale*(r**2/(1-r**2))**.5
        mr=self.mass_frac**(1/3.)
        maxr=self.rscale*(mr**2/(1-mr**2))**.5        
        mass=numpy.ones_like(x)*self.mass/self.actualN
        internal_energy=self.internal_energy/(1+(rtarget/self.rscale)**2)**(1./2)
        r=r.clip(1.e-8,maxr)
        x=rtarget*x/r
        y=rtarget*y/r
        z=rtarget*z/r
        vx=numpy.zeros_like(x)
        vy=numpy.zeros_like(x)
        vz=numpy.zeros_like(x)
        return (mass,x,y,z,vx,vy,vz,internal_energy)
    
    @property
    def result(self):
        mass,x,y,z,vx,vy,vz,u = self.new_model()
        result = datamodel.Particles(self.actualN)
        result.mass = nbody_system.mass.new_quantity(mass)
        result.x = nbody_system.length.new_quantity(x)
        result.y = nbody_system.length.new_quantity(y)
        result.z = nbody_system.length.new_quantity(z)
        result.vx = nbody_system.speed.new_quantity(vx)
        result.vy = nbody_system.speed.new_quantity(vy)
        result.vz = nbody_system.speed.new_quantity(vz)
        result.u = (nbody_system.speed**2).new_quantity(u)

        if not self.convert_nbody is None:
            result = datamodel.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_nbody())
            result = result.copy_to_memory()
            
        return result

if __name__=="__main__":
    convert_nbody = nbody_system.nbody_to_si(100. | units.MSun, 1.0 | units.parsec)
    sphere=MakePlummerGasModel(10000,convert_nbody)
    parts=sphere.result
    print parts[0].internal_energy**0.5
    print len(parts)*parts[0].mass.in_(units.MSun)
    
