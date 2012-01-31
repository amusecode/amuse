"""
initial conditions for the SPH evrard collapse test
"""
import numpy

from math import *

from amuse.units import nbody_system
from amuse.units import units

from amuse.datamodel import Particles
from amuse.datamodel import ParticlesWithUnitsConverted
class uniform_random_unit_cube(object):
    def __init__(self,targetN):
        self.targetN=targetN
        self.par=long(targetN)
    def make_xyz(self):
        x=numpy.random.uniform(-1.,1.,self.par)
        y=numpy.random.uniform(-1.,1.,self.par)
        z=numpy.random.uniform(-1.,1.,self.par)
        return x,y,z

class regular_grid_unit_cube(object):
    def __init__(self,targetN):
        self.targetN=targetN
        self.par=long(float(targetN)**(1./3.)+1.5) 
    def make_xyz(self):
        nf=self.par
        x,y,z=numpy.mgrid[-1.:1.:nf*1j,-1.:1.:nf*1j,-1.:1.:nf*1j] 
        x=x.flatten()
        y=y.flatten()
        z=z.flatten()
        return x,y,z

class body_centered_grid_unit_cube(object):
    def __init__(self,targetN):
        self.targetN=targetN
        self.par=long(float(targetN/2.)**(1./3.)+1.5)   

    def make_xyz(self):
        nf=self.par
        x1,y1,z1=numpy.mgrid[-1.:1.:nf*1j,-1.:1.:nf*1j,-1.:1.:nf*1j] 
        x2,y2,z2=numpy.mgrid[-1.+1./2/nf:1.-1./2/nf:(nf-1)*1j,  
          -1.+1./2/nf:1.-1./2/nf:(nf-1)*1j,-1.+1./2/nf:1.-1./2/nf:(nf-1)*1j]                        
        x=numpy.concatenate( (x1.flatten(),x2.flatten()) )
        y=numpy.concatenate( (y1.flatten(),y2.flatten()) )
        z=numpy.concatenate( (z1.flatten(),z2.flatten()) )
        a=numpy.where((x>=-1) & (y>=-1) & (z>=-1) & (x<1) & (y<1) & (z<1) )
        return x[a],y[a],z[a]


class glass_unit_cube(object):
    def __init__(self,targetN,target_rms=0.01):
        self.targetN=targetN
        self.target_rms=target_rms
        if target_rms < 0.0001:
            print "warning: target_rms may not succeed"
        if targetN < 1000:
            print "warning: not enough particles"  
          
    def make_xyz(self):
        from amuse.community.fi.interface import Fi

        N=self.targetN
        target_rms=self.target_rms

        L=1| nbody_system.length
        dt=0.01 | nbody_system.time
        x,y,z=uniform_random_unit_cube(N).make_xyz()
        vx,vy,vz=uniform_unit_sphere(N).make_xyz()

        p=Particles(N)
        p.x=L*x
        p.y=L*y
        p.z=L*z
        p.h_smooth=0. * L
        p.vx= 0.1*vx | (nbody_system.speed)
        p.vy= 0.1*vy | (nbody_system.speed)
        p.vz= 0.1*vz | (nbody_system.speed)
        p.u= (0.1*0.1) | nbody_system.speed**2 
        p.mass=(8./N) | nbody_system.mass

        sph=Fi(use_gl=False,mode='periodic',redirection='none')   
        sph.initialize_code()

        sph.parameters.use_hydro_flag=True
        sph.parameters.radiation_flag=False
        sph.parameters.self_gravity_flag=False
        sph.parameters.gamma=1.
        sph.parameters.isothermal_flag=True
        sph.parameters.integrate_entropy_flag=False
        sph.parameters.timestep=dt  
        sph.parameters.verbosity=0
        sph.parameters.pboxsize=2*L
        sph.parameters.artificial_viscosity_alpha = 1.
        sph.parameters.beta = 2.
        sph.commit_parameters()
        sph.gas_particles.add_particles(p)
        sph.commit_particles()

#        sph.start_viewer()

        t=0. | nbody_system.time
        rms=1.
        minrms=1.
        i=0
        while rms > target_rms:
            i+=1
            t=t+(0.25 | nbody_system.time)
            sph.evolve_model(t)
            rho=sph.particles.rho.value_in(nbody_system.density)
            rms=rho.std()/rho.mean()
            minrms=min(minrms,rms)
            if rms>2.*minrms or i>300:
                print " RMS(rho) convergence warning:", i, rms,minrms
            if i>100000:
                print "i> 100k steps - not sure about this..."
                print " rms:", rms
                break


        x=sph.particles.x.value_in(nbody_system.length)
        y=sph.particles.y.value_in(nbody_system.length)
        z=sph.particles.z.value_in(nbody_system.length)

        del sph  
        return x,y,z

def uniform_unit_cube(targetN, base_grid=None):
    if base_grid is None:
        return body_centered_grid_unit_cube(targetN)
    else:
        return base_grid(targetN)
    
class uniform_unit_sphere(object):
    def __init__(self,targetN, base_grid=None):
        cube_sphere_ratio=4/3.*numpy.pi*0.5**3
        self.targetN=targetN
        self.estimatedN=targetN/cube_sphere_ratio
        if base_grid is None:
            self.base_grid=uniform_random_unit_cube
        else:
            self.base_grid=base_grid
   
    def cutout_sphere(self,x,y,z):
        r=x**2+y**2+z**2
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
                x,y,z=self.cutout_sphere(*(self.base_grid(estimatedN)).make_xyz())
            return x[0:self.targetN],y[0:self.targetN],z[0:self.targetN]  
        else:
            return self.cutout_sphere(*(self.base_grid(self.estimatedN)).make_xyz())
        
class MakeEvrardTest(object):
    def __init__(self, targetN, base_grid=None, size=1.,
                   mass=1.,internal_energy=0.05,seed=345672):
        numpy.random.seed(seed)
        self.targetN = targetN
        self.size=size
        self.mass=mass
        self.internal_energy=internal_energy
        self.base_sphere=uniform_unit_sphere(targetN,base_grid)   
           
    def new_model(self):
        x,y,z=self.base_sphere.make_xyz()
        self.actualN=len(x)
        r=numpy.sqrt(x**2+y**2+z**2)
        rtarget=self.size*r**1.5
        mass=numpy.ones_like(x)/self.actualN
        internal_energy=numpy.ones_like(x)*self.internal_energy
        r=r.clip(1.e-8,2*self.size)
        x=rtarget*x/r
        y=rtarget*y/r
        z=rtarget*z/r
        vx=numpy.zeros_like(x)
        vy=numpy.zeros_like(x)
        vz=numpy.zeros_like(x)
        return (mass,x,y,z,vx,vy,vz,internal_energy)
    

class MakeEvrardModel(object):
    
    def __init__(self, target_number_of_particles, convert_nbody = None, base_grid = None, 
            internal_energy = 0.05, do_scale = False, seed = None,size=1.):
        self.target_number_of_particles = target_number_of_particles
        self.convert_nbody = convert_nbody
        self.internal_energy = internal_energy
        self.size=size
        self.do_scale = do_scale
        self.base_sphere = uniform_unit_sphere(target_number_of_particles, base_grid)   
        numpy.random.seed(seed)
    
    def new_model(self):
        x, y, z = self.base_sphere.make_xyz()
        self.actual_number_of_particles = len(x)
        r = numpy.sqrt(x**2+y**2+z**2)
        rtarget = self.size*r**1.5
        mass = numpy.ones_like(x)/self.actual_number_of_particles
        internal_energy = numpy.ones_like(x)*self.internal_energy
        r = r.clip(1.e-8, 2.0*self.size)
        x = rtarget*x/r
        y = rtarget*y/r
        z = rtarget*z/r
        vx = numpy.zeros_like(x)
        vy = numpy.zeros_like(x)
        vz = numpy.zeros_like(x)
        return (mass, numpy.hstack((x, y, z)), numpy.hstack((vx, vy, vz)), internal_energy)
        
    @property
    def result(self):
        masses, positions, velocities, internal_energies = self.new_model()
        result = Particles(self.actual_number_of_particles)
        result.mass = nbody_system.mass.new_quantity(masses)
        result.position = nbody_system.length.new_quantity(positions)
        result.velocity = nbody_system.speed.new_quantity(velocities)
        result.u = nbody_system.specific_energy.new_quantity(internal_energies)
        
        result.position -= result.center_of_mass()
        if self.do_scale:
            scale_factor = (result.potential_energy(G=nbody_system.G)) / (-0.5 | nbody_system.energy)
            result.position *= scale_factor
        
        if not self.convert_nbody is None:
            result = ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_nbody())
            result = result.copy_to_memory()
            
        return result
    
"""
Create an evrard gas sphere with approximately the given number of particles. 
Returns a set of particles with equal mass and specific internal energy. 
Positions are randomly distributed to fit an evrard gas distribution model 
(density proportional to r^-1). Velocities are set to zero initially. The 
model is centered around the origin. Positions are optionally scaled such 
that the potential energy is -0.5 in nbody-units.

:argument target_number_of_particles: Target number of particles to include in the model
:argument convert_nbody:  When given will convert the resulting set to SI units
:argument internal_energy: The specific internal energy of each particle (defaults to 0.05)
:argument do_scale: scale the positions to exact nbody units (U=-0.5)
:argument seed:  Seed for the random number generator
"""
def new_evrard_gas_sphere(target_number_of_particles, *list_arguments, **keyword_arguments):
    uc = MakeEvrardModel(target_number_of_particles, *list_arguments, **keyword_arguments)
    return uc.result


if __name__=="__main__":
    x,y,z=uniform_unit_sphere(10000).make_xyz()
    print len(x)
