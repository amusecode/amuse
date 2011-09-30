"""
Plummer gas model generator

This module contains the function new_plummer_gas_model, used to create Plummer 
models consisting of (SPH) gas particles.
"""

import numpy

from amuse.ext.evrard_test import uniform_unit_sphere
from amuse.units import nbody_system
from amuse.units import units
from amuse import datamodel

class _MakePlummerGasModel(object):
    
    def __init__(self, targetN, convert_nbody = None, base_grid=None, rscale=1/1.695,
                   mass_cutoff = 0.999, do_scale = False):
        self.targetN = targetN
        self.convert_nbody = convert_nbody
        self.rscale=rscale
        self.mass_frac=mass_cutoff
        self.do_scale = do_scale
        self.internal_energy = 0.25 / self.rscale
        self.base_sphere=uniform_unit_sphere(targetN,base_grid)   
           
    def new_model(self):
        x,y,z=self.base_sphere.make_xyz()
        self.actualN=len(x)
        r=numpy.sqrt(x**2+y**2+z**2)*self.mass_frac**(1/3.)
        rtarget=self.rscale*(r**2/(1-r**2))**.5
        mr=self.mass_frac**(1/3.)
        maxr=self.rscale*(mr**2/(1-mr**2))**.5        
        mass=numpy.ones_like(x)/self.actualN
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

        result.move_to_center()
        if self.do_scale:
            potential_energy = result.potential_energy(G = nbody_system.G)
            result.position *= potential_energy / (-0.5 | nbody_system.energy)
            
            internal_energy = result.thermal_energy()
            result.u *= ((0.25 | nbody_system.energy) / internal_energy)
        
        if not self.convert_nbody is None:
            result = datamodel.ParticlesWithUnitsConverted(result, self.convert_nbody.as_converter_from_si_to_nbody())
            result = result.copy_to_memory()
            
        return result


def new_plummer_gas_model(number_of_particles, *list_arguments, **keyword_arguments):
    """
    Create a plummer gas model with the given number of particles. Returns
    a set of SPH particles with equal masses and positions distributed
    to fit a plummer distribution model. The model is centered around the
    origin. Velocities are set to zero, and internal energies are set to 
    balance the gravitational forces between the gas particles.

    :argument number_of_particles: Number of particles to include in the plummer sphere
    :argument convert_nbody:  When given will convert the resulting set to SI units
    :argument mass_cutoff: Mass percentage inside radius of 1
    :argument do_scale: scale the result, similar to true nbody units (M=1, Q=0.25, U=-0.5)
    """
    uc = _MakePlummerGasModel(number_of_particles, *list_arguments, **keyword_arguments)
    return uc.result

