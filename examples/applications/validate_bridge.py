"""
first test of Fujii et al 2007

(not yet fully checked)

"""

import numpy 
import os

from amuse.io import store
from amuse.units import nbody_system
from amuse.units import units
from amuse.units import constants
from amuse.community.hermite0.interface import Hermite
from amuse.community.phiGRAPE.interface import PhiGRAPE
from amuse.community.bhtree.interface import BHTree
from amuse.community.fi.interface import Fi
from amuse.ext.bridge import bridge

from amuse.ext.derived_grav_systems import copycat
from amuse.ext.kingmodel import new_king_model

def sys_from_parts(base_class,parts,converter,eps=None,timestep=None,usegl=False,mode=None):
    if mode is None:
        if usegl:
            interface=base_class(converter,use_gl=True)
        else:
            interface=base_class(converter)
    else:
        if usegl:
            interface=base_class(converter,use_gl=True,mode=mode)
        else:
            interface=base_class(converter,mode=mode)
    
    interface.initialize_code()
    if eps is not None:
        interface.parameters.epsilon_squared = eps**2 
    if timestep is not None:
        interface.parameters.timestep = timestep 
    interface.particles.add_particles(parts)
    return interface

def sys_from_king(base_class,N,W, converter, eps=None,\
                            timestep=None,usegl=False,mode=None):
    parts = new_king_model(N, W0=W, convert_nbody=converter)
    parts.radius=0.| nbody_system.length
    return sys_from_parts(base_class,parts,converter,eps,timestep,usegl,mode)

def shift_sys(system,dx,dy,dz,dvx,dvy,dvz):
    parts=system.particles.copy()
    parts.x=parts.x+dx
    parts.y=parts.y+dy
    parts.z=parts.z+dz
    parts.vx=parts.vx+dvx
    parts.vy=parts.vy+dvy
    parts.vz=parts.vz+dvz
    parts.copy_values_of_all_attributes_to(system.particles)

def validate_bridge():
   
    convert_clu = nbody_system.nbody_to_si(1.e5 | units.MSun, 1.0 | units.parsec)
    convert_gal = nbody_system.nbody_to_si(1.e7 | units.MSun, 10. | units.parsec)

    convert_gal_to_clu=lambda x: convert_clu.to_nbody(convert_gal.to_si(x))
    convert_clu_to_gal=lambda x: convert_gal.to_nbody(convert_clu.to_si(x))

    total_sim_time=convert_clu_to_gal(200 | nbody_system.time)
    total_sim_time=200 | nbody_system.time

    sim_step=1./128 | nbody_system.time
    dt_diag=1./64 | nbody_system.time    

    eps_clu= convert_gal_to_clu(2.e-4 | nbody_system.length)
    eps_gal= 6.25e-3 | nbody_system.length

    cluster=sys_from_king(PhiGRAPE, 1000,7, convert_clu, eps_clu,usegl=False,mode="gpu")
    galaxy=sys_from_king(Fi,10000,9,convert_gal,eps_gal,sim_step,usegl=False)

    shift_sys(cluster,
      convert_clu.to_si( convert_gal_to_clu(2.5 | nbody_system.length) ),0 | units.m,0 | units.m,
      0| units.m/units.s,convert_clu.to_si( convert_gal_to_clu(0.65 | nbody_system.speed) ),0| units.m/units.s )
    
    if hasattr(cluster,"start_viewer"): cluster.start_viewer()
    if hasattr(galaxy,"viewer"): galaxy.viewer()

    clustercopy=copycat(BHTree,cluster,convert_clu)
      
    bridgesys=bridge(verbose=False)
    bridgesys.add_system(galaxy, (clustercopy,), False)
    bridgesys.add_system(cluster, (galaxy,), True )

    t=0. | nbody_system.time

    if os.path.exists("cluster.hdf"):
        os.remove("cluster.hdf")
    storage=store.StoreHDF("cluster.hdf")

    bridgesys.synchronize_model()
    Ep0=bridgesys.potential_energy
    Ek0=bridgesys.kinetic_energy
    Esp0=cluster.kinetic_energy
    Esk0=cluster.potential_energy
    tc=cluster.model_time
    tg=galaxy.model_time

    print convert_gal.to_nbody(tg)
    print 0.
#  print Ep0,Ek0
#  print Esp0,Esk0
    for x in cluster.get_center_of_mass_position():
        print convert_gal.to_nbody(x),
    print    

    part=bridgesys.particles.copy()
    part.savepoint(tg)
    storage.store(part.previous_state())    

    while( t < total_sim_time):
        t=t+dt_diag
        bridgesys.evolve_model( convert_gal.to_si(t),
                        timestep=convert_gal.to_si( sim_step ))     

        bridgesys.synchronize_model()
        Ep=bridgesys.potential_energy
        Ek=bridgesys.kinetic_energy
        Esp=cluster.kinetic_energy
        Esk=cluster.potential_energy
        tc=cluster.model_time
        tg=galaxy.model_time
        print convert_gal.to_nbody(tg)
        print (Ep+Ek-Ep0-Ek0)/(Ep0+Ek0)
#    print (Esp+Esk-Esp0-Esk0)/(Esp0+Esk0)
        for x in cluster.get_center_of_mass_position():
            print convert_gal.to_nbody(x),
        print  
        
        part=bridgesys.particles.copy()
        part.savepoint(tg)
        storage.store(part.previous_state())  
  
if __name__ == '__main__':
    validate_bridge()
