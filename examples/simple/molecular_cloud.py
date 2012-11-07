"""
Evolves a molecular cloud with explictly split gravtiy evolution

The gas is evolved in a sph code, and it's gravity is determined by a tree code.

Initial condition is a smooth spherical cloud with random velocities as in Bonnell et al. (2003)  
"""  

import numpy
  
from matplotlib import pyplot 

from amuse.units import nbody_system
from amuse.units import units

from amuse.community.fi.interface import Fi
from amuse.community.bhtree.interface import BHTree

from amuse.ext.molecular_cloud import molecular_cloud
from amuse.ext.evrard_test import body_centered_grid_unit_cube
from amuse.couple import bridge

def make_map(sph, N=100, grid_size= 1 | units.parsec):
    
    cell_center_positions_1d = numpy.linspace(-0.5, 0.5, N)
    x,y = numpy.meshgrid(cell_center_positions_1d, cell_center_positions_1d)
    x = x.flatten()
    y = y.flatten()
    
    x *= grid_size
    y *= grid_size
    z = x * 0.0
    
    vx=0.0 * x / (1.0 | units.s)
    vy=0.0 * x / (1.0 | units.s)
    vz=0.0 * x / (1.0 | units.s)
    
    rho,rhovx,rhovy,rhovz,rhoe=sph.get_hydro_state_at_point(x,y,z,vx,vy,vz)
    rho=rho.reshape((N,N))

    return rho
    
def run_mc(N=5000,Mcloud=10000. | units.MSun,Rcloud=1. | units.parsec):

    conv = nbody_system.nbody_to_si(Mcloud,Rcloud)

    interaction_timestep = 0.01 | units.Myr
    time_end=0.36 | units.Myr

    parts=molecular_cloud(targetN=N,convert_nbody=conv,
            base_grid=body_centered_grid_unit_cube).result

    sph=Fi(conv)

    # need to turn off self gravity, the bridge will calculate this
    sph.parameters.self_gravity_flag=False
    
    # some typical Fi flags (just for reference, most are default values)
    sph.parameters.use_hydro_flag=True
    sph.parameters.isothermal_flag=True
    sph.parameters.integrate_entropy_flag=False
    sph.parameters.gamma=1
    sph.parameters.verbosity = 0
    
    # setting the hydro timestep is important
    # the sph code will take 2 timesteps every interaction timestep
    sph.parameters.timestep=interaction_timestep/2  

    sph.gas_particles.add_particles(parts)

    def new_code_to_calculate_gravity_of_gas_particles():
        result = BHTree(conv)
        return result
        
    calculate_gravity_code=bridge.CalculateFieldForCodes(
        new_code_to_calculate_gravity_of_gas_particles,  # the code that calculates the gravity
        input_codes = [sph] # the codes to calculate the gravity for
    )
    
    bridged_system = bridge.Bridge()
    bridged_system.timestep=interaction_timestep
    bridged_system.add_system(
        sph,
        [calculate_gravity_code]
    )

    fig=pyplot.figure(figsize=(12,12))

    ncolumn=2
    nrow=2
    nplot=ncolumn*nrow

    grid_size = 3 | units.parsec
    extent = (grid_size * (-0.5, 0.5, -0.5, 0.5)).value_in(units.parsec)
    
    if nplot > 1:
        plot_timestep = time_end / (nplot - 1)
    else:
        plot_timestep = time_end
        
    for i in range(nplot):
        ttarget=i*plot_timestep
        print "evolving to time:", ttarget.as_quantity_in(units.Myr)
        bridged_system.evolve_model(ttarget) 
        
        rho=make_map(sph,N=200,grid_size=grid_size)
        subplot=fig.add_subplot(ncolumn,nrow,i+1)
        subplot.imshow(
            numpy.log10(1.e-5+rho.value_in(units.amu/units.cm**3)),
            extent = extent,
            vmin=1,
            vmax=5
        )
        subplot.set_title(ttarget.as_quantity_in(units.Myr))
    
    sph.stop()
    pyplot.show()
  
if __name__ in ("__main__","__plot__"):
    run_mc(10000, Mcloud=1000. | units.MSun, Rcloud=1. | units.parsec)
  
